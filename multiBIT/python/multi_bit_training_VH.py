#!/usr/bin/env python

# Standard imports
import ROOT
import numpy as np
import random
import cProfile
import time
import os, sys
from math import log, exp, sin, cos, sqrt, pi
import copy
import pickle
import itertools

# RootTools
from RootTools.core.standard   import *
ROOT.gROOT.LoadMacro("$CMSSW_BASE/src/Analysis/Tools/scripts/tdrstyle.C")
ROOT.setTDRStyle()

import helpers

# Analysis
import Analysis.Tools.syncer as syncer

# BIT
# Boosting
import sys, os, time
sys.path.insert(0,os.path.expandvars("$CMSSW_BASE/src/BIT"))
from MultiBoostedInformationTree import MultiBoostedInformationTree

# User
import TMB.Tools.user as user

# Parser
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument("--plot_directory",     action="store",      default="multiBIT_VH",                 help="plot sub-directory")
argParser.add_argument('--config_module',      action='store', type=str, default = "TMB.multiBIT.configs", help = "config directory")
argParser.add_argument('--config',             action='store', type=str, default = "WH_delphes", help="config")
argParser.add_argument('--name',               action='store', type=str,   default='default', help="Name of the training")
argParser.add_argument('--output_directory',   action='store', type=str,   default=os.path.expandvars('/groups/hephy/cms/$USER/BIT/'))
argParser.add_argument('--input_directory',    action='store', type=str,   default=os.path.expandvars("/groups/hephy/cms/$USER/multiBIT/training-ntuple-WH"))
argParser.add_argument("--coefficients",       action="store",      default=None,       nargs="*", help="Which coefficients?")
argParser.add_argument('--variable_set',       action='store', type=str,   default='mva_variables', help="List of variables for training")
argParser.add_argument("--plot_iterations",    action="store",      default=None,          nargs="*", type=int, help="Certain iterations to plot? If first iteration is -1, plot only list provided.")
argParser.add_argument('--overwrite',          action='store_true', help="overwrite training?")
argParser.add_argument('--bias',               action='store',      default=None, nargs = "*",  help="Bias training? Example:  --bias 'pT' '10**(({}-200)/200) ")
argParser.add_argument('--small',              action='store_true', help="small?")
argParser.add_argument('--debug',              action='store_true', help="Make debug plots?")
argParser.add_argument('--feature_plots',      action='store_true', help="Feature plots?")
argParser.add_argument('--lumi_norm',          action='store_true', help="Normalize the events according to lumiweight1fb?")
argParser.add_argument('--maxEvents',          action='store', default = None, type=int, help="Maximum number of training events")

args, extra = argParser.parse_known_args(sys.argv[1:])

def parse_value( s ):
    try:
        r = int( s )
    except ValueError:
        try:
            r = float(s)
        except ValueError:
            r = s
    return r

extra_args = {}
key        = None
for arg in extra:
    if arg.startswith('--'):
        # previous no value? -> Interpret as flag
        #if key is not None and extra_args[key] is None:
        #    extra_args[key]=True
        key = arg.lstrip('-')
        extra_args[key] = True # without values, interpret as flag
        continue
    else:
        if type(extra_args[key])==type([]):
            extra_args[key].append( parse_value(arg) )
        else:
            extra_args[key] = [parse_value(arg)]
for key, val in extra_args.iteritems():
    if type(val)==type([]) and len(val)==1:
        extra_args[key]=val[0]

#Logger
import Analysis.Tools.logger as logger
logger = logger.get_logger("INFO", logFile = None )

# MVA configuration
import importlib
configs = importlib.import_module(args.config_module)
config  = getattr( configs, args.config)
config.multi_bit_cfg.update( extra_args )

import uproot
import awkward
import numpy as np
import pandas as pd
#import h5py

#########################################################################################
if args.small:
    args.name+='_small'

# directories
plot_directory   = os.path.join( user. plot_directory, 'multiBit', args.config, args.name)
output_directory = os.path.join( args.output_directory, 'models', args.config, args.name)
# saving
if not os.path.exists(output_directory):
    try:
        os.makedirs(output_directory)
    except OSError:
        pass

# fix random seed for reproducibility
np.random.seed(1)

# get the training variable names
feature_names = [ mva_variable[0] for mva_variable in getattr(config, args.variable_set) ]

n_var_flat   = len(feature_names)

training_features = {}
training_weights = {}
lumi_weights       = {}
for i_training_sample, training_sample in enumerate(config.training_samples):
    upfile_name = os.path.join(os.path.expandvars(args.input_directory), training_sample.name, training_sample.name+'.root')
    logger.info( "Loading upfile %i: %s from %s", i_training_sample, training_sample.name, upfile_name)
    upfile = uproot.open(upfile_name)
    training_features[training_sample.name]            = upfile["Events"].pandas.df(branches = feature_names )
    training_weights[training_sample.name]  = upfile["Events"].pandas.df(branches = ["weight_derivatives"] )
    if args.lumi_norm:
        lumi_weights[training_sample.name]    = upfile["Events"].pandas.df(branches = ["lumiweight1fb"] )

training_features     = pd.concat([training_features[training_sample.name] for training_sample in config.training_samples])

## Quick fix in training
len_features = len(training_features)
if args.config == "WH_nlo_delphes":
    map_ = training_features.mva_WH_deltaR<5.5
else:
    map_ = np.ones(len(training_features)).astype('bool')

training_features = training_features[map_].values

training_weights = pd.concat([training_weights[training_sample.name] for training_sample in config.training_samples])
if args.lumi_norm:
    lumi_weights = pd.concat([lumi_weights[training_sample.name] for training_sample in config.training_samples])

training_weights = training_weights.values.reshape((len_features,-1))[map_][:args.maxEvents]
training_features           = training_features[:args.maxEvents]
if args.lumi_norm:
    lumi_weights       = lumi_weights.values[:,0][map_][:args.maxEvents]
    lumi_weights       = lumi_weights/np.mean(lumi_weights)

# split dataset into Input and output data

# number of samples with 'small'
n_small_samples = 10000
# small
if args.small:
    training_features = training_features[:n_small_samples]
    training_weights = training_weights[:n_small_samples]
    if args.lumi_norm:
        lumi_weights       = lumi_weights[:n_small_samples]

training_features  = training_features[:,0:n_var_flat]

if args.lumi_norm:
    training_weights = training_weights * lumi_weights[:,np.newaxis]

## Clip the most extreme scores 
#for derivative in config.combinations:
#
#    if config.bit_cfg[derivative]['clip_score_quantile'] is None: continue
#    i_derivative          = config.weight_derivative_combinations.index(derivative)
#    clip_score_quantile = config.bit_cfg[derivative]['clip_score_quantile']
#
#    clip_both_sides = True #(len(derivative)%2==0)
#
#    derivative = config.weight_derivative_combinations[i_derivative]
#    training_scores          = np.divide( weight_derivatives[:,i_derivative], weight_derivatives[:,0], np.zeros_like(weight_derivatives[:,i_derivative]), where=weight_derivatives[:,0]!=0 )
#    training_score_quantiles = np.quantile( training_scores, [clip_score_quantile, 1-clip_score_quantile] if clip_both_sides else [ 1-clip_score_quantile] )
#
#    ## Always clip pos & neg because we're dealing with nlo
#    #if len(derivative)%2==0:
#    if False: #clip_both_sides:
#        # quadratic: Clip only positive values
#        print "Clip score at percentile %3.2f for %r: s>%3.2f" % (1-clip_score_quantile, derivative, training_score_quantiles[0] )
#        to_clip                         = training_scores>training_score_quantiles[0]
#        weight_derivatives[:,i_derivative][to_clip] = training_score_quantiles[0]*weight_derivatives[:,0][to_clip]
#
#    else:
#        # linear: Clip positive and negative values
#        print "Clip score at percentiles %3.2f and %3.2f for %r: s<%3.2f and s>%3.2f" % (clip_score_quantile, 1-clip_score_quantile, derivative, training_score_quantiles[0], training_score_quantiles[1] )
#        to_clip                  = training_scores>training_score_quantiles[1]
#        weight_derivatives[:,i_derivative][to_clip] = training_score_quantiles[1]*weight_derivatives[:,0][to_clip]
#        to_clip                  = training_scores<training_score_quantiles[0]
#        weight_derivatives[:,i_derivative][to_clip] = training_score_quantiles[0]*weight_derivatives[:,0][to_clip]

# dictionary format
training_weights = {comb:training_weights[:,i_comb] for i_comb, comb in enumerate(config.combinations)}

###############
## Plot Model #
###############

if args.feature_plots and hasattr( config, "eft_plot_points"):
    h    = {}
    h_lin= {}
    h_rw = {}
    h_rw_lin = {}
    for i_eft, eft_plot_point in enumerate(config.eft_plot_points):
        eft = eft_plot_point['eft']

        if i_eft == 0:
            eft_sm     = eft
        name = ''
        name= '_'.join( [ (wc+'_%3.2f'%eft[wc]).replace('.','p').replace('-','m') for wc in config.coefficients if eft.has_key(wc) ])
        tex_name = eft_plot_point['tex'] 

        if i_eft==0: name='SM'
        h[name] = {}
        eft['name']=name
        
        for i_feature, feature in enumerate(feature_names):
            h[name][feature]        = ROOT.TH1F(name+'_'+feature+'_nom',    name+'_'+feature, *config.plot_options[feature]['binning'] )

        # make reweights for x-check
        reweight = copy.deepcopy(training_weights[()])
        # linear term
        for param1 in config.coefficients:
            reweight += (eft[param1]-eft_sm[param1])*training_weights[(param1,)] 
        # quadratic term
        for param1 in config.coefficients:
            if eft[param1]-eft_sm[param1] ==0: continue
            for param2 in config.coefficients:
                if eft[param2]-eft_sm[param2] ==0: continue
                reweight += .5*(eft[param1]-eft_sm[param1])*(eft[param2]-eft_sm[param2])*training_weights[tuple(sorted((param1,param2)))]

        sign_postfix = ""
        if False:
            reweight_sign = np.sign(np.sin(2*np.arccos(training_features[:,feature_names.index('cos_theta')]))*np.sin(2*np.arccos(training_features[:,feature_names.index('cos_theta_hat')])))
            reweight     *= reweight_sign
            #reweight_lin_sign = reweight_sign*reweight_lin
            sign_postfix    = " weighted with sgn(sin(2#theta)sin(2#hat{#theta}))"

        for i_feature, feature in enumerate(feature_names):
            binning = config.plot_options[feature]['binning']

            h[name][feature] = helpers.make_TH1F( np.histogram(training_features[:,i_feature], np.linspace(binning[1], binning[2], binning[0]+1), weights=reweight) )
            h[name][feature].style      = styles.lineStyle( eft_plot_point['color'], width=2, dashed=False )
            h[name][feature].legendText = tex_name

    for i_feature, feature in enumerate(feature_names):

        for _h in [h]:
            norm = _h[config.eft_plot_points[0]['eft']['name']][feature].Integral()
            if norm>0:
                for eft_plot_point in config.eft_plot_points:
                    _h[eft_plot_point['eft']['name']][feature].Scale(1./norm) 

        histos = [[h[eft_plot_point['eft']['name']][feature]] for eft_plot_point in reversed(config.eft_plot_points)]
        plot   = Plot.fromHisto( feature+'_nom',  histos, texX=config.plot_options[feature]['tex'], texY="1/#sigma_{SM}d#sigma/d%s"%config.plot_options[feature]['tex'] )

        for log in [True, False]:

            # Add subdirectory for lin/log plots
            plot_directory_ = os.path.join( plot_directory, "feature_plots", "log" if log else "lin" )
            for p in [plot] :#, plot_lin, plot_rw_lin]:
                plotting.draw( p,
                               plot_directory = plot_directory_,
                               logX = False, logY = log, sorting = False,
                               yRange = "auto" if not log else (0.002,"auto"),
                               ratio = None,
        #                       drawObjects = drawObjects( lumi, offset=titleOffset ),
                                legend=[(0.2,0.68,0.9,0.91),2],
                               #histModifications = histModifications,
                               copyIndexPHP = True,
                               )
print ("Done with plots")
syncer.sync()

# if not specified, take all
if args.coefficients is None:
    args.coefficients = config.coefficients

# delete coefficients we don't need (the BIT coefficients are determined from the training weight keys)
for key in training_weights.keys():
    if not all( [k in args.coefficients for k in key]):
        del training_weights[key]

bit_name = "multiBit_%s_%s_%s_nTrees_%i"%(args.config, args.name, "_".join(args.coefficients), config.multi_bit_cfg["n_trees"])

filename = os.path.join(user.model_directory, bit_name)+'.pkl'
try:
    print ("Loading %s for %s"%(bit_name, filename))
    bit = MultiBoostedInformationTree.load(filename)
except IOError:
    bit = None

# reweight training data according to bias
if args.bias is not None:
    if len(args.bias)!=2: raise RuntimeError ("Bias is defined by <var> <function>, i.e. 'x' '10**(({}-200)/200). Got instead %r"%args.bias)
    function     = eval( 'lambda x:'+args.bias[1].replace('{}','x') ) 
    bias_weights = np.array(map( function, training_features[:, feature_names.index(args.bias[0])] ))
    bias_weights /= np.mean(bias_weights)
    training_weights = {k:v*bias_weights for k,v in training_weights.iteritems()} 

# construct base_points
base_points = []
for comb in list(itertools.combinations_with_replacement(args.coefficients,1))+list(itertools.combinations_with_replacement(args.coefficients,2)):
    base_points.append( {c:comb.count(c) for c in args.coefficients} )

if bit is None or args.overwrite:
    time1 = time.time()
    bit = MultiBoostedInformationTree(
            training_features     = training_features,
            training_weights      = training_weights,
            base_points           = base_points,
            feature_names         = feature_names,
            **config.multi_bit_cfg
                )
    bit.boost()
    bit.save(filename)
    print ("Written %s"%( filename ))

    time2 = time.time()
    boosting_time = time2 - time1
    print ("Boosting time: %.2f seconds" % boosting_time)

if args.debug:

    # Loss plot
    training_losses = helpers.make_TH1F((bit.losses(training_features, training_weights),None), ignore_binning = True)
    #test_losses     = helpers.make_TH1F((bit.losses(test_features, test_weights),None),         ignore_binning = True)

    c1 = ROOT.TCanvas("c1");

    l = ROOT.TLegend(0.2,0.8,0.9,0.85)
    l.SetNColumns(2)
    l.SetFillStyle(0)
    l.SetShadowColor(ROOT.kWhite)
    l.SetBorderSize(0)


    training_losses.GetXaxis().SetTitle("N_{B}")
    training_losses.GetYaxis().SetTitle("Loss")
    l.AddEntry(training_losses, "train")
    #l.AddEntry(test_losses, "test")


    #test_losses.SetLineWidth(2)
    #test_losses.SetLineColor(ROOT.kRed+2)
    #test_losses.SetMarkerColor(ROOT.kRed+2)
    #test_losses.SetMarkerStyle(0)
    training_losses.SetLineWidth(2)
    training_losses.SetLineColor(ROOT.kRed+2)
    training_losses.SetMarkerColor(ROOT.kRed+2)
    training_losses.SetMarkerStyle(0)

    training_losses.Draw("hist") 
    #test_losses.Draw("histsame")

    for logY in [True, False]:
        plot_directory_ = os.path.join( plot_directory, "training_plots", bit_name, "log" if logY else "lin" )
        c1.Print(os.path.join(plot_directory_, "loss.png"))

    # GIF animation
    tex = ROOT.TLatex()
    tex.SetNDC()
    tex.SetTextSize(0.06)

    # Which iterations to plot
    plot_iterations = range(1,10)+range(10,bit.n_trees+1,10)
    # Add plot iterations from command line, if provided
    if type(args.plot_iterations)==type([]):
        if args.plot_iterations[0]<0:
            plot_iterations+=args.plot_iterations[1:]
        else:
            plot_iterations = args.plot_iterations
        plot_iterations.sort()

    for max_n_tree in plot_iterations:
        if max_n_tree==0: max_n_tree=1
        stuff = []
        training_predictions = bit.vectorized_predict(training_features, max_n_tree = max_n_tree)

        # colors
        color = {}
        i_lin, i_diag, i_mixed = 0,0,0
        for i_der, der in enumerate(bit.derivatives):
            if len(der)==1:
                color[der] = ROOT.kAzure + i_lin
                i_lin+=1
            elif len(der)==2 and len(set(der))==1:
                color[der] = ROOT.kRed + i_diag
                i_diag+=1
            elif len(der)==2 and len(set(der))==2:
                color[der] = ROOT.kGreen + i_mixed
                i_mixed+=1

        w0 = training_weights[()]
        h_w0, h_ratio_prediction, h_ratio_truth, lin_binning = {}, {}, {}, {}
        for i_feature, feature in enumerate(feature_names):
            # root style binning
            binning     = config.plot_options[feature]['binning']
            # linspace binning
            lin_binning[feature] = np.linspace(binning[1], binning[2], binning[0]+1)
            #digitize feature
            binned      = np.digitize(training_features[:,i_feature], lin_binning[feature] )
            # for each digit, create a mask to select the corresponding event in the bin (e.g. training_features[mask[0]] selects features in the first bin
            mask        = np.transpose( binned.reshape(-1,1)==range(1,len(lin_binning[feature])) )

            h_w0[feature]           = np.array([  w0[m].sum() for m in mask])
            h_derivative_prediction = np.array([ (w0.reshape(-1,1)*training_predictions)[m].sum(axis=0) for m in mask])
            h_derivative_truth      = np.array([ (np.transpose(np.array([training_weights[der] for der in bit.derivatives])))[m].sum(axis=0) for m in mask])

            h_ratio_prediction[feature] = h_derivative_prediction/h_w0[feature].reshape(-1,1) 
            h_ratio_truth[feature]      = h_derivative_truth/h_w0[feature].reshape(-1,1)

        n_pads = len(feature_names)+1
        n_col  = int(1.5*sqrt(n_pads))
        n_rows = n_pads/n_col
        if n_rows*n_col<n_pads: n_rows+=1

        for logY in [False, True]:
            c1 = ROOT.TCanvas("c1","multipads",500*n_col,500*n_rows);
            c1.Divide(n_col,n_rows)

            l = ROOT.TLegend(0.2,0.1,0.9,0.85)
            stuff.append(l)
            l.SetNColumns(int(len(config.combinations)/12))
            l.SetFillStyle(0)
            l.SetShadowColor(ROOT.kWhite)
            l.SetBorderSize(0)

            for i_feature, feature in enumerate(feature_names):

                th1d_yield       = helpers.make_TH1F( (h_w0[feature], lin_binning[feature]) )
                c1.cd(i_feature+1)
                ROOT.gStyle.SetOptStat(0)
                th1d_ratio_pred  = { der: helpers.make_TH1F( (h_ratio_prediction[feature][:,i_der], lin_binning[feature])) for i_der, der in enumerate( bit.derivatives ) }
                th1d_ratio_truth = { der: helpers.make_TH1F( (h_ratio_truth[feature][:,i_der], lin_binning[feature])) for i_der, der in enumerate( bit.derivatives ) }
                stuff.append(th1d_yield)
                stuff.append(th1d_ratio_truth)
                stuff.append(th1d_ratio_pred)

                th1d_yield.SetLineColor(ROOT.kGray+2)
                th1d_yield.SetMarkerColor(ROOT.kGray+2)
                th1d_yield.SetMarkerStyle(0)
                th1d_yield.GetXaxis().SetTitle(config.plot_options[feature]['tex'])
                th1d_yield.SetTitle("")

                th1d_yield.Draw("hist")

                for i_der, der in enumerate(bit.derivatives):
                    th1d_ratio_truth[der].SetTitle("")
                    th1d_ratio_truth[der].SetLineColor(color[der])
                    th1d_ratio_truth[der].SetMarkerColor(color[der])
                    th1d_ratio_truth[der].SetMarkerStyle(0)
                    th1d_ratio_truth[der].SetLineWidth(2)
                    th1d_ratio_truth[der].SetLineStyle(ROOT.kDashed)
                    th1d_ratio_truth[der].GetXaxis().SetTitle(config.plot_options[feature]['tex'])

                    th1d_ratio_pred[der].SetTitle("")
                    th1d_ratio_pred[der].SetLineColor(color[der])
                    th1d_ratio_pred[der].SetMarkerColor(color[der])
                    th1d_ratio_pred[der].SetMarkerStyle(0)
                    th1d_ratio_pred[der].SetLineWidth(2)
                    th1d_ratio_pred[der].GetXaxis().SetTitle(config.plot_options[feature]['tex'])

                    tex_name = "_{%s}"%(",".join([config.tex[c].lstrip("C_{")[:-1] if config.tex[c].startswith('C_') else config.tex[c] for c in der]))

                    if i_feature==0:
                        l.AddEntry( th1d_ratio_truth[der], "R"+tex_name)
                        l.AddEntry( th1d_ratio_pred[der],  "#hat{R}"+tex_name)

                if i_feature==0:
                    l.AddEntry( th1d_yield, "yield (SM)")

                max_ = max( map( lambda h:h.GetMaximum(), th1d_ratio_truth.values() ))
                min_ = min( map( lambda h:h.GetMinimum(), th1d_ratio_truth.values() ))
                th1d_yield.Scale(max_/th1d_yield.GetMaximum())
                th1d_yield   .Draw("hist")
                ROOT.gPad.SetLogy(logY)
                th1d_yield   .GetYaxis().SetRangeUser(0.1 if logY else (1.5*min_ if min_<0 else 0.75*min_), 10**(1.5)*max_ if logY else 1.5*max_)
                th1d_yield   .Draw("hist")
                for h in th1d_ratio_truth.values()+th1d_ratio_pred.values():
                    h .Draw("hsame")

            c1.cd(len(feature_names)+1)
            l.Draw()

            lines = [ (0.29, 0.9, 'N_{B} =%5i'%( max_n_tree )) ]
            drawObjects = [ tex.DrawLatex(*line) for line in lines ]
            for o in drawObjects:
                o.Draw()

            plot_directory_ = os.path.join( plot_directory, "training_plots", bit_name, "log" if logY else "lin" )
            if not os.path.isdir(plot_directory_):
                try:
                    os.makedirs( plot_directory_ )
                except IOError:
                    pass
            from RootTools.plot.helpers import copyIndexPHP
            copyIndexPHP( plot_directory_ )
            c1.Print( os.path.join( plot_directory_, "epoch_%05i.png"%(max_n_tree) ) )
            syncer.makeRemoteGif(plot_directory_, pattern="epoch_*.png", name="epoch" )
        syncer.sync()
