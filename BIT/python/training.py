#!/usr/bin/env python

import ROOT, os, sys
from RootTools.core.standard import *
import Analysis.Tools.syncer as syncer

import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--config_module',      action='store', type=str, default = "TMB.BIT.configs", help = "config directory")
argParser.add_argument('--config',             action='store', type=str, default = "ttZ_3l_flavor", help="config")
argParser.add_argument('--name',               action='store', type=str,   default='default', help="Name of the training")
argParser.add_argument('--variable_set',       action='store', type=str,   default='mva_variables', help="List of variables for training")
argParser.add_argument('--output_directory',   action='store', type=str,   default=os.path.expandvars('/groups/hephy/cms/$USER/BIT/'))
argParser.add_argument('--input_directory',    action='store', type=str,   default=os.path.expandvars("/scratch-cbe/users/$USER/BIT/training-ntuples-TTZ-flavor/MVA-training"))
argParser.add_argument('--small',              action='store_true', help="small?")
argParser.add_argument('--debug',              action='store_true', help="Make debug plots?")
#argParser.add_argument('--calibrated',         action='store_true', help="Calibrate output?")
argParser.add_argument('--maxEvents',          action='store', default = None, type=int, help="Maximum number of training events")
argParser.add_argument('--bagging_fraction',   action='store', default = 1., type=float, help="Bagging fraction")
argParser.add_argument('--overwrite',          action='store_true', help="Overwrite output?")
argParser.add_argument('--derivative',         action='store', nargs = '*', default = [], help="What to train?")
argParser.add_argument('--lumi_norm',          action='store_true', help="Normalize the events according to lumiweight1fb?")

#args = argParser.parse_args()
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
        
#extra_args = {key.lstrip('-'):parse_value(value) for key, value in zip(extra[::2], extra[1::2])}

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
for derivative in config.bit_derivatives:
    config.bit_cfg[derivative].update( extra_args )

import uproot
import awkward
import numpy as np
import pandas as pd
#import h5py

#########################################################################################
# variable definitions
if args.small:
    args.name+='_small'
import Analysis.Tools.user as user
# directories
plot_directory   = os.path.join( user. plot_directory, 'MVA', args.config, args.name)
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
mva_variables = [ mva_variable[0] for mva_variable in getattr(config, args.variable_set) ]

n_var_flat   = len(mva_variables)

features = {}
weight_derivatives = {}
lumi_weights       = {}
for i_training_sample, training_sample in enumerate(config.training_samples):
    upfile_name = os.path.join(os.path.expandvars(args.input_directory), training_sample.name, training_sample.name+'.root')
    logger.info( "Loading upfile %i: %s from %s", i_training_sample, training_sample.name, upfile_name)
    upfile = uproot.open(upfile_name)
    features[training_sample.name]            = upfile["Events"].pandas.df(branches = mva_variables )
    weight_derivatives[training_sample.name]  = upfile["Events"].pandas.df(branches = ["weight_derivatives"] )
    if args.lumi_norm:
        lumi_weights[training_sample.name]    = upfile["Events"].pandas.df(branches = ["lumiweight1fb"] )    

features     = pd.concat([features[training_sample.name] for training_sample in config.training_samples])

# Quick fix in training
len_features = len(features)
if args.config == "WH_nlo_delphes":
    map_ = features.mva_WH_deltaR<5.5 
else:
    map_ = np.ones(len(features)).astype('bool')

features = features[map_].values

weight_derivatives = pd.concat([weight_derivatives[training_sample.name] for training_sample in config.training_samples])
if args.lumi_norm:
    lumi_weights = pd.concat([lumi_weights[training_sample.name] for training_sample in config.training_samples])

weight_derivatives = weight_derivatives.values.reshape((len_features,-1))[map_][:args.maxEvents]
features           = features[:args.maxEvents]
if args.lumi_norm:
    lumi_weights       = lumi_weights.values[:,0][map_][:args.maxEvents]
    lumi_weights       = lumi_weights/np.mean(lumi_weights)

# split dataset into Input and output data

# number of samples with 'small'
n_small_samples = 10000

# small
if args.small:
    features = features[:n_small_samples]
    weight_derivatives = weight_derivatives[:n_small_samples]
    if args.lumi_norm:
        lumi_weights       = lumi_weights[:n_small_samples]

features  = features[:,0:n_var_flat]

if args.lumi_norm:
    weight_derivatives = weight_derivatives * lumi_weights[:,np.newaxis]

# Clip the most extreme scores 
for derivative in config.bit_derivatives:

    if config.bit_cfg[derivative]['clip_score_quantile'] is None: continue
    i_derivative          = config.weight_derivative_combinations.index(derivative)
    clip_score_quantile = config.bit_cfg[derivative]['clip_score_quantile']
    
    clip_both_sides = True #(len(derivative)%2==0)

    derivative = config.weight_derivative_combinations[i_derivative]
    training_scores          = np.divide( weight_derivatives[:,i_derivative], weight_derivatives[:,0], np.zeros_like(weight_derivatives[:,i_derivative]), where=weight_derivatives[:,0]!=0 )
    training_score_quantiles = np.quantile( training_scores, [clip_score_quantile, 1-clip_score_quantile] if clip_both_sides else [ 1-clip_score_quantile] )

    ## Always clip pos & neg because we're dealing with nlo
    #if len(derivative)%2==0:
    if False: #clip_both_sides:
        # quadratic: Clip only positive values
        print "Clip score at percentile %3.2f for %r: s>%3.2f" % (1-clip_score_quantile, derivative, training_score_quantiles[0] )
        to_clip                         = training_scores>training_score_quantiles[0]
        weight_derivatives[:,i_derivative][to_clip] = training_score_quantiles[0]*weight_derivatives[:,0][to_clip]

    else:
        # linear: Clip positive and negative values
        print "Clip score at percentiles %3.2f and %3.2f for %r: s<%3.2f and s>%3.2f" % (clip_score_quantile, 1-clip_score_quantile, derivative, training_score_quantiles[0], training_score_quantiles[1] )
        to_clip                  = training_scores>training_score_quantiles[1]
        weight_derivatives[:,i_derivative][to_clip] = training_score_quantiles[1]*weight_derivatives[:,0][to_clip]
        to_clip                  = training_scores<training_score_quantiles[0]
        weight_derivatives[:,i_derivative][to_clip] = training_score_quantiles[0]*weight_derivatives[:,0][to_clip]

## split data into train and test, test_size = 0.2 is quite standard for this
from sklearn.model_selection import train_test_split
train_test_options = {'test_size':0.5, 'random_state':7, 'shuffle':True}
training_features, test_features, training_weights, test_weights = train_test_split(features, weight_derivatives, **train_test_options)

# Boosting
import sys, os, time
sys.path.insert(0,os.path.expandvars("$CMSSW_BASE/src/BIT"))
from BoostedInformationTree import BoostedInformationTree 

# Plot a 1D histogram
def plot1DHist( plot, plot_directory, yRange=(0.3,"auto"), ratio={'yRange':(0.1,1.9)}, legend=(0.2,0.7,0.9,0.9), lumi=137, plotLog=True, histModifications=[], titleOffset=0 ):

    for log in [True, False] if plotLog else [False]:

        if yRange[0] == 0 and log:
            yRange = list(yRange)
            yRange[0] = 0.0003
            yRange = tuple(yRange)

        # Add subdirectory for lin/log plots
        plot_directory_ = os.path.join( plot_directory, "log" if log else "lin" )

        plotting.draw( plot,
                       plot_directory = plot_directory_,
                       logX = False, logY = log, sorting = False,
                       yRange = yRange,
                       ratio = ratio,
#                       drawObjects = drawObjects( lumi, offset=titleOffset ),
                       legend = legend,
                       histModifications = histModifications,
                       copyIndexPHP = True,
                       )

bits = {}
for derivative in config.bit_derivatives:
    if derivative == tuple(): continue
    if args.derivative!=[] and derivative !=tuple(args.derivative): continue
    filename = os.path.join(output_directory, "bit_derivative_%s"% ('_'.join(derivative))) + '.pkl'

    try:
        print ("Loading %s for %r"%( filename, derivative))
        bits[derivative] = BoostedInformationTree.load(filename)
    except IOError:
        args.overwrite = True

    if args.overwrite:

        n_trees = config.bit_cfg[derivative]['n_trees']

        time1 = time.time()
        print ("Learning %s"%( str(derivative)))
        bits[derivative]= BoostedInformationTree(
                training_features     = training_features,
                training_weights      = training_weights[:,0],
                training_diff_weights = training_weights[:, config.weight_derivative_combinations.index(derivative)],
                split_method          = 'vectorized_split_and_weight_sums',
                weights_update_method = 'vectorized',
                bagging_fraction      = args.bagging_fraction,
                **config.bit_cfg[derivative]
                    )

        bits[derivative].boost(debug=args.debug)
        bits[derivative].save(filename)
        print ("Written %s"%( filename ))

        time2 = time.time()
        boosting_time = time2 - time1
        print ("Boosting time: %.2f seconds" % boosting_time)

        # plot loss
        test_scores     = bits[derivative].vectorized_predict(test_features)
        #training_scores = bits[derivative].vectorized_predict(test_features)
        max_score = max(test_scores)
        min_score = min(test_scores)

        test_FIs            = np.zeros(n_trees)
        training_FIs        = np.zeros(n_trees)

        test_weights_     = test_weights    [:, config.weight_derivative_combinations.index(derivative)]
        training_weights_ = training_weights[:, config.weight_derivative_combinations.index(derivative)]

        if args.debug:
            from debug import make_debug_plots
            make_debug_plots( bits[derivative], 
                              training_features, training_weights[:,0],  
                              training_weights[:, config.weight_derivative_combinations.index(derivative)], 
                              test_features, 
                              test_weights[:,0], 
                              test_weights_, 
                              os.path.join(plot_directory, ('_'.join(derivative))),
                              mva_variables = [ config.plot_options[mva_variable[0]]['tex'] for mva_variable in  config.mva_variables ],
                            )

        test_train_factor = (1.-train_test_options['test_size'])/train_test_options['test_size']
        
        for i in range(len(test_features)): 
            test_scores = bits[derivative].predict( test_features[i],     summed = False)
            test_score  = sum( test_scores )
            test_FIs   += test_weights_[i]*test_scores*test_train_factor
        for i in range(len(training_features)): 
            training_scores = bits[derivative].predict( training_features[i],     summed = False)
            training_score  = sum( training_scores )
            training_FIs   += training_weights_[i]*training_scores

        training_FI_histo     = ROOT.TH1D("trainFI", "trainFI",          n_trees, 1, n_trees+1 )
        test_FI_histo         = ROOT.TH1D("testFI",  "testFI",           n_trees, 1, n_trees+1 )

        for i_tree in range(n_trees):
            test_FI_histo    .SetBinContent( i_tree+1, -sum(test_FIs[:i_tree]) )
            training_FI_histo.SetBinContent( i_tree+1, -sum(training_FIs[:i_tree]) )

        # Histo style
        test_FI_histo    .style = styles.lineStyle( ROOT.kBlue, width=2 )
        training_FI_histo.style = styles.lineStyle( ROOT.kRed, width=2 )
        test_FI_histo    .legendText = "Test"
        training_FI_histo.legendText = "Training"

        minY   = 0.01 * min( test_FI_histo.GetBinContent(test_FI_histo.GetMaximumBin()), training_FI_histo.GetBinContent(training_FI_histo.GetMaximumBin()))
        maxY   = 1.5  * max( test_FI_histo.GetBinContent(test_FI_histo.GetMaximumBin()), training_FI_histo.GetBinContent(training_FI_histo.GetMaximumBin()))

        histos = [ [test_FI_histo], [training_FI_histo] ]
        plot   = Plot.fromHisto( "bit_derivative_%s_evolution"% ('_'.join(derivative)), histos, texX="b", texY="L(D,b)" )

        # Plot Style
        histModifications      = []
        histModifications      += [ lambda h: h.GetYaxis().SetTitleOffset(1.4) ]
        histModifications += [ lambda h: h.GetXaxis().SetTitleSize(26) ]
        histModifications += [ lambda h: h.GetYaxis().SetTitleSize(26) ]
        histModifications += [ lambda h: h.GetXaxis().SetLabelSize(22)  ]
        histModifications += [ lambda h: h.GetYaxis().SetLabelSize(22)  ]

        ratioHistModifications = []
        ratio                  = None
        legend                 = (0.6,0.75,0.9,0.88)
        yRange                 = "auto" #( minY, maxY )
        plot1DHist( plot, plot_directory, yRange=yRange, ratio=ratio, legend=legend, histModifications=histModifications )

# make score plot and score-binned feature plots
for derivative in config.bit_derivatives:
    if derivative == tuple(): continue
    if args.derivative!=[] and derivative !=tuple(args.derivative): continue

    # plot loss
    test_scores     = bits[derivative].vectorized_predict(test_features)
    #training_scores = bits[derivative].vectorized_predict(test_features)

    n_digi    = 10
    quantiles = np.quantile( test_scores, np.arange(0,1.1,.1))
    try:
        digi      = np.digitize( test_scores, quantiles)
    except ValueError:
        continue

    #h_score   = ROOT.TH1F("score","score", len(quantiles)-1, quantiles)
    h_feature = { var[0]: {i_bin:ROOT.TH1F("%s_%i"%(var[0], i_bin), "%s_%i"%(var[0], i_bin), 50, min(test_features[:,i_var]), max(test_features[:,i_var])) for i_bin in range(n_digi)} for i_var, var in enumerate(config.mva_variables) }

    test_weights_ = test_weights[:,0]
    for i_event in range(len(test_scores)):
        if i_event%10000==0: print "At",i_event
        #h_score.Fill(test_scores[i_event], test_weights[i_event][0])
        for i_var in range(len(config.mva_variables)):
            h_feature[config.mva_variables[i_var][0]][min(n_digi,digi[i_event])-1].Fill(test_features[i_event][i_var], test_weights_[i_event])

    # Plot Style
    histModifications      = []
    histModifications      += [ lambda h: h.GetYaxis().SetTitleOffset(1.4) ]
    histModifications += [ lambda h: h.GetXaxis().SetTitleSize(26) ]
    histModifications += [ lambda h: h.GetYaxis().SetTitleSize(26) ]
    histModifications += [ lambda h: h.GetXaxis().SetLabelSize(22)  ]
    histModifications += [ lambda h: h.GetYaxis().SetLabelSize(22)  ]

    #plot = Plot.fromHisto( 
    #    "score", 
    #    [[h_score]],
    #    texX = "d_"+"_".join(derivative),
    #    )

    #ratio                  = None
    #legend                 = None #[ (0.15,0.75,0.9,0.88), 3]
    #yRange                 = "auto" #( minY, maxY )
    #plot1DHist( plot, os.path.join(plot_directory, ('_'.join(derivative))), yRange=yRange, ratio=ratio, legend=legend, histModifications=histModifications )

    #assert False, ""
    for i_var in range(len(config.mva_variables)):
        for i_digi in range(n_digi):
            h_feature[config.mva_variables[i_var][0]][i_digi].style     = styles.fillStyle(ROOT.kMagenta-10+i_digi)
            h_feature[config.mva_variables[i_var][0]][i_digi].legendText= "%4.3f #leq d(%s) < %4.3f" % ( quantiles[i_digi], "_".join(derivative), quantiles[i_digi+1] )
        plot = Plot.fromHisto( 
            config.mva_variables[i_var][0]+"_binned", 
            [ [h_feature[config.mva_variables[i_var][0]][i_digi] for i_digi in range(n_digi)]],
            texX = config.mva_variables[i_var][0],
            )

        ratio                  = None
        legend                 = [ (0.15,0.75,0.9,0.88), 2]
        yRange                 = "auto" #( minY, maxY )
        plot1DHist( plot, os.path.join(plot_directory, ('_'.join(derivative))), yRange=yRange, ratio=ratio, legend=legend, histModifications=histModifications )

