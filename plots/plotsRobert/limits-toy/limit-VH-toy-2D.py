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
import array

# RootTools
from RootTools.core.standard   import *

# Analysis
import Analysis.Tools.syncer as syncer

#TMB 
import  TMB.Tools.stat_helpers as stat_helpers
ROOT.gROOT.LoadMacro("$CMSSW_BASE/src/TMB/Tools/scripts/tdrstyle.C")
ROOT.setTDRStyle()

# BIT
sys.path.insert(0,os.path.expandvars("$CMSSW_BASE/src/BIT"))
from BoostedInformationTree import BoostedInformationTree

# User
import user

# Parser
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument("--plot_directory",     action="store",      default="BIT_VH_12",                 help="plot sub-directory")
argParser.add_argument("--prefix",             action="store",      default="bit_ZH_Nakamura_nTraining_5000000",                 help="plot sub-directory")
argParser.add_argument("--model",              action="store",      default="ZH_Nakamura",                 help="plot sub-directory")
argParser.add_argument("--WCs",                action="store",      nargs='*', default=["cHQ3", -.08, .08, "cHW", 0,.5],                 help="Wilson coefficients")
argParser.add_argument("--nBins",              action="store",      type=int, default=30,                 help="Number of bins in each dimension")
argParser.add_argument("--nEvents",            action="store",      type=int, default=200000,             help="Number of events")
argParser.add_argument('--truth',              action='store_true', help="Truth?" )
argParser.add_argument('--lumi_factor',        action='store',      type=float, default=1.0, help="Lumi factor" )

args = argParser.parse_args()

WC1, theta1_min, theta1_max, WC2, theta2_min, theta2_max = map( lambda x:float(x[1]) if x[0] in [1,2,4,5] else x[1], enumerate(args.WCs) )

# import the VH model
import VH_models
model = getattr(VH_models, args.model)

feature_names = model.feature_names

# directory for plots
plot_directory = os.path.join( user.plot_directory, args.plot_directory, args.model )

if not os.path.isdir(plot_directory):
    try:
        os.makedirs( plot_directory )
    except IOError:
        pass

# initiate plot
Plot.setDefaults()

features = model.getEvents(args.nEvents)

if args.model.startswith("ZH"):
    adhoc = (features[:,model.feature_names.index('fLL')]>0.2) & (features[:,model.feature_names.index('f2TT')]<3) & (features[:,model.feature_names.index('cos_theta')]>-0.9) & (features[:,model.feature_names.index('cos_theta')]<0.9) & (features[:,model.feature_names.index('f1TT')]>-0.9) & (features[:,model.feature_names.index('f1TT')]<0.9)  & (features[:,model.feature_names.index('f2TT')]<3.5) 
    features = features[adhoc]

nEvents  = len(features)
weights  = model.getWeights(features, eft=model.default_eft_parameters)
print ("Created data set of size %i" % nEvents )

# normalize to SM event yield
lambda_expected_sm      = 90.13 if args.model.startswith("ZH") else 599.87 #Delphes, ptZ>200
lambda_current          = np.sum(weights[tuple()])
for key in weights.keys():
    weights[key] = lambda_expected_sm/lambda_current*weights[key]
    
# predict total yields
sigma_coefficients  = {key:np.sum(weights[key]) for key in model.derivatives}
def lambda_tot( lin=False, **kwargs):
    result =  sigma_coefficients[tuple()]
    result += sum( [ (kwargs[coeff] - model.default_eft_parameters[coeff])*sigma_coefficients[(coeff,)] for coeff in kwargs.keys() ])
    if not lin:
        result += sum( [ .5*(kwargs[coeff1] - model.default_eft_parameters[coeff1])*(kwargs[coeff2] - model.default_eft_parameters[coeff2])*sigma_coefficients[tuple(sorted((coeff1,coeff2)))] for coeff1 in kwargs.keys()  for coeff2 in kwargs.keys()])
    return result 

# xsec ratio
def lambda_ratio( lin=False, **kwargs):
    return lambda_tot( lin=lin, **kwargs ) / lambda_expected_sm

# compute weights for arbitrary WC
def make_weights( lin=False, **kwargs):
    result =  copy.deepcopy(weights[tuple()])
    result += sum( [ (kwargs[coeff] - model.default_eft_parameters[coeff])*weights[(coeff,)] for coeff in kwargs.keys() ])
    if not lin:
        result += sum( [ .5*(kwargs[coeff1] - model.default_eft_parameters[coeff1])*(kwargs[coeff2] - model.default_eft_parameters[coeff2])*weights[tuple(sorted((coeff1,coeff2)))] for coeff1 in kwargs.keys()  for coeff2 in kwargs.keys()])
    return result 

# precompute BITS
bits = model.load(prefix=args.prefix)
predictions = { der:bits[der].vectorized_predict(features) for der in bits.keys() } 

#def make_q( order, eft_norm=None, truth=False, **kwargs ):
#    if eft_norm is not None:
#        # we interpret kwargs as the direction in parameter space, i.e., q = n_theta.t + eft_norm*n_theta^T.s.n_theta
#        eft     = copy.deepcopy(kwargs)
#        norm    = sqrt(sum( np.array(kwargs.values())**2 ) )
#        eft     = { k:eft[k]/norm for k in eft.keys() }
#    else:
#        # we interpret kwargs as the eft parameter point, i.e., q = theta.t + theta^T.s.theta
#        eft_norm = 1
#        eft      = kwargs
#    
#    if order not in ["lin", "quad", "total"]:
#        raise RuntimeError("Order %s not known" % order )
#    result = np.zeros(nEvents)    
#    if order in ["lin", "total"]:
#        for coeff in eft.keys():
#            result += (eft[coeff] - model.default_eft_parameters[coeff])*( weights[(coeff,)]/weights[tuple()] if truth else predictions[(coeff,)]) 
#    if order in ["quad", "total"]:
#        for coeff1 in eft.keys():
#            for coeff2 in eft.keys():
#                prefac = eft_norm if order=="total" else 1
#                result += prefac* .5*(eft[coeff1] - model.default_eft_parameters[coeff1])*(eft[coeff2] - model.default_eft_parameters[coeff2])*( weights[tuple(sorted((coeff1,coeff2)))]/weights[tuple()] if truth else predictions[tuple(sorted((coeff1,coeff2)))]) 
#    return result

def make_q( order, truth=False, predictions = predictions, **kwargs ):
    eft      = model.make_eft(**kwargs)
    if order not in ["lin", "quad", "total"]:
        raise RuntimeError("Order %s not known" % order )
    result = np.zeros(nEvents)
    if order in ["lin", "total"]:
        for coeff in model.wilson_coefficients:
            if eft[coeff] == model.default_eft_parameters[coeff]: continue
            result += (eft[coeff] - model.default_eft_parameters[coeff])*( weights[(coeff,)]/weights[tuple()] if truth else predictions[(coeff,)])
    if order in ["quad", "total"]:
        for coeff1 in model.wilson_coefficients:
            if eft[coeff1] == model.default_eft_parameters[coeff1]: continue
            for coeff2 in model.wilson_coefficients:
                if eft[coeff2] == model.default_eft_parameters[coeff2]: continue
                result += .5*(eft[coeff1] - model.default_eft_parameters[coeff1])*(eft[coeff2] - model.default_eft_parameters[coeff2])*( weights[tuple(sorted((coeff1,coeff2)))]/weights[tuple()] if truth else predictions[tuple(sorted((coeff1,coeff2)))])
    result += 1
    neg_frac = len(result[result<0])/float(len(result))
    if neg_frac>10**-3:
        print "Fraction of negative test statistics for %s: %3.2f"% ( order, neg_frac )
    return 0.5*np.log( result**2 )

    #if order == "lin":
    #    return np.log( (1. + result)**2 )
    #else:
    #    return np.log( (1. + result) )

event_indices = np.arange(nEvents)
def make_toys( yield_per_toy, n_toys, lin=False, **kwargs):
    weights_      = make_weights(lin=lin, **kwargs) 
    biased_sample = np.random.choice( event_indices, size=10*nEvents,  p = weights_/np.sum(weights_) )

    return np.array( [ np.random.choice( biased_sample, size=n_observed ) for n_observed in np.random.poisson(yield_per_toy, n_toys) ])

tex = ROOT.TLatex()
tex.SetNDC()
tex.SetTextSize(0.04)

colors   = [ ROOT.kRed, ROOT.kBlue, ROOT.kBlack, ROOT.kBlue, ROOT.kRed]
#extended = True
n_toys   = 50000

# do not make the following inconsistent
levels          = [ 0.95, 0.68]
quantile_levels = [0.025, 0.16, .5, 1-0.16, 1-0.025]
#quantile_levels = [0.05, 0.32]

def getContours( h, level):
    _h     = h.Clone()
    _h.Smooth(1, "k5b")
    ctmp = ROOT.TCanvas()
    _h.SetContour(1,array.array('d', [level]))
    _h.Draw("contzlist")
    _h.GetZaxis().SetRangeUser(0.0001,1)
    ctmp.Update()
    contours = ROOT.gROOT.GetListOfSpecials().FindObject("contours")
    return contours.At(0).Clone()

##################################
#          Training plot         #
##################################

predictions_iterations = { der:np.cumsum(bits[der].vectorized_predict(features, summed = False), axis=0) for der in bits.keys() } 
n_iterations   = len(predictions_iterations.values()[0])
h_power = ROOT.TH1F("power", "power", n_iterations, 0 ,n_iterations )
h_power.style = styles.lineStyle( ROOT.kBlack, width = 2)
#eft = {'cHWtil':0.35, 'cHW':-0.4}
for eft in [ 
        #{'cHWtil':0.1,  'cHW':0,     'cHQ3':0},
        #{'cHWtil':0.25,  'cHW':0,     'cHQ3':0},
        {'cHWtil':0.20,  'cHW':0,     'cHQ3':0},
        #{'cHWtil':0.15,  'cHW':0,     'cHQ3':0},
        #{'cHW':0.3,  'cHWtil':0,     'cHQ3':0},
        #{'cHW':0.25,  'cHWtil':0,     'cHQ3':0},
        #{'cHW':0.20,  'cHWtil':0,     'cHQ3':0},
        #{'cHW':0.15,  'cHWtil':0,     'cHQ3':0},

        #{'cHQ3':0.03,  'cHWtil':0,      'cHW':0},
        #{'cHQ3':0.025,  'cHWtil':0,     'cHW':0},
        #{'cHQ3':0.020,  'cHWtil':0,     'cHW':0},
        #{'cHQ3':0.015,  'cHWtil':0,     'cHW':0},
        #{'cHWtil':0.15, 'cHW':-0.15, 'cHQ3':0},
        #{'cHWtil':0.35, 'cHW':-0.4,  'cHQ3':0},
        ]:
    sm_toys = make_toys( args.lumi_factor*lambda_tot(), n_toys ) 
    eft_toys= make_toys( args.lumi_factor*lambda_tot(**eft), n_toys, **eft)

    for iteration in range( -1, n_iterations ):
        if iteration<0:
            q_event = make_q( "total", truth=True,  **eft )
        else:
            predictions_i = {key:(value[iteration] if iteration<bits[key].n_trees else value[bits[key].n_trees-1]) for key, value in predictions_iterations.iteritems()}
            q_event = make_q( "total", truth=False, predictions = predictions_i, **eft )

        #log_sigma_tot_ratio_subtraction = np.log(sigma_tot_ratio(theta)) if not extended else 0
        q_theta_given_SM    = np.array([np.sum( q_event[toy_] ) for toy_ in sm_toys ])
        q_theta_given_theta = np.array([np.sum( q_event[toy_] ) for toy_ in eft_toys ])


        assert False, ""
        # calibration according to SM
        if True:
            n = float(len(q_theta_given_SM))
            mean_q_theta_given_SM     = np.sum(q_theta_given_SM)/n
            sigma_q_theta_given_SM    = sqrt( np.sum((q_theta_given_SM - mean_q_theta_given_SM)**2)/(n-1) ) 
            q_theta_given_SM    = (q_theta_given_SM - mean_q_theta_given_SM)/sigma_q_theta_given_SM
            q_theta_given_theta = (q_theta_given_theta - mean_q_theta_given_SM)/sigma_q_theta_given_SM

        # Exclusion: The null hypothesis is the BSM point, the alternate is the SM.
        quantiles_theta = np.quantile( q_theta_given_theta, quantile_levels )
        quantiles_SM    = np.quantile( q_theta_given_SM, quantile_levels )
        size_           = np.sum(np.histogram( q_theta_given_SM, quantiles_SM)[0])/float(n_toys)
        #power_histo     = np.histogram( q_theta_given_theta, quantiles_SM)
        for i_level, level in enumerate(levels[-1:]):
            #if level != 0.68: continue
            power_toy_count = np.count_nonzero((q_theta_given_SM>=quantiles_theta[i_level]) & (q_theta_given_SM<quantiles_theta[-1-i_level]))
            power_ = 1. - power_toy_count/float(n_toys)
            #power[test_statistic][level].SetBinContent( power[test_statistic][level].FindBin( theta1, theta2 ), power_ )
            #power_ = 1-np.sum(np.histogram( q_theta_given_theta, quantiles_SM)[0])/float(n_toys)
            if iteration>=0:
                h_power.SetBinContent( h_power.FindBin( iteration ), power_ )
            else:
                truth = power_
            print "iteration",iteration, "size", quantile_levels[-1-i_level] - quantile_levels[i_level], "power", round(power_,3), WC1, WC2

    plot = Plot.fromHisto(name = "power_evolution_cHW_%3.2f_cHWtil_%3.2f_cHQ3_%3.2f"%( eft['cHW'], eft['cHWtil'], eft['cHQ3']), histos = [[h_power]], texX = "Boosting iteration", texY = "power" )
    line = ROOT.TLine(0,truth,n_iterations,truth)
    line.SetLineWidth(2)
    plotting.draw(plot, plot_directory = os.path.join( plot_directory, "unbinned"), logY = False, logX = False, copyIndexPHP=True, drawObjects = [ROOT.TLine(0,truth,n_iterations,truth)], yRange = (0,1))

syncer.sync()

##################################
#            2D plot             #
##################################

step1 = (theta1_max-theta1_min)/args.nBins
step2 = (theta2_max-theta2_min)/args.nBins
theta1_vals = np.arange(theta1_min, theta1_max+step1, (theta1_max-theta1_min)/args.nBins) 
theta2_vals = np.arange(theta2_min, theta2_max+step2, (theta2_max-theta2_min)/args.nBins) 

test_statistics = ["total"]

power = {}
for test_statistic in test_statistics: 
#for test_statistic in ["total"]: 

    truth_txt = "truth" if args.truth else "predicted"
    print "Test statistic", test_statistic, "truth?", args.truth

    power[test_statistic] = {level:ROOT.TH2D("power_"+test_statistic, "power_"+test_statistic, len(theta1_vals)-1, array.array('d', theta1_vals), len(theta2_vals)-1, array.array('d', theta2_vals)) for level in levels}

    min_, max_ = float('inf'), -float('inf')

    sm_toys = make_toys( args.lumi_factor*lambda_tot(), n_toys ) 
    for i_theta1, theta1 in enumerate( theta1_vals ):
        #for i_theta1, theta1 in enumerate( [0] ):
        for i_theta2, theta2 in enumerate( theta2_vals ):
            #for i_theta2, theta2 in enumerate( [.01, .05, .1, .15, .2] ):

            if theta1==theta2==0: continue

            eft     = {WC1:theta1, WC2:theta2}
            q_event = make_q( test_statistic, truth=args.truth, **eft )

            #log_sigma_tot_ratio_subtraction = np.log(sigma_tot_ratio(theta)) if not extended else 0
            q_theta_given_SM    = np.array([np.sum( q_event[toy_] ) for toy_ in sm_toys ])
            q_theta_given_theta = np.array([np.sum( q_event[toy_] ) for toy_ in make_toys( args.lumi_factor*lambda_tot(**eft), n_toys, **eft) ])

            # calibration according to SM
            if True:
                n = float(len(q_theta_given_SM))
                mean_q_theta_given_SM     = np.sum(q_theta_given_SM)/n
                sigma_q_theta_given_SM    = sqrt( np.sum((q_theta_given_SM - mean_q_theta_given_SM)**2)/(n-1) ) 
                q_theta_given_SM    = (q_theta_given_SM - mean_q_theta_given_SM)/sigma_q_theta_given_SM
                q_theta_given_theta = (q_theta_given_theta - mean_q_theta_given_SM)/sigma_q_theta_given_SM

            # Exclusion: The null hypothesis is the BSM point, the alternate is the SM.
            quantiles_theta = np.quantile( q_theta_given_theta, quantile_levels )
            quantiles_SM    = np.quantile( q_theta_given_SM, quantile_levels )
            size_           = np.sum(np.histogram( q_theta_given_SM, quantiles_SM)[0])/float(n_toys)
            #power_histo     = np.histogram( q_theta_given_theta, quantiles_SM)
            for i_level, level in enumerate(levels):
                #if level != 0.68: continue
                power_toy_count = np.count_nonzero((q_theta_given_SM>=quantiles_theta[i_level]) & (q_theta_given_SM<quantiles_theta[-1-i_level]))
                power_ = 1. - power_toy_count/float(n_toys)
                power[test_statistic][level].SetBinContent( power[test_statistic][level].FindBin( theta1, theta2 ), power_ )
                #power_ = 1-np.sum(np.histogram( q_theta_given_theta, quantiles_SM)[0])/float(n_toys)
                print "theta", round(theta1,3), round(theta2,3), "size", quantile_levels[-1-i_level] - quantile_levels[i_level], "power", round(power_,3), test_statistic, WC1, WC2,  "truth", args.truth

colors   = { 'quad':ROOT.kRed, 'lin':ROOT.kBlue, 'total':ROOT.kBlack}

contours = { key:{level:getContours( power[key][level], 0.5 ) for level in levels } for key in power.keys()} 
contour_objects = []
for test_statistic in contours.keys():
    for level in contours[test_statistic].keys():
        for i_tgraph in range(contours[test_statistic][level].GetSize()):
            tgr = contours[test_statistic][level].At(i_tgraph)
            print i_tgraph, tgr, test_statistic, level

            tgr.SetLineColor(colors[test_statistic])
            tgr.SetLineWidth(2)
            tgr.SetLineStyle(ROOT.kDashed if level!=0.95 else 1)
            tgr.SetMarkerStyle(0)
            contour_objects.append( tgr )

for test_statistic in test_statistics:
    for level in levels: 
        plot2D = Plot2D.fromHisto(name = "power_%s_%s_vs_%s_%s_lumi_factor_%3.2f_level_%3.2f"%(test_statistic, WC1, WC2, ("truth" if args.truth else "predicted"), args.lumi_factor, level), histos = [[power[test_statistic][level]]], texX = WC1, texY = WC2 )
        plotting.draw2D(plot2D, plot_directory = os.path.join( plot_directory, "unbinned"), histModifications = [lambda h:ROOT.gStyle.SetPalette(58)], logY = False, logX = False, logZ = True, copyIndexPHP=True, drawObjects = contour_objects, zRange = (0.01,1))

syncer.sync()
