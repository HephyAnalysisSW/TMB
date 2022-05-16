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
argParser.add_argument("--plot_directory",     action="store",      default="BIT_VH_14",                 help="plot sub-directory")
argParser.add_argument("--model",              action="store",      default="ZH_Nakamura",                 help="Name of the model")
argParser.add_argument("--nEvents",            action="store",      type=int, default=200000,             help="Number of events")
argParser.add_argument("--nToys",              action="store",      type=int, default=100,             help="Number of toys")
argParser.add_argument('--lumi_factor',        action='store',      type=float, default=1.0, help="Lumi factor" )

args = argParser.parse_args()

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

def make_logR_to_SM( order, truth=False, predictions=None, **kwargs ):

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

event_indices = np.arange(nEvents)
def make_toys( yield_per_toy, n_toys, lin=False, **kwargs):
    weights_      = make_weights(lin=lin, **kwargs) 
    biased_sample = np.random.choice( event_indices, size=50*nEvents,  p = weights_/np.sum(weights_) )

    return np.array( [ np.random.choice( biased_sample, size=n_observed ) for n_observed in np.random.poisson(yield_per_toy, n_toys) ])

tex = ROOT.TLatex()
tex.SetNDC()
tex.SetTextSize(0.04)

colors   = [ ROOT.kRed, ROOT.kBlue, ROOT.kBlack, ROOT.kBlue, ROOT.kRed]
#extended = True
n_toys   = 50000

# do not make the following inconsistent
level          = 0.95

# precompute BITS
h_power     = {}

cfgs = [
        ( "bit_%s_MD3_nTraining_2000000"%args.model, "D = 2", ROOT.kBlue),
        ( "bit_%s_nTraining_2000000"    %args.model ,"D = 3", ROOT.kRed),
        ( "bit_%s_MD5_nTraining_2000000"%args.model, "D = 4", ROOT.kGreen+2),
        ( "bit_%s_MD6_nTraining_2000000"%args.model, "D = 5", ROOT.kBlack),
    ]

#eft = {'cHWtil':0.35, 'cHW':-0.4}
#        {'cHWtil':0.1,  'cHW':0,     'cHQ3':0},
#        {'cHWtil':0.25,  'cHW':0,     'cHQ3':0},
#        {'cHWtil':0.15,  'cHW':0,     'cHQ3':0},
#        {'cHW':0.3,  'cHWtil':0,     'cHQ3':0},
#        {'cHW':0.25,  'cHWtil':0,     'cHQ3':0},
#        {'cHW':0.20,  'cHWtil':0,     'cHQ3':0},
#        {'cHW':0.15,  'cHWtil':0,     'cHQ3':0},

        #{'cHQ3':0.03,  'cHWtil':0,      'cHW':0},
        #{'cHQ3':0.025,  'cHWtil':0,     'cHW':0},
        #{'cHQ3':0.020,  'cHWtil':0,     'cHW':0},
        #{'cHQ3':0.015,  'cHWtil':0,     'cHW':0},
        #{'cHWtil':0.15, 'cHW':-0.15, 'cHQ3':0},
        #{'cHWtil':0.35, 'cHW':-0.4,  'cHQ3':0},
eft = {'cHWtil':0.20,  'cHW':0,     'cHQ3':0}
const     = args.lumi_factor*(lambda_tot(**eft) - lambda_tot())
bits        = {prefix:model.load(prefix=prefix) for prefix,_,_ in cfgs}
predictions = {prefix:{ der:bits[prefix][der].vectorized_predict(features) for der in bits[prefix].keys() } for prefix,_,_ in cfgs} 

def make_plot(i_plot):  
#if True:
#    i_plot=0
    for prefix, tex, color in cfgs:

        predictions_iterations = { der:np.cumsum(bits[prefix][der].vectorized_predict(features, summed = False), axis=0) for der in bits[prefix].keys() } 
        n_iterations   = len(predictions_iterations.values()[0])
        h_power[prefix] = ROOT.TH1F("power", "power", n_iterations, 0 ,n_iterations )
        h_power[prefix].style = styles.lineStyle( color, width = 2)
        h_power[prefix].legendText = tex

        sm_toys = make_toys( args.lumi_factor*lambda_tot(), n_toys ) 
        eft_toys= make_toys( args.lumi_factor*lambda_tot(**eft), n_toys, **eft)

        for iteration in range( -1, n_iterations ):
            if iteration<0:
                event_logR_to_SM = make_logR_to_SM( "total", truth=True, **eft)
            else:
                predictions_i = {key:(value[iteration] if iteration<bits[prefix][key].n_trees else value[bits[prefix][key].n_trees-1]) for key, value in predictions_iterations.iteritems()}
                event_logR_to_SM = make_logR_to_SM( "total", truth=False, predictions=predictions_i, **eft)

            q_null = const-np.array([np.sum( event_logR_to_SM[toy_] ) for toy_ in eft_toys ]) #NULL
            q_alt  = const-np.array([np.sum( event_logR_to_SM[toy_] ) for toy_ in sm_toys ])    #ALT
            # calibration according to the null
            if True:
                n = float(len(q_null))
                mean_q_null     = np.sum(q_null)/n
                sigma_q_null    = sqrt( np.sum((q_null - mean_q_null)**2)/(n-1) )
                q_alt  = (q_alt  - mean_q_null)/sigma_q_null
                q_null = (q_null - mean_q_null)/sigma_q_null

            # Exclusion: The null hypothesis is the BSM point, the alternate is the SM.
            quantile_alt  = np.quantile( q_alt, level )
            quantile_null = np.quantile( q_null, level )
            size_         = np.count_nonzero(q_null>=quantile_null)/float(n_toys)
            #power_histo     = np.histogram( q_theta_given_theta, quantile_SM)

            power_ = np.count_nonzero(q_alt>=quantile_null)/float(n_toys)
            if iteration>=0:
                h_power[prefix].SetBinContent( h_power[prefix].FindBin( iteration ), truth-power_ )
            else:
                truth = power_
            print "iteration",iteration, "size", size_, "power", round(power_,3)
#assert False, ""
    return [h_power[prefix] for prefix, _, _ in cfgs ]

from multiprocessing import Pool
p = Pool(5)
results = p.map(make_plot, range(args.nToys))

h_power = {}
boxes = []
for i_cfg, (prefix, tex, color) in enumerate(cfgs):
    histos = [ results[i][i_cfg] for i in range(len(results)) ]
    h_power[prefix] = {'low':histos[0].Clone(), 'med':histos[0].Clone(), 'high':histos[0].Clone()}
    h_power[prefix]['low'].style = styles.lineStyle( color, width = 1)
    h_power[prefix]['low'].legendText = tex+" (-1#sigma)"
    h_power[prefix]['med'].style = styles.lineStyle( color, width = 2)
    h_power[prefix]['med'].legendText = tex
    h_power[prefix]['high'].style = styles.lineStyle( color, width = 1)
    h_power[prefix]['high'].legendText = tex+" (+1#sigma)"
    
    for i_bin in range( 1, histos[0].GetNbinsX()+1 ):
        vals = np.array([ histos[i_histo].GetBinContent(i_bin) for i_histo in range(len(histos)) ])
        q = np.quantile( vals, np.array([0.16, 0.5, 1-0.16]) )
        h_power[prefix]['low']. SetBinContent( i_bin, q[0] )
        h_power[prefix]['med']. SetBinContent( i_bin, q[1] )
        h_power[prefix]['high'].SetBinContent( i_bin, q[2] )
        boxes.append( ROOT.TBox(h_power[prefix]['low'].GetBinLowEdge(i_bin), q[0], h_power[prefix]['low'].GetBinLowEdge(i_bin) + h_power[prefix]['low'].GetBinWidth(i_bin), q[2]) )
        boxes[-1].SetFillStyle(3004)
        boxes[-1].SetFillColor(color)  

#histos = [[h_power[prefix]['low']] for prefix, _, _ in cfgs ] + [[h_power[prefix]['med']] for prefix, _, _ in cfgs ] + [[h_power[prefix]['high']] for prefix, _, _ in cfgs ]
histos = [[h_power[prefix]['med']] for prefix, _, _ in cfgs ]
plot = Plot.fromHisto(name = "power_evolution_cHW_%3.2f_cHWtil_%3.2f_cHQ3_%3.2f"%( eft['cHW'], eft['cHWtil'], eft['cHQ3']), 
    histos = histos, 
    texX = "Boosting iteration", texY = "#beta-#beta_{opt}" )
plotting.draw(plot, 
    plot_directory = os.path.join( plot_directory, "unbinned"), 
    logY = False, logX = False, copyIndexPHP=True, yRange = (0,.8),
    legend=(0.6,0.7,0.9,0.9),
    drawObjects = boxes,
)
