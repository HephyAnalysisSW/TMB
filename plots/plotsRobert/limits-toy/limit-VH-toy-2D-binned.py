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
argParser.add_argument("--plot_directory",     action="store",      default="BIT_VH_8",                 help="plot sub-directory")
argParser.add_argument("--model",              action="store",      default="ZH_Nakamura",                 help="plot sub-directory")
argParser.add_argument("--WCs",                action="store",      nargs='*', default=["cHQ3", -.08, .08, "cHW", 0,.5],                 help="Wilson coefficients")
argParser.add_argument("--nBins",              action="store",      type=int, default=30,                 help="Number of bins in each dimension")
argParser.add_argument("--nBinsTestStat",      action="store",      type=int, default=20,                 help="Number of bins for the test statistic")
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

adhoc = (features[:,model.feature_names.index('fLL')]>0.2) & (features[:,model.feature_names.index('f2TT')]<3) & (features[:,model.feature_names.index('cos_theta')]>-0.9) & (features[:,model.feature_names.index('cos_theta')]<0.9) & (features[:,model.feature_names.index('f1TT')]>-0.9) & (features[:,model.feature_names.index('f1TT')]<0.9)  & (features[:,model.feature_names.index('f2TT')]<3.5) 
 
features = features[adhoc]

nEvents  = len(features)
weights  = model.getWeights(features, eft=model.default_eft_parameters)
print ("Created data set of size %i" % nEvents )

# normalize to SM event yield
lambda_expected_sm      = 90.13 #Delphes, ptZ>200
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
bits = model.load(prefix="bit_ZH_Nakamura_MD6_nTraining_2000000")
predictions = { der:bits[der].vectorized_predict(features) for der in bits.keys() } 

def make_q( order, truth=False, **kwargs ):
    eft      = kwargs
    if order not in ["lin", "quad", "total"]:
        raise RuntimeError("Order %s not known" % order )
    result = np.zeros(nEvents)
    if order in ["lin", "total"]:
        for coeff in eft.keys():
            result += (eft[coeff] - model.default_eft_parameters[coeff])*( weights[(coeff,)]/weights[tuple()] if truth else predictions[(coeff,)])
    if order in ["quad", "total"]:
        for coeff1 in eft.keys():
            for coeff2 in eft.keys():
                result += .5*(eft[coeff1] - model.default_eft_parameters[coeff1])*(eft[coeff2] - model.default_eft_parameters[coeff2])*( weights[tuple(sorted((coeff1,coeff2)))]/weights[tuple()] if truth else predictions[tuple(sorted((coeff1,coeff2)))])
    result+=1
    neg_frac = len(result[result<0])/float(len(result))
    if neg_frac>10**-3:
        print "Fraction of negative test statistics for %s: %3.2f"% ( order, neg_frac )
    return 0.5*np.log( result**2 )
#event_indices = np.arange(nEvents)
#def make_toys( yield_per_toy, n_toys, lin=False, **kwargs):
#    weights_      = make_weights(lin=lin, **kwargs)
#    biased_sample = np.random.choice( event_indices, size=10*nEvents,  p = weights_/np.sum(weights_) )
#
#    return np.array( [ np.random.choice( biased_sample, size=n_observed ) for n_observed in np.random.poisson(yield_per_toy, n_toys) ])

tex = ROOT.TLatex()
tex.SetNDC()
tex.SetTextSize(0.04)


#colors   = [ ROOT.kRed, ROOT.kBlue, ROOT.kBlack, ROOT.kBlue, ROOT.kRed]
#levels   = [ .05,       .32,        .5,          .68,        .95      ]
#extended = True
#n_toys   = 50000

def getContours( h, level):
    _h     = h.Clone()
    ctmp = ROOT.TCanvas()
    _h.SetContour(1,array.array('d', [level]))
    _h.Draw("contzlist")
    _h.GetZaxis().SetRangeUser(0.0001,1)
    ctmp.Update()
    contours = ROOT.gROOT.GetListOfSpecials().FindObject("contours")
    return contours.At(0).Clone()

colors     = [ROOT.kBlack, ROOT.kBlue, ROOT.kGreen, ROOT.kMagenta, ROOT.kCyan, ROOT.kRed]
for WC, theta_vals in [
            ("cHW",    [-.5, -.2, -.1, .1, .2, .5]), 
            ("cHWtil", [-.5, -.2, -.1, .1, .2, .5]),
            #("cHQ3",   [-.05, -.02, -.01, .01, .02, .05]),
            ("cHQ3",   [-1.2, -.08, -.02, .01, .02, .05]),
        ]: 

    q_theta_given_theta = {}
    q_theta_given_SM    = {}
    for test_statistic in ["lin", "quad", "total"]: 
    #for test_statistic in ["lin"]: 
        q_theta_given_theta [test_statistic]= {}
        q_theta_given_SM    [test_statistic]= {}

        truth_txt = "truth" if args.truth else "predicted"
        print "Test statistic", test_statistic, "truth?", args.truth

        histos = []
        for i_theta, theta in enumerate(theta_vals):
            print "theta", theta
            q_event = make_q( test_statistic, truth=args.truth, **{WC:theta} )

            q_event_argsort     = np.argsort(q_event)
            q_event_argsort_inv = np.argsort(q_event_argsort)
            cdf_sm = np.cumsum(weights[()][q_event_argsort])
            cdf_sm/=cdf_sm[-1]

            q_event_cdf = cdf_sm[q_event_argsort_inv] #uniformly distributed under the SM hypothesis

            #min_, max_ = 0, 1 
            binning = np.linspace(0, 1, args.nBinsTestStat+1)

            np_histo_SM    = np.histogram(q_event_cdf, bins=binning, weights = args.lumi_factor*weights[()])
            np_histo_theta = np.histogram(q_event_cdf, bins=binning, weights = args.lumi_factor*make_weights(lin=False, **{WC:theta})  )
            #print np_histo_SM
            #print np_histo_theta
            histo_SM       = stat_helpers.make_TH1F(np_histo_SM)
            histo_theta    = stat_helpers.make_TH1F(np_histo_theta)

            histo_SM.legendText    = "#color[%i]{p(q_{%s=%3.2f}|SM)}" % ( colors[i_theta], WC, theta)
            histo_theta.legendText = "#color[%i]{p(q_{%s=%3.2f}|%s=%3.2f)}" % ( colors[i_theta], WC, theta, WC, theta )
            histo_SM.style         = styles.lineStyle( colors[i_theta], dashed = True)
            histo_theta.style      = styles.lineStyle( colors[i_theta] ) 

            histos.append( histo_SM )
            histos.append( histo_theta )

        ## Text on the plots
        #lines = [ 
        #        #  (0.25, 0.88, "#color[4]{%i%% qu. q_{BSM} = %3.2f}" % ( 100*(1-CL), q_theta_given_theta_1mCL[theta]) ),
        #        #  (0.25, 0.83, "#color[2]{q_{SM} = %3.2f}" % ( q_theta_SM ) ),
        #        #  (0.25, 0.78, "#theta_{current} = %5.4f" % theta ),
        #        ]
        #drawObjects = [ tex.DrawLatex(*line) for line in lines ]
        drawObjects = [ ]

        plot = Plot.fromHisto( "test_stat_%s_%s_%s_lumi_factor_%3.2f_nBinsTestStat_%i"%(test_statistic, truth_txt, WC, args.lumi_factor, args.nBinsTestStat), [[h] for h in histos], texX = "q (%s)"%truth_txt, texY = "Entries" )
        plotting.draw( plot,
            plot_directory = os.path.join( plot_directory, "binned" ),
            #ratio          = {'yRange':(0.6,1.4)} if len(plot.stack)>=2 else None,
            logX = False, sorting = False,
            legend         = ( (0.15,0.7,0.9,0.92),2),
            drawObjects    = drawObjects,
            copyIndexPHP   = True,
            extensions     = ["png"], 
          )  

n_toys = 50000

# do not make the following inconsistent
levels          = [ 0.95, 0.68]
quantile_levels = [ 0.025, 0.16, .5, 1-0.16, 1-0.025 ]

exp_nll_ratio = {}

step1 = (theta1_max-theta1_min)/args.nBins
step2 = (theta2_max-theta2_min)/args.nBins
theta1_vals = np.arange(theta1_min, theta1_max+step1, (theta1_max-theta1_min)/args.nBins) 
theta2_vals = np.arange(theta2_min, theta2_max+step2, (theta2_max-theta2_min)/args.nBins) 

test_statistics = ["total", "lin", "quad"]
#test_statistics = ["total"]
#test_statistics = ["quad"]
exp_nll_ratio = {}
power         = {}
for test_statistic in test_statistics: 
#for test_statistic in ["total"]: 
    truth_txt = "truth" if args.truth else "predicted"
    print "Test statistic", test_statistic, "truth?", args.truth
    power[test_statistic] = {level:ROOT.TH2D("power_"+test_statistic, "power_"+test_statistic, len(theta1_vals)-1, array.array('d', theta1_vals), len(theta2_vals)-1, array.array('d', theta2_vals)) for level in levels}

    exp_nll_ratio[test_statistic] = ROOT.TH2D("exp_nll_ratio_"+test_statistic, "exp_nll_ratio_"+test_statistic, len(theta1_vals)-1, array.array('d', theta1_vals), len(theta2_vals)-1, array.array('d', theta2_vals))

    for i_theta1, theta1 in enumerate( theta1_vals ):
    #for i_theta1, theta1 in enumerate( [0] ):
        for i_theta2, theta2 in enumerate( theta2_vals ):
        #for i_theta2, theta2 in enumerate( [.01, .05, .1, .15, .2] ):

            if theta1==theta2==0: continue

            eft     = {WC1:theta1, WC2:theta2}
            q_event = make_q( test_statistic, truth=args.truth, **eft )


            q_event_argsort     = np.argsort(q_event)
            q_event_argsort_inv = np.argsort(q_event_argsort)
            cdf_sm = np.cumsum(weights[()][q_event_argsort])
            cdf_sm/=cdf_sm[-1]

            q_event_cdf = cdf_sm[q_event_argsort_inv] #uniformly distributed under the SM hypothesis

            #min_, max_ = min( q_event_cdf ), max( q_event_cdf )
            binning    = np.linspace(0, 1, args.nBinsTestStat+1) 

            np_histo_SM    = np.histogram(q_event_cdf, bins=binning, weights = args.lumi_factor*weights[()])[0]
            np_histo_theta = np.histogram(q_event_cdf, bins=binning, weights = args.lumi_factor*make_weights(lin=False, **eft)  )[0]

            # Expectation_BSM( -Log( Prod_i( Pois_i( n_i, lambda_i(theta))/Pois_i( n_i, lambda_i(0)) ) ))
            exp_nll_ratio_ = 2*np.sum(np_histo_SM - np_histo_theta - np_histo_theta*np.log(np_histo_SM/np_histo_theta))
            exp_nll_ratio[test_statistic].SetBinContent( exp_nll_ratio[test_statistic].FindBin( theta1, theta2 ), exp_nll_ratio_)

            binned_toys_SM    = np.random.poisson(lam=np_histo_SM, size=(n_toys, len(np_histo_SM)))
            binned_toys_theta = np.random.poisson(lam=np_histo_theta, size=(n_toys, len(np_histo_theta)))

            q_theta_given_SM    = [ np.sum( toy_ll ) for toy_ll in 2*(np_histo_SM - np_histo_theta  - binned_toys_SM*np.log(np_histo_SM/np_histo_theta))]
            q_theta_given_theta = [ np.sum( toy_ll ) for toy_ll in 2*(np_histo_SM - np_histo_theta  - binned_toys_theta*np.log(np_histo_SM/np_histo_theta))]

            #toys_theta = make_toys( args.lumi_factor*lambda_tot(**eft), n_toys, **eft)
            #toys_sm    = make_toys( args.lumi_factor*lambda_tot(), n_toys )

            if True:
                n = float(len(q_theta_given_SM))
                mean_q_theta_given_SM     = np.sum(q_theta_given_SM)/n
                sigma_q_theta_given_SM    = sqrt( np.sum((q_theta_given_SM - mean_q_theta_given_SM)**2)/(n-1) )
                q_theta_given_SM    = (q_theta_given_SM - mean_q_theta_given_SM)/sigma_q_theta_given_SM
                q_theta_given_theta = (q_theta_given_theta - mean_q_theta_given_SM)/sigma_q_theta_given_SM

            print i_theta1, theta1, i_theta2, theta2, "sqrt(2NLL)", sqrt(abs(exp_nll_ratio_))

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
                print "theta", round(theta1,3), round(theta2,3), "level", level, "size", quantile_levels[-1-i_level] - quantile_levels[i_level], "power", round(power_,3), test_statistic, WC1, WC2,  "truth", args.truth


colors   = { 'quad':ROOT.kRed, 'lin':ROOT.kBlue, 'total':ROOT.kBlack}
nll_levels = [2.27, 5.99]
contours = { key:{level:getContours( exp_nll_ratio[key], level) for level in nll_levels} for key in exp_nll_ratio.keys() }
contour_objects = []
for test_statistic in contours.keys():
    for level in contours[test_statistic].keys():
        for i_tgraph in range(contours[test_statistic][level].GetSize()):
            tgr = contours[test_statistic][level].At(i_tgraph)
            print i_tgraph, tgr, test_statistic, level

            tgr.SetLineColor(colors[test_statistic])
            tgr.SetLineWidth(2)
            tgr.SetLineStyle(ROOT.kDashed if level!=nll_levels[-1] else 1)
            tgr.SetMarkerStyle(0)
            contour_objects.append( tgr )

for test_statistic in test_statistics:    
    plot2D = Plot2D.fromHisto(name = "exp_nll_ratio_%s_%s_vs_%s_%s_lumi_factor_%3.2f_nBinsTestStat_%i"%(test_statistic, WC1, WC2, ("truth" if args.truth else "predicted"), args.lumi_factor, args.nBinsTestStat), histos = [[exp_nll_ratio[test_statistic]]], texX = WC1, texY = WC2 )
    plotting.draw2D(plot2D, plot_directory = os.path.join( plot_directory, "binned"), histModifications = [lambda h:ROOT.gStyle.SetPalette(58)], logY = False, logX = False, logZ = True, copyIndexPHP=True, drawObjects = contour_objects, zRange = (0.05,25))


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
        plot2D = Plot2D.fromHisto(name = "power_%s_%s_vs_%s_%s_lumi_factor_%3.2f_level_%3.2f_nBinsTestStat_%i"%(test_statistic, WC1, WC2, ("truth" if args.truth else "predicted"), args.lumi_factor, level, args.nBinsTestStat), histos = [[power[test_statistic][level]]], texX = WC1, texY = WC2 )
        plotting.draw2D(plot2D, plot_directory = os.path.join( plot_directory, "binned"), histModifications = [lambda h:ROOT.gStyle.SetPalette(58)], logY = False, logX = False, logZ = True, copyIndexPHP=True, drawObjects = contour_objects, zRange = (0.005,1))

syncer.sync()
