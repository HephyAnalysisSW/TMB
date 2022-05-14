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
if args.model == 'ZH_Nakamura':
    bits = model.load(prefix="bit_ZH_Nakamura_MD6_nTraining_2000000")
elif args.model == 'WH_Nakamura':
    bits = model.load(prefix="bit_WH_Nakamura_MD6_nTraining_2000000")
predictions = { der:bits[der].vectorized_predict(features) for der in bits.keys() } 

def make_logR_to_SM( order, truth=False, **kwargs ):

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

tex = {'cHW':"C_{HW}", 'cHWtil':"C_{H#tilde{W}}", "cHQ3":"C_{HQ}^{(3)}"}

#colors   = [ ROOT.kRed, ROOT.kBlue, ROOT.kBlack, ROOT.kBlue, ROOT.kRed]
#levels   = [ .05,       .32,        .5,          .68,        .95      ]
#extended = True
#n_toys   = 50000

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
            #q_event = make_q( test_statistic, truth=args.truth, **{WC:theta} )

            eft              = {WC:theta}
            event_logR_to_SM = make_logR_to_SM( test_statistic, truth=args.truth, **eft )
            const            = args.lumi_factor*(lambda_tot(**eft) - lambda_tot())

            q_event = const - event_logR_to_SM

            q_event_argsort     = np.argsort(q_event)
            q_event_argsort_inv = np.argsort(q_event_argsort)
            cdf_sm = np.cumsum(weights[()][q_event_argsort])
            cdf_sm/=cdf_sm[-1]

            q_event_cdf = cdf_sm[q_event_argsort_inv] #uniformly distributed under the SM hypothesis

            #min_, max_ = 0, 1 
            binning = np.linspace(0, 1, args.nBinsTestStat+1)

            np_histo_SM    = np.histogram(q_event_cdf, bins=binning, weights = args.lumi_factor*weights[()])
            np_histo_theta = np.histogram(q_event_cdf, bins=binning, weights = args.lumi_factor*make_weights(lin=False, **eft)  )
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

levels          = [ 0.68, 0.95]
sigmas          = { 0.68:1, 0.95:2}         
n_sigma_gauss = {}

step1 = (theta1_max-theta1_min)/args.nBins
step2 = (theta2_max-theta2_min)/args.nBins
theta1_vals = np.arange(theta1_min, theta1_max+step1, (theta1_max-theta1_min)/args.nBins) 
theta2_vals = np.arange(theta2_min, theta2_max+step2, (theta2_max-theta2_min)/args.nBins) 

test_statistics = ["total", "lin", "quad"]
#test_statistics = ["total"]
#test_statistics = ["quad"]
n_sigma_gauss = {}
power         = {}
for test_statistic in test_statistics: 
#for test_statistic in ["total"]: 
    truth_txt = "truth" if args.truth else "predicted"
    print "Test statistic", test_statistic, "truth?", args.truth
    power[test_statistic] = {level:ROOT.TH2D("power_"+test_statistic, "power_"+test_statistic, len(theta1_vals)-1, array.array('d', theta1_vals), len(theta2_vals)-1, array.array('d', theta2_vals)) for level in levels}

    n_sigma_gauss[test_statistic] = ROOT.TH2D("n_sigma_gauss_"+test_statistic, "n_sigma_gauss_"+test_statistic, len(theta1_vals)-1, array.array('d', theta1_vals), len(theta2_vals)-1, array.array('d', theta2_vals))

    for i_theta1, theta1 in enumerate( theta1_vals ):
    #for i_theta1, theta1 in enumerate( [0] ):
        for i_theta2, theta2 in enumerate( theta2_vals ):
        #for i_theta2, theta2 in enumerate( [.01, .05, .1, .15, .2] ):

            if theta1==theta2==0: continue

            eft     = {WC1:theta1, WC2:theta2}
            #q_event = make_q( test_statistic, truth=args.truth, **eft )

            event_logR_to_SM = make_logR_to_SM( test_statistic, truth=args.truth, **eft )
            const            = args.lumi_factor*(lambda_tot(**eft) - lambda_tot())

            q_event = const - event_logR_to_SM

            q_event_argsort     = np.argsort(q_event)
            q_event_argsort_inv = np.argsort(q_event_argsort)
            cdf_sm = np.cumsum(weights[()][q_event_argsort])
            cdf_sm/=cdf_sm[-1]

            q_event_cdf = cdf_sm[q_event_argsort_inv] #uniformly distributed under the SM hypothesis

            #min_, max_ = min( q_event_cdf ), max( q_event_cdf )
            binning    = np.linspace(0, 1, args.nBinsTestStat+1) 

            np_histo_SM    = np.histogram(q_event_cdf, bins=binning, weights = args.lumi_factor*weights[()])[0]
            np_histo_theta = np.histogram(q_event_cdf, bins=binning, weights = args.lumi_factor*make_weights(lin=False, **eft)  )[0]

            n_sigma_gauss_ = np.sum( (np_histo_theta - np_histo_SM)*np.log(np_histo_theta/np_histo_SM) )/sqrt( np.sum(np_histo_theta*np.log(np_histo_theta/np_histo_SM)**2) )
            n_sigma_gauss[test_statistic].SetBinContent( n_sigma_gauss[test_statistic].FindBin( theta1, theta2 ), n_sigma_gauss_)
            print i_theta1, theta1, i_theta2, theta2, "n_sigma_gauss", round(n_sigma_gauss_,3)

            binned_toys_SM    = np.random.poisson(lam=np_histo_SM, size=(n_toys, len(np_histo_SM)))
            binned_toys_theta = np.random.poisson(lam=np_histo_theta, size=(n_toys, len(np_histo_theta)))

            q_theta_given_SM    = np.array([ np.sum( toy_ll ) for toy_ll in (np_histo_theta - np_histo_SM  - binned_toys_SM*np.log(np_histo_theta/np_histo_SM))])
            q_theta_given_theta = np.array([ np.sum( toy_ll ) for toy_ll in (np_histo_theta - np_histo_SM  - binned_toys_theta*np.log(np_histo_theta/np_histo_SM))])

            toy_null_argsort = np.argsort( q_theta_given_theta )
            p_value_alt  = np.interp( q_theta_given_SM,    q_theta_given_theta[toy_null_argsort],1-np.insert(np.cumsum(np.ones(len(q_theta_given_theta)-1))/float(len(q_theta_given_theta)-1),0,0) )
            p_value_null = np.interp( q_theta_given_theta, q_theta_given_theta[toy_null_argsort],1-np.insert(np.cumsum(np.ones(len(q_theta_given_theta)-1))/float(len(q_theta_given_theta)-1),0,0) )

            sizes  = {level: np.count_nonzero(p_value_null<=1-level)/float(n_toys) for level in levels}
            powers = {level: np.count_nonzero(p_value_alt <=1-level)/float(n_toys) for level in levels}
            #if True:
            #    n = float(len(q_theta_given_SM))
            #    mean_q_theta_given_SM     = np.sum(q_theta_given_SM)/n
            #    sigma_q_theta_given_SM    = sqrt( np.sum((q_theta_given_SM - mean_q_theta_given_SM)**2)/(n-1) )
            #    q_theta_given_SM    = (q_theta_given_SM - mean_q_theta_given_SM)/sigma_q_theta_given_SM
            #    q_theta_given_theta = (q_theta_given_theta - mean_q_theta_given_SM)/sigma_q_theta_given_SM


            ## Exclusion: The null hypothesis is the BSM point, the alternate is the SM.
            #quantiles_theta = np.quantile( q_theta_given_theta, levels )
            #quantiles_SM    = np.quantile( q_theta_given_SM, levels )
            ##power_histo     = np.histogram( q_theta_given_theta, quantiles_SM)
            #sizes  = {level:np.count_nonzero(q_theta_given_theta<=quantiles_theta[i_level])/float(n_toys) for i_level, level in enumerate(levels)}
            #powers = {level:np.count_nonzero(q_theta_given_SM>quantiles_theta[i_level])/float(n_toys) for i_level, level in enumerate(levels)}

            if True: #debug histo
                drawObjects = []
                binning = np.linspace(np.min(np.concatenate((q_theta_given_theta, q_theta_given_SM))), np.max(np.concatenate((q_theta_given_theta, q_theta_given_SM))),20)
                h_null = stat_helpers.make_TH1F( np.histogram( q_theta_given_theta, binning) )
                h_null.style = styles.lineStyle(ROOT.kBlack)
                h_alt  = stat_helpers.make_TH1F( np.histogram( q_theta_given_SM, binning) )
                h_alt.style = styles.lineStyle(ROOT.kRed)

                min_, max_ = min([h_null.GetMinimum(), h_alt.GetMinimum()]), max([h_null.GetMaximum(), h_alt.GetMaximum()])

                null_mean  = np.mean( q_theta_given_theta )
                null_sigma = sqrt(np.var( q_theta_given_theta ))
                null_median = np.median( q_theta_given_theta )
                null_quantiles = np.quantile( q_theta_given_theta, levels )

                drawObjects.append( ROOT.TLine( null_mean, min_, null_mean, max_ ) )
                drawObjects[-1].SetLineColor(ROOT.kBlack)
                drawObjects.append( ROOT.TLine( null_median, min_, null_median, max_ ) )
                drawObjects[-1].SetLineColor(ROOT.kBlack)
                drawObjects[-1].SetLineStyle(ROOT.kDashed)

                alt_mean  = np.mean( q_theta_given_SM )
                alt_sigma = sqrt(np.var( q_theta_given_SM ))
                alt_median = np.median( q_theta_given_SM )
                alt_quantiles = np.quantile( q_theta_given_SM, levels )
                drawObjects.append( ROOT.TLine( alt_mean, min_, alt_mean, max_ ) )
                drawObjects[-1].SetLineColor(ROOT.kRed)
                drawObjects.append( ROOT.TLine( alt_median, min_, alt_median, max_ ) )
                drawObjects[-1].SetLineColor(ROOT.kRed)
                drawObjects[-1].SetLineStyle(ROOT.kDashed)


                plot = Plot.fromHisto( "debug_%s_%s_%s_lumi_factor_%3.2f_nBinsTestStat_%i"%(test_statistic, truth_txt, WC, args.lumi_factor, args.nBinsTestStat), [[h] for h in [h_null,h_alt]], texX = "q", texY = "Entries" )
                plotting.draw( plot,
                    plot_directory = os.path.join( plot_directory, "binned", "debug" ),
                    logX = False, sorting = False,
                    #legend         = ( (0.15,0.7,0.9,0.92),2),
                    drawObjects    = drawObjects,
                    copyIndexPHP   = True,
                    extensions     = ["png","pdf"], 
                  )  
            
            for i_level, level in enumerate(levels):
                power[test_statistic][level].SetBinContent( power[test_statistic][level].FindBin( theta1, theta2 ), powers[level] )
                print "theta", round(theta1,3), round(theta2,3), "level", level, "size", round(sizes[level],3), "power", round(powers[level],3), test_statistic, WC1, WC2,  "truth", args.truth

            assert False, ""

colors   = { 'quad':ROOT.kRed, 'lin':ROOT.kBlue, 'total':ROOT.kBlack}
contours = { key:{level:getContours( n_sigma_gauss[key], sigmas[level]) for level in levels} for key in n_sigma_gauss.keys() }
contour_objects = []
for test_statistic in contours.keys():
    for level in contours[test_statistic].keys():
        for i_tgraph in range(contours[test_statistic][level].GetSize()):
            tgr = contours[test_statistic][level].At(i_tgraph)
            print i_tgraph, tgr, test_statistic, level

            tgr.SetLineColor(colors[test_statistic])
            tgr.SetLineWidth(2)
            tgr.SetLineStyle(ROOT.kDashed if sigmas[level]==1 else 1)
            tgr.SetMarkerStyle(0)
            contour_objects.append( tgr )

for test_statistic in test_statistics:    
    plot2D = Plot2D.fromHisto(name = "n_sigma_gauss_%s_%s_vs_%s_%s_lumi_factor_%3.2f_nBinsTestStat_%i"%(test_statistic, WC1, WC2, ("truth" if args.truth else "predicted"), args.lumi_factor, args.nBinsTestStat), histos = [[n_sigma_gauss[test_statistic]]], texX = tex[WC1], texY = tex[WC2] )
    plotting.draw2D(plot2D, plot_directory = os.path.join( plot_directory, "binned"), histModifications = [lambda h:ROOT.gStyle.SetPalette(58)], logY = False, logX = False, logZ = True, copyIndexPHP=True, drawObjects = contour_objects, zRange = (0.01,25))


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
        plot2D = Plot2D.fromHisto(name = "power_%s_%s_vs_%s_%s_lumi_factor_%3.2f_level_%3.2f_nBinsTestStat_%i"%(test_statistic, WC1, WC2, ("truth" if args.truth else "predicted"), args.lumi_factor, level, args.nBinsTestStat), histos = [[power[test_statistic][level]]], texX = tex[WC1], texY = tex[WC2] )
        plotting.draw2D(plot2D, plot_directory = os.path.join( plot_directory, "binned"), histModifications = [lambda h:ROOT.gStyle.SetPalette(58)], logY = False, logX = False, logZ = True, copyIndexPHP=True, drawObjects = contour_objects, zRange = (0.01,1))

syncer.sync()
