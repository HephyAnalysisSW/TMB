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
argParser.add_argument("--plot_directory",     action="store",      default="BIT_VH_4",                 help="plot sub-directory")
argParser.add_argument("--model",              action="store",      default="WH_Spannowsky",                 help="plot sub-directory")
#argParser.add_argument('--debug',              action='store_true', help="Make debug plots?")
#argParser.add_argument('--feature_plots',      action='store_true', help="Feature plots?")
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

nEvents = 200000

features = model.getEvents(nEvents)
nEvents  = len(features)
weights  = model.getWeights(features, eft=model.default_eft_parameters)
print ("Created data set of size %i" % nEvents )

# normalize to SM event yield
lambda_expected_sm      = 100
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
bits = model.load()
predictions = { der:bits[der].vectorized_predict(features) for der in bits.keys() } 
def make_q( order, eft_norm=None, truth=False, **kwargs ):
    if eft_norm is not None:
        # we interpret kwargs as the direction in parameter space, i.e., q = n_theta.t + eft_norm*n_theta^T.s.n_theta
        eft     = copy.deepcopy(kwargs)
        norm    = sqrt(sum( np.array(kwargs.values())**2 ) )
        eft     = { k:eft[k]/norm for k in eft.keys() }
    else:
        # we interpret kwargs as the eft parameter point, i.e., q = theta.t + theta^T.s.theta
        eft_norm = 1
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
                prefac = eft_norm if order=="total" else 1
                result += prefac* .5*(eft[coeff1] - model.default_eft_parameters[coeff1])*(eft[coeff2] - model.default_eft_parameters[coeff2])*( weights[(coeff1,coeff2)]/weights[tuple()] if truth else predictions[tuple(sorted((coeff1,coeff2)))]) 
    return result

event_indices = np.arange(nEvents)
def make_toys( yield_per_toy, n_toys, lin=False, **kwargs):
    weights_      = make_weights(lin=lin, **kwargs) 
    biased_sample = np.random.choice( event_indices, size=10*nEvents,  p = weights_/np.sum(weights_) )

    return np.array( [ np.random.choice( biased_sample, size=n_observed ) for n_observed in np.random.poisson(yield_per_toy, n_toys) ])

tex = ROOT.TLatex()
tex.SetNDC()
tex.SetTextSize(0.04)

################################
## Plot test statistics  (toys)#
################################
#colors     = [ROOT.kBlack, ROOT.kBlue, ROOT.kGreen, ROOT.kMagenta, ROOT.kCyan, ROOT.kRed]
##extended   = True
#lumi_factor= 1
#n_toys     = 50000
#
#for WC, theta_vals in [
##            ("cHW",    [-.5, -.2, -.1, .1, .2, .5]), 
#            ("cHWtil", [-.5, -.2, -.1, .1, .2, .5]),
##            ("cHQ3",   [-.05, -.02, -.01, .01, .02, .05]),
#        ]: 
#
#    q_theta_given_theta = {}
#    q_theta_given_SM    = {}
#    for truth in [True, False]:
#    #for truth in [True]:
#        for test_statistic in ["lin", "quad", "total"]: 
#        #for test_statistic in ["lin"]: 
#            q_theta_given_theta [test_statistic]= {}
#            q_theta_given_SM    [test_statistic]= {}
#
#            truth_txt = "truth" if truth else "predicted"
#            print "Test statistic", test_statistic, "truth?", truth
#
#            for i_theta, theta in enumerate(theta_vals):
#                print "theta", theta
#                q_event = make_q( test_statistic, theta, truth=truth, **{WC:theta} )
#
#                #log_sigma_tot_ratio_subtraction = np.log(sigma_tot_ratio(theta)) if not extended else 0
#                q_theta_given_theta[test_statistic][theta] = np.array([np.sum( q_event[toy_] ) for toy_ in make_toys( lumi_factor*lambda_tot(**{WC:theta}), n_toys, **{WC:theta}) ])
#                q_theta_given_SM   [test_statistic][theta] = np.array([np.sum( q_event[toy_] ) for toy_ in make_toys( lumi_factor*lambda_tot(), n_toys ) ])
#
#                # calibration according to SM
#                if True:
#                    n = float(len(q_theta_given_SM[test_statistic][theta]))
#                    mean_q_theta_given_SM     = np.sum(q_theta_given_SM[test_statistic][theta])/n
#                    sigma_q_theta_given_SM    = sqrt( np.sum((q_theta_given_SM[test_statistic][theta] - mean_q_theta_given_SM)**2)/(n-1) ) 
#                    q_theta_given_SM[test_statistic][theta]    = (q_theta_given_SM[test_statistic][theta] - mean_q_theta_given_SM)/sigma_q_theta_given_SM
#                    q_theta_given_theta[test_statistic][theta] = (q_theta_given_theta[test_statistic][theta] - mean_q_theta_given_SM)/sigma_q_theta_given_SM
#
#            histos = []
#            quantile_lines  = []
#
#            all_vals = sum( [list(q_theta_given_theta[test_statistic][theta])+list(q_theta_given_SM[test_statistic][theta]) for theta in theta_vals], [] )
#            min_, max_ = min( all_vals ), max( all_vals )
#            binning = np.arange(min_, max_, (max_-min_)/100.)
#
#            for i_theta, theta in enumerate(theta_vals):
#
#                np_histo_SM    = np.histogram(q_theta_given_SM   [test_statistic][theta], bins=binning)
#                np_histo_theta = np.histogram(q_theta_given_theta[test_statistic][theta], bins=binning)
#                histo_SM       = stat_helpers.make_TH1F(np_histo_SM)
#                histo_theta    = stat_helpers.make_TH1F(np_histo_theta)
#
#                histo_SM.legendText    = "#color[%i]{p(q_{%s=%3.2f}|SM)}" % ( colors[i_theta], WC, theta)
#                histo_theta.legendText = "#color[%i]{p(q_{%s=%3.2f}|%s=%3.2f)}" % ( colors[i_theta], WC, theta, WC, theta )
#                histo_SM.style         = styles.lineStyle( colors[i_theta], dashed = True)
#                histo_theta.style      = styles.lineStyle( colors[i_theta] ) 
#
#                histos.append( histo_SM )
#                histos.append( histo_theta )
#
#                for x in np.quantile( q_theta_given_SM[test_statistic][theta], [.05, .95 ] ):
#                    quantile_lines.append( ROOT.TLine(x, 0, x, histo_SM.GetBinContent(histo_SM.FindBin(x))) )
#                    quantile_lines[-1].SetLineColor( colors[i_theta] )
#                    quantile_lines[-1].SetLineStyle( 7 )
#                for x in np.quantile( q_theta_given_theta[test_statistic][theta], [.05, .95 ] ):
#                    quantile_lines.append( ROOT.TLine(x, 0, x, histo_theta.GetBinContent(histo_theta.FindBin(x))) )
#                    quantile_lines[-1].SetLineColor( colors[i_theta] )
#
#            # Text on the plots
#            lines = [ 
#                    #  (0.25, 0.88, "#color[4]{%i%% qu. q_{BSM} = %3.2f}" % ( 100*(1-CL), q_theta_given_theta_1mCL[theta]) ),
#                    #  (0.25, 0.83, "#color[2]{q_{SM} = %3.2f}" % ( q_theta_SM ) ),
#                    #  (0.25, 0.78, "#theta_{current} = %5.4f" % theta ),
#                    ]
#            drawObjects = [ tex.DrawLatex(*line) for line in lines ]
#
#            plot = Plot.fromHisto( "test_stat_%s_%s_%s_lumi_factor_%3.2f"%(test_statistic, truth_txt, WC, lumi_factor), [[h] for h in histos], texX = "q (%s)"%truth_txt, texY = "Entries" )
#            plotting.draw( plot,
#                plot_directory = os.path.join( plot_directory, "test_statistics" ),
#                #ratio          = {'yRange':(0.6,1.4)} if len(plot.stack)>=2 else None,
#                logX = False, sorting = False,
#                legend         = ( (0.15,0.7,0.9,0.92),2),
#                drawObjects    =  quantile_lines + drawObjects,
#                copyIndexPHP   = True,
#                extensions     = ["png"], 
#              )            

###########################
## 1D+2D Plot derivatives #
###########################
#
#min_, max_ = {},{}
#for key in predictions.keys():
#    min_[key], max_[key] = stat_helpers.weighted_quantile( predictions[key], [.01,.99], sample_weight=weights[tuple()])
#
#h_1D_predicted = {}
#h_2D_predicted = {}
#h_1D_predicted = {}
#h_2D_predicted = {}
#h_1D_truth = {}
#h_2D_truth = {}
#h_1D_truth = {}
#h_2D_truth = {}
#h_2D = {}
#for i_der1, der1 in enumerate(model.derivatives):
#    if der1==tuple(): continue
#    h_1D_predicted[der1] = ROOT.TH1F("h_1D_predicted_%s"%("_".join(der1)), "h_1D_predicted_%s"%("_".join(der1)), 50, min_[der1], max_[der1])
#    h_1D_predicted[der1].legendText = "predicted %s"%("_".join(der1))
#    h_1D_predicted[der1].GetXaxis().SetTitle(",".join(der1))
#    h_1D_truth[der1] = ROOT.TH1F("h_1D_truth_%s"%("_".join(der1)), "h_1D_truth_%s"%("_".join(der1)), 50, min_[der1], max_[der1])
#    h_1D_truth[der1].legendText = "truth %s"%("_".join(der1))
#    h_1D_truth[der1].GetXaxis().SetTitle(",".join(der1))
#    h_2D[der1] = ROOT.TH2F("h_2D_truth_vs_predicted_%s"%("_".join(der1)), "h_2D_truth_vs_predicted_%s"%("_".join(der1)), 30, min_[der1], max_[der1], 30, min_[der1], max_[der1])
#    for i_der2, der2 in enumerate(model.derivatives):
#        if der2==tuple(): continue
#        if i_der1>=i_der2: continue
#        h_2D_predicted[(der1,der2)] = ROOT.TH2F("h_2D_predicted_%s_VS_%s"%("_".join(der1), "_".join(der2)), "h_2D_predicted_%s_VS_%s"%("_".join(der1), "_".join(der2)), 30, min_[der1], max_[der1], 30, min_[der2], max_[der2])
#        h_2D_predicted[(der1,der2)].GetXaxis().SetTitle(",".join(der1))
#        h_2D_predicted[(der1,der2)].GetYaxis().SetTitle(",".join(der2))
#        h_2D_truth[(der1,der2)] = ROOT.TH2F("h_2D_truth_%s_VS_%s"%("_".join(der1), "_".join(der2)), "h_2D_truth_%s_VS_%s"%("_".join(der1), "_".join(der2)), 30, min_[der1], max_[der1], 30, min_[der2], max_[der2])
#        h_2D_truth[(der1,der2)].GetXaxis().SetTitle(",".join(der1))
#        h_2D_truth[(der1,der2)].GetYaxis().SetTitle(",".join(der2))
#
#for i_event in range(nEvents):
#    for i_der1, der1 in enumerate(model.derivatives):
#        if der1==tuple(): continue
#        h_1D_predicted[der1].Fill(predictions[der1][i_event], weights[()][i_event])
#        h_1D_truth[der1].Fill(weights[der1][i_event]/weights[()][i_event], weights[()][i_event])
#        h_2D[der1].Fill(weights[der1][i_event]/weights[()][i_event], predictions[der1][i_event], weights[()][i_event])
#        for i_der2, der2 in enumerate(model.derivatives):
#            if der2==tuple(): continue
#            if i_der1>=i_der2: continue
#            h_2D_predicted[(der1,der2)].Fill(predictions[der1][i_event], predictions[der2][i_event], weights[()][i_event])
#            h_2D_truth[(der1,der2)].Fill(weights[der1][i_event]/weights[()][i_event], weights[der2][i_event]/weights[()][i_event], weights[()][i_event])
#
#for i_der1, der1 in enumerate(model.derivatives):
#    if der1==tuple(): continue
#
#    h_1D_predicted[der1].style = styles.lineStyle( ROOT.kRed )
#    plot = Plot.fromHisto( "derivative_1D_%s"%("_".join(der1)), [[h_1D_truth[der1]], [h_1D_predicted[der1]]], texX = ",".join(der1), texY = "Events" )
#    plotting.draw( plot,
#        plot_directory = os.path.join( plot_directory, "derivatives" ),
#        #ratio          = {'yRange':(0.6,1.4)} if len(plot.stack)>=2 else None,
#        logX = False, sorting = False,
#        #legend         = ( (0.15,0.7,0.9,0.92),2),
#        copyIndexPHP   = True,
#        extensions     = ["png"], 
#      )            
#    plot = Plot2D.fromHisto( "derivative_2D_truth_VS_predicted_%s"%("_".join(der1)), [[h_2D[der1]]], texX = "truth", texY = "predicted" )
#    plotting.draw2D( plot,
#        plot_directory = os.path.join( plot_directory, "derivatives" ),
#        #ratio          = {'yRange':(0.6,1.4)} if len(plot.stack)>=2 else None,
#        logX = False,
#        copyIndexPHP   = True,
#        extensions     = ["png"], 
#      )            
#
#    for i_der2, der2 in enumerate(model.derivatives):
#        if der2==tuple(): continue
#        if i_der1>=i_der2: continue
#        plot = Plot2D.fromHisto( "derivative_2D_predicted_%s_VS_%s"%("_".join(der1), "_".join(der2)), [[h_2D_predicted[(der1,der2)]]], texX = ",".join(der1), texY = ",".join(der2) )
#        plotting.draw2D( plot,
#            plot_directory = os.path.join( plot_directory, "derivatives" ),
#            #ratio          = {'yRange':(0.6,1.4)} if len(plot.stack)>=2 else None,
#            logX = False,
#            copyIndexPHP   = True,
#            extensions     = ["png"], 
#          )            
#        plot = Plot2D.fromHisto( "derivative_2D_truth_%s_VS_%s"%("_".join(der1), "_".join(der2)), [[h_2D_truth[(der1,der2)]]], texX = ",".join(der1), texY = ",".join(der2) )
#        plotting.draw2D( plot,
#            plot_directory = os.path.join( plot_directory, "derivatives" ),
#            #ratio          = {'yRange':(0.6,1.4)} if len(plot.stack)>=2 else None,
#            logX = False,
#            copyIndexPHP   = True,
#            extensions     = ["png"], 
#          )            

###############################################
# Plot quantiles of test statistics (TGraphs) #
###############################################

colors   = [ ROOT.kRed, ROOT.kBlue, ROOT.kBlack, ROOT.kBlue, ROOT.kRed]
levels   = [ .05,       .32,        .5,          .68,        .95      ]
#extended = True
n_toys   = 50000

#WC       = "cHQ3"
#theta_min   = -.03
#theta_max   =  .03

#WC       = "cHW"
#theta_min   = -.5
#theta_max   =  .3

WC       = "cHWtil"
theta_min   = -.3
theta_max   =  .3

Nbins       = 60
UL = {}

for WC, theta_min, theta_max in [
            ("cHWtil", -1.5,  .5),
            ("cHW",    -1.5,  .5), 
            ("cHQ3",   -.06, .06),
        ]: 

    power = {}
    #lumi_factors = [1., .5, 2.]
    lumi_factors = [.05, .1, .2, .5, 1.]
    for lumi_factor in lumi_factors:

        UL[lumi_factor] = {level:{} for level in levels if level<0.5}
        power[lumi_factor] = {}

        theta_vals      = theta_vals = np.arange(theta_min, theta_max, (theta_max-theta_min)/Nbins) 

        #for truth in [True, False]:
        for truth in [True]:

            power[lumi_factor][truth] = {}
            test_statistics = ["total", "lin", "quad"]
            for test_statistic in test_statistics: 

                truth_txt = "truth" if truth else "predicted"
                print "Test statistic", test_statistic, "truth?", truth

                power[lumi_factor][truth][test_statistic] = ROOT.TGraph(len(theta_vals))

                tgraphs_theta   = { level: ROOT.TGraph(len(theta_vals)) for level in levels }
                tgraphs_SM      = { level: ROOT.TGraph(len(theta_vals)) for level in levels }

                min_, max_ = float('inf'), -float('inf')

                sm_toys = make_toys( lumi_factor*lambda_tot(), n_toys ) 
                for i_theta, theta in enumerate( theta_vals ):

                    q_event = make_q( test_statistic, theta, truth=truth, **{WC:1} )

                    #log_sigma_tot_ratio_subtraction = np.log(sigma_tot_ratio(theta)) if not extended else 0
                    q_theta_given_SM    = np.array([np.sum( q_event[toy_] ) for toy_ in sm_toys ])
                    q_theta_given_theta = np.array([np.sum( q_event[toy_] ) for toy_ in make_toys( lumi_factor*lambda_tot(**{WC:theta}), n_toys, **{WC:theta}) ])

                    # calibration according to SM
                    if True:
                        n = float(len(q_theta_given_SM))
                        mean_q_theta_given_SM     = np.sum(q_theta_given_SM)/n
                        sigma_q_theta_given_SM    = sqrt( np.sum((q_theta_given_SM - mean_q_theta_given_SM)**2)/(n-1) ) 
                        q_theta_given_SM    = (q_theta_given_SM - mean_q_theta_given_SM)/sigma_q_theta_given_SM
                        q_theta_given_theta = (q_theta_given_theta - mean_q_theta_given_SM)/sigma_q_theta_given_SM

                    quantiles_theta = np.quantile( q_theta_given_theta, levels )
                    quantiles_SM    = np.quantile( q_theta_given_SM, levels )
                    size_  = np.sum(np.histogram( q_theta_given_SM, quantiles_SM)[0])/float(n_toys)
                    power_ = 1-np.sum(np.histogram( q_theta_given_theta, quantiles_SM)[0])/float(n_toys)
                    print "theta", round(theta,3), "power", round(power_,3), "size", size_, test_statistic, WC, "truth",truth
                    power[lumi_factor][truth][test_statistic].SetPoint( i_theta, theta, power_ ) 

                    for quantile, level in zip( quantiles_theta, levels ):
                        if not level<0.5: continue
                        sm_toy_fraction_below_level = np.count_nonzero( q_theta_given_SM<=quantile )/float(len(sm_toys))
                        if not UL[lumi_factor][level].has_key(test_statistic) and sm_toy_fraction_below_level>.5:
                            UL[lumi_factor][level][test_statistic] = {"theta":theta, 'frac':sm_toy_fraction_below_level}

                        #print "SM toys below level Q(%i)=%3.2f: %3.2f at theta %3.2f"% ( 100*level, quantile, sm_toy_fraction_below_level, theta)

                    min__ = min(list(quantiles_theta)+list(quantiles_SM))
                    max__ = max(list(quantiles_theta)+list(quantiles_SM))
                    if min__<min_: min_=min__
                    if max__>max_: max_=max__

                    [ tgraphs_theta[level].SetPoint( i_theta, quantile, theta) for level, quantile in zip( levels, quantiles_theta ) ]
                    [ tgraphs_SM   [level].SetPoint( i_theta, quantile, theta) for level, quantile in zip( levels, quantiles_SM ) ]

                # TGraph contour plots
                c1 = ROOT.TCanvas()
                ROOT.gStyle.SetOptStat(0)
                c1.SetTitle("")
                h_empty = ROOT.TH1F("x","",1,min_ - 0.1*(max_-min_),max_+0.1*(max_-min_))
                h_empty.Draw()
                h_empty.GetYaxis().SetRangeUser(min(theta_vals), 1.7*max(theta_vals))
                h_empty.Draw()
                l1 = ROOT.TLegend(0.25, 0.77, 0.65, 0.92)
                l1.SetFillStyle(0)
                l1.SetShadowColor(ROOT.kWhite)
                l1.SetBorderSize(0)

                l2 = ROOT.TLegend(0.6,  0.77, 0.95, 0.92)
                l2.SetFillStyle(0)
                l2.SetShadowColor(ROOT.kWhite)
                l2.SetBorderSize(0)
                for i_level, level in enumerate(levels):
                    l1.AddEntry( tgraphs_SM[level], "Q(%i%%) SM"%(100*level) )
                    tgraphs_SM[level].SetTitle("")
                    tgraphs_SM[level].GetXaxis().SetTitle("q")
                    tgraphs_SM[level].GetYaxis().SetTitle("#theta")
                    tgraphs_SM[level].SetLineColor( colors[i_level] )
                    tgraphs_SM[level].SetLineStyle( 7 )
                    tgraphs_SM[level].SetMarkerStyle(0 )
                    tgraphs_SM[level].SetMarkerColor( colors[i_level] )
                    tgraphs_SM[level].Draw("L") 

                    l2.AddEntry( tgraphs_theta[level], "Q(%i%%) BSM"%(100*level) )
                    tgraphs_theta[level].SetTitle("")
                    tgraphs_theta[level].GetXaxis().SetTitle("q")
                    tgraphs_theta[level].GetYaxis().SetTitle("#theta")
                    tgraphs_theta[level].SetLineColor( colors[i_level] )
                    tgraphs_theta[level].SetMarkerStyle(0 )
                    tgraphs_theta[level].SetMarkerColor( colors[i_level] )
                    tgraphs_theta[level].Draw("L")

                l1.Draw()
                l2.Draw()

                c1.RedrawAxis()
                c1.Update()
                c1.Print(os.path.join( plot_directory, "test_statistics/contours_%s_%s_%s_lumi_factor_%3.2f.png"% ( WC, test_statistic, ("truth" if truth else "predicted"), lumi_factor) ))

        # power plots
        c1 = ROOT.TCanvas()
        ROOT.gStyle.SetOptStat(0)
        c1.SetTitle("")
        h_empty = ROOT.TH1F("x","",1,theta_min,theta_max)
        h_empty.Draw()
        h_empty.GetYaxis().SetRangeUser(0,1)
        h_empty.GetXaxis().SetTitle("#theta")
        h_empty.GetYaxis().SetTitle("power")
        h_empty.Draw()

        l = ROOT.TLegend(0.2, 0.8-0.05*(len(test_statistics))*len(power[lumi_factor].keys()), 0.4, 0.8)
        l.SetFillStyle(0)
        l.SetShadowColor(ROOT.kWhite)
        l.SetBorderSize(0)
        #legends = []
        for i_test_statistic, test_statistic in enumerate(reversed(test_statistics)):
            #l = ROOT.TLegend(0.2 + i_test_statistic*(0.65/3.), 0.77, (i_test_statistic+1)*(0.65/3.), 0.92)
            #legends.append(l)
            for truth in [True, False]:
                if not power[lumi_factor].has_key(truth): continue
                tgr = power[lumi_factor][truth][test_statistic]
                tgr.SetTitle("")
                tgr.GetXaxis().SetTitle("#theta")
                tgr.GetYaxis().SetTitle("power")
                tgr.SetLineColor( colors[i_test_statistic] )
                tgr.SetLineStyle( 7 if not truth else 1 )
                tgr.SetMarkerStyle(0 )
                tgr.SetMarkerColor( colors[i_test_statistic] )
                l.AddEntry( tgr, test_statistic+" (%s)"%( "truth" if truth else "pred.") )
                tgr.Draw("L") 

        #for l in legends:
        l.Draw()                

        c1.RedrawAxis()
        c1.Update()
        c1.Print(os.path.join( plot_directory, "test_statistics/power_%s_lumi_factor_%3.2f.png"%(WC, lumi_factor) ))
        syncer.sync()
