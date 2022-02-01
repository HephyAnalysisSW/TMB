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

# BIT
sys.path.insert(0,os.path.expandvars("$CMSSW_BASE/src/BIT"))
from BoostedInformationTree import BoostedInformationTree

# User
import user

# Parser
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument("--plot_directory",     action="store",      default="BIT_VH",                 help="plot sub-directory")
argParser.add_argument("--model",              action="store",      default="WH_Spannowsky",                 help="plot sub-directory")
#argParser.add_argument('--overwrite',          action='store_true', help="Overwrite output?")
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


################
### Plot model #
################
##
##colors = [ROOT.kBlack, ROOT.kBlue, ROOT.kGreen, ROOT.kMagenta, ROOT.kCyan, ROOT.kRed]
##plot_features = [ 
##    {'name':'x',   'texX':'x',    'binning':[50,0,10],   },
##    {'name':'phi', 'texX':'#phi', 'binning':[50,-pi,pi], },
##    ]
##
##h = {}
##thetas = [.05, .1, .15, .2]
##n_events_plot = 10**5
##features, weights = get_sampled_dataset(n_events_plot)
##for i_theta, theta in enumerate(thetas):
##    name = ("theta_%3.2f"%theta ).replace('.','p').replace('-','m') 
##    h[theta] = {}
##    for i_feature, feature in enumerate(plot_features):
##        h[theta][i_feature]              = ROOT.TH1F(name+'_'+feature['name'], name+'_'+feature['name'], *plot_features[i_feature]['binning'] )
##        h[theta][i_feature].style        = styles.lineStyle( colors[i_theta], width=2, dashed=False )
##        h[theta][i_feature].legendText   = "#theta = %3.2f"%theta 
##
##for i_event, event in enumerate(features):
##    if i_event%10**4==0: print "At %i/%i"%( i_event, n_events_plot )
##    for i_feature, feature in enumerate(plot_features):
##        for i_theta, theta in enumerate(thetas):
##            h[theta][i_feature].Fill(event[i_feature], weights[()][i_event] + theta*weights[(0,)][i_event] + 0.5*theta**2*weights[(0,0)][i_event])
##
##for i_feature, feature in enumerate(plot_features):
##    histos = [[h[theta][i_feature]] for theta in thetas]
##    plot   = Plot.fromHisto( feature['name'],  histos, texX=feature['texX'], texY="a.u." )
##
##    for log in [True, False]:
##
##        # Add subdirectory for lin/log plots
##        plot_directory_ = os.path.join( plot_directory, "newman", "log" if log else "lin" )
##        plotting.draw( plot,
##                       plot_directory = plot_directory_,
##                       logX = False, logY = log, sorting = False,
##                       yRange = "auto",
##                       ratio = None,
###                       drawObjects = drawObjects( lumi, offset=titleOffset ),
##                        legend=(0.2,0.7,0.9,0.9),
##                       #histModifications = histModifications,
##                       copyIndexPHP = True,
##                       )
#
## we simulate a large number of events
#dataset_size      = 10**6
#features, weights = get_sampled_dataset(dataset_size)
#event_indices     = np.arange(dataset_size)
#
## normalize to SM xsec
#sigma_tot_sm      = 1
##weights[(0,)]     = weights[(0,)]*sigma_tot_sm/np.sum(weights[tuple()])
##weights[(0,0)]    = weights[(0,0)]*sigma_tot_sm/np.sum(weights[tuple()])
##weights[tuple()]  = weights[tuple()]*sigma_tot_sm/np.sum(weights[tuple()])
#
## total xsec ratio
#sigma_tot_ratio_lin  = np.sum(weights[(0,)])/np.sum(weights[tuple()])
#sigma_tot_ratio_quad = np.sum(weights[(0,0)])/np.sum(weights[tuple()])
#
## xsec ratio
#def sigma_tot_ratio( theta, lin=False):
#    return 1 + theta*sigma_tot_ratio_lin + ( 0.5*theta**2*sigma_tot_ratio_quad if not lin else 0 )
## total xsec
#def sigma_tot( theta, lin=False):
#    return sigma_tot_sm*sigma_tot_ratio( theta, lin=lin ) 
#
#
## expansion: (1 + theta*q_event_lin)**2 + (theta*q_event_quad)**2
#q_event_lin  = np.cos(n*features[:,1])*np.exp(features[:,0]/(2*x12m))
#q_event_quad = np.sin(n*features[:,1])*np.exp(features[:,0]/(2*x12m))
#
##def wrapper( i_toy, toys, plot = plot, test_statistic = "total", extended = False, verbose = True):
##for i_toy in range(5):
#
#def make_toys( yield_per_toy, theta, n_toys):
#    weights_ =  weights[tuple()]+theta*weights[(0,)]+.5*theta**2*weights[(0,0)]
#    biased_sample = np.random.choice( event_indices, size=dataset_size,  p = weights_/np.sum(weights_) )
#
#    return np.array( [ np.random.choice( biased_sample, size=n_observed ) for n_observed in np.random.poisson(yield_per_toy, n_toys) ])
#
#tex = ROOT.TLatex()
#tex.SetNDC()
#tex.SetTextSize(0.04)
#
#################################
### Plot test statistics  (toys)#
#################################
##colors   = [ROOT.kBlack, ROOT.kBlue, ROOT.kGreen, ROOT.kMagenta, ROOT.kCyan, ROOT.kRed]
##extended = True
##lumi     = 137
##n_toys   = 50000
##theta_SM = 0
##theta_vals = [1., .5, .1, 0.05]
##
###for lumi in [ 13.70/5., 13.70/2, 13.70 , 2*13.70, 5*13.70, 10*13.70 ]: 
##for lumi in [ 13.70/2 ]: 
##    q_theta_given_theta = {}
##    q_theta_given_SM    = {}
##    for test_statistic in ["lin", "quad", "total"]: 
##        q_theta_given_theta [test_statistic]= {}
##        q_theta_given_SM    [test_statistic]= {}
##
##        print "Test statistic", test_statistic
##
##        histos = []
##        quantile_lines  = []
##        for i_theta, theta in enumerate(theta_vals):
##            print "theta", theta
##            if test_statistic == "quad":
##                q_event = np.log( (q_event_lin**2+q_event_quad**2) )
##            elif test_statistic == "lin":
##                q_event = 1./theta * np.log( (1 + theta*q_event_lin)**2 )
##            elif test_statistic == "total":
##                q_event = 1./theta * np.log( (1 + theta*q_event_lin)**2 + (theta*q_event_quad)**2)
##            else:
##                raise RuntimeError( "Unknwon test statistc %s" % test_statistic )
##
##            log_sigma_tot_ratio_subtraction = np.log(sigma_tot_ratio(theta)) if not extended else 0
##            q_theta_given_theta[test_statistic][theta] = np.array([np.sum( q_event[toy_] - log_sigma_tot_ratio_subtraction ) for toy_ in make_toys( lumi*sigma_tot(theta),    theta,    n_toys ) ])
##            q_theta_given_SM   [test_statistic][theta] = np.array([np.sum( q_event[toy_] - log_sigma_tot_ratio_subtraction ) for toy_ in make_toys( lumi*sigma_tot(theta_SM), theta_SM, n_toys ) ])
##
##        for i_theta, theta in enumerate(theta_vals):
##
##            if i_theta == 0:
##                all_vals = sum( [list(q_theta_given_theta[k])+list(q_theta_given_SM[k]) for k in q_theta_given_theta.keys()], [])
##                min_, max_ = min( all_vals ), max( all_vals )
##                binning = np.arange(min_, max_, (max_-min_)/100.)
##
##            #np_histo_all   = np.histogram(q_theta_given_theta+q_theta_given_SM, 100)
##            #histo_all      = stat_helpers.make_TH1F(np_histo_all)
##            #binning        = np_histo_all[1] 
##
##            np_histo_SM    = np.histogram(q_theta_given_SM   [test_statistic][theta],    bins=binning)
##            np_histo_theta = np.histogram(q_theta_given_theta[test_statistic][theta], bins=binning)
##            histo_SM       = stat_helpers.make_TH1F(np_histo_SM)
##            histo_theta    = stat_helpers.make_TH1F(np_histo_theta)
##
##            histo_SM.legendText    = "#color[%i]{p(q_{#theta}|#theta)}, #theta =%3.2f" % ( colors[i_theta], theta_SM )
##            histo_theta.legendText = "#color[%i]{p(q_{#theta}|#theta_{SM})}, #theta = %3.2f" % ( colors[i_theta], theta )
##            histo_SM.style         = styles.lineStyle( colors[i_theta], dashed = True)
##            histo_theta.style      = styles.lineStyle( colors[i_theta] ) 
##
##            #histo_all = histo_all
##            #histo_all.legendText = None
##            #histo_all.style      = styles.invisibleStyle()
##
##            #histos.append( histo_all )
##            histos.append( histo_SM )
##            histos.append( histo_theta )
##
##            for x in np.quantile( q_theta_given_SM[test_statistic][theta], [.05, .95 ] ):
##                quantile_lines.append( ROOT.TLine(x, 0, x, histo_SM.GetBinContent(histo_SM.FindBin(x))) )
##                quantile_lines[-1].SetLineColor( colors[i_theta] )
##                quantile_lines[-1].SetLineStyle( 7 )
##            for x in np.quantile( q_theta_given_theta[test_statistic][theta], [.05, .95 ] ):
##                quantile_lines.append( ROOT.TLine(x, 0, x, histo_theta.GetBinContent(histo_theta.FindBin(x))) )
##                quantile_lines[-1].SetLineColor( colors[i_theta] )
##
##            # Text on the plots
##            lines = [ 
##                    #  (0.25, 0.88, "#color[4]{%i%% qu. q_{BSM} = %3.2f}" % ( 100*(1-CL), q_theta_given_theta_1mCL[theta]) ),
##                    #  (0.25, 0.83, "#color[2]{q_{SM} = %3.2f}" % ( q_theta_SM ) ),
##                    #  (0.25, 0.78, "#theta_{current} = %5.4f" % theta ),
##                    ]
##            drawObjects = [ tex.DrawLatex(*line) for line in lines ]
##
##        plot = Plot.fromHisto( "test_stat_%s_n%i_nEvents_%3.2f"%(test_statistic, n, lumi*sigma_tot(theta_SM)), [[h] for h in histos], texX = "q", texY = "Entries" )
##        plotting.draw( plot,
##            plot_directory = os.path.join( plot_directory, "newman_new2_q" ),
##            #ratio          = {'yRange':(0.6,1.4)} if len(plot.stack)>=2 else None,
##            logX = False, sorting = False,
##            legend         = ( (0.15,0.7,0.9,0.92),2),
##            drawObjects    =  quantile_lines + drawObjects,
##            copyIndexPHP   = True,
##            extensions     = ["png"], 
##          )            
#
################################################
## Plot quantiles of test statistics (TGraphs) #
################################################
#
#colors   = [ ROOT.kRed, ROOT.kBlue, ROOT.kBlack, ROOT.kBlue, ROOT.kRed]
#levels   = [ .05,       .32,        .5,          .68,        .95      ]
#extended = True
#n_toys   = 50000
#theta_SM = 0
#theta_max= .5
#Nbins    = 20
#UL = {}
#for lumi in 137*np.array([ 1/50., 1/20. , 1/10., 1/5., 1/2., 1., 2., 5., 10., 20., 50. ]): 
#    UL[lumi] = {level:{} for level in levels if level<0.5}
#    for test_statistic in ["total", "lin", "quad"]: 
#
#        print "Test statistic", test_statistic
#
#        theta_vals      = np.arange( theta_SM, theta_max, (theta_max-theta_SM)/Nbins)
#        tgraphs_theta   = { level: ROOT.TGraph(len(theta_vals)) for level in levels }
#        tgraphs_SM      = { level: ROOT.TGraph(len(theta_vals)) for level in levels }
#
#        min_, max_ = float('inf'), -float('inf')
#
#        sm_toys = make_toys( lumi*sigma_tot(theta_SM), theta_SM, n_toys ) 
#        for i_theta, theta in enumerate( theta_vals ):
#
#            if test_statistic == "quad":
#                q_event = np.log( (q_event_lin**2+q_event_quad**2) )
#            elif test_statistic == "lin":
#                q_event = 1./theta * np.log( (1 + theta*q_event_lin)**2 )                           if theta>0.001 else 2*q_event_lin - theta*q_event_lin**2
#            elif test_statistic == "total":
#                q_event = 1./theta * np.log( (1 + theta*q_event_lin)**2 + (theta*q_event_quad)**2)  if theta>0.001 else 2*q_event_lin + theta*(q_event_lin**2+q_event_quad**2)
#            else:
#                raise RuntimeError( "Unknwon test statistc %s" % test_statistic )
#
#            log_sigma_tot_ratio_subtraction = np.log(sigma_tot_ratio(theta)) if not extended else 0
#            q_theta_given_SM    = np.array([np.sum( q_event[toy_] - log_sigma_tot_ratio_subtraction ) for toy_ in sm_toys ])
#            q_theta_given_theta = np.array([np.sum( q_event[toy_] - log_sigma_tot_ratio_subtraction ) for toy_ in make_toys( lumi*sigma_tot(theta), theta, n_toys ) ])
#
#            quantiles_theta = np.quantile( q_theta_given_theta, levels )
#            quantiles_SM    = np.quantile( q_theta_given_SM, levels )
#
#            for quantile, level in zip( quantiles_theta, levels ):
#                if not level<0.5: continue
#                sm_toy_fraction_below_level = np.count_nonzero( q_theta_given_SM<=quantile )/float(len(sm_toys))
#                if not UL[lumi][level].has_key(test_statistic) and sm_toy_fraction_below_level>.5:
#                    UL[lumi][level][test_statistic] = {"theta":theta, 'frac':sm_toy_fraction_below_level}
#
#                #print "SM toys below level Q(%i)=%3.2f: %3.2f at theta %3.2f"% ( 100*level, quantile, sm_toy_fraction_below_level, theta)
#
#            min__ = min(list(quantiles_theta)+list(quantiles_SM))
#            max__ = max(list(quantiles_theta)+list(quantiles_SM))
#            if min__<min_: min_=min__
#            if max__>max_: max_=max__
#
#            [ tgraphs_theta[level].SetPoint( i_theta, quantile, theta) for level, quantile in zip( levels, quantiles_theta ) ]
#            [ tgraphs_SM   [level].SetPoint( i_theta, quantile, theta) for level, quantile in zip( levels, quantiles_SM ) ]
#
#        c1 = ROOT.TCanvas()
#        ROOT.gStyle.SetOptStat(0)
#        c1.SetTitle("")
#        h_empty = ROOT.TH1F("x","",1,min_ - 0.1*(max_-min_),max_+0.1*(max_-min_))
#        l1 = ROOT.TLegend(0.15, 0.7, 0.65, 0.85)
#        l1.SetFillStyle(0)
#        l1.SetShadowColor(ROOT.kWhite)
#        l1.SetBorderSize(0)
#
#        l2 = ROOT.TLegend(0.5, 0.70, 0.75, 0.85)
#        l2.SetFillStyle(0)
#        l2.SetShadowColor(ROOT.kWhite)
#        l2.SetBorderSize(0)
#        h_empty.Draw()
#        h_empty.GetYaxis().SetRangeUser(min(theta_vals), 1.4*max(theta_vals))
#        h_empty.Draw()
#        for i_level, level in enumerate(levels):
#            l1.AddEntry( tgraphs_SM[level], "Q(%i%%) SM"%(100*level) )
#            tgraphs_SM[level].SetTitle("")
#            tgraphs_SM[level].GetXaxis().SetTitle("q")
#            tgraphs_SM[level].GetYaxis().SetTitle("#theta")
#            tgraphs_SM[level].SetLineColor( colors[i_level] )
#            tgraphs_SM[level].SetLineStyle( 7 )
#            tgraphs_SM[level].Draw("L") 
#
#            l2.AddEntry( tgraphs_theta[level], "Q(%i%%) BSM"%(100*level) )
#            tgraphs_theta[level].SetTitle("")
#            tgraphs_theta[level].GetXaxis().SetTitle("q")
#            tgraphs_theta[level].GetYaxis().SetTitle("#theta")
#            tgraphs_theta[level].SetLineColor( colors[i_level] )
#            tgraphs_theta[level].Draw("L")
#
#        l1.Draw()
#        l2.Draw()
#
#        c1.RedrawAxis()
#        c1.Update()
#        c1.Print(os.path.join( plot_directory, "newman_new2_q/contours_%s_extended_%i_%3.2f.png"% ( test_statistic, extended, lumi*sigma_tot(theta_SM) )))
#
###for i_toy, toy in enumerate( toys_SM[:5] ):
##def wrapper( i_toy, toys, plot = True, test_statistic = "total", extended = True, verbose = True):
##    q_theta_given_theta_1mCL = {}
##
##    import  Analysis.Tools.syncer as syncer
##    theta_current    = theta_UL_start 
##
##    #print "Toy", i_toy
##    toy = toys[i_toy]
##    i_iter    = 0
##    while True:
##
##        # This is the second toy -> again, we guess theta
##        if i_iter==1:
##            delta_q_previous = q_theta_given_theta_1mCL[theta_current] - q_theta_SM 
##            theta_previous   = theta_UL_start
##            theta_current    = theta_UL_start/2. 
##
##        # We are at the third iteration. Time to update.
##        elif i_iter>1:
##
##            delta_q_current = q_theta_given_theta_1mCL[theta_current] - q_theta_SM
##             
##            #print q_theta_SM,  q_theta_given_theta_1mCL[theta_current]
##            if abs(theta_current - theta_previous)<=tolerance:
##                # Converged!
##                if verbose: 
##                    print "Toy", i_toy, "done."
##                    print "theta_current",  theta_current, "delta_q_current", delta_q_current, "i_toy", i_toy, "i_iter", i_iter
##                    print
##                break
##                #return [ theta_current, delta_q_current, i_toy, i_iter]
##            elif i_iter>=max_iter:
##                if verbose: 
##                    print "Toy", i_toy, "max_iter %i reached" % max_iter
##                    print "theta_current", theta_current, "theta_previous", theta_previous, "delta_q_current", delta_q_current
##                    print
##                break
##
##            # Newton step
##            k = (theta_current-theta_previous)/(delta_q_current-delta_q_previous)
##            if k>2:
##                if verbose: print "k=%3.2f too large, set it to 2." % k
##                k=2
##            if k<-2:
##                if verbose: print "k=%3.2f too small, set it to -2." % k
##                k=-2
##            theta_current, theta_previous, delta_q_previous = theta_current - delta_q_current*k, theta_current, delta_q_current 
##            # If the predicted value is negative, cool down and half the distance (to zero)
##            if theta_current<0:
##                if verbose: print "Predicted negative, cooling down."
##                theta_current = theta_previous/2.
##
##        # compute value of q_theta test statistic for all events 
##        #q_event = 1/theta_current * np.log( (1 + theta_current*q_event_lin)**2 )
##        if test_statistic == "quad":
##            q_event = np.log( theta_current**2*(q_event_lin**2+q_event_quad**2) )
##        elif test_statistic == "lin":
##            q_event = np.log( (1 + theta_current*q_event_lin)**2 )
##        elif test_statistic == "total":
##            q_event = np.log( (1 + theta_current*q_event_lin)**2 + (theta_current*q_event_quad)**2)
##        else:
##            raise RuntimeError( "Unknwon test statistc %s" % test_statistic )
##
##        log_sigma_tot_ratio_subtraction = np.log(sigma_tot_ratio(theta_current)) if not extended else 0
##        q_theta_SM          = np.sum( q_event[toy] - log_sigma_tot_ratio_subtraction )
##        #print i_toy, "subtracting", theta_current, log_sigma_tot_ratio_subtraction, q_theta_SM, np.sum( q_event[toy])
##
##        # compute 1-CL quantile of the q_theta test statistic under the theta hypothesis
##        if not q_theta_given_theta_1mCL.has_key(theta_current) or plot:
##            # compute distribution of test statistic q_theta under theta 
##            q_theta_given_theta = np.array([np.sum( q_event[toy_] - log_sigma_tot_ratio_subtraction ) for toy_ in make_toys( lumi*sigma_tot(theta_current), theta_current, n_toys ) ])
##            q_theta_given_theta_1mCL[theta_current] = np.quantile( q_theta_given_theta, 1-CL)
##            if plot:
##
##                # Text on the plots
##                lines = [ 
##                          (0.25, 0.88, "#color[4]{%i%% qu. q_{BSM} = %3.2f}" % ( 100*(1-CL), q_theta_given_theta_1mCL[theta_current]) ),
##                          (0.25, 0.83, "#color[2]{q_{SM} = %3.2f}" % ( q_theta_SM ) ),
##                          (0.25, 0.78, "#theta_{current} = %5.4f" % theta_current ),
##                        ]
##                drawObjects = [ tex.DrawLatex(*line) for line in lines ]
##
##                np_histo      = np.histogram(q_theta_given_theta, 100)
##                histo         = stat_helpers.make_TH1F(np_histo)
##                quantile_line = ROOT.TLine(q_theta_given_theta_1mCL[theta_current], 0, q_theta_given_theta_1mCL[theta_current], histo.GetMaximum())
##                quantile_line.SetLineColor(ROOT.kRed)
##                histo.style = styles.lineStyle( ROOT.kRed )
##                histo.legendText = "#color[2]{q_{BSM}}" 
##
##                q_theta_SM_line = ROOT.TLine(q_theta_SM, 0, q_theta_SM, histo.GetMaximum())
##                q_theta_SM_line.SetLineColor(ROOT.kBlue)
##
##                q_theta_given_SM = np.array([np.sum( q_event[toy_]-log_sigma_tot_ratio_subtraction ) for toy_ in toys_SM])
##                histo_SM         = stat_helpers.make_TH1F(np.histogram(q_theta_given_SM, bins=np_histo[1]))
##                histo_SM.legendText = "#color[4]{q_{SM}}" 
##                histo_SM.style = styles.lineStyle( ROOT.kBlue )
##
##                plot = Plot.fromHisto( "itoy_%03i_iter_%02i"%( i_toy, i_iter ),[[histo], [histo_SM]], texX = "q", texY = "Entries" )
##                plotting.draw( plot,
##                    plot_directory = os.path.join( plot_directory, "newman_new2_%3.2f"% ( lumi*sigma_tot(theta_SM)) ),
##                    #ratio          = {'yRange':(0.6,1.4)} if len(plot.stack)>=2 else None,
##                    logX = False, sorting = False,
##                    legend         = (0.65,0.78,0.9,0.9),
##                    drawObjects    =  [quantile_line, q_theta_SM_line] + drawObjects,
##                    copyIndexPHP   = True,
##                    extensions     = ["png"], 
##                  )            
##        i_iter +=1
##    return [ theta_current, delta_q_current, i_toy, i_iter]
##
##
##UL                       = {}
##
###UL = [ wrapper(i_toy, toys_SM, plot=True) for i_toy in range(100) ]
##
##
##from multiprocessing import Pool
##
##theta_SM    = 0
##
##CL          = 0.95
##n_toys      = 50000
##theta_UL_start = 1
##tolerance      = 0.001
##max_iter       = 50
##test_statistic = "total"
##
##plot           = True
##extended       = False
##verbose        = False
##
##n_toys_UL      = 20
##UL             = {}
###assert False, ""
##
##for lumi in [ 13.70 , 2*13.70, 5*13.70, 10*13.70, 20*13.70, 50*13.70 ]: 
###for lumi in [  5*13.70, 10*13.70, 20*13.70, 50*13.70 ]: 
###for lumi in [137]:
##    toys_SM = make_toys( lumi*sigma_tot(theta_SM), theta_SM, n_toys )
##
##    UL[lumi] = {}
##    for test_statistic in ["total", "lin", "quad"]:
##
##        #wrapper_(0)
##
##        class wrapper_(object):
##            def __init__(self, *args, **kwargs ):
##                self.args = args
##                self.kwargs = kwargs
##            def __call__(self, i_toy):
##                return wrapper( i_toy, *self.args, **self.kwargs )
##
##        p  = Pool(20)
##
##        UL_ = np.array(p.map( wrapper_(verbose=verbose, toys=toys_SM, extended=extended, plot=False, test_statistic=test_statistic), range(n_toys_UL)))
##        UL[lumi][test_statistic] = UL_
##
##        print "lumi %5.1f extended %i test_statistic %10s median expected %i%% CL UL: %5.4f coverage %3.2f%%" % ( lumi, extended, test_statistic, 100*CL, np.quantile( UL_[:,0], 0.5 ), 100*np.count_nonzero(UL_[:,0]>tolerance)/float(len(UL_)) ) 
