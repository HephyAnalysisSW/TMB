import ROOT
import os
import argparse
import array
from   RootTools.core.Sample import Sample
from   math import log
from   TMB.Tools.delphesCutInterpreter  import cutInterpreter
from   Analysis.Tools                   import u_float
import Analysis.Tools.syncer            as syncer 
from   RootTools.core.standard          import *
import TMB.Tools.helpers                as helpers
import TMB.Tools.stat_helpers           as stat_helpers

argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',       action='store',      default='INFO',         nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'],             help="Log level for logging")
argParser.add_argument("--lumi",               action='store',      type=float,             default=137, help='Which lumi?')
argParser.add_argument('--config',             action='store', type=str, default = "ZH_delphes_bkgs", help="config")
argParser.add_argument('--config_module',      action='store', type=str, default = "TMB.BIT.configs", help = "config directory")
argParser.add_argument('--output_directory',   action='store', type=str,   default=os.path.expandvars('/mnt/hephy/cms/$USER/BIT/'))
argParser.add_argument('--plot_directory',     action='store', type=str,   default="BIT_VH_9")
argParser.add_argument('--flavor',     action='store', type=str,   default="nom")
argParser.add_argument('--small',              action='store_true', help="small?")
argParser.add_argument('--name',               action='store', type=str,   default='v2', help="Name of the training")
argParser.add_argument('--input_directory',    action='store', type=str,   default=os.path.expandvars("/groups/hephy/cms/$USER/BIT/training-ntuple-ZH/MVA-training"))
argParser.add_argument("--nBinsTestStat",      action="store",      type=int, default=20,                 help="Number of bins for the test statistic")

args = argParser.parse_args()

# Logging
import Analysis.Tools.logger as logger
logger = logger.get_logger(args.logLevel, logFile = None )
import RootTools.core.logger as logger_rt
logger_rt = logger_rt.get_logger(args.logLevel, logFile = None )

# MVA configuration
import importlib
configs = importlib.import_module(args.config_module)
config  = getattr( configs, args.config)

if args.small:
    args.name+='_small'

import TMB.Tools.user as user
# directories
plot_directory      = os.path.join( user. plot_directory, args.plot_directory, args.config )

import uproot
import awkward
import numpy as np
import pandas as pd

if args.small:
    args.name+='_small'

# get the training variable names
mva_variables = [ mva_variable[0] for mva_variable in config.mva_variables]

n_var_flat   = len(mva_variables)

features           = {}
weight_derivatives = {}

bit_branches = [ 'BIT_%s_%s'%( args.flavor, "_".join(der) ) for der in config.bit_derivatives]
n_der        = len(config.weight_derivative_combinations)
n_bits       = len(config.bit_derivatives)

# only signal for now
training_sample = config.training_samples[0]

n_small_samples = 10000 if args.small else None
bit          = {}
#for i_training_sample, training_sample in enumerate(config.training_samples):
upfile_name = os.path.join(os.path.expandvars(args.input_directory), args.config, training_sample.name, training_sample.name+'.root')
logger.info( "Loading upfile: %s from %s", training_sample.name, upfile_name)
upfile = uproot.open(upfile_name)

features            = upfile["Events"].pandas.df(branches = mva_variables ).values[:n_small_samples]
weight_derivatives  = upfile["Events"].pandas.df(branches = ["weight_derivatives"] ).values.reshape((-1,n_der))[:n_small_samples]
bit                 = upfile["Events"].pandas.df(branches = bit_branches).values.reshape((-1,n_bits))[:n_small_samples]

##############################################
# Histogram of test statistic (total sample) #
##############################################
colors     = [ROOT.kBlack, ROOT.kBlue, ROOT.kGreen, ROOT.kMagenta, ROOT.kCyan, ROOT.kRed]
for WC, WC_vals in [
            ("cHW",    [-.5, -.2, -.1, .1, .2, .5]),
            ("cHWtil", [-.5, -.2, -.1, .1, .2, .5]),
            #("cHQ3",   [-.05, -.02, -.01, .01, .02, .05]),
            ("cHj3",   [-0.12, -.08, -.02, .01, .02, .05]),
        ]:

    # find position of linear and quadratic estimated coefficient 
    i_lin       = config.bit_derivatives.index((WC,))
    i_quad      = config.bit_derivatives.index((WC, WC))
    # find linear and quadratic weights
    i_der_lin   = config.weight_derivative_combinations.index((WC,))
    i_der_quad  = config.weight_derivative_combinations.index((WC, WC))

    for test_statistic in [ "lin", "quad", "total"] :

        histos = []
        for i_WC_val, WC_val in enumerate(WC_vals):
            # calculate test statistic
            if test_statistic == "lin":
                q_event = bit[:,i_lin]
            elif test_statistic == "quad":
                q_event =  bit[:,i_quad]
            elif test_statistic == "total":
                q_event =  bit[:,i_lin] + 0.5*WC_val*bit[:,i_quad]
            else:
                raise RuntimeError

            # compute BSM weights
            w_sm   = args.lumi/float(config.scale_weight)*(weight_derivatives[:,0]) 
            w_bsm  = args.lumi/float(config.scale_weight)*(weight_derivatives[:,0] + WC_val*weight_derivatives[:,i_der_lin]+0.5*WC_val**2*weight_derivatives[:,i_der_quad])

            q_event_argsort     = np.argsort(q_event)
            q_event_argsort_inv = np.argsort(q_event_argsort)
            cdf_sm = np.cumsum(w_sm[q_event_argsort])
            cdf_sm/=cdf_sm[-1]

            q_event_cdf = cdf_sm[q_event_argsort_inv] #uniformly distributed under the SM hypothesis

            #min_, max_ = 0, 1 
            binning = np.linspace(0, 1, args.nBinsTestStat+1)

            np_histo_sm  = np.histogram(q_event_cdf, bins=binning, weights = w_sm ) 
            np_histo_bsm = np.histogram(q_event_cdf, bins=binning, weights = w_bsm )
            histo_sm     = stat_helpers.make_TH1F(np_histo_sm)
            histo_bsm    = stat_helpers.make_TH1F(np_histo_bsm)

            histo_sm.legendText  = "#color[%i]{p(q_{%s=%3.2f}|SM)}" % ( colors[i_WC_val], WC, WC_val)
            histo_bsm.legendText = "#color[%i]{p(q_{%s=%3.2f}|%s=%3.2f)}" % ( colors[i_WC_val], WC, WC_val, WC, WC_val )
            histo_sm.style       = styles.lineStyle( colors[i_WC_val], dashed = True)
            histo_bsm.style      = styles.lineStyle( colors[i_WC_val] )

            histos.append( histo_sm )
            histos.append( histo_bsm )

        drawObjects = [ ]

        plot = Plot.fromHisto( "test_stat_%s_%s_lumi_%i_nBinsTestStat_%i"%(test_statistic, WC, args.lumi, args.nBinsTestStat), [[h] for h in histos], texX = "q", texY = "Entries" )
        plotting.draw( plot,
            plot_directory = os.path.join( plot_directory, "binned" ),
            #ratio          = {'yRange':(0.6,1.4)} if len(plot.stack)>=2 else None,
            logX = False, sorting = False,
            legend         = ( (0.15,0.7,0.9,0.92),2),
            drawObjects    = drawObjects,
            copyIndexPHP   = True,
            extensions     = ["png"],
          )

#for q_name, q in qs.iteritems():
#    plot = Plot.fromHisto(name = "q_%s_%s"%(q_name,args.WC), 
#            histos =  [ [h] for h in q['histos'] ],
#            texX = "q_{%s} %s"%(args.WC, q_name) , texY = "Number of Events" )
#
#    for log_ in [True, False]:
#        plot_directory_ = os.path.join(plot_directory, ("log" if log_ else "lin"))
#        plotting.draw(plot, plot_directory = plot_directory_, ratio = {'histos':[(1,0)], 'texY': 'Ratio'}, logY = log_, logX = False, yRange = (10**-2,"auto"), 
#            #yRange = (0, 0.5) if 'Mult' not in var else (0, 15 ),  
#            legend = ([0.20,0.75,0.9,0.88],3), copyIndexPHP=True, )

#######################################
## Histogram of test statistic (toys) #
#######################################
#
#colors   = [ROOT.kBlack, ROOT.kBlue, ROOT.kGreen, ROOT.kMagenta, ROOT.kCyan, ROOT.kRed]
#extended = True
#lumi     = 137
#n_toys   = 5000
#theta_SM = 0
#theta_vals = [0.05, 0.04, 0.03, 0.02, 0.01]
#
##for lumi in [ 13.70/5., 13.70/2, 13.70 , 2*13.70, 5*13.70, 10*13.70 ]: 
#for lumi in [ 1370 ]: 
#    q_theta_given_theta = {}
#    q_theta_given_SM    = {}
#    for test_statistic in ["lin", "quad", "total"]: 
#        q_theta_given_theta [test_statistic]= {}
#        q_theta_given_SM    [test_statistic]= {}
#
#        print "Test statistic", test_statistic
#
#        histos = []
#        quantile_lines  = []
#        for i_theta, theta in enumerate(theta_vals):
#            print "theta", theta
#            if test_statistic == "quad":
#                q_event = np.nan_to_num( np.log( 1 +  0.5*bit[:,i_quad] ) )
#            elif test_statistic == "lin":
#                q_event = 1./theta * np.nan_to_num( np.log( (1 + theta*bit[:,i_lin])**2 ))  if theta>0.001 else 2*bit[:,i_lin]
#            elif test_statistic == "total":
#                q_event = 1./theta * np.nan_to_num(np.log( 1 + theta*bit[:,i_lin] + 0.5*theta**2*bit[:,i_quad]))  if theta>0.001 else bit[:,i_lin] + 0.5*theta*(bit[:,i_quad]-bit[:,i_lin]**2) 
#            else:
#                raise RuntimeError( "Unknwon test statistc %s" % test_statistic )
#
#            log_sigma_tot_ratio_subtraction = np.log(sigma_tot_ratio(theta)) if not extended else 0
#            q_theta_given_theta[test_statistic][theta] = np.array([np.sum( q_event[toy_] - log_sigma_tot_ratio_subtraction ) for toy_ in make_toys( lumi*sigma_tot(theta),    theta,    n_toys ) ])
#            q_theta_given_SM   [test_statistic][theta] = np.array([np.sum( q_event[toy_] - log_sigma_tot_ratio_subtraction ) for toy_ in make_toys( lumi*sigma_tot(theta_SM), theta_SM, n_toys ) ])
#
#        for i_theta, theta in enumerate(theta_vals):
#
#            if i_theta == 0:
#                all_vals = sum( [list(q_theta_given_theta[test_statistic][k])+list(q_theta_given_SM[test_statistic][k]) for k in q_theta_given_theta[test_statistic].keys()], [])
#                min_, max_ = min( all_vals ), max( all_vals )
#                binning = np.arange(min_, max_, (max_-min_)/100.)
#
#            #np_histo_all   = np.histogram(q_theta_given_theta+q_theta_given_SM, 100)
#            #histo_all      = stat_helpers.make_TH1F(np_histo_all)
#            #binning        = np_histo_all[1] 
#
#            np_histo_SM    = np.histogram(q_theta_given_SM   [test_statistic][theta],    bins=binning)
#            np_histo_theta = np.histogram(q_theta_given_theta[test_statistic][theta], bins=binning)
#            histo_SM       = stat_helpers.make_TH1F(np_histo_SM)
#            histo_theta    = stat_helpers.make_TH1F(np_histo_theta)
#
#            histo_SM.legendText    = "#color[%i]{p(q_{#theta}|#theta)}, #theta =%3.2f" % ( colors[i_theta], theta_SM )
#            histo_theta.legendText = "#color[%i]{p(q_{#theta}|#theta_{SM})}, #theta = %3.2f" % ( colors[i_theta], theta )
#            histo_SM.style         = styles.lineStyle( colors[i_theta], dashed = True)
#            histo_theta.style      = styles.lineStyle( colors[i_theta] ) 
#            #histo_all = histo_all
#            #histo_all.legendText = None
#            #histo_all.style      = styles.invisibleStyle()
#
#            #histos.append( histo_all )
#            histos.append( histo_SM )
#            histos.append( histo_theta )
#
#            for x in np.quantile( q_theta_given_SM[test_statistic][theta], [.05, .95 ] ):
#                quantile_lines.append( ROOT.TLine(x, 0, x, histo_SM.GetBinContent(histo_SM.FindBin(x))) )
#                quantile_lines[-1].SetLineColor( colors[i_theta] )
#                quantile_lines[-1].SetLineStyle( 7 )
#            for x in np.quantile( q_theta_given_theta[test_statistic][theta], [.05, .95 ] ):
#                quantile_lines.append( ROOT.TLine(x, 0, x, histo_theta.GetBinContent(histo_theta.FindBin(x))) )
#                quantile_lines[-1].SetLineColor( colors[i_theta] )
#
#            # Text on the plots
#            lines = [ 
#                    #  (0.25, 0.88, "#color[4]{%i%% qu. q_{BSM} = %3.2f}" % ( 100*(1-CL), q_theta_given_theta_1mCL[theta]) ),
#                    #  (0.25, 0.83, "#color[2]{q_{SM} = %3.2f}" % ( q_theta_SM ) ),
#                    #  (0.25, 0.78, "#theta_{current} = %5.4f" % theta ),
#                    ]
#            drawObjects = [ tex.DrawLatex(*line) for line in lines ]
#
#        plot = Plot.fromHisto( "test_stat_%s_nEvents_%3.2f"%(test_statistic, lumi*sigma_tot(theta_SM)), [[h] for h in histos], texX = "q", texY = "Entries" )
#        plotting.draw( plot,
#            plot_directory = os.path.join( plot_directory, ),
#            #ratio          = {'yRange':(0.6,1.4)} if len(plot.stack)>=2 else None,
#            logX = False, sorting = False,
#            legend         = ( (0.15,0.7,0.9,0.92),2),
#            drawObjects    =  quantile_lines + drawObjects,
#            copyIndexPHP   = True,
#            extensions     = ["png"], 
#          )            


################################################
## Plot quantiles of test statistics (TGraphs) #
################################################
#
#colors   = [ ROOT.kRed, ROOT.kBlue, ROOT.kBlack, ROOT.kBlue, ROOT.kRed]
#levels   = [ .05,       .32,        .5,          .68,        .95      ]
#extended = True
#n_toys   = 50000
#
#theta_SM = 0
#
#theta_min= -.7
#theta_max= .5
#Nbins    = 50
#UL = {}
#for lumi in 137*np.array([ 1/50., 1/20. , 1/10., 1/5., 1/2., 1., 2., 5. ]):
##for lumi in 137*np.array([ 1/10. ]):
#    UL[lumi] = {level:{} for level in levels if level<0.5}
#    for test_statistic in ["total", "lin", "quad"]:
#        print "Test statistic", test_statistic
#
#        theta_vals      = np.arange( theta_min, theta_max, (theta_max-theta_min)/Nbins)
#        tgraphs_theta   = { level: ROOT.TGraph(len(theta_vals)) for level in levels }
#        tgraphs_SM      = { level: ROOT.TGraph(len(theta_vals)) for level in levels }
#
#        min_, max_ = float('inf'), -float('inf')
#
#        sm_toys = make_toys( lumi*sigma_tot(theta_SM), theta_SM, n_toys )
#
#        for i_theta, theta in enumerate(theta_vals):
#            print "theta", theta
#            if test_statistic == "quad":
#                q_event = np.nan_to_num( np.log( 1 +  0.5*bit[:,i_quad] ) )
#            elif test_statistic == "lin":
#                q_event = 1./theta * np.nan_to_num( np.log( (1 + theta*bit[:,i_lin])**2 ))  if abs(theta)>0.001 else 2*bit[:,i_lin]
#            elif test_statistic == "total":
#                q_event = 1./theta * np.nan_to_num(np.log( 1 + theta*bit[:,i_lin] + 0.5*theta**2*bit[:,i_quad]))  if theta>0.001 else bit[:,i_lin] + 0.5*theta*(bit[:,i_quad]-bit[:,i_lin]**2) 
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
#        l1 = ROOT.TLegend(0.2, 0.8, 0.7, 0.93)
#        l1.SetFillStyle(0)
#        l1.SetShadowColor(ROOT.kWhite)
#        l1.SetBorderSize(0)
#
#        l2 = ROOT.TLegend(0.55, 0.8, 0.8, 0.93)
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
#            tgraphs_SM[level].SetMarkerStyle(0)
#            tgraphs_SM[level].SetMarkerColor(colors[i_level])
#            tgraphs_SM[level].Draw("L")
#
#            l2.AddEntry( tgraphs_theta[level], "Q(%i%%) BSM"%(100*level) )
#            tgraphs_theta[level].SetTitle("")
#            tgraphs_theta[level].GetXaxis().SetTitle("q")
#            tgraphs_theta[level].GetYaxis().SetTitle("#theta")
#            tgraphs_theta[level].SetLineColor( colors[i_level] )
#            tgraphs_theta[level].SetMarkerStyle(0)
#            tgraphs_theta[level].SetMarkerColor(colors[i_level])
#            tgraphs_theta[level].Draw("L")
#
#        l1.Draw()
#        l2.Draw()
#
#        c1.RedrawAxis()
#        c1.Update()
#        c1.Print(os.path.join( plot_directory, "test_stat_quantiles_%s_extended_%i_nEvents_%3.2f.png"% ( test_statistic, extended, lumi*sigma_tot(theta_SM) )))


syncer.sync()

