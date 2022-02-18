import ROOT
import os
import argparse
import array
import copy
from   RootTools.core.Sample import Sample
from   math import log, sqrt
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
argParser.add_argument('--flavor',             action='store', type=str,   choices = ["nom", "bkgs"], default="bkgs")
argParser.add_argument('--small',              action='store_true', help="small?")
argParser.add_argument('--name',               action='store', type=str,   default='v2', help="Name of the training")
argParser.add_argument('--input_directory',    action='store', type=str,   default=os.path.expandvars("/groups/hephy/cms/$USER/BIT/training-ntuple-ZH/MVA-training"))
argParser.add_argument("--nBins",              action="store",      type=int, default=30,                 help="Number of bins in each dimension")
argParser.add_argument("--nBinsTestStat",      action="store",      type=int, default=20,                 help="Number of bins for the test statistic")
argParser.add_argument("--WCs",                action="store",      nargs='*', default=["cHj3", -.08, .08, "cHW", 0,.5],                 help="Wilson coefficients")

args = argParser.parse_args()
WC1, WC1_min, WC1_max, WC2, WC2_min, WC2_max = map( lambda x:float(x[1]) if x[0] in [1,2,4,5] else x[1], enumerate(args.WCs) )

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

# load signal and background 
features          = {}
weight_derivatives= {}
bit               = {}
for i_sample, sample in enumerate(config.training_samples):
    n_small_samples = 10000 if args.small else None
    #for i_sample, sample in enumerate(config.samples):
    upfile_name = os.path.join(os.path.expandvars(args.input_directory), args.config, sample.name, sample.name+'.root')
    logger.info( "Loading upfile: %s from %s", sample.name, upfile_name)
    upfile = uproot.open(upfile_name)
    key = "sig" if i_sample==0 else "bkg"
    if features.has_key(key):
        features          [key] = np.append(features          [key],  upfile["Events"].pandas.df(branches = mva_variables ).values[:n_small_samples], axis = 0)
        weight_derivatives[key] = np.append(weight_derivatives[key],  upfile["Events"].pandas.df(branches = ["weight_derivatives"] ).values.reshape((-1,n_der))[:n_small_samples], axis = 0)
        bit               [key] = np.append(bit               [key],  upfile["Events"].pandas.df(branches = bit_branches).values.reshape((-1,n_bits))[:n_small_samples], axis = 0)
    else :
        features          [key] =                                     upfile["Events"].pandas.df(branches = mva_variables ).values[:n_small_samples]
        weight_derivatives[key] =                                     upfile["Events"].pandas.df(branches = ["weight_derivatives"] ).values.reshape((-1,n_der))[:n_small_samples]
        bit               [key] =                                     upfile["Events"].pandas.df(branches = bit_branches).values.reshape((-1,n_bits))[:n_small_samples]

def make_cdf_map( x, y ):
    import scipy.interpolate
    map__ = scipy.interpolate.interp1d(x, y, kind='linear')
    max_x, min_x = max(x), min(x)
    max_y, min_y = max(y), min(y)
    def map_( x_ ):
        x__ = np.array(x_)
        result = np.zeros_like(x__).astype('float')
        result[x__>max_x] = max_y
        result[x__<min_x] = min_y
        vals = (x__>=min_x) & (x__<=max_x)
        result[vals] = map__(x__[vals]) 
        return result 

    return map_ 

def make_q( order, bit, **kwargs ):
    eft      = kwargs
    if order not in ["lin", "quad", "total"]:
        raise RuntimeError("Order %s not known" % order )
    index = {}
    for i_comb, comb in enumerate(config.bit_derivatives):
        index[comb] = i_comb
        index[tuple(reversed(comb))] = i_comb
    result = np.zeros(len(bit))
    if order in ["lin", "total"]:
        for coeff in eft.keys():
            result += eft[coeff]*bit[:,index[(coeff,)]]
    if order in ["quad", "total"]:
        for coeff1 in eft.keys():
            for coeff2 in eft.keys():
                result += .5*eft[coeff1]*eft[coeff2]*bit[:, index[(coeff1, coeff2)]]
    return result 

def make_weights( weight_derivatives, **kwargs ):
    eft      = kwargs
    index = {}
    for i_comb, comb in enumerate(config.weight_derivative_combinations):
        index[comb] = i_comb
        index[tuple(reversed(comb))] = i_comb

    result = copy.deepcopy(weight_derivatives[:,0])
    for coeff1 in eft.keys():
        result += eft[coeff1]*weight_derivatives[:,index[(coeff1,)]]
        for coeff2 in eft.keys():
            result += .5*eft[coeff1]*eft[coeff2]*weight_derivatives[:, index[(coeff1, coeff2)]]
    return result

def getContours( h, level):
    _h     = h.Clone()
    ctmp = ROOT.TCanvas()
    _h.SetContour(1,array.array('d', [level]))
    _h.Draw("contzlist")
    _h.GetZaxis().SetRangeUser(0.0001,1)
    ctmp.Update()
    contours = ROOT.gROOT.GetListOfSpecials().FindObject("contours")
    return contours.At(0).Clone()

################################################
## Histogram of test statistic (signal sample) #
################################################
#colors     = [ROOT.kBlack, ROOT.kBlue, ROOT.kGreen, ROOT.kMagenta, ROOT.kCyan, ROOT.kRed]
#for WC, WC_vals in [
#            ("cHj3",   [-0.08, -.05, -.02, .01, .02, .05]),
#            ("cHW",    [-.5, -.2, -.1, .1, .2, .5]),
#            ("cHWtil", [-.5, -.2, -.1, .1, .2, .5]),
#            #("cHQ3",   [-.05, -.02, -.01, .01, .02, .05]),
#        ]:
#
#    for test_statistic in [ "total", "lin", "quad"] :
#        for key in ["sig", "bkg"]:
#            histos = []
#            for i_WC_val, WC_val in enumerate(WC_vals):
#
#                eft = {WC:WC_val}
#                q_event = make_q( test_statistic, bit[key], **eft )
#
#                # compute BSM weights
#                w_sm   = args.lumi/float(config.scale_weight)*make_weights( weight_derivatives[key] )
#                w_bsm  = args.lumi/float(config.scale_weight)*make_weights( weight_derivatives[key], **eft)
#
#                q_event_argsort     = np.argsort(q_event)
#                q_event_argsort_inv = np.argsort(q_event_argsort)
#                cdf_sm = np.cumsum(w_sm[q_event_argsort])
#                cdf_sm/=cdf_sm[-1]
#
#                # map to the SM CDF of q
#                if key=="sig":
#                    cdf_map = make_cdf_map( q_event[q_event_argsort], cdf_sm )
#
#                #q_event_cdf = cdf_sm[q_event_argsort_inv] #uniformly distributed under the SM hypothesis
#                q_event_cdf = cdf_map( q_event )
#
#                #min_, max_ = 0, 1 
#                binning = np.linspace(0, 1, args.nBinsTestStat+1)
#
#                np_histo_sm  = np.histogram(q_event_cdf, bins=binning, weights = w_sm ) 
#                np_histo_bsm = np.histogram(q_event_cdf, bins=binning, weights = w_bsm )
#
#                histo_sm     = stat_helpers.make_TH1F(np_histo_sm)
#                histo_bsm    = stat_helpers.make_TH1F(np_histo_bsm)
#
#                histo_sm.legendText  = "#color[%i]{p(q_{%s=%3.2f}|SM)}" % ( colors[i_WC_val], WC, WC_val)
#                histo_bsm.legendText = "#color[%i]{p(q_{%s=%3.2f}|%s=%3.2f)}" % ( colors[i_WC_val], WC, WC_val, WC, WC_val )
#                histo_sm.style       = styles.lineStyle( colors[i_WC_val], dashed = True)
#                histo_bsm.style      = styles.lineStyle( colors[i_WC_val] )
#
#                histos.append( histo_sm )
#                histos.append( histo_bsm )
#
#            drawObjects = [ ]
#
#            plot = Plot.fromHisto( "test_stat_%s_%s_%s_lumi_%i_nBinsTestStat_%i"%(test_statistic, WC, key, args.lumi, args.nBinsTestStat), [[h] for h in histos], texX = "q", texY = "Entries" )
#            plotting.draw( plot,
#                plot_directory = os.path.join( plot_directory, "binned" ),
#                #ratio          = {'yRange':(0.6,1.4)} if len(plot.stack)>=2 else None,
#                logX = False, sorting = False,
#                legend         = ( (0.15,0.7,0.9,0.92),2),
#                drawObjects    = drawObjects,
#                copyIndexPHP   = True,
#                extensions     = ["png"],
#              )

n_toys = 50000

# do not make the following inconsistent
levels          = [ 0.95, 0.68]
quantile_levels = [ 0.025, 0.16, .5, 1-0.16, 1-0.025 ]

exp_nll_ratio = {}

step1 = (WC1_max-WC1_min)/args.nBins
step2 = (WC2_max-WC2_min)/args.nBins
WC1_vals = np.arange(WC1_min, WC1_max+step1, (WC1_max-WC1_min)/args.nBins)
WC2_vals = np.arange(WC2_min, WC2_max+step2, (WC2_max-WC2_min)/args.nBins)

test_statistics = ["total", "lin", "quad"]
exp_nll_ratio = {}
power         = {}
for test_statistic in test_statistics:
    print "Test statistic", test_statistic
    power[test_statistic] = {level:ROOT.TH2D("power_"+test_statistic, "power_"+test_statistic, len(WC1_vals)-1, array.array('d', WC1_vals), len(WC2_vals)-1, array.array('d', WC2_vals)) for level in levels}

    exp_nll_ratio[test_statistic] = ROOT.TH2D("exp_nll_ratio_"+test_statistic, "exp_nll_ratio_"+test_statistic, len(WC1_vals)-1, array.array('d', WC1_vals), len(WC2_vals)-1, array.array('d', WC2_vals))

    for i_WC1_val, WC1_val in enumerate( WC1_vals ):
    #for i_WC1_val, WC1_val in enumerate( [-.08] ):
        for i_WC2_val, WC2_val in enumerate( WC2_vals ):
        #for i_WC2_val, WC2_val in enumerate( [0.] ):
            if WC1_val==WC2_val==0: continue

            eft = {WC1:WC1_val, WC2:WC2_val}

            for key in ["sig", ]:#"bkg"]:

                q_event = make_q( test_statistic, bit[key], **eft )

                # compute BSM weights
                w_sm   = args.lumi/float(config.scale_weight)*make_weights( weight_derivatives[key] )
                w_bsm  = args.lumi/float(config.scale_weight)*make_weights( weight_derivatives[key], **eft)

                q_event_argsort     = np.argsort(q_event)
                q_event_argsort_inv = np.argsort(q_event_argsort)
                cdf_sm = np.cumsum(w_sm[q_event_argsort])
                cdf_sm/=cdf_sm[-1]

                # map to the SM CDF of q
                if key=="sig":
                    cdf_map = make_cdf_map( q_event[q_event_argsort], cdf_sm )

                #q_event_cdf = cdf_sm[q_event_argsort_inv] #uniformly distributed under the SM hypothesis
                q_event_cdf = cdf_map( q_event )

                #min_, max_ = 0, 1 
                binning = np.linspace(0, 1, args.nBinsTestStat+1)

                np_histo_sm  = np.histogram(q_event_cdf, bins=binning, weights = w_sm )[0] 
                np_histo_bsm = np.histogram(q_event_cdf, bins=binning, weights = w_bsm )[0]

                # Expectation_BSM( -Log( Prod_i( Pois_i( n_i, lambda_i(theta))/Pois_i( n_i, lambda_i(0)) ) ))
                exp_nll_ratio_ = 2*np.sum(np_histo_sm - np_histo_bsm - np_histo_bsm*np.log(np_histo_sm/np_histo_bsm))
                exp_nll_ratio[test_statistic].SetBinContent( exp_nll_ratio[test_statistic].FindBin( WC1_val, WC2_val ), exp_nll_ratio_)

                binned_toys_sm    = np.random.poisson(lam=np_histo_sm, size=(n_toys, len(np_histo_sm)))
                binned_toys_theta = np.random.poisson(lam=np_histo_bsm, size=(n_toys, len(np_histo_bsm)))

                q_theta_given_sm    = [ np.sum( toy_ll ) for toy_ll in 2*(np_histo_sm - np_histo_bsm  - binned_toys_sm*np.log(np_histo_sm/np_histo_bsm))]
                q_theta_given_theta = [ np.sum( toy_ll ) for toy_ll in 2*(np_histo_sm - np_histo_bsm  - binned_toys_theta*np.log(np_histo_sm/np_histo_bsm))]

                if True:
                    n = float(len(q_theta_given_sm))
                    mean_q_theta_given_sm     = np.sum(q_theta_given_sm)/n
                    sigma_q_theta_given_sm    = sqrt( np.sum((q_theta_given_sm - mean_q_theta_given_sm)**2)/(n-1) )
                    q_theta_given_sm    = (q_theta_given_sm - mean_q_theta_given_sm)/sigma_q_theta_given_sm
                    q_theta_given_theta = (q_theta_given_theta - mean_q_theta_given_sm)/sigma_q_theta_given_sm

                print i_WC1_val, WC1_val,  i_WC2_val, WC2_val,  "sqrt(2NLL)", sqrt(exp_nll_ratio_)

               # Exclusion: The null hypothesis is the BSM point, the alternate is the SM.
                quantiles_theta = np.quantile( q_theta_given_theta, quantile_levels )
                quantiles_sm    = np.quantile( q_theta_given_sm, quantile_levels )
                size_           = np.sum(np.histogram( q_theta_given_sm, quantiles_sm)[0])/float(n_toys)
                #power_histo     = np.histogram( q_theta_given_theta, quantiles_sm)
                for i_level, level in enumerate(levels):
                    #if level != 0.68: continue
                    power_toy_count = np.count_nonzero((q_theta_given_sm>=quantiles_theta[i_level]) & (q_theta_given_sm<quantiles_theta[-1-i_level]))
                    power_ = 1. - power_toy_count/float(n_toys)
                    power[test_statistic][level].SetBinContent( power[test_statistic][level].FindBin( WC1_val, WC2_val ), power_ )
                #power_ = 1-np.sum(np.histogram( q_theta_given_theta, quantiles_sm)[0])/float(n_toys)
                    print "theta", round(WC1_val,3), round(WC2_val,3), "level", level, "size", quantile_levels[-1-i_level] - quantile_levels[i_level], "power", round(power_,3), test_statistic, WC1, WC2

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
    plot2D = Plot2D.fromHisto(name = "exp_nll_ratio_%s_%s_vs_%s_lumi_%3.2f_nBinsTestStat_%i"%(test_statistic, WC1, WC2, args.lumi, args.nBinsTestStat), histos = [[exp_nll_ratio[test_statistic]]], texX = WC1, texY = WC2 )
    plotting.draw2D(plot2D, plot_directory = os.path.join( plot_directory, "binned"), logY = False, logX = False, logZ = True, copyIndexPHP=True, drawObjects = contour_objects, zRange = (0.05,25))


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
        plot2D = Plot2D.fromHisto(name = "power_%s_%s_vs_%s_lumi_%3.2f_level_%3.2f_nBinsTestStat_%i"%(test_statistic, WC1, WC2, args.lumi, level, args.nBinsTestStat), histos = [[power[test_statistic][level]]], texX = WC1, texY = WC2 )
        plotting.draw2D(plot2D, plot_directory = os.path.join( plot_directory, "binned"), logY = False, logX = False, logZ = True, copyIndexPHP=True, drawObjects = contour_objects, zRange = (0.005,1))


syncer.sync()

