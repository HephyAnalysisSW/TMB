#!/usr/bin/env python

import ROOT, os, sys
from ROOT import TLegend
from RootTools.core.standard import *
import Analysis.Tools.syncer as syncer
from TMB.Tools.user import plot_directory

import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--config_module',      action='store', type=str, default = "TMB.BIT.configs", help = "config directory")
argParser.add_argument('--config',             action='store', type=str, default = "ZH_delphes", help="config")
argParser.add_argument('--name',               action='store', type=str,   default='default', help="Name of the training")
argParser.add_argument('--variable_set',       action='store', type=str,   default='mva_variables', help="List of variables for training")
argParser.add_argument('--output_directory',   action='store', type=str,   default=os.path.expandvars('/mnt/hephy/cms/$USER/BIT/'))
argParser.add_argument('--input_directory',    action='store', type=str,   default=os.path.expandvars("/scratch-cbe/users/$USER/BIT/training-ntuples-TTZ-flavor/MVA-training"))
argParser.add_argument('--small',              action='store_true', help="small?")
argParser.add_argument('--debug',              action='store_true', help="Make debug plots?")
#argParser.add_argument('--calibrated',         action='store_true', help="Calibrate output?")
argParser.add_argument('--maxEvents',          action='store', default = None, type=int, help="Maximum number of training events")
argParser.add_argument('--bagging_fraction',   action='store', default = 1., type=float, help="Bagging fraction")
argParser.add_argument('--overwrite',          action='store_true', help="Overwrite output?")
argParser.add_argument('--derivative',         action='store', nargs = '*', default = [], help="What to train?")
argParser.add_argument('--lumi_norm',          action='store_true', help="Normalize the events according to lumiweight1fb?")
#argParser.add_argument('--max_local_score',    action='store', type=float, default = None, help="Maximum local score")
#argParser.add_argument('--rel_max_local_score',action='store', type=float, default = None, help="Relative maximum local score - share of highest scores capped")
#argParser.add_argument('--rel_max_local_score_test',action='store', type=float, default = None, help="Relative maximum local score of test Data - share of highest scores capped")

argParser.add_argument('--name_dir',           action='store', type=str,   default='default', help="Name of Plot Directory")
argParser.add_argument('--deri',               action='store', type=str,   default='cHW', help="Name of derivative for LLR")
argParser.add_argument('--bin_number',         action='store', type=int, default = 20, help ="Number of Bins")
argParser.add_argument('--binning_var',        action='store', type=str, default ='mva_Z_pt', help="Binning Variable")
argParser.add_argument('--min_bin',            action='store', type=float, default=0., help="Minimum Bin")
argParser.add_argument('--max_bin',            action='store', type=float, default=400., help="Maximum Bin")

#args = argParser.parse_args()
args, extra = argParser.parse_known_args(sys.argv[1:])

def parse_value( s ):
    try:
        r = int( s )
    except ValueError:
        r = float(s)
    return r
        
extra_args = {key.lstrip('-'):parse_value(value) for key, value in zip(extra[::2], extra[1::2])}

#Logger
import Analysis.Tools.logger as logger
logger = logger.get_logger("INFO", logFile = None )

# MVA configuration
import importlib
configs = importlib.import_module(args.config_module)
config  = getattr( configs, args.config)
config.bit_cfg.update( extra_args )

#if args.calibrated:
#    args.name+="_calibrated"
#    config.bit_cfg['calibrated'] = True
#
import uproot
import awkward
import numpy as np
import pandas as pd
#import h5py
#########################################################################################
# variable definitions
######Be aware - for other plots (histos of BIT training - plot_directory different, since next three lines outcommented
#import Analysis.Tools.user as user
# directories
#plot_directory   = os.path.join( user. plot_directory, 'MVA', args.config, args.name)
output_directory = os.path.join( args.output_directory, 'models', args.config, args.name)
# saving
if not os.path.exists(output_directory):
    os.makedirs(output_directory)

# fix random seed for reproducibility
np.random.seed(1)

# get the training variable names
mva_variables = [ mva_variable[0] for mva_variable in getattr(config, args.variable_set) ]
print "##########################################"
print mva_variables
print "##########################################"
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

features = pd.concat([features[training_sample.name] for training_sample in config.training_samples])
features = features.values

weight_derivatives = pd.concat([weight_derivatives[training_sample.name] for training_sample in config.training_samples])
if args.lumi_norm:
    lumi_weights = pd.concat([lumi_weights[training_sample.name] for training_sample in config.training_samples])

weight_derivatives = weight_derivatives.values.reshape((len(features),-1))[:args.maxEvents]
features           = features[:args.maxEvents]
if args.lumi_norm:
    lumi_weights       = lumi_weights.values[:,0][:args.maxEvents]
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
        #Maximum Local Score
        if args.max_local_score is not None and args.rel_max_local_score is not None:
                raise NameError('max_local_score and rel_max_local_score defined')
        if args.max_local_score is not None:
                mask = np.divide(bits[derivative].training_diff_weights, bits[derivative].training_weights,out = np.zeros_like(bits[derivative].training_diff_weights),where=bits[derivative].training_weights!=0) > args.max_local_score
                bits[derivative].training_diff_weights[mask] = args.max_local_score * bits[derivative].training_weights[mask]
        if args.rel_max_local_score is not None:
                score_store = np.divide(bits[derivative].training_diff_weights, bits[derivative].training_weights,out = np.zeros_like(bits[derivative].training_diff_weights),where=bits[derivative].training_weights!=0) 
                threshold = np.quantile(score_store,1.0-args.rel_max_local_score)
                bits[derivative].training_diff_weights[score_store > threshold] = threshold * bits[derivative].training_weights[score_store > threshold]

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
print "#########################################################TEST#######################################################################"
test_weights_capped = np.copy(test_weights)

#for derivative in config.bit_derivatives:
        #if args.rel_max_local_score_test is not None:
        #        score_store = np.divide(test_weights_capped[:, config.weight_derivative_combinations.index(derivative)], test_weights_capped[:,0],out = np.zeros_like(test_weights_capped[:, config.weight_derivative_combinations.index(derivative)]),where=test_weights_capped[:,0]!=0) 
        #        threshold = np.quantile(score_store,1.0-args.rel_max_local_score_test)
        #        test_weights_capped[:, config.weight_derivative_combinations.index(derivative)][score_store > threshold] = threshold * test_weights_capped[:,0][score_store > threshold]


bsm_colors = [ROOT.kBlack, ROOT.kBlue, ROOT.kGreen, ROOT.kMagenta, ROOT.kCyan]
bit_predictions  = { key:bits[key].vectorized_predict(test_features) for key in  config.bit_derivatives if key!=tuple() }

directory = os.path.join(plot_directory,'VH',args.name_dir,args.deri,args.binning_var)
if not os.path.exists(directory):
        try:
            os.makedirs(directory)
        except IOError:
            pass
directory_log = os.path.join(directory,"log")
if not os.path.exists(directory_log):
        try:
            os.makedirs(directory_log)
        except IOError:
            pass
directory_lin = os.path.join(directory,"lin")
if not os.path.exists(directory_lin):
        try:
            os.makedirs(directory_lin)
        except IOError:
            pass

c1 = ROOT.TCanvas()

deri = args.deri
h_cdf_sm = ROOT.TH1D("h_cdf","Integral Transform of Distribution of Test Statistic SM",20,0.0,1.0)
h_cdf_bsm_05_neg = ROOT.TH1D("h_cdf","Integral Transform of Distribution of Test Statistic BSM -0.5",20,0.0,1.0)
h_cdf_bsm_02_neg = ROOT.TH1D("h_cdf","Integral Transform of Distribution of Test Statistic BSM -0.2",20,0.0,1.0)
h_cdf_bsm_01_neg = ROOT.TH1D("h_cdf","Integral Transform of Distribution of Test Statistic BSM -0.1",20,0.0,1.0)
h_cdf_bsm_01 = ROOT.TH1D("h_cdf","Integral Transform of Distribution of Test Statistic BSM 0.1",20,0.0,1.0)
h_cdf_bsm_02 = ROOT.TH1D("h_cdf","Integral Transform of Distribution of Test Statistic BSM 0.2",20,0.0,1.0)
h_cdf_bsm_05 = ROOT.TH1D("h_cdf","Integral Transform of Distribution of Test Statistic BSM 0.5",20,0.0,1.0)
#h_cdf_bsm_10 = ROOT.TH1D("h_cdf","Integral Transform of Distribution of Test Statistic BSM 1.0",20,0.0,1.0)
sum_sm = np.sum(test_weights[:,0]) #need later for normalization of BSM
normalization = sum_sm/20/10**5*137 #20 bins, 10**5 scaling factor in delphes config, 137 luminosity
print "Sum SM: ", sum_sm
weights = test_weights[:,0]
for i in range(20):
    h_cdf_sm.Fill((i+0.9)/20,normalization) #normalization to number of events (sum of weights)
for i_plot_theta, plot_theta in enumerate(np.arange(-1.0,1.01,0.05)):
    h_cdf_bsm = ROOT.TH1D("h_cdf","Integral Transform of Distribution of Test Statistic BSM",20,0.0,1.0)
    q_bsm = 1.0 + plot_theta*bit_predictions[(deri,)] + 0.5*plot_theta**2 *bit_predictions[(deri,deri)]
    weights_sorted = weights[np.argsort(q_bsm)] #weights sorted by test statistic 
    q_bsm_sorted = np.sort(q_bsm)
    weights_cumsum = np.cumsum(weights_sorted)/sum_sm #normalized cumulative sum
    for bin_val in np.arange(0.05,1.01,0.05):
        q_vals = q_bsm_sorted[(weights_cumsum < bin_val)*(weights_cumsum > bin_val-0.05)]
        weight_vals = weights_sorted[(weights_cumsum < bin_val)*(weights_cumsum > bin_val-0.05)]
        #q_vals = q_bsm_sorted[weights_cumsum < bin_val]
        #weight_vals = weights_sorted[weights_cumsum < bin_val]
        qs = q_vals*weight_vals/10**5*137
        for j in range(len(qs)):
            h_cdf_bsm.Fill(bin_val-0.015,qs[j])
            if plot_theta < -0.49 and plot_theta > -0.51:
                h_cdf_bsm_05_neg.Fill(bin_val-0.015,qs[j])
            elif plot_theta <-0.19 and plot_theta > -0.21:
                h_cdf_bsm_02_neg.Fill(bin_val-0.015,qs[j])
            elif plot_theta <-0.09 and plot_theta > -0.11:
                h_cdf_bsm_01_neg.Fill(bin_val-0.015,qs[j])
            elif plot_theta <0.11 and plot_theta > 0.09:
                h_cdf_bsm_01.Fill(bin_val-0.015,qs[j])
            elif plot_theta <0.21 and plot_theta > 0.19:
                h_cdf_bsm_02.Fill(bin_val-0.015,qs[j])
            elif plot_theta < 0.51 and plot_theta > 0.49:
                h_cdf_bsm_05.Fill(bin_val-0.015,qs[j])
            #elif plot_theta > 0.99:
            #    h_cdf_bsm_10.Fill(bin_val-0.015,qs[j])
    filename = "Test_Statistic_%s"% ('_'.join([str(round(plot_theta,2))])) 
    h_cdf_sm.style = styles.lineStyle(ROOT.kRed)
    h_cdf_bsm.style = styles.lineStyle(ROOT.kBlue)
    h_cdf_sm.legendText = "SM"
    h_cdf_bsm.legendText = "BSM_%s"% ('_'.join([deri]))
    histos = [[h_cdf_sm], [h_cdf_bsm]]
    plot   = Plot.fromHisto( filename,  histos, texX = "F(q(#theta)) (cdf)", texY = "q(#theta)*Weights" )
    plotting.draw(plot,
            plot_directory = directory_log,
            yRange = 'auto',
            copyIndexPHP = True
            #ratio = {'yRange':(0.5,1.5)}, logY = False, logX = False,
            )
    plotting.draw(plot,
            plot_directory = directory_lin,
            yRange = 'auto', logY = False,
            copyIndexPHP = True
            #ratio = {'yRange':(0.5,1.5)}, logY = False, logX = False,
            )

filename = "Test_Statistic" 
h_cdf_sm.style = styles.lineStyle(ROOT.kOrange)
h_cdf_bsm_05_neg.style = styles.lineStyle(ROOT.kBlack)
h_cdf_bsm_02_neg.style = styles.lineStyle(ROOT.kBlue)
h_cdf_bsm_01_neg.style = styles.lineStyle(ROOT.kGreen)
h_cdf_bsm_01.style = styles.lineStyle(ROOT.kViolet)
h_cdf_bsm_02.style = styles.lineStyle(ROOT.kCyan)
h_cdf_bsm_05.style = styles.lineStyle(ROOT.kRed)
#h_cdf_bsm_10.style = styles.lineStyle(ROOT.kViolet)
h_cdf_sm.legendText = "SM"
h_cdf_bsm_05_neg.legendText = "BSM -0.5_%s"% ('_'.join([deri]))
h_cdf_bsm_02_neg.legendText = "BSM -0.2_%s"% ('_'.join([deri]))
h_cdf_bsm_01_neg.legendText = "BSM -0.1_%s"% ('_'.join([deri]))
h_cdf_bsm_01.legendText = "BSM 0.1_%s"% ('_'.join([deri]))
h_cdf_bsm_02.legendText = "BSM 0.2_%s"% ('_'.join([deri]))
h_cdf_bsm_05.legendText = "BSM 0.5_%s"% ('_'.join([deri]))
#h_cdf_bsm_10.legendText = "BSM 1.0_%s"% ('_'.join([deri]))
histos = [[h_cdf_sm], [h_cdf_bsm_05_neg], [h_cdf_bsm_02_neg], [h_cdf_bsm_01_neg], [h_cdf_bsm_01], [h_cdf_bsm_02], [h_cdf_bsm_05]]
plot   = Plot.fromHisto( filename,  histos, texX = "F(q(#theta)) (cdf)", texY = "q(#theta)*Weights" )
plotting.draw(plot,
        plot_directory = directory_log,
        yRange = 'auto',
        copyIndexPHP = True
        #ratio = {'yRange':(0.5,1.5)}, logY = False, logX = False,
        )
plotting.draw(plot,
        plot_directory = directory_lin,
        yRange = 'auto', logY = False,
        copyIndexPHP = True
        #ratio = {'yRange':(0.5,1.5)}, logY = False, logX = False,
        )



print "Finished"
