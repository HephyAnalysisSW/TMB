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

directory = os.path.join(plot_directory,'VH',args.name_dir,args.deri,args.binning_var,"Bins_%s"% ('_'.join([str(args.bin_number)])))
if not os.path.exists(directory):
        try:
            os.makedirs(directory)
        except IOError:
            pass



c1 = ROOT.TCanvas()
#def compute_weights(weights, theta, theta_ref):
#        return weights[:,0] \
#                        + np.array( [ (theta[i]-theta_ref[i])*weights[:,config.weight_derivative_combinations.index(('ctZ',))] for i in range(len(theta)) ] ).sum(axis=0)\
#                        + np.array( [ 0.5*(theta[i]-theta_ref[i])*(theta[j]-theta_ref[j])*weights[:,config.weight_derivative_combinations.index(('ctZ','ctZ'))] for i in range(len(theta)) for j in range(len(theta)) ] ).sum(axis=0)

#def predict_weight_ratio(bit_predictions, theta, theta_ref):
#        return 1.+ \
#                        + np.array( [ (theta[i]-theta_ref[i])*bit_predictions[('ctZ',)] for i in range(len(theta)) ] ).sum(axis=0)\
#                        + np.array( [ 0.5*(theta[i]-theta_ref[i])*(theta[j]-theta_ref[j])*bit_predictions[('ctZ','ctZ')] for i in range(len(theta)) for j in range(len(theta)) ] ).sum(axis=0)

print "#######################################################"
deri = args.deri
h_LLR_true = ROOT.TH1D("h_LLR_true","LLR Weights Summed",40,-1.0,1.0)
h_LLR_true_event = ROOT.TH1D("h_LLR_true_event","LLR Eventwise",40,-1.0,1.0)
#h_LLR_true_capped = ROOT.TH1D("h_LLR_true_capped","LLR Weights Summed",40,-1.0,1.0)
#h_LLR_true_capped_event = ROOT.TH1D("h_LLR_true_capped_event","LLR Eventwise",40,-1.0,1.0)
h_LLR_pred_lin = ROOT.TH1D("h_LLR_pred_lin","LLR Weights Summed",40,-1.0,1.0)
h_LLR_pred_lin_event = ROOT.TH1D("h_LLR_pred_lin","LLR Eventwise",40,-1.0,1.0)
h_LLR_pred_quad = ROOT.TH1D("h_LLR_pred_quad","LLR Weights Summed",40,-1.0,1.0)
h_LLR_pred_quad_event = ROOT.TH1D("h_LLR_pred_quad","LLR Eventwise",40,-1.0,1.0)
for i_plot_theta, plot_theta in enumerate(np.arange(-1.0,1.0,0.05)):
        w0_sum = 0
#        w0_capped_sum = 0
        w1_true_sum = 0
#        w1_true_capped_sum = 0
        w1_pred_quad_sum = 0
        w1_pred_lin_sum = 0
        LLR_true_event = 0
#        LLR_true_capped_event = 0
        LLR_pred_lin_event = 0
        LLR_pred_quad_event = 0
        for j in range(test_features[:,17].size):
            w1_true_inter = test_weights[j,0] + plot_theta*test_weights[j,config.weight_derivative_combinations.index((deri,))] + 0.5 * plot_theta**2 * test_weights[j,config.weight_derivative_combinations.index((deri,deri))]
#            w1_true_capped_inter = test_weights_capped[j,0] + plot_theta*test_weights_capped[j,config.weight_derivative_combinations.index(('cHW',))] + 0.5 * plot_theta**2 * test_weights_capped[j,config.weight_derivative_combinations.index(('cHW','cHW'))]
            w0_inter = test_weights[j,0]
#            w0_capped_inter = test_weights_capped[j,0]
            div = np.divide(w1_true_inter,w0_inter,out = np.zeros_like(w1_true_inter),where=w0_inter!=0)
            LLR_true_event += w1_true_inter - w0_inter -w0_inter*np.log(div,out = np.zeros_like(div),where=div>0)
            w1_true_sum += w1_true_inter
#            w1_true_capped_sum += w1_true_capped_inter
#            div = np.divide(w1_true_capped_inter,w0_capped_inter,out = np.zeros_like(w1_true_capped_inter),where=w0_capped_inter!=0)
#            LLR_true_capped_event += w1_true_capped_inter - w0_capped_inter -w0_capped_inter*np.log(div,out = np.zeros_like(div),where=div>0) #changed - set log to zero if encounter negative div value
            w0_sum += w0_inter
#            w0_capped_sum += w0_capped_inter
            w1_pred_quad_inter = test_weights[j,0] + plot_theta * bit_predictions[(deri,)][j]*test_weights[j,0] + 0.5*plot_theta**2 * bit_predictions[(deri,deri)][j]*test_weights[j,0]
            w1_pred_quad_sum += w1_pred_quad_inter
            div = np.divide(w1_pred_quad_inter,w0_inter,out = np.zeros_like(w1_pred_quad_inter),where=w0_inter!=0)
            LLR_pred_quad_event += w1_pred_quad_inter - w0_inter -w0_inter*np.log(div,out = np.zeros_like(div),where=div>0)
            w1_pred_lin_inter = test_weights[j,0] + plot_theta * bit_predictions[(deri,)][j]*test_weights[j,0] 
            w1_pred_lin_sum += w1_pred_lin_inter 
            div = np.divide(w1_pred_lin_inter,w0_inter,out = np.zeros_like(w1_pred_lin_inter),where=w0_inter!=0)
            LLR_pred_lin_event += w1_pred_lin_inter - w0_inter -w0_inter*np.log(div,out = np.zeros_like(div),where=div>0)
        LLR_true = w1_true_sum - w0_sum - w0_sum*np.log(w1_true_sum/w0_sum)
#        LLR_true_capped = w1_true_capped_sum - w0_capped_sum - w0_capped_sum*np.log(w1_true_capped_sum/w0_capped_sum)
        LLR_pred_quad = w1_pred_quad_sum - w0_sum - w0_sum*np.log(w1_pred_quad_sum/w0_sum)
        LLR_pred_lin = w1_pred_lin_sum - w0_sum - w0_sum*np.log(w1_pred_lin_sum/w0_sum)
        print "Theta: ", plot_theta, " LLR True: ", LLR_true, " LLR Pred Quad: ", LLR_pred_quad, " LLR Pred Lin: ", LLR_pred_lin
        h_LLR_true.Fill(plot_theta+0.01,LLR_true)
        h_LLR_true_event.Fill(plot_theta+0.01,LLR_true_event)
#        h_LLR_true_capped.Fill(plot_theta+0.01,LLR_true_capped)
#        h_LLR_true_capped_event.Fill(plot_theta+0.01,LLR_true_capped_event)
        h_LLR_pred_quad_event.Fill(plot_theta+0.01,LLR_pred_quad_event)
        h_LLR_pred_lin_event.Fill(plot_theta+0.01,LLR_pred_lin_event)
        h_LLR_pred_quad.Fill(plot_theta+0.01,LLR_pred_quad)
        h_LLR_pred_lin.Fill(plot_theta+0.01,LLR_pred_lin)

print "########################################################"
h_LLR_true.style = styles.lineStyle(ROOT.kRed)
h_LLR_pred_lin.style = styles.lineStyle(ROOT.kGreen)
h_LLR_pred_quad.style = styles.lineStyle(ROOT.kBlue)
h_LLR_true.legendText = "LLR true"
h_LLR_pred_lin.legendText = "LLR linear prediction"
h_LLR_pred_quad.legendText = "LLR quadratic prediction"
histos = [[h_LLR_true], [h_LLR_pred_quad], [h_LLR_pred_lin]]
plot   = Plot.fromHisto( "LLR",  histos, texX = deri, texY = "LLR" )
plotting.draw(plot,
                plot_directory = directory,
                yRange = 'auto',
                #ratio = {'yRange':(0.5,1.5)}, logY = False, logX = False,
                )

h_LLR_true_event.style = styles.lineStyle(ROOT.kRed)
h_LLR_pred_lin_event.style = styles.lineStyle(ROOT.kGreen)
h_LLR_pred_quad_event.style = styles.lineStyle(ROOT.kBlue)
h_LLR_true_event.legendText = "LLR true"
h_LLR_pred_lin_event.legendText = "LLR linear prediction"
h_LLR_pred_quad_event.legendText = "LLR quadratic prediction"
histos = [[h_LLR_true_event], [h_LLR_pred_quad_event], [h_LLR_pred_lin_event]]
plot   = Plot.fromHisto( "LLR_eventwise",  histos, texX = deri, texY = "LLR" )
plotting.draw(plot,
                plot_directory = directory,
                yRange = 'auto',
                #ratio = {'yRange':(0.5,1.5)}, logY = False, logX = False,
                )

h_LLR_true = ROOT.TH1D("h_LLR_true","LLR True",40,-1.0,1.0)
#h_LLR_true_capped = ROOT.TH1D("h_LLR_true_capped","LLR True Capped",40,-1.0,1.0)
h_LLR_pred_lin = ROOT.TH1D("h_LLR_pred_lin","LLR Linear Prediction",40,-1.0,1.0)
h_LLR_pred_quad = ROOT.TH1D("h_LLR_pred_quad","LLR Quadratic Prediction",40,-1.0,1.0)
print "##########################################################"
print "BINNED:"
bin_number = args.bin_number
binning_var = args.binning_var
min_bin = args.min_bin
max_bin = args.max_bin
bin_index = mva_variables.index(binning_var)
for i_plot_theta, plot_theta in enumerate(np.arange(-1.0,1.0,0.05)):
        h_w1_pred_lin = ROOT.TH1D("h_w1_pred_lin","Linear Prediction",bin_number,min_bin,max_bin)
        h_w1_pred_quad = ROOT.TH1D("h_w1_pred_quad","Quadratic Prediction",bin_number,min_bin,max_bin)
        h_w1_true = ROOT.TH1D("h_w1_true","True",bin_number,min_bin,max_bin)
#        h_w1_true_capped = ROOT.TH1D("h_w1_true_capped","True Capped",bin_number,min_bin,max_bin)
        h_w0 = ROOT.TH1D("h_w0", "Weight",bin_number,min_bin,max_bin)
        w0 = 0
        w1_true = 0
#        w1_true_capped = 0
        w1_pred_quad = 0
        w1_pred_lin = 0
        for j in range(test_features[:,bin_index].size):
            w1_true = test_weights[j,0] + plot_theta*test_weights[j,config.weight_derivative_combinations.index((deri,))] + 0.5 * plot_theta**2 * test_weights[j,config.weight_derivative_combinations.index((deri,deri))]
#            w1_true_capped = test_weights_capped[j,0] + plot_theta*test_weights_capped[j,config.weight_derivative_combinations.index(('ctZ',))] + 0.5 * plot_theta**2 * test_weights_capped[j,config.weight_derivative_combinations.index(('ctZ','ctZ'))]
            w0 = test_weights[j,0]
            w1_pred_quad = test_weights[j,0] + plot_theta * bit_predictions[(deri,)][j]*test_weights[j,0] + 0.5*plot_theta**2 * bit_predictions[(deri,deri)][j]*test_weights[j,0]
            w1_pred_lin = test_weights[j,0] + plot_theta * bit_predictions[(deri,)][j]*test_weights[j,0] 
            h_w1_pred_lin.Fill(test_features[j,bin_index],w1_pred_lin)
            h_w1_pred_quad.Fill(test_features[j,bin_index],w1_pred_quad)
            h_w1_true.Fill(test_features[j,bin_index],w1_true)
#            h_w1_true_capped.Fill(test_features[j,bin_index],w1_true_capped)
            h_w0.Fill(test_features[j,bin_index],w0)
        h_w1_pred_lin_clone = h_w1_pred_lin.Clone("h_w1_pred_lin_clone")
        h_w1_pred_quad_clone = h_w1_pred_quad.Clone("h_w1_pred_quad_clone")
        h_w1_true_clone = h_w1_true.Clone("h_w1_true_clone")
#        h_w1_true_capped_clone = h_w1_true_capped.Clone("h_w1_true_capped_clone")
        h_w0_clone = h_w0.Clone("h_w0_clone")
        h_w0_clone_pred_lin = h_w0.Clone("h_w0_clone_pred_lin")
        h_w0_clone_pred_quad = h_w0.Clone("h_w0_clone_pred_quad")
        h_w0_clone_true = h_w0.Clone("h_w0_clone_true")
#        h_w0_clone_true_capped = h_w0.Clone("h_w0_clone_true_capped")
        h_w1_pred_lin_clone.Divide(h_w0_clone)
        h_w1_pred_quad_clone.Divide(h_w0_clone)
        h_w1_true_clone.Divide(h_w0_clone)
#        h_w1_true_capped_clone.Divide(h_w0_clone)
        for k in range(bin_number):
                inter = h_w1_pred_lin_clone.GetBinContent(k+1)
                h_w1_pred_lin_clone.SetBinContent(k+1,np.log(inter,out = np.zeros_like(inter),where=inter!=0))
                inter = h_w1_pred_quad_clone.GetBinContent(k+1)
                h_w1_pred_quad_clone.SetBinContent(k+1,np.log(inter,out = np.zeros_like(inter),where=inter!=0))
                inter = h_w1_true_clone.GetBinContent(k+1)
                h_w1_true_clone.SetBinContent(k+1,np.log(inter,out = np.zeros_like(inter),where=inter!=0))
#                inter = h_w1_true_capped_clone.GetBinContent(k+1)
#                h_w1_true_capped_clone.SetBinContent(k+1,np.log(inter,out = np.zeros_like(inter),where=inter!=0))
        h_w0_clone_pred_lin.Multiply(h_w1_pred_lin_clone)
        h_w0_clone_pred_quad.Multiply(h_w1_pred_quad_clone)
        h_w0_clone_true.Multiply(h_w1_true_clone)
#        h_w0_clone_true_capped.Multiply(h_w1_true_capped_clone)
        h_w0_clone = h_w0.Clone("h_w0_clone")
        h_w0_clone.Add(h_w0_clone_pred_lin)
        h_w1_pred_lin.Add(h_w0_clone,-1)
        h_w0_clone = h_w0.Clone("h_w0_clone")
        h_w0_clone.Add(h_w0_clone_pred_quad)
        h_w1_pred_quad.Add(h_w0_clone,-1)
        h_w0_clone = h_w0.Clone("h_w0_clone")
        h_w0_clone.Add(h_w0_clone_true)
        h_w1_true.Add(h_w0_clone,-1)
 #       h_w0_clone = h_w0.Clone("h_w0_clone")
 #       h_w0_clone.Add(h_w0_clone_true_capped)
 #       h_w1_true_capped.Add(h_w0_clone,-1)
        LLR_pred_lin = h_w1_pred_lin.Clone("LLR_pred_lin")
        LLR_pred_quad = h_w1_pred_quad.Clone("LLR_pred_quad")
        LLR_true = h_w1_true.Clone("LLR_true")
 #       LLR_true_capped = h_w1_true_capped.Clone("LLR_true_capped")

        true_inter = 0
#        true_capped_inter = 0
        lin_inter = 0
        quad_inter = 0
        for k in range(bin_number):
                true_inter += LLR_true.GetBinContent(k+1)
#                true_capped_inter += LLR_true_capped.GetBinContent(k+1)
                lin_inter += LLR_pred_lin.GetBinContent(k+1)
                quad_inter += LLR_pred_quad.GetBinContent(k+1)
        h_LLR_true.Fill(plot_theta+0.01,true_inter)
#        h_LLR_true_capped.Fill(plot_theta+0.01,true_capped_inter)
        h_LLR_pred_quad.Fill(plot_theta+0.01,quad_inter)
        h_LLR_pred_lin.Fill(plot_theta+0.01,lin_inter)
        if binning_var == "mva_Z_pt":
                axis_label= "p_{T}(Z) (GeV)"
        elif binning_var == "mva_W_pt":
                axis_label = "p_{T}(W) (GeV)"
        else:
                axis_label = binning_var
        filename = "LLR_%s"% ('_'.join([str(round(plot_theta,2))])) 
        LLR_true.style = styles.lineStyle(ROOT.kRed)
        LLR_pred_lin.style = styles.lineStyle(ROOT.kGreen)
        LLR_pred_quad.style = styles.lineStyle(ROOT.kBlue)
        LLR_true.legendText = "LLR true"
        LLR_pred_lin.legendText = "LLR linear prediction"
        LLR_pred_quad.legendText = "LLR quadratic prediction"
        histos = [[LLR_true], [LLR_pred_quad], [LLR_pred_lin]]
        plot   = Plot.fromHisto( filename,  histos, texX = axis_label, texY = "LLR" )
        plotting.draw(plot,
                plot_directory = directory,
                yRange = 'auto',
                #ratio = {'yRange':(0.5,1.5)}, logY = False, logX = False,
                )


h_LLR_true.style = styles.lineStyle(ROOT.kRed)
h_LLR_pred_lin.style = styles.lineStyle(ROOT.kGreen)
h_LLR_pred_quad.style = styles.lineStyle(ROOT.kBlue)
h_LLR_true.legendText = "LLR true"
h_LLR_pred_lin.legendText = "LLR linear prediction"
h_LLR_pred_quad.legendText = "LLR quadratic prediction"
histos = [[h_LLR_true], [h_LLR_pred_quad], [h_LLR_pred_lin]]
plot   = Plot.fromHisto( "LLR_binned",  histos, texX = deri, texY = "LLR" )
plotting.draw(plot,
                plot_directory = directory,
                yRange = 'auto',
                #ratio = {'yRange':(0.5,1.5)}, logY = False, logX = False,
                )
print "Finished"
