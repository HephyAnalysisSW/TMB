#!/usr/bin/env python

import ROOT, os, sys
from ROOT import TLegend
from RootTools.core.standard import *
import Analysis.Tools.syncer as syncer
from TMB.Tools.user import plot_directory

import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--config_module',      action='store', type=str, default = "TMB.BIT.configs", help = "config directory")
argParser.add_argument('--config',             action='store', type=str, default = "ttZ_3l_flavor", help="config")
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
argParser.add_argument('--max_local_score',    action='store', type=float, default = None, help="Maximum local score")
argParser.add_argument('--rel_max_local_score',action='store', type=float, default = None, help="Relative maximum local score - share of highest scores capped")
argParser.add_argument('--rel_max_local_score_test',action='store', type=float, default = None, help="Relative maximum local score of test Data - share of highest scores capped")

argParser.add_argument('--name_dir',           action='store', type=str,   default='default', help="Name of Plot Directory")
argParser.add_argument('--bin_number',         action='store', type=int, default = 20, help ="Number of Bins")
argParser.add_argument('--binning_var',        action='store', type=str, default ='mva_photon_pt', help="Binning Variable")
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

#        if args.debug:
#            from debug import make_debug_plots
#            make_debug_plots( bits[derivative], 
#                              training_features, training_weights[:,0],  
#                              training_weights[:, config.weight_derivative_combinations.index(derivative)], 
#                              test_features, 
#                              test_weights[:,0], 
#                              test_weights_, 
#                              os.path.join(plot_directory, ('_'.join(derivative))),
#                              mva_variables = config.mva_variables) 

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

        # Histo style
#        test_FI_histo    .style = styles.lineStyle( ROOT.kBlue, width=2 )
#        training_FI_histo.style = styles.lineStyle( ROOT.kRed, width=2 )
#        test_FI_histo    .legendText = "Test"
#        training_FI_histo.legendText = "Training"
#
#        minY   = 0.01 * min( test_FI_histo.GetBinContent(test_FI_histo.GetMaximumBin()), training_FI_histo.GetBinContent(training_FI_histo.GetMaximumBin()))
#        maxY   = 1.5  * max( test_FI_histo.GetBinContent(test_FI_histo.GetMaximumBin()), training_FI_histo.GetBinContent(training_FI_histo.GetMaximumBin()))

#        histos = [ [test_FI_histo], [training_FI_histo] ]
#        plot   = Plot.fromHisto( "bit_derivative_%s_evolution"% ('_'.join(derivative)), histos, texX="b", texY="L(D,b)" )

        # Plot Style
#        histModifications      = []
#        histModifications      += [ lambda h: h.GetYaxis().SetTitleOffset(1.4) ]
#        histModifications += [ lambda h: h.GetXaxis().SetTitleSize(26) ]
#        histModifications += [ lambda h: h.GetYaxis().SetTitleSize(26) ]
#        histModifications += [ lambda h: h.GetXaxis().SetLabelSize(22)  ]
#        histModifications += [ lambda h: h.GetYaxis().SetLabelSize(22)  ]

#        ratioHistModifications = []
#        ratio                  = None
#        legend                 = (0.6,0.75,0.9,0.88)
#        yRange                 = "auto" #( minY, maxY )
#        plot1DHist( plot, plot_directory, yRange=yRange, ratio=ratio, legend=legend, histModifications=histModifications )
test_weights_capped = np.copy(test_weights)

for derivative in config.bit_derivatives:
        if args.rel_max_local_score_test is not None:
                score_store = np.divide(test_weights_capped[:, config.weight_derivative_combinations.index(derivative)], test_weights_capped[:,0],out = np.zeros_like(test_weights_capped[:, config.weight_derivative_combinations.index(derivative)]),where=test_weights_capped[:,0]!=0) 
                threshold = np.quantile(score_store,1.0-args.rel_max_local_score_test)
                test_weights_capped[:, config.weight_derivative_combinations.index(derivative)][score_store > threshold] = threshold * test_weights_capped[:,0][score_store > threshold]


bsm_colors = [ROOT.kBlack, ROOT.kBlue, ROOT.kGreen, ROOT.kMagenta, ROOT.kCyan]
bit_predictions  = { key:bits[key].vectorized_predict(test_features) for key in  config.bit_derivatives if key!=tuple() }

directory = os.path.join(plot_directory, 'Plot_LLR',args.name_dir,args.binning_var,"Bins_%s"% ('_'.join([str(args.bin_number)])))
if not os.path.exists(directory):
        try:
            os.makedirs(directory)
        except IOError:
            pass



c1 = ROOT.TCanvas()
def compute_weights(weights, theta, theta_ref):
        return weights[:,0] \
                        + np.array( [ (theta[i]-theta_ref[i])*weights[:,config.weight_derivative_combinations.index(('ctZ',))] for i in range(len(theta)) ] ).sum(axis=0)\
                        + np.array( [ 0.5*(theta[i]-theta_ref[i])*(theta[j]-theta_ref[j])*weights[:,config.weight_derivative_combinations.index(('ctZ','ctZ'))] for i in range(len(theta)) for j in range(len(theta)) ] ).sum(axis=0)

def predict_weight_ratio(bit_predictions, theta, theta_ref):
        return 1.+ \
                        + np.array( [ (theta[i]-theta_ref[i])*bit_predictions[('ctZ',)] for i in range(len(theta)) ] ).sum(axis=0)\
                        + np.array( [ 0.5*(theta[i]-theta_ref[i])*(theta[j]-theta_ref[j])*bit_predictions[('ctZ','ctZ')] for i in range(len(theta)) for j in range(len(theta)) ] ).sum(axis=0)

#theta_vec = np.array([])
print "#######################################################"

h_LLR_true = ROOT.TH1D("h_LLR_true","LLR Weights Summed",40,-1.0,1.0)
h_LLR_true_event = ROOT.TH1D("h_LLR_true_event","LLR Eventwise",40,-1.0,1.0)
h_LLR_true_capped = ROOT.TH1D("h_LLR_true_capped","LLR Weights Summed",40,-1.0,1.0)
h_LLR_true_capped_event = ROOT.TH1D("h_LLR_true_capped_event","LLR Eventwise",40,-1.0,1.0)
h_LLR_pred_lin = ROOT.TH1D("h_LLR_pred_lin","LLR Weights Summed",40,-1.0,1.0)
h_LLR_pred_lin_event = ROOT.TH1D("h_LLR_pred_lin","LLR Eventwise",40,-1.0,1.0)
h_LLR_pred_quad = ROOT.TH1D("h_LLR_pred_quad","LLR Weights Summed",40,-1.0,1.0)
h_LLR_pred_quad_event = ROOT.TH1D("h_LLR_pred_quad","LLR Eventwise",40,-1.0,1.0)
#h_LLR_true.GetXAxis().SetTitle(args.binning_var)
#h_LLR_true.GetYAxis().SetTitle("LLR")
for i_plot_theta, plot_theta in enumerate(np.arange(-1.0,1.0,0.05)):
        w0_sum = 0
        w0_capped_sum = 0
        w1_true_sum = 0
        w1_true_capped_sum = 0
        w1_pred_quad_sum = 0
        w1_pred_lin_sum = 0
        LLR_true_event = 0
        LLR_true_capped_event = 0
        LLR_pred_lin_event = 0
        LLR_pred_quad_event = 0
        for j in range(test_features[:,17].size):
            w1_true_inter = test_weights[j,0] + plot_theta*test_weights[j,config.weight_derivative_combinations.index(('ctZ',))] + 0.5 * plot_theta**2 * test_weights[j,config.weight_derivative_combinations.index(('ctZ','ctZ'))]
            w1_true_capped_inter = test_weights_capped[j,0] + plot_theta*test_weights_capped[j,config.weight_derivative_combinations.index(('ctZ',))] + 0.5 * plot_theta**2 * test_weights_capped[j,config.weight_derivative_combinations.index(('ctZ','ctZ'))]
            w0_inter = test_weights[j,0]
            w0_capped_inter = test_weights_capped[j,0]
            div = np.divide(w1_true_inter,w0_inter,out = np.zeros_like(w1_true_inter),where=w0_inter!=0)
#            div_capped = np.divide(w1_true_capped_inter,w0_inter,out = np.zeros_like(w1_true_capped_inter),where=w0_inter!=0)
            LLR_true_event += w1_true_inter - w0_inter -w0_inter*np.log(div,out = np.zeros_like(div),where=div!=0)
            w1_true_sum += w1_true_inter
#            if plot_theta < 0.11 and plot_theta > -0.06 and w1_true_capped_inter <= 0:
#                    print "J: ", j
#                    print "Log: ",np.log(div,out = np.zeros_like(div),where=div!=0), ", Div: ", div
#                    print "Weight Diff: ",w1_true_inter, " Weight: ", w0_inter, " Added:",w1_true_inter - w0_inter -w0_inter*np.log(div,out = np.zeros_like(div),where=div!=0), "LLR: ", LLR_true_event
#            LLR_true_capped_event += w1_true_capped_inter - w0_inter -w0_inter*np.log(div_capped,out = np.zeros_like(div_capped),where=div_capped!=0)
            w1_true_capped_sum += w1_true_capped_inter
            div = np.divide(w1_true_capped_inter,w0_capped_inter,out = np.zeros_like(w1_true_capped_inter),where=w0_capped_inter!=0)
            LLR_true_capped_event += w1_true_capped_inter - w0_capped_inter -w0_capped_inter*np.log(div,out = np.zeros_like(div),where=div>0) #changed - set log to zero if encounter negative div value
#            if plot_theta < 0.11 and plot_theta > -0.06 and w1_true_capped_inter <= 0:
#                    print "Cap Log: ",np.log(div,out = np.zeros_like(div),where=div>0), ", Div: ", div
#                    print "Cap Weight Diff: ",w1_true_capped_inter, " Cap Weight: ", w0_capped_inter, " Added: ", w1_true_capped_inter - w0_capped_inter -w0_capped_inter*np.log(div,out = np.zeros_like(div),where=div>0), "LLR: ", LLR_true_capped_event
            w0_sum += w0_inter
            w0_capped_sum += w0_capped_inter
            w1_pred_quad_inter = test_weights[j,0] + plot_theta * bit_predictions[('ctZ',)][j]*test_weights[j,0] + 0.5*plot_theta**2 * bit_predictions[('ctZ','ctZ')][j]*test_weights[j,0]
            w1_pred_quad_sum += w1_pred_quad_inter
            div = np.divide(w1_pred_quad_inter,w0_inter,out = np.zeros_like(w1_pred_quad_inter),where=w0_inter!=0)
            LLR_pred_quad_event += w1_pred_quad_inter - w0_inter -w0_inter*np.log(div,out = np.zeros_like(div),where=div!=0)
            w1_pred_lin_inter = test_weights[j,0] + plot_theta * bit_predictions[('ctZ',)][j]*test_weights[j,0] 
            w1_pred_lin_sum += w1_pred_lin_inter 
            div = np.divide(w1_pred_lin_inter,w0_inter,out = np.zeros_like(w1_pred_lin_inter),where=w0_inter!=0)
            LLR_pred_lin_event += w1_pred_lin_inter - w0_inter -w0_inter*np.log(div,out = np.zeros_like(div),where=div!=0)
        LLR_true = w1_true_sum - w0_sum - w0_sum*np.log(w1_true_sum/w0_sum)
        LLR_true_capped = w1_true_capped_sum - w0_capped_sum - w0_capped_sum*np.log(w1_true_capped_sum/w0_capped_sum)
        LLR_pred_quad = w1_pred_quad_sum - w0_sum - w0_sum*np.log(w1_pred_quad_sum/w0_sum)
        LLR_pred_lin = w1_pred_lin_sum - w0_sum - w0_sum*np.log(w1_pred_lin_sum/w0_sum)
        print "Theta: ", plot_theta, " LLR True: ", LLR_true, " LLR Pred Quad: ", LLR_pred_quad, " LLR Pred Lin: ", LLR_pred_lin, " LLR Cap: ", LLR_true_capped
        h_LLR_true.Fill(plot_theta+0.01,LLR_true)
        h_LLR_true_event.Fill(plot_theta+0.01,LLR_true_event)
        h_LLR_true_capped.Fill(plot_theta+0.01,LLR_true_capped)
        h_LLR_true_capped_event.Fill(plot_theta+0.01,LLR_true_capped_event)
        h_LLR_pred_quad_event.Fill(plot_theta+0.01,LLR_pred_quad_event)
        h_LLR_pred_lin_event.Fill(plot_theta+0.01,LLR_pred_lin_event)
        h_LLR_pred_quad.Fill(plot_theta+0.01,LLR_pred_quad)
        h_LLR_pred_lin.Fill(plot_theta+0.01,LLR_pred_lin)

print "########################################################"
h_LLR_true.SetTitle("LLR Weights Summed;ctZ;LLR")
h_LLR_true.GetYaxis().SetTitleOffset(1.4)
#h_LLR_true.GetXaxis().SetTitleSize(0.15)
#h_LLR_true.GetYaxis().SetTitleSize(0.15)
h_LLR_true.SetLineColor(ROOT.kRed)
h_LLR_true_event.SetLineColor(ROOT.kRed)
h_LLR_true_capped.SetLineColor(ROOT.kRed)
h_LLR_true_capped_event.SetLineColor(ROOT.kRed)
h_LLR_true_capped.SetLineStyle(2)
h_LLR_true_capped_event.SetLineStyle(2)
h_LLR_pred_lin.SetLineColor(ROOT.kGreen)
h_LLR_pred_lin_event.SetLineColor(ROOT.kGreen)
h_LLR_pred_quad.SetLineColor(ROOT.kBlue)
h_LLR_pred_quad_event.SetLineColor(ROOT.kBlue)
legend = TLegend(0.65,0.82,0.95,0.92)
legend.SetFillStyle(0)
legend.SetShadowColor(ROOT.kWhite)
legend.SetBorderSize(0)
legend.AddEntry(h_LLR_true,"LLR true","l")
legend.AddEntry(h_LLR_true_capped,"LLR true capped","l")
legend.AddEntry(h_LLR_pred_lin,"LLR linear prediction","l")
legend.AddEntry(h_LLR_pred_quad,"LLR quadratic prediction","l")
h_LLR_true.Draw("hist")
h_LLR_true_capped.Draw("histsame")
h_LLR_pred_lin.Draw("histsame")
h_LLR_pred_quad.Draw("histsame")
legend.Draw()
filename = "LLR.png"
c1.Print(os.path.join(directory,filename))
h_LLR_true_capped.Draw("hist")
filename = "LLR_Capped.png"
c1.Print(os.path.join(directory,filename))

h_LLR_true_event.SetTitle("LLR Eventwise;ctZ;LLR")
legend = TLegend(0.65,0.82,0.95,0.92)
legend.SetFillStyle(0)
legend.SetShadowColor(ROOT.kWhite)
legend.SetBorderSize(0)
legend.AddEntry(h_LLR_true_event,"LLR true","l")
legend.AddEntry(h_LLR_true_capped_event,"LLR true capped","l")
legend.AddEntry(h_LLR_pred_lin_event,"LLR linear prediction","l")
legend.AddEntry(h_LLR_pred_quad_event,"LLR quadratic prediction","l")
h_LLR_true_event.Draw("hist")
h_LLR_true_capped_event.Draw("histsame")
h_LLR_pred_lin_event.Draw("histsame")
h_LLR_pred_quad_event.Draw("histsame")
legend.Draw()
filename = "LLR_eventwise.png"
c1.Print(os.path.join(directory,filename))
h_LLR_true_capped_event.Draw("hist")
filename = "LLR_eventwise_capped.png"
c1.Print(os.path.join(directory,filename))

h_LLR_true = ROOT.TH1D("h_LLR_true","LLR True",40,-1.0,1.0)
h_LLR_true_capped = ROOT.TH1D("h_LLR_true_capped","LLR True Capped",40,-1.0,1.0)
h_LLR_pred_lin = ROOT.TH1D("h_LLR_pred_lin","LLR Linear Prediction",40,-1.0,1.0)
h_LLR_pred_quad = ROOT.TH1D("h_LLR_pred_quad","LLR Quadratic Prediction",40,-1.0,1.0)
print "##########################################################"
print "BINNED:"
bin_number = args.bin_number
binning_var = args.binning_var
min_bin = args.min_bin
max_bin = args.max_bin
#print np.where(mva_variables == binning_var)
#bin_index = np.where(mva_variables == binning_var)[0][0]
bin_index = mva_variables.index(binning_var)
for i_plot_theta, plot_theta in enumerate(np.arange(-1.0,1.0,0.05)):
        h_w1_pred_lin = ROOT.TH1D("h_w1_pred_lin","Linear Prediction",bin_number,min_bin,max_bin)
        h_w1_pred_quad = ROOT.TH1D("h_w1_pred_quad","Quadratic Prediction",bin_number,min_bin,max_bin)
        h_w1_true = ROOT.TH1D("h_w1_true","True",bin_number,min_bin,max_bin)
        h_w1_true_capped = ROOT.TH1D("h_w1_true_capped","True Capped",bin_number,min_bin,max_bin)
        h_w0 = ROOT.TH1D("h_w0", "Weight",bin_number,min_bin,max_bin)
        w0 = 0
        w1_true = 0
        w1_true_capped = 0
        w1_pred_quad = 0
        w1_pred_lin = 0
        for j in range(test_features[:,bin_index].size):
            w1_true = test_weights[j,0] + plot_theta*test_weights[j,config.weight_derivative_combinations.index(('ctZ',))] + 0.5 * plot_theta**2 * test_weights[j,config.weight_derivative_combinations.index(('ctZ','ctZ'))]
            w1_true_capped = test_weights_capped[j,0] + plot_theta*test_weights_capped[j,config.weight_derivative_combinations.index(('ctZ',))] + 0.5 * plot_theta**2 * test_weights_capped[j,config.weight_derivative_combinations.index(('ctZ','ctZ'))]
            w0 = test_weights[j,0]
            w1_pred_quad = test_weights[j,0] + plot_theta * bit_predictions[('ctZ',)][j]*test_weights[j,0] + 0.5*plot_theta**2 * bit_predictions[('ctZ','ctZ')][j]*test_weights[j,0]
            w1_pred_lin = test_weights[j,0] + plot_theta * bit_predictions[('ctZ',)][j]*test_weights[j,0] 
            h_w1_pred_lin.Fill(test_features[j,bin_index],w1_pred_lin)
            h_w1_pred_quad.Fill(test_features[j,bin_index],w1_pred_quad)
            h_w1_true.Fill(test_features[j,bin_index],w1_true)
            h_w1_true_capped.Fill(test_features[j,bin_index],w1_true_capped)
            h_w0.Fill(test_features[j,bin_index],w0)
        h_w1_pred_lin_clone = h_w1_pred_lin.Clone("h_w1_pred_lin_clone")
        h_w1_pred_quad_clone = h_w1_pred_quad.Clone("h_w1_pred_quad_clone")
        h_w1_true_clone = h_w1_true.Clone("h_w1_true_clone")
        h_w1_true_capped_clone = h_w1_true_capped.Clone("h_w1_true_capped_clone")
        h_w0_clone = h_w0.Clone("h_w0_clone")
        h_w0_clone_pred_lin = h_w0.Clone("h_w0_clone_pred_lin")
        h_w0_clone_pred_quad = h_w0.Clone("h_w0_clone_pred_quad")
        h_w0_clone_true = h_w0.Clone("h_w0_clone_true")
        h_w0_clone_true_capped = h_w0.Clone("h_w0_clone_true_capped")
        h_w1_pred_lin_clone.Divide(h_w0_clone)
        h_w1_pred_quad_clone.Divide(h_w0_clone)
        h_w1_true_clone.Divide(h_w0_clone)
        h_w1_true_capped_clone.Divide(h_w0_clone)
        for k in range(bin_number):
                inter = h_w1_pred_lin_clone.GetBinContent(k+1)
                h_w1_pred_lin_clone.SetBinContent(k+1,np.log(inter,out = np.zeros_like(inter),where=inter!=0))
                inter = h_w1_pred_quad_clone.GetBinContent(k+1)
                h_w1_pred_quad_clone.SetBinContent(k+1,np.log(inter,out = np.zeros_like(inter),where=inter!=0))
                inter = h_w1_true_clone.GetBinContent(k+1)
                h_w1_true_clone.SetBinContent(k+1,np.log(inter,out = np.zeros_like(inter),where=inter!=0))
                inter = h_w1_true_capped_clone.GetBinContent(k+1)
                h_w1_true_capped_clone.SetBinContent(k+1,np.log(inter,out = np.zeros_like(inter),where=inter!=0))
        h_w0_clone_pred_lin.Multiply(h_w1_pred_lin_clone)
        h_w0_clone_pred_quad.Multiply(h_w1_pred_quad_clone)
        h_w0_clone_true.Multiply(h_w1_true_clone)
        h_w0_clone_true_capped.Multiply(h_w1_true_capped_clone)
        h_w0_clone = h_w0.Clone("h_w0_clone")
        h_w0_clone.Add(h_w0_clone_pred_lin)
        h_w1_pred_lin.Add(h_w0_clone,-1)
        h_w0_clone = h_w0.Clone("h_w0_clone")
        h_w0_clone.Add(h_w0_clone_pred_quad)
        h_w1_pred_quad.Add(h_w0_clone,-1)
        h_w0_clone = h_w0.Clone("h_w0_clone")
        h_w0_clone.Add(h_w0_clone_true)
        h_w1_true.Add(h_w0_clone,-1)
        h_w0_clone = h_w0.Clone("h_w0_clone")
        h_w0_clone.Add(h_w0_clone_true_capped)
        h_w1_true_capped.Add(h_w0_clone,-1)
        LLR_pred_lin = h_w1_pred_lin.Clone("LLR_pred_lin")
        LLR_pred_quad = h_w1_pred_quad.Clone("LLR_pred_quad")
        LLR_true = h_w1_true.Clone("LLR_true")
        LLR_true_capped = h_w1_true_capped.Clone("LLR_true_capped")

        #LLR_true.legendText = "LLR (true)"
        #LLR_pred_lin.legendText = "LLR lin (pred)"
        #LLR_pred_quad.legendText = "LLR quad (pred)"
        #LLR_true.style = styles.lineStyle(ROOT.kRed)
        #LLR_pred_lin.style = styles.lineStyle(ROOT.kBlue)
        #LLR_pred_quad.style = styles.lineStyle(ROOT.kGreen)
        LLR_true.SetLineColor(ROOT.kRed)
        LLR_true_capped.SetLineColor(ROOT.kRed)
        LLR_true_capped.SetLineStyle(2)
        #LLR_true.UseCurrentStyle()
        LLR_pred_lin.SetLineColor(ROOT.kGreen)
        #LLR_pred_lin.UseCurrentStyle()
        LLR_pred_quad.SetLineColor(ROOT.kBlue)
        #LLR_pred_quad.UseCurrentStyle()
        true_inter = 0
        true_capped_inter = 0
        lin_inter = 0
        quad_inter = 0
        for k in range(bin_number):
                true_inter += LLR_true.GetBinContent(k+1)
                true_capped_inter += LLR_true_capped.GetBinContent(k+1)
                lin_inter += LLR_pred_lin.GetBinContent(k+1)
                quad_inter += LLR_pred_quad.GetBinContent(k+1)
        h_LLR_true.Fill(plot_theta+0.01,true_inter)
        h_LLR_true_capped.Fill(plot_theta+0.01,true_capped_inter)
        h_LLR_pred_quad.Fill(plot_theta+0.01,quad_inter)
        h_LLR_pred_lin.Fill(plot_theta+0.01,lin_inter)
        if binning_var == "mva_photon_pt":
                LLR_true.GetXaxis().SetTitle("p_{T}(#gamma) (GeV)")
        elif binning_var == "mva_l1_pt":
                LLR_true.GetXaxis().SetTitle("p_{T}(l_{1}) (GeV)")
        elif binning_var == "mva_nBTag":
                LLR_true.GetXaxis().SetTitle("N_{b-tag}")
        elif binning_var == "mva_nJetGood":
                LLR_true.GetXaxis().SetTitle("N_{jets}")
        else:
                LLR_true.GetXaxis().SetTitle(binning_var)
        LLR_true.GetYaxis().SetTitle("LLR")
        legend = TLegend(0.65,0.82,0.95,0.92)
        legend.SetFillStyle(0)
        legend.SetShadowColor(ROOT.kWhite)
        legend.SetBorderSize(0)
        legend.AddEntry(LLR_true,"LLR true","l")
        legend.AddEntry(LLR_true_capped,"LLR true capped","l")
        legend.AddEntry(LLR_pred_lin,"LLR linear prediction","l")
        legend.AddEntry(LLR_pred_quad,"LLR quadratic prediction","l")

        LLR_true.Draw("hist")
        LLR_true_capped.Draw("histsame")
        LLR_pred_lin.Draw("histsame")
        LLR_pred_quad.Draw("histsame")
        legend.Draw()
        filename = "LLR_ctZ_%s"% ('_'.join([str(round(plot_theta,2))])) + ".png"
        c1.Print(os.path.join(directory,filename))


h_LLR_true.GetXaxis().SetTitle("ctZ")
h_LLR_true.GetYaxis().SetTitle("LLR")
legend = TLegend(0.65,0.82,0.95,0.92)
legend.SetFillStyle(0)
legend.SetShadowColor(ROOT.kWhite)
legend.SetBorderSize(0)
legend.AddEntry(h_LLR_true,"LLR true","l")
legend.AddEntry(h_LLR_true_capped,"LLR true capped","l")
legend.AddEntry(h_LLR_pred_lin,"LLR linear prediction","l")
legend.AddEntry(h_LLR_pred_quad,"LLR quadratic prediction","l")
h_LLR_true.SetLineColor(ROOT.kRed)
h_LLR_true_capped.SetLineColor(ROOT.kRed)
h_LLR_true_capped.SetLineStyle(2)
h_LLR_pred_lin.SetLineColor(ROOT.kGreen)
h_LLR_pred_quad.SetLineColor(ROOT.kBlue)
h_LLR_true.Draw("hist")
h_LLR_true_capped.Draw("histsame")
h_LLR_pred_lin.Draw("histsame")
h_LLR_pred_quad.Draw("histsame")
legend.Draw()
filename = "LLR_binned.png"
c1.Print(os.path.join(directory,filename))
#print binning_var
#print bin_index
#print mva_variables
#        histos = [[LLR_true], [LLR_pred_lin], [LLR_pred_quad]]
#        plot = Plot.fromHisto( "LLR_Theta_%s"% ('_'.join([str(plot_theta)])), histos, texX="p_{T}(#gamma) (GeV)", texY="LLR" )
#        plotting.draw(plot, 
#            plot_directory = plot_directory,
#            yRange = 'auto',
            #legend = None,
            #ratio = {'yRange':(0.5,1.5)}, logY = False, logX = False, 
#            drawObjects = [h_nom]
#        )


        

#        ext_Delta_NLL[pd.isna(ext_Delta_NLL)] = 0.0
#        ext_Delta_NLL_pred = w1 - w0 - w0*np.log(predict_weight_ratio(bit_predictions, plot_theta, theta_ref))
#        print plot_theta, "true", ext_Delta_NLL.sum(), "pred", ext_Delta_NLL_pred.sum()
h_lin_pred = ROOT.TH1D("h_lin_pred","Linear Prediction",20,0,400)
h_quad_pred = ROOT.TH1D("h_quad_pred","Quadratic Prediction",20,0,400)
h_lin_true = ROOT.TH1D("h_lin_true","Linear True",20,0,400)
h_quad_true = ROOT.TH1D("h_quad_true","Quadratic True",20,0,400)
h_weight_0 = ROOT.TH1D("h_quad_true","Quadratic True",20,0,400)

pred_lin_total = 0.0
pred_quad_total = 0.0
for i in range(test_features[:,17].size):
        pred_lin_inter =bit_predictions[('ctZ',)][i]*test_weights[i,0]
        pred_quad_inter =bit_predictions[('ctZ','ctZ')][i]*test_weights[i,0]
        true_lin_inter =test_weights[i,config.weight_derivative_combinations.index(('ctZ',))]
        true_quad_inter =test_weights[i,config.weight_derivative_combinations.index(('ctZ','ctZ'))] 
        h_lin_pred.Fill(test_features[i,17],pred_lin_inter)
        h_quad_pred.Fill(test_features[i,17],pred_quad_inter)
        h_lin_true.Fill(test_features[i,17],true_lin_inter)
        h_quad_true.Fill(test_features[i,17],true_quad_inter)
        h_weight_0.Fill(test_features[i,17],test_weights[i,0])

#cloning histogramm
#h_new = h.Clone("h_new")

#for i in range(bin_weight_0.size):
#        h_lin_pred.Fill(test_features[i,17],bin_pred_lin[i]/bin_weight_0[i])
#        h_quad_pred.Fill(test_features[i,17],bin_pred_quad[i]/bin_weight_0[i])
#        h_lin_true.Fill(test_features[i,17],bin_true_lin[i]/bin_weight_0[i])
#        h_quad_true.Fill(test_features[i,17],bin_true_quad[i]/bin_weight_0[i])

h_quad_pred.legendText = "s_{tZ} (pred.)"
h_quad_true.legendText = "s_{tZ} (true)"
h_quad_pred.style = styles.lineStyle(ROOT.kRed)
h_quad_true.style = styles.lineStyle(ROOT.kBlue)

histos = [[h_quad_true], [h_quad_pred]]
plot   = Plot.fromHisto( "closure",  histos, texX = "p_{T}(#gamma) (GeV)", texY = "s_{tZ}" )
plotting.draw(plot, 
        plot_directory = directory,
        yRange = 'auto',
        #legend = None,
        ratio = {'yRange':(0.5,1.5)}, logY = False, logX = False, 
        #drawObjects = [h_nom]
    )
#h_lin_pred.Draw("hist")
#filename2 = "h_lin_pred_without.png"
#c1.Print(os.path.join(plot_directory, 'Plot_LLR',args.name_dir,filename2))
#h_quad_pred.Draw("hist")
#filename2 = "h_quad_pred_without.png"
#c1.Print(os.path.join(plot_directory, 'Plot_LLR',args.name_dir,filename2))
#h_lin_true.Draw("hist")
#filename2 = "h_lin_true_without.png"
#c1.Print(os.path.join(plot_directory, 'Plot_LLR',args.name_dir,filename2))
#h_quad_true.Draw("hist")
#filename2 = "h_quad_true_without.png"
#c1.Print(os.path.join(plot_directory, 'Plot_LLR',args.name_dir,filename2))
#
#h_lin_pred.Divide(h_weight_0)
#h_quad_pred.Divide(h_weight_0)
#h_lin_true.Divide(h_weight_0)
#h_quad_true.Divide(h_weight_0)
#
#
#h_lin_pred.Draw("hist")
#filename2 = "h_lin_pred.png"
#c1.Print(os.path.join(plot_directory, 'Plot_LLR',args.name_dir,filename2))
#h_quad_pred.Draw("hist")
#filename2 = "h_quad_pred.png"
#c1.Print(os.path.join(plot_directory, 'Plot_LLR',args.name_dir,filename2))
#h_lin_true.Draw("hist")
#filename2 = "h_lin_true.png"
#c1.Print(os.path.join(plot_directory, 'Plot_LLR',args.name_dir,filename2))
#h_quad_true.Draw("hist")
#filename2 = "h_quad_true.png"
#c1.Print(os.path.join(plot_directory, 'Plot_LLR',args.name_dir,filename2))
#
#h_lin_pred.Divide(h_lin_true)
#h_quad_pred.Divide(h_quad_true)
#
#h_lin_pred.Draw("hist")
#filename2 = "h_lin_pred_vs_true.png"
#c1.Print(os.path.join(plot_directory, 'Plot_LLR',args.name_dir,filename2))
#h_quad_pred.Draw("hist")
#filename2 = "h_quad_pred_vs_true.png"
#c1.Print(os.path.join(plot_directory, 'Plot_LLR',args.name_dir,filename2))
#
#
##h2 = ROOT.TH1D("h2", "Histogramm New", 10,0,100)
##x_list = np.random.uniform(0,100,1000)
##w_list = np.random.uniform(0,1,1000)
##for i in range(x_list.size):
##    h2.Fill(x_list[i],w_list[i])
##h2.Draw("hist") #Draw Histogramm
##filename = "test_histo.png"
##c1.Print(os.path.join(plot_directory, 'Plot_LLR',filename))
#
#print "Success"

#def compute_weights(weights, theta, theta_ref, param):
#        return weights[()] \
#                        + np.array( [ (theta-theta_ref)*weights[(param,)]]).sum(axis=0)\
#                        #+ np.array( [ 0.5*(theta-theta_ref)*(theta-theta_ref)*weights[(param, param)]] ).sum(axis=0)

#def predict_weight_ratio(bit_predictions, theta, theta_ref, param):
#        lin       = np.array( [ (theta-theta_ref)*bit_predictions[(param,)] ] ).sum(axis=0)
#        #quadratic = np.array( [ 0.5*(theta-theta_ref)*(theta-theta_ref)*bit_predictions[(param, param)] ] ).sum(axis=0)
#        return 1.+lin #, 1.+lin+quadratic

#for i_plot_theta, plot_theta in enumerate(np.arange(-.03,.05,.002)):
#        param = 'ctZ'
#        w1 = compute_weights( test_weights, plot_theta, 0, param)
#        w0 = compute_weights( test_weights, 0,0, param)

#        ext_Delta_NLL      = w1 - w0 - w0*np.log(w1/w0)

#        lin = np.log(predict_weight_ratio(bit_predictions, plot_theta, 0, param))
        #ext_Delta_NLL_pred     = w1 - w0 - w0*quad
#        ext_Delta_NLL_pred_lin = w1 - w0 - w0*lin

#        print plot_theta, "true", ext_Delta_NLL.sum(), "lin", ext_Delta_NLL_pred_lin.sum()

#        # Histo style
#        test_FI_histo    .style = styles.lineStyle( ROOT.kBlue, width=2 )
#        training_FI_histo.style = styles.lineStyle( ROOT.kRed, width=2 )
#        test_FI_histo    .legendText = "Test"
#        training_FI_histo.legendText = "Training"
#
#        minY   = 0.01 * min( test_FI_histo.GetBinContent(test_FI_histo.GetMaximumBin()), training_FI_histo.GetBinContent(training_FI_histo.GetMaximumBin()))
#        maxY   = 1.5  * max( test_FI_histo.GetBinContent(test_FI_histo.GetMaximumBin()), training_FI_histo.GetBinContent(training_FI_histo.GetMaximumBin()))
#
#        histos = [ [test_FI_histo], [training_FI_histo] ]
#        plot   = Plot.fromHisto( "bit_derivative_%s_evolution"% ('_'.join(derivative)), histos, texX="b", texY="L(D,b)" )
#
#        # Plot Style
#        histModifications      = []
#        histModifications      += [ lambda h: h.GetYaxis().SetTitleOffset(1.4) ]
#        histModifications += [ lambda h: h.GetXaxis().SetTitleSize(26) ]
#        histModifications += [ lambda h: h.GetYaxis().SetTitleSize(26) ]
#        histModifications += [ lambda h: h.GetXaxis().SetLabelSize(22)  ]
#        histModifications += [ lambda h: h.GetYaxis().SetLabelSize(22)  ]
#
#        ratioHistModifications = []
#        ratio                  = None
#        legend                 = (0.6,0.75,0.9,0.88)
#        yRange                 = "auto" #( minY, maxY )
#        plot1DHist( plot, plot_directory, yRange=yRange, ratio=ratio, legend=legend, histModifications=histModifications )
#
