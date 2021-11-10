#!/usr/bin/env python

import ROOT, os, sys
from RootTools.core.standard import *
import Analysis.Tools.syncer as syncer

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

argParser.add_argument('--name_dir',               action='store', type=str,   default='default', help="Name of Plot Directory")

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

import Analysis.Tools.user as user
# directories
plot_directory   = os.path.join( user. plot_directory, 'MVA', args.config, args.name)
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

        n_trees = config.bit_cfg['n_trees']

        time1 = time.time()
        print ("Learning %s"%( str(derivative)))
        bits[derivative]= BoostedInformationTree(
                training_features     = training_features,
                training_weights      = training_weights[:,0],
                training_diff_weights = training_weights[:, config.weight_derivative_combinations.index(derivative)],
                split_method          = 'vectorized_split_and_weight_sums',
                weights_update_method = 'vectorized',
                bagging_fraction      = args.bagging_fraction,
                **config.bit_cfg
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

        if args.debug:
            from debug import make_debug_plots
            make_debug_plots( bits[derivative], 
                              training_features, training_weights[:,0],  
                              training_weights[:, config.weight_derivative_combinations.index(derivative)], 
                              test_features, 
                              test_weights[:,0], 
                              test_weights_, 
                              os.path.join(plot_directory, ('_'.join(derivative))),
                              mva_variables = config.mva_variables) 

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
        test_FI_histo    .style = styles.lineStyle( ROOT.kBlue, width=2 )
        training_FI_histo.style = styles.lineStyle( ROOT.kRed, width=2 )
        test_FI_histo    .legendText = "Test"
        training_FI_histo.legendText = "Training"

        minY   = 0.01 * min( test_FI_histo.GetBinContent(test_FI_histo.GetMaximumBin()), training_FI_histo.GetBinContent(training_FI_histo.GetMaximumBin()))
        maxY   = 1.5  * max( test_FI_histo.GetBinContent(test_FI_histo.GetMaximumBin()), training_FI_histo.GetBinContent(training_FI_histo.GetMaximumBin()))

        histos = [ [test_FI_histo], [training_FI_histo] ]
        plot   = Plot.fromHisto( "bit_derivative_%s_evolution"% ('_'.join(derivative)), histos, texX="b", texY="L(D,b)" )

        # Plot Style
        histModifications      = []
        histModifications      += [ lambda h: h.GetYaxis().SetTitleOffset(1.4) ]
        histModifications += [ lambda h: h.GetXaxis().SetTitleSize(26) ]
        histModifications += [ lambda h: h.GetYaxis().SetTitleSize(26) ]
        histModifications += [ lambda h: h.GetXaxis().SetLabelSize(22)  ]
        histModifications += [ lambda h: h.GetYaxis().SetLabelSize(22)  ]

        ratioHistModifications = []
        ratio                  = None
        legend                 = (0.6,0.75,0.9,0.88)
        yRange                 = "auto" #( minY, maxY )
        plot1DHist( plot, plot_directory, yRange=yRange, ratio=ratio, legend=legend, histModifications=histModifications )

# make score plot and score-binned feature plots
#for derivative in config.bit_derivatives:
#    if derivative == tuple(): continue
#    if args.derivative!=[] and derivative !=tuple(args.derivative): continue

    # plot loss
#    test_scores     = bits[derivative].vectorized_predict(test_features)
    #training_scores = bits[derivative].vectorized_predict(test_features)

#    n_digi    = 10
#    quantiles = np.quantile( test_scores,np.arange(0,1.1,.1))
#    try:
#        digi      = np.digitize( test_scores, quantiles)
#    except ValueError:
#        continue

#    h_score   = ROOT.TH1F("score","score", len(quantiles)-1, quantiles)
#    h_feature = { var[0]: {i_bin:ROOT.TH1F("%s_%i"%(var[0], i_bin), "%s_%i"%(var[0], i_bin), 50, min(test_features[:,i_var]), max(test_features[:,i_var])) for i_bin in range(n_digi)} for i_var, var in enumerate(config.mva_variables) }

#    test_weights_ = test_weights[:,config.weight_derivative_combinations.index(derivative)]
#    for i_event in range(len(test_scores)):
#        if i_event%10000==0: print "At",i_event
#        h_score.Fill(test_scores[i_event], test_weights_[i_event])
#        for i_var in range(len(config.mva_variables)):
#            h_feature[config.mva_variables[i_var][0]][min(n_digi,digi[i_event])-1].Fill(test_features[i_event][i_var], test_weights_[i_event])

    # Plot Style
#    histModifications      = []
#    histModifications      += [ lambda h: h.GetYaxis().SetTitleOffset(1.4) ]
#    histModifications += [ lambda h: h.GetXaxis().SetTitleSize(26) ]
#    histModifications += [ lambda h: h.GetYaxis().SetTitleSize(26) ]
#    histModifications += [ lambda h: h.GetXaxis().SetLabelSize(22)  ]
#    histModifications += [ lambda h: h.GetYaxis().SetLabelSize(22)  ]

#    plot = Plot.fromHisto( 
#        "score", 
#        [[h_score]],
#        texX = "d_"+"_".join(derivative),
#        )

#    ratio                  = None
#    legend                 = None #[ (0.15,0.75,0.9,0.88), 3]
#    yRange                 = "auto" #( minY, maxY )
#    plot1DHist( plot, os.path.join(plot_directory, ('_'.join(derivative))), yRange=yRange, ratio=ratio, legend=legend, histModifications=histModifications )

#    for i_var in range(len(config.mva_variables)):
#        for i_digi in range(n_digi):
#            h_feature[config.mva_variables[i_var][0]][i_digi].style     = styles.fillStyle(ROOT.kMagenta-10+i_digi)
#            h_feature[config.mva_variables[i_var][0]][i_digi].legendText= "%4.3f #leq d(%s) < %4.3f" % ( quantiles[i_digi], "_".join(derivative), quantiles[i_digi+1] )
#        plot = Plot.fromHisto( 
#            config.mva_variables[i_var][0]+"_binned", 
#            [ [h_feature[config.mva_variables[i_var][0]][i_digi] for i_digi in range(n_digi)]],
#            texX = config.mva_variables[i_var][0],
#            )

#        ratio                  = None
#        legend                 = [ (0.15,0.75,0.9,0.88), 2]
#        yRange                 = "auto" #( minY, maxY )
#        plot1DHist( plot, os.path.join(plot_directory, ('_'.join(derivative))), yRange=yRange, ratio=ratio, legend=legend, histModifications=histModifications )

bsm_colors = [ROOT.kBlack, ROOT.kBlue, ROOT.kGreen, ROOT.kMagenta, ROOT.kCyan]
#for derivative in config.bit_derivatives:
#        print "Derivative: ", derivative
#bit_predictions = [bits[('ctZ',)].vectorized_predict(test_features), bits[('ctZ','ctZ')].vectorized_predict(test_features)]
bit_predictions  = { key:bits[key].vectorized_predict(test_features) for key in  config.bit_derivatives if key!=tuple() }
#print "#################################################################"
#print "Bit Pred", bit_predictions
#print "Bit Pred ctZ", bit_predictions[('ctZ',)]
#print "Bit Pred ctZ, ctZ", bit_predictions[('ctZ','ctZ')]
#print "Bit Pred ctZ Sum", bit_predictions[('ctZ',)].sum(axis = 0)
#print "Bit Pred ctZ, ctZ sum", bit_predictions[('ctZ','ctZ')].sum(axis=0)
#print "Bit Pred ctZ Length", bit_predictions[('ctZ',)].size
#print "Bit Pred ctZ, ctZ Length", bit_predictions[('ctZ','ctZ')].size
#print "test weight size", test_weights.size
#print "test_weights", test_weights
#print "test_weights ctZ", test_weights[:, config.weight_derivative_combinations.index(('ctZ',))]
#print "test_weights ctZ Index", config.weight_derivative_combinations.index(('ctZ',))
#print "test_weights ctZ ctZ", test_weights[:, config.weight_derivative_combinations.index(('ctZ','ctZ'))]
#print "test_weights ctZ ctZ Index", config.weight_derivative_combinations.index(('ctZ','ctZ'))
#print "test_weights ctZ Sum", test_weights[:, config.weight_derivative_combinations.index(('ctZ',))].sum()
#print "test_weights ctZ ctZ Sum", test_weights[:, config.weight_derivative_combinations.index(('ctZ','ctZ'))].sum()
#print "test features", test_features
#print "test features - variables", mva_variables
#print "test features - variables length", len(mva_variables)
#print "test features length (1 event)", test_features[1,:].size
#print "test features length (1 variable)", test_features[:,1].size

#bit_predictions  = { key:bits[key].vectorized_predict(test_features) for key in  [('ctZ',),('ctZ','ctZ')] if key!=tuple() }
def compute_weights(weights, theta, theta_ref):
        return weights[:,0] \
                        + np.array( [ (theta[i]-theta_ref[i])*weights[:,config.weight_derivative_combinations.index(('ctZ',))] for i in range(len(theta)) ] ).sum(axis=0)\
                        + np.array( [ 0.5*(theta[i]-theta_ref[i])*(theta[j]-theta_ref[j])*weights[:,config.weight_derivative_combinations.index(('ctZ','ctZ'))] for i in range(len(theta)) for j in range(len(theta)) ] ).sum(axis=0)

def predict_weight_ratio(bit_predictions, theta, theta_ref):
        return 1.+ \
                        + np.array( [ (theta[i]-theta_ref[i])*bit_predictions[('ctZ',)] for i in range(len(theta)) ] ).sum(axis=0)\
                        + np.array( [ 0.5*(theta[i]-theta_ref[i])*(theta[j]-theta_ref[j])*bit_predictions[('ctZ','ctZ')] for i in range(len(theta)) for j in range(len(theta)) ] ).sum(axis=0)

theta_ref = np.array([0.])

true = np.array([])
pred = np.array([])
theta_vec = np.array([])
for i_plot_theta, plot_theta in enumerate(np.arange(-1.0,1.0,0.05).reshape(-1,1)):
        w1 = compute_weights( test_weights, plot_theta, theta_ref )
        w0 = compute_weights( test_weights, theta_ref , theta_ref )
        #w0[np.absolute(w0) < 1.*10**-14] = 0.0
        #print "weight 0: ", w0, "weight 1: ", w1
        div      = np.divide(w1,w0,out = np.zeros_like(w1),where=w0!=0)

        ext_Delta_NLL      = w1 - w0 - w0*np.log(div,out = np.zeros_like(div),where=div!=0)
        ext_Delta_NLL[pd.isna(ext_Delta_NLL)] = 0.0
#        print "#################################################################################################"
#        print "ext Delta: ", ext_Delta_NLL
#        print ", Div: ",np.divide(w1,w0,out = np.zeros_like(w1),where = w0!=0)
#        print "weight 0: ", w0
#        print "weight 1: ", w1
        ext_Delta_NLL_pred = w1 - w0 - w0*np.log(predict_weight_ratio(bit_predictions, plot_theta, theta_ref))
        print plot_theta, "true", ext_Delta_NLL.sum(), "pred", ext_Delta_NLL_pred.sum()
#        print "true ratio (first value): ", div[1], "pred ratio: ", predict_weight_ratio(bit_predictions,plot_theta,theta_ref)[1]
#        print "true ratio (second value): ", div[2], "pred ratio: ", predict_weight_ratio(bit_predictions,plot_theta,theta_ref)[2]
#        print "true ratio (third value): ", div[3], "pred ratio: ", predict_weight_ratio(bit_predictions,plot_theta,theta_ref)[3]
#print "#########################################################################################"
#for i in range(len(bit_predictions[('ctZ',)])):
#                print "ctZ predict: ", bit_predictions[('ctZ',)][i], " weights: ", test_weights[i,config.weight_derivative_combinations.index(('ctZ',))]
#                print "ctZ ctZ predict: ", bit_predictions[('ctZ','ctZ')][i], " weights: ", test_weights[i,config.weight_derivative_combinations.index(('ctZ','ctZ'))]
#print "ctZ predict: ", bit_predictions[('ctZ',)], " weights: ", test_weights[:,config.weight_derivative_combinations.index(('ctZ',))]/ test_weights[:,0]
#print "ctZ ctZ predict: ", bit_predictions[('ctZ','ctZ')], " weights: ", test_weights[:,config.weight_derivative_combinations.index(('ctZ','ctZ'))]/ test_weights[:,0]
#print "#####################################################################"
#print "ctZ predict one: ", bit_predictions[('ctZ',)][1], " weights: ", test_weights[1,config.weight_derivative_combinations.index(('ctZ',))]/ test_weights[:,0]
#print "ctZ ctZ predict one: ", bit_predictions[('ctZ','ctZ')][1], " weights: ", test_weights[1,config.weight_derivative_combinations.index(('ctZ','ctZ'))]/ test_weights[:,0]
#print "ctZ predict two: ", bit_predictions[('ctZ',)][2], " weights: ", test_weights[2,config.weight_derivative_combinations.index(('ctZ',))]/ test_weights[:,0]
#print "ctZ ctZ predict two: ", bit_predictions[('ctZ','ctZ')][2], " weights: ", test_weights[2,config.weight_derivative_combinations.index(('ctZ','ctZ'))]/ test_weights[:,0]
#print "##########################################################################################################"
#counter = 0
#for i in range(test_features[:,17].size):
#        if test_features[i,17] > 390.:
#                print "Event: ", i, "Photon pt: ", test_features[i,17], " Bit-Pred Lin: ", bit_predictions[('ctZ',)][i]," Bit-Pred Quad: ", bit_predictions[('ctZ','ctZ')][i], " weights lin: ", test_weights[i,config.weight_derivative_combinations.index(('ctZ',))]/ test_weights[i,0], " weights quad: ", test_weights[i,config.weight_derivative_combinations.index(('ctZ','ctZ'))]/ test_weights[i,0]
#                print "########################"
#                counter += 1
#                if counter >= 20:
#                        break
#bin_150_200_pred_lin = np.array([0.,0.,0.,0.,0.])
#bin_150_200_true_lin = np.array([0.,0.,0.,0.,0.])
#bin_150_200_pred_quad = np.array([0.,0.,0.,0.,0.])
#bin_150_200_true_quad = np.array([0.,0.,0.,0.,0.])
#bin_150_200_weight_0 = np.array([0.,0.,0.,0.,0.])
#for i in range(test_features[:,17].size):
#        if test_features[i,17] >= 150. and test_features[i,17] < 200.:
#                pred_lin_inter =bit_predictions[('ctZ',)][i]*test_weights[i,0]
#                pred_quad_inter =bit_predictions[('ctZ','ctZ')][i]*test_weights[i,0]
#                true_lin_inter =test_weights[i,config.weight_derivative_combinations.index(('ctZ',))]#/ test_weights[i,0]
#                true_quad_inter =test_weights[i,config.weight_derivative_combinations.index(('ctZ','ctZ'))]#/ test_weights[i,0] 
#                index_var = 0
#                if test_features[i,17] < 160:
#                        index_var = 0
#                elif test_features[i,17] < 170:
#                        index_var = 1
#                elif test_features[i,17] < 180:
#                        index_var = 2
#                elif test_features[i,17] < 190:
#                        index_var = 3
#                else:
#                        index_var = 4
#                bin_150_200_pred_lin[index_var] += pred_lin_inter
#                bin_150_200_pred_quad[index_var] += pred_quad_inter
#                bin_150_200_true_lin[index_var] += true_lin_inter
#                bin_150_200_true_quad[index_var] += true_quad_inter
#                bin_150_200_weight_0[index_var] += test_weights[i,0]

#print "##################################################################################################"
#print "##################################################################################################"
#print "Sums in 5 bins from 150 to 200 in ptGamma"
#print "Linear Pred: ", bin_150_200_pred_lin/bin_150_200_weight_0, " Linear True: ", bin_150_200_true_lin/bin_150_200_weight_0, " Quad Pred: ", bin_150_200_pred_quad/bin_150_200_weight_0, " Quad True: ", bin_150_200_true_quad/bin_150_200_weight_0

bin_pred_lin = np.zeros(40)
bin_true_lin = np.zeros(40) 
bin_pred_quad = np.zeros(40) 
bin_true_quad = np.zeros(40) 
bin_weight_0 =np.zeros(40)  
for i in range(test_features[:,17].size):
                pred_lin_inter =bit_predictions[('ctZ',)][i]*test_weights[i,0]
                pred_quad_inter =bit_predictions[('ctZ','ctZ')][i]*test_weights[i,0]
                true_lin_inter =test_weights[i,config.weight_derivative_combinations.index(('ctZ',))]#/ test_weights[i,0]
                true_quad_inter =test_weights[i,config.weight_derivative_combinations.index(('ctZ','ctZ'))]#/ test_weights[i,0] 
                index_var = 0
                for j in range(40):
                        if test_features[i,17] > j*10 and test_features[i,17] <= (j+1)*10:
                                bin_pred_lin[j] += pred_lin_inter
                                bin_pred_quad[j] += pred_quad_inter
                                bin_true_lin[j] += true_lin_inter
                                bin_true_quad[j] += true_quad_inter
                                bin_weight_0[j] += test_weights[i,0]

print "##################################################################################################"
print "##################################################################################################"
print "Sums in ptGamma bins (10, from 0 to 400)"
print "Linear Pred: ", bin_pred_lin/bin_weight_0, " Linear True: ", bin_true_lin/bin_weight_0, " Quad Pred: ", bin_pred_quad/bin_weight_0, " Quad True: ", bin_true_quad/bin_weight_0

for i in range(bin_weight_0.size):
        print "#################################"
        print "Bin: from ", i*10, " to ", (i+1)*10
        print "Linear Pred: ", bin_pred_lin[i]/bin_weight_0[i], " Linear True: ", bin_true_lin[i]/bin_weight_0[i], " Quad Pred: ", bin_pred_quad[i]/bin_weight_0[i], " Quad True: ", bin_true_quad[i]/bin_weight_0[i] 

plot_directory = "/mnt/hephy/cms/stefan.rohshap/www/TMB/"
directory = os.path.join(plot_directory, 'Plot_LLR',args.name_dir)
if not os.path.exists(directory):
        try:
            os.makedirs(directory)
        except IOError:
            pass

c1 = ROOT.TCanvas()
h_lin_pred = ROOT.TH1D("h_lin_pred","Linear Prediction",40,0,400)
h_quad_pred = ROOT.TH1D("h_quad_pred","Quadratic Prediction",40,0,400)
h_lin_true = ROOT.TH1D("h_lin_true","Linear True",40,0,400)
h_quad_true = ROOT.TH1D("h_quad_true","Quadratic True",40,0,400)
h_weight_0 = ROOT.TH1D("h_quad_true","Quadratic True",40,0,400)
for i in range(test_features[:,17].size):
        pred_lin_inter =bit_predictions[('ctZ',)][i]*test_weights[i,0]
        pred_quad_inter =bit_predictions[('ctZ','ctZ')][i]*test_weights[i,0]
        true_lin_inter =test_weights[i,config.weight_derivative_combinations.index(('ctZ',))]#/ test_weights[i,0]
        true_quad_inter =test_weights[i,config.weight_derivative_combinations.index(('ctZ','ctZ'))]#/ test_weights[i,0] 
        h_lin_pred.Fill(test_features[i,17],pred_lin_inter)
        h_quad_pred.Fill(test_features[i,17],pred_quad_inter)
        h_lin_true.Fill(test_features[i,17],true_lin_inter)
        h_quad_true.Fill(test_features[i,17],true_quad_inter)
        h_weight_0.Fill(test_features[i,17],test_weights[i,0])

#for i in range(bin_weight_0.size):
#        h_lin_pred.Fill(test_features[i,17],bin_pred_lin[i]/bin_weight_0[i])
#        h_quad_pred.Fill(test_features[i,17],bin_pred_quad[i]/bin_weight_0[i])
#        h_lin_true.Fill(test_features[i,17],bin_true_lin[i]/bin_weight_0[i])
#        h_quad_true.Fill(test_features[i,17],bin_true_quad[i]/bin_weight_0[i])


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
