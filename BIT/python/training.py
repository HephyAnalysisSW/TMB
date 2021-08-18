#!/usr/bin/env python

import ROOT, os
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

args = argParser.parse_args()

input_postfix = ''
if args.small:
    args.name    +="_small"
    input_postfix = '_small'

#Logger
import Analysis.Tools.logger as logger
logger = logger.get_logger("INFO", logFile = None )

# MVA configuration
import importlib
configs = importlib.import_module(args.config_module)
config  = getattr( configs, args.config)

import uproot
import awkward
import numpy as np
import pandas as pd
#import h5py

#########################################################################################
# variable definitions

import Analysis.Tools.user as user

# directories
plot_directory   = os.path.join( user. plot_directory, 'MVA', args.name, args.config )
output_directory = os.path.join( args.output_directory, 'models', args.name, args.config)
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
for i_training_sample, training_sample in enumerate(config.training_samples):
    upfile_name = os.path.join(os.path.expandvars(args.input_directory), args.config+input_postfix, training_sample.name, training_sample.name+'.root')
    logger.info( "Loading upfile %i: %s from %s", i_training_sample, training_sample.name, upfile_name)
    upfile = uproot.open(upfile_name)
    features[training_sample.name]  = upfile["Events"].pandas.df(branches = mva_variables )
    weight_derivatives[training_sample.name]  = upfile["Events"].pandas.df(branches = ["weight_derivatives"] )

features = pd.concat([features[training_sample.name] for training_sample in config.training_samples])
features = features.values

weight_derivatives = pd.concat([weight_derivatives[training_sample.name] for training_sample in config.training_samples])
weight_derivatives = weight_derivatives.values.reshape((len(features),-1))

# split dataset into Input and output data

# number of samples with 'small'
n_small_samples = 10000

# small
if args.small:
    features = features[:n_small_samples]
    weight_derivatives = weight_derivatives[:n_small_samples]

features  = features[:,0:n_var_flat]

## split data into train and test, test_size = 0.2 is quite standard for this
#from sklearn.model_selection import train_test_split
#
#options = {'test_size':0.2, 'random_state':7, 'shuffle':True}
#X_train, X_test, Y_train, Y_test                  = train_test_split(X, Y, **options)
#validation_data = ( X_test,  Y_test)
#training_data   =   X_train

# Boosting
import sys, os, time
sys.path.insert(0,os.path.expandvars("$CMSSW_BASE/src/BIT"))
from BoostedInformationTree import BoostedInformationTree 

bits = {}
for derivative in config.bit_derivatives:
    if derivative == tuple(): continue

    filename = os.path.join(output_directory, "bit_derivative_%s"% ('_'.join(derivative))) + '.pkl'
    try:
        print ("Loading %s for %r"%( filename, derivative))
        bits[derivative] = BoostedInformationTree.load(filename)
    except IOError:
        time1 = time.time()
        print ("Learning %s"%( str(derivative)))
        bits[derivative]= BoostedInformationTree(
                training_features     = features,
                training_weights      = weight_derivatives[:,0],
                training_diff_weights = weight_derivatives[:, config.weight_derivative_combinations.index(derivative)],
                split_method          = 'vectorized_split_and_weight_sums',
                weights_update_method = 'vectorized',
                calibrated            = False,
                **config.bit_cfg
                    )
        bits[derivative].boost()
        bits[derivative].save(filename)
        print ("Written %s"%( filename ))

        time2 = time.time()
        boosting_time = time2 - time1
        print ("Boosting time: %.2f seconds" % boosting_time)

        ## plot loss
        #test_scores     = bits[derivative].vectorized_predict(test_features)
        #training_scores = bits[derivative].vectorized_predict(test_features)
        #max_score = max(test_scores)
        #min_score = min(test_scores)

        #test_FIs            = np.zeros(n_trees)
        #training_FIs        = np.zeros(n_trees)

        #for i in range(args.nTraining):
        #    test_scores     = bits[derivative].predict( test_features[i],     summed = False)
        #    training_scores = bits[derivative].predict( training_features[i], summed = False)

        #    test_score  = sum( test_scores )
        #    train_score = sum( training_scores )

        #    # compute test and training FI evolution during training
        #    test_FIs     += test_weights[derivative][i]*test_scores
        #    training_FIs += training_weights[derivative][i]*training_scores

        #training_FI_histo     = ROOT.TH1D("trainFI", "trainFI",          n_trees, 1, n_trees+1 )
        #test_FI_histo         = ROOT.TH1D("testFI",  "testFI",           n_trees, 1, n_trees+1 )

        #for i_tree in range(n_trees):
        #    test_FI_histo    .SetBinContent( i_tree+1, -sum(test_FIs[:i_tree]) )
        #    training_FI_histo.SetBinContent( i_tree+1, -sum(training_FIs[:i_tree]) )

        ## Histo style
        #test_FI_histo    .style = styles.lineStyle( ROOT.kBlue, width=2 )
        #training_FI_histo.style = styles.lineStyle( ROOT.kRed, width=2 )
        #test_FI_histo    .legendText = "Test"
        #training_FI_histo.legendText = "Training"

        #minY   = 0.01 * min( test_FI_histo.GetBinContent(test_FI_histo.GetMaximumBin()), training_FI_histo.GetBinContent(training_FI_histo.GetMaximumBin()))
        #maxY   = 1.5  * max( test_FI_histo.GetBinContent(test_FI_histo.GetMaximumBin()), training_FI_histo.GetBinContent(training_FI_histo.GetMaximumBin()))

        #histos = [ [test_FI_histo], [training_FI_histo] ]
        #plot   = Plot.fromHisto( filename+"_evolution", histos, texX="b", texY="L(D,b)" )

        ## Plot Style
        #histModifications      = []
        #histModifications      += [ lambda h: h.GetYaxis().SetTitleOffset(1.4) ]
        #histModifications += [ lambda h: h.GetXaxis().SetTitleSize(26) ]
        #histModifications += [ lambda h: h.GetYaxis().SetTitleSize(26) ]
        #histModifications += [ lambda h: h.GetXaxis().SetLabelSize(22)  ]
        #histModifications += [ lambda h: h.GetYaxis().SetLabelSize(22)  ]

        #ratioHistModifications = []
        #ratio                  = None
        #legend                 = (0.6,0.75,0.9,0.88)
        #yRange                 = "auto" #( minY, maxY )

        #plot1DHist( plot, plot_directory, yRange=yRange, ratio=ratio, legend=legend, plotLog=False, titleOffset=0.08, histModifications=histModifications )

