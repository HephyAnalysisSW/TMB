import ROOT, os
from RootTools.core.standard import *
import Analysis.Tools.syncer as syncer
c1 = ROOT.TCanvas() # do this to avoid version conflict in png.h with keras import ...
c1.Draw()
c1.Print('/tmp/delete.png')

import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--trainingfiles',      action='store', type=str,   nargs = '*', default=['input.root'], help="Input files for training")
argParser.add_argument('--FI_branch',          action='store', type=str,   default='FI_ctZ_BSM', help="Regression target")
argParser.add_argument('--config',             action='store', type=str,   default='ttG_WG', help="Name of the config file")
argParser.add_argument('--name',               action='store', type=str,   default='default', help="Name of the training")
argParser.add_argument('--variable_set',       action='store', type=str,   default='mva_variables', help="List of variables for training")
argParser.add_argument('--output_directory',   action='store', type=str,   default='/mnt/hephy/cms/robert.schoefbeck/TMB/models/')
argParser.add_argument('--truth_input',        action='store_true', help="Include truth in training input")

args = argParser.parse_args()

#Logger
import TMB.Tools.logger as logger
logger = logger.get_logger("INFO", logFile = None )

# MVA configuration
import TMB.MVA.configs  as configs 

#config
config = getattr( configs, args.config)

import uproot
import numpy as np
import pandas as pd
#import h5py

#########################################################################################
# variable definitions

import TMB.Tools.user as user 

# directories
plot_directory   = os.path.join( user. plot_directory, 'MVA', args.name, args.config, args.FI_branch )
output_directory = os.path.join( args.output_directory, args.name, args.config, args.FI_branch) 

# input samples
import TMB.Samples.pp_TTGammaEFT as samples
sample = samples.ttG_noFullyHad_fast

# fix random seed for reproducibility
np.random.seed(1)

# get the training variable names
mva_variables = [ mva_variable[0] for mva_variable in getattr(config, args.variable_set) ]

n_var_input   = len(mva_variables)

df_file = {}
for trainingfile in args.trainingfiles:
    upfile = uproot.open(os.path.join("/eos/vbc/user/robert.schoefbeck/TMB/MVA-training/", args.config, trainingfile))
    df_file[trainingfile]  = upfile["Events"].pandas.df(branches = mva_variables+[args.FI_branch])

df = pd.concat([df_file[trainingfile] for trainingfile in args.trainingfiles])

#branches = upfile['Events'].arrays(namedecode='utf-8')
#basic    = (branches['mva_m3'] >= 0 ) 

df = df.dropna() # removes all Events with nan -> amounts to M3 cut

# split dataset into Input and output data
dataset = df.values
X = dataset[:,0:n_var_input]
#Y = np.log(dataset[:,n_var_input])
log_FI = np.log10(dataset[:,n_var_input])
# cap negative infinity
min_log_FI = -30

log_FI[log_FI<min_log_FI] = min_log_FI

#Y[Y < -15] = -15
#quantiles = [.2, .5, .8, .9, .95]
#quantiles = [.9, .95]
quantiles = [.05, .10, .15, .20, .25, .30, .35, .40, .45, .50, .55, .60, .65, .70, .75, .80, .85, .90, .95]
# compute quantiles for values above the minimum
q=np.quantile(log_FI[log_FI>min_log_FI], quantiles)
q=np.concatenate(([-np.inf], q, [np.inf]))

Y = np.digitize(log_FI, q)

## train with truth!! FIXME!!

if args.truth_input:
    logger.waring("Training with truth input!")
    X = np.concatenate( (X, Y.reshape(len(Y),1)), axis=1)
    n_var_input+=1

# split data into train and test, test_size = 0.2 is quite standard for this
from sklearn.model_selection import train_test_split
X_train_val, X_test, Y_train_val, Y_test = train_test_split(X, Y, test_size=0.2, random_state=7, shuffle = True)

#assert False, ""

#########################################################################################
# define model (neural network)
from keras.models import Sequential, Model
from keras.optimizers import SGD
from keras.layers import Input, Activation, Dense, Convolution2D, MaxPooling2D, Dropout, Flatten
from keras.layers import BatchNormalization
from keras.utils import np_utils

model = Sequential([
                    #Flatten(input_shape=(n_var_input, 1)), # instead of (n_var_input,1)
                    BatchNormalization(input_shape=(n_var_input, )),
                    Dense(n_var_input*2, activation='sigmoid'),
                    Dense(n_var_input+5, activation='sigmoid'),
                    #Dense(n_var_input*5, activation='sigmoid'),
                    Dense(1),#,kernel_initializer='normal')
                    ])

#model.compile(optimizer='adam', loss='binary_crossentropy', metrics=['accuracy'])
model.compile(optimizer='adam', loss='mean_squared_error', metrics=['mean_absolute_percentage_error'])
model.summary()

# define callback for early stopping
import tensorflow as tf
callback = tf.keras.callbacks.EarlyStopping(monitor='loss', patience=5) # patience can be higher if a more accurate result is preferred
                                                                        # I would recommmend at least 3, otherwise it might cancel too early

# train the model
batch_size = 1024*6
history = model.fit(X_train_val, 
                    Y_train_val, 
                    epochs=500, 
                    batch_size=batch_size,
                    #verbose=0, # switch to 1 for more verbosity, 'silences' the output
                    callbacks=[callback],
                    #validation_split=0.1
                    validation_data=(X_test,Y_test) # use either validation_split or validation_data
                   )
print('training finished')

# saving
if not os.path.exists(output_directory):
    os.makedirs(output_directory)

output_file = os.path.join(output_directory, 'regression_model.h5')
model.save(output_file)
logger.info("Written model to: %s", output_file)

#########################################################################################
# Apply the model

import array
c = ROOT.TCanvas()

n = len(quantiles)+1
h = ROOT.TH2F("scatter", "scatter", n,.5,n+.5,n,.5,n+.5)
pred_test = list(model.predict(X_test))
for x,y in zip( list(Y_test), pred_test):
    h.Fill( x,y )


plot2D = Plot2D.fromHisto(name = "scatter", histos = [[h]], texX = "truth", texY = "prediction" )
plotting.draw2D(plot2D, plot_directory = plot_directory, logY = False, logX = False, logZ = True, copyIndexPHP=True)

for var, binning in [ 
        ('mva_photon_pt',   [30,0,300]),
        ('mva_photonJetdR', [30,0,5]),
        ('mva_photonLepdR', [30,0,5]),
        ('mva_photon_eta',  [30,-2.5,2.5]),
        ]:

    if var not in mva_variables: continue

    h = ROOT.TH2F("truth_"+var, "truth_"+var, n,.5,n+.5, *binning)
    for x,y in zip( list(Y_test), list(X_test[:,mva_variables.index(var)]) ):
        h.Fill( x,y )

    plot2D = Plot2D.fromHisto(name = "truth_vs_"+var, histos = [[h]], texX = "truth", texY = var )
    plotting.draw2D(plot2D, plot_directory = plot_directory, logY = False, logX = False, logZ = True)

    h = ROOT.TH2F("pred_"+var, "pred_"+var, n,.5,n+.5, *binning)
    for x,y in zip( pred_test, list(X_test[:,mva_variables.index(var)]) ):
        h.Fill( x,y )

    plot2D = Plot2D.fromHisto(name = "pred_vs_"+var, histos = [[h]], texX = "pred", texY = var )
    plotting.draw2D(plot2D, plot_directory = plot_directory, logY = False, logX = False, logZ = True)

syncer.sync()
