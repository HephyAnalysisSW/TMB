import ROOT, os
from RootTools.core.standard import *
import Analysis.Tools.syncer as syncer
c1 = ROOT.TCanvas() # do this to avoid version conflict in png.h with keras import ...
c1.Draw()
c1.Print('/tmp/delete.png')

import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--config',             action='store', type=str,   default='tttt_3l', help="Name of the config file")
argParser.add_argument('--name',               action='store', type=str,   default='default', help="Name of the training")
argParser.add_argument('--variable_set',       action='store', type=str,   default='mva_variables', help="List of variables for training")
argParser.add_argument('--output_directory',   action='store', type=str,   default='/mnt/hephy/cms/robert.schoefbeck/TMB/models/')
argParser.add_argument('--input_directory',    action='store', type=str,   default=os.path.expandvars("/eos/vbc/user/$USER/TMB/training-ntuples-tttt-v1/MVA-training/") )
argParser.add_argument('--small',              action='store_true', help="small?")

args = argParser.parse_args()

if args.small: args.name+="_small"

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
plot_directory   = os.path.join( user. plot_directory, 'MVA', args.name, args.config )
output_directory = os.path.join( args.output_directory, args.name, args.config) 

# fix random seed for reproducibility
np.random.seed(1)

# get the training variable names
mva_variables = [ mva_variable[0] for mva_variable in getattr(config, args.variable_set) ]

n_var_input   = len(mva_variables)

df_file = {}
for i_training_sample, training_sample in enumerate(config.training_samples):
    upfile_name = os.path.join(os.path.expandvars(args.input_directory), args.config, training_sample.name, training_sample.name+'.root')
    logger.info( "Loading upfile %i: %s from %s", i_training_sample, training_sample.name, upfile_name)
    upfile = uproot.open(os.path.join(os.path.expandvars(args.input_directory), args.config, training_sample.name, training_sample.name+'.root'))
    df_file[training_sample.name]  = upfile["Events"].pandas.df(branches = mva_variables )
    # enumerate
    df_file[training_sample.name]['signal_type'] =  np.ones(len(df_file[training_sample.name])) * i_training_sample

df = pd.concat([df_file[training_sample.name] for training_sample in config.training_samples])

#df = df.dropna() # removes all Events with nan -> amounts to M3 cut

# split dataset into Input and output data
dataset = df.values

# small
if args.small:
    dataset = dataset[:10000]

X  = dataset[:,0:n_var_input]

# regress FI
Y = dataset[:, n_var_input] 

from sklearn.preprocessing import label_binarize
classes = range(len(config.training_samples))
Y = label_binarize(Y, classes=classes)

# split data into train and test, test_size = 0.2 is quite standard for this
from sklearn.model_selection import train_test_split

X_train, X_test, Y_train, Y_test                  = train_test_split(X, Y,    test_size=0.2, random_state=7, shuffle = True)

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
                    Dense(len(config.training_samples), kernel_initializer='normal', activation='sigmoid'),
                    ])

#model.compile(optimizer='adam', loss='binary_crossentropy', metrics=['accuracy'])
model.compile(optimizer='adam', loss='mean_squared_error', metrics=['mean_absolute_percentage_error'])
model.summary()

# define callback for early stopping
import tensorflow as tf
callback = tf.keras.callbacks.EarlyStopping(monitor='loss', patience=3) # patience can be higher if a more accurate result is preferred
                                                                        # I would recommmend at least 3, otherwise it might cancel too early

# train the model
batch_size = 1024*6
history = model.fit(X_train, 
                    Y_train, 
                    sample_weight = None,
                    epochs=100, 
                    batch_size=batch_size,
                    #verbose=0, # switch to 1 for more verbosity, 'silences' the output
                    callbacks=[callback],
                    #validation_split=0.1
                    validation_data= ( (X_test,Y_test) ),
                   )
print('training finished')

# saving
if not os.path.exists(output_directory):
    os.makedirs(output_directory)

output_file = os.path.join(output_directory, 'regression_model.h5')
model.save(output_file)
logger.info("Written model to: %s", output_file)
