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
argParser.add_argument('--output_directory',   action='store', type=str,   default='/mnt/hephy/cms/rosmarie.schoefbeck/TMB/models/')
argParser.add_argument('--input_directory',    action='store', type=str,   default=os.path.expandvars("/eos/vbc/user/rosmarie.schoefbeck/TMB/training-ntuples-tttt-v1/MVA-training") )
argParser.add_argument('--small',              action='store_true', help="small?")
argParser.add_argument('--add_LSTM',           action='store_true', help="add LSTM?")

args = argParser.parse_args()

if args.add_LSTM:   args.name+="_LSTM"
if args.small:      args.name+="_small"

#Logger
import TMB.Tools.logger as logger
logger = logger.get_logger("INFO", logFile = None )

# MVA configuration
import TMB.MVA.configs  as configs 

#config
config = getattr( configs, args.config)

import uproot
import awkward
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

n_var_flat   = len(mva_variables)

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

# number of samples with 'small'
n_small_samples = 10000

# small
if args.small:
    dataset = dataset[:n_small_samples]

X  = dataset[:,0:n_var_flat]

# regress FI
Y = dataset[:, n_var_flat] 

from sklearn.preprocessing import label_binarize
classes = range(len(config.training_samples))
Y = label_binarize(Y, classes=classes)

# loading vector branches for LSTM
if args.add_LSTM:
    vector_branches = ["mva_JetGood_%s" % varname for varname in config.jetVarNames]
    max_timestep = 10 # for LSTM

    vec_br_f  = {}

    for i_training_sample, training_sample in enumerate(config.training_samples):
        upfile_name = os.path.join(os.path.expandvars(args.input_directory), args.config, training_sample.name, training_sample.name+'.root')
        logger.info( "Loading vector branches %i: %s from %s", i_training_sample, training_sample.name, upfile_name)
        with uproot.open(upfile_name) as upfile:
            vec_br_f[i_training_sample]   = {}
            for name, branch in upfile["Events"].arrays(vector_branches).iteritems():
                vec_br_f[i_training_sample][name] = branch.pad(max_timestep)[:,:max_timestep].fillna(0)

    vec_br = {name: awkward.JaggedArray.concatenate( [vec_br_f[i_training_sample][name] for i_training_sample in range(len(config.training_samples))] ) for name in vector_branches}
    if args.small:
        for key, branch in vec_br.iteritems():
            vec_br[key] = branch[:n_small_samples]

    # put columns side by side and transpose the innermost two axis
    len_samples = len(vec_br.values()[0])
    V           = np.column_stack( [vec_br[name] for name in vector_branches] ).reshape( len_samples, len(vector_branches), max_timestep).transpose((0,2,1))

# split data into train and test, test_size = 0.2 is quite standard for this
from sklearn.model_selection import train_test_split

options = {'test_size':0.2, 'random_state':7, 'shuffle':True}
if args.add_LSTM:
    X_train, X_test, Y_train, Y_test, V_train, V_test = train_test_split(X, Y, V, **options)
    validation_data = ( [X_test,  V_test], Y_test )
    training_data   =   [X_train, V_train]
else:
    X_train, X_test, Y_train, Y_test                  = train_test_split(X, Y, **options)
    validation_data = ( X_test,  Y_test)
    training_data   =   X_train

#########################################################################################
# define model (neural network)
from keras.models import Sequential, Model
from keras.optimizers import SGD
from keras.layers import Input, Activation, Dense, Convolution2D, MaxPooling2D, Dropout, Flatten, LSTM, Concatenate
from keras.layers import BatchNormalization
from keras.utils import np_utils

#model = Sequential([
#                    #Flatten(input_shape=(n_var_flat, 1)), # instead of (n_var_flat,1)
#                    BatchNormalization(input_shape=(n_var_flat, )),
#                    Dense(n_var_flat*2, activation='sigmoid'),
#                    Dense(n_var_flat+5, activation='sigmoid'),
#                    #Dense(n_var_flat*5, activation='sigmoid'),
#                    Dense(len(config.training_samples), kernel_initializer='normal', activation='sigmoid'),
#                    ])

# flat layers
flat_inputs = Input(shape=(n_var_flat, ))
x = BatchNormalization(input_shape=(n_var_flat, ))(flat_inputs)
x = Dense(n_var_flat*2, activation='sigmoid')(x)
x = Dense(n_var_flat+5, activation='sigmoid')(x)

inputs = flat_inputs

# LSTMs
if args.add_LSTM:
    vec_inputs = Input(shape=(max_timestep, len(vector_branches),) )
    v = LSTM(10, activation='relu', input_shape=(max_timestep, len(vector_branches)))( vec_inputs )
    x = Concatenate()( [x, v])
    inputs = ( flat_inputs, vec_inputs)

outputs = Dense(len(config.training_samples), kernel_initializer='normal', activation='sigmoid')(x)
model = Model( inputs, outputs )

#model.compile(optimizer='adam', loss='binary_crossentropy', metrics=['accuracy'])
model.compile(optimizer='adam', loss='mean_squared_error', metrics=['mean_absolute_percentage_error'])
model.summary()

# define callback for early stopping
import tensorflow as tf
callback = tf.keras.callbacks.EarlyStopping(monitor='loss', patience=3) # patience can be higher if a more accurate result is preferred
                                                                        # I would recommmend at least 3, otherwise it might cancel too early

# train the model
batch_size = 1024*6
history = model.fit(training_data, 
                    Y_train, 
                    sample_weight = None,
                    epochs=100, 
                    batch_size=batch_size,
                    #verbose=0, # switch to 1 for more verbosity, 'silences' the output
                    callbacks=[callback],
                    #validation_split=0.1
                    validation_data = validation_data,
                   )
print('training finished')

# saving
if not os.path.exists(output_directory):
    os.makedirs(output_directory)

output_file = os.path.join(output_directory, 'regression_model.h5')
model.save(output_file)
logger.info("Written model to: %s", output_file)

#plot roc curves 
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc

Y_predict = model.predict(X_test)
#if args.add_LSTM:
#    Y_predict = model.predict( [X_test,  V_test] )
fpr = dict()
tpr = dict()
roc_auc = dict()
for i in classes:
    fpr[i], tpr[i], _ = roc_curve(Y_test[:, i], Y_predict[:, i])
    roc_auc[i] = auc(fpr[i], tpr[i])

#from itertools import cycle
#
#colors = cycle(['aqua', 'darkorange', 'cornflowerblue', 'indigo', 'darkgreen', 'crimson',
#    'sienna', 'darkmagenta', 'darkslatergrey', 'maroon', 'olive', 'purple'])
#
#
#for i, color in zip(classes, colors):
#    ra = int(roc_auc[i] * 100)
#    ra = ra / 100
#    plt.plot(fpr[i], tpr[i], color=color, lw=2,
#             label=key_list[i] + ' (area = ' + str(ra) +')')
#
lw = 2
plt.plot(fpr[0], tpr[0], color='darkorange', lw=lw, label='ROC curve (area = %0.2f)' %roc_auc[0])
plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
#plt.plot([0, 1], [0, 1], 'k--', lw=2)
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('ROC_curve')
plt.legend(loc="lower right")
plt.savefig(plot_directory + 'ROC_curve.png')

#signal vs all 
Y_predict_s = np.zeros((len(Y_predict), 2))
Y_predict_s[:,0] = Y_predict[:,0]
Y_predict_s[:,1] = Y_predict[:,1] + Y_predict[:,2] + Y_predict[:,3]
for i, y in enumerate(Y_predict_s):
    Y_predict_s[i] = Y_predict_s[i] / np.sum(y)

fpr, tpr, _ = roc_curve(Y_test[:, 0], Y_predict_s[:, 0])
plt.clf()
plt.plot([0, 1], [0, 1], 'k--')
plt.plot(fpr, tpr)
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('ROC_signalvsall')
plt.legend(loc="lower right")
plt.savefig(plot_directory+'ROC_s.png')
