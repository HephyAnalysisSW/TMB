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
argParser.add_argument('--weighted_training',  action='store_true', help="Train weighted?")
argParser.add_argument('--add_weighted_quantiles', action='store_true', help="Make quantiles from weighted FI?")
argParser.add_argument('--small',              action='store_true', help="add LSTM?")
argParser.add_argument('--add_LSTM',           action='store_true', help="add LSTM?")

args = argParser.parse_args()

if args.truth_input:        args.name+="_truth"
if args.weighted_training:  args.name+="_wtr"
if args.add_weighted_quantiles: args.name+="_wq"
if args.add_LSTM:           args.name+="_LSTM"
if args.small:              args.name+="_small"

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
plot_directory   = os.path.join( user. plot_directory, 'MVA', args.name, args.config, args.FI_branch )
output_directory = os.path.join( args.output_directory, args.name, args.config, args.FI_branch) 

# fix random seed for reproducibility
np.random.seed(1)

# number of samples with 'small'
n_small_samples = 10000

# load flat input 
mva_variables = [ mva_variable[0] for mva_variable in getattr(config, args.variable_set) ]
n_var_flat   = len(mva_variables)
df_file = {}
for trainingfile in args.trainingfiles:
    with uproot.open(os.path.join("/eos/vbc/user/robert.schoefbeck/TMB/training-ntuples-v6/MVA-training/", args.config, trainingfile)) as upfile:
        df_file[trainingfile]  = upfile["Events"].pandas.df(branches = mva_variables+[args.FI_branch])
df = pd.concat([df_file[trainingfile] for trainingfile in args.trainingfiles])
df = df.dropna() # removes all Events with nan -> amounts to M3 cut
dataset = df.values
if args.small:
    dataset = dataset[:n_small_samples]

X  = dataset[:,0:n_var_flat]

# deal with regression target
FI     = dataset[:,n_var_flat]
# make positive
FI[FI<0] = 0

log_FI = np.log10(dataset[:,n_var_flat])

min_log_FI = -30
log_FI[np.isnan(log_FI)]=-float('inf')
log_FI[log_FI<min_log_FI] = min_log_FI

quantiles = np.array([ i/100. for i in range(1,100)] )

# compute quantiles for values above the minimum
q=np.quantile(log_FI[log_FI>min_log_FI], quantiles)

# add FI weighted quantils
if args.add_weighted_quantiles:
    sorter = np.argsort(log_FI)
    log_FI_sorted = log_FI[sorter]
    #sample_weight = np.ones(FI.size)
    FI_weight =     FI[sorter]

    weighted_quantiles = np.cumsum(FI_weight) - 0.5 * FI_weight #https://stackoverflow.com/questions/21844024/weighted-percentile-using-numpy
    weighted_quantiles /= np.sum(FI_weight)

    weighted_quantiles = np.interp(quantiles, weighted_quantiles, log_FI_sorted)

    q = np.concatenate( (q[q<weighted_quantiles.min()], weighted_quantiles[1:] ))

q = np.unique(q)
q.sort()

logger.info("Using these quantiles: %r", list(q))

# regress FI
q = np.concatenate(([-np.inf], q, [np.inf]))
Y_dig = np.digitize(log_FI, q)
Y = (Y_dig-Y_dig.min())/float(Y_dig.max()-Y_dig.min())

if args.truth_input:
    logger.waring("Training with truth input!")
    X= np.concatenate( (X, Y.reshape(len(Y),1)), axis=1)
    n_var_flat+=1

# loading vector branches for LSTM
if args.add_LSTM:
    vector_branches = ["mva_JetGood_%s" % varname for varname in config.jetVarNames]
    max_timestep = 10 # for LSTM

    vec_br_f  = {}
    for trainingfile in args.trainingfiles:
        with uproot.open(os.path.join("/eos/vbc/user/robert.schoefbeck/TMB/training-ntuples-v6/MVA-training/", args.config, trainingfile)) as upfile:
            vec_br_f[trainingfile]   = {}
            for name, branch in upfile["Events"].arrays(vector_branches).iteritems():
                vec_br_f[trainingfile][name] = branch.pad(max_timestep)[:,:max_timestep].fillna(0)

    vec_br = {name: awkward.JaggedArray.concatenate( [vec_br_f[trainingfile][name] for trainingfile in args.trainingfiles] ) for name in vector_branches}
    if args.small:
        for key, branch in vec_br.iteritems():
            vec_br[key] = branch[:n_small_samples]

    # put columns side by side and transpose the innermost two axis
    len_samples = len(vec_br.values()[0])
    V           = np.column_stack( [vec_br[name] for name in vector_branches] ).reshape( len_samples, len(vector_branches), max_timestep).transpose((0,2,1))
# split data into train and test, test_size = 0.2 is quite standard for this
from sklearn.model_selection import train_test_split

# now compute weights

options = {'test_size':0.2, 'random_state':7, 'shuffle':True}
if args.add_LSTM:
    if args.weighted_training:
        counts  = np.unique(Y,return_counts=True)[1]
        weights = float(counts.min())/counts
        W       = np.array(map(weights.__getitem__, Y_dig-1))
        X_train, X_test, Y_train, Y_test, V_train, V_test, W_train, W_test = train_test_split(X, Y, V, W, **options)

        validation_data = ( [X_test,  V_test],  Y_test,  W_test )
        training_data   =   [X_train, V_train]
    else:
        X_train, X_test, Y_train, Y_test, V_train, V_test = train_test_split(X, Y, V, **options)

        W_train = None        
        validation_data = ( [X_test,  V_test], Y_test )
        training_data   =   [X_train, V_train]
else:
    if args.weighted_training:
        counts  = np.unique(Y,return_counts=True)[1]
        weights = float(counts.min())/counts
        W       = np.array(map(weights.__getitem__, Y_dig-1))
        X_train, X_test, Y_train, Y_test, W_train, W_test = train_test_split(X, Y, W, **options)

        validation_data = ( X_test,  Y_test,  W_test )
        training_data   =   X_train

    else:
        X_train, X_test, Y_train, Y_test = train_test_split(X, Y, **options)

        W_train = None        
        validation_data = ( X_test,  Y_test)
        training_data   =   X_train

# define model (neural network)
from keras.models import Sequential, Model
from keras.optimizers import SGD
from keras.layers import Input, Activation, Dense, Convolution2D, MaxPooling2D, Dropout, Flatten, LSTM, Concatenate
from keras.layers import BatchNormalization
from keras.utils import np_utils

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

outputs = Dense(1)(x)

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
history = model.fit( 
                    training_data,
                    Y_train, 
                    sample_weight = W_train,
                    epochs = 500, 
                    batch_size = batch_size,
                    verbose=0, # switch to 1 for more verbosity, 'silences' the output
                    callbacks = [callback],
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

#########################################################################################
# Apply the model

import array
c = ROOT.TCanvas()

h = ROOT.TH2F("scatter", "scatter", 100,0,1,100,0,1)
if args.add_LSTM:
    pred_test = list(model.predict([X_test,  V_test]))
else:
    pred_test = list(model.predict(X_test))

for x,y in zip( list(Y_test), pred_test):
    h.Fill( x,y )

plot2D = Plot2D.fromHisto(name = "scatter", histos = [[h]], texX = "truth", texY = "prediction" )
plotting.draw2D(plot2D, plot_directory = plot_directory, logY = False, logX = False, logZ = True, copyIndexPHP=True)

for var, binning in [ 
        ('mva_photon_pt',   [25,0,500]),
        ('mva_photonJetdR', [30,0,5]),
        ('mva_photonLepdR', [30,0,5]),
        ('mva_photon_eta',  [30,-2.5,2.5]),
        ]:

    if var not in mva_variables: continue

    h = ROOT.TH2F("truth_"+var, "truth_"+var, 100,0,1, *binning)
    for x,y in zip( list(Y_test), list(X_test[:,mva_variables.index(var)]) ):
        h.Fill( x,y )

    plot2D = Plot2D.fromHisto(name = "truth_vs_"+var, histos = [[h]], texX = "truth", texY = var )
    plotting.draw2D(plot2D, plot_directory = plot_directory, logY = False, logX = False, logZ = True)

    h = ROOT.TH2F("pred_"+var, "pred_"+var, 100,0,1, *binning)
    for x,y in zip( pred_test, list(X_test[:,mva_variables.index(var)]) ):
        h.Fill( x,y )

    plot2D = Plot2D.fromHisto(name = "pred_vs_"+var, histos = [[h]], texX = "pred", texY = var )
    plotting.draw2D(plot2D, plot_directory = plot_directory, logY = False, logX = False, logZ = True)

syncer.sync()
