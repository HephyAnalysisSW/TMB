import ROOT
import Analysis.Tools.syncer as syncer
c1 = ROOT.TCanvas() # do this to avoid version conflict in png.h with keras import ...
c1.Draw()
c1.Print('/tmp/delete.png')

import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
#argParser.add_argument('--sample',             action='store', type=str,   default='ttG_noFullyHad_fast')
argParser.add_argument('--config',             action='store', type=str,   default='ttG')
#argParser.add_argument('--output_directory',   action='store', type=str,   default='.')
#argParser.add_argument('--small',              action='store_true')

args = argParser.parse_args()

#Logger
import TMB.Tools.logger as logger
logger = logger.get_logger("INFO", logFile = None )

filename = "/eos/vbc/user/robert.schoefbeck/TMB/ttG_noFullyHad_fast.root"

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

# model savepath:
model_path = '.'

# for the plots
save_path = '.'


# input samples
import TMB.Samples.pp_TTGammaEFT as samples
sample = samples.ttG_noFullyHad_fast

# fix random seed for reproducibility
np.random.seed(1)

regression_target = "FI_ctZ_SM"

#mva_variables = config.mva_variables.keys()

mva_variables = [
    #'mva_photonJetdR',
    #'mva_photonLepdR',
    #'mva_photon_eta',
    'mva_photon_pt']

mva_variables.sort()
n_var_input   = len(mva_variables)

upfile = uproot.open(filename)
df     = upfile["Events"].pandas.df(branches = mva_variables+[regression_target])

#branches = upfile['Events'].arrays(namedecode='utf-8')
#basic    = (branches['mva_m3'] >= 0 ) 

df = df.dropna() # removes all Events with nan -> amounts to M3 cut

# split dataset into Input and output data
dataset = df.values
X = dataset[:,0:n_var_input]
#Y = np.log(dataset[:,n_var_input])
Y = np.log(dataset[:,n_var_input])
Y[Y < -15] = -15

# split data into train and test, test_size = 0.2 is quite standard for this
from sklearn.model_selection import train_test_split
X_train_val, X_test, Y_train_val, Y_test = train_test_split(X, Y, test_size=0.2, random_state=7, shuffle = True)


#########################################################################################
# define model (neural network)
from keras.models import Sequential, Model
from keras.optimizers import SGD
from keras.layers import Input, Activation, Dense, Convolution2D, MaxPooling2D, Dropout, Flatten
from keras.layers import BatchNormalization
from keras.utils import np_utils

model = Sequential()
model.add(BatchNormalization(input_shape=(n_var_input, )))

layers = [n_var_input+20]
for dim in layers:
    model.add(Dense(dim, activation='sigmoid'))

model.add(Dense(1))

model.compile(optimizer='adam', loss='mean_squared_error', metrics=['mean_absolute_percentage_error'])
model.summary()

# define callback for early stopping
import tensorflow as tf
callback = tf.keras.callbacks.EarlyStopping(monitor='loss', patience=10) # patience can be higher if a more accurate result is preferred
                                                                        # I would recommmend at least 3, otherwise it might cancel too early

# train the model
batch_size = 1024*4
history = model.fit(X_train_val, 
                    Y_train_val, 
                    epochs=100, 
                    batch_size=batch_size,
                    #verbose=0, # switch to 1 for more verbosity, 'silences' the output
                    callbacks=[callback],
                    #validation_split=0.1
                    validation_data=(X_test,Y_test) # use either validation_split or validation_data
                   )
print('trainig finished')

# saving
model.save(model_path + '_regression_model.h5')

#########################################################################################
# Apply the model

import array
x,y = pickle.load(file('tmp.pkl'))
c = ROOT.TCanvas()
g = ROOT.TGraph(len(Y_test), Y_test, model.predict(X_test)) 
g.Draw('ap')
c.Print("/mnt/hephy/cms/robert.schoefbeck/www/TMB/MVA/scatter.png")
