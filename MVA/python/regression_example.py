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

n_var_input = len(variables)

# input samples
import TMB.Samples.pp_TTGammaEFT as samples
sample = samples.ttG_noFullyHad_fast

# fix random seed for reproducibility
seed = 7
np.random.seed(seed)


    #branches = uproot['Events'].arrays(namedecode='utf-8')
    #basic    = (branches['PhotonGood0_pt'] >= 0 ) & (branches['LeptonTight0_pt'] >= 0)

regression_target = "FI_ctZ_SM"

mva_variables = config.mva_variables.keys()
mva_variables.sort()

upfile = uproot.open(filename)
df     = upfile[key]["Events"].pandas.df(branches = mva_variables+[regression_target])

## check for NaN in the dataframe, .root file might be slightly broken
#for key in key_list:
#    for var in config.mva_variables:
#        if df[key][var].isnull().values.any() == True:
#            print(key,':', var, ' has some entries as nan:')
#            print(df[key][var].isnull().sum(), ' are nan')
#            print('nan Events will be removed')

df_all = df_all.dropna() # removes all Events with nan

# split dataset into Input and output data
dataset = df_all.values
X = dataset[:,0:len(mva_variables)]
Y = dataset[:,len(mva_variables)]

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

layers = [n_var_input*5, n_var_input*5, n_var_input*5]
for dim in layers:
    model.add(Dense(dim, activation='sigmoid'))

model.add(Dense(1))

model.compile(optimizer='adam', loss='mean_squared_error', metrics=['mean_absolute_percentage_error'])
model.summary()

# define callback for early stopping
import tensorflow as tf
callback = tf.keras.callbacks.EarlyStopping(monitor='loss', patience=1) # patience can be higher if a more accurate result is preferred
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
print(model.predict(X_test[9:10]), Y_test[9:10])
