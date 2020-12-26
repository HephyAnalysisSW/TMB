''' Benchmark samples for TopEFT (EDM)'''

# standard imports
import os

# RootTools
from RootTools.core.standard import *

# Logging
import logging
logger = logging.getLogger(__name__)

# Logging
if __name__ == "__main__":
    import TMB.Tools.logger as logger
    logger = logger.get_logger('DEBUG')
    import RootTools.core.logger as logger_rt
    logger_rt = logger_rt.get_logger('DEBUG')

import glob

gridpack_dir = "/eos/vbc/user/robert.schoefbeck/gridpacks/"
pp_dir      = "/scratch-cbe/users/robert.schoefbeck/TMB/postprocessed/gen/v1/"

# no reference point samples 8/3
WW               = Sample.fromDirectory("WW",  texName = "WW", directory = [os.path.join( pp_dir, "WW" )]) 
WW.reweight_pkl  = os.path.join( pp_dir, "boson_v1_reweight_card.pkl" ) 

WZ               = Sample.fromDirectory("WZ",  texName = "WZ", directory = [os.path.join( pp_dir, "WZ" )]) 
WZ.reweight_pkl  = os.path.join( pp_dir, "boson_v1_reweight_card.pkl" ) 

ZZ               = Sample.fromDirectory("ZZ",  texName = "ZZ", directory = [os.path.join( pp_dir, "ZZ" )]) 
ZZ.reweight_pkl  = os.path.join( pp_dir, "boson_v1_reweight_card.pkl" ) 

WWW               = Sample.fromDirectory("WWW",  texName = "WWW", directory = [os.path.join( pp_dir, "WWW" )]) 
WWW.reweight_pkl  = os.path.join( pp_dir, "boson_v1_reweight_card.pkl" ) 

WWZ               = Sample.fromDirectory("WWZ",  texName = "WWZ", directory = [os.path.join( pp_dir, "WWZ" )]) 
WWZ.reweight_pkl  = os.path.join( pp_dir, "boson_v1_reweight_card.pkl" ) 

WZZ               = Sample.fromDirectory("WZZ",  texName = "WZZ", directory = [os.path.join( pp_dir, "WZZ" )]) 
WZZ.reweight_pkl  = os.path.join( pp_dir, "boson_v1_reweight_card.pkl" ) 

ZZZ               = Sample.fromDirectory("ZZZ",  texName = "ZZZ", directory = [os.path.join( pp_dir, "ZZZ" )]) 
ZZZ.reweight_pkl  = os.path.join( pp_dir, "boson_v1_reweight_card.pkl" ) 
