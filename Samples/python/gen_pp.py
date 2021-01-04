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
pp_dir      = "/scratch-cbe/users/robert.schoefbeck/TMB/postprocessed/gen/v3/"

# no reference point samples 8/3
WW               = Sample.fromDirectory("WW",  texName = "WW", directory = [os.path.join( pp_dir, "WW" )]) 
WW.reweight_pkl  = os.path.join( gridpack_dir, "boson_v1_reweight_card.pkl" ) 

WZ               = Sample.fromDirectory("WZ",  texName = "WZ", directory = [os.path.join( pp_dir, "WZ" )]) 
WZ.reweight_pkl  = os.path.join( gridpack_dir, "boson_v1_reweight_card.pkl" ) 

ZZ               = Sample.fromDirectory("ZZ",  texName = "ZZ", directory = [os.path.join( pp_dir, "ZZ" )]) 
ZZ.reweight_pkl  = os.path.join( gridpack_dir, "boson_v1_reweight_card.pkl" ) 

WWW               = Sample.fromDirectory("WWW",  texName = "WWW", directory = [os.path.join( pp_dir, "WWW" )]) 
WWW.reweight_pkl  = os.path.join( gridpack_dir, "boson_v1_reweight_card.pkl" ) 

WWZ               = Sample.fromDirectory("WWZ",  texName = "WWZ", directory = [os.path.join( pp_dir, "WWZ" )]) 
WWZ.reweight_pkl  = os.path.join( gridpack_dir, "boson_v1_reweight_card.pkl" ) 

WZZ               = Sample.fromDirectory("WZZ",  texName = "WZZ", directory = [os.path.join( pp_dir, "WZZ" )]) 
WZZ.reweight_pkl  = os.path.join( gridpack_dir, "boson_v1_reweight_card.pkl" ) 

ZZZ               = Sample.fromDirectory("ZZZ",  texName = "ZZZ", directory = [os.path.join( pp_dir, "ZZZ" )]) 
ZZZ.reweight_pkl  = os.path.join( gridpack_dir, "boson_v1_reweight_card.pkl" )

WA                = Sample.fromDirectory("WA", texName = "WA", directory = [os.path.join( pp_dir, "WA" )])
WA.reweight_pkl  = os.path.join( gridpack_dir, "boson_v1_reweight_card.pkl" )

WAjj              = Sample.fromDirectory("WAjj", texName = "WAjj", directory = [os.path.join( pp_dir, "WAjj" )])
WAjj.reweight_pkl  = os.path.join( gridpack_dir, "boson_v1_reweight_card.pkl" )

WWjj_OS           = Sample.fromDirectory("WWjj_OS", texName = "WWjj_OS", directory = [os.path.join( pp_dir, "WWjj_OS" )])
WWjj_OS.reweight_pkl  = os.path.join( gridpack_dir, "boson_v1_reweight_card.pkl" )

WWjj_SS           = Sample.fromDirectory("WWjj_SS", texName = "WWjj_SS", directory = [os.path.join( pp_dir, "WWjj_SS" )])
WWjj_SS.reweight_pkl  = os.path.join( gridpack_dir, "boson_v1_reweight_card.pkl" )

WZjj              = Sample.fromDirectory("WZjj", texName = "WZjj", directory = [os.path.join( pp_dir, "WZjj" )])
WZjj.reweight_pkl  = os.path.join( gridpack_dir, "boson_v1_reweight_card.pkl" )

ZAjj              = Sample.fromDirectory("ZAjj", texName = "ZAjj", directory = [os.path.join( pp_dir, "ZAjj" )])
ZAjj.reweight_pkl  = os.path.join( gridpack_dir, "boson_v1_reweight_card.pkl" )

ZZjj              = Sample.fromDirectory("ZZjj", texName = "ZZjj", directory = [os.path.join( pp_dir, "ZZjj" )])
ZZjj.reweight_pkl  = os.path.join( gridpack_dir, "boson_v1_reweight_card.pkl" )

ttg_noFullyHad    = Sample.fromDirectory("ttg_noFullyHad", texName = "ttg_noFullyHad", directory = [os.path.join( pp_dir, "ttg_noFullyHad" )])
ttg_noFullyHad.reweight_pkl  = os.path.join( gridpack_dir, "top_boson_v1_reweight_card.pkl" )

ttW01j            = Sample.fromDirectory("ttW01j", texName = "ttW01j", directory = [os.path.join( pp_dir, "ttW01j" )])
ttW01j.reweight_pkl  = os.path.join( gridpack_dir, "top_boson_v1_reweight_card.pkl" )

ttZ01j            = Sample.fromDirectory("ttZ01j", texName = "ttZ01j", directory = [os.path.join( pp_dir, "ttZ01j" )])
ttZ01j.reweight_pkl  = os.path.join( gridpack_dir, "top_boson_v1_reweight_card.pkl" )
