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

gridpack_dir = "/eos/vbc/user/robert.schoefbeck/gridpacks/v6/"
#pp_dir      = "cratch-cbe/users/robert.schoefbeck/TMB/postprocessed/gen/v6_small/ttG_noFullyHad_test/"

# no reference point samples 8/3
#WW               = Sample.fromDirectory("WW",  texName = "WW", directory = [os.path.join( pp_dir, "WW" )]) 
#WW.reweight_pkl  = os.path.join( gridpack_dir, "boson_v1_reweight_card.pkl" ) 
#WW.objects       = ['W']
#
#WZ               = Sample.fromDirectory("WZ",  texName = "WZ", directory = [os.path.join( pp_dir, "WZ" )]) 
#WZ.reweight_pkl  = os.path.join( gridpack_dir, "boson_v1_reweight_card.pkl" ) 
#WZ.objects       = ['W','Z']
#
#ZZ               = Sample.fromDirectory("ZZ",  texName = "ZZ", directory = [os.path.join( pp_dir, "ZZ" )]) 
#ZZ.reweight_pkl  = os.path.join( gridpack_dir, "boson_v1_reweight_card.pkl" ) 
#ZZ.objects       = ['Z']
#
#WWW               = Sample.fromDirectory("WWW",  texName = "WWW", directory = [os.path.join( pp_dir, "WWW" )]) 
#WWW.reweight_pkl  = os.path.join( gridpack_dir, "boson_v1_reweight_card.pkl" ) 
#WWW.objects       = ['W']
#
#WWZ               = Sample.fromDirectory("WWZ",  texName = "WWZ", directory = [os.path.join( pp_dir, "WWZ" )]) 
#WWZ.reweight_pkl  = os.path.join( gridpack_dir, "boson_v1_reweight_card.pkl" ) 
#WWZ.objects       = ['W', 'Z']
#
#WZZ               = Sample.fromDirectory("WZZ",  texName = "WZZ", directory = [os.path.join( pp_dir, "WZZ" )]) 
#WZZ.reweight_pkl  = os.path.join( gridpack_dir, "boson_v1_reweight_card.pkl" ) 
#WZZ.objects       = ['W', 'Z']
#
#ZZZ               = Sample.fromDirectory("ZZZ",  texName = "ZZZ", directory = [os.path.join( pp_dir, "ZZZ" )]) 
#ZZZ.reweight_pkl  = os.path.join( gridpack_dir, "boson_v1_reweight_card.pkl" )
#ZZZ.objects       = ['Z']
#
#ZA                = Sample.fromDirectory("ZA", texName = "ZA", directory = [os.path.join( pp_dir, "ZA" )])
#ZA.reweight_pkl   = os.path.join( gridpack_dir, "boson_v1_reweight_card.pkl" )
#ZA.objects        = ['Z', 'g']
#
#WA                = Sample.fromDirectory("WA", texName = "WA", directory = [os.path.join( pp_dir, "WA" )])
#WA.reweight_pkl   = os.path.join( gridpack_dir, "boson_v1_reweight_card.pkl" )
#WA.objects        = ['W', 'g']
#
#WA_LO             = Sample.fromDirectory("WA_LO", texName = "WA_LO", directory = [os.path.join( pp_dir, "WA_LO" )])
#WA_LO.reweight_pkl= os.path.join( gridpack_dir, "boson_v1_reweight_card.pkl" )
#WA_LO.objects     = ['W', 'g']
#
#WAjj              = Sample.fromDirectory("WAjj", texName = "WAjj", directory = [os.path.join( pp_dir, "WAjj" )])
#WAjj.reweight_pkl = os.path.join( gridpack_dir, "boson_v1_reweight_card.pkl" )
#WAjj.objects      = ['W', 'g']
#
#WWjj_OS           = Sample.fromDirectory("WWjj_OS", texName = "WWjj_OS", directory = [os.path.join( pp_dir, "WWjj_OS" )])
#WWjj_OS.reweight_pkl = os.path.join( gridpack_dir, "boson_v1_reweight_card.pkl" )
#WWjj_OS.objects   = ['W']
#
#WWjj_SS           = Sample.fromDirectory("WWjj_SS", texName = "WWjj_SS", directory = [os.path.join( pp_dir, "WWjj_SS" )])
#WWjj_SS.reweight_pkl = os.path.join( gridpack_dir, "boson_v1_reweight_card.pkl" )
#WWjj_SS.objects   = ['W']
#
#WZjj              = Sample.fromDirectory("WZjj", texName = "WZjj", directory = [os.path.join( pp_dir, "WZjj" )])
#WZjj.reweight_pkl = os.path.join( gridpack_dir, "boson_v1_reweight_card.pkl" )
#WZjj.objects      = ['W', 'Z']
#
#ZAjj              = Sample.fromDirectory("ZAjj", texName = "ZAjj", directory = [os.path.join( pp_dir, "ZAjj" )])
#ZAjj.reweight_pkl = os.path.join( gridpack_dir, "boson_v1_reweight_card.pkl" )
#ZAjj.objects      = ['Z', 'g']
#
#ZZjj              = Sample.fromDirectory("ZZjj", texName = "ZZjj", directory = [os.path.join( pp_dir, "ZZjj" )])
#ZZjj.reweight_pkl = os.path.join( gridpack_dir, "boson_v1_reweight_card.pkl" )
#ZZjj.objects      = ['Z']

pp_dir      = "/scratch-cbe/users/robert.schoefbeck/TMB/postprocessed/gen/v6/"
ttG_noFullyHad_test    = Sample.fromDirectory("ttG_noFullyHad_test", texName = "ttG_noFullyHad_test", directory = [os.path.join( pp_dir, "ttG_noFullyHad_test" )])
ttG_noFullyHad_test.reweight_pkl = os.path.join( gridpack_dir, "ttG_noFullyHad_test_ctZ_slc7_amd64_gcc700_CMSSW_10_6_0_tarball.pkl" )
ttG_noFullyHad_test.objects      = ['t', 'W', 'g']

ttg_noFullyHad_test2    = Sample.fromDirectory("ttg_noFullyHad_test2", texName = "ttg_noFullyHad_test2", directory = [os.path.join( pp_dir, "ttg_noFullyHad_test2" )])
ttg_noFullyHad_test2.reweight_pkl = os.path.join( gridpack_dir, "ttg_noFullyHad_slc7_amd64_gcc700_CMSSW_10_6_0_tarball.pkl" )
ttg_noFullyHad_test2.objects      = ['t', 'W', 'g']

#ttW01j            = Sample.fromDirectory("ttW01j", texName = "ttW01j", directory = [os.path.join( pp_dir, "ttW01j" )])
#ttW01j.reweight_pkl = os.path.join( gridpack_dir, "top_boson_v1_reweight_card.pkl" )
#ttW01j.objects    = ['t', 'W']
#
#ttZ01j            = Sample.fromDirectory("ttZ01j", texName = "ttZ01j", directory = [os.path.join( pp_dir, "ttZ01j" )])
#ttZ01j.reweight_pkl = os.path.join( gridpack_dir, "top_boson_reweight_card.pkl" )
#ttZ01j.objects    = ['t', 'W', 'Z']
