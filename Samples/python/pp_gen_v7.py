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

gridpack_dir = "/eos/vbc/user/robert.schoefbeck/gridpacks/flavor/order_2"
pp_dir       = "/scratch-cbe/users/robert.schoefbeck/TMB/postprocessed/gen/v7"

ttG_noFullyHad    = Sample.fromDirectory("ttG_noFullyHad", texName = "ttG_noFullyHad", directory = [os.path.join( pp_dir, "ttG_noFullyHad" )])
ttG_noFullyHad.reweight_pkl = os.path.join( gridpack_dir, ttG_noFullyHad.name+'_reweight_card.pkl' )
ttG_noFullyHad.objects      = ['t', 'W', 'g']

ttW01j    = Sample.fromDirectory("ttW01j", texName = "ttW01j", directory = [os.path.join( pp_dir, "ttW01j" )])
ttW01j.reweight_pkl = os.path.join( gridpack_dir, ttW01j.name+'_reweight_card.pkl' )
ttW01j.objects      = ['t', 'W']

ttZ01j    = Sample.fromDirectory("ttZ01j", texName = "ttZ01j", directory = [os.path.join( pp_dir, "ttZ01j" )])
ttZ01j.reweight_pkl = os.path.join( gridpack_dir, ttZ01j.name+'_reweight_card.pkl' )
ttZ01j.objects      = ['t', 'W', 'Z']

WGToLNu    = Sample.fromDirectory("WGToLNu", texName = "WGToLNu", directory = [os.path.join( pp_dir, "WGToLNu" )])
WGToLNu.reweight_pkl = os.path.join( gridpack_dir, WGToLNu.name+'_reweight_card.pkl' )
WGToLNu.objects      = ['W', 'g']

WW    = Sample.fromDirectory("WW", texName = "WW", directory = [os.path.join( pp_dir, "WW" )])
WW.reweight_pkl = os.path.join( gridpack_dir, WW.name+'_reweight_card.pkl' )
WW.objects      = ['W']

WZTo3L1Nu    = Sample.fromDirectory("WZTo3L1Nu", texName = "WZTo3L1Nu", directory = [os.path.join( pp_dir, "WZTo3L1Nu" )])
WZTo3L1Nu.reweight_pkl = os.path.join( gridpack_dir, WZTo3L1Nu.name+'_reweight_card.pkl' )
WZTo3L1Nu.objects      = [ 'W', 'Z']

ZGTo2L    = Sample.fromDirectory("ZGTo2L", texName = "ZGTo2L", directory = [os.path.join( pp_dir, "ZGTo2L" )])
ZGTo2L.reweight_pkl = os.path.join( gridpack_dir, ZGTo2L.name+'_reweight_card.pkl' )
ZGTo2L.objects      = ['Z', 'g']

ZZ    = Sample.fromDirectory("ZZ", texName = "ZZ", directory = [os.path.join( pp_dir, "ZZ" )])
ZZ.reweight_pkl = os.path.join( gridpack_dir, ZZ.name+'_reweight_card.pkl' )
ZZ.objects      = ['Z']
