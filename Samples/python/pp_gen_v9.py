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

gridpack_directory = "/eos/vbc/user/robert.schoefbeck/gridpacks/flavor/vec/"
pp_dir       = "/scratch-cbe/users/robert.schoefbeck/TMB/postprocessed/gen/v9"

#ttZ01j    = Sample.fromDirectory("ttZ01j", texName = "ttZ01j", directory = [os.path.join( pp_dir, "ttZ01j" )])
#ttZ01j.reweight_pkl = os.path.join(gridpack_directory, "ttZ01j-vec_reweight_card.pkl")
#ttZ01j.objects      = ['t', 'W', 'Z']

WZTo3L1Nu    = Sample.fromDirectory("WZTo3L1Nu", texName = "WZTo3L1Nu", directory = [os.path.join( pp_dir, "WZTo3L1Nu" )])
WZTo3L1Nu.reweight_pkl = os.path.join(gridpack_directory, "WZTo3L1Nu-vec-cw_reweight_card.pkl")
WZTo3L1Nu.objects      = [ 'W', 'Z']

#ZZ    = Sample.fromDirectory("ZZ", texName = "ZZ", directory = [os.path.join( pp_dir, "ZZ" )])
#ZZ.reweight_pkl = os.path.join(gridpack_directory, "ZZ-vec_reweight_card.pkl")
#ZZ.objects      = ['Z']
