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

gridpack_directory = "/eos/vbc/user/robert.schoefbeck/gridpacks/flavor/vec/"
pp_dir       = "/scratch-cbe/users/robert.schoefbeck/TMB/postprocessed/gen/v10"


gridpack_directory = "/eos/vbc/user/robert.schoefbeck/gridpacks/v7/"

WGjj = Sample.fromDirectory("WGjj", texName = "WGjj", directory = [os.path.join( pp_dir, "WGjj" )]) 
WGjj.reweight_pkl = os.path.join(gridpack_directory, "WGjj_VBF_reweight_card.pkl")
WGjj.objects      = ['g', 'W']

WGToLNu = Sample.fromDirectory("WGToLNu", texName = "WGToLNu", directory = [os.path.join( pp_dir, "WGToLNu" )]) 
WGToLNu.reweight_pkl = os.path.join(gridpack_directory, "WGToLNu_reweight_card.pkl")
WGToLNu.objects      = ['g', 'W']

WG20To130ToLNu = Sample.fromDirectory("WG20To130ToLNu", texName = "WG20To130ToLNu", directory = [os.path.join( pp_dir, "WG20To130ToLNu" )]) 
WG20To130ToLNu.reweight_pkl = os.path.join(gridpack_directory, "WG20To130ToLNu_reweight_card.pkl")
WG20To130ToLNu.objects      = ['g', 'W']

WG130To300ToLNu = Sample.fromDirectory("WG130To300ToLNu", texName = "WG130To300ToLNu", directory = [os.path.join( pp_dir, "WG130To300ToLNu" )]) 
WG130To300ToLNu.reweight_pkl = os.path.join(gridpack_directory, "WG130To300ToLNu_reweight_card.pkl")
WG130To300ToLNu.objects      = ['g', 'W']

WG300To500ToLNu = Sample.fromDirectory("WG300To500ToLNu", texName = "WG300To500ToLNu", directory = [os.path.join( pp_dir, "WG300To500ToLNu" )]) 
WG300To500ToLNu.reweight_pkl = os.path.join(gridpack_directory, "WG300To500ToLNu_reweight_card.pkl")
WG300To500ToLNu.objects      = ['g', 'W']

WG500ToLNu = Sample.fromDirectory("WG500ToLNu", texName = "WG500ToLNu", directory = [os.path.join( pp_dir, "WG500ToLNu" )]) 
WG500ToLNu.reweight_pkl = os.path.join(gridpack_directory, "WG500ToLNu_reweight_card.pkl")
WG500ToLNu.objects      = ['g', 'W']

ZZ = Sample.fromDirectory("ZZ", texName = "ZZ", directory = [os.path.join( pp_dir, "ZZ" )]) 
ZZ.reweight_pkl = os.path.join(gridpack_directory, "ZZ_reweight_card.pkl")
ZZ.objects      = ['Z']

ZZjj = Sample.fromDirectory("ZZjj", texName = "ZZjj", directory = [os.path.join( pp_dir, "ZZjj" )]) 
ZZjj.reweight_pkl = os.path.join(gridpack_directory, "ZZjj_VBF_reweight_card.pkl")
ZZjj.objects      = ['Z']
