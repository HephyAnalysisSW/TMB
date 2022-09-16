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

gridpack_directory = "/eos/vbc/group/cms/robert.schoefbeck/gridpacks/ParticleNet/" 
pp_dir             = "/groups/hephy/cms/robert.schoefbeck/TMB/postprocessed/gen/v2/"

#tt1LepHad    = Sample.fromDirectory("tt1LepHad", texName = "tt1LepHad", directory = [os.path.join( pp_dir, "tt1LepHad" )])
#tt1LepHad.reweight_pkl = os.path.join(gridpack_directory, "tt01j-1l-NPtHad_HT800_slc7_amd64_gcc700_CMSSW_10_6_19_tarball.pkl")

tschRefPoint    = Sample.fromDirectory("tschRefPoint", texName = "single-t (s-ch.)", directory = [os.path.join( pp_dir, "tschRefPoint" )])
tschRefPoint.reweight_pkl = os.path.join(pp_dir, "tschRefPoint/t-sch-RefPoint_reweight_card.pkl")

tschRefPointNoWidthRW    = Sample.fromDirectory("tschRefPointNoWidthRW", texName = "single-t (s-ch.)", directory = [os.path.join( pp_dir, "tschRefPointNoWidthRW" )])
tschRefPointNoWidthRW.reweight_pkl = os.path.join(pp_dir, "tschRefPoint/t-sch-RefPoint_reweight_card.pkl")
