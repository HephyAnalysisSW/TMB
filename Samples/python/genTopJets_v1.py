''' GEN samples for TMB'''

# standard imports
import os

# RootTools
from RootTools.core.standard import *

# Logging
if __name__ == "__main__":
    import TMB.Tools.logger as logger
    logger = logger.get_logger('DEBUG')
    import RootTools.core.logger as logger_rt
    logger_rt = logger_rt.get_logger('DEBUG')

# for flavor analysis 

gridpack_directory = "/eos/vbc/group/cms/robert.schoefbeck/gridpacks/ParticleNet/"

#tt1LepHad = FWLiteSample.fromFiles("tt1LepHad", ["/eos/vbc/group/cms/robert.schoefbeck/ParticleNet/test/GEN_LO_0j_102X.root"])
#tt1LepHad.reweight_pkl = os.path.join(gridpack_directory, "tt01j-1l-NPtHad_HT800_slc7_amd64_gcc700_CMSSW_10_6_19_tarball.pkl")
tmp = FWLiteSample.fromFiles("GEN", ["/users/robert.schoefbeck/CMS/CMSSW_10_6_27/src/Samples/cfg/GEN_LO_0j_102X.root"]) 
tmp.reweight_pkl = "/users/robert.schoefbeck/tt01j-1l-NPtHad_reweight_card.pkl" 
