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

# TMB
from TMB.Tools.user import cache_directory, gridpack_directory 

# sqlite3 sample cache file
dbFile = os.path.join( cache_directory, 'sample_cache', 'genTopJets_v1.db')
overwrite = False

# for flavor analysis 

gridpack_directory = "/eos/vbc/group/cms/robert.schoefbeck/gridpacks/ParticleNet/"

tt1LepHad_test = FWLiteSample.fromFiles("tt1LepHad_test", ["/eos/vbc/group/cms/robert.schoefbeck/ParticleNet/test/GEN_LO_0j_102X.root"])
tt1LepHad_test.reweight_pkl = os.path.join(gridpack_directory, "tt01j-1l-NPtHad_HT800_slc7_amd64_gcc700_CMSSW_10_6_19_tarball.pkl")
#tmp = FWLiteSample.fromFiles("GEN", ["/users/robert.schoefbeck/CMS/CMSSW_10_6_27/src/Samples/cfg/GEN_LO_0j_102X.root"]) 
#tmp.reweight_pkl = "/users/robert.schoefbeck/tt01j-1l-NPtHad_reweight_card.pkl" 

tt1LepHad = FWLiteSample.fromDAS("tt1LepHad", "/PNet/schoef-PNet-54f3d650012e99eaac4060c8f33d741f/USER", "phys03", dbFile = dbFile, overwrite=overwrite, skipCheck=True)
tt1LepHad.reweight_pkl = os.path.join(gridpack_directory, "tt01j-1l-NPtHad_HT800_slc7_amd64_gcc700_CMSSW_10_6_19_tarball.pkl") 
