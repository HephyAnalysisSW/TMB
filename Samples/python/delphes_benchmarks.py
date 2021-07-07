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

gridpack_dir = "/mnt/hephy/cms/lukas.lechner/gridpacks/TTXPheno/gridpacks/"
data_dir     = "/mnt/hephy/cms/robert.schoefbeck/afs/TTXPheno/"

ttZ_ll_LO_order2_15weights_ref_ext_delphes_RunII               = Sample.fromDirectory("ttZ_ll_LO_order2_15weights_ref_ext",  texName = "ttZ",      directory = [os.path.join( data_dir, "fwlite_ttZ_ll_LO_order2_15weights_ref" )]) 
ttZ_ll_LO_order2_15weights_ref_ext_delphes_RunII.reweight_pkl  = gridpack_dir + "18052018_ref/ttZ/order2/ttZ0j_rwgt_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl"
ttZ_ll_LO_order2_15weights_ref_ext_delphes_RunII.nEvents       = 1994*5000
ttZ_ll_LO_order2_15weights_ref_ext_delphes_RunII.xsec          = 0.5205 * (0.0915 / 0.0565) #pb ttZ, Z->ll, ttZ gridpack * ttZ NLO Daniel / ttZ LO run.py UFO
ttZ_ll_LO_order2_15weights_ref_ext_delphes_RunII.xsec14        = 0.5205 * (0.0915 / 0.0565) * (0.7152 / 0.616) #pb ttZ, Z->ll, ttZ gridpack * ttZ NLO Daniel / ttZ LO run.py UFO * ttZ jets 14 TeV / ttZ jets 13 TeV
