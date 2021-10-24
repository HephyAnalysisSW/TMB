''' GEN samples for TMB'''

# standard imports
import os

# RootTools
from RootTools.core.standard import *

# TMB
from TMB.Tools.user import cache_directory, gridpack_directory 

# sqlite3 sample cache file
dbFile = os.path.join( cache_directory, 'sample_cache', 'gen_v8.db')
overwrite = False

# Logging
if __name__ == "__main__":
    import TMB.Tools.logger as logger
    logger = logger.get_logger('DEBUG')
    import RootTools.core.logger as logger_rt
    logger_rt = logger_rt.get_logger('DEBUG')

# for flavor analysis 

gridpack_directory = "/eos/vbc/user/robert.schoefbeck/gridpacks/flavor/vec/"

ttZ01j = FWLiteSample.fromDAS("ttZ01j", "/flavor_vec_gen_ttZ01j/schoef-flavor_vec_gen_ttZ01j-72dd436ff602e070e806a1a87ef88c24/USER", "phys03", dbFile = dbFile, overwrite=overwrite, skipCheck=True)
ttZ01j.reweight_pkl = os.path.join(gridpack_directory, "ttZ01j-vec_reweight_card.pkl")
ttZ01j.xsec         = 5.852e-02 
ttZ01j.nEvents      = 8717754

WZTo3L1Nu = FWLiteSample.fromDAS("WZTo3L1Nu", "/flavor_vec_gen_WZTo3L1Nu/schoef-flavor_vec_gen_WZTo3L1Nu-63a6e745b000ec3e2146517f610f18f1/USER", "phys03", dbFile = dbFile, overwrite=overwrite, skipCheck=True)
WZTo3L1Nu.reweight_pkl = os.path.join(gridpack_directory, "WZTo3L1Nu-vec_reweight_card.pkl")
WZTo3L1Nu.xsec    = 1.146e+00 
WZTo3L1Nu.nEvents = 15646293

ZZ = FWLiteSample.fromDAS("ZZ", "/flavor_vec_gen_ZZ/schoef-flavor_vec_gen_ZZ-351f9bda1d3fa7d3617998f8fee76c02/USER", "phys03", dbFile = dbFile, overwrite=overwrite, skipCheck=True)
ZZ.reweight_pkl = os.path.join(gridpack_directory, "ZZ-vec_reweight_card.pkl")
ZZ.xsec    = 1.146e-01
ZZ.nEvents = 15484325
