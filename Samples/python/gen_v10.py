''' GEN samples for TMB'''

# standard imports
import os

# RootTools
from RootTools.core.standard import *

# TMB
from TMB.Tools.user import cache_directory, gridpack_directory 

# sqlite3 sample cache file
dbFile = os.path.join( cache_directory, 'sample_cache', 'gen_v10.db')
overwrite = False

# Logging
if __name__ == "__main__":
    import TMB.Tools.logger as logger
    logger = logger.get_logger('DEBUG')
    import RootTools.core.logger as logger_rt
    logger_rt = logger_rt.get_logger('DEBUG')

# for flavor analysis 

gridpack_directory = "/eos/vbc/user/robert.schoefbeck/gridpacks/v7/"
WGjj = FWLiteSample.fromDAS("WGjj", "/WGjj_VBF/schoef-WGjj_VBF-ded83965d0d897daece5550e2bd76aff/USER", "phys03", dbFile = dbFile, overwrite=overwrite, skipCheck=True)
WGjj.xsec         = 5.245e-01 
WGjj.reweight_pkl = os.path.join(gridpack_directory, "WGjj_VBF_reweight_card.pkl")
WGjj.nEvents      = 29450000 

WGToLNu  = FWLiteSample.fromDAS("WGToLNu", "/WGToLNu/schoef-WGToLNu-3396e5d19c87e3124f7e9b553fcaca2c/USER", "phys03", dbFile = dbFile, overwrite=overwrite, skipCheck=True)
WGToLNu.xsec = 6.275e+01  
WGToLNu.reweight_pkl = os.path.join(gridpack_directory, "WGToLNu_reweight_card.pkl")
WGToLNu.nEvents      = 11298801

#ZZ  = FWLiteSample.fromDAS("ZZ", "/ZZ/schoef-ZZ-4d8081e3495990e85b27c44ea3ccd438/USER", "phys03", dbFile = dbFile, overwrite=overwrite, skipCheck=True)
#ZZ.xsec = 2.245e-01

ZZjj  = FWLiteSample.fromDAS("ZZjj", "/ZZjj_VBF/schoef-ZZjj_VBF-fb8d234640c8186dca4b1e7cc5fc4553/USER", "phys03", dbFile = dbFile, overwrite=overwrite, skipCheck=True)
ZZjj.xsec=1.096e-03
ZZjj.reweight_pkl = os.path.join(gridpack_directory, "ZZjj_VBF_reweight_card.pkl")
ZZjj.nEvents      = 25650000
