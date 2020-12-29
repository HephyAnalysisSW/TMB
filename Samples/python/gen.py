''' Benchmark samples for TTXPheno (EDM)'''

# standard imports
import os

# RootTools
from RootTools.core.standard import *

# TTXPheno
from TMB.Tools.user import cache_directory, gridpack_directory 

# sqlite3 sample cache file
dbFile = os.path.join( cache_directory, 'sample_cache', 'gen.db')
overwrite = False

# Logging
if __name__ == "__main__":
    import TTXPheno.Tools.logger as logger
    logger = logger.get_logger('DEBUG')
    import RootTools.core.logger as logger_rt
    logger_rt = logger_rt.get_logger('DEBUG')

# VV and VVV with boson reweighting
# Only leptonic decays for all V!!

test = FWLiteSample.fromFiles("test", ["/eos/vbc/experiments/cms/store/user/schoef/WW-test/WW-test/201223_221102/0000/GEN_LO_0j_93X_1.root"] )
test.reweight_pkl = os.path.join(gridpack_directory, "boson_v1_reweight_card.pkl")
test.xsec         = 1 
test.nEvents      = 100 


WW  = FWLiteSample.fromDAS("WW", "/WW-v4/schoef-WW-v4-238a1f3c8105c56c183c394f2927064f/USER", "phys03", dbFile = dbFile, overwrite=overwrite)
WW.reweight_pkl = os.path.join(gridpack_directory, "boson_v1_reweight_card.pkl")
WW.xsec         = 1.751 
WW.nEvents      = 4990816 

WZ  = FWLiteSample.fromDAS("WZ", "/WZ-v4/schoef-WZ-v4-093f53fdc754836a25887f03d0613010/USER", "phys03", dbFile = dbFile, overwrite=overwrite)
WZ.reweight_pkl = os.path.join(gridpack_directory, "boson_v1_reweight_card.pkl")
WZ.xsec         = 1.751 
WZ.nEvents      = 5209260

ZZ  = FWLiteSample.fromDAS("ZZ", "/ZZ-v4/schoef-ZZ-v4-f26891a5e446d02810e35e2cfeda4b23/USER", "phys03", dbFile = dbFile, overwrite=overwrite)
ZZ.reweight_pkl = os.path.join(gridpack_directory, "boson_v1_reweight_card.pkl")
ZZ.xsec         = 110.74415084
ZZ.nEvents      = 5144574

WWW = FWLiteSample.fromDAS("WWW", "/WWW-v4/schoef-WWW-v4-2d24f84d57d61aa0ed9149b9c9610b37/USER", "phys03", dbFile = dbFile, overwrite=overwrite)
WWW.reweight_pkl = os.path.join(gridpack_directory, "boson_v1_reweight_card.pkl")
WWW.xsec         = 0.01834 
WWW.nEvents      = 3991770 

WWZ = FWLiteSample.fromDAS("WWZ", "/WWZ-v4/schoef-WWZ-v4-d2861d31cc6c25778d9def6ca333eb5d/USER", "phys03", dbFile = dbFile, overwrite=overwrite)
WWZ.reweight_pkl = os.path.join(gridpack_directory, "boson_v1_reweight_card.pkl")
WWZ.xsec         = 0.004744 
WWZ.nEvents      = 4162978

WZZ = FWLiteSample.fromDAS("WZZ", "/WZZ-v4/schoef-WZZ-v4-808cecf77e16e2367f66c7edcf59b94d/USER", "phys03", dbFile = dbFile, overwrite=overwrite)
WZZ.reweight_pkl = os.path.join(gridpack_directory, "boson_v1_reweight_card.pkl")
WZZ.xsec         = 0.0003355 
WZZ.nEvents      = 4399466

ZZZ = FWLiteSample.fromDAS("ZZZ", "/ZZZ-v4/schoef-ZZZ-v4-fdf8e8138c23af78baac6d1fc1ca4a27/USER", "phys03", dbFile = dbFile, overwrite=overwrite)
ZZZ.reweight_pkl = os.path.join(gridpack_directory, "boson_v1_reweight_card.pkl")
ZZZ.xsec         = 0.00003884 
ZZZ.nEvents      = 4272023
