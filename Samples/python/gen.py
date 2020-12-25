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

# WZ with boson reweighting
WZ              = FWLiteSample.fromDAS("WZ", "/WZ-v2/schoef-WZ-v2-093f53fdc754836a25887f03d0613010/USER", "phys03", dbFile = dbFile, overwrite=overwrite, )
WZ.reweight_pkl = os.path.join(gridpack_directory, "boson_v1_reweight_card.pkl")
WZ.xsec         = 110.74415084
WZ.nEvents      = 1000 
# ZZ with boson reweighting
ZZ              = FWLiteSample.fromDAS("ZZ", "/ZZ-v2/schoef-ZZ-v2-f26891a5e446d02810e35e2cfeda4b23/USER", "phys03", dbFile = dbFile, overwrite=overwrite )
ZZ.reweight_pkl = os.path.join(gridpack_directory, "boson_v1_reweight_card.pkl")
ZZ.xsec         = 33.59376163
ZZ.nEvents      = 1000 

