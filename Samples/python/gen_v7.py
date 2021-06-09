''' GEN samples for TMB'''

# standard imports
import os

# RootTools
from RootTools.core.standard import *

# TMB
from TMB.Tools.user import cache_directory, gridpack_directory 

# sqlite3 sample cache file
dbFile = os.path.join( cache_directory, 'sample_cache', 'gen_v7.db')
overwrite = False

# Logging
if __name__ == "__main__":
    import TMB.Tools.logger as logger
    logger = logger.get_logger('DEBUG')
    import RootTools.core.logger as logger_rt
    logger_rt = logger_rt.get_logger('DEBUG')

# for flavor analysis 
# https://cmsweb.cern.ch/das/request?view=list&limit=50&instance=prod%2Fphys03&input=%2F*flavor_o2*%2F*%2F*

WGToLNu = FWLiteSample.fromDirectory("WGToLNu", "/eos/vbc/experiments/cms/store/user/schoef/flavor_o2_WGToLNu/flavor_o2_WGToLNu/210608_200214/0000")
WGToLNu.reweight_pkl = os.path.join(gridpack_directory, "flavor/order_2", "WGToLNu_reweight_card.pkl")
WGToLNu.xsec         = 6.235e+01 
WGToLNu.nEvents      = 378275

ZGTo2L = FWLiteSample.fromDirectory("ZGTo2L", "/eos/vbc/experiments/cms/store/user/schoef/flavor_o2_ZGTo2L/flavor_o2_ZGTo2L/210608_200746/0000")
ZGTo2L.reweight_pkl = os.path.join(gridpack_directory, "flavor/order_2", "ZGTo2L_reweight_card.pkl")
ZGTo2L.xsec    = 1.027e+01
ZGTo2L.nEvents = 426937

WW = FWLiteSample.fromDirectory("WW", "/eos/vbc/experiments/cms/store/user/schoef/flavor_o2_WW/flavor_o2_WW/210608_200334/0000")
WW.reweight_pkl = os.path.join(gridpack_directory, "flavor/order_2", "WW_reweight_card.pkl")
WW.xsec    = 8.648e+00
WW.nEvents = 533994

WZTo3L1Nu = FWLiteSample.fromDirectory("WZTo3L1Nu", "/eos/vbc/experiments/cms/store/user/schoef/flavor_o2_WZTo3L1Nu/flavor_o2_WZTo3L1Nu/210608_200615/0000/")
WZTo3L1Nu.reweight_pkl = os.path.join(gridpack_directory, "flavor/order_2", "WZTo3L1Nu_reweight_card.pkl")
WZTo3L1Nu.xsec    = 1.086e+00
WZTo3L1Nu.nEvents = 537247 

ZZ = FWLiteSample.fromDirectory("ZZ", "/eos/vbc/experiments/cms/store/user/schoef/flavor_o2_ZZ/flavor_o2_ZZ/210608_200906/0000/")
ZZ.reweight_pkl = os.path.join(gridpack_directory, "flavor/order_2", "ZZ_reweight_card.pkl")
ZZ.xsec    = 1.137e-01
ZZ.nEvents = 526969 
