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

WG20To130ToLNu  = FWLiteSample.fromDAS("WG20To130ToLNu", "/WG20To130ToLNu/schoef-WG20To130ToLNu-383dc2402bb0fa0c195be4a5319db7c8/USER", "phys03", dbFile = dbFile, overwrite=overwrite, skipCheck=True)
WG20To130ToLNu.xsec = 2.548e+01 
WG20To130ToLNu.reweight_pkl = os.path.join(gridpack_directory, "WG20To130ToLNu_reweight_card.pkl")
WG20To130ToLNu.nEvents      = 10485159

WG130To300ToLNu  = FWLiteSample.fromDAS("WG130To300ToLNu", "/WG130To300ToLNu/schoef-WG130To300ToLNu-a18c6b1e17bf5a0e6bed9fd71b7ca83c/USER", "phys03", dbFile = dbFile, overwrite=overwrite, skipCheck=True)
WG130To300ToLNu.xsec = 6.481e-01  
WG130To300ToLNu.reweight_pkl = os.path.join(gridpack_directory, "WG130To300ToLNu_reweight_card.pkl")
WG130To300ToLNu.nEvents      = 9701591

WG300To500ToLNu  = FWLiteSample.fromDAS("WG300To500ToLNu", "/WG300To500ToLNu/schoef-WG300To500ToLNu-ed17bc4e8efeba8bccab6b5d5fc7f5dd/USER", "phys03", dbFile = dbFile, overwrite=overwrite, skipCheck=True)
WG300To500ToLNu.xsec = 4.782e-02  
WG300To500ToLNu.reweight_pkl = os.path.join(gridpack_directory, "WG300To500ToLNu_reweight_card.pkl")
WG300To500ToLNu.nEvents      = 9098674

WG500ToLNu  = FWLiteSample.fromDAS("WG500ToLNu", "/WG500ToLNu/schoef-WG500ToLNu-809d9b04747d0bfff5a6caabc1a59bf0/USER", "phys03", dbFile = dbFile, overwrite=overwrite, skipCheck=True)
WG500ToLNu.xsec = 7.962e-03
WG500ToLNu.reweight_pkl = os.path.join(gridpack_directory, "WG500ToLNu_reweight_card.pkl")
WG500ToLNu.nEvents      = 8810781

ZZ  = FWLiteSample.fromDAS("ZZ", "/ZZ/schoef-ZZ-3d151550d5cebaa869431ac7c227b1f2/USER", "phys03", dbFile = dbFile, overwrite=overwrite, skipCheck=True)
ZZ.xsec = 1.187e-01
ZZ.reweight_pkl = os.path.join(gridpack_directory, "ZZ_reweight_card.pkl")
ZZ.nEvents      = 15139425

ZZjj  = FWLiteSample.fromDAS("ZZjj", "/ZZjj_VBF/schoef-ZZjj_VBF-fb8d234640c8186dca4b1e7cc5fc4553/USER", "phys03", dbFile = dbFile, overwrite=overwrite, skipCheck=True)
ZZjj.xsec=1.096e-03
ZZjj.reweight_pkl = os.path.join(gridpack_directory, "ZZjj_VBF_reweight_card.pkl")
ZZjj.nEvents      = 25650000

DYJets = FWLiteSample.fromDAS("DYJets", "/DYJets/schoef-DYJets-ceb522b2459205cbd7e9c77e14159897/USER", "phys03", dbFile = dbFile, overwrite=overwrite, skipCheck=True)
DYJets.xsec    = 5.475e+03
DYJets.nEvents = 13472406

WJetsToLNu = FWLiteSample.fromDAS("WJetsToLNu", "/WJetsToLNu/schoef-WJetsToLNu-14ef8bf820f268bd24130b312e7bbd9b/USER", "phys03", dbFile = dbFile, overwrite=overwrite, skipCheck=True)
WJetsToLNu.xsec     = 5.408e+04
WJetsToLNu.nEvents  = 14287034

TTJets = FWLiteSample.fromDAS("TTJets", "/TTJets_LO/schoef-TTJets_LO-2354f74db9d92f5749622da9bb4dfc72/USER", "phys03", dbFile = dbFile, overwrite=overwrite, skipCheck=True)
TTJets.xsec    = 4.833e+02
TTJets.nEvents = 2324192

WH = FWLiteSample.fromDAS("WH", "/WH/schoef-WH-bf335b0285fb1dbca5c66b3278d630c1/USER", "phys03", dbFile = dbFile, overwrite=overwrite, skipCheck=True)
WH.xsec     = 3.779e-01
WH.nEvents  = 15066155
WH.reweight_pkl = "/eos/vbc/user/robert.schoefbeck/gridpacks/VH/SMEFTsim_VH_reweight_card.pkl"

ZH = FWLiteSample.fromDAS("ZH", "/ZH/schoef-ZH-fa0e05e99c7eb8089d963cc72f0f9faf/USER", "phys03", dbFile = dbFile, overwrite=overwrite, skipCheck=True)
ZH.xsec    = 5.846e-02
ZH.nEvents = 14769665
ZH.reweight_pkl = "/eos/vbc/user/robert.schoefbeck/gridpacks/VH/SMEFTsim_VH_reweight_card.pkl"
