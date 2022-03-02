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

gridpack_directory = "/eos/vbc/group/cms/robert.schoefbeck/gridpacks/v7/"
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

DYJets_HT_70to100_AODSIM = FWLiteSample.fromDAS("DYJets_HT_70to100_AODSIM", "/DYJetsToLL_M-50_HT-70to100_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIIAutumn18DRPremix-102X_upgrade2018_realistic_v15-v1/AODSIM", "phys03", dbFile = dbFile, overwrite=overwrite, skipCheck=True)
DYJets_HT_70to100_AODSIM.nEvents = 10077604
DYJets_HT_70to100_AODSIM.xsec    = 1.235e+02

#DYJets_HT_70to100_MINIAODSIM = FWLiteSample.fromDAS("DYJets_HT_70to100_MINIAODSIM", "/DYJetsToLL_M-50_HT-70to100_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM", "phys03", dbFile = dbFile, overwrite=overwrite, skipCheck=True)
#DYJets_HT_70to100_MINIAODSIM.nEvents = 10019684
#DYJets_HT_70to100_MINIAODSIM.xsec    = 1.235e+02

#498
DYJets_HT_70to100 = FWLiteSample.fromDAS("DYJets_HT_70to100", "/DYJets_HT-70to100/schoef-DYJets_HT-70to100-6b303b663d22933412cd7b1ff512a296/USER", "phys03", dbFile = dbFile, overwrite=overwrite, skipCheck=True)
DYJets_HT_70to100.nEvents = 735956
DYJets_HT_70to100.xsec    = 1.235e+02

#500
DYJets_HT_100to200 = FWLiteSample.fromDAS("DYJets_HT_100to200", "/DYJets_HT-100to200/schoef-DYJets_HT-100to200-4f8c750d548dd65ec9aa8edea87a4a02/USER", "phys03", dbFile = dbFile, overwrite=overwrite, skipCheck=True)
DYJets_HT_100to200.nEvents = 620456
DYJets_HT_100to200.xsec    = 1.353e+02

#499
DYJets_HT_200to400 = FWLiteSample.fromDAS("DYJets_HT_200to400", "/DYJets_HT-200to400/schoef-DYJets_HT-200to400-4bf38cc0584e5abf92fd7a9baaf94731/USER", "phys03", dbFile = dbFile, overwrite=overwrite, skipCheck=True)
DYJets_HT_200to400.nEvents = 550996
DYJets_HT_200to400.xsec    = 4.399e+01

#200
DYJets_HT_400to600 = FWLiteSample.fromDAS("DYJets_HT_400to600", "/DYJets_HT-400to600/schoef-DYJets_HT-400to600-9971d8727dcb8732f6d2b6aa95c4cdce/USER", "phys03", dbFile = dbFile, overwrite=overwrite, skipCheck=True)
DYJets_HT_400to600.nEvents = 103084
DYJets_HT_400to600.xsec    = 6.040e+00

#195
DYJets_HT_600to800 = FWLiteSample.fromDAS("DYJets_HT_600to800", "/DYJets_HT-600to800/schoef-DYJets_HT-600to800-62c1f7c295f35b0b5b9baa81c061382d/USER", "phys03", dbFile = dbFile, overwrite=overwrite, skipCheck=True)
DYJets_HT_600to800.nEvents = 97558
DYJets_HT_600to800.xsec    = 1.524e+00

#200
DYJets_HT_800to1200 = FWLiteSample.fromDAS("DYJets_HT_800to1200", "/DYJets_HT-800to1200/schoef-DYJets_HT-800to1200-be0e2b9e34c3092d6c43d77c4b7cb0e0/USER", "phys03", dbFile = dbFile, overwrite=overwrite, skipCheck=True)
DYJets_HT_800to1200.nEvents = 96916
DYJets_HT_800to1200.xsec    = 7.665e-01

#1962
DYJets_HT_1200to2500 = FWLiteSample.fromDAS("DYJets_HT_1200to2500", "/DYJets_HT-1200to2500/schoef-DYJets_HT-1200to2500-12e6a3c76a833f7cd33046fc1529d076/USER", "phys03", dbFile = dbFile, overwrite=overwrite, skipCheck=True)
DYJets_HT_1200to2500.nEvents = 929027
DYJets_HT_1200to2500.xsec    = 1.872e-01

#128
DYJets_HT_2500toInf = FWLiteSample.fromDAS("DYJets_HT_2500toInf", "/DYJets_HT-2500toInf/schoef-DYJets_HT-2500toInf-d0684908fdf6bb6d1c1498a95ca6871a/USER", "phys03", dbFile = dbFile, overwrite=overwrite, skipCheck=True)
DYJets_HT_2500toInf.nEvents = 63444
DYJets_HT_2500toInf.xsec    = 3.431e-03

WJetsToLNu_HT_70to100 = FWLiteSample.fromDAS("WJetsToLNu_HT_70to100", "/WJetsToLNu_HT-70to100/schoef-WJetsToLNu_HT-70to100-bd4a8d8de733317826409edb551d8252/USER", "phys03", dbFile = dbFile, overwrite=overwrite, skipCheck=True)
WJetsToLNu_HT_70to100.nEvents = 659918
WJetsToLNu_HT_70to100.xsec    = 1.086e+03
#498

WJetsToLNu_HT_100to200 = FWLiteSample.fromDAS("WJetsToLNu_HT_100to200", "/WJetsToLNu_HT-100to200/schoef-WJetsToLNu_HT-100to200-dec4781bb005a405a35a19ef42d5f1dd/USER", "phys03", dbFile = dbFile, overwrite=overwrite, skipCheck=True)
WJetsToLNu_HT_100to200.xsec     = 1.173e+03
WJetsToLNu_HT_100to200.nEvents  = 575604
#499

WJetsToLNu_HT_200to400 = FWLiteSample.fromDAS("WJetsToLNu_HT_200to400", "/WJetsToLNu_HT-200to400/schoef-WJetsToLNu_HT-200to400-25f879b5f8e82f7f10cc9e604590572f/USER", "phys03", dbFile = dbFile, overwrite=overwrite, skipCheck=True)
WJetsToLNu_HT_200to400.xsec     = 3.680e+02
WJetsToLNu_HT_200to400.nEvents  = 532790
#500

WJetsToLNu_HT_400to600 = FWLiteSample.fromDAS("WJetsToLNu_HT_400to600", "/WJetsToLNu_HT-400to600/schoef-WJetsToLNu_HT-400to600-502fbcea228579016f4127fbe81015a3/USER", "phys03", dbFile = dbFile, overwrite=overwrite, skipCheck=True)
WJetsToLNu_HT_400to600.xsec     = 5.600e+01
WJetsToLNu_HT_400to600.nEvents  = 101111
#200

WJetsToLNu_HT_600to800 = FWLiteSample.fromDAS("WJetsToLNu_HT_600to800", "/WJetsToLNu_HT-600to800/schoef-WJetsToLNu_HT-600to800-aa2c3a2e3e720b3232f3b592514c0fba/USER", "phys03", dbFile = dbFile, overwrite=overwrite, skipCheck=True)
WJetsToLNu_HT_600to800.xsec    = 1.471e+01
WJetsToLNu_HT_600to800.nEvents = 96116
#198

WJetsToLNu_HT_800to1200 = FWLiteSample.fromDAS("WJetsToLNu_HT_800to1200", "/WJetsToLNu_HT-800to1200/schoef-WJetsToLNu_HT-800to1200-e6efc9c649846793d727b0698f13dc5a/USER", "phys03", dbFile = dbFile, overwrite=overwrite, skipCheck=True)
WJetsToLNu_HT_800to1200.xsec    = 6.807e+00
WJetsToLNu_HT_800to1200.nEvents = 92675
#196

WJetsToLNu_HT_1200to2500 = FWLiteSample.fromDAS("WJetsToLNu_HT_1200to2500", "/WJetsToLNu_HT-1200to2500/schoef-WJetsToLNu_HT-1200to2500-c898df7427c9ecc9d386d96f8da1d67b/USER", "phys03", dbFile = dbFile, overwrite=overwrite, skipCheck=True)
WJetsToLNu_HT_1200to2500.xsec   = 1.627e+00
WJetsToLNu_HT_1200to2500.nEvents= 91634
#200

WJetsToLNu_HT_2500toInf = FWLiteSample.fromDAS("WJetsToLNu_HT_2500toInf", "/WJetsToLNu_HT-2500toInf/schoef-WJetsToLNu_HT-2500toInf-adf0c94d1fa7de330f987ac552330ea0/USER", "phys03", dbFile = dbFile, overwrite=overwrite, skipCheck=True)
WJetsToLNu_HT_2500toInf.xsec    = 3.844e-02
WJetsToLNu_HT_2500toInf.nEvents = 87154
#196


WJetsToLNu = FWLiteSample.fromDAS("WJetsToLNu", "/WJetsToLNu/schoef-WJetsToLNu-14ef8bf820f268bd24130b312e7bbd9b/USER", "phys03", dbFile = dbFile, overwrite=overwrite, skipCheck=True)
WJetsToLNu.xsec     = 5.408e+04
WJetsToLNu.nEvents  = 14287034

TTJets = FWLiteSample.fromDAS("TTJets", "/TTJets_LO/schoef-TTJets_LO-2354f74db9d92f5749622da9bb4dfc72/USER", "phys03", dbFile = dbFile, overwrite=overwrite, skipCheck=True)
TTJets.xsec    = 4.833e+02
TTJets.nEvents = 2324192

WH = FWLiteSample.fromDAS("WH", "/WH/schoef-WH-bf335b0285fb1dbca5c66b3278d630c1/USER", "phys03", dbFile = dbFile, overwrite=overwrite, skipCheck=True)
WH.xsec     = 3.779e-01
WH.nEvents  = 15066155
WH.reweight_pkl = "/eos/vbc/group/cms/robert.schoefbeck/gridpacks/VH/SMEFTsim_VH_reweight_card.pkl"

ZH = FWLiteSample.fromDAS("ZH", "/ZH/schoef-ZH-fa0e05e99c7eb8089d963cc72f0f9faf/USER", "phys03", dbFile = dbFile, overwrite=overwrite, skipCheck=True)
ZH.xsec    = 5.846e-02
ZH.nEvents = 14769665
ZH.reweight_pkl = "/eos/vbc/group/cms/robert.schoefbeck/gridpacks/VH/SMEFTsim_VH_reweight_card.pkl"
