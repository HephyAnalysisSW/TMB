''' miniAOD samples for TMB'''

# standard imports
import os

# RootTools
from RootTools.core.standard import *

# TMB
from TMB.Tools.user import cache_directory, gridpack_directory 

# sqlite3 sample cache file
dbFile = os.path.join( cache_directory, 'sample_cache', 'gen_TMB_miniAOD.db')
ov = False

# Logging
if __name__ == "__main__":
    import TMB.Tools.logger as logger
    logger = logger.get_logger('DEBUG')
    import RootTools.core.logger as logger_rt
    logger_rt = logger_rt.get_logger('DEBUG')

from Samples.Tools.config import redirector, redirector_global

DYBBJetsToLL_M50_LO                 = FWLiteSample.fromDAS("DYBBJetsToLL_M50_LO",         "/DYBBJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM", dbFile=dbFile, overwrite=ov, prefix=redirector_global, skipCheck=True)
DYBBJetsToLL_M50_LO.xsec            = 14.49*1.23
DYBBJetsToLL_M50_LO.nEvents         = 5039926 

DYJetsToLL_M50_HT70to100_LO         = FWLiteSample.fromDAS("DYJetsToLL_M50_HT70to100_LO",         "/DYJetsToLL_M-50_HT-70to100_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM", dbFile=dbFile, overwrite=ov, prefix=redirector_global, skipCheck=True)
DYJetsToLL_M50_HT70to100_LO.xsec    = 146.5*1.23
DYJetsToLL_M50_HT70to100_LO.nEvents = 10019684

DYJetsToLL_M50_HT100to200_LO        = FWLiteSample.fromDAS("DYJetsToLL_M50_HT100to200_LO",        "/DYJetsToLL_M-50_HT-100to200_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM", dbFile=dbFile, overwrite=ov, prefix=redirector_global, skipCheck=True)
DYJetsToLL_M50_HT100to200_LO.xsec   = 160.1*1.23
DYJetsToLL_M50_HT100to200_LO.nEvents = 11530510 

DYJetsToLL_M50_HT200to400_LO        = FWLiteSample.fromDAS("DYJetsToLL_M50_HT200to400_LO",        "/DYJetsToLL_M-50_HT-200to400_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM", dbFile=dbFile, overwrite=ov, prefix=redirector_global, skipCheck=True)
DYJetsToLL_M50_HT200to400_LO.xsec   = 49.22*1.23
DYJetsToLL_M50_HT200to400_LO.nEvents = 11225887

DYJetsToLL_M50_HT400to600_LO        = FWLiteSample.fromDAS("DYJetsToLL_M50_HT400to600_LO",        "/DYJetsToLL_M-50_HT-400to600_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v7/MINIAODSIM", dbFile=dbFile, overwrite=ov, prefix=redirector_global, skipCheck=True)
DYJetsToLL_M50_HT400to600_LO.xsec   = 6.996*1.23
DYJetsToLL_M50_HT400to600_LO.nEvents = 9697098 

DYJetsToLL_M50_HT400to600_LO_ext2   = FWLiteSample.fromDAS("DYJetsToLL_M50_HT400to600_LO_ext2",   "/DYJetsToLL_M-50_HT-400to600_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext2-v3/MINIAODSIM", dbFile=dbFile, overwrite=ov, prefix=redirector_global, skipCheck=True)
DYJetsToLL_M50_HT400to600_LO_ext2.xsec = 6.996*1.23
DYJetsToLL_M50_HT400to600_LO_ext2.nEvents = 9840466

DYJetsToLL_M50_HT600to800_LO        = FWLiteSample.fromDAS("DYJetsToLL_M50_HT600to800_LO",        "/DYJetsToLL_M-50_HT-600to800_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM", dbFile=dbFile, overwrite=ov, prefix=redirector_global, skipCheck=True)
DYJetsToLL_M50_HT600to800_LO.xsec   = 1.754*1.23
DYJetsToLL_M50_HT600to800_LO.nEvents = 8862104

DYJetsToLL_M50_HT800to1200_LO       = FWLiteSample.fromDAS("DYJetsToLL_M50_HT800to1200_LO",       "/DYJetsToLL_M-50_HT-800to1200_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM", dbFile=dbFile, overwrite=ov, prefix=redirector_global, skipCheck=True)
DYJetsToLL_M50_HT800to1200_LO.xsec  = 0.8718*1.23
DYJetsToLL_M50_HT800to1200_LO.nEvents = 3138129

DYJetsToLL_M50_HT1200to2500_LO      = FWLiteSample.fromDAS("DYJetsToLL_M50_HT1200to2500_LO",      "/DYJetsToLL_M-50_HT-1200to2500_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM", dbFile=dbFile, overwrite=ov, prefix=redirector_global, skipCheck=True)
DYJetsToLL_M50_HT1200to2500_LO.xsec = 0.1938*1.23
DYJetsToLL_M50_HT1200to2500_LO.nEvents = 536416

DYJetsToLL_M50_HT2500toInf_LO       = FWLiteSample.fromDAS("DYJetsToLL_M50_HT2500toInf_LO",       "/DYJetsToLL_M-50_HT-2500toInf_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM", dbFile=dbFile, overwrite=ov, prefix=redirector_global, skipCheck=True)
DYJetsToLL_M50_HT2500toInf_LO.xsec  = 0.003511*1.23
DYJetsToLL_M50_HT2500toInf_LO.nEvents = 427051

WJetsToLNu_HT_70To100   = FWLiteSample.fromDAS("WJetsToLNu_HT_70To100", "/WJetsToLNu_HT-70To100_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM", dbFile=dbFile, overwrite=ov, prefix=redirector_global, skipCheck=True)
WJetsToLNu_HT_70To100.nEvents = 28084244
WJetsToLNu_HT_70To100.xsec    = 1264.0

WJetsToLNu_HT_100To200 = FWLiteSample.fromDAS("WJetsToLNu_HT_100To200", "/WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM", dbFile=dbFile, overwrite=ov, prefix=redirector_global, skipCheck=True)
WJetsToLNu_HT_100To200.nEvents = 29521158
WJetsToLNu_HT_100To200.xsec    = 1256.0

WJetsToLNu_HT_200To400 = FWLiteSample.fromDAS("WJetsToLNu_HT_200To400", "/WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM", dbFile=dbFile, overwrite=ov, prefix=redirector_global, skipCheck=True)
WJetsToLNu_HT_200To400.nEvents = 25468933
WJetsToLNu_HT_200To400.xsec    = 335.5  

WJetsToLNu_HT_400To600 = FWLiteSample.fromDAS("WJetsToLNu_HT_400To600", "/WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM", dbFile=dbFile, overwrite=ov, prefix=redirector_global, skipCheck=True)
WJetsToLNu_HT_400To600.nEvents = 5932701
WJetsToLNu_HT_400To600.xsec    = 45.25

WJetsToLNu_HT_600To800 = FWLiteSample.fromDAS("WJetsToLNu_HT_600To800", "/WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM", dbFile=dbFile, overwrite=ov, prefix=redirector_global, skipCheck=True)
WJetsToLNu_HT_600To800.nEvents = 19771294
WJetsToLNu_HT_600To800.xsec    = 10.97

WJetsToLNu_HT_800To1200 = FWLiteSample.fromDAS("WJetsToLNu_HT_800To1200", "/WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM", dbFile=dbFile, overwrite=ov, prefix=redirector_global, skipCheck=True)
WJetsToLNu_HT_800To1200.nEvents = 8402687
WJetsToLNu_HT_800To1200.xsec    = 4.933

WJetsToLNu_HT_1200To2500 = FWLiteSample.fromDAS("WJetsToLNu_HT_1200To2500", "/WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM", dbFile=dbFile, overwrite=ov, prefix=redirector_global, skipCheck=True)
WJetsToLNu_HT_1200To2500.nEvents = 7633949
WJetsToLNu_HT_1200To2500.xsec    = 1.16

WJetsToLNu_HT_2500ToInf = FWLiteSample.fromDAS("WJetsToLNu_HT_2500ToInf", "/WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM", dbFile=dbFile, overwrite=ov, prefix=redirector_global, skipCheck=True)
WJetsToLNu_HT_2500ToInf.nEvents = 3273980
WJetsToLNu_HT_2500ToInf.xsec    = 0.008001

TTJets_DiLept = FWLiteSample.fromDAS("TTJets_DiLept", "/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM", dbFile=dbFile, overwrite=ov, prefix=redirector_global, skipCheck=True)
TTJets_DiLept.nEvents = 28701360
TTJets_DiLept.xsec    = 54.17
# 841

TTJets_DiLept_genMET80 = FWLiteSample.fromDAS("TTJets_DiLept_genMET80", "/TTJets_DiLept_genMET-80_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM", dbFile=dbFile, overwrite=ov, prefix=redirector_global, skipCheck=True)
TTJets_DiLept_genMET80.nEvents = 58442000
TTJets_DiLept_genMET80.xsec    = 22.45   
# 1578

TTJets_HT_1200to2500 = FWLiteSample.fromDAS("TTJets_HT_1200to2500", "/TTJets_HT-1200to2500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM", dbFile=dbFile, overwrite=ov, prefix=redirector_global, skipCheck=True)
TTJets_HT_1200to2500.nEvents = 2779427
TTJets_HT_1200to2500.xsec    = 0.1316
# 123

TTJets_HT_2500toInf = FWLiteSample.fromDAS("TTJets_HT_2500toInf", "/TTJets_HT-2500toInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM", dbFile=dbFile, overwrite=ov, prefix=redirector_global, skipCheck=True)
TTJets_HT_2500toInf.nEvents = 1451104
TTJets_HT_2500toInf.xsec    = 0.001407
# 63

TTJets_HT_600to800 = FWLiteSample.fromDAS("TTJets_HT_600to800", "/TTJets_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM", dbFile=dbFile, overwrite=ov, prefix=redirector_global, skipCheck=True)
TTJets_HT_600to800.nEvents = 14149394
TTJets_HT_600to800.xsec    = 1.821
# 503

TTJets_HT_800to1200 = FWLiteSample.fromDAS("TTJets_HT_800to1200", "/TTJets_HT-800to1200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM", dbFile=dbFile, overwrite=ov, prefix=redirector_global, skipCheck=True)
TTJets_HT_800to1200.nEvents = 10372802
TTJets_HT_800to1200.xsec    = 0.7532
# 399

TTJets_SingleLeptFromT = FWLiteSample.fromDAS("TTJets_SingleLeptFromT", "/TTJets_SingleLeptFromT_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM", dbFile=dbFile, overwrite=ov, prefix=redirector_global, skipCheck=True)
TTJets_SingleLeptFromT.nEvents = 57259880
TTJets_SingleLeptFromT.xsec    = 109.6
# 1505

TTJets_SingleLeptFromT_genMET80 = FWLiteSample.fromDAS("TTJets_SingleLeptFromT_genMET80", "/TTJets_SingleLeptFromT_genMET-80_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM", dbFile=dbFile, overwrite=ov, prefix=redirector_global, skipCheck=True)
TTJets_SingleLeptFromT_genMET80.nEvents = 84200008
TTJets_SingleLeptFromT_genMET80.xsec    = 32.27
# 2230

TTJets_SingleLeptFromTbar = FWLiteSample.fromDAS("TTJets_SingleLeptFromTbar", "/TTJets_SingleLeptFromTbar_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM", dbFile=dbFile, overwrite=ov, prefix=redirector_global, skipCheck=True)
TTJets_SingleLeptFromTbar.nEvents = 59929205
TTJets_SingleLeptFromTbar.xsec    = 108.7
# 1540

TTJets_SingleLeptFromTbar_genMET80 = FWLiteSample.fromDAS("TTJets_SingleLeptFromTbar_genMET80", "/TTJets_SingleLeptFromTbar_genMET-80_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM", dbFile=dbFile, overwrite=ov, prefix=redirector_global, skipCheck=True)
TTJets_SingleLeptFromTbar_genMET80.nEvents = 79533863
TTJets_SingleLeptFromTbar_genMET80.xsec    = 31.68
# 2055

TTJets = FWLiteSample.fromDAS("TTJets", "/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM", dbFile=dbFile, overwrite=ov, prefix=redirector_global, skipCheck=True)
TTJets.nEvents = 10244307
TTJets.xsec    = 496.1
# 309
