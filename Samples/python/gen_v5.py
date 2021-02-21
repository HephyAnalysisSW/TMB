''' GEN samples for TMB'''

# standard imports
import os

# RootTools
from RootTools.core.standard import *

# TMB
from TMB.Tools.user import cache_directory, gridpack_directory 

# sqlite3 sample cache file
dbFile = os.path.join( cache_directory, 'sample_cache', 'gen_v5.db')
overwrite = False

# Logging
if __name__ == "__main__":
    import TMB.Tools.logger as logger
    logger = logger.get_logger('DEBUG')
    import RootTools.core.logger as logger_rt
    logger_rt = logger_rt.get_logger('DEBUG')

# VV and VVV with boson reweighting
# Only leptonic decays for all V!!

#test = FWLiteSample.fromFiles("test", ["/eos/vbc/experiments/cms/store/user/schoef/WW-test/WW-test/201223_221102/0000/GEN_LO_0j_93X_1.root"] )
#test.reweight_pkl = os.path.join(gridpack_directory, "v4", "boson_v1_reweight_card.pkl")
#test.xsec         = 1 
#test.nEvents      = 100 
#
#test2 = FWLiteSample.fromFiles("test2", ["/users/robert.schoefbeck/CMS/test/CMSSW_10_2_22/src/reco_FastSim_LO_0j_102X_CP5_FastSim.root "] )
#test2.reweight_pkl = os.path.join(gridpack_directory, "v4", "boson_v1_reweight_card.pkl")
#test2.xsec         = 1 
#test2.nEvents      = 100 
#
#WW  = FWLiteSample.fromDAS("WW", "/WW-v4/schoef-WW-v4-238a1f3c8105c56c183c394f2927064f/USER", "phys03", dbFile = dbFile, overwrite=overwrite)
#WW.reweight_pkl = os.path.join(gridpack_directory, "v4", "boson_v1_reweight_card.pkl")
#WW.xsec         = 1.751 
#WW.nEvents      = 5118391
#
WZTo3L1Nu  = FWLiteSample.fromDAS("WZTo3L1Nu", "/WZTo3L1Nu-v5/schoef-WZTo3L1Nu-v5-151ec2871e2ed107def2a3aa630d2b24/USER", "phys03", dbFile = dbFile, overwrite=overwrite, skipCheck=True)
WZTo3L1Nu.reweight_pkl = os.path.join(gridpack_directory, "v5", "boson_reweight_card.pkl")
WZTo3L1Nu.xsec         = 1.729e+00 
WZTo3L1Nu.nEvents      = 1547385.0
#
#ZZ  = FWLiteSample.fromDAS("ZZ", "/ZZ-v4/schoef-ZZ-v4-f26891a5e446d02810e35e2cfeda4b23/USER", "phys03", dbFile = dbFile, overwrite=overwrite)
#ZZ.reweight_pkl = os.path.join(gridpack_directory, "v4", "boson_v1_reweight_card.pkl")
#ZZ.xsec         = 110.74415084
#ZZ.nEvents      = 5144574
#
#ZA  = FWLiteSample.fromDAS("ZA", "/ZA-v4/schoef-ZA-v4-c271581693d6eda368761f0699a46798/USER", "phys03", dbFile = dbFile, overwrite=overwrite)
#ZA.reweight_pkl = os.path.join(gridpack_directory, "v4", "boson_v1_reweight_card.pkl")
#ZA.xsec         = 5.794e+00
#ZA.nEvents      = 4361967
#
#WWW = FWLiteSample.fromDAS("WWW", "/WWW-v4/schoef-WWW-v4-2d24f84d57d61aa0ed9149b9c9610b37/USER", "phys03", dbFile = dbFile, overwrite=overwrite)
#WWW.reweight_pkl = os.path.join(gridpack_directory, "v4", "boson_v1_reweight_card.pkl")
#WWW.xsec         = 0.01834 
#WWW.nEvents      = 4130664 
#
#WWZ = FWLiteSample.fromDAS("WWZ", "/WWZ-v4/schoef-WWZ-v4-d2861d31cc6c25778d9def6ca333eb5d/USER", "phys03", dbFile = dbFile, overwrite=overwrite)
#WWZ.reweight_pkl = os.path.join(gridpack_directory, "v4", "boson_v1_reweight_card.pkl")
#WWZ.xsec         = 0.004744 
#WWZ.nEvents      = 4205193
#
#WZZ = FWLiteSample.fromDAS("WZZ", "/WZZ-v4/schoef-WZZ-v4-808cecf77e16e2367f66c7edcf59b94d/USER", "phys03", dbFile = dbFile, overwrite=overwrite)
#WZZ.reweight_pkl = os.path.join(gridpack_directory, "v4", "boson_v1_reweight_card.pkl")
#WZZ.xsec         = 0.0003355 
#WZZ.nEvents      = 4399466
#
#ZZZ = FWLiteSample.fromDAS("ZZZ", "/ZZZ-v4/schoef-ZZZ-v4-fdf8e8138c23af78baac6d1fc1ca4a27/USER", "phys03", dbFile = dbFile, overwrite=overwrite)
#ZZZ.reweight_pkl = os.path.join(gridpack_directory, "v4", "boson_v1_reweight_card.pkl")
#ZZZ.xsec         = 0.00003884 
#ZZZ.nEvents      = 4272023
#
#WAjj = FWLiteSample.fromDAS("WAjj", "/WAjj-v4/schoef-WAjj-v4-e218d7828d2ff40ab0212e5db8c3c3d7/USER", "phys03", dbFile = dbFile, overwrite=overwrite)
#WAjj.reweight_pkl = os.path.join(gridpack_directory, "v4", "boson_v1_reweight_card.pkl")
#WAjj.xsec       = 3.872e-01
#WAjj.nEvents    = 10000000
#
WGToLNu = FWLiteSample.fromDAS("WGToLNu", "/WGToLNu/schoef-WGToLNu-888b7a86c2f3c15fead55bb8986384d5/USER", "phys03", dbFile = dbFile, overwrite=overwrite, skipCheck=True)
WGToLNu.reweight_pkl = os.path.join(gridpack_directory, "v5", "boson_reweight_card.pkl")
WGToLNu.xsec     = 5.453e+01
WGToLNu.nEvents  = 1515120
#
#WA_LO = FWLiteSample.fromDAS("WA_LO", "/WA_LO-v4/schoef-WA_LO-v4-7615e614e93751a66f1db38b71b2c3c8/USER", "phys03", dbFile = dbFile, overwrite=overwrite)
#WA_LO.reweight_pkl = os.path.join(gridpack_directory, "v4", "boson_v1_reweight_card.pkl")
#WA_LO.xsec    = 2.352e+01
#WA_LO.nEvents = 9850000
#
#WWjj_OS = FWLiteSample.fromDAS("WWjj_OS", "/WWjj_OS-v4/schoef-WWjj_OS-v4-f92080c6914bdc8d26ba686c331a8c58/USER", "phys03", dbFile = dbFile, overwrite=overwrite)
#WWjj_OS.reweight_pkl = os.path.join(gridpack_directory, "v4", "boson_v1_reweight_card.pkl")
#WWjj_OS.xsec    = 1.327e+00
#WWjj_OS.nEvents = 8290000
#
#WWjj_SS = FWLiteSample.fromDAS("WWjj_SS", "/WWjj_SS-v4/schoef-WWjj_SS-v4-d4f9a470ff7eca744278e830c3726bfc/USER", "phys03", dbFile = dbFile, overwrite=overwrite)
#WWjj_SS.reweight_pkl = os.path.join(gridpack_directory, "v4", "boson_v1_reweight_card.pkl")
#WWjj_SS.xsec    = 3.720e-02
#WWjj_SS.nEvents = 9995000
#
#WZjj = FWLiteSample.fromDAS("WZjj", "/WZjj-v4/schoef-WZjj-v4-11f9d5415a2912befe8a7b02fdf13990/USER", "phys03", dbFile = dbFile, overwrite=overwrite)
#WZjj.reweight_pkl = os.path.join(gridpack_directory, "v4", "boson_v1_reweight_card.pkl")
#WZjj.xsec    = 1.863e-02
#WZjj.nEvents = 9995000
#
#ZAjj = FWLiteSample.fromDAS("ZAjj", "/ZAjj-v4/schoef-ZAjj-v4-60aadee52814bd0f27220d33d89e99f3/USER", "phys03", dbFile = dbFile, overwrite=overwrite)
#ZAjj.reweight_pkl = os.path.join(gridpack_directory, "v4", "boson_v1_reweight_card.pkl")
#ZAjj.xsec    = 3.671e-02
#ZAjj.nEvents = 10000000
#
#ZZjj = FWLiteSample.fromDAS("ZZjj", "/ZZjj-v4/schoef-ZZjj-v4-8fdc3294d133ba1fc99a45ab33fb5565/USER", "phys03", dbFile = dbFile, overwrite=overwrite)
#ZZjj.reweight_pkl = os.path.join(gridpack_directory, "v4", "boson_v1_reweight_card.pkl")
#ZZjj.xsec    = 1.435e-03
#ZZjj.nEvents = 10000000

ttG_noFullyHad = FWLiteSample.fromDAS("ttG_noFullyHad", "/ttG_noFullyHad-v5/schoef-ttG_noFullyHad-v5-fa75f6d0d64d97494f28d49b68b6ba8c/USER", "phys03", dbFile = dbFile, overwrite=overwrite, skipCheck=True)
ttG_noFullyHad.reweight_pkl = os.path.join(gridpack_directory, "v5", "top_boson_reweight_card.pkl")
ttG_noFullyHad.xsec    = 1.947
ttG_noFullyHad.nEvents = 3972000

#ttW01j = FWLiteSample.fromDAS("ttW01j", "/ttW01j-v4/schoef-ttW01j-v4-b005ba525445f45dfb6340de320c29c1/USER", "phys03", dbFile = dbFile, overwrite=overwrite)
#ttW01j.reweight_pkl = os.path.join(gridpack_directory, "v4", "top_boson_v1_reweight_card.pkl")
#ttW01j.xsec    = 1.052
#ttW01j.nEvents = 3734478
#
ttZ01j = FWLiteSample.fromDAS("ttZ01j", "/ttZ01j-v5/schoef-ttZ01j-v5-5941b9a9132a067bf58859850614a6bf/USER", "phys03", dbFile = dbFile, overwrite=overwrite, skipCheck=True)
ttZ01j.reweight_pkl = os.path.join(gridpack_directory, "v5", "top_boson_reweight_card.pkl")
ttZ01j.xsec    = 1.369e-01
ttZ01j.nEvents = 3345859
#
#WH = FWLiteSample.fromDAS("WH", "/WH-v1/chatterj-WH-v1-13982263318d64186461fb17c44d9cef/USER", "phys03", dbFile = dbFile, overwrite=overwrite)                                                              
#WH.reweight_pkl = os.path.join(gridpack_directory, "v4", "WH_LeptonicW_LO_reweight_card.pkl")
#WH.xsec    = 7.505
#WH.nEvents = 3017257
