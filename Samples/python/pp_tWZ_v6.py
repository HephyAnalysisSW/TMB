# Standard Imports
import os, sys
import ROOT

# RootTools Imports
from RootTools.core.Sample import Sample

# Colors
from TMB.Samples.color import color

# Logging
if __name__=="__main__":
    import Analysis.Tools.logger as logger
    logger = logger.get_logger("INFO", logFile = None )
    import RootTools.core.logger as logger_rt
    logger_rt = logger_rt.get_logger("INFO", logFile = None )
else:
    import logging
    logger = logging.getLogger(__name__)

data_path           = "/scratch-cbe/users/robert.schoefbeck/tWZ/nanoTuples/tWZ_nAODv6_private_v6/2018/" 

ttG_noFullyHad_fast = Sample.fromDirectory("ttG_noFullyHad_fast", directory = os.path.join( data_path, "singlelep-photon", "ttG_noFullyHad_fast") )
ttG_noFullyHad_fast.reweight_pkl = "/eos/vbc/user/robert.schoefbeck/gridpacks/v6/ttG_noFullyHad_reweight_card.pkl"
#ttZ01j_fast         = Sample.fromDirectory("ttZ01j_fast", directory = os.path.join( data_path,"ttZ01j_fast") )
#ttZ01j_fast.reweight_pkl = "/eos/vbc/user/robert.schoefbeck/gridpacks/v5/top_boson_reweight_card.pkl"
#ttW01j_fast         = Sample.fromDirectory("ttW01j_fast", directory = os.path.join( data_path,"ttW01j_fast") )
#ttW01j_fast.reweight_pkl = "/eos/vbc/user/robert.schoefbeck/gridpacks/v5/top_boson_reweight_card.pkl"
#WZTo3L1Nu_fast      = Sample.fromDirectory("WZTo3L1Nu_fast", directory = os.path.join( data_path,"WZTo3L1Nu_fast") )
#WZTo3L1Nu_fast.reweight_pkl = "/eos/vbc/user/robert.schoefbeck/gridpacks/v5/boson_reweight_card.pkl"
#WZTojj2L_fast       = Sample.fromDirectory("WZTojj2L_fast", directory = os.path.join( data_path,"WZTojj2L_fast") )
#WZTojj2L_fast.reweight_pkl = "/eos/vbc/user/robert.schoefbeck/gridpacks/v5/boson_reweight_card.pkl"
#WZToLNujj_fast      = Sample.fromDirectory("WZToLNujj_fast", directory = os.path.join( data_path,"WZToLNujj_fast") )
#WZToLNujj_fast.reweight_pkl = "/eos/vbc/user/robert.schoefbeck/gridpacks/v5/boson_reweight_card.pkl"
WGToLNu_fast        = Sample.fromDirectory("WGToLNu_fast", directory = os.path.join( data_path, "singlelep-photon", "WGToLNu_fast") )
WGToLNu_fast.reweight_pkl = "/eos/vbc/user/robert.schoefbeck/gridpacks/v6/WGToLNu_reweight_card.pkl"
#ZGTo2L_fast         = Sample.fromDirectory("ZGTo2L_fast", directory = os.path.join( data_path, "dilep-photon", "ZGTo2L_fast") )
#ZGTo2L_fast.reweight_pkl = "/eos/vbc/user/robert.schoefbeck/gridpacks/v5/boson_reweight_card.pkl"
#WW_fast             = Sample.fromDirectory("WW_fast", directory = os.path.join( data_path,"WW_fast") )
#WW_fast.reweight_pkl     = "/eos/vbc/user/robert.schoefbeck/gridpacks/v5/boson_reweight_card.pkl"

allSamples = [
    ttG_noFullyHad_fast,
#    ttZ01j_fast,
#    ttW01j_fast,
#    WZTo3L1Nu_fast,
#    WZTojj2L_fast,
#    WZToLNujj_fast,
    WGToLNu_fast,
#    ZGTo2L_fast,
#    WW_fast,

]

