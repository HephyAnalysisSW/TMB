
#/scratch-cbe/users/robert.schoefbeck/tWZ/nanoTuples/tWZ_nAODv6_private_v3/2016/dilep:
#DYJetsToLL_M10to50             DYJetsToLL_M50_HT100to200_ext   DYJetsToLL_M50_HT400to600_comb  DYJetsToLL_M50_LO_ext1_comb          DYJetsToLL_M5to50_HT400to600_comb  TTTT
#DYJetsToLL_M10to50_LO          DYJetsToLL_M50_HT1200to2500     DYJetsToLL_M50_HT600to800       DYJetsToLL_M50_LO_ext1_comb_lheHT70  DYJetsToLL_M5to50_HT600toInf       TTWW
#DYJetsToLL_M10to50_LO_lheHT70  DYJetsToLL_M50_HT200to400_comb  DYJetsToLL_M50_HT70to100        DYJetsToLL_M5to50_HT100to200_comb    DYJetsToLL_M5to50_HT70to100        TTWZ
#DYJetsToLL_M50_ext2            DYJetsToLL_M50_HT2500toInf      DYJetsToLL_M50_HT800to1200      DYJetsToLL_M5to50_HT200to400_comb    TTLep_pow                          TTZZ
#
#/scratch-cbe/users/robert.schoefbeck/tWZ/nanoTuples/tWZ_nAODv6_private_v3/2017/dilep:
#DYJetsToLL_M10to50_LO           DYJetsToLL_M4to50_HT200to400_comb  DYJetsToLL_M50_HT100to200       DYJetsToLL_M50_HT2500toInf  DYJetsToLL_M50_HT70to100    DYJetsToLL_M50_LO_comb_lheHT70  TTWW
#DYJetsToLL_M10to50_LO_lheHT100  DYJetsToLL_M4to50_HT400to600_comb  DYJetsToLL_M50_HT1200to2500     DYJetsToLL_M50_HT400to600   DYJetsToLL_M50_HT800to1200  TTLep_pow                       TTWZ
#DYJetsToLL_M4to50_HT100to200    DYJetsToLL_M4to50_HT600toInf_comb  DYJetsToLL_M50_HT200to400_comb  DYJetsToLL_M50_HT600to800   DYJetsToLL_M50_LO_comb      TTTT                            TTZZ



#DYJetsToLL_M10to50_LO
#DYJetsToLL_M50_HT1200to2500
#DYJetsToLL_M50_HT2500toInf
#DYJetsToLL_M50_HT600to800
#DYJetsToLL_M50_HT800to1200
#DYJetsToLL_M50_HT100to200
#DYJetsToLL_M50_HT200to400
#DYJetsToLL_M50_HT400to600
#DYJetsToLL_M50_HT70to100
#DYJetsToLL_M50_LO_lheHT70

import os
from RootTools.core.Sample import Sample
import ROOT

# Logging
import logging
logger = logging.getLogger(__name__)

from TMB.Samples.color import color

directory_ = "/scratch-cbe/users/robert.schoefbeck/tWZ/nanoTuples/tWZ_nAODv6_private_v4/2018/dilep"

def make_dirs( dirs ):
    return [ os.path.join( directory_, dir_ ) for dir_ in dirs ]

dirs = {}

dirs['TTLep']  = ["TTLep_pow"]
TTLep = Sample.fromDirectory(name="TTLep", treeName="Events", isData=False, color=color.TT, texName="TTLep", directory=make_dirs( dirs['TTLep'] ))

dirs['TTH']  = ["TTHbbLep", "TTHnobb_pow"]
TTH = Sample.fromDirectory(name="TTH", treeName="Events", isData=False, color=color.TTH, texName="TTH", directory=make_dirs( dirs['TTH'] ))

dirs['TTW']  = ["TTWToLNu", "TTWToQQ"]
TTW = Sample.fromDirectory(name="TTW", treeName="Events", isData=False, color=color.TTW, texName="TTW", directory=make_dirs( dirs['TTW'] ))

dirs['TTZ']  = ["TTZToLLNuNu", "TTZToLLNuNu_m1to10", "TTZToQQ"]
TTZ = Sample.fromDirectory(name="TTZ", treeName="Events", isData=False, color=color.TTZ, texName="TTZ", directory=make_dirs( dirs['TTZ'] ))

dirs['DYJetsToLL']  = ["DYJetsToLL_M50_LO"]
DYJetsToLL = Sample.fromDirectory(name="DYJetsToLL", treeName="Events", isData=False, color=color.DY, texName="DYJetsToLL", directory=make_dirs( dirs['DYJetsToLL'] ))

dirs['TTTT']  = ["TTTT"]
TTTT = Sample.fromDirectory(name="TTTT", treeName="Events", isData=False, color=color.TTTT, texName="TTTT", directory=make_dirs( dirs['TTTT'] ))

dirs['TTWZ']  = ["TTWZ"]
TTWZ = Sample.fromDirectory(name="TTWZ", treeName="Events", isData=False, color=color.TTWZ, texName="TTWZ", directory=make_dirs( dirs['TTWZ'] ))

dirs['TTWW']  = ["TTWW"]
TTWW = Sample.fromDirectory(name="TTWW", treeName="Events", isData=False, color=color.TTWW, texName="TTWW", directory=make_dirs( dirs['TTWW'] ))

dirs['TTZZ']  = ["TTZZ"]
TTZZ = Sample.fromDirectory(name="TTZZ", treeName="Events", isData=False, color=color.TTZZ, texName="TTZZ", directory=make_dirs( dirs['TTZZ'] ))
