from RootTools.core.standard import *

import TMB.Samples.nanoTuples_Summer16_nanoAODv6_dilep_pp as Summer16
import TMB.Samples.nanoTuples_Fall17_nanoAODv6_dilep_pp   as Fall17
import TMB.Samples.nanoTuples_Autumn18_nanoAODv6_dilep_pp as Autumn18

#DYJetsToLL  = Sample.combine( "DYJetsToLL", [Summer16.DYJetsToLL, Fall17.DYJetsToLL, Autumn18.DYJetsToLL])
TTLep       = Sample.combine( "TTLep", [Summer16.TTLep, Fall17.TTLep, Autumn18.TTLep])

TTTT       = Sample.combine( "TTTT", [Summer16.TTTT, Fall17.TTTT, Autumn18.TTTT])
TTWW       = Sample.combine( "TTWW", [Summer16.TTWW, Fall17.TTWW, Autumn18.TTWW])
TTWZ       = Sample.combine( "TTWZ", [Summer16.TTWZ, Fall17.TTWZ, Autumn18.TTWZ])
TTZZ       = Sample.combine( "TTZZ", [Summer16.TTZZ, Fall17.TTZZ, Autumn18.TTZZ])

TTW       = Sample.combine( "TTW", [Summer16.TTW, Fall17.TTW, Autumn18.TTW])
TTZ       = Sample.combine( "TTZ", [Summer16.TTZ, Fall17.TTZ, Autumn18.TTZ])
TTH       = Sample.combine( "TTH", [Summer16.TTH, Fall17.TTH, Autumn18.TTH])

TTWZH     = Sample.combine( "TTWZH", [Summer16.TTH, Fall17.TTH, Autumn18.TTH])
TTWZH.legendText = "TT+W/Z/H"

DY        = Summer16.DY
