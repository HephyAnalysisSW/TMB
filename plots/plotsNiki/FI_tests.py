#!/usr/bin/env python
''' Analysis script for standard plots
'''
#
# Standard imports and batch mode
#
import ROOT, os
ROOT.gROOT.SetBatch(True)
c1 = ROOT.TCanvas() # do this to avoid version conflict in png.h with keras import ...
c1.Draw()
c1.Print('/tmp/delete.png')

import itertools
import copy
import array
import operator
from math                                import sqrt, cos, sin, pi, atan2, cosh, log

# RootTools
from RootTools.core.standard             import *

# tWZ
from TMB.Tools.user                      import plot_directory
from TMB.Tools.cutInterpreter            import cutInterpreter
from tWZ.Tools.objectSelection           import lepString # probably will merge TMB and tWZ repos 
# Analysis
from Analysis.Tools.WeightInfo                import WeightInfo
import Analysis.Tools.syncer             as     syncer
import numpy as np

# Arguments
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',       action='store',      default='INFO', nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'], help="Log level for logging")
argParser.add_argument('--small',                             action='store_true', help='Run only on a small subset of the data?', )
#argParser.add_argument('--sorting',                           action='store', default=None, choices=[None, "forDYMB"],  help='Sort histos?', )
argParser.add_argument('--plot_directory', action='store', default='FI-test-v6')
argParser.add_argument('--sample',        action='store', type=str, default="WGToLNu_fast")
argParser.add_argument('--selection',      action='store', default='singlelep-photon')
args = argParser.parse_args()

# Logger
import TMB.Tools.logger as logger
import RootTools.core.logger as logger_rt
logger    = logger.get_logger(   args.logLevel, logFile = None)
logger_rt = logger_rt.get_logger(args.logLevel, logFile = None)


import TMB.Samples.pp_tWZ_v6 as samples

sample = getattr( samples, args.sample )

lumi_scale = 137

#lumi_weight = lambda event, sample: lumi_scale*event.weight 

def lumi_weight( event, sample):
    return 1.# lumi_scale*event.weight

if args.small:
    sample.reduceFiles( to = 1 )
    args.plot_directory += "_small"

# WeightInfo
w = WeightInfo(sample.reweight_pkl)
w.set_order(2)

# selection
selectionString = cutInterpreter.cutString(args.selection)

# Read variables and sequences

read_variables = [
    "weight/F", "year/I", "met_pt/F", "met_phi/F", "nBTag/I", "nJetGood/I", "PV_npvsGood/I",
    "l1_pt/F", "l1_eta/F" , "l1_phi/F", "l1_mvaTOP/F", "l1_mvaTOPWP/I", "l1_index/I", 
    #"l2_pt/F", "l2_eta/F" , "l2_phi/F", "l2_mvaTOP/F", "l2_mvaTOPWP/I", "l2_index/I",
    #"l3_pt/F", "l3_eta/F" , "l3_phi/F", "l3_mvaTOP/F", "l3_mvaTOPWP/I", "l3_index/I",
    "JetGood[pt/F,eta/F,phi/F]",
    "lep[pt/F,eta/F,phi/F,pdgId/I,muIndex/I,eleIndex/I]",
#    "Z1_l1_index/I", "Z1_l2_index/I", "nonZ1_l1_index/I", "nonZ1_l2_index/I", 
    "Z1_phi/F", "Z1_pt/F", "Z1_mass/F", "Z1_cosThetaStar/F", "Z1_eta/F", "Z1_lldPhi/F", "Z1_lldR/F",
    "Muon[pt/F,eta/F,phi/F,dxy/F,dz/F,ip3d/F,sip3d/F,jetRelIso/F,miniPFRelIso_all/F,pfRelIso03_all/F,mvaTOP/F,mvaTTH/F,pdgId/I,segmentComp/F,nStations/I,nTrackerLayers/I]",
    "Electron[pt/F,eta/F,phi/F,dxy/F,dz/F,ip3d/F,sip3d/F,jetRelIso/F,miniPFRelIso_all/F,pfRelIso03_all/F,mvaTOP/F,mvaTTH/F,pdgId/I,vidNestedWPBitmap/I]",
    "np/I", VectorTreeVariable.fromString("p[C/F]", nMax=500),
    "photon_pt/F",
    "photon_eta/F",
    "photon_phi/F",
    "photonJetdR/F", "photonLepdR/F",
]

WC    = 'cWWW'
WCval = 1.0

# string to compute Fisher information. 
FI_string = w.get_fisher_weight_string(WC,WC, **{WC:WCval})
#Try: sample.chain.Scan('photon_pt:'+FI_string, "photon_pt>400") 

# Get the coefficient list for a selection
coeffList = w.getCoeffListFromDraw( sample, 'photon_pt>400' ) 
# compute FI from coefficient list:
print ("FI for %s and %r"%(WC, {WC:WCval}), w.get_fisherInformation_matrix(coeffList, [WC], **{WC:WCval}))

# read coefficients of individual events:
coeffList_events = w.getCoeffListFromEvents( sample, 'photon_pt>400' )
# accumulate events, then compute FI
print w.get_fisherInformation_matrix( map(sum, zip(*coeffList_events)), [WC], **{WC:WCval} )
# (watch out ... compute per-event FI)
for coeffList in coeffList_events:
    print w.get_fisherInformation_matrix( coeffList, [WC], **{WC:WCval} )

