#!/usr/bin/env python
''' Optimize 1D FI over single photon_pt cut
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

import time
from cut_opt import *

# Arguments
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',       action='store',      default='INFO', nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'], help="Log level for logging")
argParser.add_argument('--small',                             action='store_true', help='Run only on a small subset of the data?', )
#argParser.add_argument('--sorting',                           action='store', default=None, choices=[None, "forDYMB"],  help='Sort histos?', )
argParser.add_argument('--plot_directory', action='store', default='FI-test-v6')
argParser.add_argument('--sample',        action='store', type=str, default="WGToLNu_fast")
argParser.add_argument('--selection',      action='store', default='singlelep-photon')
argParser.add_argument('--first_photon_pt_cut',      action='store', type=float, default=400.0)
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

photon_pt_initial_cut = args.first_photon_pt_cut

t_start = time.time()

coeffs_with_photon_pt_list_events = get_coeff_list_with_photon_pt_from_events(sample, "photon_pt>%f" % photon_pt_initial_cut)
number_of_events = len(coeffs_with_photon_pt_list_events)

sorted_coeffs_with_photon_pt_list_events = sorted(coeffs_with_photon_pt_list_events, key=lambda x: x['photon_pt'])
sorted_coeffs_list_events = map(lambda x: x['coeffs'], sorted_coeffs_with_photon_pt_list_events)
sorted_coeffs_matrix = np.array(zip(*sorted_coeffs_list_events))

t_draw = time.time() - t_start
print "time needed for data selection: %f seconds" % t_draw

t_start = time.time()
fi_calc_time = 0.0

second_photon_pt_cut_with_fi = []
max_fi_sum = 0
max_photon_pt = photon_pt_initial_cut
total_fi_sum = w.get_fisherInformation_matrix(np.sum(sorted_coeffs_matrix, axis=1), [WC], **{WC:WCval})[1][0][0]

for i, event_data in enumerate(sorted_coeffs_with_photon_pt_list_events):
    fi_sum = 0
    # can be implemented more efficiently by incremental changes
    left_coeffs = np.sum(sorted_coeffs_matrix[:, 0:i+1], axis=1)
    right_coeffs = np.sum(sorted_coeffs_matrix[:, i+1:], axis=1)
    assert left_coeffs[0] + right_coeffs[0] == number_of_events
    t_start_2 = time.time()
    fi_sum += w.get_fisherInformation_matrix(left_coeffs, [WC], **{WC:WCval})[1][0][0]
    fi_sum += w.get_fisherInformation_matrix(right_coeffs, [WC], **{WC:WCval})[1][0][0]
    fi_calc_time += time.time() -t_start_2
    if fi_sum > max_fi_sum:
        max_photon_pt = event_data['photon_pt']
        max_fi_sum = fi_sum
    second_photon_pt_cut_with_fi.append([event_data['photon_pt'], fi_sum/total_fi_sum-1])

t_cut_optimization = time.time() - t_start

print "time needed for cut optimization: %f seconds" % t_cut_optimization
print "of this time, fi calc time: %f second" % fi_calc_time
print "total fi sum without cut %f, max fi sum: %f at second cut on photon_pt: %f, max relative fi increase in: %.02f%%" % (total_fi_sum, max_fi_sum, max_photon_pt, (max_fi_sum/total_fi_sum-1)*100)

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

second_photon_pt_cut_with_fi_matrix = np.array(second_photon_pt_cut_with_fi)

plt.plot(second_photon_pt_cut_with_fi_matrix[:, 0], second_photon_pt_cut_with_fi_matrix[:, 1]*100)
plt.grid(True)
plt.xlabel('second photon pt cut')
plt.ylabel('relative fisher information increase in %')
plt.savefig('FI_over_second_photon_pt_cut.pdf', dpi=1000)