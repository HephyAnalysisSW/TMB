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
from Analysis.Tools.helpers              import deltaPhi, deltaR
from Analysis.Tools.puProfileCache       import *
from Analysis.Tools.puReweighting        import getReweightingFunction
from TMB.Tools.WeightInfo                import WeightInfo
import Analysis.Tools.syncer             as     syncer
import numpy as np

# Arguments
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',       action='store',      default='INFO', nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'], help="Log level for logging")
argParser.add_argument('--small',                             action='store_true', help='Run only on a small subset of the data?', )
#argParser.add_argument('--sorting',                           action='store', default=None, choices=[None, "forDYMB"],  help='Sort histos?', )
argParser.add_argument('--plot_directory', action='store', default='FI-test')
argParser.add_argument('--WC',                 action='store',      default='cWWW')
argParser.add_argument('--WCval',              action='store',   type=float,    default=1.0,  help='Value of the Wilson coefficient for the distribution.')
argParser.add_argument('--era',            action='store', type=str, default="Autumn18")
argParser.add_argument('--sample',        action='store', type=str, default="WGToLNu_fast")
argParser.add_argument('--selection',      action='store', default='singlelep-photon')
args = argParser.parse_args()

# Logger
import TMB.Tools.logger as logger
import RootTools.core.logger as logger_rt
logger    = logger.get_logger(   args.logLevel, logFile = None)
logger_rt = logger_rt.get_logger(args.logLevel, logFile = None)


import TMB.Samples.pp_tWZ as samples

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

# define which Wilson coefficients to plot
FIs = []

param =  {'legendText':'%s %2.1f'%(args.WC, args.WCval), 'WC':{args.WC:args.WCval} } 

# selection
selectionString = cutInterpreter.cutString(args.selection)

# get quantiles of log10(FI)
h_FI = sample.get1DHistoFromDraw("TMath::Log10(%s)" %w.get_fisher_weight_string(args.WC,args.WC, **{args.WC:0}), binning = [500,-30,20], selectionString=selectionString)
for i_bin in range(1, h_FI.GetNbinsX()+1):
    h_FI.SetBinContent( i_bin, 10**h_FI.GetBinLowEdge(i_bin)*h_FI.GetBinContent(i_bin) )
quantiles = [0, 0.2, 0.4, 0.6, 0.8, 1.0]
prob  =   array.array('d', quantiles)
q     =   array.array('d', [0]*len(quantiles) )
h_FI.GetQuantiles(len(q),q, prob)

colors = [ ROOT.kViolet-9, ROOT.kViolet-6, ROOT.kViolet-4, ROOT.kViolet-3 , ROOT.kViolet-2, ROOT.kViolet]

var = w.get_fisher_weight_string(args.WC,args.WC, **{args.WC:args.WCval})

params = [{ 'legendText':'log_{10}(I_{F})<%2.1f'%(q[1]), 'selection':"{var}<%f"%(q[1]) , 'style':styles.fillStyle(colors[0])}]
for i_th, th in enumerate( q[1:-2] ):
    params.append( { 'legendText':'%2.1f<=log_{10}(I_{F})<%2.1f'%(q[i_th+1], q[i_th+2]), 'selection':"%f<={var}&&{var}<%f"%(q[i_th+1], q[i_th+2]) , 'style':styles.fillStyle(colors[i_th+1])})
params.append( { 'legendText':'%2.1f<=log_{10}(I_{F})'%(q[-2]), 'selection':"%f<={var}"%(q[-2]) , 'style':styles.fillStyle(colors[-1])})

stack = Stack([ copy.deepcopy(sample) for param in params ] )
for i_s_, s_ in enumerate(stack[0]):
    s_.setSelectionString( params[i_s_]['selection'].format(var= "TMath::Log10(%s)"%var) )
    s_.name+='_%i'%i_s_
    s_.weight = lumi_weight
weight = [ [ w.get_weight_func(**{args.WC:args.WCval}) for param in params ] ]

# Read variables and sequences
sequence       = []

if args.selection.count('dilep'):
    def make_dilep_angles(event, sample):

        event.photon_Z1_deltaR   = float('nan')
        event.photon_Z1_deltaPhi = float('nan')
        event.jet0_Z1_deltaR = float('nan')
        event.jet1_Z1_deltaR = float('nan')

        if event.Z1_pt>=0:
            event.photon_Z1_deltaR   = deltaR({'eta':event.photon_eta, 'phi':event.photon_phi}, {'eta':event.Z1_eta, 'phi':event.Z1_phi})
            event.photon_Z1_deltaPhi = deltaPhi(event.photon_phi, event.Z1_phi)
            if event.nJetGood>=1:
                event.jet0_Z1_deltaR      = deltaR({'eta':event.JetGood_eta[0], 'phi':event.JetGood_phi[0]}, {'eta':event.Z1_eta, 'phi':event.Z1_phi})
            if event.nJetGood>=1:
                event.jet1_Z1_deltaR      = deltaR({'eta':event.JetGood_eta[1], 'phi':event.JetGood_phi[1]}, {'eta':event.Z1_eta, 'phi':event.Z1_phi})

    sequence.append( make_dilep_angles )

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
    "np/I", "p[C/F]",
    "photon_pt/F",
    "photon_eta/F",
    "photon_phi/F",
    "photonJetdR/F", "photonLepdR/F",
]

read_variables_MC = ['reweightBTag_SF/F', 'reweightPU/F', 'reweightL1Prefire/F', 'reweightLeptonSF/F'] #'reweightTrigger/F']

yields     = {}
allPlots   = {}

#yt_TWZ_filter.scale = lumi_scale * 1.07314

# Use some defaults
Plot.setDefaults(stack = stack, weight = weight, selectionString = selectionString)

plots        = []

plots.append(Plot(
  name = 'nVtxs', texX = 'vertex multiplicity', texY = 'Number of Events',
  attribute = TreeVariable.fromString( "PV_npvsGood/I" ),
  binning=[50,0,50],
  addOverFlowBin='upper',
))

plots.append(Plot(
    name = 'M3',
    texX = 'M_{3} (GeV)', texY = 'Number of Events / 10 GeV',
    attribute = TreeVariable.fromString( "m3/F" ),
    binning=[30,0,300],
    addOverFlowBin='upper',
))

plots.append(Plot(
    name = 'log10FI',
    texX = 'log_{10}(FI(%s=%2.1f))' % (args.WC, args.WCval), texY = 'Number of Events',
    attribute = lambda event, sample: ROOT.TMath.Log10(w.get_fisherInformation_matrix([event.p_C[i] for i in range(event.np)], variables=[args.WC,args.WC], **{args.WC:args.WCval})[1][0][0]),
    binning=[300,-20,10],
    addOverFlowBin='upper',
))

plots.append(Plot(
    name = 'photonJetdR',
    texX = '#Delta R(#gamma, jets)', texY = 'Number of Events',
    attribute = TreeVariable.fromString( "photonJetdR/F" ),
    binning=[30,0,6],
    addOverFlowBin='upper',
))

plots.append(Plot(
    name = 'photonLepdR',
    texX = '#Delta R(#gamma, leptons)', texY = 'Number of Events',
    attribute = TreeVariable.fromString( "photonLepdR/F" ),
    binning=[30,0,6],
    addOverFlowBin='upper',
))

plots.append(Plot(
    name = 'mT',
    texX = 'M_{T} (GeV)', texY = 'Number of Events / 10 GeV',
    attribute = lambda event, sample: sqrt(2.*event.met_pt*event.l1_pt*(1-cos(event.met_phi-event.l1_phi))),
    binning=[30,0,300],
    addOverFlowBin='upper',
))

plots.append(Plot(
    name = 'photon_pt',
    texX = 'p_{T}(#gamma) (GeV)', texY = 'Number of Events / 20 GeV',
    attribute = lambda event, sample:event.photon_pt,
    binning=[15,0,300],
    addOverFlowBin='upper',
))

plots.append(Plot(
    name = 'l1_pt',
    texX = 'p_{T}(l_{1}) (GeV)', texY = 'Number of Events / 20 GeV',
    attribute = lambda event, sample:event.l1_pt,
    binning=[15,0,300],
    addOverFlowBin='upper',
))

plots.append(Plot(
    name = 'l1_eta',
    texX = '#eta(l_{1})', texY = 'Number of Events',
    attribute = lambda event, sample: event.l1_eta,
    binning=[20,-3,3],
))

plots.append(Plot(
    name = 'l1_mvaTOP',
    texX = 'MVA_{TOP}(l_{1})', texY = 'Number of Events',
    attribute = lambda event, sample: event.l1_mvaTOP,
    binning=[20,-1,1],
))

plots.append(Plot(
    name = 'l1_mvaTOPWP',
    texX = 'MVA_{TOP}(l_{1}) WP', texY = 'Number of Events',
    attribute = lambda event, sample: event.l1_mvaTOPWP,
    binning=[5,0,5],
))

plots.append(Plot(
    texX = 'E_{T}^{miss} (GeV)', texY = 'Number of Events / 20 GeV',
    attribute = TreeVariable.fromString( "met_pt/F" ),
    binning=[400/20,0,400],
    addOverFlowBin='upper',
))

plots.append(Plot(
  texX = 'N_{jets}', texY = 'Number of Events',
  attribute = TreeVariable.fromString( "nJetGood/I" ), #nJetSelected
  binning=[8,-0.5,7.5],
))

plots.append(Plot(
  texX = 'N_{b-tag}', texY = 'Number of Events',
  attribute = TreeVariable.fromString( "nBTag/I" ), #nJetSelected
  binning=[4,-0.5,3.5],
))

plots.append(Plot(
  texX = 'p_{T}(leading jet) (GeV)', texY = 'Number of Events / 30 GeV',
  name = 'jet0_pt', attribute = lambda event, sample: event.JetGood_pt[0],
  binning=[600/30,0,600],
))

plots.append(Plot(
  texX = 'p_{T}(subleading jet) (GeV)', texY = 'Number of Events / 30 GeV',
  name = 'jet1_pt', attribute = lambda event, sample: event.JetGood_pt[1],
  binning=[600/30,0,600],
))

if args.selection.count('dilep'):

    plots.append(Plot(
        name = "Z1_pt",
        texX = 'p_{T}(Z_{1}) (GeV)', texY = 'Number of Events / 20 GeV',
        attribute = TreeVariable.fromString( "Z1_pt/F" ),
        binning=[20,0,400],
        addOverFlowBin='upper',
    ))

    plots.append(Plot(
        name = 'Z1_pt_coarse', texX = 'p_{T}(Z_{1}) (GeV)', texY = 'Number of Events / 50 GeV',
        attribute = TreeVariable.fromString( "Z1_pt/F" ),
        binning=[16,0,800],
        addOverFlowBin='upper',
    ))

    plots.append(Plot(
        texX = '#Delta#phi(Z_{1}(ll))', texY = 'Number of Events',
        attribute = TreeVariable.fromString( "Z1_lldPhi/F" ),
        binning=[10,0,pi],
    ))

    plots.append(Plot(
        texX = '#Delta R(Z_{1}(ll))', texY = 'Number of Events',
        attribute = TreeVariable.fromString( "Z1_lldR/F" ),
        binning=[10,0,6],
    ))

    plots.append(Plot(name="photon_Z1_deltaR",
        texX = '#Delta R(Z_{1}, #gamma)', texY = 'Number of Events',
        attribute = lambda event, sample: event.photon_Z1_deltaR,
        binning=[10,0,6],
    ))

    plots.append(Plot(name="photon_Z1_deltaPhi",
        texX = '#Delta #phi(Z_{1}, #gamma)', texY = 'Number of Events',
        attribute = lambda event, sample: event.photon_Z1_deltaPhi,
        binning=[10,0,pi],
    ))

    plots.append(Plot(name="jet0_Z1_deltaR",
        texX = '#Delta R(Z_{1}, j_{0})', texY = 'Number of Events',
        attribute = lambda event, sample: event.jet0_Z1_deltaR,
        binning=[10,0,6],
    ))

    plots.append(Plot(name="jet1_Z1_deltaR",
        texX = '#Delta R(Z_{1}, j_{1})', texY = 'Number of Events',
        attribute = lambda event, sample: event.jet1_Z1_deltaR,
        binning=[10,0,6],
    ))

for plot in plots:
    plot.name = 'FI_binned_'+plot.name

# Text on the plots
def drawObjects( hasData = False ):
    tex = ROOT.TLatex()
    tex.SetNDC()
    tex.SetTextSize(0.04)
    tex.SetTextAlign(11) # align right
    lines = [
      (0.15, 0.95, 'CMS Preliminary' if hasData else "CMS Simulation"),
      (0.45, 0.95, 'L=%3.1f fb^{-1} (13 TeV)' % lumi_scale),
    ]
    return [tex.DrawLatex(*l) for l in lines]

# draw function for plots
def drawPlots(plots):
  for log in [False, True]:
    plot_directory_ = os.path.join(plot_directory, args.plot_directory, args.sample, args.selection, ("log" if log else "lin") )
    for plot in plots:
      if not max(l[0].GetMaximum() for l in plot.histos): continue # Empty plot

      plotting.draw(plot,
        plot_directory = plot_directory_,
        ratio = None, #{'histos':[(i,0) for i in range(1,len(plot.histos))], 'yRange':(0.1,1.9)} ,
        logX = False, logY = log, sorting = False,
        yRange = (0.03, "auto") if log else "auto",
        scaling = {},
        legend =  ( (0.17,0.9-0.05*sum(map(len, plot.histos))/2,1.,0.9), 2),
        drawObjects = drawObjects( ),
        copyIndexPHP = True,
      )

plotting.fill(plots, read_variables = read_variables, sequence = sequence, max_events = -1 if args.small else -1)

for plot in plots:
    for i_h, h in enumerate(plot.histos[0]):
        # dress up
        h.legendText = params[i_h]['legendText']
        h.style = params[i_h]['style']

drawPlots(plots)

syncer.sync()

logger.info( "Done with prefix %s and selectionString %s", args.selection, cutInterpreter.cutString(args.selection) )
