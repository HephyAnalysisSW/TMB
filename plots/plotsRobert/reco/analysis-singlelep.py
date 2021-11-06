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
from Analysis.Tools.WeightInfo           import WeightInfo
import Analysis.Tools.syncer             as     syncer
import numpy as np

# Arguments
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',       action='store',      default='INFO', nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'], help="Log level for logging")
argParser.add_argument('--small',                             action='store_true', help='Run only on a small subset of the data?', )
#argParser.add_argument('--sorting',                           action='store', default=None, choices=[None, "forDYMB"],  help='Sort histos?', )
argParser.add_argument('--plot_directory', action='store', default='analysis-v6')
argParser.add_argument('--era',            action='store', type=str, default="Autumn18")
argParser.add_argument('--sample',        action='store', type=str, default="ttG_noFullyHad_fast")
argParser.add_argument('--WC',            action='store', type=str, default="cWWW")
argParser.add_argument('--selection',      action='store', default='singlelep-photon')
argParser.add_argument('--onlyMVA',       action='store', default=None, help='Plot only this MVA')
args = argParser.parse_args()

# Logger
import TMB.Tools.logger as logger
import RootTools.core.logger as logger_rt
logger    = logger.get_logger(   args.logLevel, logFile = None)
logger_rt = logger_rt.get_logger(args.logLevel, logFile = None)

if args.small: args.plot_directory += "_small"

import TMB.Samples.pp_tWZ_v6 as samples

lumi_scale  = 137
lumi_weight = lambda event, sample: lumi_scale*event.weight

if args.WC == 'cWWW':
    mc = [ samples.WGToLNu_fast, samples.ttG_noFullyHad_fast ]
else:
    mc = [ samples.ttG_noFullyHad_fast, samples.WGToLNu_fast ]

for sample in mc:
    sample.style  = styles.fillStyle( sample.color )
    sample.weightInfo = WeightInfo( sample.reweight_pkl )
    sample.weightInfo.set_order(2)

    total_norm = len(sample.files) # ...well
    if args.small:
        sample.reduceFiles( to = 10 )
    sample.scale = total_norm/float(len(sample.files))

    sample.weight = sample.weightInfo.get_weight_func( ) 
        

BSM = copy.deepcopy(mc)
for sample in BSM:
    sample.style = styles.invisibleStyle()
    sample.notInLegend = True
    if args.WC in  sample.weightInfo.variables:
        sample.weight = sample.weightInfo.get_weight_func( **{args.WC:1} )
        logger.info( "Set %s=1 for sample %s", args.WC, sample.name )
    else:
        sample.weight = sample.weightInfo.get_weight_func( )
        logger.info( "DO NOT set %s=1 for sample %s", args.WC, sample.name )

BSM[0].notInLegend = False
BSM[0].style = styles.lineStyle(ROOT.kBlack, dashed = True)
BSM[0].texName = '%s=1'%args.WC

stack = Stack( mc, BSM)

# Read variables and sequences
sequence       = []

read_variables = [
    "weight/F", "year/I", "met_pt/F", "met_phi/F", "nBTag/I", "nJetGood/I", "PV_npvsGood/I",
    "l1_pt/F", "l1_eta/F" , "l1_phi/F", "l1_pdgId/I", "l1_mvaTOP/F", "l1_mvaTOPWP/I", "l1_index/I", 
    #"l2_pt/F", "l2_eta/F" , "l2_phi/F", "l2_mvaTOP/F", "l2_mvaTOPWP/I", "l2_index/I",
    #"l3_pt/F", "l3_eta/F" , "l3_phi/F", "l3_mvaTOP/F", "l3_mvaTOPWP/I", "l3_index/I",
    "JetGood[pt/F,eta/F,phi/F]",
    "lep[pt/F,eta/F,phi/F,pdgId/I,muIndex/I,eleIndex/I]",
#    "Z1_l1_index/I", "Z1_l2_index/I", "nonZ1_l1_index/I", "nonZ1_l2_index/I", 
#    "Z1_phi/F", "Z1_pt/F", "Z1_mass/F", "Z1_cosThetaStar/F", "Z1_eta/F", "Z1_lldPhi/F", "Z1_lldR/F",
    "Muon[pt/F,eta/F,phi/F,dxy/F,dz/F,ip3d/F,sip3d/F,jetRelIso/F,miniPFRelIso_all/F,pfRelIso03_all/F,mvaTOP/F,mvaTTH/F,pdgId/I,segmentComp/F,nStations/I,nTrackerLayers/I]",
    "Electron[pt/F,eta/F,phi/F,dxy/F,dz/F,ip3d/F,sip3d/F,jetRelIso/F,miniPFRelIso_all/F,pfRelIso03_all/F,mvaTOP/F,mvaTTH/F,pdgId/I,vidNestedWPBitmap/I]",
    "np/I", VectorTreeVariable.fromString("p[C/F]",nMax=500),
    "photon_pt/F",
    "photon_eta/F",
    "photon_phi/F",
    "photonJetdR/F", "photonLepdR/F", "m3/F",
]

# define 3l selections

import TMB.Tools.VV_angles          as VV_angles
import random
def make_VV( event, sample ):
    ##AN2019_059_v8 p22
    #mW      = 80.4
    #mt2     = 2*event.l1_pt*event.met_pt*(1-cos(event.met_phi-event.l1_phi))
    #delta   = sqrt((mW**2-mt2)/(2.*event.l1_pt*event.met_pt)) 
    #if mW**2<mt2: 
    #    random_sign  = 2*(-.5+(random.random()>0.5))
    #    eta_neu = l1_eta + random_sign*log(1.+delta*sqrt(2+delta**2)+delta**2)
    #else:
    #    eta_neu = l1_eta

    # A. Wulzer lep decay angles in 2007.10356:  The latter angles are in the rest frame of each boson and they are 
    # defined as those of the final fermion or anti-fermion with helicity +1/2 (e.g. the l+ in the case
    # of a W+ and the nu-bar for a W-), denoted as f_+ in the Fig. 6

    lep_4 = ROOT.TLorentzVector()
    lep_4.SetPtEtaPhiM(event.l1_pt, event.l1_eta, event.l1_phi, 0)

    random_number = ROOT.gRandom.Uniform() 
    neu_4         = VV_angles.neutrino_mom( lep_4, event.met_pt, event.met_phi, random_number ) 

    # the helicity+ fermion is the l+ (from a W+), otherwise it's the neutrino
    lep_m, lep_p  = (neu_4, lep_4) if event.l1_pdgId<0 else (lep_4, neu_4)

    gamma_4 = ROOT.TLorentzVector()
    gamma_4.SetPtEtaPhiM(event.photon_pt, event.photon_eta, event.photon_phi, 0)

    event.thetaW = VV_angles.gettheta(lep_m, lep_p, gamma_4)
    event.Theta  = VV_angles.getTheta(lep_m, lep_p, gamma_4)
    event.phiW   = VV_angles.getphi(lep_m, lep_p, gamma_4)

    #print "MT", sqrt(2*event.l1_pt*event.met_pt*(1-cos(event.l1_phi-event.met_phi)))
    #lep_4.Print()
    #neu_4.Print()
    #gamma_4.Print()

    #print event.thetaW, event.Theta, event.phiW
    #print

sequence.append( make_VV )

#MVA
import TMB.MVA.configs as configs
config = configs.ttG_WG
read_variables += config.read_variables

# Add sequence that computes the MVA inputs
def make_mva_inputs( event, sample ):
    for mva_variable, func in config.mva_variables:
        setattr( event, mva_variable, func(event, sample) )
sequence.append( make_mva_inputs )

## load models
#from keras.models import load_model
#
#if args.onlyMVA is not None:
#    has_lstm = ('LSTM' in args.onlyMVA)
#    name = args.onlyMVA.split('/')[-4]
#    models = [ (name, has_lstm, load_model(args.onlyMVA) ), ]
#else:
#    if args.WC == 'ctZ':
#        models = [
#    #        ("FI_ctZ_BSM_TTG",      False, load_model("/mnt/hephy/cms/robert.schoefbeck/TMB/models/ctZ_BSM_TTG/ttG_WG/FI_ctZ_BSM/regression_model.h5")),
#    #        ("FI_ctZ_BSM_TTG_wq",   False, load_model("/mnt/hephy/cms/robert.schoefbeck/TMB/models/ctZ_BSM_TTG_wq/ttG_WG/FI_ctZ_BSM/regression_model.h5")),
#            ("FI_ctZ_BSM_TTGWG_wq", False, load_model("/mnt/hephy/cms/robert.schoefbeck/TMB/models/ctZ_BSM_TTGWG_wq/ttG_WG/FI_ctZ_BSM/regression_model.h5")),
#            ("FI_ctZ_BSM_TTGWG",    False, load_model("/mnt/hephy/cms/robert.schoefbeck/TMB/models/ctZ_BSM_TTGWG/ttG_WG/FI_ctZ_BSM/regression_model.h5")),
#            ("FI_ctZ_BSM_TTGWG_LSTM_wq",  True, load_model("/mnt/hephy/cms/robert.schoefbeck/TMB/models/ctZ_TTGWG_wq_LSTM/ttG_WG/FI_ctZ_BSM/regression_model.h5")),
#        ]
#    elif args.WC == 'cWWW':
#        models = [
#            ("FI_cWWW_regression_ttGWG", False, load_model("/mnt/hephy/cms/robert.schoefbeck/TMB/models/train_ctZ_BSM_in_ttG_noFullyHad_WGToLNu/ttG_WG/FI_cWWW_BSM/regression_model.h5")),
#            ("FI_cWWW_regression_WG",    False, load_model("/mnt/hephy/cms/robert.schoefbeck/TMB/models/train_ctZ_BSM_in_WGToLNu/ttG_WG/FI_cWWW_BSM/regression_model.h5")),
#        ]
#
#def keras_predict( event, sample ):
#
#    # get model inputs assuming lstm
#    flat_variables, lstm_jets = config.predict_inputs( event, sample, jet_lstm = True)
#    for name, has_lstm, model in models:
#        #print has_lstm, flat_variables, lstm_jets
#        prediction = model.predict( flat_variables if not has_lstm else [flat_variables, lstm_jets] )
#        setattr( event, name, prediction )
#        if not prediction>-float('inf'):
#            print name, prediction, [[getattr( event, mva_variable) for mva_variable, _ in config.mva_variables]]
#            print "mva_m3", event.mva_m3, "m3", event.m3, "event.nJetGood", event.nJetGood
#            raise RuntimeError("Found NAN prediction?")
#
#sequence.append( keras_predict )

#BITs
import sys, os, time
sys.path.insert(0,os.path.expandvars("$CMSSW_BASE/src/BIT"))
from BoostedInformationTree import BoostedInformationTree 
import TMB.BIT.configs.ttG_WG as config

bits        = config.load("/mnt/hephy/cms/robert.schoefbeck/BIT/models/ttG_WG/clipScore/")
#bit_ctZ_ctZ_clipped   = BoostedInformationTree.load('/mnt/hephy/cms/robert.schoefbeck/BIT/models/ttG_WG/cpSfix/bit_derivative_ctZ_ctZ.pkl')
#bit_cWWW_cWWW_clipped = BoostedInformationTree.load('/mnt/hephy/cms/robert.schoefbeck/BIT/models/ttG_WG/cpSfix/bit_derivative_cWWW_cWWW.pkl')
models = [
    ("BIT_ctZ", bits[('ctZ',)], [50, -.4, .6,]),
    ("BIT_ctZ_ctZ", bits[('ctZ','ctZ')], [50, -1, 9,]),
#    ("BIT_ctZ_ctZ_clipped", bit_ctZ_ctZ_clipped, [50, -1, 9,]),
    ("BIT_cWWW",    bits[('cWWW',)], [50, -.1, .1,]),
    ("BIT_cWWW_cWWW",  bits[('cWWW','cWWW')], [50, -.02, .04,]),
#    ("BIT_cWWW_cWWW_clipped", bit_cWWW_cWWW_clipped, [50, -0.02, 0.04,]),
]

def bit_predict( event, sample ):

    # get model inputs assuming lstm
    features = config.predict_inputs( event, sample)
    for name, model, _ in models:
        #print has_lstm, flat_variables, lstm_jets
        prediction = model.predict( features )
        setattr( event, name, prediction )
#        if not prediction>-float('inf'):
#            print name, prediction, [[getattr( event, mva_variable) for mva_variable, _ in config.mva_variables]]
#            print "mva_m3", event.mva_m3, "m3", event.m3, "event.nJetGood", event.nJetGood
#            raise RuntimeError("Found NAN prediction?")

sequence.append( bit_predict )


mu_string  = lepString('mu','VL')
ele_string = lepString('ele','VL')
def getLeptonSelection():
    return "Sum$({mu_string})+Sum$({ele_string})==1".format(mu_string=mu_string,ele_string=ele_string)

yields     = {}
allPlots   = {}

#yt_TWZ_filter.scale = lumi_scale * 1.07314

selectionString = cutInterpreter.cutString(args.selection) 
logger.info("selectionString: %s", selectionString)

# Use some defaults
Plot.setDefaults(stack = stack, weight = staticmethod(lumi_weight), selectionString = selectionString)

plots        = []
fisher_plots = []

for model_name, _, binning in models:
    plots.append(Plot(
        name = model_name,
        texX = model_name, texY = 'Number of Events / 10 GeV',
        attribute = lambda event, sample, model_name=model_name: getattr(event, model_name),
        binning=binning,
        addOverFlowBin='upper',
    ))

if args.onlyMVA is None:

    plots.append(Plot(
      name = 'nVtxs', texX = 'vertex multiplicity', texY = 'Number of Events',
      attribute = TreeVariable.fromString( "PV_npvsGood/I" ),
      binning=[50,0,50],
      addOverFlowBin='upper',
    ))

    plots.append(Plot(
        name = 'M3',
        texX = 'M_{3} (GeV)', texY = 'Number of Events / 10 GeV',
        attribute = lambda event, sample: event.m3 if event.nJetGood >=3 else 0,
        binning=[30,0,600],
        addOverFlowBin='both',
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
        texX = 'p_{T}(#gamma) (GeV)', texY = 'Number of Events / 40 GeV',
        attribute = lambda event, sample:event.photon_pt,
        binning=[25,0,500],
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

    plots.append(Plot(
        name = 'phiW',
        texX = '#phi(W)', texY = 'Number of Events',
        attribute = lambda event, sample: event.phiW,
        binning=[30,-pi,pi],
        addOverFlowBin='upper',
    ))
    plots.append(Plot(
        name = 'thetaW',
        texX = '#theta(W)', texY = 'Number of Events',
        attribute = lambda event, sample: event.thetaW,
        binning=[30,0,pi],
        addOverFlowBin='both',
    ))
    plots.append(Plot(
        name = 'Theta',
        texX = '#Theta', texY = 'Number of Events',
        attribute = lambda event, sample: event.Theta,
        binning=[30,0,pi],
        addOverFlowBin='both',
    ))

# Text on the plots
def drawObjects( dev = None, hasData = False ):
    tex = ROOT.TLatex()
    tex.SetNDC()
    tex.SetTextSize(0.04)
    tex.SetTextAlign(11) # align right
    lines = [
      (0.15, 0.95, 'CMS Preliminary' if hasData else "CMS Simulation"),
      (0.45, 0.95, 'L=%3.1f fb{}^{-1} (13 TeV)' % lumi_scale),
    ]
    if dev is not None:
        lines.append( (0.65, 0.8, 'dev=%3.1f'%dev) )
    return [tex.DrawLatex(*l) for l in lines]

# draw function for plots
def drawPlots(plots):
  for log in [False, True]:
    plot_directory_ = os.path.join(plot_directory, args.plot_directory, args.WC, args.selection, ("log" if log else "lin") )
    for plot in plots:
      if not max(l[0].GetMaximum() for l in plot.histos): continue # Empty plot

      plotting.draw(plot,
        plot_directory = plot_directory_,
        ratio = {'histos':[(1,0)], 'yRange':(0.8,1.2), 'texY':'BSM/SM'},
        logX = False, logY = log, sorting = False,
        yRange = (0.03, "auto") if log else "auto",
        scaling = {},
        legend =  ( (0.17,0.9-0.05*sum(map(len, plot.histos))/2,1.,0.9), 2),
        drawObjects = drawObjects( dev = plot.dev if hasattr(plot, "dev") else None),
        copyIndexPHP = True,
      )

plotting.fill(plots, read_variables = read_variables, sequence = sequence, max_events = -1 if args.small else -1)

for p in plots:
    h = p.histos_added
    p.dev = 0
    for i_bin in range(1,h[0][0].GetNbinsX()+1):
        if h[0][0].GetBinContent(i_bin)>0:
            p.dev += abs(h[1][0].GetBinContent(i_bin) - h[0][0].GetBinContent(i_bin))/sqrt(h[0][0].GetBinContent(i_bin))

#for plot in plots:
#    for i_h, hl in enumerate(plot.histos):
#        # dress up
#        hl[0].legendText = params[i_h]['legendText']
#        hl[0].style = params[i_h]['style']
#for p in plots:
#    if p.histos[0][0].Integral()<p.histos[0][1].Integral():
#        p.histos[0][0], p.histos[0][1] = p.histos[0][1], p.histos[0][0]

drawPlots(plots)

syncer.sync()

logger.info( "Done with prefix %s and selectionString %s", args.selection, cutInterpreter.cutString(args.selection) )
