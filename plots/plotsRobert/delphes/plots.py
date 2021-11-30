#!/usr/bin/env python
''' Analysis script for standard plots
'''

# Standard imports and batch mode
import ROOT, os, itertools
#ROOT.gROOT.SetBatch(True)
import copy
import random
from math                           import sqrt, cos, sin, pi, isnan, sinh, cosh, log

# Analysis
import Analysis.Tools.syncer        as syncer
from   Analysis.Tools.WeightInfo    import WeightInfo
from   Analysis.Tools.helpers       import deltaPhi, deltaR, getObjDict

# RootTools
from RootTools.core.standard        import *

# TMB
from TMB.Tools.user                 import plot_directory
from TMB.Tools.helpers              import deltaPhi, getCollection, deltaR, mZ
from TMB.Tools.delphesCutInterpreter import cutInterpreter
import TMB.Tools.VV_angles          as VV_angles
from TMB.Tools.genObjectSelection   import isGoodGenJet

# Arguments
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',           action='store',      default='INFO',          nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'], help="Log level for logging")
argParser.add_argument('--plot_directory',     action='store',      default='delphes')
argParser.add_argument('--selection',          action='store',      default='singlelep-WHJet-onH')
argParser.add_argument('--signal',             action='store',      default='WH', choices = ['WH', 'ZH'])
argParser.add_argument('--small',                                   action='store_true',     help='Run only on a small subset of the data?')
args = argParser.parse_args()

# Logger'singlelep-WHJet' if sample.name=='WH' else 'dilep-ZHJet-onZ'
import TMB.Tools.logger as logger
import RootTools.core.logger as logger_rt
logger    = logger.get_logger(   args.logLevel, logFile = None)
logger_rt = logger_rt.get_logger(args.logLevel, logFile = None)

plot_directory = os.path.join(plot_directory, args.plot_directory,  args.signal )
if args.small: plot_directory += "_small"

# Import samples
import TMB.Samples.pp_gen_v10 as samples
    
signal = getattr( samples, args.signal)
 
# WeightInfo
signal.weightInfo = WeightInfo(signal.reweight_pkl)
signal.weightInfo.set_order(2)
signal.read_variables = [VectorTreeVariable.fromString( "p[C/F]", nMax=200 )]

eft_configs = [
    (ROOT.kBlack, {}, "SM"), 
    (ROOT.kGreen+1, {'cHW':1}, "c_{W}=1"), 
    (ROOT.kGreen+3, {'cHWtil':1}, "c_{#tilde{W}}=1"), 
    (ROOT.kOrange-1, {'cHj3':1}, "c_{Hq3}=1"),
    (ROOT.kOrange-2, {'cHj3':2}, "c_{Hq3}=2"),
    ]

def get_eft_reweight( eft, weightInfo_):
    func1 = signal.weightInfo.get_weight_func(**eft)
    func2 = signal.weightInfo.get_weight_func()
    
    def eft_reweight_( event, sample, func1=func1, func2=func2):
        return func1(event, sample)/func2(event, sample)
    return eft_reweight_

if args.signal == "WH":
    stack       = Stack( [samples.TTJets, samples.WJetsToLNu_HT] )
elif args.signal == "ZH":
    stack       = Stack( [samples.DYBBJets]) 

eft_weights = [[]]
for sample in stack.samples:
    sample.style = styles.fillStyle(sample.color)
    eft_weights[0].append( None )

for _, eft, _ in eft_configs:
    stack.append( [signal] )
    eft_weights.append( [get_eft_reweight(eft, signal.weightInfo)] )


lumi  = 59.7
lumi_weight = lambda event, sample: lumi*event.lumiweight1fb#*sin(2*event.VV_angles['Theta'])*sin(2*event.VV_angles['theta_V1'])

for sample in stack.samples:
    sample.weight = lumi_weight

# Read variables and sequences
jetVars          = ['pt/F', 'eta/F', 'phi/F', 'bTag/F', 'bTagPhys/I']

jetVarNames      = [x.split('/')[0] for x in jetVars]

lepVars          = ['pt/F','eta/F','phi/F','pdgId/I','isolationVar/F', 'isolationVarRhoCorr/F']
lepVarNames      = [x.split('/')[0] for x in lepVars]

# Training variables
read_variables = [\
#    "nBTag/I", 
    "nBTag_loose/I",
    "recoMet_pt/F", "recoMet_phi/F",
    "genMet_pt/F", "genMet_phi/F",
    "recoZ_pt/F", "recoZ_eta/F", "recoZ_phi/F", "recoZ_mass/F", "recoZ_cosThetaStar/F", "recoZ_lldPhi/F", "recoZ_lldR/F", "recoZ_l1_index/I", "recoZ_l2_index/I",
    "nrecoJet/I",
    "recoJet[%s]"%(",".join(jetVars)),
    "nrecoLep/I",
    "recoLep[%s]"%(",".join(lepVars)),
    "ngenLep/I", "genLep[pt/F,eta/F,phi/F,pdgId/I,mother_pdgId/I]", 
    "lumiweight1fb/F",
    "genW[pt/F,eta/F,phi/F,l1_index/I,l2_index/I]", "ngenW/I",
    "evt/l", "run/I", "lumi/I",
    "H_dijet_mass/F", "H_pt/F", "H_j1_index/I", "H_j2_index/I", 
    "H_j1_index/I", "H_j2_index/I", 
    "WH_W_pt/F", "WH_dPhiMetLep/F", "WH_MT/F", "WH_nu_pt/F", "WH_nu_eta/F", "WH_nu_phi/F", "WH_nu_E/F", "WH_Theta/F", "WH_theta/F", "WH_phi/F", 
    "ZH_Theta/F", "ZH_theta/F", "ZH_phi/F",
]


preselection = [ 
    #("debug", "(evt==25857178)") 
    #("debug", "(genW_pt[0]>0)") 
]

selectionString  = "&&".join( [ c[1] for c in preselection] + ([cutInterpreter.cutString(args.selection)] if args.selection is not None else []))
subDirectory     =  '-'.join( [ c[0] for c in preselection] + ([args.selection] if args.selection is not None else []))
if subDirectory  == '': 
    subDirectory = 'inc'

for sample in stack.samples:
    if selectionString != "":
        sample.addSelectionString( selectionString )
    if args.small:
        #sample.reduceFiles( factor = 30 )
        sample.reduceFiles( to = 15 )

#sequence functions
sequence = []

#BITs
import sys, os, time
sys.path.insert(0,os.path.expandvars("$CMSSW_BASE/src/BIT"))
from BoostedInformationTree import BoostedInformationTree
if signal.name == 'WH':
    import TMB.BIT.configs.WH_delphes as config
    bits        = config.load("/mnt/hephy/cms/robert.schoefbeck/BIT/models/WH_delphes/default/")
elif signal.name == 'ZH':
    import TMB.BIT.configs.ZH_delphes as config
    bits        = config.load("/mnt/hephy/cms/robert.schoefbeck/BIT/models/ZH_delphes/default/")

models = [
    ("BIT_cHW",             bits[('cHW',)],             [20,-5,5]), 
    ("BIT_cHW_cHW",         bits[('cHW','cHW')],        [20,-5,15]), 
    ("BIT_cHWtil",          bits[('cHWtil',)],          [20,-5,5]), 
    ("BIT_cHWtil_cHWtil",   bits[('cHWtil','cHWtil')],  [20,-5,15]), 

#        ("BIT_cHW_cHW_coarse", bits[('cHW','cHW')],              [50,-5,5]), 
#        ("BIT_cHWtil_cHWtil_coarse",  bits[('cHWtil','cHWtil')], [50,-5,5]),
]
sequence.extend( config.sequence )

def bit_predict( event, sample ):

    for var, func in config.mva_variables:
        setattr( event, var, func(event, sample) )
    
    # get model inputs assuming lstm
    event.features = config.predict_inputs( event, sample)
    for name, model, _ in models:
        #print has_lstm, flat_variables, lstm_jets
        prediction = model.predict( event.features )
        setattr( event, name, prediction )
#        if not prediction>-float('inf'):
#            print name, prediction, [[getattr( event, mva_variable) for mva_variable, _ in config.mva_variables]]
#            print "mva_m3", event.mva_m3, "m3", event.m3, "event.nJetGood", event.nJetGood
#            raise RuntimeError("Found NAN prediction?")

sequence.append( bit_predict )

### Helpers
def addTransverseVector( p_dict ):
    ''' add a transverse vector for further calculations
    '''
    p_dict['vec2D'] = ROOT.TVector2( p_dict['pt']*cos(p_dict['phi']), p_dict['pt']*sin(p_dict['phi']) )

def addTLorentzVector( p_dict ):
    ''' add a TLorentz 4D Vector for further calculations
    '''
    p_dict['vecP4'] = ROOT.TLorentzVector( p_dict['pt']*cos(p_dict['phi']), p_dict['pt']*sin(p_dict['phi']),  p_dict['pt']*sinh(p_dict['eta']), p_dict['pt']*cosh(p_dict['eta']) )

# Use some defaults
Plot.setDefaults(stack = stack, weight = eft_weights, addOverFlowBin="upper")
  
plots        = []

postfix = "" 

for model_name, _, binning in models:
    plots.append(Plot(
        name = model_name+postfix,
        texX = model_name, texY = 'Number of Events / 10 GeV',
        attribute = lambda event, sample, model_name=model_name: getattr(event, model_name),
        #binning=Binning.fromThresholds([0, 0.5, 1, 2,3,4,10]),
        binning   = binning,
        addOverFlowBin = 'upper',
    ))

for i_key, (key, _) in enumerate( config.mva_variables ):
    plots.append(Plot( name = key.replace("mva_", "")+postfix,
      texX = config.plot_options[key]['tex'], texY = 'Number of Events',
      attribute = lambda event, sample, i_key=i_key: event.features[i_key],
      binning   =  config.plot_options[key]['binning'],
    ))


#plots.append(Plot( name = "j0_pt"+postfix,
#  texX = 'p_{T}(j_{0}) (GeV)', texY = 'Number of Events',
#  attribute = lambda event, sample: event.jets[0]['pt'] if len(event.jets)>0 else -float('inf'),
#  binning=[600/20,0,600],
#))
#
#plots.append(Plot( name = "j1_pt"+postfix,
#  texX = 'p_{T}(j_{1}) (GeV)', texY = 'Number of Events',
#  attribute = lambda event, sample: event.jets[1]['pt'] if len(event.jets)>1 else -float('inf'),
#  binning=[600/20,0,600],
#))
#
#plots.append(Plot( name = "j2_pt"+postfix,
#  texX = 'p_{T}(j_{2}) (GeV)', texY = 'Number of Events',
#  attribute = lambda event, sample: event.jets[2]['pt'] if len(event.jets)>2 else -float('inf'),
#  binning=[600/20,0,600],
#))
#
#plots.append(Plot( name = "j0_eta"+postfix,
#  texX = '#eta(j_{0}) (GeV)', texY = 'Number of Events',
#  attribute = lambda event, sample: event.jets[0]['eta'] if len(event.jets)>0 else -float('inf'),
#  binning=[30,-3,3],
#))
#
#plots.append(Plot( name = "j1_eta"+postfix,
#  texX = '#eta(j_{1}) (GeV)', texY = 'Number of Events',
#  attribute = lambda event, sample: event.jets[1]['eta'] if len(event.jets)>1 else -float('inf'),
#  binning=[30,-3,3],
#))
#
#plots.append(Plot( name = "j2_eta"+postfix,
#  texX = '#eta(j_{2}) (GeV)', texY = 'Number of Events',
#  attribute = lambda event, sample: event.jets[2]['eta'] if len(event.jets)>2 else -float('inf'),
#  binning=[30,-3,3],
#))
#
#plots.append(Plot( name = "b0_pt"+postfix,
#  texX = 'p_{T}(b_{0}) (GeV)', texY = 'Number of Events',
#  attribute = lambda event, sample: event.bJets[0]['pt'] if len(event.bJets)>0 else -float('inf'),
#  binning=[600/20,0,600],
#))
#
#plots.append(Plot( name = "b1_pt"+postfix,
#  texX = 'p_{T}(b_{1}) (GeV)', texY = 'Number of Events',
#  attribute = lambda event, sample: event.bJets[1]['pt'] if len(event.bJets)>1 else -float('inf'),
#  binning=[600/20,0,600],
#))
#
#plots.append(Plot( name = "b01PtRatio"+postfix,
#  texX = 'p_{T}(b_{1})/p_{T}(b_{0})', texY = 'Number of Events',
#  attribute = lambda event, sample: event.bJets[1]['pt']/event.bJets[0]['pt'] if len(event.bJets)>1 else -float('inf'),
#  binning=[20,0,1],
#))
#
#plots.append(Plot( name = "deltaPhib01"+postfix,
#  texX = '#Delta#Phi(b_{0},b_{1})', texY = 'Number of Events',
#  attribute = lambda event, sample: deltaPhi(event.bJets[0]['phi'], event.bJets[1]['phi']) if len(event.bJets)>1 else -float('inf'),
#  binning=[20,0,pi],
#))
#
#plots.append(Plot( name = "b0_eta"+postfix,
#  texX = '#eta(b_{0}) (GeV)', texY = 'Number of Events',
#  attribute = lambda event, sample: event.bJets[0]['eta'] if len(event.bJets)>0 else -float('inf'),
#  binning=[30,-3,3],
#))
#
#plots.append(Plot( name = "b1_eta"+postfix,
#  texX = '#eta(b_{1}) (GeV)', texY = 'Number of Events',
#  attribute = lambda event, sample: event.bJets[1]['eta'] if len(event.bJets)>1 else -float('inf'),
#  binning=[30,-3,3],
#))
#
#plots.append(Plot( name = 'Met_pt'+postfix,
#  texX = 'E_{T}^{miss} (GeV)', texY = 'Number of Events / 20 GeV',
#  attribute = lambda event, sample: event.recoMet_pt,
#  binning=[400/20,0,400],
#))
#
#plots.append(Plot( name = 'thrust'+postfix,
#  texX = 'Thrust (GeV)', texY = 'Number of Events / 20 GeV',
#  attribute = lambda event, sample: event.thrust,
#  binning=[20,0,1.2],
#))
#
#plots.append(Plot( name = 'thrust_min'+postfix,
#  texX = 'min(Thrust) (GeV)', texY = 'Number of Events / 20 GeV',
#  attribute = lambda event, sample: event.thrust_min,
#  binning=[20,0,1.2],
#))
#
#plots.append(Plot( name = 'nJet'+postfix,
#  texX = 'jet multiplicity', texY = 'Number of Events / 20 GeV',
#  attribute = lambda event, sample: len(event.jets),
#  binning=[8,0,8],
#))
#
#plots.append(Plot( name = 'H_pt'+postfix,
#  texX = ' p_{T}(H) (GeV)', texY = 'Number of Events / 20 GeV',
#  attribute = lambda event, sample: event.H_pt,
#  binning=[600/20,0,600],
#))
#
#plots.append(Plot( name = 'H_mass_noCut'+postfix,
#  texX = ' M(b_{1},b_{2}) (GeV)', texY = 'Number of Events / 20 GeV',
#  attribute = lambda event, sample: event.H_vecP4.M(),
#  binning=[600/20,0,600],
#))
#
#plots.append(Plot( name = 'H_mass'+postfix,
#  texX = ' M(b_{1},b_{2}) (GeV)', texY = 'Number of Events / 20 GeV',
#  attribute = lambda event, sample: event.H_vecP4.M(),
#  binning=[600/20,0,600],
#))
#
#plots.append(Plot( name = "deltaPhiVH"+postfix,
#  texX = '#Delta#Phi(V,H)', texY = 'Number of Events',
#  attribute = lambda event, sample: deltaPhi(event.V_vecP4.Phi(), event.H_vecP4.Phi()),
#  binning=[20,0,pi],
#))
#
#plots.append(Plot( name = "Theta"+postfix,
#  texX = '#Theta', texY = 'Number of Events',
#  attribute = lambda event, sample: event.Theta,
#  binning=[20,0,pi],
#))
#
#plots.append(Plot( name = "theta"+postfix,
#  texX = '#theta', texY = 'Number of Events',
#  attribute = lambda event, sample: event.theta,
#  binning=[20,0,pi],
#))
#
#plots.append(Plot( name = 'V_pt'+postfix,
#  texX = 'p_{T}(V) (GeV)', texY = 'Number of Events / 20 GeV',
#  attribute = lambda event, sample: event.V_pt,
#  binning=[600/20,0,600],
#))
#
#plots.append(Plot( name = 'V_pt_coarse'+postfix,
#  texX = 'p_{T}(V) (GeV)', texY = 'Number of Events',
#  attribute = lambda event, sample: event.V_pt,
#  binning=Binning.fromThresholds([150, 200, 250, 300, 360, 430, 510, 590, 690, 800, 900, 1100]),
#))
#
#plots.append(Plot( name = 'V_eta'+postfix,
#  texX = '#eta(V) ', texY = 'Number of Events / 20 GeV',
#  attribute = lambda event, sample: event.V_vecP4.Eta(),
#  binning=[20,-3,3],
#))
#
#
#if args.signal == 'ZH':
#
#    plots.append(Plot( name = 'Z_mass'+postfix,
#      texX = 'm(l_{1},l_{2}) (GeV)', texY = 'Number of Events / 20 GeV',
#      attribute = lambda event, sample: event.recoZ_mass if event.selection else -float('inf'),
#      binning=[20,70,110],
#    ))
#
#    plots.append(Plot( name = "Z_cosThetaStar"+postfix,
#      texX = 'cos(#theta^{*})', texY = 'Number of Events',
#      attribute = lambda event, sample: event.recoZ_cosThetaStar if event.selection else -float('inf'),
#      binning=[20,-1,1],
#    ))
#
#    plots.append(Plot( name = "Z_lldPhi"+postfix,
#      texX = '#Delta#phi(l_{1},l_{2})', texY = 'Number of Events',
#      attribute = lambda event, sample: event.recoZ_lldPhi if event.selection else -float('inf'),
#      binning=[20,0,pi],
#    ))
#
#    plots.append(Plot( name = "Z_lldR"+postfix,
#      texX = '#Delta R(l_{1},l_{2})', texY = 'Number of Events',
#      attribute = lambda event, sample: event.recoZ_lldR if event.selection else -float('inf'),
#      binning=[20,0,7],
#    ))
#
#    plots.append(Plot( name = "lep1_pt"+postfix,
#      texX = 'leading p_{T}(l) (GeV)', texY = 'Number of Events',
#      attribute = lambda event, sample: event.lepton1['pt'] if event.selection else -float('inf'),
#      binning=[300/20,0,300],
#    ))
#
#    plots.append(Plot( name = "lep1_eta"+postfix,
#      texX = 'leading l. #eta(l)', texY = 'Number of Events',
#      attribute = lambda event, sample: event.lepton1['eta'] if event.selection else -float('inf'),
#      binning=[20,-3,3],
#    ))
#
#    plots.append(Plot( name = "lep2_pt"+postfix,
#      texX = 'leading p_{T}(l) (GeV)', texY = 'Number of Events',
#      attribute = lambda event, sample: event.lepton2['pt'] if event.selection else -float('inf'),
#      binning=[300/20,0,300],
#    ))
#
#    plots.append(Plot( name = "lep2_eta"+postfix,
#      texX = 'leading l. #eta(l)', texY = 'Number of Events',
#      attribute = lambda event, sample: event.lepton2['eta'] if event.selection else -float('inf'),
#      binning=[20,-3,3],
#    ))
#
#if args.signal == 'WH':
#
#    plots.append(Plot( name = 'genW_pt'+postfix,
#      texX = 'gen p_{T}(W) (GeV)', texY = 'Number of Events / 20 GeV',
#      attribute = lambda event, sample: event.genW_pt[0] if event.ngenW>0 and event.selection else -float('inf'),
#      binning=[600/20,0,600],
#    ))
#
#    plots.append(Plot( name = 'genW_pt_noHighPtV'+postfix,
#      texX = 'gen p_{T}(W) (GeV)', texY = 'Number of Events / 20 GeV',
#      attribute = lambda event, sample: event.genW_pt[0] if event.ngenW>0 and event.selection_noHighPtV else -float('inf'),
#      binning=[600/20,0,600],
#    ))
#
#    plots.append(Plot( name = "lep_pt"+postfix,
#      texX = 'p_{T}(l) (GeV)', texY = 'Number of Events',
#      attribute = lambda event, sample: event.lepton['pt'] if event.selection else -float('inf'),
#      binning=[300/20,0,300],
#    ))
#
#    plots.append(Plot( name = "lep_eta"+postfix,
#      texX = '#eta(l)', texY = 'Number of Events',
#      attribute = lambda event, sample: event.lepton['eta'] if event.selection else -float('inf'),
#      binning=[20,-3,3],
#    ))
#
#    plots.append(Plot( name = 'MT'+postfix,
#      texX = 'M_{T} (GeV)', texY = 'Number of Events / 20 GeV',
#      attribute = lambda event, sample: event.MT if event.selection else -float('inf'),
#      binning=[400/20,0,400],
#    ))
#
#    plots.append(Plot( name = 'deltaPhiMetLep_noCut'+postfix,
#      texX = '#Delta#Phi(E_{T}^{miss},l)', texY = 'Number of Events / 20 GeV',
#      attribute = lambda event, sample: event.dPhiMetLep if event.selection_noPhiMetLep else -float('inf'),
#      binning=[20,0,pi],
#    ))
#
#    plots.append(Plot( name = 'deltaPhiMetLep'+postfix,
#      texX = '#Delta#Phi(E_{T}^{miss},l)', texY = 'Number of Events / 20 GeV',
#      attribute = lambda event, sample: event.dPhiMetLep if event.selection else -float('inf'),
#      binning=[20,0,pi],
#    ))
#
#    plots.append(Plot( name = 'W_mass'+postfix,
#      texX = ' m(l,#nu) (GeV)', texY = 'Number of Events / 20 GeV',
#      attribute = lambda event, sample: event.V_vecP4.M() if event.selection else -float('inf'),
#      binning=[600/20,0,600],
#    ))


# Text on the plots
def drawObjects( hasData = False ):
    tex = ROOT.TLatex()
    tex.SetNDC()
    tex.SetTextSize(0.04)
    tex.SetTextAlign(11) # align right
    lines = [
      (0.15, 0.95, 'CMS Preliminary' if hasData else "CMS Simulation"), 
      #(0.45, 0.95, 'L=%3.1f fb{}^{-1} (13 TeV) Scale %3.2f'% ( lumi_scale, dataMCScale ) ) if plotData else (0.45, 0.95, 'L=%3.1f fb{}^{-1} (13 TeV)' % lumi_scale)
    ]
    return [tex.DrawLatex(*l) for l in lines] 

# draw function for plots
def drawPlots(plots, subDirectory=''):
  for log in [False, True]:
    plot_directory_ = os.path.join(plot_directory, subDirectory)
    plot_directory_ = os.path.join(plot_directory_, "log") if log else os.path.join(plot_directory_, "lin")
    for plot in plots:
      if not max(l[0].GetMaximum() for l in plot.histos): continue # Empty plot

      plotting.draw(plot,
	    plot_directory = plot_directory_,
	    ratio =  None,
	    logX = False, logY = log, sorting = False,
	    yRange = (0.03, "auto") if log else "auto",
	    scaling = {},
	    legend =  ( (0.17,0.9-0.05*sum(map(len, plot.histos))/2,1.,0.9), 2),
	    drawObjects = drawObjects( ),
        copyIndexPHP = True,
      )

plotting.fill(plots, read_variables = read_variables, sequence = sequence, max_events = -1 if args.small else -1)

#plot_phi_subtr = copy.deepcopy(plots[subtr_phi_ind])
#plot_phi_subtr.name = plot_phi_subtr.name.replace("pos_","subtr_")
#plot_phi_subtr.texX = plot_phi_subtr.texX.replace("(pos.)","(subtr)")
#for i_h_, h_ in enumerate(plots[subtr_phi_ind].histos):
#    for i_h, h in enumerate(h_):
#        h_sub = plots[subtr_phi_ind+1].histos[i_h_][i_h].Clone()
#        h_sub.Scale(-1)
#        plot_phi_subtr.histos[i_h_][i_h].Add(h_sub)
#
#plots.append( plot_phi_subtr )

#color EFT
for plot in plots:
    for i, (color, _, texName) in enumerate(eft_configs):
        plot.histos[i+1][0].legendText = texName
        plot.histos[i+1][0].style = styles.lineStyle(color,width=2)

#plot_phi_subtr.histos = plot_phi_subtr.histos[1:]

drawPlots(plots, subDirectory = subDirectory)

logger.info( "Done with prefix %s and selectionString %s", args.selection, cutInterpreter.cutString(args.selection) )

syncer.sync()

