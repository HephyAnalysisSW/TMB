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
argParser.add_argument('--selection',          action='store',      default='singlelep-WHJet')
argParser.add_argument('--signal',             action='store',      default='WH', choices = ['WH', 'ZH'])
argParser.add_argument('--small',                                   action='store_true',     help='Run only on a small subset of the data?')
args = argParser.parse_args()

# Logger
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

lumi  = 137
lumi_weight = lambda event, sample: lumi*event.lumiweight1fb

eft_configs = [
    (ROOT.kBlack, {}, "SM"), 
    (ROOT.kGreen+1, {'cHW':1}, "c_{W}=1"), 
    (ROOT.kOrange-1, {'cHj3':1}, "c_{Hq3}=1"),
    (ROOT.kOrange-2, {'cHj3':.7}, "c_{Hq3}=.7"),
    ]

def get_eft_reweight( eft, weightInfo_):
    func1 = signal.weightInfo.get_weight_func(**eft)
    func2 = signal.weightInfo.get_weight_func()
    
    def eft_reweight_( event, sample, func1=func1, func2=func2):
        return func1(event, sample)/func2(event, sample)
    return eft_reweight_

if args.signal == "WH":
    stack       = Stack( [samples.TTJets, samples.WJetsToLNu] )
elif args.signal == "ZH":
    stack       = Stack( [samples.DYJets, samples.WJetsToLNu] )

eft_weights = [[]]
for sample in stack.samples:
    sample.style = styles.fillStyle(sample.color)
    eft_weights[0].append( None )

for _, eft, _ in eft_configs:
    stack.append( [signal] )
    eft_weights.append( [get_eft_reweight(eft, signal.weightInfo)] )

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
    "recoZ_pt/F", "recoZ_eta/F", "recoZ_phi/F", "recoZ_mass/F", "recoZ_cosThetaStar/F", "recoZ_lldPhi/F", "recoZ_lldR/F", "recoZ_l1_index/I", "recoZ_l2_index/I",
    "nrecoJet/I",
    "recoJet[%s]"%(",".join(jetVars)),
    "nrecoLep/I",
    "recoLep[%s]"%(",".join(lepVars)),
    "lumiweight1fb/F",
]

preselection = [ 
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
        sample.reduceFiles( to = 2 )

### Helpers
def addTransverseVector( p_dict ):
    ''' add a transverse vector for further calculations
    '''
    p_dict['vec2D'] = ROOT.TVector2( p_dict['pt']*cos(p_dict['phi']), p_dict['pt']*sin(p_dict['phi']) )

def addTLorentzVector( p_dict ):
    ''' add a TLorentz 4D Vector for further calculations
    '''
    p_dict['vecP4'] = ROOT.TLorentzVector( p_dict['pt']*cos(p_dict['phi']), p_dict['pt']*sin(p_dict['phi']),  p_dict['pt']*sinh(p_dict['eta']), p_dict['pt']*cosh(p_dict['eta']) )

#sequence functions
sequence = []

from TMB.Tools.objectSelection import isBJet
def makeJets( event, sample ):
    event.jets     = [getObjDict(event, 'recoJet_', jetVarNames, i) for i in range(int(event.nrecoJet))]
    for p in event.jets:
        #addTransverseVector( p )
        addTLorentzVector( p )
    event.bJets    = filter(lambda j:j['bTag']>=1 and abs(j['eta'])<=2.4    , event.jets)

sequence.append( makeJets )

gRandom = ROOT.TRandom3() #TRandom3(time.gmtime(0)) # for changing seed
def makeLeptonic( event, sample ):

    # Read leptons, do not yet filter
    event.all_leps = getCollection( event, 'recoLep', ['pt', 'eta', 'phi', 'pdgId', 'isolationVar', 'isolationVarRhoCorr'], 'nrecoLep' )
    # Add extra vectors
    for p in event.all_leps:
        #addTransverseVector( p )
        addTLorentzVector( p )

    # Sort
    event.leps = sorted( event.all_leps, key=lambda k: -k['pt'] )

    # Cross-cleaning: remove leptons that overlap with a jet within 0.4
    event.leps = list(filter( lambda l: min( [ deltaR(l, j) for j in event.jets ] + [999] ) > 0.4 , event.leps ))
    if args.signal == "WH":
        if len(event.leps)>0:
            event.lepton = event.leps[0]
            random_no      = gRandom.Uniform(0,1)
            event.neutrino_vecP4 = VV_angles.neutrino_mom(event.lepton['vecP4'], event.recoMet_pt, event.recoMet_phi, random_no)
            event.V_vecP4  = event.neutrino_vecP4 + event.lepton['vecP4']  

            event.V_pt          = event.V_vecP4.Pt()
            event.has_highPt_V  = event.V_pt > 150
            event.dPhiMetLep    = abs(deltaPhi( event.recoMet_phi, event.lepton['phi'] ))
            event.MT            = sqrt(2.*event.recoMet_pt*event.lepton['pt']*(1-cos(event.dPhiMetLep)))
        else:
            event.lepton            = None 
            event.neutrino_vecP4    = None
            event.V_vecP4           = None
            event.V_pt              = float('nan')
            event.has_highPt_V      = False
            event.dPhiMetLep        = float('nan')
            event.MT                = float('nan')

    elif args.signal == "ZH":

        #"recoZ_pt/F", "recoZ_eta/F", "recoZ_phi/F", "recoZ_mass/F", "recoZ_cosThetaStar/F", "recoZ_lldPhi/F", "recoZ_lldR/F", "recoZ_l1_index/I", "recoZ_l2_index/I",
        if event.recoZ_pt>0:
            event.V_pt          =   event.recoZ_pt
            event.V_vecP4       =   event.all_leps[event.recoZ_l1_index]['vecP4'] + event.all_leps[event.recoZ_l2_index]['vecP4'] 
            event.lepton1       =   event.all_leps[event.recoZ_l1_index]
            event.lepton2       =   event.all_leps[event.recoZ_l1_index]
            event.has_highPt_V  =   event.recoZ_pt > 75
        else:
            event.V_pt          =   float('nan')
            event.V_vecP4       =   False 
            event.lepton1       =   False 
            event.lepton2       =   False 
            event.has_highPt_V  =   False 

sequence.append( makeLeptonic )

def makeH( event, sample ):
    event.dijet_mass = (event.bJets[0]['vecP4'] + event.bJets[1]['vecP4']).M() 
    event.H_vecP4 = event.bJets[0]['vecP4'] + event.bJets[1]['vecP4']
    event.H_pt    = event.H_vecP4.Pt()
    if event.dijet_mass>90 and event.dijet_mass<150:
        event.has_H = 1
    else:
        event.has_H = 0  

sequence.append( makeH )

# thrust
from TMB.Tools.Thrust import Thrust
def makeThrust( event, sample ):
    if event.V_vecP4 is not None:
        (event.thrust, event.thrust_min) = Thrust(event.V_vecP4,[j['vecP4'] for j in event.jets] )
    else:
        (event.thrust, event.thrust_min) = -1, -1 
        
sequence.append( makeThrust )

# selection bools
if args.signal == 'WH':
    def makeSelection( event, sample):
        event.selection     = event.has_highPt_V and event.has_H and event.dPhiMetLep < 2
        event.selection_noH = event.has_highPt_V and event.dPhiMetLep < 2 
        event.selection_noPhiMetLep = event.has_highPt_V and event.has_H
elif args.signal == 'ZH':
    def makeSelection( event, sample):
        event.selection     = event.has_highPt_V and event.has_H 
        event.selection_noH = event.has_highPt_V 

sequence.append( makeSelection )

# Use some defaults
Plot.setDefaults(stack = stack, weight = eft_weights, addOverFlowBin="upper")
  
plots        = []
fisher_plots = []

postfix = "" 

plots.append(Plot( name = "j0_pt"+postfix,
  texX = 'p_{T}(j_{0}) (GeV)', texY = 'Number of Events',
  attribute = lambda event, sample: event.jets[0]['pt'] if event.selection and len(event.jets)>0 else -float('inf'),
  binning=[600/20,0,600],
))

plots.append(Plot( name = "j1_pt"+postfix,
  texX = 'p_{T}(j_{1}) (GeV)', texY = 'Number of Events',
  attribute = lambda event, sample: event.jets[1]['pt'] if event.selection and len(event.jets)>1 else -float('inf'),
  binning=[600/20,0,600],
))

plots.append(Plot( name = "j2_pt"+postfix,
  texX = 'p_{T}(j_{2}) (GeV)', texY = 'Number of Events',
  attribute = lambda event, sample: event.jets[2]['pt'] if event.selection and len(event.jets)>2 else -float('inf'),
  binning=[600/20,0,600],
))

plots.append(Plot( name = "j0_eta"+postfix,
  texX = '#eta(j_{0}) (GeV)', texY = 'Number of Events',
  attribute = lambda event, sample: event.jets[0]['eta'] if event.selection and len(event.jets)>0 else -float('inf'),
  binning=[30,-3,3],
))

plots.append(Plot( name = "j1_eta"+postfix,
  texX = '#eta(j_{1}) (GeV)', texY = 'Number of Events',
  attribute = lambda event, sample: event.jets[1]['eta'] if event.selection and len(event.jets)>1 else -float('inf'),
  binning=[30,-3,3],
))

plots.append(Plot( name = "j2_eta"+postfix,
  texX = '#eta(j_{2}) (GeV)', texY = 'Number of Events',
  attribute = lambda event, sample: event.jets[2]['eta'] if event.selection and len(event.jets)>2 else -float('inf'),
  binning=[30,-3,3],
))

plots.append(Plot( name = "b0_pt"+postfix,
  texX = 'p_{T}(b_{0}) (GeV)', texY = 'Number of Events',
  attribute = lambda event, sample: event.bJets[0]['pt'] if event.selection and len(event.bJets)>0 else -float('inf'),
  binning=[600/20,0,600],
))

plots.append(Plot( name = "b1_pt"+postfix,
  texX = 'p_{T}(b_{1}) (GeV)', texY = 'Number of Events',
  attribute = lambda event, sample: event.bJets[1]['pt'] if event.selection and len(event.bJets)>1 else -float('inf'),
  binning=[600/20,0,600],
))

plots.append(Plot( name = "b01PtRatio"+postfix,
  texX = 'p_{T}(b_{1})/p_{T}(b_{0})', texY = 'Number of Events',
  attribute = lambda event, sample: event.bJets[1]['pt']/event.bJets[0]['pt'] if event.selection and len(event.bJets)>1 else -float('inf'),
  binning=[20,0,1],
))

plots.append(Plot( name = "deltaPhib01"+postfix,
  texX = '#Delta#Phi(b_{0},b_{1})', texY = 'Number of Events',
  attribute = lambda event, sample: deltaPhi(event.bJets[0]['phi'], event.bJets[1]['phi']) if event.selection and len(event.bJets)>1 else -float('inf'),
  binning=[20,0,pi],
))

plots.append(Plot( name = "b0_eta"+postfix,
  texX = '#eta(b_{0}) (GeV)', texY = 'Number of Events',
  attribute = lambda event, sample: event.bJets[0]['eta'] if event.selection and len(event.bJets)>0 else -float('inf'),
  binning=[30,-3,3],
))

plots.append(Plot( name = "b1_eta"+postfix,
  texX = '#eta(b_{1}) (GeV)', texY = 'Number of Events',
  attribute = lambda event, sample: event.bJets[1]['eta'] if event.selection and len(event.bJets)>1 else -float('inf'),
  binning=[30,-3,3],
))

plots.append(Plot( name = 'Met_pt'+postfix,
  texX = 'E_{T}^{miss} (GeV)', texY = 'Number of Events / 20 GeV',
  attribute = lambda event, sample: event.recoMet_pt if event.selection else -float('inf'),
  binning=[400/20,0,400],
))

plots.append(Plot( name = 'thrust'+postfix,
  texX = 'Thrust (GeV)', texY = 'Number of Events / 20 GeV',
  attribute = lambda event, sample: event.thrust if event.selection else -float('inf'),
  binning=[20,0,1.2],
))

plots.append(Plot( name = 'thrust_min'+postfix,
  texX = 'min(Thrust) (GeV)', texY = 'Number of Events / 20 GeV',
  attribute = lambda event, sample: event.thrust_min if event.selection else -float('inf'),
  binning=[20,0,1.2],
))

plots.append(Plot( name = 'nJet'+postfix,
  texX = 'jet multiplicity', texY = 'Number of Events / 20 GeV',
  attribute = lambda event, sample: len(event.jets) if event.selection else -float('inf'),
  binning=[8,0,8],
))

plots.append(Plot( name = 'H_pt'+postfix,
  texX = ' p_{T}(H) (GeV)', texY = 'Number of Events / 20 GeV',
  attribute = lambda event, sample: event.H_pt if event.selection else -float('inf'),
  binning=[600/20,0,600],
))

plots.append(Plot( name = 'H_mass_noCut'+postfix,
  texX = ' M(b_{1},b_{2}) (GeV)', texY = 'Number of Events / 20 GeV',
  attribute = lambda event, sample: event.H_vecP4.M() if event.selection_noH else -float('inf'),
  binning=[600/20,0,600],
))

plots.append(Plot( name = 'H_mass'+postfix,
  texX = ' M(b_{1},b_{2}) (GeV)', texY = 'Number of Events / 20 GeV',
  attribute = lambda event, sample: event.H_vecP4.M() if event.selection else -float('inf'),
  binning=[600/20,0,600],
))

plots.append(Plot( name = "deltaPhiVH"+postfix,
  texX = '#Delta#Phi(V,H)', texY = 'Number of Events',
  attribute = lambda event, sample: deltaPhi(event.V_vecP4.Phi(), event.H_vecP4.Phi()) if event.selection else -float('inf'),
  binning=[20,0,pi],
))

plots.append(Plot( name = 'V_pt'+postfix,
  texX = 'p_{T}(V) (GeV)', texY = 'Number of Events / 20 GeV',
  attribute = lambda event, sample: event.V_pt if event.selection else -float('inf'),
  binning=[600/20,0,600],
))

plots.append(Plot( name = 'V_pt'+postfix,
  texX = '#eta(V) ', texY = 'Number of Events / 20 GeV',
  attribute = lambda event, sample: event.V_vecP4.Eta() if event.selection else -float('inf'),
  binning=[20,-3,3],
))


if args.signal == 'ZH':

    plots.append(Plot( name = 'Z_mass'+postfix,
      texX = 'm(l_{1},l_{2}) (GeV)', texY = 'Number of Events / 20 GeV',
      attribute = lambda event, sample: event.recoZ_mass if event.selection else -float('inf'),
      binning=[20,70,110],
    ))

    plots.append(Plot( name = "Z_cosThetaStar"+postfix,
      texX = 'cos(#theta^{*})', texY = 'Number of Events',
      attribute = lambda event, sample: event.recoZ_cosThetaStar if event.selection else -float('inf'),
      binning=[20,-1,1],
    ))

    plots.append(Plot( name = "Z_lldPhi"+postfix,
      texX = '#Delta#phi(l_{1},l_{2})', texY = 'Number of Events',
      attribute = lambda event, sample: event.recoZ_lldPhi if event.selection else -float('inf'),
      binning=[20,0,pi],
    ))

    plots.append(Plot( name = "Z_lldR"+postfix,
      texX = '#Delta R(l_{1},l_{2})', texY = 'Number of Events',
      attribute = lambda event, sample: event.recoZ_lldR if event.selection else -float('inf'),
      binning=[20,0,7],
    ))

    plots.append(Plot( name = "lep1_pt"+postfix,
      texX = 'leading p_{T}(l) (GeV)', texY = 'Number of Events',
      attribute = lambda event, sample: event.lepton1['pt'] if event.selection else -float('inf'),
      binning=[300/20,0,300],
    ))

    plots.append(Plot( name = "lep1_eta"+postfix,
      texX = 'leading l. #eta(l)', texY = 'Number of Events',
      attribute = lambda event, sample: event.lepton1['eta'] if event.selection else -float('inf'),
      binning=[20,-3,3],
    ))

    plots.append(Plot( name = "lep2_pt"+postfix,
      texX = 'leading p_{T}(l) (GeV)', texY = 'Number of Events',
      attribute = lambda event, sample: event.lepton2['pt'] if event.selection else -float('inf'),
      binning=[300/20,0,300],
    ))

    plots.append(Plot( name = "lep2_eta"+postfix,
      texX = 'leading l. #eta(l)', texY = 'Number of Events',
      attribute = lambda event, sample: event.lepton2['eta'] if event.selection else -float('inf'),
      binning=[20,-3,3],
    ))

if args.signal == 'WH':
    plots.append(Plot( name = "lep_pt"+postfix,
      texX = 'p_{T}(l) (GeV)', texY = 'Number of Events',
      attribute = lambda event, sample: event.lepton['pt'] if event.selection else -float('inf'),
      binning=[300/20,0,300],
    ))

    plots.append(Plot( name = "lep_eta"+postfix,
      texX = '#eta(l)', texY = 'Number of Events',
      attribute = lambda event, sample: event.lepton['eta'] if event.selection else -float('inf'),
      binning=[20,-3,3],
    ))

    plots.append(Plot( name = 'MT'+postfix,
      texX = 'M_{T} (GeV)', texY = 'Number of Events / 20 GeV',
      attribute = lambda event, sample: event.MT if event.selection else -float('inf'),
      binning=[400/20,0,400],
    ))

    plots.append(Plot( name = 'deltaPhiMetLep_noCut'+postfix,
      texX = '#Delta#Phi(E_{T}^{miss},l)', texY = 'Number of Events / 20 GeV',
      attribute = lambda event, sample: event.dPhiMetLep if event.selection_noPhiMetLep else -float('inf'),
      binning=[20,0,pi],
    ))

    plots.append(Plot( name = 'deltaPhiMetLep'+postfix,
      texX = '#Delta#Phi(E_{T}^{miss},l)', texY = 'Number of Events / 20 GeV',
      attribute = lambda event, sample: event.dPhiMetLep if event.selection else -float('inf'),
      binning=[20,0,pi],
    ))

    plots.append(Plot( name = 'W_mass'+postfix,
      texX = ' m(l,#nu) (GeV)', texY = 'Number of Events / 20 GeV',
      attribute = lambda event, sample: event.V_vecP4.M() if event.selection else -float('inf'),
      binning=[600/20,0,600],
    ))


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

      len_FI = len(plot.fisher_plots) if hasattr(plot, "fisher_plots") else 0
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

plotting.fill(plots+fisher_plots, read_variables = read_variables, sequence = sequence, max_events = -1 if args.small else -1)

#color EFT
for plot in plots:
    for i, (color, _, texName) in enumerate(eft_configs):

        plot.histos[i+1][0].legendText = texName
        plot.histos[i+1][0].style = styles.lineStyle(color,width=2)

drawPlots(plots, subDirectory = subDirectory)

logger.info( "Done with prefix %s and selectionString %s", args.selection, cutInterpreter.cutString(args.selection) )

syncer.sync()

