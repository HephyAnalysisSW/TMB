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
argParser.add_argument('--plot_directory',     action='store',      default='delphes/WH')
argParser.add_argument('--selection',          action='store',      default='singlelep-WHJet')
argParser.add_argument('--WC',                 action='store',      default='cHj3')
argParser.add_argument('--WCval',              action='store',      nargs = '*',             type=float,    default=[1.0],  help='Values of the Wilson coefficient')
argParser.add_argument('--small',                                   action='store_true',     help='Run only on a small subset of the data?')
args = argParser.parse_args()

# Logger
import TMB.Tools.logger as logger
import RootTools.core.logger as logger_rt
logger    = logger.get_logger(   args.logLevel, logFile = None)
logger_rt = logger_rt.get_logger(args.logLevel, logFile = None)

if args.small: args.plot_directory += "_small"

# Import samples
import TMB.Samples.pp_gen_v10 as samples
#from TMB.Samples.pp_gen_v4 import *

lumi  = 137
stack = Stack( [samples.TTJets, samples.WJetsToLNu, samples.WH] )
weight= lambda event, sample: lumi*event.lumiweight1fb

# WeightInfo
for sample in stack.samples:
    sample.style = styles.fillStyle(sample.color)

    if hasattr( sample, "reweight_pkl" ):
        sample.w = WeightInfo(sample.reweight_pkl)
        sample.w.set_order(2)
        sample.read_variables = [VectorTreeVariable.fromString( "p[C/F]", nMax=200 )]

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
    "recoZ_pt/F", "recoZ_eta/F", "recoZ_phi/F", "recoZ_mass/F", "recoZ_cosThetaStar/F", "recoZ_lldPhi/F", "recoZ_lldR/F",
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
    p_dict['vec4D'] = ROOT.TLorentzVector( p_dict['pt']*cos(p_dict['phi']), p_dict['pt']*sin(p_dict['phi']),  p_dict['pt']*sinh(p_dict['eta']), p_dict['pt']*cosh(p_dict['eta']) )

#sequence functions
sequence = []

from TMB.Tools.objectSelection import isBJet
def make_jets( event, sample ):
    event.jets     = [getObjDict(event, 'recoJet_', jetVarNames, i) for i in range(int(event.nrecoJet))]
    for p in event.jets:
        #addTransverseVector( p )
        addTLorentzVector( p )
    event.bJets    = filter(lambda j:j['bTag']>=1 and abs(j['eta'])<=2.4    , event.jets)

sequence.append( make_jets )

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
    if len(event.leps)>0:
        event.lepton = event.leps[0]
        random_no      = gRandom.Uniform(0,1)
        event.neutrino_vec4D = VV_angles.neutrino_mom(event.lepton['vec4D'], event.recoMet_pt, event.recoMet_phi, random_no)
        event.W_vec4D  = event.neutrino_vec4D + event.lepton['vec4D']  

        event.W_pt         = event.W_vec4D.Pt()
        event.has_highPt_W = event.W_pt > 150
    else:
        event.lepton = None 
        event.neutrino_vec4D = None
        event.W_vec4D = None
        event.W_pt = -1
        event.has_highPt_W = False

sequence.append( makeLeptonic )

def makeH( event, sample ):
    event.dijet_mass = (event.bJets[0]['vec4D'] + event.bJets[1]['vec4D']).M() 
    event.H_vec4D = event.bJets[0]['vec4D'] + event.bJets[1]['vec4D']
    event.H_pt    = event.H_vec4D.Pt()
    if event.dijet_mass>90 and event.dijet_mass<150:
        event.has_H =  1
    else:
        event.has_H = -1 

sequence.append( makeH )

def makeSelection( event, sample):
    event.selection     = event.has_highPt_W and event.has_H 
    event.selection_noH = event.has_highPt_W 

sequence.append( makeSelection )

# Use some defaults
Plot.setDefaults(stack = stack, weight = staticmethod(weight), addOverFlowBin="upper")
  
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
  attribute = lambda event, sample: event.recoMet_pt,
  binning=[400/20,0,400],
))

plots.append(Plot( name = 'nJet'+postfix,
  texX = 'jet multiplicity', texY = 'Number of Events / 20 GeV',
  attribute = lambda event, sample: len(event.jets),
  binning=[8,0,8],
))

plots.append(Plot( name = 'W_pt'+postfix,
  texX = ' p_{T}(W) (GeV)', texY = 'Number of Events / 20 GeV',
  attribute = lambda event, sample: event.W_pt if event.selection else -float('inf'),
  binning=[600/20,0,600],
))

plots.append(Plot( name = 'H_pt'+postfix,
  texX = ' p_{T}(H) (GeV)', texY = 'Number of Events / 20 GeV',
  attribute = lambda event, sample: event.H_pt if event.selection else -float('inf'),
  binning=[600/20,0,600],
))

plots.append(Plot( name = 'W_mass'+postfix,
  texX = ' m(l,#nu) (GeV)', texY = 'Number of Events / 20 GeV',
  attribute = lambda event, sample: event.W_vec4D.M() if event.selection else -float('inf'),
  binning=[600/20,0,600],
))

plots.append(Plot( name = 'H_mass_noCut'+postfix,
  texX = ' M(b_{1},b_{2}) (GeV)', texY = 'Number of Events / 20 GeV',
  attribute = lambda event, sample: event.H_vec4D.M() if event.selection_noH else -float('inf'),
  binning=[600/20,0,600],
))

plots.append(Plot( name = 'H_mass'+postfix,
  texX = ' M(b_{1},b_{2}) (GeV)', texY = 'Number of Events / 20 GeV',
  attribute = lambda event, sample: event.H_vec4D.M() if event.selection else -float('inf'),
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
    plot_directory_ = os.path.join(plot_directory, args.plot_directory, subDirectory)
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

drawPlots(plots, subDirectory = subDirectory)

logger.info( "Done with prefix %s and selectionString %s", args.selection, cutInterpreter.cutString(args.selection) )
