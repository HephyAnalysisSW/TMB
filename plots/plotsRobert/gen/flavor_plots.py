#!/usr/bin/env python
''' Analysis script for standard plots
'''

# Standard imports and batch mode
import ROOT, os, itertools
#ROOT.gROOT.SetBatch(True)
import copy
from math                           import sqrt, cos, sin, pi, isnan, sinh, cosh, log

# Analysis
import Analysis.Tools.syncer        as syncer
from   Analysis.Tools.WeightInfo    import WeightInfo

# RootTools
from RootTools.core.standard        import *

# TMB
from TMB.Tools.user                 import plot_directory
from TMB.Tools.helpers              import deltaPhi, getCollection, deltaR, mZ
from TMB.Tools.genCutInterpreter    import cutInterpreter
import TMB.Tools.VV_angles          as VV_angles
from TMB.Tools.genObjectSelection   import isGoodGenJet

# Arguments
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',           action='store',      default='INFO',          nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'], help="Log level for logging")
argParser.add_argument('--plot_directory',     action='store',      default='gen_flavor')
argParser.add_argument('--selection',          action='store',      default=None)
argParser.add_argument('--sample',             action='store',      default='ttG_noFullyHad')
argParser.add_argument('--small',                                   action='store_true',     help='Run only on a small subset of the data?')
argParser.add_argument('--ZZ',                                      action='store_true',     help='Add diboson ZZ angles?')
argParser.add_argument('--WG',                                      action='store_true',     help='Add diboson WG angles?')
args = argParser.parse_args()

# Logger
import TMB.Tools.logger as logger
import RootTools.core.logger as logger_rt
logger    = logger.get_logger(   args.logLevel, logFile = None)
logger_rt = logger_rt.get_logger(args.logLevel, logFile = None)

if args.small: args.plot_directory += "_small"

# Import samples
from TMB.Samples.pp_gen_v7 import *
sample = eval(args.sample) 
# objects to plot
objects = sample.objects if hasattr( sample, "objects") else []

# WeightInfo
w = WeightInfo(sample.reweight_pkl)
w.set_order(2)

# define which Wilson coefficients to plot

WCs = [
    ('cWWW', 1.0, ROOT.kRed),
    ('cpq3i', 1.0,ROOT.kGreen),
    ('cpqMi', 1.0, ROOT.kOrange),
    ('cpu', 1.0, ROOT.kCyan),
    ('cpt', 1.0, ROOT.kMagenta),
    ('cpQ3', 1.0, ROOT.kBlue),
    ('cpQM', 1.0, ROOT.kViolet-9),
    ('cpb', 1.0, ROOT.kGray),
    ('cpd', 1.0, ROOT.kAzure+10),
]

params =  [ ] 

for i_wc, (WC, WCval, color) in enumerate(WCs):
    params.append ({'legendText':'%s=%3.2f'%(WC, WCval), 'color':color,  'WC':{WC:WCval} })
params +=  [ {'legendText':'SM',  'color':ROOT.kBlack, 'WC':{}} ] 

for i_param, param in enumerate(params):
    param['sample'] = sample
    param['style']  = styles.lineStyle( param['color'] )

stack = Stack(*[ [ param['sample'] ] for param in params ] )
weight= [ [ w.get_weight_func(**param['WC']) ] for param in params ]
    

# Read variables and sequences
read_variables = [
    "genMet_pt/F", "genMet_phi/F", 
    "ngenJet/I", "genJet[pt/F,eta/F,phi/F,matchBParton/I]", 
    "ngenLep/I", "genLep[pt/F,eta/F,phi/F,pdgId/I,mother_pdgId/I]", 
    "ngenTop/I", "genTop[pt/F,eta/F,phi/F]", 
    "ngenZ/I", "genZ[pt/F,phi/F,eta/F,daughter_pdgId/I,l1_index/I,l2_index/I,mother_pdgId/I,grandmother_pdgId/I]",
    "ngenW/I", "genW[pt/F,phi/F,eta/F,daughter_pdgId/I,l1_index/I,l2_index/I,mother_pdgId/I,grandmother_pdgId/I]",
    "ngenPhoton/I", "genPhoton[pt/F,phi/F,eta/F,mother_pdgId/I,grandmother_pdgId/I]",
]
read_variables.append( VectorTreeVariable.fromString('p[C/F]', nMax=2000) )

preselection = [ 
    #( "SMP-20-005",  "genPhoton_pt[0]>300&&genMet_pt>80&&genLep_pt[genW_l1_index[0]]>80&&sqrt(acos(cos(genLep_phi[genW_l1_index[0]]-genPhoton_phi[0]))**2+(genLep_eta[genW_l1_index[0]]-genPhoton_eta[0])**2)>3.0"),
    #( "SMP-20-005-light",  "genPhoton_pt[0]>150&genMet_pt>30&&genLep_pt[genW_l1_index[0]]>30&&sqrt(acos(cos(genLep_phi[genW_l1_index[0]]-genPhoton_phi[0]))**2+(genLep_eta[genW_l1_index[0]]-genPhoton_eta[0])**2)>3.0"),
]

selectionString  = "&&".join( [ c[1] for c in preselection] + ([cutInterpreter.cutString(args.selection)] if args.selection is not None else []))
subDirectory     =  '-'.join( [ c[0] for c in preselection] + ([args.selection] if args.selection is not None else []))
if subDirectory  == '': 
    subDirectory = 'inc'

for sample in stack.samples:
    if selectionString != "":
        sample.addSelectionString( selectionString )
    if args.small:
        sample.reduceFiles( to = 1 )

## Helpers
def addTransverseVector( p_dict ):
    ''' add a transverse vector for further calculations
    '''
    p_dict['vec2D'] = ROOT.TVector2( p_dict['pt']*cos(p_dict['phi']), p_dict['pt']*sin(p_dict['phi']) )

def addTLorentzVector( p_dict ):
    ''' add a TLorentz 4D Vector for further calculations
    '''
    p_dict['vec4D'] = ROOT.TLorentzVector( p_dict['pt']*cos(p_dict['phi']), p_dict['pt']*sin(p_dict['phi']),  p_dict['pt']*sinh(p_dict['eta']), 0 )

#sequence functions
sequence = []

def makeJets( event, sample ):
    ''' Add a list of filtered all jets to the event
    '''

    # Retrieve & filter
    event.jets = getCollection( event, 'genJet', ['pt', 'eta', 'phi', 'matchBParton'], 'ngenJet' )
    event.jets = list(filter( lambda j:isGoodGenJet( j ), event.jets))
    # sort
    event.jets = sorted( event.jets, key=lambda k: -k['pt'] )

    # Add extra vectors
    for p in event.jets:
        addTransverseVector( p )
        addTLorentzVector( p )

    # True B's
    event.trueBjets    = list( filter( lambda j: j['matchBParton'], event.jets ) )
    event.trueNonBjets = list( filter( lambda j: not j['matchBParton'], event.jets ) )

    # Mimick b reconstruction ( if the trailing b fails acceptance, we supplement with the leading non-b jet ) 
    # event.bj0, event.bj1 = (event.trueBjets + event.trueNonBjets + [None, None])[:2]
    
sequence.append( makeJets )

def makeMET( event, sample ):
    ''' Make a MET vector to facilitate further calculations
    '''
    event.MET = {'pt':event.genMet_pt, 'phi':event.genMet_phi}
    addTransverseVector( event.MET )

sequence.append( makeMET )

def makeLeps( event, sample ):

    # Read leptons, do not yet filter
    event.all_leps = getCollection( event, 'genLep', ['pt', 'eta', 'phi', 'pdgId', 'mother_pdgId'], 'ngenLep' )
    # Add extra vectors
    for p in event.all_leps:
        addTransverseVector( p )
        addTLorentzVector( p )

    # Sort
    event.leps = sorted( event.all_leps, key=lambda k: -k['pt'] )

    # Cross-cleaning: remove leptons that overlap with a jet within 0.4
    event.leps = list(filter( lambda l: min( [ deltaR(l, j) for j in event.jets ] + [999] ) > 0.4 , event.leps ))

    # find leptons from Z
    event.lepsFromZ = list( filter( lambda j: j['mother_pdgId'] == 23 , event.leps ) )
    event.foundZ    = len( event.lepsFromZ )==2 and event.lepsFromZ[0]['pdgId'] * event.lepsFromZ[1]['pdgId'] < 0
    event.Z_deltaPhi_ll = deltaPhi( event.lepsFromZ[0]['phi'], event.lepsFromZ[1]['phi']) if event.foundZ else float('nan')
    event.Z_deltaR_ll   = deltaR( *event.lepsFromZ) if event.foundZ else float('nan')
 
    # convinience
    #event.Z_unitVec2D = UnitVectorT2( event.Z_phi )
    #event.Z_vec4D     = ROOT.TLorentzVector()
    #event.Z_vec4D.SetPtEtaPhiM( event.Z_pt, event.Z_eta, event.Z_phi, event.Z_mass )
    #event.Z_unitVec3D = event.Z_vec4D.Vect()
    #event.Z_unitVec3D /= event.Z_unitVec3D.Mag() 

    # find leptons that are NOT from Z 
    event.leptonsNotFromZ = [lepton for lepton in event.leps if lepton not in event.lepsFromZ] 

sequence.append( makeLeps )

# interference resurrection
if args.ZZ:
    def make_ZZ_angles( event, sample ):
        # nothing to resurrect here
        event.genZs = filter( lambda z: abs(z['daughter_pdgId']) in [11, 13, 15], getCollection( event, 'genZ', ['pt', 'eta', 'phi', 'daughter_pdgId', 'l1_index', 'l2_index'], 'ngenZ' ) )
        if len(event.genZs)>1:
            id1 = event.genZs[0]['l1_index']
            id2 = event.genZs[0]['l2_index']
            id3 = event.genZs[1]['l1_index']
            id4 = event.genZs[1]['l2_index']

            l = event.all_leps
            if l[id1]['pdgId']<0:
                id1, id2 = id2, id1 
            if l[id3]['pdgId']<0:
                id3, id4 = id4, id3 

            v = [ ROOT.TLorentzVector() for i in range(4) ]
            for i_id, id_ in enumerate( [id1, id2, id3, id4] ):
                #print "Lepton",i_id, l[id_]
                v[i_id].SetPtEtaPhiM(l[id_]['pt'], l[id_]['eta'], l[id_]['phi'], 0)

            event.VV_Theta = VV_angles.getTheta( v[0], v[1], v[2]+v[3] )
            event.VV_theta1= VV_angles.gettheta( v[0], v[1], v[2]+v[3] )
            event.VV_theta2= VV_angles.gettheta( v[2], v[3], v[0]+v[1] )
            event.VV_phi1  = VV_angles.getphi( v[0], v[1], v[2]+v[3] )
            event.VV_phi2  = VV_angles.getphi( v[2], v[3], v[0]+v[1] )
            event.VV_dphi  = event.VV_phi1 + event.VV_phi2 #ZZ_phis differ by a sign
        else:
            event.VV_Theta = float('nan') 
            event.VV_theta1= float('nan') 
            event.VV_theta2= float('nan') 
            event.VV_phi1  = float('nan') 
            event.VV_phi2  = float('nan') 
            event.VV_dphi  = float('nan') 
        #print event.ZZ_angles
        #print
        #assert False, "" 
    sequence.append( make_ZZ_angles )

if args.WG:
    def make_WG_angles( event, sample ):
        # SMP-20-005 (gen-level version) 
        #print  getCollection( event, 'genW', ['pt', 'eta', 'phi', 'daughter_pdgId', 'l1_index', 'l2_index'], 'ngenW' )
        event.genWs      = filter( lambda w: abs(w['daughter_pdgId']) in [11, 13, 15], getCollection( event, 'genW', ['pt', 'eta', 'phi', 'daughter_pdgId', 'l1_index', 'l2_index'], 'ngenW' ) )
        event.genPhotons = getCollection( event, 'genPhoton', ['pt', 'eta', 'phi'], 'ngenPhoton' )
        if len(event.genWs)>0 and len(event.genPhotons)>0:
            genW = event.genWs[0]
            genPhoton = event.genPhotons[0]

            l = event.all_leps

            v = [ ROOT.TLorentzVector() for i in range(2) ]
            id1, id2 = (event.genWs[0]['l1_index'], event.genWs[0]['l2_index']) if abs(event.genLep_pdgId[event.genWs[0]['l1_index']]) in [ 11, 13, 15 ] else (event.genWs[0]['l2_index'], event.genWs[0]['l1_index'])

            for i_id, id_ in enumerate( [id1, id2] ):
                v[i_id].SetPtEtaPhiM(l[id_]['pt'], l[id_]['eta'], l[id_]['phi'], 0)
            v_photon = ROOT.TLorentzVector()
            v_photon.SetPtEtaPhiM(genPhoton['pt'], genPhoton['eta'],genPhoton['phi'], 0)

            event.VV_Theta = VV_angles.getTheta( v[0], v[1], v_photon )
            event.VV_theta1= VV_angles.gettheta( v[0], v[1], v_photon )
            event.VV_phi1  = VV_angles.getphi( v[0], v[1], v_photon )
        else:
            event.VV_Theta = float('nan') 
            event.VV_theta1= float('nan') 
            event.VV_phi1  = float('nan') 

        #print event.WG_angles
        #print 
    sequence.append( make_WG_angles )

# Use some defaults
Plot.setDefaults(stack = stack, weight = weight, addOverFlowBin=None)
  
plots        = []
fisher_plots = []

postfix = "" #'_'+args.WC + ("_FIth" if args.FI_thresholds is not None else "") 

if args.ZZ:
    plots.append(Plot( name = "phi1"+postfix,
      texX = '#phi_{1}', texY = 'Number of Events',
      attribute = lambda event, sample: event.VV_phi1,
      binning=[40, -pi, pi], addOverFlowBin="both",
    ))

    plots.append(Plot( name = "phi2"+postfix,
      texX = '#phi_{2}', texY = 'Number of Events',
      attribute = lambda event, sample: event.VV_phi2,
      binning=[40, -pi, pi], addOverFlowBin="both",
    ))

    plots.append(Plot( name = "theta1"+postfix,
      texX = '#theta_{1}', texY = 'Number of Events',
      attribute = lambda event, sample: event.VV_theta1,
      binning=[40, 0, pi], addOverFlowBin="both",
    ))

    plots.append(Plot( name = "theta2"+postfix,
      texX = '#theta_{2}', texY = 'Number of Events',
      attribute = lambda event, sample: event.VV_theta2,
      binning=[40, 0, pi], addOverFlowBin="both",
    ))

    plots.append(Plot( name = "Theta"+postfix,
      texX = '#Theta', texY = 'Number of Events',
      attribute = lambda event, sample: event.VV_Theta,
      binning=[40, 0, pi], addOverFlowBin="both",
    ))

    plots.append(Plot( name = "deltaPhi"+postfix,
      texX = '#Delta #Phi', texY = 'Number of Events',
      attribute = lambda event, sample: event.VV_dphi,
      binning=[40, 0, 2*pi], addOverFlowBin="both",
    ))


if 'g' in objects:

    plots.append(Plot( name = "Photon0_pt"+postfix,
      texX = 'p_{T}(#gamma_{0}) (GeV)', texY = 'Number of Events',
      attribute = lambda event, sample: event.genPhoton_pt[0] if event.ngenPhoton>0 else float('nan'),
      binning=[600/20,0,600],
    ))

    plots.append(Plot( name = "Photon0_eta"+postfix,
      texX = '#eta(#gamma_{0}) (GeV)', texY = 'Number of Events',
      attribute = lambda event, sample: event.genPhoton_eta[0] if event.ngenPhoton>0 else float('nan'),
      binning=[30,-3,3],
    ))

    plots.append(Plot( name = "Photon0_mother_pdgId"+postfix,
      texX = 'id(#gamma_{0}) (GeV)', texY = 'Number of Events',
      attribute = lambda event, sample: event.genPhoton_mother_pdgId[0] if event.ngenPhoton>0 else float('nan'),
      binning=[50,-25,25], 
    ))

    plots.append(Plot( name = "Photon0_pt_ljm"+postfix,
      texX = 'p_{T}(#gamma_{0}) (GeV)', texY = 'Number of Events',
      attribute = lambda event, sample: event.genPhoton_pt[0] if event.ngenPhoton>0 and abs(event.genPhoton_mother_pdgId[0]) in [1,2,3,4] else float('nan'),
      binning=[600/20,0,600],
    ))

    plots.append(Plot( name = "Photon0_eta_ljm"+postfix,
      texX = '#eta(#gamma_{0}) (GeV)', texY = 'Number of Events',
      attribute = lambda event, sample: event.genPhoton_eta[0] if event.ngenPhoton>0 and abs(event.genPhoton_mother_pdgId[0]) in [1,2,3,4] else float('nan'),
      binning=[30,-3,3],
    ))

    plots.append(Plot( name = "Photon0_pt_nljm"+postfix,
      texX = 'p_{T}(#gamma_{0}) (GeV)', texY = 'Number of Events',
      attribute = lambda event, sample: event.genPhoton_pt[0] if event.ngenPhoton>0 and abs(event.genPhoton_mother_pdgId[0]) not in [1,2,3,4] else float('nan'),
      binning=[600/20,0,600],
    ))

    plots.append(Plot( name = "Photon0_eta_nljm"+postfix,
      texX = '#eta(#gamma_{0}) (GeV)', texY = 'Number of Events',
      attribute = lambda event, sample: event.genPhoton_eta[0] if event.ngenPhoton>0 and abs(event.genPhoton_mother_pdgId[0]) not in [1,2,3,4] else float('nan'),
      binning=[30,-3,3],
    ))

if args.WG:
    plots.append(Plot( name = "phi1"+postfix,
      texX = '#phi_{1}', texY = 'Number of Events',
      attribute = lambda event, sample: event.VV_phi1,
      binning=[40, -pi, pi], addOverFlowBin="both",
    ))

    plots.append(Plot( name = "theta1"+postfix,
      texX = '#theta_{1}', texY = 'Number of Events',
      attribute = lambda event, sample: event.VV_theta1,
      binning=[40, 0, pi], addOverFlowBin="both",
    ))

    plots.append(Plot( name = "Theta"+postfix,
      texX = '#Theta', texY = 'Number of Events',
      attribute = lambda event, sample: event.VV_Theta,
      binning=[40, 0, pi], addOverFlowBin="both",
    ))


if 'Z' in objects:
    plots.append(Plot( name = "Z0_pt"+postfix,
      texX = 'p_{T}(Z_{0}) (GeV)', texY = 'Number of Events',
      attribute = lambda event, sample: event.genZ_pt[0] if event.ngenZ>0 else float('nan'),
      binning=[600/20,0,600],
    ))

    plots.append(Plot( name = "Z0_eta"+postfix,
      texX = '#eta(Z_{0}) (GeV)', texY = 'Number of Events',
      attribute = lambda event, sample: event.genZ_eta[0] if event.ngenZ>0 else float('nan'),
      binning=[30,-3,3],
    ))

    plots.append(Plot( name = "Z1_pt"+postfix,
      texX = 'p_{T}(Z_{1}) (GeV)', texY = 'Number of Events',
      attribute = lambda event, sample: event.genZ_pt[1] if event.ngenZ>1 else float('nan'),
      binning=[600/20,0,600],
    ))

    plots.append(Plot( name = "Z1_eta"+postfix,
      texX = '#eta(Z_{1}) (GeV)', texY = 'Number of Events',
      attribute = lambda event, sample: event.genZ_eta[1] if event.ngenZ>1 else float('nan'),
      binning=[30,-3,3],
    ))

    plots.append(Plot( name = "Z0_mother_pdgId"+postfix,
      texX = 'id(Z_{0}) (GeV)', texY = 'Number of Events',
      attribute = lambda event, sample: event.genZ_mother_pdgId[0] if event.ngenZ>0 else float('nan'),
      binning=[50,-25,25], 
    ))

    plots.append(Plot( name = "Z0_pt_ljm"+postfix,
      texX = 'p_{T}(Z_{0}) (GeV)', texY = 'Number of Events',
      attribute = lambda event, sample: event.genZ_pt[0] if event.ngenZ>0 and abs(event.genZ_mother_pdgId[0]) in [1,2,3,4] else float('nan'),
      binning=[600/20,0,600],
    ))

    plots.append(Plot( name = "Z0_eta_ljm"+postfix,
      texX = '#eta(Z_{0}) (GeV)', texY = 'Number of Events',
      attribute = lambda event, sample: event.genZ_eta[0] if event.ngenZ>0 and abs(event.genZ_mother_pdgId[0]) in [1,2,3,4] else float('nan'),
      binning=[30,-3,3],
    ))

    plots.append(Plot( name = "Z0_pt_nljm"+postfix,
      texX = 'p_{T}(Z_{0}) (GeV)', texY = 'Number of Events',
      attribute = lambda event, sample: event.genZ_pt[0] if event.ngenZ>0 and abs(event.genZ_mother_pdgId[0]) not in [1,2,3,4] else float('nan'),
      binning=[600/20,0,600],
    ))

    plots.append(Plot( name = "Z0_eta_nljm"+postfix,
      texX = '#eta(Z_{0}) (GeV)', texY = 'Number of Events',
      attribute = lambda event, sample: event.genZ_eta[0] if event.ngenZ>0 and abs(event.genZ_mother_pdgId[0]) not in [1,2,3,4] else float('nan'),
      binning=[30,-3,3],
    ))

if 'W' in objects:
    plots.append(Plot( name = "W0_pt"+postfix,
      texX = 'p_{T}(W_{0}) (GeV)', texY = 'Number of Events',
      attribute = lambda event, sample: event.genW_pt[0] if event.ngenW>0 else float('nan'),
      binning=[600/20,0,600],
    ))

    plots.append(Plot( name = "W0_eta"+postfix,
      texX = '#eta(W_{0}) (GeV)', texY = 'Number of Events',
      attribute = lambda event, sample: event.genW_eta[0] if event.ngenW>0 else float('nan'),
      binning=[30,-3,3],
    ))

    plots.append(Plot( name = "W1_pt"+postfix,
      texX = 'p_{T}(W_{1}) (GeV)', texY = 'Number of Events',
      attribute = lambda event, sample: event.genW_pt[1] if event.ngenW>1 else float('nan'),
      binning=[600/20,0,600],
    ))

    plots.append(Plot( name = "W1_eta"+postfix,
      texX = '#eta(W_{1}) (GeV)', texY = 'Number of Events',
      attribute = lambda event, sample: event.genW_eta[1] if event.ngenW>1 else float('nan'),
      binning=[30,-3,3],
    ))

    plots.append(Plot( name = "W0_mother_pdgId"+postfix,
      texX = 'id(W_{0}) (GeV)', texY = 'Number of Events',
      attribute = lambda event, sample: event.genW_mother_pdgId[0] if event.ngenW>0 else float('nan'),
      binning=[50,-25,25],
    ))

    plots.append(Plot( name = "W0_pt_ljm"+postfix,
      texX = 'p_{T}(W_{0}) (GeV)', texY = 'Number of Events',
      attribute = lambda event, sample: event.genW_pt[0] if event.ngenW>0 and abs(event.genW_mother_pdgId[0]) in [1,2,3,4] else float('nan'),
      binning=[600/20,0,600],
    ))

    plots.append(Plot( name = "W0_eta_ljm"+postfix,
      texX = '#eta(W_{0}) (GeV)', texY = 'Number of Events',
      attribute = lambda event, sample: event.genW_eta[0] if event.ngenW>0 and abs(event.genW_mother_pdgId[0]) in [1,2,3,4] else float('nan'),
      binning=[30,-3,3],
    ))

    plots.append(Plot( name = "W0_pt_nljm"+postfix,
      texX = 'p_{T}(W_{0}) (GeV)', texY = 'Number of Events',
      attribute = lambda event, sample: event.genW_pt[0] if event.ngenW>0 and abs(event.genW_mother_pdgId[0]) not in [1,2,3,4] else float('nan'),
      binning=[600/20,0,600],
    ))

    plots.append(Plot( name = "W0_eta_nljm"+postfix,
      texX = '#eta(W_{0}) (GeV)', texY = 'Number of Events',
      attribute = lambda event, sample: event.genW_eta[0] if event.ngenW>0 and abs(event.genW_mother_pdgId[0]) not in [1,2,3,4] else float('nan'),
      binning=[30,-3,3],
    ))

if 't' in objects:
    plots.append(Plot( name = "top0_pt"+postfix,
      texX = 'p_{T}(top_{0}) (GeV)', texY = 'Number of Events',
      attribute = lambda event, sample: event.genTop_pt[0] if event.ngenTop>0 else float('nan'),
      binning=[600/20,0,600],
    ))

    plots.append(Plot( name = "top0_eta"+postfix,
      texX = '#eta(top_{0}) (GeV)', texY = 'Number of Events',
      attribute = lambda event, sample: event.genTop_eta[0] if event.ngenTop>0 else float('nan'),
      binning=[30,-3,3],
    ))

    plots.append(Plot( name = "top1_pt"+postfix,
      texX = 'p_{T}(top_{1}) (GeV)', texY = 'Number of Events',
      attribute = lambda event, sample: event.genTop_pt[1] if event.ngenTop>1 else float('nan'),
      binning=[600/20,0,600],
    ))

    plots.append(Plot( name = "top1_eta"+postfix,
      texX = '#eta(top_{1}) (GeV)', texY = 'Number of Events',
      attribute = lambda event, sample: event.genTop_eta[1] if event.ngenTop>1 else float('nan'),
      binning=[30,-3,3],
    ))

plots.append(Plot( name = "j0_pt"+postfix,
  texX = 'p_{T}(j_{0}) (GeV)', texY = 'Number of Events',
  attribute = lambda event, sample: event.jets[0]['pt'] if len(event.jets)>0 else float('nan'),
  binning=[600/20,0,600],
))

plots.append(Plot( name = "j1_pt"+postfix,
  texX = 'p_{T}(j_{1}) (GeV)', texY = 'Number of Events',
  attribute = lambda event, sample: event.jets[1]['pt'] if len(event.jets)>1 else float('nan'),
  binning=[600/20,0,600],
))

plots.append(Plot( name = "j2_pt"+postfix,
  texX = 'p_{T}(j_{2}) (GeV)', texY = 'Number of Events',
  attribute = lambda event, sample: event.jets[2]['pt'] if len(event.jets)>2 else float('nan'),
  binning=[600/20,0,600],
))

plots.append(Plot( name = "j0_eta"+postfix,
  texX = '#eta(j_{0}) (GeV)', texY = 'Number of Events',
  attribute = lambda event, sample: event.jets[0]['eta'] if len(event.jets)>0 else float('nan'),
  binning=[30,-3,3],
))

plots.append(Plot( name = "j1_eta"+postfix,
  texX = '#eta(j_{1}) (GeV)', texY = 'Number of Events',
  attribute = lambda event, sample: event.jets[1]['eta'] if len(event.jets)>1 else float('nan'),
  binning=[30,-3,3],
))

plots.append(Plot( name = "j2_eta"+postfix,
  texX = '#eta(j_{2}) (GeV)', texY = 'Number of Events',
  attribute = lambda event, sample: event.jets[2]['eta'] if len(event.jets)>2 else float('nan'),
  binning=[30,-3,3],
))

plots.append(Plot( name = "b0_pt"+postfix,
  texX = 'p_{T}(b_{0}) (GeV)', texY = 'Number of Events',
  attribute = lambda event, sample: event.trueBjets[0]['pt'] if len(event.trueBjets)>0 else float('nan'),
  binning=[600/20,0,600],
))

plots.append(Plot( name = "b1_pt"+postfix,
  texX = 'p_{T}(b_{1}) (GeV)', texY = 'Number of Events',
  attribute = lambda event, sample: event.trueBjets[1]['pt'] if len(event.trueBjets)>1 else float('nan'),
  binning=[600/20,0,600],
))

plots.append(Plot( name = "b0_eta"+postfix,
  texX = '#eta(b_{0}) (GeV)', texY = 'Number of Events',
  attribute = lambda event, sample: event.trueBjets[0]['eta'] if len(event.trueBjets)>0 else float('nan'),
  binning=[30,-3,3],
))

plots.append(Plot( name = "b1_eta"+postfix,
  texX = '#eta(b_{1}) (GeV)', texY = 'Number of Events',
  attribute = lambda event, sample: event.trueBjets[1]['eta'] if len(event.trueBjets)>1 else float('nan'),
  binning=[30,-3,3],
))

plots.append(Plot( name = 'Met_pt'+postfix,
  texX = 'E_{T}^{miss} (GeV)', texY = 'Number of Events / 20 GeV',
  attribute = lambda event, sample: event.genMet_pt,
  binning=[400/20,0,400],
))

plots.append(Plot( name = 'nJet'+postfix,
  texX = 'jet multiplicity', texY = 'Number of Events / 20 GeV',
  attribute = lambda event, sample: len(event.jets),
  binning=[8,0,8],
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
    plot_directory_ = os.path.join(plot_directory, args.plot_directory, args.sample, subDirectory)
    plot_directory_ = os.path.join(plot_directory_, "log") if log else os.path.join(plot_directory_, "lin")
    for plot in plots:
      if not max(l[0].GetMaximum() for l in plot.histos): continue # Empty plot

      len_FI = len(plot.fisher_plots) if hasattr(plot, "fisher_plots") else 0
      plotting.draw(plot,
	    plot_directory = plot_directory_,
	    ratio = {'histos':[(i,-1) for i in range(0,-1+len(plot.histos))], 'yRange':(0.1,1.9)},
	    logX = False, logY = log, sorting = False,
	    yRange = (0.03, "auto") if log else "auto",
	    #scaling = {i:len(plot.histos) for i in range(len(plot.histos)-1)},
	    legend =  ( (0.17,0.9-0.05*sum(map(len, plot.histos))/2,1.,0.9), 2),
	    drawObjects = drawObjects( ),
        copyIndexPHP = True,
      )

plotting.fill(plots+fisher_plots, read_variables = read_variables, sequence = sequence, max_events = -1 if args.small else -1)

for plot in plots:
    for i_h, hl in enumerate(plot.histos):
        # dress up
        hl[0].legendText = params[i_h]['legendText'] 
        hl[0].style = params[i_h]['style']

drawPlots(plots, subDirectory = subDirectory)

logger.info( "Done with prefix %s and selectionString %s", args.selection, cutInterpreter.cutString(args.selection) )

syncer.sync()
