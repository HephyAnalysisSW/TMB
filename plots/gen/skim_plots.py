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

# RootTools
from RootTools.core.standard        import *

# TMB
from TMB.Tools.user                 import plot_directory
from TMB.Tools.helpers              import deltaPhi, getCollection, deltaR, mZ
from TMB.Tools.WeightInfo           import WeightInfo
from TMB.Tools.genCutInterpreter       import genCutInterpreter
from TMB.Tools.VV_angles            import VV_angles
from TMB.Tools.objectSelection      import isGoodGenJet

# Arguments
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',           action='store',      default='INFO',          nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'], help="Log level for logging")
argParser.add_argument('--plot_directory',     action='store',      default='gen')
argParser.add_argument('--selection',          action='store',      default=None)
argParser.add_argument('--sample',             action='store',      default='WZ')
argParser.add_argument('--WC',                 action='store',      default='cWWW')
argParser.add_argument('--WCval',              action='store',      nargs = '*',             type=float,    default=[1.0],  help='Values of the Wilson coefficient')
argParser.add_argument('--WCval_FI',           action='store',      nargs = '*',             type=float,    default=[0.0],  help='Values of the Wilson coefficient to show FI for.')
argParser.add_argument('--FI_thresholds',      action='store',      nargs = '*',             type=float,    default=None,   help='Bins in FI?')
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
from TMB.Samples.gen_pp import *
sample = eval(args.sample) 
# objects to plot
objects = sample.objects if hasattr( sample, "objects") else []

# WeightInfo
w = WeightInfo(sample.reweight_pkl)
w.set_order(2)

# define which Wilson coefficients to plot
FIs = []
if args.FI_thresholds is None:
    params =  [ {'legendText':'SM',  'color':ROOT.kBlue, 'WC':{}} ] 
    params += [ {'legendText':'%s %3.2f'%(args.WC, wc), 'color':ROOT.kOrange+i_wc,  'WC':{args.WC:wc} } for i_wc, wc in enumerate(args.WCval)]

    for i_param, param in enumerate(params):
        param['sample'] = sample
        param['style']  = styles.lineStyle( param['color'] )

    stack = Stack(*[ [ param['sample'] ] for param in params ] )
    weight= [ [ w.get_weight_func(**param['WC']) ] for param in params ]
    FIs.append( ( ROOT.kGray+1, "WC@SM",   w.get_fisher_weight_string(args.WC,args.WC, **{args.WC:0})) )
    for i_WCval, WCval in enumerate(args.WCval_FI):
        FIs.append( ( ROOT.kGray+i_WCval,   "WC@WC=%3.2f"%WCval, w.get_fisher_weight_string(args.WC,args.WC, **{args.WC:WCval})) )
else:

    colors = [ ROOT.kViolet-9, ROOT.kViolet-6, ROOT.kViolet-4, ROOT.kViolet-3 , ROOT.kViolet-2, ROOT.kViolet]
    WCval = args.WCval[0]
    if len(args.WCval)>1:
        raise RuntimeError, "Call with exactly one WC!"
    var = w.get_fisher_weight_string(args.WC,args.WC, **{args.WC:WCval})

    params = [ {'legendText':'I_{F}<%3.2f'%args.FI_thresholds[0],'selection':"{var}<%f"%(args.FI_thresholds[0]), 'style':styles.fillStyle(colors[0])} ]
    for i_th, th in enumerate( args.FI_thresholds[:-1] ):
        params.append( { 'legendText':'%3.2f<=I_{F}<%3.2f'%(args.FI_thresholds[i_th], args.FI_thresholds[i_th+1]), 'selection':"%f<={var}&&{var}<%f"%(args.FI_thresholds[i_th], args.FI_thresholds[i_th+1]) , 'style':styles.fillStyle(colors[i_th+1])})
    params.append( {'legendText':'I_{F}>=%3.2f'%args.FI_thresholds[-1], 'selection': "%f<={var}"%(args.FI_thresholds[-1]), 'style':styles.fillStyle(colors[len(args.FI_thresholds)])} )
     
    stack = Stack([ copy.deepcopy(sample) for param in params ] )
    for i_s_, s_ in enumerate(stack[0]):
        s_.setSelectionString( params[i_s_]['selection'].format(var=var) )
        s_.name+='_%i'%i_s_
    weight= [ [ w.get_weight_func(**{args.WC:WCval}) for param in params ] ]
    

# Read variables and sequences
read_variables = [
    "genMet_pt/F", "genMet_phi/F", 
    "ngenJet/I", "genJet[pt/F,eta/F,phi/F,matchBParton/I]", 
    "ngenLep/I", "genLep[pt/F,eta/F,phi/F,pdgId/I,mother_pdgId/I]", 
    "ngenTop/I", "genTop[pt/F,eta/F,phi/F]", 
    "ngenZ/I", "genZ[pt/F,phi/F,eta/F,daughter_pdgId/I,l1_index/I,l2_index/I]",
    "ngenW/I", "genW[pt/F,phi/F,eta/F,daughter_pdgId/I,l1_index/I,l2_index/I]",
    "ngenPhoton/I", "genPhoton[pt/F,phi/F,eta/F]",
]
read_variables.append( VectorTreeVariable.fromString('p[C/F]', nMax=2000) )

preselection = [ 
    #( "SMP-20-005",  "genPhoton_pt[0]>300&&genMet_pt>80&&genLep_pt[genW_l1_index[0]]>80&&sqrt(acos(cos(genLep_phi[genW_l1_index[0]]-genPhoton_phi[0]))**2+(genLep_eta[genW_l1_index[0]]-genPhoton_eta[0])**2)>3.0"),
    #( "SMP-20-005-light",  "genPhoton_pt[0]>150&genMet_pt>30&&genLep_pt[genW_l1_index[0]]>30&&sqrt(acos(cos(genLep_phi[genW_l1_index[0]]-genPhoton_phi[0]))**2+(genLep_eta[genW_l1_index[0]]-genPhoton_eta[0])**2)>3.0"),
]

selectionString  = "&&".join( [ c[1] for c in preselection] + ([genCutInterpreter.cutString(args.selection)] if args.selection is not None else []))
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
    def make_VV_angles( event, sample ):
        # nothing to resurrect here
        event.genZs = filter( lambda z: abs(z['daughter_pdgId']) in [11, 13, 15], getCollection( event, 'genZ', ['pt', 'eta', 'phi', 'daughter_pdgId', 'l1_index', 'l2_index'], 'ngenZ' ) )
        if len(event.genZs)>1:
            id1 = event.genZs[0]['l1_index']
            id2 = event.genZs[0]['l2_index']
            id3 = event.genZs[1]['l1_index']
            id4 = event.genZs[1]['l2_index']
            v = [ ROOT.TLorentzVector() for i in range(4) ]
            l = event.all_leps
            for i_id, id_ in enumerate( [id1, id2, id3, id4] ):
                #print "Lepton",i_id, l[id_]
                v[i_id].SetPtEtaPhiM(l[id_]['pt'], l[id_]['eta'], l[id_]['phi'], 0)
            #print
            event.ZZ_angles = VV_angles( *v, debug=False)
        else:
            event.ZZ_angles = None
        #print event.ZZ_angles
        #print
        #assert False, "" 
    sequence.append( make_VV_angles )

if args.WG:
    def make_VV_angles( event, sample ):
        # SMP-20-005 (gen-level version) 
        #print  getCollection( event, 'genW', ['pt', 'eta', 'phi', 'daughter_pdgId', 'l1_index', 'l2_index'], 'ngenW' )
        event.genWs      = filter( lambda w: abs(w['daughter_pdgId']) in [11, 13, 15], getCollection( event, 'genW', ['pt', 'eta', 'phi', 'daughter_pdgId', 'l1_index', 'l2_index'], 'ngenW' ) )
        event.genPhotons = getCollection( event, 'genPhoton', ['pt', 'eta', 'phi'], 'ngenPhoton' )
        if len(event.genWs)>0 and len(event.genPhotons)>0:
            genW = event.genWs[0]
            genPhoton = event.genPhotons[0]

            l = event.all_leps

            v = [ ROOT.TLorentzVector() for i in range(2) ]
            id1, id2 = (event.genWs[0]['l1_index'], event.genWs[0]['l2_index']) if event.genLep_pdgId[event.genWs[0]['l1_index']]>0 else (event.genWs[0]['l2_index'], event.genWs[0]['l1_index'])

            for i_id, id_ in enumerate( [id1, id2] ):
                v[i_id].SetPtEtaPhiM(l[id_]['pt'], l[id_]['eta'], l[id_]['phi'], 0)
            v_photon = ROOT.TLorentzVector()
            v_photon.SetPtEtaPhiM(genPhoton['pt'], genPhoton['eta'],genPhoton['phi'], 0)
            v.append( v_photon )
            v.append( ROOT.TLorentzVector(0,0,0,0) )

            event.WG_angles = VV_angles( *v, debug=False)
        else:
            event.WG_angles = None
        #print event.WG_angles
        #print 
    sequence.append( make_VV_angles )

import re
# make stack and construct weights
def coeff_getter( coeff):
    def getter_( event, sample):
        #print "Filling with", coeff, event.p_C[coeff]
        return event.p_C[coeff]
    return getter_

def add_fisher_plot( plot, fisher_string, color, legendText):
    # the coefficients we need
    required_coefficients = map( lambda s: int(s.replace('[','').replace(']','')), list(set(re.findall(r'\[[0-9][0-9]*\]', fisher_string ) )) )
    # copy the plot
    fisher_plot = copy.deepcopy( plot )
    # modify weights etc
    fisher_plot.name+="_fisher_coeffs"
    fisher_plot.weight = [ [ coeff_getter(coeff) ] for coeff in required_coefficients ]
    fisher_plot.stack  = Stack( *[[sample] for coeff in required_coefficients] )
    # for computing the FI
    fisher_plot.coefficients = required_coefficients
    fisher_plot.formula      = re.sub(r'p_C\[([0-9][0-9]*)\]', r'{lambda_\1}', fisher_string)
    # pass through information for plotting
    fisher_plot.color        = color
    fisher_plot.legendText   = legendText
    # add the fisher_plot to the original plot so we know which is which
    if hasattr( plot, "fisher_plots"):
        plot.fisher_plots.append( fisher_plot )
    else: 
        plot.fisher_plots = [fisher_plot]
    return fisher_plot 
    
# Use some defaults
Plot.setDefaults(stack = stack, weight = weight, addOverFlowBin=None)
  
plots        = []
fisher_plots = []

postfix = '_'+args.WC + ("_FIth" if args.FI_thresholds is not None else "") 

if args.ZZ:
    plots.append(Plot( name = "phi1"+postfix,
      texX = '#phi_{1}', texY = 'Number of Events',
      attribute = lambda event, sample: event.ZZ_angles['phi1'] if event.ZZ_angles is not None else float('nan'),
      binning=[40, -pi, pi], addOverFlowBin="both",
    ))
    for color, legendText, fisher_string in FIs:
        fisher_plots.append( add_fisher_plot( plots[-1], fisher_string, color, legendText ) )

    plots.append(Plot( name = "phi2"+postfix,
      texX = '#phi_{2}', texY = 'Number of Events',
      attribute = lambda event, sample: event.ZZ_angles['phi2'] if event.ZZ_angles is not None else float('nan'),
      binning=[40, -pi, pi], addOverFlowBin="both",
    ))
    for color, legendText, fisher_string in FIs:
        fisher_plots.append( add_fisher_plot( plots[-1], fisher_string, color, legendText ) )

    plots.append(Plot( name = "deltaPhi"+postfix,
      texX = '#Delta #Phi', texY = 'Number of Events',
      attribute = lambda event, sample: event.ZZ_angles['deltaPhi'] if event.ZZ_angles is not None else float('nan'),
      binning=[40, 0, 2*pi], addOverFlowBin="both",
    ))
    for color, legendText, fisher_string in FIs:
        fisher_plots.append( add_fisher_plot( plots[-1], fisher_string, color, legendText ) )

    plots.append(Plot( name = "theta_V1"+postfix,
      texX = '#theta(V1)', texY = 'Number of Events',
      attribute = lambda event, sample: event.ZZ_angles['theta_V1'] if event.ZZ_angles is not None else float('nan'),
      binning=[20, 0, pi], addOverFlowBin="both",
    ))
    for color, legendText, fisher_string in FIs:
        fisher_plots.append( add_fisher_plot( plots[-1], fisher_string, color, legendText ) )

    plots.append(Plot( name = "theta_V2"+postfix,
      texX = '(#theta(V2)', texY = 'Number of Events',
      attribute = lambda event, sample: event.ZZ_angles['theta_V2'] if event.ZZ_angles is not None else float('nan'),
      binning=[20, 0, pi], addOverFlowBin="both",
    ))
    for color, legendText, fisher_string in FIs:
        fisher_plots.append( add_fisher_plot( plots[-1], fisher_string, color, legendText ) )

    plots.append(Plot( name = "theta"+postfix,
      texX = '#theta', texY = 'Number of Events',
      attribute = lambda event, sample: event.ZZ_angles['theta'] if event.ZZ_angles is not None else float('nan'),
      binning=[20, 0,pi], addOverFlowBin="both",
    ))
    for color, legendText, fisher_string in FIs:
        fisher_plots.append( add_fisher_plot( plots[-1], fisher_string, color, legendText ) )

if 'g' in objects:
    plots.append(Plot( name = "Photon0_pt"+postfix,
      texX = 'p_{T}(#gamma_{0}) (GeV)', texY = 'Number of Events',
      attribute = lambda event, sample: event.genPhoton_pt[0] if event.ngenPhoton>0 else float('nan'),
      binning=[600/20,0,600],
    ))
    for color, legendText, fisher_string in FIs:
        fisher_plots.append( add_fisher_plot( plots[-1], fisher_string, color, legendText ) )

    plots.append(Plot( name = "Photon0_eta"+postfix,
      texX = '#eta(#gamma_{0}) (GeV)', texY = 'Number of Events',
      attribute = lambda event, sample: event.genPhoton_eta[0] if event.ngenPhoton>0 else float('nan'),
      binning=[30,-3,3],
    ))
    for color, legendText, fisher_string in FIs:
        fisher_plots.append( add_fisher_plot( plots[-1], fisher_string, color, legendText ) )

if args.WG:
    plots.append(Plot( name = "phi1"+postfix,
      texX = '#phi_{1}', texY = 'Number of Events',
      attribute = lambda event, sample: event.WG_angles['phi1'] if event.WG_angles is not None else float('nan'),
      binning=[20, -pi, pi], addOverFlowBin="both",
    ))
    for color, legendText, fisher_string in FIs:
        fisher_plots.append( add_fisher_plot( plots[-1], fisher_string, color, legendText ) )

    plots.append(Plot( name = "phi2"+postfix,
      texX = '#phi_{2}', texY = 'Number of Events',
      attribute = lambda event, sample: event.WG_angles['phi2'] if event.WG_angles is not None else float('nan'),
      binning=[20, -pi, pi], addOverFlowBin="both",
    ))
    for color, legendText, fisher_string in FIs:
        fisher_plots.append( add_fisher_plot( plots[-1], fisher_string, color, legendText ) )

    plots.append(Plot( name = "deltaPhi"+postfix,
      texX = '#Delta #Phi', texY = 'Number of Events',
      attribute = lambda event, sample: event.WG_angles['deltaPhi'] if event.WG_angles is not None else float('nan'),
      binning=[20, 0, 2*pi], addOverFlowBin="both",
    ))
    for color, legendText, fisher_string in FIs:
        fisher_plots.append( add_fisher_plot( plots[-1], fisher_string, color, legendText ) )

    plots.append(Plot( name = "theta_V1"+postfix,
      texX = '#theta(V1)', texY = 'Number of Events',
      attribute = lambda event, sample: event.WG_angles['theta_V1'] if event.WG_angles is not None else float('nan'),
      binning=[20, 0, pi], addOverFlowBin="both",
    ))
    for color, legendText, fisher_string in FIs:
        fisher_plots.append( add_fisher_plot( plots[-1], fisher_string, color, legendText ) )

    plots.append(Plot( name = "theta_V2"+postfix,
      texX = '(#theta(V2)', texY = 'Number of Events',
      attribute = lambda event, sample: event.WG_angles['theta_V2'] if event.WG_angles is not None else float('nan'),
      binning=[20, 0, pi], addOverFlowBin="both",
    ))
    for color, legendText, fisher_string in FIs:
        fisher_plots.append( add_fisher_plot( plots[-1], fisher_string, color, legendText ) )

    plots.append(Plot( name = "theta"+postfix,
      texX = '#theta', texY = 'Number of Events',
      attribute = lambda event, sample: event.WG_angles['theta'] if event.WG_angles is not None else float('nan'),
      binning=[20, 0, pi], addOverFlowBin="both",
    ))
    for color, legendText, fisher_string in FIs:
        fisher_plots.append( add_fisher_plot( plots[-1], fisher_string, color, legendText ) )

if 'Z' in objects:
    plots.append(Plot( name = "Z0_pt"+postfix,
      texX = 'p_{T}(Z_{0}) (GeV)', texY = 'Number of Events',
      attribute = lambda event, sample: event.genZ_pt[0] if event.ngenZ>0 else float('nan'),
      binning=[600/20,0,600],
    ))
    for color, legendText, fisher_string in FIs:
        fisher_plots.append( add_fisher_plot( plots[-1], fisher_string, color, legendText ) )

    plots.append(Plot( name = "Z0_eta"+postfix,
      texX = '#eta(Z_{0}) (GeV)', texY = 'Number of Events',
      attribute = lambda event, sample: event.genZ_eta[0] if event.ngenZ>0 else float('nan'),
      binning=[30,-3,3],
    ))
    for color, legendText, fisher_string in FIs:
        fisher_plots.append( add_fisher_plot( plots[-1], fisher_string, color, legendText ) )

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

if 'W' in objects:
    plots.append(Plot( name = "W0_pt"+postfix,
      texX = 'p_{T}(W_{0}) (GeV)', texY = 'Number of Events',
      attribute = lambda event, sample: event.genW_pt[0] if event.ngenW>0 else float('nan'),
      binning=[600/20,0,600],
    ))
    for color, legendText, fisher_string in FIs:
        fisher_plots.append( add_fisher_plot( plots[-1], fisher_string, color, legendText ) )

    plots.append(Plot( name = "W0_eta"+postfix,
      texX = '#eta(W_{0}) (GeV)', texY = 'Number of Events',
      attribute = lambda event, sample: event.genW_eta[0] if event.ngenW>0 else float('nan'),
      binning=[30,-3,3],
    ))
    for color, legendText, fisher_string in FIs:
        fisher_plots.append( add_fisher_plot( plots[-1], fisher_string, color, legendText ) )

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

if 't' in objects:
    plots.append(Plot( name = "top0_pt"+postfix,
      texX = 'p_{T}(top_{0}) (GeV)', texY = 'Number of Events',
      attribute = lambda event, sample: event.genTop_pt[0] if event.ngenTop>0 else float('nan'),
      binning=[600/20,0,600],
    ))
    for color, legendText, fisher_string in FIs:
        fisher_plots.append( add_fisher_plot( plots[-1], fisher_string, color, legendText ) )

    plots.append(Plot( name = "top0_eta"+postfix,
      texX = '#eta(top_{0}) (GeV)', texY = 'Number of Events',
      attribute = lambda event, sample: event.genTop_eta[0] if event.ngenTop>0 else float('nan'),
      binning=[30,-3,3],
    ))
    for color, legendText, fisher_string in FIs:
        fisher_plots.append( add_fisher_plot( plots[-1], fisher_string, color, legendText ) )

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
	    ratio = {'histos':[(i,0) for i in range(1,len(plot.histos)-len_FI)], 'yRange':(0.1,1.9)} if args.FI_thresholds is None else None,
	    logX = False, logY = log, sorting = False,
	    yRange = (0.03, "auto") if log else "auto",
	    scaling = {},
	    legend =  ( (0.17,0.9-0.05*sum(map(len, plot.histos))/2,1.,0.9), 2),
	    drawObjects = drawObjects( ),
        copyIndexPHP = True,
      )

plotting.fill(plots+fisher_plots, read_variables = read_variables, sequence = sequence, max_events = -1 if args.small else -1)

for plot in plots:
    if args.FI_thresholds is None:
        for i_h, hl in enumerate(plot.histos):
            # dress up
            hl[0].legendText = params[i_h]['legendText'] 
            hl[0].style = params[i_h]['style']
    else: 
        for i_h, hl in enumerate(plot.histos[0]):
            # dress up
            hl.legendText = params[i_h]['legendText'] 
            hl.style = params[i_h]['style']

    # calculate & add FI histo
    if hasattr(plot, "fisher_plots" ):
        for fisher_plot in plot.fisher_plots:

            # make empty histo
            FI = plot.histos[0][0].Clone()
            FI.Reset()

            # calculate FI in each bin
            for i_bin in range(1, FI.GetNbinsX()+1):
                yields = {'lambda_%i'%coeff : fisher_plot.histos[i_coeff][0].GetBinContent(i_bin) for i_coeff, coeff in enumerate(fisher_plot.coefficients)}
                try:
                    FI.SetBinContent( i_bin,  eval(fisher_plot.formula.format(**yields)) )
                except ZeroDivisionError:
                    pass

            # dress up
            FI.legendText = fisher_plot.legendText
            FI.style = styles.lineStyle( fisher_plot.color )

            # scale the FI & indicate log(I) in the plot
            if FI.Integral()>0:
                FI.legendText += "(%3.2f)" % log(FI.Integral(),10)
                FI.Scale(plot.histos[0][0].Integral()/FI.Integral())

            # add the FI histos back to the plot and fake the stack so that the draw command does not get confused
            plot.histos.append( [FI] )
            stack.append( [sample] )

drawPlots(plots, subDirectory = subDirectory)

logger.info( "Done with prefix %s and selectionString %s", args.selection, genCutInterpreter.cutString(args.selection) )
