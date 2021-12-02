#!/usr/bin/env python

# Standard imports
from operator                   import attrgetter
from math                       import pi, sqrt, cosh, cos, acos, sin, sinh, copysign
import ROOT, os, copy

# RootTools
from RootTools.core.standard     import *

# helpers
from Analysis.Tools.helpers              import deltaPhi, deltaR, getObjDict
from TMB.Tools.helpers          import getCollection

from Analysis.Tools.WeightInfo       import WeightInfo

import logging
logger = logging.getLogger(__name__)

# read_variables
jetVars          = ['pt/F', 'eta/F', 'phi/F', 'bTag/F', 'bTagPhys/I']
jetVarNames      = [x.split('/')[0] for x in jetVars]
lepVars          = ['pt/F','eta/F','phi/F','pdgId/I','isolationVar/F', 'isolationVarRhoCorr/F']
lepVarNames      = [x.split('/')[0] for x in lepVars]
read_variables = [\
    "nBTag_loose/I",
    "recoMet_pt/F", "recoMet_phi/F",
    "genMet_pt/F", "genMet_phi/F",
    "nrecoJet/I",
    "recoJet[%s]"%(",".join(jetVars)),
    "nrecoLep/I",
    "recoLep[%s]"%(",".join(lepVars)),
    "ngenLep/I", "genLep[pt/F,eta/F,phi/F,pdgId/I,mother_pdgId/I]",
    "lumiweight1fb/F",
    "genW[pt/F,eta/F,phi/F,l1_index/I,l2_index/I]", "ngenW/I",
    "evt/l", "run/I", "lumi/I", 
    #"np/I", "nweight/I",
    "H_dijet_mass/F", "H_pt/F", "H_j1_index/I", "H_j2_index/I", 
    "H_j1_index/I", "H_j2_index/I", 
    "WH_W_pt/F", "WH_dPhiMetLep/F", "WH_MT/F", "WH_nu_pt/F", "WH_nu_eta/F", "WH_nu_phi/F", "WH_nu_E/F", "WH_Theta/F", "WH_theta/F", "WH_phi/F", 
]

def addTLorentzVector( p_dict ):
    ''' add a TLorentz 4D Vector for further calculations
    '''
    p_dict['vecP4'] = ROOT.TLorentzVector( p_dict['pt']*cos(p_dict['phi']), p_dict['pt']*sin(p_dict['phi']),  p_dict['pt']*sinh(p_dict['eta']), p_dict['pt']*cosh(p_dict['eta']) )

sequence = []

from TMB.Tools.objectSelection import isBJet
def makeJets( event, sample ):
    event.jets     = [getObjDict(event, 'recoJet_', jetVarNames, i) for i in range(int(event.nrecoJet))]
    for p in event.jets:
        #addTransverseVector( p )
        addTLorentzVector( p )
    event.bJets      = filter(lambda j:j['bTag']>=1 and abs(j['eta'])<=2.4, event.jets)
    event.extraJets  = [ j for i_j, j in enumerate( event.jets ) if i_j not in [event.H_j1_index, event.H_j2_index] ] 
    event.nextraJet = len(event.extraJets)

sequence.append( makeJets )

def makeLeptonic( event, sample):
    all_leps = getCollection( event, 'recoLep', ['pt', 'eta', 'phi', 'pdgId', 'isolationVar', 'isolationVarRhoCorr'], 'nrecoLep' )
    # Add extra vectors
    for p in all_leps:
        #addTransverseVector( p )
        addTLorentzVector( p )
    event.lepton = all_leps[0]

sequence.append( makeLeptonic )

def makeAngles(event, sample):

    if event.WH_theta < float('inf'):

        #print event.WH_Theta, event.WH_theta, event.WH_phi, event.WH_nu_pt, event.WH_nu_eta
        event.sin2thetaSin2Theta = sin(2*event.WH_theta)*sin(2*event.WH_Theta)

        event.fLL         = sin(event.WH_Theta)**2*sin(event.WH_theta)**2
        event.f1TT        = cos(event.WH_theta)*cos(event.WH_Theta)
        event.f2TT        = (1+cos(event.WH_theta)**2)*(1+cos(event.WH_Theta)**2)
        event.f1LT        = cos(event.WH_phi)*sin(event.WH_theta)*sin(event.WH_Theta)
        event.f2LT        = event.f1LT*cos(event.WH_theta)*cos(event.WH_Theta)
        event.f1tildeLT   = sin(event.WH_phi)*sin(event.WH_theta)*sin(event.WH_Theta)
        event.f2tildeLT   = event.f1tildeLT*cos(event.WH_theta)*cos(event.WH_Theta)
        event.fTTprime    = cos(2*event.WH_phi)*event.fLL
        event.ftildeTTprime=sin(2*event.WH_phi)*event.fLL
    else:
        event.WH_theta    = -2
        event.WH_phi      = -4

        event.fLL         = -2 
        event.f1TT        = -2 
        event.f2TT        = -2 
        event.f1LT        = -2 
        event.f2LT        = -2 
        event.f1tildeLT   = -2 
        event.f2tildeLT   = -2 
        event.fTTprime    = -2 
        event.ftildeTTprime=-2

sequence.append( makeAngles )

def makeH( event, sample ):
    event.dijet_mass = (event.bJets[0]['vecP4'] + event.bJets[1]['vecP4']).M()
    event.H_j1    = event.jets[event.H_j1_index]
    event.H_j2    = event.jets[event.H_j2_index]
    event.H_vecP4 = event.H_j1['vecP4'] + event.H_j2['vecP4']
    event.H_pt    = event.H_vecP4.Pt()
    event.H_eta   = event.H_vecP4.Eta()
    event.H_phi   = event.H_vecP4.Phi()

sequence.append( makeH )

from TMB.Tools.Thrust import Thrust
def makeThrust( event, sample ):

    event.W_nu_vecP4 = ROOT.TLorentzVector()
    event.W_nu_vecP4.SetPtEtaPhiM(event.WH_nu_pt, event.WH_nu_eta, event.WH_nu_phi, 0 )

    event.W_vecP4 = event.W_nu_vecP4 + event.lepton['vecP4'] 
    (event.thrust, event.thrust_min) = Thrust(event.W_vecP4,[j['vecP4'] for j in event.jets] )

sequence.append( makeThrust )


weight_variables = ['cHj3', 'cHW', 'cHWtil']
max_order        = 2

from TMB.Samples.pp_gen_v10 import *
training_samples = [ WH, TTJets, WJetsToLNu_HT]

assert len(training_samples)==len(set([s.name for s in training_samples])), "training_samples names are not unique!"

# Add weight infos
for sample in training_samples:
    if hasattr( sample, "reweight_pkl" ):
        sample.weightInfo = WeightInfo( sample.reweight_pkl )
        sample.weightInfo.set_order(2)
        sample.read_variables = [
                VectorTreeVariable.fromString( "p[C/F]", nMax=200 ),
                VectorTreeVariable.fromString( "weight[base/F]", nMax=200 )
            ]

# make all combinations
weight_derivative_combinations = []
for i_comb, comb in enumerate(WH.weightInfo.make_combinations(weight_variables, max_order)):
    weight_derivative_combinations.append(comb)

scale_weight = 10**5

for sample in training_samples:
    sample.weight_derivatives = []
    for i_comb, comb in enumerate(weight_derivative_combinations):
        #print name, i_comb, comb, weightInfo.get_diff_weight_string(comb)
        weight = {}
        if hasattr( sample, "weightInfo"):
            if all( v in sample.weightInfo.variables for v in comb ):
                # func_ takes care of p_C. We also normalize with the lumi-weight 'weight' from the original sample
                #weight['string'] = "lumiweight1fb*("+sample.weightInfo.get_diff_weight_string(comb)+")"
                func_            = sample.weightInfo.get_diff_weight_func(comb)
                # evaluate weight at BSM point, divide by first base weight (SM), norm to lumi
                weight['func']   = lambda event, sample, func_=func_: func_(event,sample)*event.lumiweight1fb/event.weight_base[0]*scale_weight
            else:
                print "Warning! Derivative %r put to zero in sample %s because some WC are missing."%( comb, sample.name )
                weight = {}
                #weight['string'] = '(0.)'
                weight['func']   = lambda event, sample: 0.
        else:
            const = 1 if len(comb)==0 else 0
            print "Warning. Sample %s has no weightInfo. Therefore, weight %r will be %f." % ( sample.name, comb, const )
            weight['func']   = lambda event, sample, const=const: const*event.lumiweight1fb*scale_weight
        weight['name']   = '_'.join(comb)
        weight['comb']   = comb
        sample.weight_derivatives.append( weight )

def compute_weight_derivatives( event, sample ):
    vector = [{'derivatives':weight['func'](event, sample)} for weight in sample.weight_derivatives]

    #print vector[115]
    #m =  [abs(weight['func'](event, sample)/nominal) for weight in weights[1:]]
    #if max(m)>100:
    #    print max(m)

    return vector
    #for weight in weights[1:]:
    #    setattr( event, "locsc_"+weight['name'], 0 if nominal==0 else weight['func'](event, sample)/nominal )  

#sequence.append( compute_weight_derivatives )

all_mva_variables = {

# global event properties
     "mva_ht"                    :(lambda event, sample: sum( [event.recoJet_pt[i] for i in range(event.nrecoJet) ])),
     "mva_met_pt"                :(lambda event, sample: event.recoMet_pt),
     "mva_nextraJet"            :(lambda event, sample: event.nextraJet),
     "mva_nBTag"                 :(lambda event, sample: event.nBTag_loose),

# jet kinmatics
     "mva_jet0_pt"               :(lambda event, sample: event.recoJet_pt[0]          if event.nrecoJet >=1 else 0),
     "mva_jet0_eta"              :(lambda event, sample: event.recoJet_eta[0]         if event.nrecoJet >=1 else -10),
     "mva_jet0_btag"             :(lambda event, sample: event.recoJet_bTag[0]        if event.nrecoJet >=1 else -1),
     "mva_jet1_pt"               :(lambda event, sample: event.recoJet_pt[1]          if event.nrecoJet >=2 else 0),
     "mva_jet1_eta"              :(lambda event, sample: event.recoJet_eta[1]         if event.nrecoJet >=2 else -10),
     "mva_jet1_btag"             :(lambda event, sample: event.recoJet_bTag[1]        if event.nrecoJet >=2 else -1),
     "mva_jet2_pt"               :(lambda event, sample: event.recoJet_pt[2]          if event.nrecoJet >=3 else 0),
     "mva_jet2_eta"              :(lambda event, sample: event.recoJet_eta[2]         if event.nrecoJet >=3 else -10),
     "mva_jet2_btag"             :(lambda event, sample: event.recoJet_bTag[2]        if event.nrecoJet >=3 else -1),

# W kinematics
     "mva_W_pt"                  :(lambda event, sample: event.WH_W_pt),
     "mva_W_eta"                 :(lambda event, sample: event.W_vecP4.Eta()),
     #"mva_W_cosThetaStar"       :(lambda event, sample: event.recoZ_cosThetaStar),
     "mva_WH_dPhiMetLep"         :(lambda event, sample: event.WH_dPhiMetLep), 
     "mva_WH_MT"                 :(lambda event, sample: event.WH_MT),
     "mva_W_RatioPtlPtnu"        :(lambda event, sample: event.WH_nu_pt/event.lepton['pt']) ,
     "mva_W_lepton_pt"           :(lambda event, sample: event.lepton['pt']) ,
     "mva_W_lepton_eta"          :(lambda event, sample: event.lepton['eta']) ,
     "mva_WH_nu_pt"              :(lambda event, sample: event.WH_nu_pt) ,
     "mva_WH_nu_eta"             :(lambda event, sample: event.WH_nu_eta) ,

# W decay angles & signed
    "mva_WH_Theta"               :(lambda event, sample: event.WH_Theta),
    "mva_WH_theta"               :(lambda event, sample: event.WH_theta),
    "mva_WH_phi"                 :(lambda event, sample: event.WH_phi),
    "mva_WH_fLL"                 :(lambda event, sample: event.fLL),
    "mva_WH_f1TT"                :(lambda event, sample: event.f1TT),
    "mva_WH_f2TT"                :(lambda event, sample: event.f2TT),
    "mva_WH_f1LT"                :(lambda event, sample: event.f1LT),
    "mva_WH_f2LT"                :(lambda event, sample: event.f2LT),
    "mva_WH_f1tildeLT"           :(lambda event, sample: event.f1tildeLT),
    "mva_WH_f2tildeLT"           :(lambda event, sample: event.f2tildeLT),
    "mva_WH_fTTprime"            :(lambda event, sample: event.fTTprime),
    "mva_WH_ftildeTTprime"       :(lambda event, sample: event.ftildeTTprime),

# H kinematics
     "mva_H_pt"                  :(lambda event, sample: event.H_pt),
     "mva_H_dijet_mass"          :(lambda event, sample: event.H_dijet_mass),
     "mva_H_eta"                 :(lambda event, sample: event.H_eta ),
     "mva_H_DeltaPhib1b2"        :(lambda event, sample: deltaPhi(event.H_j1['phi'], event.H_j2['phi'])) ,
     "mva_H_RatioPtb1b2"         :(lambda event, sample: event.H_j2['pt']/event.H_j1['pt']) ,
     "mva_H_maxPtb1b2"           :(lambda event, sample: max([event.H_j1['pt'], event.H_j2['pt']])) ,
     "mva_H_minPtb1b2"           :(lambda event, sample: min([event.H_j1['pt'], event.H_j2['pt']])) ,
     "mva_H_DeltaRb1b2"          :(lambda event, sample: deltaR(event.H_j1, event.H_j2)) ,

# W vs. H
     "mva_WH_deltaPhi"           :(lambda event, sample: deltaPhi(event.H_phi, event.W_vecP4.Phi())),
     "mva_WH_deltaR"             :(lambda event, sample: deltaR({'phi':event.H_phi, 'eta':event.H_eta},{'phi': event.W_vecP4.Phi(), 'eta':event.W_vecP4.Eta()})),
     "mva_mWH"                   :(lambda event, sample: (event.W_vecP4+event.H_vecP4).M()),

# thrust
    "mva_thrust"                 :(lambda event, sample: event.thrust),
    "mva_thrust_min"             :(lambda event, sample: event.thrust_min),

# W vs. other objects
     "mva_extraJet0_W_deltaR"    :(lambda event, sample: deltaR({'phi':event.extraJets[0]['phi'], 'eta':event.extraJets[0]['eta']}, {'phi':event.W_vecP4.Phi(), 'eta':event.W_vecP4.Eta()} )  if event.nextraJet >=1 else -1),
     "mva_extraJet1_W_deltaR"    :(lambda event, sample: deltaR({'phi':event.extraJets[1]['phi'], 'eta':event.extraJets[1]['eta']}, {'phi':event.W_vecP4.Phi(), 'eta':event.W_vecP4.Eta()} )  if event.nextraJet >=2 else -1),
     #"mva_extraJet2_W_deltaR"    :(lambda event, sample: deltaR({'phi':event.extraJets[2]['phi'], 'eta':event.extraJets[2]['eta']}, {'phi':event.W_vecP4.Phi(), 'eta':event.W_vecP4.Eta()} )  if event.nextraJet >=3 else -1),

# H vs. other objects
     "mva_extraJet0_H_deltaR"    :(lambda event, sample: deltaR({'phi':event.extraJets[0]['phi'], 'eta':event.extraJets[0]['eta']}, {'phi':event.H_phi, 'eta':event.H_eta} )  if event.nextraJet >=1 else -1),
     "mva_extraJet1_H_deltaR"    :(lambda event, sample: deltaR({'phi':event.extraJets[1]['phi'], 'eta':event.extraJets[1]['eta']}, {'phi':event.H_phi, 'eta':event.H_eta} )  if event.nextraJet >=2 else -1),
     #"mva_extraJet2_H_deltaR"    :(lambda event, sample: deltaR({'phi':event.extraJets[2]['phi'], 'eta':event.extraJets[2]['eta']}, {'phi':event.H_phi, 'eta':event.H_eta} )  if event.nextraJet >=3 else -1),
}
# global event properties

Nbins = 30
plot_options = {
     "mva_ht"              :{'tex':'H_{T} (GeV)',           'binning':[20,0,1000]},
     "mva_met_pt"          :{'tex':'p_{T}^{miss}',          'binning':[20,0,400]},
     "mva_nextraJet"       :{'tex':'N_{jet, extra}',        'binning':[4,0,4]},
     "mva_nBTag"           :{'tex':'N_{b-tag}',             'binning':[4,0,4]},
     "mva_jet0_pt"         :{'tex':'p_{T}(j_{0}) (GeV)',    'binning':[600/20,0,600]},
     "mva_jet0_eta"        :{'tex':'#eta (j_{0})',          'binning':[30,-3,3]},     
     "mva_jet0_btag"       :{'tex':'b-jet disc. of j_{0}',  'binning':[5,-2,3]},
     "mva_jet1_pt"         :{'tex':'p_{T}(j_{1}) (GeV)',    'binning':[600/20,0,600]},
     "mva_jet1_eta"        :{'tex':'#eta (j_{1})',          'binning':[30,-3,3]},
     "mva_jet1_btag"       :{'tex':'b-jet disc. of j_{1}',  'binning':[5,-2,3]},
     "mva_jet2_pt"         :{'tex':'p_{T}(j_{2}) (GeV)',    'binning':[600/20,0,600]},
     "mva_jet2_eta"        :{'tex':'#eta (j_{2})',          'binning':[30,-3,3]},
     "mva_jet2_btag"       :{'tex':'b-jet disc. of j_{2}',  'binning':[5,-2,3]},
     "mva_W_pt"            :{'tex':'p_{T}(W) (GeV)',        'binning':[600/20,0,600]},       
     "mva_W_eta"           :{'tex':'#eta(W)',               'binning':[30,-3,3]},
     "mva_WH_dPhiMetLep"   :{'tex':'#Delta#phi(l,#nu) from W',   'binning':[30,0,pi]},
     "mva_WH_MT"           :{'tex':'M_{T}(W) (GeV)',        'binning':[300/20,0,300]},

     "mva_W_RatioPtlPtnu"  :{'tex':'p_{T}(l)/p_{T}(#nu)',   'binning':[30,0,3]},
     "mva_W_lepton_pt"     :{'tex':'p_{T}(l) from W',       'binning':[300/20,0,300]},
     "mva_W_lepton_eta"    :{'tex':'#eta(l) from W',        'binning':[30,-3,3]},
     "mva_WH_nu_pt"        :{'tex':'p_{T}(#nu) from W',     'binning':[300/20,0,300]},
     "mva_WH_nu_eta"       :{'tex':'#eta(l) from W',        'binning':[30,-3,3]},

     "mva_WH_Theta"        :{'tex':'#Theta',                'binning':[30,0,pi]},
     "mva_WH_theta"        :{'tex':'#theta',                'binning':[30,0,pi]},
     "mva_WH_phi"          :{'tex':'#phi',                  'binning':[30,-pi,pi]},

     'mva_WH_fLL'          :{ 'tex':'f_{LL}'         ,      'binning':[Nbins,0,1]},
     'mva_WH_f1TT'         :{ 'tex':'f_{1TT}'        ,      'binning':[Nbins,-1,1]},
     'mva_WH_f2TT'         :{ 'tex':'f_{2TT}'        ,      'binning':[Nbins, 0,4]},
     'mva_WH_f1LT'         :{ 'tex':'f_{1LT}'        ,      'binning':[Nbins,-1,1]},
     'mva_WH_f2LT'         :{ 'tex':'f_{2LT}'        ,      'binning':[Nbins,-1,1]},
     'mva_WH_f1tildeLT'    :{ 'tex':'#tilde{f}_{1LT}',      'binning':[Nbins,-1,1]},
     'mva_WH_f2tildeLT'    :{ 'tex':'#tilde{f}_{2LT}',      'binning':[Nbins,-1,1]},
     'mva_WH_fTTprime'     :{ 'tex':'f_{TT}'         ,      'binning':[Nbins,-1,1]},
     'mva_WH_ftildeTTprime':{ 'tex':'#tilde{f}_{TT}' ,      'binning':[Nbins,-1,1]},

     "mva_H_pt"            :{'tex':'p_{T}(H) (GeV)',        'binning':[1000/40,0,1000]},
     "mva_H_dijet_mass"    :{'tex':'M(H) (GeV)',            'binning':[30,90,150]},
     "mva_H_eta"           :{'tex':'#eta(H)',               'binning':[40,-4,4]},
     "mva_H_DeltaPhib1b2"  :{'tex':'#Delta#Phi(b_{1},b_{2}) from H',            'binning':[30,0,pi]},
     "mva_H_RatioPtb1b2"   :{'tex':'p_{T}(b_{2})/p_{T}(b_{2})',                 'binning':[30,0,1]},
     "mva_H_maxPtb1b2"     :{'tex':'max(p_{T}(b_{1}), p_{T}(b_{2})) from H',    'binning':[800/20,0,800]},
     "mva_H_minPtb1b2"     :{'tex':'min(p_{T}(b_{1}), p_{T}(b_{2})) from H',    'binning':[800/20,0,800]},
     "mva_H_DeltaRb1b2"    :{'tex':'#DeltaR(b_{1}, b_{2}) from H',              'binning':[30,0,6]},
     "mva_WH_deltaPhi"     :{'tex':'#Delta#Phi(W,H)',                           'binning':[30,0,pi]},
     "mva_WH_deltaR"       :{'tex':'#Delta R(W,H)',                             'binning':[30,0,6]},
     "mva_mWH"             :{'tex':'M(W,H)',                                    'binning':[50,0,1000]},
     "mva_thrust"          :{'tex':'Thrust',                                    'binning':[Nbins,0.5,1]},
     "mva_thrust_min"      :{'tex':'Thrust minor',                              'binning':[Nbins,0,.8]},
     "mva_extraJet0_W_deltaR"   :{'tex':'#Delta R( j_{0, extra}, W)',           'binning':[30,0,6]},
     "mva_extraJet1_W_deltaR"   :{'tex':'#Delta R( j_{1, extra}, W)',           'binning':[30,0,6]},
     #"mva_extraJet2_W_deltaR"   :{'tex':'#Delta R( j_{2, extra}, W)',          'binning':[30,0,6]},
     "mva_extraJet0_H_deltaR"   :{'tex':'#Delta R( j_{0, extra}, H)',           'binning':[30,0,6]},
     "mva_extraJet1_H_deltaR"   :{'tex':'#Delta R( j_{1, extra}, H)',           'binning':[30,0,6]},
     #"mva_extraJet2_H_deltaR"   :{'tex':'#Delta R( j_{2, extra}, H)',          'binning':[30,0,6]},
}

mva_vector_variables    =   {
    "weight":  { "name":"weight", "func":compute_weight_derivatives, "vars":["derivatives/F"], "varnames":["derivatives"]} 
}

## Using all variables
mva_variables_ = all_mva_variables.keys()
mva_variables_.sort()
mva_variables  = [ (key, value) for key, value in all_mva_variables.iteritems() if key in mva_variables_ ]

all_mva_variables['lumiweight1fb'] = (lambda event, sample: event.lumiweight1fb)

import numpy as np
import operator

# make predictions to be used with keras.predict
def predict_inputs( event, sample):
    return np.array([getattr( event, mva_variable) for mva_variable, _ in mva_variables])

# training selection
from TMB.Tools.delphesCutInterpreter import cutInterpreter

selectionString = cutInterpreter.cutString( 'singlelep-WHJet-onH' )
# selectionString = cutInterpreter.cutString( 'trilepT-minDLmass12-onZ1-njet4p-btag1' )

bit_derivatives  = weight_derivative_combinations[1:] 

bit_cfg = { derivative : { 
            'n_trees': 200,
            'max_depth'     : 3,
            'learning_rate' : 0.25,
            'min_size'      : 50,
            'clip_score_quantile': None,
            'calibrated'    : False,} for derivative in bit_derivatives }

def load(directory = '/mnt/hephy/cms/$USER/BIT/models/default/WH/', bit_derivatives=bit_derivatives):
    import sys, os
    sys.path.insert(0,os.path.expandvars("$CMSSW_BASE/src/BIT"))
    from BoostedInformationTree import BoostedInformationTree
    bits = {} 
    for derivative in bit_derivatives:
        if derivative == tuple(): continue

        filename = os.path.expandvars(os.path.join(directory, "bit_derivative_%s"% ('_'.join(derivative))) + '.pkl')
        try:
            print ("Loading %s for %r"%( filename, derivative))
            bits[derivative] = BoostedInformationTree.load(filename) 
        except IOError:
            print ("Could not load %s for %r"%( filename, derivative))

    return bits

#BITs -> Add pre-calculated BITs for limits, comparisons, etc.
if True:
    import sys, os, time
    sys.path.insert(0,os.path.expandvars("$CMSSW_BASE/src/BIT"))
    from BoostedInformationTree import BoostedInformationTree

    save_derivatives = [ ('cHW',), ('cHW','cHW'), ('cHWtil',), ('cHWtil','cHWtil'), ('cHW','cHWtil'), ('cHj3',), ('cHj3','cHj3'), ('cHj3', 'cHW'), ('cHj3', 'cHWtil')]
    bits             = {}

    directory = "/mnt/hephy/cms/robert.schoefbeck/BIT/models/"
    save_cfgs = [
        ( "nom",  "WH_delphes/first_try"),
        ( "bkgs", "WH_delphes_bkgs/first_try"),
        ]

    for name, subdir in save_cfgs:
        bits[name] = load(os.path.join( directory, subdir), bit_derivatives=save_derivatives)
        for der in save_derivatives:
            der_name = "BIT_"+name+"_"+"_".join(list(der))
            all_mva_variables[der_name] = (lambda event, sample, der_name=der_name: getattr( event, der_name))
            #all_mva_variables[der_name] = (lambda event, sample, der_name=der_name: bits[name][der].predict( event ))

    def bit_save_predict( event, sample ):

        for name, func in mva_variables:
            setattr( event, name, func(event, sample) )

        features = [ getattr( event, var[0] ) for var in mva_variables ]

        for name, _ in save_cfgs:
            for der in save_derivatives:
                prediction = bits[name][der].predict( features )
                setattr( event, "BIT_"+name+"_"+"_".join(der), prediction )

    sequence.append( bit_save_predict )

