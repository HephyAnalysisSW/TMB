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
    "recoZ_pt/F", "recoZ_eta/F", "recoZ_phi/F", "recoZ_mass/F", "recoZ_cosThetaStar/F", "recoZ_lldPhi/F", "recoZ_lldR/F", "recoZ_l1_index/I", "recoZ_l2_index/I",
    "nrecoJet/I",
    "recoJet[%s]"%(",".join(jetVars)),
    "nrecoLep/I",
    "recoLep[%s]"%(",".join(lepVars)),
    "ngenLep/I", "genLep[pt/F,eta/F,phi/F,pdgId/I,mother_pdgId/I]",
    "lumiweight1fb/F",
    "genW[pt/F,eta/F,phi/F,l1_index/I,l2_index/I]", "ngenW/I",
    "evt/l", "run/I", "lumi/I", "np/I", "nweight/I",
    "H_dijet_mass/F", "H_pt/F", "H_j1_index/I", "H_j2_index/I", 
    #"WH_W_pt/F", "WH_dPhiMetLep/F", "WH_MT/F", "WH_nu_pt/F", "WH_nu_eta/F", "WH_nu_phi/F", "WH_nu_E", "WH_Theta/F", "WH_theta/F", "WH_phi/F", 
    "H_j1_index/I", "H_j2_index/I", 
    "ZH_Theta/F", "ZH_theta/F", "ZH_phi/F", 

]
read_variables += [VectorTreeVariable.fromString( "p[C/F]", nMax=200 )]

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
    event.lepton1 = all_leps[event.recoZ_l1_index]
    event.lepton2 = all_leps[event.recoZ_l2_index]

sequence.append( makeLeptonic )

def makeAngles(event, sample):

    event.sin2thetaSin2Theta = sin(2*event.ZH_theta)*sin(2*event.ZH_Theta)

    event.fLL         = sin(event.ZH_Theta)**2*sin(event.ZH_theta)**2
    event.f1TT        = cos(event.ZH_theta)*cos(event.ZH_Theta)
    event.f2TT        = (1+cos(event.ZH_theta)**2)*(1+cos(event.ZH_Theta)**2)
    event.f1LT        = cos(event.ZH_phi)*sin(event.ZH_theta)*sin(event.ZH_Theta)
    event.f2LT        = event.f1LT*cos(event.ZH_theta)*cos(event.ZH_Theta)
    event.f1tildeLT   = sin(event.ZH_phi)*sin(event.ZH_theta)*sin(event.ZH_Theta)
    event.f2tildeLT   = event.f1tildeLT*cos(event.ZH_theta)*cos(event.ZH_Theta)
    event.fTTprime    = cos(2*event.ZH_phi)*event.fLL
    event.ftildeTTprime=sin(2*event.ZH_phi)*event.fLL

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
    event.Z_vecP4 = ROOT.TLorentzVector( event.recoZ_pt*cos(event.recoZ_phi), event.recoZ_pt*sin(event.recoZ_phi),  event.recoZ_pt*sinh(event.recoZ_eta), event.recoZ_pt*cosh(event.recoZ_eta) )
    (event.thrust, event.thrust_min) = Thrust(event.Z_vecP4,[j['vecP4'] for j in event.jets] )

sequence.append( makeThrust )


weight_variables = ['cHj3', 'cHW', 'cHWtil']
max_order        = 2

from TMB.Samples.pp_gen_v10 import *
training_samples = [ ZH ]

assert len(training_samples)==len(set([s.name for s in training_samples])), "training_samples names are not unique!"

weightInfo       = WeightInfo(ZH.reweight_pkl )
weightInfo.set_order(2)

weight_derivatives = []
weight_derivative_combinations = []
for i_comb, comb in enumerate(weightInfo.make_combinations(weight_variables, max_order)):
    #print name, i_comb, comb, weightInfo.get_diff_weight_string(comb)
    weight = {}
    weight['string'] = weightInfo.get_diff_weight_string(comb)
    weight['func']   = weightInfo.get_diff_weight_func(comb)
    weight['name']   = '_'.join(comb)
    weight['comb']   = comb
    weight_derivatives.append( weight )
    weight_derivative_combinations.append(comb)

def compute_weight_derivatives( event, sample ):
    vector = [{'derivatives':weight['func'](event, sample)} for weight in weight_derivatives]

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

# Z kinematics
     "mva_Z_pt"                  :(lambda event, sample: event.recoZ_pt),
     "mva_Z_eta"                 :(lambda event, sample: event.recoZ_eta),
     "mva_Z_cosThetaStar"        :(lambda event, sample: event.recoZ_cosThetaStar),
     "mva_Z_lldPhi"              :(lambda event, sample: event.recoZ_lldPhi), 
     "mva_Z_lldR"                :(lambda event, sample: event.recoZ_lldR),
     "mva_Z_RatioPtl1l2"         :(lambda event, sample: event.lepton2['pt']/event.lepton1['pt']) ,
     "mva_Z_maxPtl1l2"           :(lambda event, sample: max([event.lepton1['pt'], event.lepton2['pt']])) ,
     "mva_Z_minPtl1l2"           :(lambda event, sample: min([event.lepton1['pt'], event.lepton2['pt']])) ,

# Z decay angles & signed
     "mva_ZH_Theta"              :(lambda event, sample: event.ZH_Theta),
     "mva_ZH_theta"              :(lambda event, sample: event.ZH_theta),
     "mva_ZH_phi"                :(lambda event, sample: event.ZH_phi),
#     "mva_ZH_signed_Theta"       :(lambda event, sample: copysign(event.ZH_Theta, event.sin2thetaSin2Theta) ),
#     "mva_ZH_signed_theta"       :(lambda event, sample: copysign(event.ZH_theta, event.sin2thetaSin2Theta) ),
#     "mva_ZH_signed_phi"         :(lambda event, sample: copysign(event.ZH_phi, event.sin2thetaSin2Theta) ),
    "mva_ZH_fLL"                 :(lambda event, sample: event.fLL),
    "mva_ZH_f1TT"                :(lambda event, sample: event.f1TT),
    "mva_ZH_f2TT"                :(lambda event, sample: event.f2TT),
    "mva_ZH_f1LT"                :(lambda event, sample: event.f1LT),
    "mva_ZH_f2LT"                :(lambda event, sample: event.f2LT),
    "mva_ZH_f1tildeLT"           :(lambda event, sample: event.f1tildeLT),
    "mva_ZH_f2tildeLT"           :(lambda event, sample: event.f2tildeLT),
    "mva_ZH_fTTprime"            :(lambda event, sample: event.fTTprime),
    "mva_ZH_ftildeTTprime"       :(lambda event, sample: event.ftildeTTprime),

# H kinematics
     "mva_H_pt"                  :(lambda event, sample: event.H_pt),
     "mva_H_dijet_mass"          :(lambda event, sample: event.H_dijet_mass),
     "mva_H_eta"                 :(lambda event, sample: event.H_eta ),
     "mva_H_DeltaPhib1b2"        :(lambda event, sample: deltaPhi(event.H_j1['phi'], event.H_j2['phi'])) ,
     "mva_H_RatioPtb1b2"         :(lambda event, sample: event.H_j2['pt']/event.H_j1['pt']) ,
     "mva_H_maxPtb1b2"           :(lambda event, sample: max([event.H_j1['pt'], event.H_j2['pt']])) ,
     "mva_H_minPtb1b2"           :(lambda event, sample: min([event.H_j1['pt'], event.H_j2['pt']])) ,
     "mva_H_DeltaRb1b2"          :(lambda event, sample: deltaR(event.H_j1, event.H_j2)) ,

# Z vs. H
     "mva_ZH_deltaPhi"           :(lambda event, sample: deltaPhi(event.H_phi, event.recoZ_phi)),
     "mva_ZH_deltaR"             :(lambda event, sample: deltaR({'phi':event.H_phi, 'eta':event.H_eta},{'phi': event.recoZ_phi, 'eta':event.recoZ_eta})),
     "mva_mZH"                   :(lambda event, sample: (event.Z_vecP4+event.H_vecP4).M()),

# thrust
    "mva_thrust"                 :(lambda event, sample: event.thrust),
    "mva_thrust_min"             :(lambda event, sample: event.thrust_min),

# Z vs. other objects
     "mva_extraJet0_Z_deltaR"    :(lambda event, sample: deltaR({'phi':event.extraJets[0]['phi'], 'eta':event.extraJets[0]['eta']}, {'phi':event.recoZ_phi, 'eta':event.recoZ_eta} )  if event.nextraJet >=1 else -1),
     "mva_extraJet1_Z_deltaR"    :(lambda event, sample: deltaR({'phi':event.extraJets[1]['phi'], 'eta':event.extraJets[1]['eta']}, {'phi':event.recoZ_phi, 'eta':event.recoZ_eta} )  if event.nextraJet >=2 else -1),
     #"mva_extraJet2_Z_deltaR"    :(lambda event, sample: deltaR({'phi':event.extraJets[2]['phi'], 'eta':event.extraJets[2]['eta']}, {'phi':event.recoZ_phi, 'eta':event.recoZ_eta} )  if event.nextraJet >=3 else -1),

# H vs. other objects
     "mva_extraJet0_H_deltaR"    :(lambda event, sample: deltaR({'phi':event.extraJets[0]['phi'], 'eta':event.extraJets[0]['eta']}, {'phi':event.H_phi, 'eta':event.H_eta} )  if event.nextraJet >=1 else -1),
     "mva_extraJet1_H_deltaR"    :(lambda event, sample: deltaR({'phi':event.extraJets[1]['phi'], 'eta':event.extraJets[1]['eta']}, {'phi':event.H_phi, 'eta':event.H_eta} )  if event.nextraJet >=2 else -1),
     #"mva_extraJet2_H_deltaR"    :(lambda event, sample: deltaR({'phi':event.extraJets[2]['phi'], 'eta':event.extraJets[2]['eta']}, {'phi':event.H_phi, 'eta':event.H_eta} )  if event.nextraJet >=3 else -1),
}
# global event properties

Nbins = 30
plot_options = {
     "mva_ht"              :{'tex':'H_{T} (GeV)'},
     "mva_met_pt"          :{'tex':'p_{T}^{miss}'},
     "mva_nextraJet"       :{'tex':'N_{jet, extra}',   'binning':[8,0,8]},
     "mva_nBTag"           :{'tex':'N_{b-tag}', 'binning':[4,0,4]},
     "mva_jet0_pt"         :{'tex':'p_{T}(j_{0}) (GeV)'},
     "mva_jet0_eta"        :{'tex':'#eta (j_{0})'},
     "mva_jet0_btag"       :{'tex':'b-jet disc. of j_{0}'},
     "mva_jet1_pt"         :{'tex':'p_{T}(j_{1}) (GeV)'},
     "mva_jet1_eta"        :{'tex':'#eta (j_{1})'},
     "mva_jet1_btag"       :{'tex':'b-jet disc. of j_{1}'},
     "mva_jet2_pt"         :{'tex':'p_{T}(j_{2}) (GeV)'},
     "mva_jet2_eta"        :{'tex':'#eta (j_{2})'},
     "mva_jet2_btag"       :{'tex':'b-jet disc. of j_{2}'},
     "mva_Z_pt"            :{'tex':'p_{T}(Z) (GeV)'},
     "mva_Z_eta"           :{'tex':'#eta(Z)'},
     "mva_Z_cosThetaStar"  :{'tex':'cos(#theta^{*})'},
     "mva_Z_lldPhi"        :{'tex':'#Delta#phi(ll) from Z'},
     "mva_Z_lldR"          :{'tex':'#Delta R(ll) from Z'},

     "mva_Z_RatioPtl1l2"   :{'tex':'p_{T}(l_{2})/p_{T}(l_{2}) from Z'},
     "mva_Z_maxPtl1l2"     :{'tex':'max(p_{T}(l_{1}), p_{T}(l_{2})) from Z'},
     "mva_Z_minPtl1l2"     :{'tex':'min(p_{T}(l_{1}), p_{T}(l_{2})) from Z'},

     "mva_ZH_Theta"        :{'tex':'#Theta'},
     "mva_ZH_theta"        :{'tex':'#theta'},
     "mva_ZH_phi"          :{'tex':'#phi'},
     #"mva_ZH_signed_Theta" :{'tex':'#Theta#times sign(sin(2#Theta) sin(2#theta)'},
     #"mva_ZH_signed_theta" :{'tex':'#theta#times sign(sin(2#Theta) sin(2#theta)'},
     #"mva_ZH_signed_phi"   :{'tex':'#phi#times sign(sin(2#Theta) sin(2#theta)'},

     'mva_ZH_fLL'                  : {'binning':[Nbins,0,1],        'tex':'f_{LL}'          ,},
     'mva_ZH_f1TT'                 : {'binning':[Nbins,-1,1],        'tex':'f_{1TT}'         ,},
     'mva_ZH_f2TT'                 : {'binning':[Nbins, 0,4],        'tex':'f_{2TT}'         ,},
     'mva_ZH_f1LT'                 : {'binning':[Nbins,-1,1],        'tex':'f_{1LT}'         ,},
     'mva_ZH_f2LT'                 : {'binning':[Nbins,-1,1],        'tex':'f_{2LT}'         ,},
     'mva_ZH_f1tildeLT'            : {'binning':[Nbins,-1,1],        'tex':'#tilde{f}_{1LT}' ,},
     'mva_ZH_f2tildeLT'            : {'binning':[Nbins,-1,1],        'tex':'#tilde{f}_{2LT}' ,},
     'mva_ZH_fTTprime'             : {'binning':[Nbins,-1,1],        'tex':'f_{TT}'     ,},
     'mva_ZH_ftildeTTprime'        : {'binning':[Nbins,-1,1],        'tex':'#tilde{f}_{TT}',},

     "mva_H_pt"            :{'tex':'p_{T}(H) (GeV)'},
     "mva_H_dijet_mass"    :{'tex':'M(H) (GeV)'},
     "mva_H_absEta"        :{'tex':'|#eta|(H)'},
     "mva_H_DeltaPhib1b2"  :{'tex':'#Delta#Phi(b_{1},b_{2}) from H'},
     "mva_H_RatioPtb1b2"   :{'tex':'p_{T}(b_{2})/p_{T}(b_{2}) from H'},
     "mva_H_maxPtb1b2"     :{'tex':'max(p_{T}(b_{1}), p_{T}(b_{2})) from H'},
     "mva_H_minPtb1b2"     :{'tex':'min(p_{T}(b_{1}), p_{T}(b_{2})) from H'},
     "mva_H_DeltaRb1b2"    :{'tex':'#DeltaR(b_{1}, b_{2}) from H'},
     "mva_ZH_deltaPhi"     :{'tex':'#Delta#Phi(Z,H)'},
     "mva_ZH_deltaR"       :{'tex':'#Delta R(Z,H)'},
     "mva_mZH"             :{'tex':'M(Z,H)'},
     "mva_thrust"          :{'tex':'Thrust'},
     "mva_thrust_min"      :{'tex':'Thrust minor'},
     "mva_extraJet0_Z_deltaR"   :{'tex':'#Delta R( j_{0, extra}, Z)'},
     "mva_extraJet1_Z_deltaR"   :{'tex':'#Delta R( j_{1, extra}, Z)'},
     #"mva_extraJet2_Z_deltaR"   :{'tex':'#Delta R( j_{2, extra}, Z)'},
     "mva_extraJet0_H_deltaR"   :{'tex':'#Delta R( j_{0, extra}, H)'},
     "mva_extraJet1_H_deltaR"   :{'tex':'#Delta R( j_{1, extra}, H)'},
     #"mva_extraJet2_H_deltaR"   :{'tex':'#Delta R( j_{2, extra}, H)'},
}

mva_vector_variables    =   {
    "weight":  { "name":"weight", "func":compute_weight_derivatives, "vars":["derivatives/F"], "varnames":["derivatives"], 'nMax':len(weight_derivatives)} 
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

selectionString = cutInterpreter.cutString( 'dilep-ZHJet-onZ-onH' )
# selectionString = cutInterpreter.cutString( 'trilepT-minDLmass12-onZ1-njet4p-btag1' )

bit_derivatives  = weight_derivative_combinations[1:] 

bit_cfg = { derivative : { 
            'n_trees': 100,
            'max_depth'     : 3,
            'learning_rate' : 0.25,
            'min_size'      : 50,
            'clip_score_quantile': None,
            'calibrated'    : False,} for derivative in bit_derivatives }

def load(directory = '/mnt/hephy/cms/$USER/BIT/models/default/ZH/', bit_derivatives=bit_derivatives):
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
