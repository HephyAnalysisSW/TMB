#!/usr/bin/env python

# Standard imports
from operator                   import attrgetter
from math                       import pi, sqrt, cosh, cos, acos
import ROOT, os, copy

# RootTools
from RootTools.core.standard     import *

# helpers
from Analysis.Tools.helpers              import deltaPhi, deltaR, getObjDict
from TMB.Tools.helpers          import getCollection

from Analysis.Tools.WeightInfo       import WeightInfo

import logging
logger = logging.getLogger(__name__)

from Analysis.Tools.leptonJetArbitration     import cleanJetsAndLeptons

jetVars          = ['pt/F', 'eta/F', 'phi/F', 'bTag/F', 'bTagPhys/I']

jetVarNames      = [x.split('/')[0] for x in jetVars]

lepVars          = ['pt/F','eta/F','phi/F','pdgId/I','isolationVar/F', 'isolationVarRhoCorr/F']
lepVarNames      = [x.split('/')[0] for x in lepVars]

# Training variables
read_variables = [\
    "nBTag/I", "nBTag_loose/I",
    "recoMet_pt/F", "recoMet_phi/F",
    "recoLep[pt/F,eta/F,phi/F]",
    "recoZ_pt/F", "recoZ_eta/F", "recoZ_phi/F", "recoZ_mass/F", "recoZ_cosThetaStar/F", "recoZ_lldPhi/F", "recoZ_lldR/F",
    "recoNonZ_l1_index/I",
    "nrecoJet/I",
    "recoJet[%s]"%(",".join(jetVars)),
    "nrecoLep/I",
    "recoLep[%s]"%(",".join(lepVars)),
    VectorTreeVariable.fromString( "p[C/F]", nMax=200 ),
    "lumiweight1fb/F",
]

sequence = []

def getBJetindex( event ):
    ''' return highest pt bjet'''
    for i in range(event.nrecoJet):
        if event.recoJet_bTag[i]>=1:
            return i
    return -1

def getWpt( event, sample=None):
    # get the lepton and met
    lepton  = ROOT.TLorentzVector()
    met     = ROOT.TLorentzVector()
    lepton.SetPtEtaPhiM(event.recoLep_pt[event.recoNonZ_l1_index], event.recoLep_eta[event.recoNonZ_l1_index], event.recoLep_phi[event.recoNonZ_l1_index], 0)
    met.SetPtEtaPhiM(event.recoMet_pt, 0, event.recoMet_phi, 0)
    # get the W boson candidate
    W   = lepton + met
    event.W_pt = W.Pt()
sequence.append( getWpt )

def getM3l( event, sample=None):
    # get the invariant mass of the 3l system
    l = []
    for i in range(3):
        l.append(ROOT.TLorentzVector())
        l[i].SetPtEtaPhiM(event.recoLep_pt[i], event.recoLep_eta[i], event.recoLep_phi[i],0)
    event.m3l = (l[0] + l[1] + l[2]).M()
sequence.append( getM3l )

from TMB.Tools.objectSelection import isBJet
def make_jets( event, sample ):
    event.jets     = [getObjDict(event, 'recoJet_', jetVarNames, i) for i in range(int(event.nrecoJet))]
    event.bJets    = filter(lambda j:j['bTag']>=1 and abs(j['eta'])<=2.4    , event.jets)
sequence.append( make_jets )

def getAngles(event, sample=None):
    event.nonZ1_l1_Z1_deltaPhi = deltaPhi(event.recoLep_phi[event.recoNonZ_l1_index], event.recoZ_phi)
    event.Z1_j1_deltaPhi       = deltaPhi(event.recoZ_phi, event.recoJet_phi[0])
    event.nonZ1_l1_Z1_deltaEta = abs(event.recoLep_eta[event.recoNonZ_l1_index] - event.recoZ_eta)
    event.nonZ1_l1_Z1_deltaR   = deltaR({'eta':event.recoLep_eta[event.recoNonZ_l1_index], 'phi':event.recoLep_phi[event.recoNonZ_l1_index]}, {'eta':event.recoZ_eta, 'phi':event.recoZ_phi})
    event.jet0_Z1_deltaR       = deltaR({'eta':event.recoJet_eta[0], 'phi':event.recoJet_phi[0]}, {'eta':event.recoZ_eta, 'phi':event.recoZ_phi})
    event.jet0_nonZ1_l1_deltaR = deltaR({'eta':event.recoJet_eta[0], 'phi':event.recoJet_phi[0]}, {'eta':event.recoLep_eta[event.recoNonZ_l1_index], 'phi':event.recoLep_phi[event.recoNonZ_l1_index]})
    event.jet1_Z1_deltaR       = deltaR({'eta':event.recoJet_eta[1], 'phi':event.recoJet_phi[1]}, {'eta':event.recoZ_eta, 'phi':event.recoZ_phi})
    event.jet1_nonZ1_l1_deltaR = deltaR({'eta':event.recoJet_eta[1], 'phi':event.recoJet_phi[1]}, {'eta':event.recoLep_eta[event.recoNonZ_l1_index], 'phi':event.recoLep_phi[event.recoNonZ_l1_index]})
    event.jet2_Z1_deltaR       = deltaR({'eta':event.recoJet_eta[2], 'phi':event.recoJet_phi[2]}, {'eta':event.recoZ_eta, 'phi':event.recoZ_phi})
    event.jet2_nonZ1_l1_deltaR = deltaR({'eta':event.recoJet_eta[2], 'phi':event.recoJet_phi[2]}, {'eta':event.recoLep_eta[event.recoNonZ_l1_index], 'phi':event.recoLep_phi[event.recoNonZ_l1_index]})
    i_bjet = getBJetindex(event)
    if i_bjet>=0:
        event.bJet_Z1_deltaR      = deltaR({'eta':event.recoJet_eta[i_bjet], 'phi':event.recoJet_phi[i_bjet]}, {'eta':event.recoZ_eta, 'phi':event.recoZ_phi})
        event.bJet_nonZ1l1_deltaR = deltaR({'eta':event.recoJet_eta[i_bjet], 'phi':event.recoJet_phi[i_bjet]}, {'eta':event.recoLep_eta[event.recoNonZ_l1_index], 'phi':event.recoLep_phi[event.recoNonZ_l1_index]})
    else:
        event.bJet_Z1_deltaR      = -1 
        event.bJet_nonZ1l1_deltaR = -1
    
sequence.append( getAngles )

weight_variables = ['cHq1Re11', 'cHq1Re22', 'cHq1Re33', 'cHq3Re11', 'cHq3Re22', 'cHq3Re33', 'cHuRe11', 'cHuRe22', 'cHuRe33', 'cHdRe11', 'cHdRe22', 'cHdRe33', 'cHudRe11', 'cHudRe22', 'cHudRe33']
max_order        = 2

from TMB.Samples.pp_gen_v8 import *
training_samples = [ ttZ01j, WZTo3L1Nu]

assert len(training_samples)==len(set([s.name for s in training_samples])), "training_samples names are not unique!"

reweight_pkl     = ttZ01j.reweight_pkl
weightInfo       = WeightInfo( reweight_pkl )
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
     "mva_m3l"                   :(lambda event, sample: event.m3l),
     "mva_nrecoJet"              :(lambda event, sample: event.nrecoJet),
     "mva_nBTag"                 :(lambda event, sample: event.nBTag),

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

# Z1 kinematics
     "mva_Z1_pt"                 :(lambda event, sample: event.recoZ_pt),
     "mva_Z1_eta"                :(lambda event, sample: event.recoZ_eta),
     "mva_Z1_cosThetaStar"       :(lambda event, sample: event.recoZ_cosThetaStar),

# extra lepton kinematics
     "mva_lnonZ1_pt"             :(lambda event, sample: event.recoLep_pt[event.recoNonZ_l1_index]),
     "mva_lnonZ1_eta"            :(lambda event, sample: event.recoLep_eta[event.recoNonZ_l1_index]),

# leptonic W
     "mva_W_pt"                  :(lambda event, sample: event.W_pt),

# Z1 vs. other objects
     "mva_nonZ1_l1_Z1_deltaPhi"  :(lambda event, sample: event.nonZ1_l1_Z1_deltaPhi     if event.nrecoLep >= 1 else -1 ),
     "mva_nonZ1_l1_Z1_deltaR"    :(lambda event, sample: event.nonZ1_l1_Z1_deltaR),

     "mva_jet0_Z1_deltaR"        :(lambda event, sample: event.jet0_Z1_deltaR         if event.nrecoJet >=1 else -1),
     "mva_jet1_Z1_deltaR"        :(lambda event, sample: event.jet1_Z1_deltaR         if event.nrecoJet >=2 else -1),
     "mva_jet2_Z1_deltaR"        :(lambda event, sample: event.jet2_Z1_deltaR         if event.nrecoJet >=3 else -1),

# nonZ1_l1 vs. other objects
     "mva_jet0_nonZl1_deltaR"    :(lambda event, sample: event.jet0_nonZ1_l1_deltaR    if event.nrecoJet >=1 else -1),
     "mva_jet1_nonZl1_deltaR"    :(lambda event, sample: event.jet1_nonZ1_l1_deltaR    if event.nrecoJet >=2 else -1),
     "mva_bJet_Z1_deltaR"        :(lambda event, sample: event.bJet_Z1_deltaR),
     "mva_bJet_non_Z1l1_deltaR"  :(lambda event, sample: event.bJet_nonZ1l1_deltaR),
}

#all_mva_variables.update( {  "locsc_"+weight['name']: (lambda event, sample, name="locsc_"+weight['name']: getattr( event, name) ) for weight in weights[1:] } )

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
selectionString = cutInterpreter.cutString( 'trilep-onZ-njet2p' )
# selectionString = cutInterpreter.cutString( 'trilepT-minDLmass12-onZ1-njet4p-btag1' )

bit_cfg = { 'n_trees': 100,
            'max_depth'     : 3,
            'learning_rate' : 0.20,
            'min_size'      : 50,
            'calibrated'    : False,
            'global_score_subtraction': True,
    }

bit_derivatives  = [ ('cHq1Re11',), ('cHq1Re22',), ('cHq1Re33',), ('cHq1Re11','cHq1Re11'), ('cHq1Re22','cHq1Re22'), ('cHq1Re33','cHq1Re33')]
bit_derivatives += [ ('cHq3Re11',), ('cHq3Re22',), ('cHq3Re33',), ('cHq3Re11','cHq3Re11'), ('cHq3Re22','cHq3Re22'), ('cHq3Re33','cHq3Re33')]
bit_derivatives += [ ('cHuRe11',), ('cHuRe22',), ('cHuRe33',), ('cHuRe11','cHuRe11'), ('cHuRe22','cHuRe22'), ('cHuRe33','cHuRe33')]
bit_derivatives += [ ('cHdRe11',), ('cHdRe22',), ('cHdRe33',), ('cHdRe11','cHdRe11'), ('cHdRe22','cHdRe22'), ('cHdRe33','cHdRe33')]
bit_derivatives += [ ('cHudRe11',), ('cHudRe22',), ('cHudRe33',), ('cHudRe11','cHudRe11'), ('cHudRe22','cHudRe22'), ('cHudRe33','cHudRe33')]

def load(directory = '/groups/hephy/cms/$USER/BIT/models/default/ttZ_3l_flavor_delphes/', bit_derivatives=bit_derivatives):
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
