#!/usr/bin/env python

# Standard imports
from operator                   import attrgetter
from math                       import pi, sqrt, cosh, cos, acos
import ROOT, os
import copy
# RootTools
from RootTools.core.standard     import *

# helpers
from tttt.Tools.helpers          import deltaPhi, deltaR2, deltaR, getCollection, getObjDict
#from tWZ.Tools.objectSelection  import isBJet, isAnalysisJet
from Analysis.Tools.WeightInfo       import WeightInfo

import logging
logger = logging.getLogger(__name__)

from Analysis.Tools.leptonJetArbitration     import cleanJetsAndLeptons

jetVars          = ['pt/F', 'eta/F', 'phi/F', 'btagDeepB/F', 'btagDeepFlavB/F', 'index/I']

jetVarNames      = [x.split('/')[0] for x in jetVars]

lstm_jets_maxN   = 10
lstm_jetVars     = ['pt/F', 'eta/F', 'phi/F', 'btagDeepFlavB/F', 'btagDeepFlavCvB/F', 'btagDeepFlavCvL/F', 'chEmEF/F', 'chHEF/F', 'neEmEF/F', 'neHEF/F', 'muEF/F', 'puId/F', 'qgl/F']
lstm_jetVarNames = [x.split('/')[0] for x in lstm_jetVars]


lepVars          = ['pt/F','eta/F','phi/F','pdgId/I','cutBased/I','miniPFRelIso_all/F','pfRelIso03_all/F','mvaFall17V2Iso_WP90/O', 'mvaTOP/F', 'sip3d/F','lostHits/I','convVeto/I','dxy/F','dz/F','charge/I','deltaEtaSC/F','mediumId/I','eleIndex/I','muIndex/I']
lepVarNames      = [x.split('/')[0] for x in lepVars]

# Training variables
read_variables = [\
                    "nBTag/I",
                    "nJetGood/I",
                    "nlep/I",
                    "m3/F",
                    "JetGood[%s]"%(",".join(jetVars)),
                    "Jet[%s]"%(",".join(lstm_jetVars)),
                    "lep[%s]"%(",".join(lepVars)),
                    "met_pt/F", "met_phi/F",
                    "l1_pt/F",
                    "l1_eta/F",
                    "l1_phi/F",
                    "l2_pt/F",
                    "l2_eta/F",
                    "l2_phi/F",
                    "year/I",
                    ]
# sequence 
sequence = []

# Fisher informations
FIs = {
}

from tttt.Tools.objectSelection import isBJet
def make_jets( event, sample ):
    event.jets     = [getObjDict(event, 'JetGood_', jetVarNames, i) for i in range(int(event.nJetGood))] 
    event.bJets    = filter(lambda j:isBJet(j, year=event.year) and abs(j['eta'])<=2.4    , event.jets)
sequence.append( make_jets )

all_mva_variables = {

# global event properties     
     "mva_nJetGood"              :(lambda event, sample: event.nJetGood),
     "mva_nBTag"                 :(lambda event, sample: event.nBTag),
     "mva_nlep"                  :(lambda event, sample: event.nlep),

     "mva_mT_l1"                 :(lambda event, sample: sqrt(2*event.l1_pt*event.met_pt*(1-cos(event.l1_phi-event.met_phi)))),
     "mva_mT_l2"                 :(lambda event, sample: sqrt(2*event.l2_pt*event.met_pt*(1-cos(event.l2_phi-event.met_phi)))),
     "mva_ml_12"                 :(lambda event, sample: sqrt(2*event.l1_pt*event.l2_pt*(cosh(event.l1_eta-event.l2_eta)-cos(event.l1_phi-event.l2_phi)))),
     "mva_met_pt"                :(lambda event, sample: event.met_pt),
     "mva_l1_pt"                 :(lambda event, sample: event.l1_pt),
     "mva_l1_eta"                :(lambda event, sample: event.l1_eta),
     "mva_l2_pt"                 :(lambda event, sample: event.l2_pt),
     "mva_l2_eta"                :(lambda event, sample: event.l2_eta),
     
     "mva_mj_12"                 :(lambda event, sample: sqrt(event.jets[0]['pt']*event.jets[1]['pt']*cosh(event.jets[0]['eta']-event.jets[1]['eta'])-cos(event.jets[0]['phi']-event.jets[1]['phi']))  if event.nJetGood >=2 else 0),
     "mva_mlj_11"                :(lambda event, sample: sqrt(event.l1_pt*event.jets[0]['pt']*cosh(event.l1_eta-event.jets[0]['eta'])-cos(event.l1_phi-event.jets[0]['phi'])) if event.nJetGood >=1 else 0),
     "mva_mlj_12"                :(lambda event, sample: sqrt(event.l1_pt*event.jets[1]['pt']*cosh(event.l1_eta-event.jets[1]['eta'])-cos(event.l1_phi-event.jets[1]['phi'])) if event.nJetGood >=2 else 0),

     "mva_dPhil_12"              :(lambda event, sample: acos(cos(event.l1_phi-event.l2_phi))),
     "mva_dPhij_12"              :(lambda event, sample: acos(cos(event.JetGood_phi[0]-event.JetGood_phi[1])) if event.nJetGood >=2 else 0),

     "mva_dEtal_12"              :(lambda event, sample: event.l1_eta-event.l2_eta),
     "mva_dEtaj_12"              :(lambda event, sample: event.JetGood_eta[0] - event.JetGood_eta[1] if event.nJetGood >=2 else -10),

     "mva_ht"                    :(lambda event, sample: sum( [j['pt'] for j in event.jets] ) ),
     "mva_htb"                   :(lambda event, sample: sum( [j['pt'] for j in event.bJets] ) ),
     "mva_ht_ratio"              :(lambda event, sample: sum( [j['pt'] for j in event.jets[:4]])/ sum( [j['pt'] for j in event.jets ]) if event.nJetGood>=4 else 1 ),

     "mva_jet0_pt"               :(lambda event, sample: event.JetGood_pt[0]          if event.nJetGood >=1 else 0),
     "mva_jet0_eta"              :(lambda event, sample: event.JetGood_eta[0]         if event.nJetGood >=1 else -10),
     "mva_jet0_btagDeepFlavB"        :(lambda event, sample: event.JetGood_btagDeepFlavB[0]   if (event.nJetGood >=1 and event.JetGood_btagDeepFlavB[0]>-10) else -10),
     "mva_jet1_pt"               :(lambda event, sample: event.JetGood_pt[1]          if event.nJetGood >=2 else 0),
     "mva_jet1_eta"              :(lambda event, sample: event.JetGood_eta[1]         if event.nJetGood >=2 else -10),
     "mva_jet1_btagDeepFlavB"        :(lambda event, sample: event.JetGood_btagDeepFlavB[1]   if (event.nJetGood >=2 and event.JetGood_btagDeepFlavB[1]>-10) else -10),
     "mva_jet2_pt"               :(lambda event, sample: event.JetGood_pt[2]          if event.nJetGood >=3 else 0),
     "mva_jet2_eta"              :(lambda event, sample: event.JetGood_eta[2]         if event.nJetGood >=3 else -10),
     "mva_jet2_btagDeepFlavB"        :(lambda event, sample: event.JetGood_btagDeepFlavB[2]   if (event.nJetGood >=3 and event.JetGood_btagDeepFlavB[2]>-10) else -10),

     "mva_jet3_pt"               :(lambda event, sample: event.JetGood_pt[3]          if event.nJetGood >=4 else 0),
     "mva_jet4_pt"               :(lambda event, sample: event.JetGood_pt[4]          if event.nJetGood >=5 else 0),
     "mva_jet5_pt"               :(lambda event, sample: event.JetGood_pt[5]          if event.nJetGood >=6 else 0),
     "mva_jet6_pt"               :(lambda event, sample: event.JetGood_pt[6]          if event.nJetGood >=7 else 0),
     "mva_jet7_pt"               :(lambda event, sample: event.JetGood_pt[7]          if event.nJetGood >=8 else 0),
}

def lstm_jets(event, sample):
    jets = [ getObjDict( event, 'Jet_', lstm_jetVarNames, event.JetGood_index[i] ) for i in range(int(event.nJetGood)) ]
    #jets = filter( jet_vector_var['selector'], jets )
    return jets 

# for the filler
mva_vector_variables    =   {
    #"mva_Jet":  {"name":"Jet", "vars":lstm_jetVars, "varnames":lstm_jetVarNames, "selector": (lambda jet: True), 'maxN':10} 
    "mva_Jet":  {"func":lstm_jets, "name":"Jet", "vars":lstm_jetVars, "varnames":lstm_jetVarNames} 
}

## Using all variables
mva_variables_ = all_mva_variables.keys()
mva_variables_.sort()
mva_variables  = [ (key, value) for key, value in all_mva_variables.iteritems() if key in mva_variables_ ]

import numpy as np
import operator

def predict_inputs( event, sample, jet_lstm = False):
    flat_variables = np.array([[getattr( event, mva_variable) for mva_variable, _ in mva_variables]])
    if jet_lstm:
        lstm_jets_maxN = 10 #remove after retraining
        jet_vector_var = mva_vector_variables["mva_Jet"]
        jets = mva_vector_variables["mva_Jet"]["func"](event,sample=None) 
        jets =  [ [ operator.itemgetter(varname)(jet) for varname in lstm_jetVarNames] for jet in jets[:lstm_jets_maxN] ]
        # zero padding
        jets += [ [0.]*len(lstm_jetVarNames)]*(max(0, lstm_jets_maxN-len(jets))) 
        jets = np.array([jets])

        return [ flat_variables, jets ]
    else:
        return   flat_variables

#define training samples for multiclassification
#import TMB.Samples.nanoTuples_RunII_nanoAODv6_dilep_pp as samples
#import tttt.samples.nano_mc_private_UL20_RunII_postProcessed_dilep as samples

## ttbar gen classification: https://github.com/cms-top/cmssw/blob/topNanoV6_from-CMSSW_10_2_18/TopQuarkAnalysis/TopTools/plugins/GenTtbarCategorizer.cc
#TTLep_bb    = copy.deepcopy( samples.TTLep )
#TTLep_bb.name = "TTLep_bb"
#TTLep_bb.texName = samples.TTLep.name+" (b#overline{b})"
#TTLep_bb.setSelectionString( "genTtbarId%100>=50" )
#TTLep_cc    = copy.deepcopy( samples.TTLep )
#TTLep_cc.name = "TTLep_cc"
#TTLep_cc.texName = samples.TTLep.name+" (c#overline{c})"
#TTLep_cc.setSelectionString( "genTtbarId%100>=40&&genTtbarId%100<50" )
#TTLep_other = copy.deepcopy( samples.TTLep )
#TTLep_other.name = "TTLep_other"
#TTLep_other.texName = samples.TTLep.name+" (other)"
#TTLep_other.setSelectionString( "genTtbarId%100<40" )
#
#training_samples = [ samples.TTTT, TTLep_bb, TTLep_cc, TTLep_other ]
#
#assert len(training_samples)==len(set([s.name for s in training_samples])), "training_samples names are not unique!"
#
## training selection
#
#from tttt.Tools.cutInterpreter import cutInterpreter
#selectionString = cutInterpreter.cutString( 'dilepVL-ht400' )
