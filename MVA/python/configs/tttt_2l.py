#!/usr/bin/env python

# Standard imports
from operator                   import attrgetter
from math                       import pi, sqrt, cosh, cos, acos
import ROOT, os

# RootTools
from RootTools.core.standard     import *

# helpers
from TMB.Tools.helpers          import deltaPhi, deltaR2, deltaR, getCollection, getObjDict
#from tWZ.Tools.objectSelection  import isBJet, isAnalysisJet
from Analysis.Tools.WeightInfo       import WeightInfo

import logging
logger = logging.getLogger(__name__)

from Analysis.Tools.leptonJetArbitration     import cleanJetsAndLeptons

jetVars          = ['pt/F', 'eta/F', 'phi/F', 'btagDeepB/F']
jetVarNames      = [x.split('/')[0] for x in jetVars]

lepVars          = ['pt/F','eta/F','phi/F','pdgId/I','cutBased/I','miniPFRelIso_all/F','pfRelIso03_all/F','mvaFall17V2Iso_WP90/O', 'mvaTOP/F', 'sip3d/F','lostHits/I','convVeto/I','dxy/F','dz/F','charge/I','deltaEtaSC/F','mediumId/I','eleIndex/I','muIndex/I']
lepVarNames      = [x.split('/')[0] for x in lepVars]

# Training variables
read_variables = [\
                    "nBTag/I",
                    "nJetGood/I",
                    "nlep/I",
                    "m3/F",
                    "JetGood[%s]"%(",".join(jetVars)),
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

from TMB.Tools.objectSelection import isBJet
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
     "mva_jet0_btagDeepB"        :(lambda event, sample: event.JetGood_btagDeepB[0]   if (event.nJetGood >=1 and event.JetGood_btagDeepB[0]>-10) else -10),
     "mva_jet1_pt"               :(lambda event, sample: event.JetGood_pt[1]          if event.nJetGood >=2 else 0),
     "mva_jet1_eta"              :(lambda event, sample: event.JetGood_eta[1]         if event.nJetGood >=2 else -10),
     "mva_jet1_btagDeepB"        :(lambda event, sample: event.JetGood_btagDeepB[1]   if (event.nJetGood >=2 and event.JetGood_btagDeepB[1]>-10) else -10),
     "mva_jet2_pt"               :(lambda event, sample: event.JetGood_pt[2]          if event.nJetGood >=3 else 0),
     "mva_jet2_eta"              :(lambda event, sample: event.JetGood_eta[2]         if event.nJetGood >=3 else -10),
     "mva_jet2_btagDeepB"        :(lambda event, sample: event.JetGood_btagDeepB[2]   if (event.nJetGood >=3 and event.JetGood_btagDeepB[2]>-10) else -10),

     "mva_jet3_pt"               :(lambda event, sample: event.JetGood_pt[2]          if event.nJetGood >=4 else 0),
     "mva_jet4_pt"               :(lambda event, sample: event.JetGood_pt[2]          if event.nJetGood >=5 else 0),
     "mva_jet5_pt"               :(lambda event, sample: event.JetGood_pt[2]          if event.nJetGood >=6 else 0),
     "mva_jet6_pt"               :(lambda event, sample: event.JetGood_pt[2]          if event.nJetGood >=7 else 0),
     "mva_jet7_pt"               :(lambda event, sample: event.JetGood_pt[2]          if event.nJetGood >=8 else 0),
                }

mva_vector_variables    =   {
    "mva_JetGood":  {"name":"JetGood", "vars":jetVars, "varnames":jetVarNames, "selector": (lambda jet: True), 'maxN':10} 
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
        jet_vector_var = mva_vector_variables["mva_JetGood"]
        jets = [ getObjDict( event, jet_vector_var['name']+'_', jet_vector_var['varnames'], i ) for i in range(int(getattr(event,  'n'+jet_vector_var['name']))) ]
        jets = filter( jet_vector_var['selector'], jets )
        jets =  [ [ operator.itemgetter(varname)(jet) for varname in jetVarNames] for jet in jets[:jet_vector_var['maxN']] ]
        # zero padding
        jets += [ [0.]*len(jetVarNames)]*(max(0, jet_vector_var['maxN']-len(jets))) 
        jets = np.array([jets])

        return [ flat_variables, jets ]
    else:
        return   flat_variables
