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

jetVars          = ['pt/F', 'eta/F', 'phi/F', 'btagDeepB/F', 'jetId/I', 'btagDeepFlavB/F', 'mass/F']

jetVarNames      = [x.split('/')[0] for x in jetVars]

lepVars          = ['pt/F','eta/F','phi/F','pdgId/I','cutBased/I','miniPFRelIso_all/F','pfRelIso03_all/F','mvaFall17V2Iso_WP90/O', 'mvaTOP/F', 'sip3d/F','lostHits/I','convVeto/I','dxy/F','dz/F','charge/I','deltaEtaSC/F','mediumId/I','eleIndex/I','muIndex/I']
lepVarNames      = [x.split('/')[0] for x in lepVars]

# Training variables
read_variables = [\
                    "nBTag/I",
                    "m3/F",
                    "nJetGood/I",
                    "JetGood[%s]"%(",".join(jetVars)),
                    "nlep/I",
                    "lep[%s]"%(",".join(lepVars)),
                    "met_pt/F", "met_phi/F",
                    "l1_pt/F",
                    "l1_eta/F",
                    "l1_phi/F",
                    "photon_pt/F",
                    "photon_eta/F",
                    "photon_phi/F",
                    "photonJetdR/F", "photonLepdR/F",
                    "np/I", VectorTreeVariable.fromString("p[C/F]", nMax=200),
                    ]

weight_variables = ['ctZ']
max_order        = 2

reweight_pkl     = "/eos/vbc/user/robert.schoefbeck/gridpacks/v6/ttG_noFullyHad_reweight_card.pkl" 
ttg    = Sample.fromDirectory("ttg",    ['/scratch-cbe/users/robert.schoefbeck/tWZ/nanoTuples/tWZ_nAODv6_private_v6/2018/singlelep-photon/ttG_noFullyHad_fast'])
ttg.reweight_pkl = reweight_pkl
training_samples = [ ttg ]

assert len(training_samples)==len(set([s.name for s in training_samples])), "training_samples names are not unique!"

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

all_mva_variables = {

# global event properties     
     "mva_ht"                    :(lambda event, sample: sum( [event.JetGood_pt[i] for i in range(event.nJetGood)]) ),
     "mva_mT"                    :(lambda event, sample: sqrt(2*event.l1_pt*event.met_pt*(1-cos(event.l1_phi-event.met_phi)))),
     "mva_m3"                    :(lambda event, sample: event.m3 if event.nJetGood >=3 else 0),
     "mva_met_pt"                :(lambda event, sample: event.met_pt),
     "mva_nJetGood"              :(lambda event, sample: event.nJetGood),
     "mva_nBTag"                 :(lambda event, sample: event.nBTag),
     "mva_l1_pt"                 :(lambda event, sample: event.l1_pt),
     "mva_l1_eta"                :(lambda event, sample: event.l1_eta),
     "mva_l1_phi"                :(lambda event, sample: event.l1_phi),
     "mva_photon_pt"             :(lambda event, sample: event.photon_pt),
     "mva_photon_eta"            :(lambda event, sample: event.photon_eta),
     "mva_photon_phi"            :(lambda event, sample: event.photon_phi),
# jet kinmatics
     "mva_jet0_pt"               :(lambda event, sample: event.JetGood_pt[0]          if event.nJetGood >=1 else 0),
     "mva_jet0_eta"              :(lambda event, sample: event.JetGood_eta[0]         if event.nJetGood >=1 else -10),
     "mva_jet0_btagDeepB"        :(lambda event, sample: event.JetGood_btagDeepB[0]   if (event.nJetGood >=1 and event.JetGood_btagDeepB[0]>-10) else -10),
     "mva_jet1_pt"               :(lambda event, sample: event.JetGood_pt[1]          if event.nJetGood >=2 else 0),
     "mva_jet1_eta"              :(lambda event, sample: event.JetGood_eta[1]         if event.nJetGood >=2 else -10),
     "mva_jet1_btagDeepB"        :(lambda event, sample: event.JetGood_btagDeepB[1]   if (event.nJetGood >=2 and event.JetGood_btagDeepB[1]>-10) else -10),
     "mva_jet2_pt"               :(lambda event, sample: event.JetGood_pt[2]          if event.nJetGood >=3 else 0),
     "mva_jet2_eta"              :(lambda event, sample: event.JetGood_eta[2]         if event.nJetGood >=3 else -10),
     "mva_jet2_btagDeepB"        :(lambda event, sample: event.JetGood_btagDeepB[2]   if (event.nJetGood >=3 and event.JetGood_btagDeepB[2]>-10) else -10),
     "mva_photonJetdR"           :(lambda event, sample: event.photonJetdR),
     "mva_photonLepdR"           :(lambda event, sample: event.photonLepdR),
                }

mva_vector_variables    =   {
    "weight":  { "name":"weight", "func":compute_weight_derivatives, "vars":["derivatives/F"], "varnames":["derivatives"], 'nMax':len(weight_derivatives)} 
}

## Using all variables
mva_variables_ = all_mva_variables.keys()
mva_variables_.sort()
mva_variables  = [ (key, value) for key, value in all_mva_variables.iteritems() if key in mva_variables_ ]

import numpy as np
import operator

# make predictions to be used with keras.predict
def predict_inputs( event, sample):
    return np.array([getattr( event, mva_variable) for mva_variable, _ in mva_variables])

# training selection
from TMB.Tools.cutInterpreter import cutInterpreter
selectionString = cutInterpreter.cutString( 'photon-singlelep' )

bit_cfg = { 'n_trees': 100,
            'max_depth'     : 2,
            'learning_rate' : 0.20,
            'min_size'      : 50,
    }

bit_derivatives  = [ ('ctZ',), ('ctZ','ctZ')]

def load(directory = '/mnt/hephy/cms/$USER/BIT/models/default/ttG_ctZ/', bit_derivatives=bit_derivatives):
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
