#!/usr/bin/env python

# Standard imports
from operator                   import attrgetter
from math                       import pi, sqrt, cosh, cos
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

jetVars          = ['pt/F', 'chEmEF/F', 'chHEF/F', 'neEmEF/F', 'neHEF/F', 'rawFactor/F', 'eta/F', 'phi/F', 'jetId/I', 'btagDeepB/F', 'btagDeepFlavB/F', 'btagCSVV2/F', 'area/F']
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
                    "np/I", "p[C/F]",
                    ]
# sequence 
sequence = []

# Fisher informations
FIs = {
    'cWWW_SM' :  { 'var': 'cWWW', 'point':{'cWWW':0}},
    'cWWW_1' :   { 'var': 'cWWW', 'point':{'cWWW':1}},
}

# initialize weight stuff
def init( event, sample ):
    if hasattr( sample, "EFT_init") and sample.EFT_init:
        return
    sample.EFT_init = True
    if hasattr( sample, "reweight_pkl"):
        sample.weightInfo = WeightInfo( sample.reweight_pkl )
        sample.weightInfo.set_order(2)
        for FI_name, FI in FIs.iteritems():
            FI['string'] = sample.weightInfo.get_fisher_weight_string( FI['var'], FI['var'], **FI['point'] )
            FI['func'] = lambda p_C: sample.weightInfo.get_fisherInformation_matrix( p_C, variables = [FI['var']], **FI['point'])
            #sample.chain.SetNotify( FI['TTreeFormula'] )
            logger.info( "Compute %s = FI(%s) at %r. The string expression is %s", FI_name, FI['var'], FI['point'], FI['string'] )

sequence.append( init )

all_mva_variables = {

# global event properties     
     "mva_ht"                    :(lambda event, sample: sum( [event.JetGood_pt[i] for i in range(event.nJetGood)]) ),
     "mva_mT"                    :(lambda event, sample: sqrt(2*event.l1_pt*event.met_pt*(1-cos(event.l1_phi-event.met_phi)))),
     "mva_m3"                    :(lambda event, sample: event.m3),
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

# Using all variables
mva_variables_ = all_mva_variables.keys()

mva_variables = {key:value for key, value in all_mva_variables.iteritems() if key in mva_variables_}
