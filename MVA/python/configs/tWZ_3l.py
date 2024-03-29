#!/usr/bin/env python

# Standard imports
from operator                   import attrgetter
from math                       import pi, sqrt, cosh
import ROOT, os

# RootTools
from RootTools.core.standard     import *

# helpers
from tWZ.Tools.helpers          import deltaPhi, deltaR2, deltaR, getCollection, getObjDict
from tWZ.Tools.objectSelection  import isBJet, isAnalysisJet

# Logger
import tWZ.Tools.logger as logger
logger = logger.get_logger("INFO", logFile = None )

from Analysis.Tools.leptonJetArbitration     import cleanJetsAndLeptons

jetVars          = ['pt/F', 'chEmEF/F', 'chHEF/F', 'neEmEF/F', 'neHEF/F', 'rawFactor/F', 'eta/F', 'phi/F', 'jetId/I', 'btagDeepB/F', 'index/I','btagDeepFlavB/F', 'btagCSVV2/F', 'area/F']
jetVarNames      = [x.split('/')[0] for x in jetVars]

lstm_jets_maxN   = 10
lstm_jetVars     = ['pt/F', 'eta/F', 'phi/F', 'btagDeepFlavB/F', 'btagDeepFlavC/F', 'chEmEF/F', 'chHEF/F', 'neEmEF/F', 'neHEF/F', 'muEF/F', 'puId/F', 'qgl/F']
lstm_jetVarNames = [x.split('/')[0] for x in lstm_jetVars]

lepVars          = ['pt/F','eta/F','phi/F','pdgId/I','cutBased/I','miniPFRelIso_all/F','pfRelIso03_all/F','mvaFall17V2Iso_WP90/O', 'mvaTTH/F', 'sip3d/F','lostHits/I','convVeto/I','dxy/F','dz/F','charge/I','deltaEtaSC/F','mediumId/I','eleIndex/I','muIndex/I']
lepVarNames      = [x.split('/')[0] for x in lepVars]

# Training variables
read_variables = [\
                    "nBTag/I",
                    "nJetGood/I",
                    "JetGood[%s]"%(",".join(jetVars)),
                    "met_pt/F", "met_phi/F",
                    "lep[pt/F,eta/F,phi/F]", 
                    "Z1_pt/F", "Z1_eta/F", "Z1_phi/F", "Z1_mass/F", "Z1_cosThetaStar/F",
                    "nonZ1_l1_index/I",
                    "Jet[%s]"%(",".join(jetVars)),
                    "nJet/I",
                    #"JetGood[pt/F,eta/F,phi/F,btagDeepB/F]",
                    "nlep/I",
                    "Z1_lldPhi/F",
                    "Z1_lldR/F",
                    "l1_mvaTOP/F", "l1_mvaTOPWP/I",
                    "l2_mvaTOP/F", "l2_mvaTOPWP/I",
                    "l3_mvaTOP/F", "l3_mvaTOPWP/I",
                    ]

#b tagger
b_tagger = "DeepJet"
#def flavorBin( event, sample=None):
#    event.flavorBin = 0
#
#    if      event.nMuons_tight_3l==3 and event.nElectrons_tight_3l==0: event.flavorBin = 1
#    elif    event.nMuons_tight_3l==2 and event.nElectrons_tight_3l==1: event.flavorBin = 2
#    elif    event.nMuons_tight_3l==1 and event.nElectrons_tight_3l==2: event.flavorBin = 3
#    elif    event.nMuons_tight_3l==0 and event.nElectrons_tight_3l==3: event.flavorBin = 4
#sequence.append( flavorBin )

# sequence 
sequence = []

def getbJets( event, sample=None ):
#    #bJets filtern( nur 2016 )
    alljets   = getCollection( event, 'Jet', jetVarNames , 'nJet')  
    goodjets = filter( isAnalysisJet, alljets ) 
    bJets = filter(lambda j:isBJet(j, tagger=b_tagger, year=2016) and abs(j['eta'])<=2.4, goodjets)

    if len(bJets)>=1:
        event.bJet = bJets[0]
        event.bJet_Z1_deltaR      = deltaR({'eta':bJets[0]['eta'], 'phi':bJets[0]['phi']}, {'eta':event.Z1_eta, 'phi':event.Z1_phi})
        event.bJet_nonZ1l1_deltaR = deltaR({'eta':bJets[0]['eta'], 'phi':bJets[0]['phi']}, {'eta':event.lep_eta[event.nonZ1_l1_index], 'phi':event.lep_phi[event.nonZ1_l1_index]})
    else:
        event.bJet_Z1_deltaR      = -1 
        event.bJet_nonZ1l1_deltaR = -1 

sequence.append( getbJets )

def getM3l( event, sample=None):
    # get the invariant mass of the 3l system
    l = []
    for i in range(3):
        l.append(ROOT.TLorentzVector())
        l[i].SetPtEtaPhiM(event.lep_pt[i], event.lep_eta[i], event.lep_phi[i],0)
    event.m3l = (l[0] + l[1] + l[2]).M()
sequence.append( getM3l )

def getAngles(event, sample=None):
    event.nonZ1_l1_Z1_deltaPhi = deltaPhi(event.lep_phi[event.nonZ1_l1_index], event.Z1_phi)
    event.Z1_j1_deltaPhi       = deltaPhi(event.Z1_phi, event.JetGood_phi[0]) 
    event.nonZ1_l1_Z1_deltaEta = abs(event.lep_eta[event.nonZ1_l1_index] - event.Z1_eta)
    event.nonZ1_l1_Z1_deltaR   = deltaR({'eta':event.lep_eta[event.nonZ1_l1_index], 'phi':event.lep_phi[event.nonZ1_l1_index]}, {'eta':event.Z1_eta, 'phi':event.Z1_phi})
    event.jet0_Z1_deltaR       = deltaR({'eta':event.JetGood_eta[0], 'phi':event.JetGood_phi[0]}, {'eta':event.Z1_eta, 'phi':event.Z1_phi})
    event.jet0_nonZ1_l1_deltaR = deltaR({'eta':event.JetGood_eta[0], 'phi':event.JetGood_phi[0]}, {'eta':event.lep_eta[event.nonZ1_l1_index], 'phi':event.lep_phi[event.nonZ1_l1_index]})
    event.jet1_Z1_deltaR       = deltaR({'eta':event.JetGood_eta[1], 'phi':event.JetGood_phi[1]}, {'eta':event.Z1_eta, 'phi':event.Z1_phi})
    event.jet1_nonZ1_l1_deltaR = deltaR({'eta':event.JetGood_eta[1], 'phi':event.JetGood_phi[1]}, {'eta':event.lep_eta[event.nonZ1_l1_index], 'phi':event.lep_phi[event.nonZ1_l1_index]})
    event.jet2_Z1_deltaR       = deltaR({'eta':event.JetGood_eta[2], 'phi':event.JetGood_phi[2]}, {'eta':event.Z1_eta, 'phi':event.Z1_phi})
    event.jet2_nonZ1_l1_deltaR = deltaR({'eta':event.JetGood_eta[2], 'phi':event.JetGood_phi[2]}, {'eta':event.lep_eta[event.nonZ1_l1_index], 'phi':event.lep_phi[event.nonZ1_l1_index]})

sequence.append( getAngles )

def getWpt( event, sample=None):
    # get the lepton and met
    lepton  = ROOT.TLorentzVector()
    met     = ROOT.TLorentzVector()
    lepton.SetPtEtaPhiM(event.lep_pt[event.nonZ1_l1_index], event.lep_eta[event.nonZ1_l1_index], event.lep_phi[event.nonZ1_l1_index], 0)
    met.SetPtEtaPhiM(event.met_pt, 0, event.met_phi, 0)
    # get the W boson candidate
    W   = lepton + met
    event.W_pt = W.Pt()
sequence.append( getWpt )

def forwardJets( event, sample=None ):
    #jets einlesen (in case of MC also reat the index of the genjet)
    alljets   = getCollection( event, 'Jet', jetVarNames, 'nJet')  
    alljets.sort( key = lambda j: -j['pt'] )
    leptons   = getCollection(event, "lep", lepVarNames, 'nlep') 
    # clean against good leptons
    clean_jets,_ = cleanJetsAndLeptons( alljets, leptons )
    # filter pt, but not eta 
    jets_no_eta         = filter(lambda j:j['pt']>30, clean_jets)
    if jets_no_eta: 
        event.maxAbsEta_of_pt30jets = max( [ abs(j['eta']) for j in jets_no_eta ])
    else:
        event.maxAbsEta_of_pt30jets = -1

sequence.append( forwardJets )

## met, ht, nonZ1_pt/eta, Z1_pt, nJet, nBTag, lep1_eta
#mva_variables =  {
#                    "met_pt"    :attrgetter("met_pt"), # copy
#                    "ht"        :attrgetter("ht"), # copy
#                    "lnonZ1_pt" :(lambda event: event.lep_pt[event.nonZ1_l1_index_4l]),
#                    "lnonZ1_eta":(lambda event: event.lep_eta[event.nonZ1_l1_index_4l]),
#                    "Z1_pt_4l"  :attrgetter("Z1_pt_4l"),
##                    "lep1_pt"   :(lambda event: event.lep_pt[0]),
##                    "lep2_pt"   :(lambda event: event.lep_pt[1]),
#                    "lep1_eta"  :(lambda event: event.lep_eta[0]),
##                    "lep2_eta"  :(lambda event: event.lep_eta[1]),
#                    "nJetGood":attrgetter("nJetGood"),
#                    "nBTag"     :attrgetter("nBTag"),      
##                    "yield"     :(lambda event: event.flavorBin),
##                    "jet1_pt"   :(lambda event: event.jet_pt[0]),
##                    "nLepLoose":(lambda event: event.nlep),
#                    #"myvar1" :(lambda event: event.nBTag), # calculate on the fly
#                    #"myvar2" :(lambda event: event.myFancyVar), # from sequence

all_mva_variables = {

# global event properties     
     "mva_ht"                    :(lambda event, sample: sum( [event.JetGood_pt[i] for i in range(event.nJetGood) ])),
     "mva_met_pt"                :(lambda event, sample: event.met_pt),
     "mva_m3l"                   :(lambda event, sample: event.m3l),
     "mva_nJetGood"              :(lambda event, sample: event.nJetGood),
     "mva_nBTag"                 :(lambda event, sample: event.nBTag),

# jet kinmatics
     "mva_jet0_pt"               :(lambda event, sample: event.JetGood_pt[0]          if event.nJetGood >=1 else 0),
     "mva_jet0_eta"              :(lambda event, sample: event.JetGood_eta[0]         if event.nJetGood >=1 else -10),
     "mva_jet0_btagDeepB"        :(lambda event, sample: event.JetGood_btagDeepB[0] if (event.nJetGood >=1 and event.JetGood_btagDeepB[0]>-10) else -10),
     "mva_jet1_pt"               :(lambda event, sample: event.JetGood_pt[1]          if event.nJetGood >=2 else 0),
     "mva_jet1_eta"              :(lambda event, sample: event.JetGood_eta[1]         if event.nJetGood >=2 else -10),
     "mva_jet1_btagDeepB"        :(lambda event, sample: event.JetGood_btagDeepB[1] if (event.nJetGood >=2 and event.JetGood_btagDeepB[1]>-10) else -10),
     "mva_jet2_pt"               :(lambda event, sample: event.JetGood_pt[2]          if event.nJetGood >=3 else 0),
     "mva_jet2_eta"              :(lambda event, sample: event.JetGood_eta[2]         if event.nJetGood >=3 else -10),
     "mva_jet2_btagDeepB"        :(lambda event, sample: event.JetGood_btagDeepB[2]   if (event.nJetGood >=3 and event.JetGood_btagDeepB[1]>-10) else -10),

# Z1 kinematics
     "mva_Z1_pt"                 :(lambda event, sample: event.Z1_pt),
     "mva_Z1_eta"                :(lambda event, sample: event.Z1_eta),
     "mva_Z1_cosThetaStar"       :(lambda event, sample: event.Z1_cosThetaStar),

# extra lepton kinematics
     "mva_lnonZ1_pt"             :(lambda event, sample: event.lep_pt[event.nonZ1_l1_index]),
     "mva_lnonZ1_eta"            :(lambda event, sample: event.lep_eta[event.nonZ1_l1_index]),

# leptonic W     
     "mva_W_pt"                  :(lambda event, sample: event.W_pt),

# Z1 vs. other objects
     "mva_nonZ1_l1_Z1_deltaPhi"  :(lambda event, sample: event.nonZ1_l1_Z1_deltaPhi     if event.nlep >= 1 else -1 ),
     "mva_nonZ1_l1_Z1_deltaR"    :(lambda event, sample: event.nonZ1_l1_Z1_deltaR),
     
     "mva_jet0_Z1_deltaR"        :(lambda event, sample: event.jet0_Z1_deltaR         if event.nJetGood >=1 else -1),
     "mva_jet1_Z1_deltaR"        :(lambda event, sample: event.jet1_Z1_deltaR         if event.nJetGood >=2 else -1),
     "mva_jet2_Z1_deltaR"        :(lambda event, sample: event.jet2_Z1_deltaR         if event.nJetGood >=3 else -1),

# nonZ1_l1 vs. other objects
     "mva_jet0_nonZl1_deltaR"    :(lambda event, sample: event.jet0_nonZ1_l1_deltaR    if event.nJetGood >=1 else -1),
     "mva_jet1_nonZl1_deltaR"    :(lambda event, sample: event.jet1_nonZ1_l1_deltaR    if event.nJetGood >=2 else -1),

     "mva_bJet_Z1_deltaR"        :(lambda event, sample: event.bJet_Z1_deltaR),
     "mva_bJet_non_Z1l1_deltaR"  :(lambda event, sample: event.bJet_nonZ1l1_deltaR),

     "mva_maxAbsEta_of_pt30jets" :(lambda event, sample: event.maxAbsEta_of_pt30jets),

#leptonMVA 
     "mva_l1_mvaTOP"             :(lambda event, sample: event.l1_mvaTOP),   
     "mva_l2_mvaTOP"             :(lambda event, sample: event.l2_mvaTOP),   
     "mva_l3_mvaTOP"             :(lambda event, sample: event.l3_mvaTOP),   

     "mva_l1_mvaTOPWP"           :(lambda event, sample: event.l1_mvaTOPWP),   
     "mva_l2_mvaTOPWP"           :(lambda event, sample: event.l2_mvaTOPWP),   
     "mva_l3_mvaTOPWP"           :(lambda event, sample: event.l3_mvaTOPWP),   
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

# make predictions to be used with keras.predict
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
from tWZ.samples.nanoTuples_Summer16_nanoAODv6_private_postProcessed import *
#use only Summer16
training_samples = [ TWZ_NLO_DR, TTZ ,WZ ]

assert len(training_samples)==len(set([s.name for s in training_samples])), "training_samples names are not unique!"

# training selection

from TMB.Tools.cutInterpreter import cutInterpreter
selectionString = cutInterpreter.cutString( 'trilepM' )


