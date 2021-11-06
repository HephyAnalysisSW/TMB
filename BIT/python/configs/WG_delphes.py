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

photonVars          = ['pt/F','eta/F','phi/F','pdgId/I','isolationVar/F', 'isolationVarRhoCorr/F', 'minLeptonDR/F', 'minJetDR/F']
photonVarNames      = [x.split('/')[0] for x in lepVars]

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
    "nrecoPhoton/I",
    "recoPhoton[%s]"%(",".join(photonVars)),
    VectorTreeVariable.fromString( "p[C/F]", nMax=200 ),
    "lumiweight1fb/F",
]


weight_variables = ['cWWW']
max_order        = 2

import TMB.Samples.pp_gen_v10 as samples

wg = samples.WGToLNu_ptG_binned 
training_samples = [wg]

assert len(training_samples)==len(set([s.name for s in training_samples])), "training_samples names are not unique!"

# Add weight infos
for sample in training_samples:
    sample.weightInfo = WeightInfo( sample.reweight_pkl )
    sample.weightInfo.set_order(2)

# make all combinations
weight_derivative_combinations = []
for i_comb, comb in enumerate(wg.weightInfo.make_combinations(weight_variables, max_order)):
    weight_derivative_combinations.append(comb)

for sample in training_samples:
    sample.weight_derivatives = []
    for i_comb, comb in enumerate(weight_derivative_combinations):
        #print name, i_comb, comb, weightInfo.get_diff_weight_string(comb)
        if all( v in sample.weightInfo.variables for v in comb ):

            # func_ takes care of p_C. We also normalize with the lumi-weight 'weight' from the original sample
            func_            = sample.weightInfo.get_diff_weight_func(comb)
            weight = {}
            weight['string'] = "weight*("+sample.weightInfo.get_diff_weight_string(comb)+")"
            weight['func']   = lambda event, sample, func_=func_: func_(event,sample)*event.weight 
            weight['name']   = '_'.join(comb)
            weight['comb']   = comb
            sample.weight_derivatives.append( weight )
        else:
            print "Warning! Derivative %r put to zero in sample %s because some WC are missing."%( comb, sample.name )
            weight = {}
            weight['string'] = '(0.)' 
            weight['func']   = lambda event, sample: 0. 
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

sequence = []

def getBJetindex( event ):
    ''' return highest pt bjet'''
    for i in range(event.nrecoJet):
        if event.recoJet_bTag[i]>=1:
            return i
    return -1

from TMB.Tools.objectSelection import isBJet
def make_jets( event, sample ):
    event.jets     = [getObjDict(event, 'recoJet_', jetVarNames, i) for i in range(int(event.nrecoJet))]
    event.bJets    = filter(lambda j:j['bTag']>=1 and abs(j['eta'])<=2.4    , event.jets)
sequence.append( make_jets )


#sequence = []
#
#import TMB.Tools.VV_angles          as VV_angles
#import random
#def make_VV( event, sample ):
#    ##AN2019_059_v8 p22
#    #mW      = 80.4
#    #mt2     = 2*event.l1_pt*event.met_pt*(1-cos(event.met_phi-event.l1_phi))
#    #delta   = sqrt((mW**2-mt2)/(2.*event.l1_pt*event.met_pt)) 
#    #if mW**2<mt2: 
#    #    random_sign  = 2*(-.5+(random.random()>0.5))
#    #    eta_neu = l1_eta + random_sign*log(1.+delta*sqrt(2+delta**2)+delta**2)
#    #else:
#    #    eta_neu = l1_eta
#
#    # A. Wulzer lep decay angles in 2007.10356:  The latter angles are in the rest frame of each boson and they are 
#    # defined as those of the final fermion or anti-fermion with helicity +1/2 (e.g. the l+ in the case
#    # of a W+ and the nu-bar for a W-), denoted as f_+ in the Fig. 6
#
#    lep_4 = ROOT.TLorentzVector()
#    lep_4.SetPtEtaPhiM(event.l1_pt, event.l1_eta, event.l1_phi, 0)
#
#    random_number = ROOT.gRandom.Uniform()
#    neu_4         = VV_angles.neutrino_mom( lep_4, event.met_pt, event.met_phi, random_number )
#
#    # the helicity+ fermion is the l+ (from a W+), otherwise it's the neutrino
#    lep_m, lep_p  = (neu_4, lep_4) if event.l1_pdgId<0 else (lep_4, neu_4)
#
#    gamma_4 = ROOT.TLorentzVector()
#    gamma_4.SetPtEtaPhiM(event.photon_pt, event.photon_eta, event.photon_phi, 0)
#
#    event.thetaW = VV_angles.gettheta(lep_m, lep_p, gamma_4)
#    event.Theta  = VV_angles.getTheta(lep_m, lep_p, gamma_4)
#    event.phiW   = VV_angles.getphi(lep_m, lep_p, gamma_4)
#
#    #print "MT", sqrt(2*event.l1_pt*event.met_pt*(1-cos(event.l1_phi-event.met_phi)))
#    #lep_4.Print()
#    #neu_4.Print()
#    #gamma_4.Print()
#
#    #print event.thetaW, event.Theta, event.phiW
#    #print
#
#sequence.append( make_VV )
#
all_mva_variables = {

## global event properties     
#     "mva_mT"                    :(lambda event, sample: min( sqrt(2*event.l1_pt*event.met_pt*(1-cos(event.l1_phi-event.met_phi))), 1000) ),
#     "mva_m3"                    :(lambda event, sample: event.m3 if event.nJetGood >=3 else 0),
#     "mva_met_pt"                :(lambda event, sample: event.met_pt),
#     "mva_nJetGood"              :(lambda event, sample: min(event.nJetGood,8) ),
#     "mva_nBTag"                 :(lambda event, sample: event.nBTag),
#     "mva_l1_pt"                 :(lambda event, sample: min(event.l1_pt, 800)),
#     "mva_l1_eta"                :(lambda event, sample: event.l1_eta),
#     "mva_l1_phi"                :(lambda event, sample: event.l1_phi),
#     "mva_photon_pt"             :(lambda event, sample: min(event.photon_pt,400)),
#     "mva_photon_eta"            :(lambda event, sample: event.photon_eta),
#     "mva_photon_phi"            :(lambda event, sample: event.photon_phi),
##VV angles
#     "mva_thetaW"                :(lambda event, sample: event.thetaW),
#     "mva_Theta"                 :(lambda event, sample: event.Theta),
#     "mva_phiW"                  :(lambda event, sample: event.phiW),
##lep vs. gamma
#     "mva_g1l1DR"                :(lambda event, sample: max(7., deltaR({'phi':event.l1_phi,'eta':event.l1_eta}, {'phi':event.photon_phi,'eta':event.photon_eta}) )),
#     "mva_g1l1DPhi"              :(lambda event, sample: deltaPhi(event.l1_phi, event.photon_phi) ),
#     "mva_g1l1DEta"              :(lambda event, sample: abs(event.l1_eta-event.photon_eta)),

# global event properties
     "mva_ht"                    :(lambda event, sample: sum( [event.recoJet_pt[i] for i in range(event.nrecoJet) ])),
     "mva_met_pt"                :(lambda event, sample: event.recoMet_pt),
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


#     "mva_photonJetdR"           :(lambda event, sample: min(event.photonJetdR,7)),
#     "mva_photonLepdR"           :(lambda event, sample: event.photonLepdR),
                }
#
#plot_options = {
#     "mva_ht"              :{'tex':'H_{T} (GeV)'}, 
#     "mva_mT"              :{'tex':'M_{T} (GeV)'}, 
#     "mva_m3"              :{'tex':'M_3 (GeV)'},
#     "mva_met_pt"          :{'tex':'p_{T}^{miss}'},
#     "mva_nJetGood"        :{'tex':'N_{jet}',   'binning':[8,0,8]},
#     "mva_nBTag"           :{'tex':'N_{b-tag}', 'binning':[4,0,4]},
#     "mva_l1_pt"           :{'tex':'p_{T}(l) (GeV)'},
#     "mva_l1_eta"          :{'tex':'#eta (l)'},
#     "mva_l1_phi"          :{'tex':'#phi (l)'},
#     "mva_thetaW"          :{'tex':'#theta(W)', 'binning':[30,0,pi]},
#     "mva_Theta"           :{'tex':'#Theta', 'binning':[30,0,pi]},
#     "mva_phiW"            :{'tex':'#phi (W)', 'binning':[30,-pi,pi]},
#     "mva_photon_pt"       :{'tex':'p_{T} (#gamma) (GeV)'},
#     "mva_photon_eta"      :{'tex':'#eta (#gamma)'},
#     "mva_photon_phi"      :{'tex':'#phi (#gamma)'},
## jet kinmatics
#     "mva_jet0_pt"         :{'tex':'p_{T}(j_{0}) (GeV)'},
#     "mva_jet0_eta"        :{'tex':'#eta (j_{0})'},
#     "mva_jet0_btagDeepB"  :{'tex':'b-jet disc. of j_{0}'}, 
#     "mva_jet1_pt"         :{'tex':'p_{T}(j_{1}) (GeV)'},
#     "mva_jet1_eta"        :{'tex':'#eta (j_{1})'},
#     "mva_jet1_btagDeepB"  :{'tex':'b-jet disc. of j_{1}'}, 
#     "mva_jet2_pt"         :{'tex':'p_{T}(j_{2}) (GeV)'},
#     "mva_jet2_eta"        :{'tex':'#eta (j_{2})'},
#     "mva_jet2_btagDeepB"  :{'tex':'b-jet disc. of j_{2}'}, 
#     "mva_photonJetdR"     :{'tex':'min #Delta R(#gamma, j)'},
#     "mva_photonLepdR"     :{'tex':'min #Delta R(#gamma, l)'}
#}
#
mva_vector_variables    =   {
    "weight":  { "name":"weight", "func":compute_weight_derivatives, "vars":["derivatives/F"], "varnames":["derivatives"],}# 'nMax':len(weight_derivatives)} 
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
from TMB.Tools.delphesCutInterpreter import cutInterpreter
selectionString = "(1)" #FIXME cutInterpreter.cutString( 'photon-singlelep' )

bit_derivatives  = [ ('cWWW',), ('cWWW','cWWW')]

bit_cfg = { derivative : { 'n_trees': 100,
            'max_depth'     : 4,
            'learning_rate' : 0.25,
            'min_size'      : 50,
            'clip_score_quantile': None,
            'calibrated'    : False,} for derivative in bit_derivatives }

def load(directory = '/mnt/hephy/cms/$USER/BIT/models/default/WG_delphes/', bit_derivatives=bit_derivatives):
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
