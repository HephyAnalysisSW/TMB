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

jetVars          = ['pt/F', 'eta/F', 'phi/F', 'btagDeepB/F', 'jetId/I', 'btagDeepFlavB/F']

jetVarNames      = [x.split('/')[0] for x in jetVars]

lepVars          = ['pt/F','eta/F','phi/F','pdgId/I','cutBased/I','miniPFRelIso_all/F','pfRelIso03_all/F','mvaFall17V2Iso_WP90/O', 'mvaTOP/F', 'sip3d/F','lostHits/I','convVeto/I','dxy/F','dz/F','charge/I','deltaEtaSC/F','mediumId/I','eleIndex/I','muIndex/I']
lepVarNames      = [x.split('/')[0] for x in lepVars]

# Training variables
read_variables = [\
                    "weight/F", # the lumi weight
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
                    "l1_pdgId/I",
                    "photon_pt/F",
                    "photon_eta/F",
                    "photon_phi/F",
                    "photonJetdR/F", "photonLepdR/F",
                    "np/I", VectorTreeVariable.fromString("p[C/F]", nMax=200),
                    ]

weight_variables = ['ctZ', 'cWWW']
max_order        = 2

directory = "/groups/hephy/cms/"
#directory = "/scratch-cbe/users/"

ttg    = Sample.fromDirectory("ttg",    [ '%s/robert.schoefbeck/tWZ/nanoTuples/tWZ_nAODv6_private_v6/2018/singlelep-photon/ttG_noFullyHad_fast'%directory])
wg    = Sample.fromDirectory("wg",      [ '%s/robert.schoefbeck/tWZ/nanoTuples/tWZ_nAODv6_private_v6/2018/singlelep-photon/WGToLNu_fast'%directory ])
ttg.reweight_pkl = "/eos/vbc/group/cms/robert.schoefbeck/gridpacks/v6/ttG_noFullyHad_reweight_card.pkl" 
wg.reweight_pkl = "/eos/vbc/group/cms/robert.schoefbeck/gridpacks/v6/WGToLNu_reweight_card.pkl" 
training_samples = [ ttg, wg]

assert len(training_samples)==len(set([s.name for s in training_samples])), "training_samples names are not unique!"

# Add weight infos
for sample in training_samples:
    sample.weightInfo = WeightInfo( sample.reweight_pkl )
    sample.weightInfo.set_order(2)

# make all combinations
weight_derivative_combinations = []
for i_comb, comb in enumerate(ttg.weightInfo.make_combinations(weight_variables, max_order)):
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

import TMB.Tools.VV_angles          as VV_angles
import random
def make_VV( event, sample ):
    ##AN2019_059_v8 p22
    #mW      = 80.4
    #mt2     = 2*event.l1_pt*event.met_pt*(1-cos(event.met_phi-event.l1_phi))
    #delta   = sqrt((mW**2-mt2)/(2.*event.l1_pt*event.met_pt)) 
    #if mW**2<mt2: 
    #    random_sign  = 2*(-.5+(random.random()>0.5))
    #    eta_neu = l1_eta + random_sign*log(1.+delta*sqrt(2+delta**2)+delta**2)
    #else:
    #    eta_neu = l1_eta

    # A. Wulzer lep decay angles in 2007.10356:  The latter angles are in the rest frame of each boson and they are 
    # defined as those of the final fermion or anti-fermion with helicity +1/2 (e.g. the l+ in the case
    # of a W+ and the nu-bar for a W-), denoted as f_+ in the Fig. 6

    lep_4 = ROOT.TLorentzVector()
    lep_4.SetPtEtaPhiM(event.l1_pt, event.l1_eta, event.l1_phi, 0)

    random_number = ROOT.gRandom.Uniform()
    neu_4         = VV_angles.neutrino_mom( lep_4, event.met_pt, event.met_phi, random_number )

    # the helicity+ fermion is the l+ (from a W+), otherwise it's the neutrino
    lep_m, lep_p  = (neu_4, lep_4) if event.l1_pdgId<0 else (lep_4, neu_4)

    gamma_4 = ROOT.TLorentzVector()
    gamma_4.SetPtEtaPhiM(event.photon_pt, event.photon_eta, event.photon_phi, 0)

    event.thetaW = VV_angles.gettheta(lep_m, lep_p, gamma_4)
    event.Theta  = VV_angles.getTheta(lep_m, lep_p, gamma_4)
    event.phiW   = VV_angles.getphi(lep_m, lep_p, gamma_4)

    #print "MT", sqrt(2*event.l1_pt*event.met_pt*(1-cos(event.l1_phi-event.met_phi)))
    #lep_4.Print()
    #neu_4.Print()
    #gamma_4.Print()

    #print event.thetaW, event.Theta, event.phiW
    #print

sequence.append( make_VV )

all_mva_variables = {

# global event properties     
     "mva_ht"                    :(lambda event, sample: min( sum( [event.JetGood_pt[i] for i in range(event.nJetGood)]), 2000) ),
     "mva_mT"                    :(lambda event, sample: min( sqrt(2*event.l1_pt*event.met_pt*(1-cos(event.l1_phi-event.met_phi))), 1000) ),
     "mva_m3"                    :(lambda event, sample: event.m3 if event.nJetGood >=3 else 0),
     "mva_met_pt"                :(lambda event, sample: event.met_pt),
     "mva_nJetGood"              :(lambda event, sample: min(event.nJetGood,8) ),
     "mva_nBTag"                 :(lambda event, sample: event.nBTag),
     "mva_l1_pt"                 :(lambda event, sample: min(event.l1_pt, 800)),
     "mva_l1_eta"                :(lambda event, sample: event.l1_eta),
     "mva_l1_phi"                :(lambda event, sample: event.l1_phi),
     "mva_photon_pt"             :(lambda event, sample: min(event.photon_pt,400)),
     "mva_photon_eta"            :(lambda event, sample: event.photon_eta),
     "mva_photon_phi"            :(lambda event, sample: event.photon_phi),
#VV angles
     "mva_thetaW"                :(lambda event, sample: event.thetaW),
     "mva_Theta"                 :(lambda event, sample: event.Theta),
     "mva_phiW"                  :(lambda event, sample: event.phiW),
#lep vs. gamma
     "mva_g1l1DR"                :(lambda event, sample: max(7., deltaR({'phi':event.l1_phi,'eta':event.l1_eta}, {'phi':event.photon_phi,'eta':event.photon_eta}) )),
     "mva_g1l1DPhi"              :(lambda event, sample: deltaPhi(event.l1_phi, event.photon_phi) ),
     "mva_g1l1DEta"              :(lambda event, sample: abs(event.l1_eta-event.photon_eta)),

# jet kinmatics
     "mva_jet0_pt"               :(lambda event, sample: min( event.JetGood_pt[0]     if event.nJetGood >=1 else 0, 1500)),
     "mva_jet0_eta"              :(lambda event, sample: event.JetGood_eta[0]         if event.nJetGood >=1 else -6),
     "mva_jet0_btagDeepB"        :(lambda event, sample: event.JetGood_btagDeepB[0]   if (event.nJetGood >=1 and event.JetGood_btagDeepB[0]>-10) else -3),
     "mva_jet1_pt"               :(lambda event, sample: min( event.JetGood_pt[1]          if event.nJetGood >=2 else 0,1000)),
     "mva_jet1_eta"              :(lambda event, sample: event.JetGood_eta[1]         if event.nJetGood >=2 else -6),
     "mva_jet1_btagDeepB"        :(lambda event, sample: event.JetGood_btagDeepB[1]   if (event.nJetGood >=2 and event.JetGood_btagDeepB[1]>-10) else -3),
     "mva_jet2_pt"               :(lambda event, sample: min( event.JetGood_pt[2]          if event.nJetGood >=3 else 0,1000)),
     "mva_jet2_eta"              :(lambda event, sample: event.JetGood_eta[2]         if event.nJetGood >=3 else -6),
     "mva_jet2_btagDeepB"        :(lambda event, sample: event.JetGood_btagDeepB[2]   if (event.nJetGood >=3 and event.JetGood_btagDeepB[2]>-10) else -3),
     "mva_photonJetdR"           :(lambda event, sample: min(event.photonJetdR,7)),
     "mva_photonLepdR"           :(lambda event, sample: event.photonLepdR),
                }

plot_options = {
     "mva_ht"              :{'tex':'H_{T} (GeV)'}, 
     "mva_mT"              :{'tex':'M_{T} (GeV)'}, 
     "mva_m3"              :{'tex':'M_3 (GeV)'},
     "mva_met_pt"          :{'tex':'p_{T}^{miss}'},
     "mva_nJetGood"        :{'tex':'N_{jet}',   'binning':[8,0,8]},
     "mva_nBTag"           :{'tex':'N_{b-tag}', 'binning':[4,0,4]},
     "mva_l1_pt"           :{'tex':'p_{T}(l) (GeV)'},
     "mva_l1_eta"          :{'tex':'#eta (l)'},
     "mva_l1_phi"          :{'tex':'#phi (l)'},
     "mva_thetaW"          :{'tex':'#theta(W)', 'binning':[30,0,pi]},
     "mva_Theta"           :{'tex':'#Theta', 'binning':[30,0,pi]},
     "mva_phiW"            :{'tex':'#phi (W)', 'binning':[30,-pi,pi]},
     "mva_photon_pt"       :{'tex':'p_{T} (#gamma) (GeV)'},
     "mva_photon_eta"      :{'tex':'#eta (#gamma)'},
     "mva_photon_phi"      :{'tex':'#phi (#gamma)'},
# jet kinmatics
     "mva_jet0_pt"         :{'tex':'p_{T}(j_{0}) (GeV)'},
     "mva_jet0_eta"        :{'tex':'#eta (j_{0})'},
     "mva_jet0_btagDeepB"  :{'tex':'b-jet disc. of j_{0}'}, 
     "mva_jet1_pt"         :{'tex':'p_{T}(j_{1}) (GeV)'},
     "mva_jet1_eta"        :{'tex':'#eta (j_{1})'},
     "mva_jet1_btagDeepB"  :{'tex':'b-jet disc. of j_{1}'}, 
     "mva_jet2_pt"         :{'tex':'p_{T}(j_{2}) (GeV)'},
     "mva_jet2_eta"        :{'tex':'#eta (j_{2})'},
     "mva_jet2_btagDeepB"  :{'tex':'b-jet disc. of j_{2}'}, 
     "mva_photonJetdR"     :{'tex':'min #Delta R(#gamma, j)'},
     "mva_photonLepdR"     :{'tex':'min #Delta R(#gamma, l)'}
}

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
from TMB.Tools.cutInterpreter import cutInterpreter
selectionString = cutInterpreter.cutString( 'photon-singlelep' )

bit_derivatives  = [ ('cWWW',), ('cWWW','cWWW'), ('ctZ',), ('ctZ','ctZ')]

bit_cfg = { derivative : { 'n_trees': 100,
            'max_depth'     : 4,
            'learning_rate' : 0.25,
            'min_size'      : 50,
            'clip_score_quantile': None,
            'calibrated'    : False,} for derivative in bit_derivatives }

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
