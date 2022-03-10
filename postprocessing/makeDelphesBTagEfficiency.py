''' compute Delphes b-tagging efficiency
'''

# TMB
import TMB.Samples.pp_gen_v10 as samples

# Standard imports
import ROOT
import array
import pickle
import os
import itertools
import operator

# binning 
ptBorders = [30, 50, 70, 100, 140, 200, 300, 600, 10**6]
absEtaBorders = [0, .5, 1., 2., 2.5]
ptBins    = [ (ptBorders[i], ptBorders[i+1]) for i in range(len(ptBorders)-1) ]
absEtaBins    = [ (absEtaBorders[i], absEtaBorders[i+1]) for i in range(len(absEtaBorders)-1) ]


filename = "/groups/hephy/cms/robert.schoefbeck/TMB/bTagEff/Delphes_pp_gen_v10.pkl"

if __name__=="__main__":

    if os.path.exists( filename ):
        print ( "File exists: %s. Not overwriting." % filename )
    else:

        # jet selection
        jet_selection = "recoJet_pt>30&&abs(recoJet_eta)<2.5&&(recoJet_nNeutrals>0&&recoJet_nCharged>=1||abs(recoJet_eta)>2.4)"

        # select sample
        trueBsample = samples.DYBBJets

        h_trueBEff = ROOT.TH2F("h_trueBEff", "h_trueBEff", len(absEtaBorders)-1, array.array('d', absEtaBorders), len(ptBorders)-1, array.array('d', ptBorders))
        h_trueBDen = ROOT.TH2F("h_trueBDen", "h_trueBDen", len(absEtaBorders)-1, array.array('d', absEtaBorders), len(ptBorders)-1, array.array('d', ptBorders))
        trueBsample.chain.Draw("recoJet_pt:abs(recoJet_eta)>>h_trueBDen", jet_selection+"&&recoJet_matchGenBJet")
        trueBsample.chain.Draw("recoJet_pt:abs(recoJet_eta)>>h_trueBEff", jet_selection+"&&recoJet_matchGenBJet&&recoJet_bTag")

        h_trueBEff.Divide(h_trueBDen)

        trueNonBsample = samples.DYJets

        h_trueNonBEff = ROOT.TH2F("h_trueNonBEff", "h_trueNonBEff", len(absEtaBorders)-1, array.array('d', absEtaBorders), len(ptBorders)-1, array.array('d', ptBorders))
        h_trueNonBDen = ROOT.TH2F("h_trueNonBDen", "h_trueNonBDen", len(absEtaBorders)-1, array.array('d', absEtaBorders), len(ptBorders)-1, array.array('d', ptBorders))
        trueNonBsample.chain.Draw("recoJet_pt:abs(recoJet_eta)>>h_trueNonBDen", jet_selection+"&&(!recoJet_matchGenBJet)")
        trueNonBsample.chain.Draw("recoJet_pt:abs(recoJet_eta)>>h_trueNonBEff", jet_selection+"&&(!recoJet_matchGenBJet)&&recoJet_bTag")

        h_trueNonBEff.Divide(h_trueNonBDen)

        trueBEff = {}
        for i_ptBin, ptBin in enumerate(ptBins):
            trueBEff[ptBin] = {}
            for i_absEtaBin, absEtaBin in enumerate(absEtaBins):
                trueBEff[ptBin][absEtaBin] =  h_trueBEff.GetBinContent(h_trueBEff.GetBin(i_absEtaBin+1,i_ptBin+1))
        trueNonBEff = {}
        for i_ptBin, ptBin in enumerate(ptBins):
            trueNonBEff[ptBin] = {}
            for i_absEtaBin, absEtaBin in enumerate(absEtaBins):
                trueNonBEff[ptBin][absEtaBin] =  h_trueNonBEff.GetBinContent(h_trueNonBEff.GetBin(i_absEtaBin+1,i_ptBin+1))

        pickle.dump( [ trueBEff, trueNonBEff], file( filename, 'w') )

trueBEff, trueNonBEff = pickle.load( file(filename) )

#https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagSFMethods
import numpy as np

trueBJets = [
    {'pt':100, 'eta':0},
    {'pt':110, 'eta':1},
    {'pt':120, 'eta':1},
]
trueNonBJets = [
    {'pt':190, 'eta':0},
    {'pt':180, 'eta':2},
]

trueBEff_    = np.array( [ [trueBEff[ptb][etab] for j, etab in enumerate(absEtaBins)] for i, ptb in enumerate(ptBins)] )
trueNonBEff_ = np.array( [ [trueNonBEff[ptb][etab] for j, etab in enumerate(absEtaBins)] for i, ptb in enumerate(ptBins)] )
#def getBTagWeights( trueBJets, trueNonBJets ):


def getBJetProbabilities( trueBJets, trueNonBJets, nbtag_min = 2):
#if True:
#    nbtag = 2

    alljets = trueBJets+trueNonBJets
    pt_sort   = np.argsort([-j['pt'] for j in alljets ])
    b_pts     = -1 + np.digitize(np.array( [j['pt'] for j in trueBJets ]), ptBorders )
    nonb_pts  = -1 + np.digitize(np.array( [j['pt'] for j in trueNonBJets ]), ptBorders )
    b_etas    = -1 + np.digitize(np.array( [abs(j['eta']) for j in trueBJets ]), absEtaBorders )
    nonb_etas = -1 + np.digitize(np.array( [abs(j['eta']) for j in trueNonBJets ]), absEtaBorders )

    tagging_efficiences   = np.concatenate( (trueBEff_[b_pts,b_etas], trueNonBEff_[nonb_pts,nonb_etas] ) )[pt_sort]
    tagging_inefficiences = 1-tagging_efficiences

    total = 0
    njet = len(tagging_efficiences)

    cumulative_leadingjet_predictions = {comb:0 for comb in itertools.combinations(range(njet), nbtag)}

    # compute the probability to select each combination of length "nbtag_" as the b-tagged jets 
    for nbtag_ in range(0, njet+1):

        tags   = list(itertools.combinations(range(njet), nbtag_))
        notags = list(reversed(list(itertools.combinations(range(njet), njet-nbtag_))))

        tag_probabilities   = map( lambda l: reduce( operator.mul, l), tagging_efficiences[np.array(tags)]) if nbtag_>0 else [1]

        notag_probabilities = map( lambda l: reduce( operator.mul, l), 1-tagging_efficiences[np.array(notags)]) if nbtag_<njet else [1]
        probabilities = np.array(tag_probabilities)*np.array(notag_probabilities)

        total+=np.sum(probabilities)
        #print "nbtag_", nbtag_, np.sum(np.array(tag_probabilities)*np.array(notag_probabilities)), tags, tag_probabilities,probabilities#, notags

        # accumulate the probabilities by restricting to the first nbtag jets because you don't care for the rest
        for i_tag, tag in enumerate(tags):
            if len(tag)<nbtag_min: continue
            cumulative_leadingjet_predictions[tag[:nbtag_min]] += probabilities[i_tag] 

    if abs(1-total)>0.01: raise RuntimeError("BTag probabilities don't add up.")
    return cumulative_leadingjet_predictions, [alljets[j] for j in pt_sort]
