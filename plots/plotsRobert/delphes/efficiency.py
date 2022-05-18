import pickle
import numpy as np

features, weight_derivatives, predictions, test_statistic_predictions = pickle.load(file('/groups/hephy/cms/robert.schoefbeck/TMB/results/delphes-v6_MD5/ZH/dilep-ZHJet-onZ-onH-even_all/data.pkl'))

for b in ['pTVMCCut_1', 'pTVMCCut_2', 'pTVMCCut_3', 'pTVMCCut_4', 'pTVMCCut_5']:
   
    bkg_yield_tot = 0. 
    bkg_yield = 0. 
    for bkg in ['DYJets_HT_comb', 'DYBBJets_comb']:
        bkg_yield_tot += np.sum(weight_derivatives[bkg][()][test_statistic_predictions[bkg][b]>-float('inf')])
        bkg_yield     += np.sum(weight_derivatives[bkg][()][test_statistic_predictions[bkg][b]>0])

    sig_yield_tot = np.sum(weight_derivatives['ZH'][()][test_statistic_predictions['ZH'][b]>-float('inf')])
    sig_yield     = np.sum(weight_derivatives['ZH'][()][test_statistic_predictions['ZH'][b]>0])

    print b, "sig-eff", sig_yield/sig_yield_tot, "bkg-eff", bkg_yield/bkg_yield_tot
