import TMB.Samples.nanoTuples_RunII_nanoAODv6_dilep_pp as samples
from RootTools.core.standard import *
import os
from TMB.Tools.cutInterpreter            import cutInterpreter
from TMB.Tools.user                      import plot_directory
import Analysis.Tools.syncer
import ROOT

selection = 'dilepL-offZ1'

sample = samples.TTLep

bGenJet    = "(GenJet_nBHadFromT+GenJet_nBHadFromTbar>=1 &&GenJet_nBHadFromW+GenJet_nBHadOther+GenJet_nCHadFromW+GenJet_nCHadOther==0)"
lightGenJet    = "(GenJet_nBHadFromT+GenJet_nBHadFromTbar+GenJet_nBHadFromW+GenJet_nBHadOther+GenJet_nCHadFromW+GenJet_nCHadOther==0)"
acceptedGenJet = "(GenJet_pt>30&&abs(GenJet_eta)<2.4)"
#acceptedGenJet = "(1)"#(GenJet_pt>30&&abs(GenJet_eta)<2.4)"
#acceptedGenJet = "(GenJet_pt>50)"
nLightGenJet   = "Sum$({gen_sel}&&{acc_sel})".format( gen_sel=lightGenJet, acc_sel=acceptedGenJet )
nBGenJet       = "Sum$({gen_sel}&&{acc_sel})".format( gen_sel=bGenJet, acc_sel=acceptedGenJet )

selectionString = "(1)"

#nlj_bt1  = sample.get1DHistoFromDraw( nLightGenJet, [7,1,8], selectionString = cutInterpreter.cutString(selection)+"&&"+nBGenJet+"==1" )
#nlj_bt2p = sample.get1DHistoFromDraw( nLightGenJet, [7,1,8], selectionString = cutInterpreter.cutString(selection)+"&&"+nBGenJet+">=2" )
#
#nlj_bt1.style = styles.lineStyle( ROOT.kBlue)
#nlj_bt2p.style = styles.lineStyle( ROOT.kRed)
#
#nlj_plot = Plot.fromHisto(name = "TTLep_lj_comp_e2p4_mc-btag", histos = [[nlj_bt1],[nlj_bt2p]], texX = "n(light jet)" , texY = "Number of Events" )
#
#for log in [True, False]:
#    plot_directory_ = os.path.join(plot_directory, 'nc_shapes',  selection, ("log" if log else "lin"))
#    plotting.draw(nlj_plot, plot_directory = plot_directory_, ratio = {1:0}, scaling={1:0}, logY = log, logX = False, 
#        #yRange = (0, 0.5) if 'Mult' not in var else (0, 15 ),  
#        #legend = [0.60,0.92-0.05*len(frac_histos),0.99,0.88]
#    )
#
#Analysis.Tools.syncer.sync()

nlj_bt1  = sample.get1DHistoFromDraw( nLightGenJet, [9,2-0.5,11-0.5], selectionString = cutInterpreter.cutString(selection+"-btag1") )
nlj_bt2p = sample.get1DHistoFromDraw( nLightGenJet, [9,2-0.5,11-0.5], selectionString = cutInterpreter.cutString(selection+"-btag2p") )

nlj_bt1.style = styles.lineStyle( ROOT.kBlue , errors=True, width=2)
nlj_bt2p.style = styles.lineStyle( ROOT.kRed, errors=True, width=2)

nlj_bt1.Scale(1./nlj_bt1.Integral())
nlj_bt2p.Scale(1./nlj_bt2p.Integral())

nlj_bt1.legendText= "n_{b-tag} = 1"
nlj_bt2p.legendText= "n_{b-tag} #geq 2"

nlj_plot = Plot.fromHisto(name = "TTLep_lj_comp_e2p4_reco-btag", histos = [[nlj_bt1],[nlj_bt2p]], texX = "Number of light jets" , texY = "Arbitrary units" )

for log in [True, False]:
    plot_directory_ = os.path.join(plot_directory, 'nc_shapes',  selection, ("log" if log else "lin"))
    plotting.draw(nlj_plot, plot_directory = plot_directory_, ratio = {1:0}, scaling={1:0}, logY = log, logX = False, yRange = (10**-6,0.9), 
        #yRange = (0, 0.5) if 'Mult' not in var else (0, 15 ),  
        legend = [0.60,0.75,0.99,0.88]
    )

Analysis.Tools.syncer.sync()
