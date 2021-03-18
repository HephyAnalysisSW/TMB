from RootTools.core.standard import *
import Analysis.Tools.syncer
import os
import ROOT
filename = "/mnt/hephy/cms/robert.schoefbeck/www/TMB/analysisPlots/TMB_4t_noData/RunII/dilepL-offZ1-njet4p-btag2p-ht500.pkl"
ROOT.gROOT.LoadMacro("$CMSSW_BASE/src/TMB/Tools/scripts/tdrstyle.C")
ROOT.setTDRStyle()

import pickle
hist_dict = pickle.load(file(filename))

mva  = 'tttt_2l_lstm'
classifier = 'TTTT_vs_TTLep_bb'
#classifier = 'TTLep_bb_vs_TTLep_cc'
bkgs = ['TTTT', 'TTLep_bb', 'TTLep_cc', 'TTLep_other'] 


niceName = {
    'TTTT':"t#bar{t}t#bar{t}",
    'TTLep_bb':"t#bar{t}b#bar{b}",
    'TTLep_cc':"t#bar{t}c#bar{c}",
    'TTLep_other':"t#bar{t} (other)",
}

ROOT.setTDRStyle()

c1 = ROOT.TCanvas()
h = {}
for bkg in bkgs:
    h[bkg] =  hist_dict[mva+'_'+bkg+'_'+classifier][0][0]
    h[bkg].Scale(1./h[bkg].Integral() )

h['TTTT']        .SetLineColor( ROOT.kOrange+1 )
h['TTTT']        .SetFillColor( ROOT.kOrange+1 )
h['TTTT']        .SetMarkerColor( ROOT.kOrange+1 )
h['TTTT']        .SetMarkerStyle(0)
h['TTLep_bb']    .SetLineColor( ROOT.kRed + 2)
h['TTLep_bb']    .SetFillColor( ROOT.kRed + 2)
h['TTLep_bb']    .SetMarkerColor( ROOT.kRed + 2)
h['TTLep_bb']    .SetMarkerStyle(0)
h['TTLep_cc']    .SetLineColor( ROOT.kRed-3 )
h['TTLep_cc']    .SetFillColor( ROOT.kRed-3 )
h['TTLep_cc']    .SetMarkerColor( ROOT.kRed-3 )
h['TTLep_cc']    .SetMarkerStyle(0)
h['TTLep_other'] .SetLineColor( ROOT.kRed-7 )
h['TTLep_other'] .SetFillColor( ROOT.kRed-7 )
h['TTLep_other'] .SetMarkerColor( ROOT.kRed-7 )
h['TTLep_other'] .SetMarkerStyle(0)

l = ROOT.TLegend(0.65, 0.65, 0.95, 0.89)
l.SetFillStyle(0)
l.SetShadowColor(ROOT.kWhite)
l.SetBorderSize(0)

l.AddEntry( h['TTLep_bb'], "t#bar{t}b#bar{b}" )
l.AddEntry( h['TTLep_cc'], "t#bar{t}c#bar{c}" )
l.AddEntry( h['TTLep_other'], "t#bar{t} (other)" )
l.AddEntry( h['TTTT'], "t#bar{t}t#bar{t}" )
c1 = ROOT.TCanvas()
h['TTLep_bb'].Draw("BOX")
h['TTLep_bb'].SetTitle("")
h['TTLep_bb'].GetXaxis().SetTitle("p(%s)"%niceName[classifier.split('_vs_')[0]])
h['TTLep_bb'].GetYaxis().SetTitle("p(%s)"%niceName[classifier.split('_vs_')[1]])

h['TTLep_cc'].Draw("BOXsame")
h['TTLep_other'].Draw("BOXsame")
h['TTTT'].Draw("BOXsame")

c1.SetLogz(0)
l.Draw()
c1.RedrawAxis()
ROOT.gStyle.SetOptStat(0)
c1.SetTitle("")
c1.Print(filename.replace('.pkl', '_'+classifier+'_mva.png'))
c1.Print(filename.replace('.pkl', '_'+classifier+'_mva.pdf'))
c1.Print(filename.replace('.pkl', '_'+classifier+'_mva.root'))
Analysis.Tools.syncer.sync()
