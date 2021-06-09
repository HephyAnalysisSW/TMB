import ROOT
import Analysis.Tools.syncer as syncer
from TMB.Tools.WeightInfo import WeightInfo

#ver = "v6"
#sample_name = "ttG_noFullyHad_test"
#varname = "genPhoton_pt"

ver = "v6"
sample_name = "ttg_noFullyHad_test2"
varname = "genPhoton_pt"

#ver = "v4"
#sample_name = "ttg_noFullyHad"
#varname = "genPhoton_pt"
#ver = "v5"
#sample_name = "ttG_noFullyHad"
#varname = "genPhoton_pt"
#ver = "v5"
#sample_name = "ttZ01j"
#varname = "genZ_pt"
#ver = "v5"
#sample_name = "ttZ01j"
#varname = "genZ_pt"

if ver == "v4":
    import pp_gen_v4 as samples
elif ver == "v5":
    import pp_gen_v5 as samples
elif ver == "v6":
    import pp_gen_v6 as samples
    
sample = getattr(samples, sample_name)

#import pp_gen_v5 as samples
#sample = samples.ttG_noFullyHad

w = WeightInfo( sample.reweight_pkl)
w.set_order(2)
sample.chain.Draw(varname+">>hSM(50,0,500)", "("+w.get_weight_string(ctZ=0)+")")
sample.chain.Draw(varname+">>hBSM(50,0,500)", "("+w.get_weight_string(ctZ=1)+")")

c1 = ROOT.TCanvas()
ROOT.hBSM.SetLineColor(ROOT.kRed)
ROOT.hSM.SetMarkerStyle(0)
ROOT.hBSM.SetMarkerStyle(0)
ROOT.hSM.Draw("h")
ROOT.hBSM.Draw("hsame")
c1.SetLogy()
c1.Print("/mnt/hephy/cms/robert.schoefbeck/www/TMB/gen/"+ver+"_"+sample_name+"_"+varname+".png")

syncer.sync()
