import ROOT
import os
from RootTools.core.standard import *
import Analysis.Tools.syncer as syncer
import array

# Unbinned
plot_dir = "/groups/hephy/cms/robert.schoefbeck/www/BIT/BIT_VH_12_MD6/ZH_Nakamura/unbinned/"

# Binned
plot_dir_binned= "/groups/hephy/cms/robert.schoefbeck/www/BIT/BIT_VH_12_MD6/ZH_Nakamura/binned/"

tex = {'cHW':"C_{HW}", 'cHWtil':"C_{H#tilde{W}}", "cHQ3":"C_{HQ}^{(3)}"}
texX = tex[pred.split('_')[2]]
texY = tex[pred.split('_')[4]]

# Smoothing
#smoother = ROOT.TGraphSmooth("normal")
def getContours( h, level):
    _h     = h.Clone()
    ctmp = ROOT.TCanvas()
    _h.SetContour(1,array.array('d', [level]))
    _h.Draw("contzlist")
    _h.GetZaxis().SetRangeUser(0.0001,1)
    ctmp.Update()
    contours = ROOT.gROOT.GetListOfSpecials().FindObject("contours")
    return contours.At(0).Clone()

stuff=[]
binned_conts = []

for pred in [
    "power_total_cHQ3_vs_cHWtil_predicted_lumi_factor_1.00_level_0.68.root",
    "power_total_cHQ3_vs_cHW_predicted_lumi_factor_1.00_level_0.68.root",
    "power_total_cHW_vs_cHWtil_predicted_lumi_factor_1.00_level_0.68.root",
        ]:

    f_pred = ROOT.TFile(os.path.join( plot_dir, pred))
    c_pred = f_pred.Get("ROOT.c1")
    pred_cont1 = c_pred.GetListOfPrimitives().At(3)
    pred_cont2 = c_pred.GetListOfPrimitives().At(4)
    pred_h_2D  = c_pred.GetListOfPrimitives().At(1)

    truth = pred.replace("predicted", "truth") 

    f_truth = ROOT.TFile(os.path.join( plot_dir, truth))
    c_truth = f_truth.Get("ROOT.c1")
    truth_cont1 = c_truth.GetListOfPrimitives().At(3)
    truth_cont2 = c_truth.GetListOfPrimitives().At(4)
    truth_h_2D  = c_truth.GetListOfPrimitives().At(1)

    truth_cont1.SetLineColor(ROOT.kRed - 2)
    truth_cont2.SetLineColor(ROOT.kRed - 2)


    for nBinsTestStat in [1, 2, 5]:
        for CL in [ "0.68", "0.95" ]:
            pred_binned = pred.replace(".root", "_nBinsTestStat_%i.root"%nBinsTestStat).replace("0.68", CL)
            f_binned    = ROOT.TFile(os.path.join( plot_dir_binned, pred_binned))
            c_binned    = f_binned.Get("ROOT.c1")
            stuff.append(f_binned)
            stuff.append(c_binned)
            th2d = c_binned.GetListOfPrimitives().At(1)
            th2d.Smooth(10)
            conts = getContours( th2d, 0.5 ) 
            for binned_cont in conts:
                binned_cont.SetLineColor(ROOT.kGray)
                if CL=="0.68":
                    binned_cont.SetLineStyle(7)
                #spline.SetLineColor(ROOT.kGray)
                #binned_conts.append( spline )
                binned_conts.append( binned_cont )

    plot2D = Plot2D.fromHisto( name=pred.replace("predicted", "comb"), histos=[[pred_h_2D]],
        texX = texX,    
        texY = texY,    
    )

    plotting.draw2D(plot2D, plot_directory = plot_dir,
        histModifications = [lambda h:ROOT.gStyle.SetPalette(58)],
        drawObjects       = [pred_cont1, pred_cont2, truth_cont1, truth_cont2] + binned_conts,
            )

syncer.sync()
