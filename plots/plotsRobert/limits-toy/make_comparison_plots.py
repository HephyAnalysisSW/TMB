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

for pred in [
    "power_total_cHQ3_vs_cHWtil_predicted_lumi_factor_1.00_level_0.68.root",
    "power_total_cHQ3_vs_cHW_predicted_lumi_factor_1.00_level_0.68.root",
    "power_total_cHW_vs_cHWtil_predicted_lumi_factor_1.00_level_0.68.root",
        ]:
    draw_conts = []

    texX = tex[pred.split('_')[2]]
    texY = tex[pred.split('_')[4]]

    for stat in ["truth", "predicted"]:
        for CL in [ "0.68", "0.95" ]:
            f = ROOT.TFile(os.path.join( plot_dir, pred.replace("predicted", stat).replace("0.68", CL)))
            c = f.Get("ROOT.c1")
            h_2D  = c.GetListOfPrimitives().At(1)
            stuff.append(h_2D)
            if CL=="0.68" and stat=="predicted":
                pred_h_2D = h_2D

            h_2D.Smooth(1, "k5b")
            h_2D.Smooth(1, "k5b")
            conts = getContours(h_2D, 0.5)
            for cont in conts:
                if stat == "truth":
                    cont.SetLineColor(ROOT.kRed - 2)
                    cont.SetLineWidth(2)
                else:
                    cont.SetLineColor(ROOT.kBlack)
                if CL == "0.68":
                    cont.SetLineStyle(2)
                cont.SetLineWidth(2)
                draw_conts.append(cont)
                stuff.append(cont)
            
    for nBinsTestStat in [1, 2, 5]:
        for CL in [ "0.68", "0.95" ]:
            pred_binned = pred.replace(".root", "_nBinsTestStat_%i.root"%nBinsTestStat).replace("0.68", CL)
            f_binned    = ROOT.TFile(os.path.join( plot_dir_binned, pred_binned))
            c_binned    = f_binned.Get("ROOT.c1")
            stuff.append(f_binned)
            stuff.append(c_binned)
            th2d = c_binned.GetListOfPrimitives().At(1)
            th2d.Smooth(1, "k5b")
            conts = getContours( th2d, 0.5 ) 
            for binned_cont in conts:
                binned_cont.SetLineColor(ROOT.kGray+2)
                if CL=="0.68":
                    binned_cont.SetLineStyle(2)
                    #binned_cont.SetLineWidth(2)
                #spline.SetLineColor(ROOT.kGray)
                #binned_conts.append( spline )
                draw_conts.append( binned_cont )

    plot2D = Plot2D.fromHisto( name=pred.replace("predicted", "comb"), histos=[[pred_h_2D]],
        texX = texX,    
        texY = texY,    
    )

    plotting.draw2D(plot2D, plot_directory = plot_dir,
        histModifications = [lambda h:ROOT.gStyle.SetPalette(58)],
        drawObjects       = draw_conts,
        zRange            = (0.095,1),
            )

syncer.sync()
