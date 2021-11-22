import ROOT

from TMB.Samples.helpers import singleton as singleton

@singleton
class color():
  pass

color.data           = ROOT.kBlack
color.Zinv           = ROOT.kGray+2
color.DY             = ROOT.kAzure+4 #ROOT.kCyan+2
color.ZGamma         = ROOT.kBlue+2
color.WGamma         = ROOT.kAzure-3
color.VGamma         = ROOT.kAzure+3 #ROOT.kAzure
color.TTH            = ROOT.kAzure+8 #ROOT.kAzure
color.TWZ            = ROOT.kRed
#color.TTZ            = ROOT.kAzure+4 #ROOT.kBlack
color.TTW            = ROOT.kGreen+2 
color.TTG            = ROOT.kOrange
color.Other          = ROOT.kOrange+2 #ROOT.kViolet+5
color.TT             = ROOT.kRed-7 #ROOT.kRed+1
color.TTTT           = ROOT.kOrange+1
color.TTWZ           = ROOT.kOrange+1
color.TTWW           = ROOT.kOrange+2
color.TTZZ           = ROOT.kOrange+3
color.T              = ROOT.kOrange+4
color.TGamma         = ROOT.kGray
color.tW             = ROOT.kCyan+3
color.W              = ROOT.kCyan+1
color.WJets          = ROOT.kCyan+1
color.QCD            = ROOT.kGreen-2 #+3
color.GJets          = ROOT.kGreen+4

color.gen            = ROOT.kOrange
color.had            = ROOT.kAzure-3
color.misID          = ROOT.kAzure+4
color.fakes          = ROOT.kRed-3
color.PU             = ROOT.kGreen+4

color.WG_misID       = ROOT.kAzure-2
color.ZG_misID       = ROOT.kAzure+2
color.DY_misID       = ROOT.kAzure-8
color.VGamma_misID   = ROOT.kAzure-2
color.TTG_had        = ROOT.kOrange-9
color.Other_misID    = ROOT.kOrange+6
color.TT_misID       = ROOT.kRed-7

color.TTG1L          = ROOT.kBlue+2
color.TTG3LBuggy     = ROOT.kRed+1
color.TTG3LPatched   = ROOT.kGreen+3
color.TTG1           = ROOT.kOrange+1
color.TTG2           = ROOT.kGray
color.TTG3           = ROOT.kCyan+1
color.TTG4           = ROOT.kViolet
color.TTG5           = ROOT.kBlue
color.TTG6           = ROOT.kBlack


color.VG1            = ROOT.kBlue+2
color.VG2            = ROOT.kRed+1
color.VG3            = ROOT.kGreen+3
color.VG4            = ROOT.kOrange+1
color.VG5            = ROOT.kBlack
color.VG6            = ROOT.kCyan+1
color.VG7            = ROOT.kViolet
color.VG8            = ROOT.kYellow


color.DY             = ROOT.kOrange-3
color.TTZ            = ROOT.kGreen+4 #ROOT.kBlack
color.TTLep          = ROOT.kAzure+6

