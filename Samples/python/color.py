import ROOT

from TMB.Tools.helpers import singleton as singleton

@singleton
class color():
  pass

color.data           = ROOT.kBlack
color.ZZ             = ROOT.kBlue + 1
color.ZZjj           = ROOT.kBlue - 9

color.WZ             = ROOT.kRed  + 1
color.WZjj           = ROOT.kRed  - 4

color.WW             = ROOT.kYellow + 1
color.WWjj_OS        = ROOT.kYellow + 2 
color.WWjj_SS        = ROOT.kYellow + 3

color.WWW            = ROOT.kMagenta
color.WWZ            = ROOT.kMagenta + 1
color.WZZ            = ROOT.kMagenta - 4
color.ZZZ            = ROOT.kMagenta - 9
color.WA             = ROOT.kGreen + 1
color.WAjj           = ROOT.kGreen - 7

color.ZA             = ROOT.kCyan  + 1
color.ZAjj           = ROOT.kCyan  + 2

color.ttg_noFullyHad = ROOT.kOrange 
color.ttW01j         = ROOT.kOrange + 10
color.ttZ01j         = ROOT.kOrange + 7
