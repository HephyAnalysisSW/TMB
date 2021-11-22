from Analysis.Tools.CutInterpreter import CutInterpreter

mZ              = 91.1876
mT              = 172.5

special_cuts = {
    'VBS':"genVBS_mjj>400&&abs(genVBS_dEta)>2.4"
  }

continous_variables  = [
    ("ptG","genPhoton_pt[0]"),
    ("Z1pt", "genZ_pt[0]"), ("Z2pt", "genZ_pt[1]")
    ]

discrete_variables  = [
    ("nGenJet", "ngenJet"),
    ("nGenBTag", "Sum$(genJet_pt>30&&abs(genJet_eta)<2.4&&genJet_matchBParton)"),
    ("nZ", "ngenZ"),
    ("nW", "ngenW"),
    ("nG", "ngenPhoton"),
    ]

cutInterpreter = CutInterpreter( continous_variables, discrete_variables, special_cuts)

if __name__ == "__main__":
    print cutInterpreter.cutString("nZ1p-nW1p")


