#python new-limit.py --signal WH --selection singlelep-WHJet-realW-onH-ptW300 --plot_directory delphes-ptV200-training --nBinsTestStat 10 --nBins 80 --WCs cHj3 -.18 .05  cHW    -0.9 0.9
#python new-limit.py --signal WH --selection singlelep-WHJet-realW-onH-ptW300 --plot_directory delphes-ptV200-training --nBinsTestStat 10 --nBins 80 --WCs cHj3 -.18 .05  cHWtil -0.9 0.9
#python new-limit.py --signal WH --selection singlelep-WHJet-realW-onH-ptW300 --plot_directory delphes-ptV200-training --nBinsTestStat 10 --nBins 80 --WCs cHW -0.75 0.35  cHWtil  -0.6 0.6 

#python new-limit.py --signal ZH --selection dilep-ZHJet-onZ-onH-ptZ300 --plot_directory delphes-v3 --DYsample all --nBinsTestStat 30 --nBins 100 --WCs cHj3 -.3 .15 cHW    -1.3 1.3
#python new-limit.py --signal ZH --selection dilep-ZHJet-onZ-onH-ptZ300 --plot_directory delphes-v3 --DYsample all --nBinsTestStat 30 --nBins 100 --WCs cHj3 -.3 .15 cHWtil -1.3 1.3
#python new-limit.py --signal ZH --selection dilep-ZHJet-onZ-onH-ptZ300 --plot_directory delphes-v3 --DYsample all --nBinsTestStat 30 --nBins 100 --WCs cHW -1.2 0.65  cHWtil  -1 1 

python new-limit.py --signal ZH --no_bkgs --selection dilep-ZHJet-onZ-onH-ptZ300 --plot_directory delphes-v3 --DYsample all --nBinsTestStat 30 --nBins 100 --WCs cHj3 -.3 .15 cHW    -1.3 1.3
python new-limit.py --signal ZH --no_bkgs --selection dilep-ZHJet-onZ-onH-ptZ300 --plot_directory delphes-v3 --DYsample all --nBinsTestStat 30 --nBins 100 --WCs cHj3 -.3 .15 cHWtil -1.3 1.3
python new-limit.py --signal ZH --no_bkgs --selection dilep-ZHJet-onZ-onH-ptZ300 --plot_directory delphes-v3 --DYsample all --nBinsTestStat 30 --nBins 100 --WCs cHW -1.2 0.65  cHWtil  -1 1 
