python new-limit.py --signal WH --selection singlelep-WHJet-realW-onH-ptW300 --plot_directory delphes --nBinsTestStat 10 --nBins 100 --WCs cHj3 -.15 .05 cHW    -0.75 0.75 
python new-limit.py --signal WH --selection singlelep-WHJet-realW-onH-ptW300 --plot_directory delphes --nBinsTestStat 10 --nBins 100 --WCs cHj3 -.12 .03 cHWtil -0.75 0.75 
python new-limit.py --signal WH --selection singlelep-WHJet-realW-onH-ptW300 --plot_directory delphes --nBinsTestStat 10 --nBins 100 --WCs cHW -0.75 0.35  cHWtil  -0.45 0.45 

python new-limit.py --signal ZH --selection dilep-ZHJet-onZ-onH-ptZ300 --plot_directory delphes --DYsample all --nBinsTestStat 10 --nBins 100 --WCs cHj3 -.2 .05 cHW    -1 1
python new-limit.py --signal ZH --selection dilep-ZHJet-onZ-onH-ptZ300 --plot_directory delphes --DYsample all --nBinsTestStat 10 --nBins 100 --WCs cHj3 -.2 .05 cHWtil -1 1
python new-limit.py --signal ZH --selection dilep-ZHJet-onZ-onH-ptZ300 --plot_directory delphes --DYsample all --nBinsTestStat 10 --nBins 100 --WCs cHW -1 0.45  cHWtil  -0.8 0.8 
