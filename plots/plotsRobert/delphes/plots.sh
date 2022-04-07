#python plots.py --signal ZH --selection dilep-ZHJet-onZ-onH-ptZ300  --combinatoricalBTags --plot_directory delphes-combWeights --DYsample DYBBJets

#python plots.py --signal ZH --selection dilep-ZHJet-onZ-onH-ptZ300  --combinatoricalBTags --plot_directory delphes-combWeights --DYsample DYJets_HT 
#python plots.py --signal ZH --selection dilep-ZHJet-onZ-onH-ptZ300  --combinatoricalBTags --plot_directory delphes-combWeights --DYsample DYJets

python plots.py --signal ZH --selection dilep-ZHJet-onZ-onH-ptZ300  --plot_directory delphes --DYsample all

#python plots.py --signal ZH --selection dilep-ZHJet-onZ-onH-ptZ300  --plot_directory delphes-combWeights --DYsample DYBBJets
#python plots.py --signal ZH --selection dilep-ZHJet-onZ-onH-ptZ300  --plot_directory delphes-combWeights --DYsample DYJets_HT 
#python plots.py --signal ZH --selection dilep-ZHJet-onZ-onH-ptZ300  --plot_directory delphes-combWeights --DYsample DYJets

#python plots.py --signal ZH_nlo --selection dilep-ZHJet-onZ-onH-ptZ300  --combinatoricalBTags --plot_directory delphes-combWeights --DYsample DYBBJets
#python plots.py --signal ZH_nlo --selection dilep-ZHJet-onZ-onH-ptZ300  --combinatoricalBTags --plot_directory delphes-combWeights --DYsample DYJets_HT 
#python plots.py --signal ZH_nlo --selection dilep-ZHJet-onZ-onH-ptZ300  --combinatoricalBTags --plot_directory delphes-combWeights --DYsample DYJets
#python plots.py --signal ZH_nlo --selection dilep-ZHJet-onZ-onH-ptZ300  --plot_directory delphes-combWeights --DYsample DYBBJets
#python plots.py --signal ZH_nlo --selection dilep-ZHJet-onZ-onH-ptZ300  --plot_directory delphes-combWeights --DYsample DYJets_HT 
#python plots.py --signal ZH_nlo --selection dilep-ZHJet-onZ-onH-ptZ300  --plot_directory delphes-combWeights --DYsample DYJets

python plots.py --signal WH --selection singlelep-WHJet-realW-onH-ptW300 --plot_directory delphes
#python plots.py --signal WH_nlo --selection singlelep-WHJet-realW-onH-ptW300 --plot_directory delphes

#python plots.py --signal ZH --selection dilep-ZHJet-onZ-onH-ptZ300       --no_bkgs --plot_directory delphes-no_bkgs
#python plots.py --signal WH --selection singlelep-WHJet-realW-onH-ptW300 --no_bkgs --plot_directory delphes-no_bkgs

#python plots.py --signal ZH --sign_reweight --selection dilep-ZHJet-onZ-onH-ptZ300       --no_bkgs --plot_directory delphes-no_bkgs
#python plots.py --signal WH --sign_reweight --selection singlelep-WHJet-realW-onH-ptW300 --no_bkgs --plot_directory delphes-no_bkgs
