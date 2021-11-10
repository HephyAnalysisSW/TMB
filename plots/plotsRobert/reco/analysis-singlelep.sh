#!/bin/sh

python analysis-singlelep.py --plot_directory clipScore --WC ctZ
python analysis-singlelep.py --plot_directory clipScore --WC cWWW
python analysis-singlelep.py --plot_directory clipScore --WC ctZ --selection singlelep-photon-btag1p-njet1p

#python analysis-singlelep.py --plot_directory clipScore --WC cWWW --selection singlelep-photon-ptG150To200
#python analysis-singlelep.py --plot_directory clipScore --WC cWWW --selection singlelep-photon-ptG200To300
#python analysis-singlelep.py --plot_directory clipScore --WC cWWW --selection singlelep-photon-ptG300To500
#python analysis-singlelep.py --plot_directory clipScore --WC cWWW --selection singlelep-photon-ptG500To800
#python analysis-singlelep.py --plot_directory clipScore --WC cWWW --selection singlelep-photon-ptG800
