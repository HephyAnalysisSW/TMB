#!/bin/sh

python analysis-singlelep2.py --WC ctZ --plot_directory analysis-v12_1
python analysis-singlelep2.py --WC cWWW --plot_directory analysis-v12_1

python analysis-singlelep21.py --WC ctZ --plot_directory analysis-v12_2
python analysis-singlelep21.py --WC cWWW --plot_directory analysis-v12_2

#python analysis-singlelep4.py --WC ctZ --plot_directory analysis-v11
#python analysis-singlelep4.py --WC cWWW --plot_directory analysis-v11

#python analysis-singlelep2.py --WC ctZ --plot_directory analysis-v9_2
#python analysis-singlelep2.py --WC cWWW --plot_directory analysis-v9_2
