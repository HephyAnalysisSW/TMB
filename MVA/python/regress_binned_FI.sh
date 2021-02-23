#!/bin/sh

python regress_binned_FI.py --name train_ctZ_BSM_in_ttG_noFullyHad_v6         --config ttG_WG --FI_branch FI_ctZ_BSM --trainingfile ttG_noFullyHad_fast/singlelep-photon/ttG_noFullyHad_fast.root
python regress_binned_FI.py --name train_ctZ_BSM_in_ttG_noFullyHad_WGToLNu_v6 --config ttG_WG --FI_branch FI_ctZ_BSM --trainingfile ttG_noFullyHad_fast/singlelep-photon/ttG_noFullyHad_fast.root WGToLNu_fast/singlelep-photon/WGToLNu_fast.root

python regress_binned_FI.py --name train_cWWW_BSM_in_WGToLNu         --config ttG_WG --FI_branch FI_cWWW_BSM --trainingfile WGToLNu_fast/singlelep-photon/WGToLNu_fast.root 
python regress_binned_FI.py --name train_cWWW_BSM_in_ttG_noFullyHad_WGToLNu --config ttG_WG --FI_branch FI_cWWW_BSM --trainingfile ttG_noFullyHad_fast/singlelep-photon/ttG_noFullyHad_fast.root WGToLNu_fast/singlelep-photon/WGToLNu_fast.root
