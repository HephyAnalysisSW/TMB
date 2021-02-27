#!/bin/sh

python regress_binned_FI_LSTM.py --name ctZ_LSTM_BSM_in_ttG_noFullyHad_v6         --config ttG_WG --FI_branch FI_ctZ_BSM --trainingfile ttG_noFullyHad_fast/singlelep-photon/ttG_noFullyHad_fast.root
python regress_binned_FI_LSTM.py --name ctZ_LSTM_BSM_in_ttG_noFullyHad_WGToLNu_v6 --config ttG_WG --FI_branch FI_ctZ_BSM --trainingfile ttG_noFullyHad_fast/singlelep-photon/ttG_noFullyHad_fast.root WGToLNu_fast/singlelep-photon/WGToLNu_fast.root
python regress_binned_FI.py      --name ctZ_BSM_in_ttG_noFullyHad_v6              --config ttG_WG --FI_branch FI_ctZ_BSM --trainingfile ttG_noFullyHad_fast/singlelep-photon/ttG_noFullyHad_fast.root
python regress_binned_FI.py      --name ctZ_BSM_in_ttG_noFullyHad_WGToLNu_v6      --config ttG_WG --FI_branch FI_ctZ_BSM --trainingfile ttG_noFullyHad_fast/singlelep-photon/ttG_noFullyHad_fast.root WGToLNu_fast/singlelep-photon/WGToLNu_fast.root

python regress_binned_FI_LSTM.py --name cWWW_LSTM_BSM_in_ttG_noFullyHad_v6         --config ttG_WG --FI_branch FI_cWWW_BSM --trainingfile ttG_noFullyHad_fast/singlelep-photon/ttG_noFullyHad_fast.root
python regress_binned_FI_LSTM.py --name cWWW_LSTM_BSM_in_ttG_noFullyHad_WGToLNu_v6 --config ttG_WG --FI_branch FI_cWWW_BSM --trainingfile ttG_noFullyHad_fast/singlelep-photon/ttG_noFullyHad_fast.root WGToLNu_fast/singlelep-photon/WGToLNu_fast.root
python regress_binned_FI.py      --name cWWW_BSM_in_ttG_noFullyHad_v6              --config ttG_WG --FI_branch FI_cWWW_BSM --trainingfile ttG_noFullyHad_fast/singlelep-photon/ttG_noFullyHad_fast.root
python regress_binned_FI.py      --name cWWW_BSM_in_ttG_noFullyHad_WGToLNu_v6      --config ttG_WG --FI_branch FI_cWWW_BSM --trainingfile ttG_noFullyHad_fast/singlelep-photon/ttG_noFullyHad_fast.root WGToLNu_fast/singlelep-photon/WGToLNu_fast.root

python regress_binned_FI_LSTM.py --name cpDC_LSTM_BSM_in_ttG_noFullyHad_v6         --config ttG_WG --FI_branch FI_cpDC_BSM --trainingfile ttG_noFullyHad_fast/singlelep-photon/ttG_noFullyHad_fast.root
python regress_binned_FI_LSTM.py --name cpDC_LSTM_BSM_in_ttG_noFullyHad_WGToLNu_v6 --config ttG_WG --FI_branch FI_cpDC_BSM --trainingfile ttG_noFullyHad_fast/singlelep-photon/ttG_noFullyHad_fast.root WGToLNu_fast/singlelep-photon/WGToLNu_fast.root
python regress_binned_FI.py      --name cpDC_BSM_in_ttG_noFullyHad_v6              --config ttG_WG --FI_branch FI_cpDC_BSM --trainingfile ttG_noFullyHad_fast/singlelep-photon/ttG_noFullyHad_fast.root
python regress_binned_FI.py      --name cpDC_BSM_in_ttG_noFullyHad_WGToLNu_v6      --config ttG_WG --FI_branch FI_cpDC_BSM --trainingfile ttG_noFullyHad_fast/singlelep-photon/ttG_noFullyHad_fast.root WGToLNu_fast/singlelep-photon/WGToLNu_fast.root

python regress_binned_FI_LSTM.py --name cpWB_LSTM_BSM_in_ttG_noFullyHad_v6         --config ttG_WG --FI_branch FI_cpWB_BSM --trainingfile ttG_noFullyHad_fast/singlelep-photon/ttG_noFullyHad_fast.root
python regress_binned_FI_LSTM.py --name cpWB_LSTM_BSM_in_ttG_noFullyHad_WGToLNu_v6 --config ttG_WG --FI_branch FI_cpWB_BSM --trainingfile ttG_noFullyHad_fast/singlelep-photon/ttG_noFullyHad_fast.root WGToLNu_fast/singlelep-photon/WGToLNu_fast.root
python regress_binned_FI.py      --name cpWB_BSM_in_ttG_noFullyHad_v6              --config ttG_WG --FI_branch FI_cpWB_BSM --trainingfile ttG_noFullyHad_fast/singlelep-photon/ttG_noFullyHad_fast.root
python regress_binned_FI.py      --name cpWB_BSM_in_ttG_noFullyHad_WGToLNu_v6      --config ttG_WG --FI_branch FI_cpWB_BSM --trainingfile ttG_noFullyHad_fast/singlelep-photon/ttG_noFullyHad_fast.root WGToLNu_fast/singlelep-photon/WGToLNu_fast.root

python regress_binned_FI_LSTM.py --name ctG_LSTM_BSM_in_ttG_noFullyHad_v6         --config ttG_WG --FI_branch FI_ctG_BSM --trainingfile ttG_noFullyHad_fast/singlelep-photon/ttG_noFullyHad_fast.root
python regress_binned_FI_LSTM.py --name ctG_LSTM_BSM_in_ttG_noFullyHad_WGToLNu_v6 --config ttG_WG --FI_branch FI_ctG_BSM --trainingfile ttG_noFullyHad_fast/singlelep-photon/ttG_noFullyHad_fast.root WGToLNu_fast/singlelep-photon/WGToLNu_fast.root
python regress_binned_FI.py      --name ctG_BSM_in_ttG_noFullyHad_v6              --config ttG_WG --FI_branch FI_ctG_BSM --trainingfile ttG_noFullyHad_fast/singlelep-photon/ttG_noFullyHad_fast.root
python regress_binned_FI.py      --name ctG_BSM_in_ttG_noFullyHad_WGToLNu_v6      --config ttG_WG --FI_branch FI_ctG_BSM --trainingfile ttG_noFullyHad_fast/singlelep-photon/ttG_noFullyHad_fast.root WGToLNu_fast/singlelep-photon/WGToLNu_fast.root

python regress_binned_FI_LSTM.py --name cpG_LSTM_BSM_in_ttG_noFullyHad_v6         --config ttG_WG --FI_branch FI_cpG_BSM --trainingfile ttG_noFullyHad_fast/singlelep-photon/ttG_noFullyHad_fast.root
python regress_binned_FI_LSTM.py --name cpG_LSTM_BSM_in_ttG_noFullyHad_WGToLNu_v6 --config ttG_WG --FI_branch FI_cpG_BSM --trainingfile ttG_noFullyHad_fast/singlelep-photon/ttG_noFullyHad_fast.root WGToLNu_fast/singlelep-photon/WGToLNu_fast.root
python regress_binned_FI.py      --name cpG_BSM_in_ttG_noFullyHad_v6              --config ttG_WG --FI_branch FI_cpG_BSM --trainingfile ttG_noFullyHad_fast/singlelep-photon/ttG_noFullyHad_fast.root
python regress_binned_FI.py      --name cpG_BSM_in_ttG_noFullyHad_WGToLNu_v6      --config ttG_WG --FI_branch FI_cpG_BSM --trainingfile ttG_noFullyHad_fast/singlelep-photon/ttG_noFullyHad_fast.root WGToLNu_fast/singlelep-photon/WGToLNu_fast.root

python regress_binned_FI_LSTM.py --name cpQ3_LSTM_BSM_in_ttG_noFullyHad_v6         --config ttG_WG --FI_branch FI_cpQ3_BSM --trainingfile ttG_noFullyHad_fast/singlelep-photon/ttG_noFullyHad_fast.root
python regress_binned_FI_LSTM.py --name cpQ3_LSTM_BSM_in_ttG_noFullyHad_WGToLNu_v6 --config ttG_WG --FI_branch FI_cpQ3_BSM --trainingfile ttG_noFullyHad_fast/singlelep-photon/ttG_noFullyHad_fast.root WGToLNu_fast/singlelep-photon/WGToLNu_fast.root
python regress_binned_FI.py      --name cpQ3_BSM_in_ttG_noFullyHad_v6              --config ttG_WG --FI_branch FI_cpQ3_BSM --trainingfile ttG_noFullyHad_fast/singlelep-photon/ttG_noFullyHad_fast.root
python regress_binned_FI.py      --name cpQ3_BSM_in_ttG_noFullyHad_WGToLNu_v6      --config ttG_WG --FI_branch FI_cpQ3_BSM --trainingfile ttG_noFullyHad_fast/singlelep-photon/ttG_noFullyHad_fast.root WGToLNu_fast/singlelep-photon/WGToLNu_fast.root
