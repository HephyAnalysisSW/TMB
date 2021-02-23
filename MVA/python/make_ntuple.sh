#!/bin/sh
python make_ntuple.py  --output /eos/vbc/user/$USER/TMB/training-ntuples-v6 --selection singlelep-photon --sample ttG_noFullyHad_fast --config ttG_WG
python make_ntuple.py  --output /eos/vbc/user/$USER/TMB/training-ntuples-v6 --selection singlelep-photon --sample WGToLNu_fast --config ttG_WG
