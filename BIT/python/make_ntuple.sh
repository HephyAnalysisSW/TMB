#!/bin/sh
python $CMSSW_BASE/src/Analysis/MVA/python/make_ntuple.py  --output /scratch-cbe/users/$USER/BIT/training-ntuples-TTZ-flavor --sample ttZ01j  --config_module TMB.BIT.configs --config ttZ_3l_flavor
python $CMSSW_BASE/src/Analysis/MVA/python/make_ntuple.py  --output /scratch-cbe/users/$USER/BIT/training-ntuples-TTZ-flavor --sample WZ      --config_module TMB.BIT.configs --config WZ_3l_flavor
