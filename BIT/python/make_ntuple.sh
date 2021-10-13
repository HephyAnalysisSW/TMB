#!/bin/sh
#python $CMSSW_BASE/src/Analysis/MVA/python/make_ntuple.py  --output /scratch-cbe/users/$USER/BIT/training-ntuples-TTZ-flavor --sample ttZ01j  --config_module TMB.BIT.configs --config ttZ_3l_flavor
#python $CMSSW_BASE/src/Analysis/MVA/python/make_ntuple.py  --output /scratch-cbe/users/$USER/BIT/training-ntuples-TTZ-flavor --sample WZ      --config_module TMB.BIT.configs --config WZ_3l_flavor
#python $CMSSW_BASE/src/Analysis/MVA/python/make_ntuple.py  --output /scratch-cbe/users/$USER/BIT/training-ntuples-ttG_ctZ --sample ttg      --config_module TMB.BIT.configs --config ttG_ctZ

#python $CMSSW_BASE/src/Analysis/MVA/python/make_ntuple.py  --output /scratch-cbe/users/$USER/BIT/training-ntuples-TTZ-flavor-delphes --sample ttZ01j --config_module TMB.BIT.configs --config ttZ_3l_flavor_delphes
#python $CMSSW_BASE/src/Analysis/MVA/python/make_ntuple.py  --output /scratch-cbe/users/$USER/BIT/training-ntuples-TTZ-flavor-delphes --sample WZTo3L1Nu --config_module TMB.BIT.configs --config WZ_3l_flavor_delphes
#python $CMSSW_BASE/src/Analysis/MVA/python/make_ntuple.py  --output /scratch-cbe/users/$USER/BIT/training-ntuples-TTZ-flavor-delphes --sample ZZ --config_module TMB.BIT.configs --config ZZ_3l_flavor_delphes

python $CMSSW_BASE/src/Analysis/MVA/python/make_ntuple.py  --output /scratch-cbe/users/$USER/BIT/training-ntuples-TTZ-flavor-delphes --sample ttZ01j --config_module TMB.BIT.configs --config ttZ_WZ_3l_flavor_delphes
python $CMSSW_BASE/src/Analysis/MVA/python/make_ntuple.py  --output /scratch-cbe/users/$USER/BIT/training-ntuples-TTZ-flavor-delphes --sample WZTo3L1Nu --config_module TMB.BIT.configs --config ttZ_WZ_3l_flavor_delphes
