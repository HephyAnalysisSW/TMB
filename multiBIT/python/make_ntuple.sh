#!/bin/sh

#python $CMSSW_BASE/src/Analysis/MVA/python/make_ntuple.py  --output /groups/hephy/cms/$USER/BIT/training-ntuple-ZH --sample ZH_nlo --config_module TMB.BIT.configs --config ZH_nlo_delphes
#python $CMSSW_BASE/src/Analysis/MVA/python/make_ntuple.py  --output /groups/hephy/cms/$USER/BIT/training-ntuple-WH --sample WH_nlo --config_module TMB.BIT.configs --config WH_nlo_delphes

python $CMSSW_BASE/src/Analysis/MVA/python/make_ntuple.py  --output /groups/hephy/cms/$USER/multiBIT/training-ntuple-WH --sample WH --config_module TMB.multiBIT.configs  --config WH_delphes 
