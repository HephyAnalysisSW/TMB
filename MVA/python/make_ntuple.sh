#!/bin/sh
#python make_ntuple.py  --output /eos/vbc/user/$USER/TMB/training-ntuples-v6 --selection singlelep-photon --sample ttG_noFullyHad_fast --config ttG_WG
#python make_ntuple.py  --output /eos/vbc/user/$USER/TMB/training-ntuples-v6 --selection singlelep-photon --sample WGToLNu_fast --config ttG_WG

#python make_ntuple.py  --output /scratch-cbe/users/$USER/TMB/training-ntuples-tttt-v3 --sample TTTT --config tttt_3l 
#python make_ntuple.py  --output /scratch-cbe/users/$USER/TMB/training-ntuples-tttt-v3 --sample TTW  --config tttt_3l 
#python make_ntuple.py  --output /scratch-cbe/users/$USER/TMB/training-ntuples-tttt-v3 --sample TTZ  --config tttt_3l 
#python make_ntuple.py  --output /scratch-cbe/users/$USER/TMB/training-ntuples-tttt-v3 --sample nonprompt_3l --config tttt_3l 
#python make_ntuple.py  --output /scratch-cbe/users/$USER/TMB/training-ntuples-tttt-v3 --sample TTH --config tttt_3l_tth 

python make_ntuple.py  --output /scratch-cbe/users/$USER/TMB/training-ntuples-ttZ --sample TTLep --config ttZ_2l 
python make_ntuple.py  --output /scratch-cbe/users/$USER/TMB/training-ntuples-ttZ --sample DY  --config ttZ_2l 
python make_ntuple.py  --output /scratch-cbe/users/$USER/TMB/training-ntuples-ttZ --sample TTZ  --config ttZ_2l 

#python make_ntuple.py  --output /scratch-cbe/users/$USER/TMB/training-ntuples-ttZ --sample TTW  --config ttZ_dy 
#python make_ntuple.py  --output /scratch-cbe/users/$USER/TMB/training-ntuples-ttZ --sample DY  --config ttZ_dy 
#python make_ntuple.py  --output /scratch-cbe/users/$USER/TMB/training-ntuples-ttZ --sample TTZ  --config ttZ_dy 
##python make_ntuple.py  --output /scratch-cbe/users/$USER/TMB/training-ntuples-tttt-v3 --sample nonprompt_3l --config tttt_3l 

#python make_ntuple.py  --output /eos/vbc/user/$USER/TMB/training-ntuples-tttt-v1 --sample TTTT         --config tttt_2l 
#python make_ntuple.py  --output /eos/vbc/user/$USER/TMB/training-ntuples-tttt-v1 --sample TTLep_bb     --config tttt_2l 
#python make_ntuple.py  --output /eos/vbc/user/$USER/TMB/training-ntuples-tttt-v1 --sample TTLep_cc     --config tttt_2l 
#python make_ntuple.py  --output /eos/vbc/user/$USER/TMB/training-ntuples-tttt-v1 --sample TTLep_other  --config tttt_2l 
#python make_ntuple.py  --output /eos/vbc/user/$USER/TMB/training-ntuples-tttt-v2 --sample TTTT         --config tttt_2l 
#python make_ntuple.py  --output /eos/vbc/user/$USER/TMB/training-ntuples-tttt-v2 --sample TTLep_bb     --config tttt_2l 
#python make_ntuple.py  --output /eos/vbc/user/$USER/TMB/training-ntuples-tttt-v2 --sample TTLep_cc     --config tttt_2l 
#python make_ntuple.py  --output /eos/vbc/user/$USER/TMB/training-ntuples-tttt-v2 --sample TTLep_other  --config tttt_2l 

#python make_ntuple.py  --output /eos/vbc/user/$USER/TMB/training-ntuples-tttt-v1 --selection dilepL --sample_file '$CMSSW_BASE/src/TMB/Samples/python/nanoTuples_RunII_nanoAODv6_dilep_pp.py' --sample TTTT --config tttt_2l 
#python make_ntuple.py  --output /eos/vbc/user/$USER/TMB/training-ntuples-tttt-v1 --selection dilepL --sample_file '$CMSSW_BASE/src/TMB/Samples/python/nanoTuples_RunII_nanoAODv6_dilep_pp.py' --sample TTLep --config tttt_2l

#python $CMSSW_BASE/src/Analysis/MVA/python/make_ntuple.py  --output /scratch-cbe/users/$USER/TMB/training-ntuples-TTZ-flavor --sample WZ  --config_module TMB.MVA.configs --config ttZ_3l_flavor
#python $CMSSW_BASE/src/Analysis/MVA/python/make_ntuple.py  --output /scratch-cbe/users/$USER/TMB/training-ntuples-TTZ-flavor --sample TTZ_ISR  --config_module TMB.MVA.configs --config ttZ_3l_flavor
#python $CMSSW_BASE/src/Analysis/MVA/python/make_ntuple.py  --output /scratch-cbe/users/$USER/TMB/training-ntuples-TTZ-flavor --sample TTZ_nonISR  --config_module TMB.MVA.configs --config ttZ_3l_flavor
 
