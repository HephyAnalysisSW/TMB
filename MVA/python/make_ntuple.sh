#!/bin/sh
#python make_ntuple.py  --output /eos/vbc/user/$USER/TMB/training-ntuples-v6 --selection singlelep-photon --sample ttG_noFullyHad_fast --config ttG_WG
#python make_ntuple.py  --output /eos/vbc/user/$USER/TMB/training-ntuples-v6 --selection singlelep-photon --sample WGToLNu_fast --config ttG_WG


python make_ntuple.py  --output /eos/vbc/user/$USER/TMB/training-ntuples-tttt-v1 --selection trilepVL --sample_file '$CMSSW_BASE/src/tWZ/samples/python/nanoTuples_RunII_nanoAODv6_private_postProcessed.py' --sample TTTT --config tttt_3l 
python make_ntuple.py  --output /eos/vbc/user/$USER/TMB/training-ntuples-tttt-v1 --selection trilepVL --sample_file '$CMSSW_BASE/src/tWZ/samples/python/nanoTuples_RunII_nanoAODv6_private_postProcessed.py' --sample TTW --config tttt_3l 
python make_ntuple.py  --output /eos/vbc/user/$USER/TMB/training-ntuples-tttt-v1 --selection trilepVL --sample_file '$CMSSW_BASE/src/tWZ/samples/python/nanoTuples_RunII_nanoAODv6_private_postProcessed.py' --sample TTZ --config tttt_3l 
