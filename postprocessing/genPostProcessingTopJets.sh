#!/bin/sh

#python genPostProcessingTopJets.py --overwrite --addReweights --sample tt1LepHad #SPLIT100
python genPostProcessingTopJets.py --overwrite --addReweights --targetDir v3 --sample tschRefPoint #SPLIT100
python genPostProcessingTopJets.py --overwrite --addReweights --targetDir v3 --sample tschRefPointNoWidthRW #SPLIT100
