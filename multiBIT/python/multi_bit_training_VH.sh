#python training.py --overwrite  --input_directory /groups/hephy/cms/robert.schoefbeck/multiBIT/training-ntuple-WH/MVA-training/WH_delphes --config WH_delphes  --name v1

python multi_bit_training_VH.py --input_directory /groups/hephy/cms/robert.schoefbeck/multiBIT/training-ntuple-WH/MVA-training/WH_delphes --config WH_delphes  --name v1 --coefficients cHW cHWtil cHj3 cHu --debug 
python multi_bit_training_VH.py --input_directory /groups/hephy/cms/robert.schoefbeck/multiBIT/training-ntuple-WH/MVA-training/WH_delphes --config WH_delphes  --name v2 --coefficients cHDD cHW cHWtil cHbox cHj3 cbHRe --debug
