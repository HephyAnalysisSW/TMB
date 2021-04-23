#python multiclass.py --name tttt_3l_ttw_ttz_nonprompt_tth_v3 --input_directory /scratch-cbe/users/$USER/TMB/training-ntuples-tttt-v3/MVA-training/ --config tttt_3l_tth --output_directory /mnt/hephy/cms/$USER/TMB/models --add_LSTM 
#python multiclass.py --name tttt_3l_ttw_ttz_nonprompt_v3 --input_directory /scratch-cbe/users/$USER/TMB/training-ntuples-tttt-v3/MVA-training/ --config tttt_3l --output_directory /mnt/hephy/cms/rosmarie.schoefbeck/TMB/models --add_LSTM 
#python multiclass.py --name tttt_2l_test --config tttt_2l 

#python multiclass.py --name ttZ_tt_dy --input_directory /scratch-cbe/users/rosmarie.schoefbeck/TMB/training-ntuples-ttZ/MVA-training/ --config ttZ_2l --output_directory /mnt/hephy/cms/rosmarie.schoefbeck/TMB/models 
#python multiclass.py --name ttZ_tt_dy --input_directory /scratch-cbe/users/rosmarie.schoefbeck/TMB/training-ntuples-ttZ/MVA-training/ --config ttZ_2l --output_directory /mnt/hephy/cms/rosmarie.schoefbeck/TMB/models --add_LSTM 

python multiclass.py --name ttZ_dy --input_directory /scratch-cbe/users/rosmarie.schoefbeck/TMB/training-ntuples-ttZ/MVA-training/ --config ttZ_dy --output_directory /mnt/hephy/cms/rosmarie.schoefbeck/TMB/models 
python multiclass.py --name ttZ_dy --input_directory /scratch-cbe/users/rosmarie.schoefbeck/TMB/training-ntuples-ttZ/MVA-training/ --config ttZ_dy --output_directory /mnt/hephy/cms/rosmarie.schoefbeck/TMB/models --add_LSTM 


