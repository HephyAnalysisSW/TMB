#python multiclass.py --name tttt_3l_ttw_ttz_nonprompt_tth_v3 --input_directory /scratch-cbe/users/$USER/TMB/training-ntuples-tttt-v3/MVA-training/ --config tttt_3l_tth --output_directory /mnt/hephy/cms/$USER/TMB/models --add_LSTM 
#python multiclass.py --name tttt_3l_ttw_ttz_nonprompt_v3 --input_directory /scratch-cbe/users/$USER/TMB/training-ntuples-tttt-v3/MVA-training/ --config tttt_3l --output_directory /mnt/hephy/cms/$USER/TMB/models --add_LSTM 
#python multiclass.py --name tttt_2l_test --config tttt_2l 

#python multiclass.py --name ttZ_tt_dy --input_directory /scratch-cbe/users/$USER/TMB/training-ntuples-ttZ/MVA-training/ --config ttZ_2l --output_directory /mnt/hephy/cms/$USER/TMB/models 
#python multiclass.py --name ttZ_tt_dy --input_directory /scratch-cbe/users/$USER/TMB/training-ntuples-ttZ/MVA-training/ --config ttZ_2l --output_directory /mnt/hephy/cms/$USER/TMB/models --add_LSTM 

#python multiclass.py --name ttZ_dy --input_directory /scratch-cbe/users/rosmarie.schoefbeck/TMB/training-ntuples-ttZ/MVA-training/ --config ttZ_dy --output_directory /mnt/hephy/cms/$USER/TMB/models --small 
#python multiclass.py --name ttZ_dy --input_directory /scratch-cbe/users/rosmarie.schoefbeck/TMB/training-ntuples-ttZ/MVA-training/ --config ttZ_dy --output_directory /mnt/hephy/cms/$USER/TMB/models --add_LSTM --small

#python multiclass.py --name ttZ_dy --input_directory /scratch-cbe/users/rosmarie.schoefbeck/TMB/training-ntuples-ttZ/MVA-training/ --config ttZ_dy --output_directory /mnt/hephy/cms/$USER/TMB/models  
#python multiclass.py --name ttZ_dy --input_directory /scratch-cbe/users/rosmarie.schoefbeck/TMB/training-ntuples-ttZ/MVA-training/ --config ttZ_dy --output_directory /mnt/hephy/cms/$USER/TMB/models --add_LSTM

#python multiclass_generator.py --name ttZ_dy_50 --input_directory /scratch-cbe/users/rosmarie.schoefbeck/TMB/training-ntuples-ttZ/MVA-training/ --config ttZ_dy --output_directory /mnt/hephy/cms/$USER/TMB/models --small
#python multiclass_generator.py --name ttZ_dy_50 --input_directory /scratch-cbe/users/rosmarie.schoefbeck/TMB/training-ntuples-ttZ/MVA-training/ --config ttZ_dy --output_directory /mnt/hephy/cms/$USER/TMB/models --add_LSTM --small

python $CMSSW_BASE/src/Analysis/MVA/python/multiclass_generator.py --name ttZ_dy_100 --input_directory /scratch-cbe/users/rosmarie.schoefbeck/TMB/training-ntuples-ttZ/MVA-training/ --config_module TMB.MVA.configs --config ttZ_dy --output_directory /mnt/hephy/cms/$USER/TMB/models 
python $CMSSW_BASE/src/Analysis/MVA/python/multiclass_generator.py --name ttZ_dy_100 --input_directory /scratch-cbe/users/rosmarie.schoefbeck/TMB/training-ntuples-ttZ/MVA-training/ --config_module TMB.MVA.configs --config ttZ_dy --output_directory /mnt/hephy/cms/$USER/TMB/models --add_LSTM 

#python multiclass_generator.py --name ttZ_dy_100_lstm7 --input_directory /scratch-cbe/users/rosmarie.schoefbeck/TMB/training-ntuples-ttZ/MVA-training/ --config ttZ_dy --output_directory /mnt/hephy/cms/$USER/TMB/models --small
#python multiclass_generator.py --name ttZ_dy_100_lstm7_2 --input_directory /scratch-cbe/users/rosmarie.schoefbeck/TMB/training-ntuples-ttZ/MVA-training/ --config ttZ_dy --output_directory /mnt/hephy/cms/$USER/TMB/models --add_LSTM --small
#python multiclass_generator.py --name ttZ_tt_dy --input_directory /scratch-cbe/users/rosmarie.schoefbeck/TMB/training-ntuples-ttZ/MVA-training/ --config ttZ_2l --output_directory /mnt/hephy/cms/$USER/TMB/models 
#python multiclass_generator.py --name ttZ_tt_dy --input_directory /scratch-cbe/users/rosmarie.schoefbeck/TMB/training-ntuples-ttZ/MVA-training/ --config ttZ_2l --output_directory /mnt/hephy/cms/$USER/TMB/models --add_LSTM
#
#python multiclass.py --name ttZ_dy --input_directory /scratch-cbe/users/rosmarie.schoefbeck/TMB/training-ntuples-ttZ/MVA-training/ --config ttZ_dy --output_directory /mnt/hephy/cms/$USER/TMB/models  
#python multiclass.py --name ttZ_dy --input_directory /scratch-cbe/users/rosmarie.schoefbeck/TMB/training-ntuples-ttZ/MVA-training/ --config ttZ_dy --output_directory /mnt/hephy/cms/$USER/TMB/models --add_LSTM

#python $CMSSW_BASE/src/Analysis/MVA/python/multiclass.py --name ttZ_3l_flavor --input_directory /scratch-cbe/users/$USER/TMB/training-ntuples-TTZ-flavor/MVA-training/ --config_module TMB.MVA.configs --config ttZ_3l_flavor --output_directory /mnt/hephy/cms/$USER/TMB
#python $CMSSW_BASE/src/Analysis/MVA/python/multiclass.py --name ttZ_3l_flavor --add_LSTM --input_directory /scratch-cbe/users/$USER/TMB/training-ntuples-TTZ-flavor/MVA-training/ --config_module TMB.MVA.configs --config ttZ_3l_flavor --output_directory /mnt/hephy/cms/$USER/TMB
