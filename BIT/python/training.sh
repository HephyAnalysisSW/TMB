#python training.py --config ttZ_3l_flavor & 
#python training.py --config ttZ_3l_flavor --name MU0.05 --max_uncertainty 0.05 &
#python training.py --config ttZ_3l_flavor --name MU0.1 --max_uncertainty 0.1 &
#python training.py --config ttZ_3l_flavor --name MU0.2  --max_uncertainty 0.2 &
#python training.py --config ttZ_3l_flavor --name MU0.5  --max_uncertainty 0.5 &
#python training.py --config ttZ_3l_flavor --name MU1.0  --max_uncertainty 1.0 &

#python training.py --overwrite --input_directory /scratch-cbe/users/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor/MVA-training/ttZ_3l_flavor --config ttZ_3l_flavor_ptZ --name UNC2 --max_uncertainty 2
#python training.py --overwrite --input_directory /scratch-cbe/users/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor/MVA-training/ttZ_3l_flavor --config ttZ_3l_flavor_ptZ --name UNC1 --max_uncertainty 1
#python training.py --overwrite --input_directory /scratch-cbe/users/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor/MVA-training/ttZ_3l_flavor --config ttZ_3l_flavor_ptZ --name UNC0p5 --max_uncertainty 0.5
#python training.py --overwrite --input_directory /scratch-cbe/users/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor/MVA-training/ttZ_3l_flavor --config ttZ_3l_flavor_ptZ --name UNC0p4 --max_uncertainty 0.4
#python training.py --overwrite --input_directory /scratch-cbe/users/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor/MVA-training/ttZ_3l_flavor --config ttZ_3l_flavor_ptZ --name UNC0p3 --max_uncertainty 0.3
#python training.py --overwrite --input_directory /scratch-cbe/users/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor/MVA-training/ttZ_3l_flavor --config ttZ_3l_flavor_ptZ --name UNC0p2 --max_uncertainty 0.2
#python training.py --overwrite --input_directory /scratch-cbe/users/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor/MVA-training/ttZ_3l_flavor --config ttZ_3l_flavor_ptZ --name UNC0p1 --max_uncertainty 0.1
#python training.py --overwrite --input_directory /scratch-cbe/users/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor/MVA-training/ttZ_3l_flavor --config ttZ_3l_flavor_ptZ --name UNC0p05 --max_uncertainty 0.05
#python training.py --overwrite --input_directory /scratch-cbe/users/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor/MVA-training/ttZ_3l_flavor --config ttZ_3l_flavor_ptZ --name UNC0p04 --max_uncertainty 0.04
#python training.py --overwrite --input_directory /scratch-cbe/users/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor/MVA-training/ttZ_3l_flavor --config ttZ_3l_flavor_ptZ --name UNC0p03 --max_uncertainty 0.03
#python training.py --overwrite --input_directory /scratch-cbe/users/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor/MVA-training/ttZ_3l_flavor --config ttZ_3l_flavor_ptZ --name UNC0p02 --max_uncertainty 0.02
#python training.py --overwrite --input_directory /scratch-cbe/users/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor/MVA-training/ttZ_3l_flavor --config ttZ_3l_flavor_ptZ --name UNC0p01 --max_uncertainty 0.01 
#python training.py --overwrite --input_directory /scratch-cbe/users/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor/MVA-training/ttZ_3l_flavor --config ttZ_3l_flavor_ptZ --name default

#python training.py --input_directory /scratch-cbe/users/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor-delphes/MVA-training/ttZ_3l_flavor_delphes --config ttZ_3l_flavor_delphes_ptZ --name default
#python training.py --input_directory /scratch-cbe/users/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor-delphes/MVA-training/ttZ_3l_flavor_delphes --config ttZ_3l_flavor_delphes_ptZ  --name BF0p6 --bagging_fraction 0.6
#python training.py --input_directory /scratch-cbe/users/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor-delphes/MVA-training/ttZ_3l_flavor_delphes --config ttZ_3l_flavor_delphes_ptZ  --name BF0p6_MNS10 --bagging_fraction 0.6 --max_n_split 10
#python training.py --input_directory /scratch-cbe/users/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor-delphes/MVA-training/ttZ_3l_flavor_delphes --config ttZ_3l_flavor_delphes_ptZ  --name MNS10  --max_n_split 10
 
#python training.py --input_directory /scratch-cbe/users/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor-delphes/MVA-training/WZ_3l_flavor_delphes --config WZ_3l_flavor_delphes --name default
#python training.py --input_directory /scratch-cbe/users/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor-delphes/MVA-training/ZZ_3l_flavor_delphes --config ZZ_3l_flavor_delphes --name default

#python training.py --input_directory /scratch-cbe/users/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor-delphes/MVA-training/ZZ_3l_flavor_delphes --config ZZ_3l_flavor_delphes --name default --derivative  cHq1Re11
#python training.py --input_directory /scratch-cbe/users/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor-delphes/MVA-training/ZZ_3l_flavor_delphes --config ZZ_3l_flavor_delphes --name default --derivative  cHq1Re22
#python training.py --input_directory /scratch-cbe/users/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor-delphes/MVA-training/ZZ_3l_flavor_delphes --config ZZ_3l_flavor_delphes --name default --derivative  cHq1Re11 cHq1Re11
#python training.py --input_directory /scratch-cbe/users/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor-delphes/MVA-training/ZZ_3l_flavor_delphes --config ZZ_3l_flavor_delphes --name default --derivative  cHq1Re22 cHq1Re22
#python training.py --input_directory /scratch-cbe/users/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor-delphes/MVA-training/ZZ_3l_flavor_delphes --config ZZ_3l_flavor_delphes --name default --derivative  cHq3Re11
#python training.py --input_directory /scratch-cbe/users/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor-delphes/MVA-training/ZZ_3l_flavor_delphes --config ZZ_3l_flavor_delphes --name default --derivative  cHq3Re22
#python training.py --input_directory /scratch-cbe/users/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor-delphes/MVA-training/ZZ_3l_flavor_delphes --config ZZ_3l_flavor_delphes --name default --derivative  cHq3Re11 cHq3Re11
#python training.py --input_directory /scratch-cbe/users/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor-delphes/MVA-training/ZZ_3l_flavor_delphes --config ZZ_3l_flavor_delphes --name default --derivative  cHq3Re22 cHq3Re22
#python training.py --input_directory /scratch-cbe/users/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor-delphes/MVA-training/ZZ_3l_flavor_delphes --config ZZ_3l_flavor_delphes --name default --derivative  cHuRe11
#python training.py --input_directory /scratch-cbe/users/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor-delphes/MVA-training/ZZ_3l_flavor_delphes --config ZZ_3l_flavor_delphes --name default --derivative  cHuRe22
#python training.py --input_directory /scratch-cbe/users/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor-delphes/MVA-training/ZZ_3l_flavor_delphes --config ZZ_3l_flavor_delphes --name default --derivative    cHuRe11 cHuRe11
#python training.py --input_directory /scratch-cbe/users/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor-delphes/MVA-training/ZZ_3l_flavor_delphes --config ZZ_3l_flavor_delphes --name default --derivative    cHuRe22 cHuRe22
#python training.py --input_directory /scratch-cbe/users/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor-delphes/MVA-training/ZZ_3l_flavor_delphes --config ZZ_3l_flavor_delphes --name default --derivative  cHdRe11
#python training.py --input_directory /scratch-cbe/users/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor-delphes/MVA-training/ZZ_3l_flavor_delphes --config ZZ_3l_flavor_delphes --name default --derivative  cHdRe22
#python training.py --input_directory /scratch-cbe/users/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor-delphes/MVA-training/ZZ_3l_flavor_delphes --config ZZ_3l_flavor_delphes --name default --derivative    cHdRe11 cHdRe11
#python training.py --input_directory /scratch-cbe/users/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor-delphes/MVA-training/ZZ_3l_flavor_delphes --config ZZ_3l_flavor_delphes --name default --derivative    cHdRe22 cHdRe22
#python training.py --input_directory /scratch-cbe/users/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor-delphes/MVA-training/ZZ_3l_flavor_delphes --config ZZ_3l_flavor_delphes --name default --derivative  cHudRe11
#python training.py --input_directory /scratch-cbe/users/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor-delphes/MVA-training/ZZ_3l_flavor_delphes --config ZZ_3l_flavor_delphes --name default --derivative  cHudRe22
#python training.py --input_directory /scratch-cbe/users/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor-delphes/MVA-training/ZZ_3l_flavor_delphes --config ZZ_3l_flavor_delphes --name default --derivative  cHudRe11 cHudRe11
#python training.py --input_directory /scratch-cbe/users/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor-delphes/MVA-training/ZZ_3l_flavor_delphes --config ZZ_3l_flavor_delphes --name default --derivative  cHudRe22 cHudRe22


#python training.py --overwrite  --input_directory /scratch-cbe/users/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor-delphes/MVA-training/ttZ_3l_flavor_delphes --config ttZ_3l_flavor_delphes  --name default --derivative  cHq1Re33 
#python training.py --overwrite  --input_directory /scratch-cbe/users/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor-delphes/MVA-training/ttZ_3l_flavor_delphes --config ttZ_3l_flavor_delphes  --name default --derivative  cHq1Re33 cHq1Re33
#python training.py --overwrite  --input_directory /scratch-cbe/users/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor-delphes/MVA-training/ttZ_3l_flavor_delphes --config ttZ_3l_flavor_delphes  --name default --derivative  cHq3Re33 
#python training.py --overwrite  --input_directory /scratch-cbe/users/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor-delphes/MVA-training/ttZ_3l_flavor_delphes --config ttZ_3l_flavor_delphes  --name default --derivative  cHq3Re33 cHq3Re33
#python training.py --overwrite  --input_directory /scratch-cbe/users/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor-delphes/MVA-training/ttZ_3l_flavor_delphes --config ttZ_3l_flavor_delphes  --name default --derivative  cHuRe33 
#python training.py --overwrite  --input_directory /scratch-cbe/users/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor-delphes/MVA-training/ttZ_3l_flavor_delphes --config ttZ_3l_flavor_delphes  --name default --derivative    cHuRe33 cHuRe33  
#python training.py --overwrite  --input_directory /scratch-cbe/users/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor-delphes/MVA-training/ttZ_3l_flavor_delphes --config ttZ_3l_flavor_delphes  --name default --derivative  cHdRe33 
#python training.py --overwrite  --input_directory /scratch-cbe/users/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor-delphes/MVA-training/ttZ_3l_flavor_delphes --config ttZ_3l_flavor_delphes  --name default --derivative    cHdRe33 cHdRe33  
#python training.py --overwrite  --input_directory /scratch-cbe/users/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor-delphes/MVA-training/ttZ_3l_flavor_delphes --config ttZ_3l_flavor_delphes  --name default --derivative  cHudRe33 
#python training.py --overwrite  --input_directory /scratch-cbe/users/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor-delphes/MVA-training/ttZ_3l_flavor_delphes --config ttZ_3l_flavor_delphes  --name default --derivative  cHudRe33 cHudRe33

#python training.py --input_directory /scratch-cbe/users/robert.schoefbeck/BIT/training-ntuples-TTZ-delphes/MVA-training/ttZ_3l_delphes --config ttZ_3l_delphes --lumi_norm --name gss 

#python training.py  --overwrite  --debug --input_directory /scratch-cbe/users/robert.schoefbeck/BIT/training-ntuples-ttG_WG/MVA-training/ttG_WG --config ttG_WG  --name capPtG --derivative  cWWW 
#python training.py --overwrite  --debug --input_directory /scratch-cbe/users/robert.schoefbeck/BIT/training-ntuples-ttG_WG/MVA-training/ttG_WG --config ttG_WG  --name capPtG --derivative  cWWW cWWW
#python training.py --overwrite  --debug --input_directory /scratch-cbe/users/robert.schoefbeck/BIT/training-ntuples-ttG_WG/MVA-training/ttG_WG --config ttG_WG  --name capPtG --derivative  ctZ
#python training.py  --overwrite  --debug --input_directory /scratch-cbe/users/robert.schoefbeck/BIT/training-ntuples-ttG_WG/MVA-training/ttG_WG --config ttG_WG  --name capPtG --derivative  ctZ ctZ 

#python training.py  --overwrite  --debug --input_directory /scratch-cbe/users/robert.schoefbeck/BIT/training-ntuples-ttG_WG/MVA-training/ttG_WG --config ttG_WG  --name capScore_100.0 --derivative  cWWW --max_local_score 100.0
#python training.py --overwrite  --debug --input_directory /scratch-cbe/users/robert.schoefbeck/BIT/training-ntuples-ttG_WG/MVA-training/ttG_WG --config ttG_WG  --name capScore_100.0 --derivative  cWWW cWWW --max_local_score 100.0
#python training.py --overwrite  --debug --input_directory /scratch-cbe/users/robert.schoefbeck/BIT/training-ntuples-ttG_WG/MVA-training/ttG_WG --config ttG_WG  --name capScore_100.0 --derivative  ctZ --max_local_score 100.0
#python training.py  --overwrite  --debug --input_directory /scratch-cbe/users/robert.schoefbeck/BIT/training-ntuples-ttG_WG/MVA-training/ttG_WG --config ttG_WG  --name capScore_100.0 --derivative  ctZ ctZ --max_local_score 100.0

#python training.py  --overwrite  --debug --input_directory /scratch-cbe/users/robert.schoefbeck/BIT/training-ntuples-ttG_WG/MVA-training/ttG_WG --config ttG_WG  --name capScore_10.0 --derivative  cWWW --max_local_score 10.0
#python training.py --overwrite  --debug --input_directory /scratch-cbe/users/robert.schoefbeck/BIT/training-ntuples-ttG_WG/MVA-training/ttG_WG --config ttG_WG  --name capScore_10.0 --derivative  cWWW cWWW --max_local_score 10.0
#python training.py --overwrite  --debug --input_directory /scratch-cbe/users/robert.schoefbeck/BIT/training-ntuples-ttG_WG/MVA-training/ttG_WG --config ttG_WG  --name capScore_10.0 --derivative  ctZ --max_local_score 10.0
#python training.py  --overwrite  --debug --input_directory /scratch-cbe/users/robert.schoefbeck/BIT/training-ntuples-ttG_WG/MVA-training/ttG_WG --config ttG_WG  --name capScore_10.0 --derivative  ctZ ctZ --max_local_score 10.0

#python training.py  --overwrite  --debug --input_directory /scratch-cbe/users/robert.schoefbeck/BIT/training-ntuples-ttG_WG/MVA-training/ttG_WG --config ttG_WG  --name capScore_1.0 --derivative  cWWW --max_local_score 1.0
#python training.py --overwrite  --debug --input_directory /scratch-cbe/users/robert.schoefbeck/BIT/training-ntuples-ttG_WG/MVA-training/ttG_WG --config ttG_WG  --name capScore_1.0 --derivative  cWWW cWWW --max_local_score 1.0
#python training.py --overwrite  --debug --input_directory /scratch-cbe/users/robert.schoefbeck/BIT/training-ntuples-ttG_WG/MVA-training/ttG_WG --config ttG_WG  --name capScore_1.0 --derivative  ctZ --max_local_score 1.0
#python training.py  --overwrite  --debug --input_directory /scratch-cbe/users/robert.schoefbeck/BIT/training-ntuples-ttG_WG/MVA-training/ttG_WG --config ttG_WG  --name capScore_1.0 --derivative  ctZ ctZ --max_local_score 1.0

#python training.py  --overwrite  --debug --input_directory /scratch-cbe/users/robert.schoefbeck/BIT/training-ntuples-ttG_WG/MVA-training/ttG_WG --config ttG_WG  --name capScore_5.0 --derivative  cWWW --max_local_score 5.0
python training.py --overwrite  --debug --input_directory /scratch-cbe/users/robert.schoefbeck/BIT/training-ntuples-ttG_WG/MVA-training/ttG_WG --config ttG_WG  --name capScore_5.0 --derivative  cWWW cWWW --max_local_score 5.0
python training.py --overwrite  --debug --input_directory /scratch-cbe/users/robert.schoefbeck/BIT/training-ntuples-ttG_WG/MVA-training/ttG_WG --config ttG_WG  --name capScore_5.0 --derivative  ctZ --max_local_score 5.0
#python training.py  --overwrite  --debug --input_directory /scratch-cbe/users/robert.schoefbeck/BIT/training-ntuples-ttG_WG/MVA-training/ttG_WG --config ttG_WG  --name capScore_5.0 --derivative  ctZ ctZ --max_local_score 5.0

#python training.py  --overwrite  --debug --input_directory /scratch-cbe/users/robert.schoefbeck/BIT/training-ntuples-ttG_WG/MVA-training/ttG_WG --config ttG_WG  --name capScore_50.0 --derivative  cWWW --max_local_score 50.0
#python training.py --overwrite  --debug --input_directory /scratch-cbe/users/robert.schoefbeck/BIT/training-ntuples-ttG_WG/MVA-training/ttG_WG --config ttG_WG  --name capScore_50.0 --derivative  cWWW cWWW --max_local_score 50.0
#python training.py --overwrite  --debug --input_directory /scratch-cbe/users/robert.schoefbeck/BIT/training-ntuples-ttG_WG/MVA-training/ttG_WG --config ttG_WG  --name capScore_50.0 --derivative  ctZ --max_local_score 50.0
#python training.py  --overwrite  --debug --input_directory /scratch-cbe/users/robert.schoefbeck/BIT/training-ntuples-ttG_WG/MVA-training/ttG_WG --config ttG_WG  --name capScore_50.0 --derivative  ctZ ctZ --max_local_score 50.0

#python training.py  --overwrite  --debug --input_directory /scratch-cbe/users/robert.schoefbeck/BIT/training-ntuples-ttG_WG/MVA-training/ttG_WG --config ttG_WG  --name capScore_0.1 --derivative  cWWW --max_local_score 0.1 
#python training.py --overwrite  --debug --input_directory /scratch-cbe/users/robert.schoefbeck/BIT/training-ntuples-ttG_WG/MVA-training/ttG_WG --config ttG_WG  --name capScore_0.1 --derivative  cWWW cWWW --max_local_score 0.1
#python training.py --overwrite  --debug --input_directory /scratch-cbe/users/robert.schoefbeck/BIT/training-ntuples-ttG_WG/MVA-training/ttG_WG --config ttG_WG  --name capScore_0.1 --derivative  ctZ --max_local_score 0.1
#python training.py  --overwrite  --debug --input_directory /scratch-cbe/users/robert.schoefbeck/BIT/training-ntuples-ttG_WG/MVA-training/ttG_WG --config ttG_WG  --name capScore_0.1 --derivative  ctZ ctZ --max_local_score 0.1

