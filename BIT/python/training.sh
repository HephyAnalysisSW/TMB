#python training.py --config ttZ_3l_flavor & 
#python training.py --config ttZ_3l_flavor --name MU0.05 --max_uncertainty 0.05 &
#python training.py --config ttZ_3l_flavor --name MU0.1 --max_uncertainty 0.1 &
#python training.py --config ttZ_3l_flavor --name MU0.2  --max_uncertainty 0.2 &
#python training.py --config ttZ_3l_flavor --name MU0.5  --max_uncertainty 0.5 &
#python training.py --config ttZ_3l_flavor --name MU1.0  --max_uncertainty 1.0 &

#python training.py --overwrite --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor/MVA-training/ttZ_3l_flavor --config ttZ_3l_flavor_ptZ --name UNC2 --max_uncertainty 2
#python training.py --overwrite --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor/MVA-training/ttZ_3l_flavor --config ttZ_3l_flavor_ptZ --name UNC1 --max_uncertainty 1
#python training.py --overwrite --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor/MVA-training/ttZ_3l_flavor --config ttZ_3l_flavor_ptZ --name UNC0p5 --max_uncertainty 0.5
#python training.py --overwrite --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor/MVA-training/ttZ_3l_flavor --config ttZ_3l_flavor_ptZ --name UNC0p4 --max_uncertainty 0.4
#python training.py --overwrite --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor/MVA-training/ttZ_3l_flavor --config ttZ_3l_flavor_ptZ --name UNC0p3 --max_uncertainty 0.3
#python training.py --overwrite --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor/MVA-training/ttZ_3l_flavor --config ttZ_3l_flavor_ptZ --name UNC0p2 --max_uncertainty 0.2
#python training.py --overwrite --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor/MVA-training/ttZ_3l_flavor --config ttZ_3l_flavor_ptZ --name UNC0p1 --max_uncertainty 0.1
#python training.py --overwrite --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor/MVA-training/ttZ_3l_flavor --config ttZ_3l_flavor_ptZ --name UNC0p05 --max_uncertainty 0.05
#python training.py --overwrite --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor/MVA-training/ttZ_3l_flavor --config ttZ_3l_flavor_ptZ --name UNC0p04 --max_uncertainty 0.04
#python training.py --overwrite --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor/MVA-training/ttZ_3l_flavor --config ttZ_3l_flavor_ptZ --name UNC0p03 --max_uncertainty 0.03
#python training.py --overwrite --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor/MVA-training/ttZ_3l_flavor --config ttZ_3l_flavor_ptZ --name UNC0p02 --max_uncertainty 0.02
#python training.py --overwrite --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor/MVA-training/ttZ_3l_flavor --config ttZ_3l_flavor_ptZ --name UNC0p01 --max_uncertainty 0.01 
#python training.py --overwrite --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor/MVA-training/ttZ_3l_flavor --config ttZ_3l_flavor_ptZ --name default

#python training.py --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor-delphes/MVA-training/ttZ_3l_flavor_delphes --config ttZ_3l_flavor_delphes_ptZ --name default
#python training.py --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor-delphes/MVA-training/ttZ_3l_flavor_delphes --config ttZ_3l_flavor_delphes_ptZ  --name BF0p6 --bagging_fraction 0.6
#python training.py --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor-delphes/MVA-training/ttZ_3l_flavor_delphes --config ttZ_3l_flavor_delphes_ptZ  --name BF0p6_MNS10 --bagging_fraction 0.6 --max_n_split 10
#python training.py --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor-delphes/MVA-training/ttZ_3l_flavor_delphes --config ttZ_3l_flavor_delphes_ptZ  --name MNS10  --max_n_split 10

#python training.py --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor-delphes/MVA-training/WZ_3l_flavor_delphes --config WZ_3l_flavor_delphes --name default
#python training.py --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor-delphes/MVA-training/ZZ_3l_flavor_delphes --config ZZ_3l_flavor_delphes --name default

#python training.py --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor-delphes/MVA-training/ZZ_3l_flavor_delphes --config ZZ_3l_flavor_delphes --name default --derivative  cHq1Re11
#python training.py --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor-delphes/MVA-training/ZZ_3l_flavor_delphes --config ZZ_3l_flavor_delphes --name default --derivative  cHq1Re22
#python training.py --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor-delphes/MVA-training/ZZ_3l_flavor_delphes --config ZZ_3l_flavor_delphes --name default --derivative  cHq1Re11 cHq1Re11
#python training.py --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor-delphes/MVA-training/ZZ_3l_flavor_delphes --config ZZ_3l_flavor_delphes --name default --derivative  cHq1Re22 cHq1Re22
#python training.py --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor-delphes/MVA-training/ZZ_3l_flavor_delphes --config ZZ_3l_flavor_delphes --name default --derivative  cHq3Re11
#python training.py --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor-delphes/MVA-training/ZZ_3l_flavor_delphes --config ZZ_3l_flavor_delphes --name default --derivative  cHq3Re22
#python training.py --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor-delphes/MVA-training/ZZ_3l_flavor_delphes --config ZZ_3l_flavor_delphes --name default --derivative  cHq3Re11 cHq3Re11
#python training.py --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor-delphes/MVA-training/ZZ_3l_flavor_delphes --config ZZ_3l_flavor_delphes --name default --derivative  cHq3Re22 cHq3Re22
#python training.py --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor-delphes/MVA-training/ZZ_3l_flavor_delphes --config ZZ_3l_flavor_delphes --name default --derivative  cHuRe11
#python training.py --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor-delphes/MVA-training/ZZ_3l_flavor_delphes --config ZZ_3l_flavor_delphes --name default --derivative  cHuRe22
#python training.py --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor-delphes/MVA-training/ZZ_3l_flavor_delphes --config ZZ_3l_flavor_delphes --name default --derivative    cHuRe11 cHuRe11
#python training.py --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor-delphes/MVA-training/ZZ_3l_flavor_delphes --config ZZ_3l_flavor_delphes --name default --derivative    cHuRe22 cHuRe22
#python training.py --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor-delphes/MVA-training/ZZ_3l_flavor_delphes --config ZZ_3l_flavor_delphes --name default --derivative  cHdRe11
#python training.py --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor-delphes/MVA-training/ZZ_3l_flavor_delphes --config ZZ_3l_flavor_delphes --name default --derivative  cHdRe22
#python training.py --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor-delphes/MVA-training/ZZ_3l_flavor_delphes --config ZZ_3l_flavor_delphes --name default --derivative    cHdRe11 cHdRe11
#python training.py --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor-delphes/MVA-training/ZZ_3l_flavor_delphes --config ZZ_3l_flavor_delphes --name default --derivative    cHdRe22 cHdRe22
#python training.py --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor-delphes/MVA-training/ZZ_3l_flavor_delphes --config ZZ_3l_flavor_delphes --name default --derivative  cHudRe11
#python training.py --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor-delphes/MVA-training/ZZ_3l_flavor_delphes --config ZZ_3l_flavor_delphes --name default --derivative  cHudRe22
#python training.py --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor-delphes/MVA-training/ZZ_3l_flavor_delphes --config ZZ_3l_flavor_delphes --name default --derivative  cHudRe11 cHudRe11
#python training.py --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor-delphes/MVA-training/ZZ_3l_flavor_delphes --config ZZ_3l_flavor_delphes --name default --derivative  cHudRe22 cHudRe22


#python training.py --overwrite  --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor-delphes/MVA-training/ttZ_3l_flavor_delphes --config ttZ_3l_flavor_delphes  --name default --derivative  cHq1Re33 
#python training.py --overwrite  --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor-delphes/MVA-training/ttZ_3l_flavor_delphes --config ttZ_3l_flavor_delphes  --name default --derivative  cHq1Re33 cHq1Re33
#python training.py --overwrite  --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor-delphes/MVA-training/ttZ_3l_flavor_delphes --config ttZ_3l_flavor_delphes  --name default --derivative  cHq3Re33 
#python training.py --overwrite  --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor-delphes/MVA-training/ttZ_3l_flavor_delphes --config ttZ_3l_flavor_delphes  --name default --derivative  cHq3Re33 cHq3Re33
#python training.py --overwrite  --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor-delphes/MVA-training/ttZ_3l_flavor_delphes --config ttZ_3l_flavor_delphes  --name default --derivative  cHuRe33 
#python training.py --overwrite  --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor-delphes/MVA-training/ttZ_3l_flavor_delphes --config ttZ_3l_flavor_delphes  --name default --derivative    cHuRe33 cHuRe33  
#python training.py --overwrite  --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor-delphes/MVA-training/ttZ_3l_flavor_delphes --config ttZ_3l_flavor_delphes  --name default --derivative  cHdRe33 
#python training.py --overwrite  --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor-delphes/MVA-training/ttZ_3l_flavor_delphes --config ttZ_3l_flavor_delphes  --name default --derivative    cHdRe33 cHdRe33  
#python training.py --overwrite  --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor-delphes/MVA-training/ttZ_3l_flavor_delphes --config ttZ_3l_flavor_delphes  --name default --derivative  cHudRe33 
#python training.py --overwrite  --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuples-TTZ-flavor-delphes/MVA-training/ttZ_3l_flavor_delphes --config ttZ_3l_flavor_delphes  --name default --derivative  cHudRe33 cHudRe33

#python training.py --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuples-TTZ-delphes/MVA-training/ttZ_3l_delphes --config ttZ_3l_delphes --lumi_norm --name gss 

#python training.py --overwrite  --debug  --clip_score_quantile 0.01 --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuples-ttG_WG/MVA-training/ttG_WG --config ttG_WG  --name clipScore2 --derivative  cWWW
#python training.py --overwrite  --debug  --clip_score_quantile 0.01 --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuples-ttG_WG/MVA-training/ttG_WG --config ttG_WG  --name cpSfix --derivative  cWWW cWWW --positive_score 
#python training.py --overwrite  --debug  --clip_score_quantile 0.01 --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuples-ttG_WG/MVA-training/ttG_WG --config ttG_WG  --name clipScore2 --derivative  ctZ
#python training.py --overwrite  --debug  --clip_score_quantile 0.01 --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuples-ttG_WG/MVA-training/ttG_WG --config ttG_WG  --name cpSfix --derivative  ctZ ctZ --positive_score 

#python training.py --overwrite  --debug  --clip_score_quantile 0.01 --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuples-ttG_WG/MVA-training/ttG_WG --config ttG_WG_VV  --name clip1p_v10 --derivative  cWWW
#python training.py --overwrite  --debug  --clip_score_quantile 0.01 --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuples-ttG_WG/MVA-training/ttG_WG --config ttG_WG_VV  --name clip1p_v10_psc --derivative  cWWW cWWW --positive_score 
#python training.py --overwrite  --debug  --clip_score_quantile 0.01 --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuples-ttG_WG/MVA-training/ttG_WG --config ttG_WG_VV  --name clip1p_v10 --derivative  cWWW cWWW 
#python training.py --overwrite  --debug  --clip_score_quantile 0.01 --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuples-ttG_WG/MVA-training/ttG_WG --config ttG_WG_VV  --name clip1p_v10 --derivative  ctZ
#python training.py --overwrite  --debug  --clip_score_quantile 0.01 --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuples-ttG_WG/MVA-training/ttG_WG --config ttG_WG_VV  --name clip1p_v10_psc --derivative  ctZ ctZ --positive_score 
#python training.py --overwrite  --debug  --clip_score_quantile 0.01 --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuples-ttG_WG/MVA-training/ttG_WG --config ttG_WG_VV  --name clip1p_v10 --derivative  ctZ ctZ 

#python training.py --overwrite  --debug  --clip_score_quantile 0.05 --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuples-WG/MVA-training/WG_delphes --config WG_delphes  --name default --derivative  cWWW
#python training.py --overwrite  --debug  --clip_score_quantile 0.05 --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuples-WG/MVA-training/WG_delphes --config WG_delphes  --name default --derivative  cWWW cWWW 

#python training.py --overwrite  --debug  --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-ZH/MVA-training/ZH_delphes_bkgs --config ZH_delphes  --name v2 --derivative cHW 
#python training.py --overwrite  --debug  --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-ZH/MVA-training/ZH_delphes_bkgs --config ZH_delphes  --name v2 --derivative cHWtil 
#python training.py --overwrite  --debug  --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-ZH/MVA-training/ZH_delphes_bkgs --config ZH_delphes  --name v2 --derivative cHWtil cHWtil
#python training.py --overwrite  --debug  --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-ZH/MVA-training/ZH_delphes_bkgs --config ZH_delphes  --name v2 --derivative cHW cHW
#
#python training.py --overwrite  --debug  --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-ZH/MVA-training/ZH_delphes_bkgs --config ZH_delphes  --name v2 --derivative cHj3 
#python training.py --overwrite  --debug  --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-ZH/MVA-training/ZH_delphes_bkgs --config ZH_delphes  --name v2 --derivative cHj3 cHj3 
#
#python training.py --overwrite  --debug  --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-ZH/MVA-training/ZH_delphes_bkgs --config ZH_delphes  --name v2 --derivative cHj3 cHW 
#python training.py --overwrite  --debug  --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-ZH/MVA-training/ZH_delphes_bkgs --config ZH_delphes  --name v2 --derivative cHj3 cHWtil
#python training.py --overwrite  --debug  --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-ZH/MVA-training/ZH_delphes_bkgs --config ZH_delphes  --name v2 --derivative cHW cHWtil


#python training.py --overwrite  --debug  --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-WH/MVA-training/WH_delphes_bkgs --config WH_delphes  --name v2 --derivative cHW 
#python training.py --overwrite  --debug  --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-WH/MVA-training/WH_delphes_bkgs --config WH_delphes  --name v2 --derivative cHWtil 
#python training.py --overwrite  --debug  --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-WH/MVA-training/WH_delphes_bkgs --config WH_delphes  --name v2 --derivative cHWtil cHWtil
#python training.py --overwrite  --debug  --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-WH/MVA-training/WH_delphes_bkgs --config WH_delphes  --name v2 --derivative cHW cHW
#
#python training.py --overwrite  --debug  --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-WH/MVA-training/WH_delphes_bkgs --config WH_delphes  --name v2 --derivative cHj3 
#python training.py --overwrite  --debug  --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-WH/MVA-training/WH_delphes_bkgs --config WH_delphes  --name v2 --derivative cHj3 cHj3 
#
#python training.py --overwrite  --debug  --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-WH/MVA-training/WH_delphes_bkgs --config WH_delphes  --name v2 --derivative cHj3 cHW 
#python training.py --overwrite  --debug  --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-WH/MVA-training/WH_delphes_bkgs --config WH_delphes  --name v2 --derivative cHj3 cHWtil
#python training.py --overwrite  --debug  --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-WH/MVA-training/WH_delphes_bkgs --config WH_delphes  --name v2 --derivative cHW cHWtil

#python training.py --debug --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-WH/MVA-training/WH_delphes_bkgs --config WH_delphes_bkgs  --name v2 --derivative cHW 
#python training.py --debug --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-WH/MVA-training/WH_delphes_bkgs --config WH_delphes_bkgs  --name v2 --derivative cHWtil 
#python training.py --debug --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-WH/MVA-training/WH_delphes_bkgs --config WH_delphes_bkgs  --name v2 --derivative cHWtil cHWtil
#python training.py --debug --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-WH/MVA-training/WH_delphes_bkgs --config WH_delphes_bkgs  --name v2 --derivative cHW cHW
#
#python training.py --debug --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-WH/MVA-training/WH_delphes_bkgs --config WH_delphes_bkgs  --name v2 --derivative cHj3 
#python training.py --debug --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-WH/MVA-training/WH_delphes_bkgs --config WH_delphes_bkgs  --name v2 --derivative cHj3 cHj3 
#
#python training.py --debug --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-WH/MVA-training/WH_delphes_bkgs --config WH_delphes_bkgs  --name v2 --derivative cHj3 cHW 
#python training.py --debug --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-WH/MVA-training/WH_delphes_bkgs --config WH_delphes_bkgs  --name v2 --derivative cHj3 cHWtil
#python training.py --debug --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-WH/MVA-training/WH_delphes_bkgs --config WH_delphes_bkgs  --name v2 --derivative cHW cHWtil

#python training.py --debug --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-WH/MVA-training/WH_delphes_bkgs_ptW200 --config WH_delphes_bkgs_ptW200  --name v2 --derivative cHW 
#python training.py --debug --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-WH/MVA-training/WH_delphes_bkgs_ptW200 --config WH_delphes_bkgs_ptW200  --name v2 --derivative cHWtil 
#python training.py --debug --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-WH/MVA-training/WH_delphes_bkgs_ptW200 --config WH_delphes_bkgs_ptW200  --name v2 --derivative cHWtil cHWtil
#python training.py --debug --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-WH/MVA-training/WH_delphes_bkgs_ptW200 --config WH_delphes_bkgs_ptW200  --name v2 --derivative cHW cHW
#
#python training.py --debug --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-WH/MVA-training/WH_delphes_bkgs_ptW200 --config WH_delphes_bkgs_ptW200  --name v2 --derivative cHj3 
#python training.py --debug --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-WH/MVA-training/WH_delphes_bkgs_ptW200 --config WH_delphes_bkgs_ptW200  --name v2 --derivative cHj3 cHj3 
#
#python training.py --debug --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-WH/MVA-training/WH_delphes_bkgs_ptW200 --config WH_delphes_bkgs_ptW200  --name v2 --derivative cHj3 cHW 
#python training.py --debug --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-WH/MVA-training/WH_delphes_bkgs_ptW200 --config WH_delphes_bkgs_ptW200  --name v2 --derivative cHj3 cHWtil
#python training.py --debug --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-WH/MVA-training/WH_delphes_bkgs_ptW200 --config WH_delphes_bkgs_ptW200  --name v2 --derivative cHW cHWtil

#python training.py  --debug  --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-ZH/MVA-training/ZH_nlo_delphes --config ZH_nlo_delphes  --name v1 --derivative cpW 
#python training.py  --debug  --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-ZH/MVA-training/ZH_nlo_delphes --config ZH_nlo_delphes  --name v1 --derivative cpqMi 
#python training.py  --debug  --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-ZH/MVA-training/ZH_nlo_delphes --config ZH_nlo_delphes  --name v1 --derivative cpqMi cpqMi
#python training.py  --debug  --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-ZH/MVA-training/ZH_nlo_delphes --config ZH_nlo_delphes  --name v1 --derivative cpW cpW
#
#python training.py  --debug  --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-ZH/MVA-training/ZH_nlo_delphes --config ZH_nlo_delphes  --name v1 --derivative cpq3i 
#python training.py  --debug  --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-ZH/MVA-training/ZH_nlo_delphes --config ZH_nlo_delphes  --name v1 --derivative cpq3i cpq3i 
#
#python training.py  --debug  --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-ZH/MVA-training/ZH_nlo_delphes --config ZH_nlo_delphes  --name v1 --derivative cpq3i cpW 
#python training.py  --debug  --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-ZH/MVA-training/ZH_nlo_delphes --config ZH_nlo_delphes  --name v1 --derivative cpq3i cpqMi
#python training.py  --debug  --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-ZH/MVA-training/ZH_nlo_delphes --config ZH_nlo_delphes  --name v1 --derivative cpW cpqMi

python training.py  --debug  --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-ZH/MVA-training/ZH_delphes_bkgs_comb --config ZH_delphes_bkgs_comb  --max_depth 5 --name v2_MD5  --derivative cHW 
python training.py  --debug  --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-ZH/MVA-training/ZH_delphes_bkgs_comb --config ZH_delphes_bkgs_comb  --max_depth 5 --name v2_MD5  --derivative cHWtil 
python training.py  --debug  --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-ZH/MVA-training/ZH_delphes_bkgs_comb --config ZH_delphes_bkgs_comb  --max_depth 5 --name v2_MD5  --derivative cHWtil cHWtil
python training.py  --debug  --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-ZH/MVA-training/ZH_delphes_bkgs_comb --config ZH_delphes_bkgs_comb  --max_depth 5 --name v2_MD5  --derivative cHW cHW
                                                                                                                                                                           --max_depth 5 --name v2_MD5                             
python training.py  --debug  --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-ZH/MVA-training/ZH_delphes_bkgs_comb --config ZH_delphes_bkgs_comb  --max_depth 5 --name v2_MD5  --derivative cHj3 
python training.py  --debug  --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-ZH/MVA-training/ZH_delphes_bkgs_comb --config ZH_delphes_bkgs_comb  --max_depth 5 --name v2_MD5  --derivative cHj3 cHj3 
                                                                                                                                                                           --max_depth 5 --name v2_MD5                             
python training.py  --debug  --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-ZH/MVA-training/ZH_delphes_bkgs_comb --config ZH_delphes_bkgs_comb  --max_depth 5 --name v2_MD5  --derivative cHj3 cHW 
python training.py  --debug  --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-ZH/MVA-training/ZH_delphes_bkgs_comb --config ZH_delphes_bkgs_comb  --max_depth 5 --name v2_MD5  --derivative cHj3 cHWtil
python training.py  --debug  --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-ZH/MVA-training/ZH_delphes_bkgs_comb --config ZH_delphes_bkgs_comb  --max_depth 5 --name v2_MD5  --derivative cHW cHWtil

#python training.py --overwrite  --debug  --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-WH/MVA-training/WH_nlo_delphes --config WH_nlo_delphes  --max_depth 2 --n_trees 300 --clip_score_quantile 0.001 --name v1_clipScore_0p001_max_depth_2_n_trees_300 --derivative cpW 
#python training.py --overwrite  --debug  --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-WH/MVA-training/WH_nlo_delphes --config WH_nlo_delphes  --max_depth 2 --n_trees 300 --clip_score_quantile 0.001 --name v1_clipScore_0p001_max_depth_2_n_trees_300 --derivative cpqMi 
#python training.py --overwrite  --debug  --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-WH/MVA-training/WH_nlo_delphes --config WH_nlo_delphes  --max_depth 2 --n_trees 300 --clip_score_quantile 0.001 --name v1_clipScore_0p001_max_depth_2_n_trees_300 --derivative cpqMi cpqMi
#python training.py --overwrite  --debug  --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-WH/MVA-training/WH_nlo_delphes --config WH_nlo_delphes  --max_depth 2 --n_trees 300 --clip_score_quantile 0.001 --name v1_clipScore_0p001_max_depth_2_n_trees_300 --derivative cpW cpW
#
#python training.py --overwrite  --debug  --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-WH/MVA-training/WH_nlo_delphes --config WH_nlo_delphes  --max_depth 2 --n_trees 300 --clip_score_quantile 0.001 --name v1_clipScore_0p001_max_depth_2_n_trees_300 --derivative cpq3i 
#python training.py --overwrite  --debug  --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-WH/MVA-training/WH_nlo_delphes --config WH_nlo_delphes  --max_depth 2 --n_trees 300 --clip_score_quantile 0.001 --name v1_clipScore_0p001_max_depth_2_n_trees_300 --derivative cpq3i cpq3i 
#
#python training.py --overwrite  --debug  --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-WH/MVA-training/WH_nlo_delphes --config WH_nlo_delphes  --max_depth 2 --n_trees 300 --clip_score_quantile 0.001 --name v1_clipScore_0p001_max_depth_2_n_trees_300 --derivative cpq3i cpW 
#python training.py --overwrite  --debug  --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-WH/MVA-training/WH_nlo_delphes --config WH_nlo_delphes  --max_depth 2 --n_trees 300 --clip_score_quantile 0.001 --name v1_clipScore_0p001_max_depth_2_n_trees_300 --derivative cpq3i cpqMi
#python training.py --overwrite  --debug  --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-WH/MVA-training/WH_nlo_delphes --config WH_nlo_delphes  --max_depth 2 --n_trees 300 --clip_score_quantile 0.001 --name v1_clipScore_0p001_max_depth_2_n_trees_300 --derivative cpW cpqMi


#python training.py  --debug  --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-ZH/MVA-training/ZH_delphes_bkgs_comb_ptZ200 --config ZH_delphes_bkgs_comb_ptZ200  --name v2 --derivative cHW 
#python training.py  --debug  --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-ZH/MVA-training/ZH_delphes_bkgs_comb_ptZ200 --config ZH_delphes_bkgs_comb_ptZ200  --name v2 --derivative cHWtil 
#python training.py  --debug  --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-ZH/MVA-training/ZH_delphes_bkgs_comb_ptZ200 --config ZH_delphes_bkgs_comb_ptZ200  --name v2 --derivative cHWtil cHWtil
#python training.py  --debug  --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-ZH/MVA-training/ZH_delphes_bkgs_comb_ptZ200 --config ZH_delphes_bkgs_comb_ptZ200  --name v2 --derivative cHW cHW
#                                                                                                                                                                                                                
#python training.py  --debug  --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-ZH/MVA-training/ZH_delphes_bkgs_comb_ptZ200 --config ZH_delphes_bkgs_comb_ptZ200  --name v2 --derivative cHj3 
#python training.py  --debug  --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-ZH/MVA-training/ZH_delphes_bkgs_comb_ptZ200 --config ZH_delphes_bkgs_comb_ptZ200  --name v2 --derivative cHj3 cHj3 
#                                                                                                                                                                                                                
#python training.py  --debug  --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-ZH/MVA-training/ZH_delphes_bkgs_comb_ptZ200 --config ZH_delphes_bkgs_comb_ptZ200  --name v2 --derivative cHj3 cHW 
#python training.py  --debug  --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-ZH/MVA-training/ZH_delphes_bkgs_comb_ptZ200 --config ZH_delphes_bkgs_comb_ptZ200  --name v2 --derivative cHj3 cHWtil
#python training.py  --debug  --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-ZH/MVA-training/ZH_delphes_bkgs_comb_ptZ200 --config ZH_delphes_bkgs_comb_ptZ200  --name v2 --derivative cHW cHWtil

#python training.py  --debug  --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-ZH/MVA-training/ZH_delphes_bkgs_comb_ptZ200 --config ZH_delphes_bkgs_comb_ptZ200  --max_depth 5 --name v2_MD5 --derivative cHW 
#python training.py  --debug  --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-ZH/MVA-training/ZH_delphes_bkgs_comb_ptZ200 --config ZH_delphes_bkgs_comb_ptZ200  --max_depth 5 --name v2_MD5 --derivative cHWtil 
#python training.py  --debug  --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-ZH/MVA-training/ZH_delphes_bkgs_comb_ptZ200 --config ZH_delphes_bkgs_comb_ptZ200  --max_depth 5 --name v2_MD5 --derivative cHWtil cHWtil
#python training.py  --debug  --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-ZH/MVA-training/ZH_delphes_bkgs_comb_ptZ200 --config ZH_delphes_bkgs_comb_ptZ200  --max_depth 5 --name v2_MD5 --derivative cHW cHW
#python training.py  --debug  --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-ZH/MVA-training/ZH_delphes_bkgs_comb_ptZ200 --config ZH_delphes_bkgs_comb_ptZ200  --max_depth 5 --name v2_MD5 --derivative cHj3 
#python training.py  --debug  --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-ZH/MVA-training/ZH_delphes_bkgs_comb_ptZ200 --config ZH_delphes_bkgs_comb_ptZ200  --max_depth 5 --name v2_MD5 --derivative cHj3 cHj3 
#python training.py  --debug  --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-ZH/MVA-training/ZH_delphes_bkgs_comb_ptZ200 --config ZH_delphes_bkgs_comb_ptZ200  --max_depth 5 --name v2_MD5 --derivative cHj3 cHW 
#python training.py  --debug  --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-ZH/MVA-training/ZH_delphes_bkgs_comb_ptZ200 --config ZH_delphes_bkgs_comb_ptZ200  --max_depth 5 --name v2_MD5 --derivative cHj3 cHWtil
#python training.py  --debug  --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-ZH/MVA-training/ZH_delphes_bkgs_comb_ptZ200 --config ZH_delphes_bkgs_comb_ptZ200  --max_depth 5 --name v2_MD5 --derivative cHW cHWtil

#python training.py  --debug  --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-ZH/MVA-training/ZH_delphes_bkgs_comb_ptZ200 --config ZH_delphes_ptZ200  --max_depth 5 --name v2_MD5 --derivative cHW 
#python training.py  --debug  --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-ZH/MVA-training/ZH_delphes_bkgs_comb_ptZ200 --config ZH_delphes_ptZ200  --max_depth 5 --name v2_MD5 --derivative cHWtil 
#python training.py  --debug  --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-ZH/MVA-training/ZH_delphes_bkgs_comb_ptZ200 --config ZH_delphes_ptZ200  --max_depth 5 --name v2_MD5 --derivative cHWtil cHWtil
#python training.py  --debug  --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-ZH/MVA-training/ZH_delphes_bkgs_comb_ptZ200 --config ZH_delphes_ptZ200  --max_depth 5 --name v2_MD5 --derivative cHW cHW
#python training.py  --debug  --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-ZH/MVA-training/ZH_delphes_bkgs_comb_ptZ200 --config ZH_delphes_ptZ200  --max_depth 5 --name v2_MD5 --derivative cHj3 
#python training.py  --debug  --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-ZH/MVA-training/ZH_delphes_bkgs_comb_ptZ200 --config ZH_delphes_ptZ200  --max_depth 5 --name v2_MD5 --derivative cHj3 cHj3 
#python training.py  --debug  --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-ZH/MVA-training/ZH_delphes_bkgs_comb_ptZ200 --config ZH_delphes_ptZ200  --max_depth 5 --name v2_MD5 --derivative cHj3 cHW 
#python training.py  --debug  --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-ZH/MVA-training/ZH_delphes_bkgs_comb_ptZ200 --config ZH_delphes_ptZ200  --max_depth 5 --name v2_MD5 --derivative cHj3 cHWtil
#python training.py  --debug  --input_directory /groups/hephy/cms/robert.schoefbeck/BIT/training-ntuple-ZH/MVA-training/ZH_delphes_bkgs_comb_ptZ200 --config ZH_delphes_ptZ200  --max_depth 5 --name v2_MD5 --derivative cHW cHWtil
