#!/bin/sh

#python genPostProcessing.py --delphesEra RunII --targetDir v10 --overwrite target --addReweights --logLevel INFO --interpolationOrder 2  --sample WGToLNu #SPLIT598
#python genPostProcessing.py --delphesEra RunII --targetDir v10 --overwrite target --addReweights --logLevel INFO --interpolationOrder 2  --sample WGjj #SPLIT589
#python genPostProcessing.py --delphesEra RunII --targetDir v10 --overwrite target --addReweights --logLevel INFO --interpolationOrder 2  --sample ZZ #SPLIT575
#python genPostProcessing.py --delphesEra RunII --targetDir v10 --overwrite target --addReweights --logLevel INFO --interpolationOrder 2  --sample ZZjj #SPLIT513
#
#python genPostProcessing.py --delphesEra RunII --targetDir v10 --overwrite target --addReweights --logLevel INFO --interpolationOrder 2  --sample WG20To130ToLNu #SPLIT581
#python genPostProcessing.py --delphesEra RunII --targetDir v10 --overwrite target --addReweights --logLevel INFO --interpolationOrder 2  --sample WG130To300ToLNu #SPLIT598
#python genPostProcessing.py --delphesEra RunII --targetDir v10 --overwrite target --addReweights --logLevel INFO --interpolationOrder 2  --sample WG300To500ToLNu #SPLIT597
#python genPostProcessing.py --delphesEra RunII --targetDir v10 --overwrite target --addReweights --logLevel INFO --interpolationOrder 2  --sample WG500ToLNu #SPLIT595

#python genPostProcessing.py --delphesEra RunII --targetDir v10 --addReweights --logLevel INFO --interpolationOrder 2  --sample ZH #SPLIT580
#python genPostProcessing.py --delphesEra RunII --targetDir v10 --addReweights --logLevel INFO --interpolationOrder 2  --sample WH #SPLIT585

#python genPostProcessing.py --delphesEra RunII --targetDir v10 --addReweights --logLevel INFO --interpolationOrder 4  --sample ZH_nlo #SPLIT599
#python genPostProcessing.py --delphesEra RunII --targetDir v10 --addReweights --logLevel INFO --interpolationOrder 4  --sample WH_nlo #SPLIT599
#
#python genPostProcessing.py --delphesEra RunII --targetDir v10 --addReweights --logLevel INFO --interpolationOrder 2  --sample ZH0jNoEFTDecay #SPLIT998
#python genPostProcessing.py --delphesEra RunII --targetDir v10 --addReweights --logLevel INFO --interpolationOrder 2  --sample ZH0jEFTDecay #SPLIT997

#python genPostProcessing.py --combinatoricalBTags --miniAOD --delphesEra RunII --targetDir v10_combinatoricalBTagWeights --logLevel INFO  --sample DYBBJetsToLL_M50_LO #SPLIT145
#python genPostProcessing.py --combinatoricalBTags --miniAOD --delphesEra RunII --targetDir v10_combinatoricalBTagWeights --logLevel INFO  --sample DYJetsToLL_M50_HT70to100_LO #SPLIT237
#python genPostProcessing.py --combinatoricalBTags --miniAOD --delphesEra RunII --targetDir v10_combinatoricalBTagWeights --logLevel INFO  --sample DYJetsToLL_M50_HT100to200_LO #SPLIT232
#python genPostProcessing.py --combinatoricalBTags --miniAOD --delphesEra RunII --targetDir v10_combinatoricalBTagWeights --logLevel INFO  --sample DYJetsToLL_M50_HT200to400_LO #SPLIT311
#python genPostProcessing.py --combinatoricalBTags --miniAOD --delphesEra RunII --targetDir v10_combinatoricalBTagWeights --logLevel INFO  --sample DYJetsToLL_M50_HT400to600_LO_ext2 #SPLIT253
#python genPostProcessing.py --combinatoricalBTags --miniAOD --delphesEra RunII --targetDir v10_combinatoricalBTagWeights --logLevel INFO  --sample DYJetsToLL_M50_HT600to800_LO #SPLIT264
#python genPostProcessing.py --combinatoricalBTags --miniAOD --delphesEra RunII --targetDir v10_combinatoricalBTagWeights --logLevel INFO  --sample DYJetsToLL_M50_HT800to1200_LO #SPLIT109
#python genPostProcessing.py --combinatoricalBTags --miniAOD --delphesEra RunII --targetDir v10_combinatoricalBTagWeights --logLevel INFO  --sample DYJetsToLL_M50_HT1200to2500_LO #SPLIT23
#python genPostProcessing.py --combinatoricalBTags --miniAOD --delphesEra RunII --targetDir v10_combinatoricalBTagWeights --logLevel INFO  --sample DYJetsToLL_M50_HT2500toInf_LO #SPLIT23
#
#python genPostProcessing.py --combinatoricalBTags           --delphesEra RunII --targetDir v10_combinatoricalBTagWeights --logLevel INFO  --sample DYJets #SPLIT595
#
#python genPostProcessing.py --miniAOD --delphesEra RunII --targetDir v10 --logLevel INFO  --sample DYBBJetsToLL_M50_LO #SPLIT145
#python genPostProcessing.py --miniAOD --delphesEra RunII --targetDir v10 --logLevel INFO  --sample DYJetsToLL_M50_HT70to100_LO #SPLIT237
#python genPostProcessing.py --miniAOD --delphesEra RunII --targetDir v10 --logLevel INFO  --sample DYJetsToLL_M50_HT100to200_LO #SPLIT232
#python genPostProcessing.py --miniAOD --delphesEra RunII --targetDir v10 --logLevel INFO  --sample DYJetsToLL_M50_HT200to400_LO #SPLIT311
#python genPostProcessing.py --miniAOD --delphesEra RunII --targetDir v10 --logLevel INFO  --sample DYJetsToLL_M50_HT400to600_LO_ext2 #SPLIT253
#python genPostProcessing.py --miniAOD --delphesEra RunII --targetDir v10 --logLevel INFO  --sample DYJetsToLL_M50_HT600to800_LO #SPLIT264
#python genPostProcessing.py --miniAOD --delphesEra RunII --targetDir v10 --logLevel INFO  --sample DYJetsToLL_M50_HT800to1200_LO #SPLIT109
#python genPostProcessing.py --miniAOD --delphesEra RunII --targetDir v10 --logLevel INFO  --sample DYJetsToLL_M50_HT1200to2500_LO #SPLIT23
#python genPostProcessing.py --miniAOD --delphesEra RunII --targetDir v10 --logLevel INFO  --sample DYJetsToLL_M50_HT2500toInf_LO #SPLIT23
#
#python genPostProcessing.py           --delphesEra RunII --targetDir v10 --logLevel INFO  --sample DYJets #SPLIT595

#python genPostProcessing.py --miniAOD --delphesEra RunII --targetDir v10 --logLevel INFO  --sample DYJetsToLL_M50_HT400to600_LO #SPLIT249
#
#python genPostProcessing.py --delphesEra RunII --targetDir v10 --logLevel INFO  --sample TTJets #SPLIT1873
#python genPostProcessing.py --delphesEra RunII --targetDir v10 --logLevel INFO  --sample WJetsToLNu #SPLIT596

##python genPostProcessing.py --delphesEra RunII --targetDir v10 --logLevel INFO  --sample WJetsToLNu_HT_70to100 #SPLIT498
##python genPostProcessing.py --delphesEra RunII --targetDir v10 --logLevel INFO  --sample WJetsToLNu_HT_100to200 #SPLIT499
##python genPostProcessing.py --delphesEra RunII --targetDir v10 --logLevel INFO  --sample WJetsToLNu_HT_200to400 #SPLIT500
##python genPostProcessing.py --delphesEra RunII --targetDir v10 --logLevel INFO  --sample WJetsToLNu_HT_400to600 #SPLIT200
##python genPostProcessing.py --delphesEra RunII --targetDir v10 --logLevel INFO  --sample WJetsToLNu_HT_600to800 #SPLIT198
##python genPostProcessing.py --delphesEra RunII --targetDir v10 --logLevel INFO  --sample WJetsToLNu_HT_800to1200 #SPLIT196
##python genPostProcessing.py --delphesEra RunII --targetDir v10 --logLevel INFO  --sample WJetsToLNu_HT_1200to2500 #SPLIT200
##python genPostProcessing.py --delphesEra RunII --targetDir v10 --logLevel INFO  --sample WJetsToLNu_HT_2500toInf #SPLIT196

#python genPostProcessing.py --miniAOD --delphesEra RunII --targetDir v10 --logLevel INFO  --sample WJetsToLNu_HT_70To100 #SPLIT514
#python genPostProcessing.py --miniAOD --delphesEra RunII --targetDir v10 --logLevel INFO  --sample WJetsToLNu_HT_100To200 #SPLIT562
#python genPostProcessing.py --miniAOD --delphesEra RunII --targetDir v10 --logLevel INFO  --sample WJetsToLNu_HT_200To400 #SPLIT522
#python genPostProcessing.py --miniAOD --delphesEra RunII --targetDir v10 --logLevel INFO  --sample WJetsToLNu_HT_400To600 #SPLIT239
#python genPostProcessing.py --miniAOD --delphesEra RunII --targetDir v10 --logLevel INFO  --sample WJetsToLNu_HT_600To800 #SPLIT595
#python genPostProcessing.py --miniAOD --delphesEra RunII --targetDir v10 --logLevel INFO  --sample WJetsToLNu_HT_800To1200 #SPLIT327
#python genPostProcessing.py --miniAOD --delphesEra RunII --targetDir v10 --logLevel INFO  --sample WJetsToLNu_HT_1200To2500 #SPLIT299
#python genPostProcessing.py --miniAOD --delphesEra RunII --targetDir v10 --logLevel INFO  --sample WJetsToLNu_HT_2500ToInf #SPLIT132
#
#
#python genPostProcessing.py --miniAOD --delphesEra RunII --targetDir v10 --logLevel INFO  --sample TTJets_DiLept                       #SPLIT841
#python genPostProcessing.py --miniAOD --delphesEra RunII --targetDir v10 --logLevel INFO  --sample TTJets_DiLept_genMET80              #SPLIT1578
#python genPostProcessing.py --miniAOD --delphesEra RunII --targetDir v10 --logLevel INFO  --sample TTJets_HT_1200to2500                #SPLIT123
#python genPostProcessing.py --miniAOD --delphesEra RunII --targetDir v10 --logLevel INFO  --sample TTJets_HT_2500toInf                 #SPLIT63
#python genPostProcessing.py --miniAOD --delphesEra RunII --targetDir v10 --logLevel INFO  --sample TTJets_HT_600to800                  #SPLIT503
#python genPostProcessing.py --miniAOD --delphesEra RunII --targetDir v10 --logLevel INFO  --sample TTJets_HT_800to1200                 #SPLIT399
#python genPostProcessing.py --miniAOD --delphesEra RunII --targetDir v10 --logLevel INFO  --sample TTJets_SingleLeptFromT              #SPLIT1505
#python genPostProcessing.py --miniAOD --delphesEra RunII --targetDir v10 --logLevel INFO  --sample TTJets_SingleLeptFromT_genMET80     #SPLIT2230
#python genPostProcessing.py --miniAOD --delphesEra RunII --targetDir v10 --logLevel INFO  --sample TTJets_SingleLeptFromTbar           #SPLIT1540
#python genPostProcessing.py --miniAOD --delphesEra RunII --targetDir v10 --logLevel INFO  --sample TTJets_SingleLeptFromTbar_genMET80  #SPLIT2055
#python genPostProcessing.py --miniAOD --delphesEra RunII --targetDir v10 --logLevel INFO  --sample TTJets                              #SPLIT309
#

#python genPostProcessing.py --delphesEra RunII --targetDir v9 --overwrite --addReweights --logLevel INFO --interpolationOrder 2  --sample WZTo3L1Nu #SPLIT599
#python genPostProcessing.py --delphesEra RunII --targetDir v9 --overwrite --addReweights --logLevel INFO --interpolationOrder 2  --sample ttZ01j #SPLIT550

#python genPostProcessing.py --delphesEra RunII --targetDir v5 --overwrite --addReweights --logLevel INFO --interpolationOrder 2  --sample ttZ01j #SPLIT600

#python genPostProcessing.py --delphesEra RunII --targetDir v8 --overwrite --addReweights --logLevel INFO --interpolationOrder 2  --sample ttZ01j #SPLIT600
#python genPostProcessing.py --delphesEra RunII --targetDir v8 --overwrite --addReweights --logLevel INFO --interpolationOrder 2  --sample WZTo3L1Nu #SPLIT600
#python genPostProcessing.py --delphesEra RunII --targetDir v8 --overwrite --addReweights --logLevel INFO --interpolationOrder 2  --sample ZZ #SPLIT600

#python genPostProcessing.py --targetDir v6 --overwrite --addReweights --logLevel INFO --interpolationOrder 2  --sample ttG_noFullyHad_test #SPLIT10
#python genPostProcessing.py --targetDir v6 --overwrite --addReweights --logLevel INFO --interpolationOrder 2  --sample ttg_noFullyHad_test2 #SPLIT10

#python genPostProcessing.py --targetDir v7 --overwrite --addReweights --logLevel INFO --interpolationOrder 2  --sample ttG_noFullyHad #SPLIT200
#python genPostProcessing.py --targetDir v7 --overwrite --addReweights --logLevel INFO --interpolationOrder 2  --sample ttW01j #SPLIT200
#python genPostProcessing.py --targetDir v7 --overwrite --addReweights --logLevel INFO --interpolationOrder 2  --sample ttZ01j #SPLIT200
#python genPostProcessing.py --targetDir v7 --overwrite --addReweights --logLevel INFO --interpolationOrder 2  --sample WGToLNu #SPLIT200
#python genPostProcessing.py --targetDir v7 --overwrite --addReweights --logLevel INFO --interpolationOrder 2  --sample ZGTo2L #SPLIT200
#python genPostProcessing.py --targetDir v7 --overwrite --addReweights --logLevel INFO --interpolationOrder 2  --sample WZTo3L1Nu #SPLIT200
#python genPostProcessing.py --targetDir v7 --overwrite --addReweights --logLevel INFO --interpolationOrder 2  --sample WW #SPLIT200
#python genPostProcessing.py --targetDir v7 --overwrite --addReweights --logLevel INFO --interpolationOrder 2  --sample ZZ #SPLIT200

#python genPostProcessing.py --targetDir v4 --overwrite --addReweights --logLevel INFO --interpolationOrder 2  --sample ttg_noFullyHad #SPLIT200
#python genPostProcessing.py --targetDir v4 --overwrite --addReweights --logLevel INFO --interpolationOrder 2  --sample ttZ01j #SPLIT200
#python genPostProcessing.py --targetDir v4 --overwrite --addReweights --logLevel INFO --interpolationOrder 2  --sample ttW01j #SPLIT200
#
#python genPostProcessing.py --targetDir v4 --overwrite --addReweights --logLevel INFO --interpolationOrder 2  --sample WW #SPLIT100
#python genPostProcessing.py --targetDir v4 --overwrite --addReweights --logLevel INFO --interpolationOrder 2  --sample WZ #SPLIT100
#python genPostProcessing.py --targetDir v4 --overwrite --addReweights --logLevel INFO --interpolationOrder 2  --sample ZZ #SPLIT100
#python genPostProcessing.py --targetDir v4 --overwrite --addReweights --logLevel INFO --interpolationOrder 2  --sample WA #SPLIT200
#python genPostProcessing.py --targetDir v4 --overwrite --addReweights --logLevel INFO --interpolationOrder 2  --sample WA_LO #SPLIT200
#python genPostProcessing.py --targetDir v4 --overwrite --addReweights --logLevel INFO --interpolationOrder 2  --sample WAjj #SPLIT200
#python genPostProcessing.py --targetDir v4 --overwrite --addReweights --logLevel INFO --interpolationOrder 2  --sample WWjj_OS #SPLIT200
#python genPostProcessing.py --targetDir v4 --overwrite --addReweights --logLevel INFO --interpolationOrder 2  --sample WWjj_SS #SPLIT200
#python genPostProcessing.py --targetDir v4 --overwrite --addReweights --logLevel INFO --interpolationOrder 2  --sample WZjj #SPLIT200
#python genPostProcessing.py --targetDir v4 --overwrite --addReweights --logLevel INFO --interpolationOrder 2  --sample ZA #SPLIT200
#python genPostProcessing.py --targetDir v4 --overwrite --addReweights --logLevel INFO --interpolationOrder 2  --sample ZAjj #SPLIT200
#python genPostProcessing.py --targetDir v4 --overwrite --addReweights --logLevel INFO --interpolationOrder 2  --sample ZZjj #SPLIT200
#
#python genPostProcessing.py --targetDir v4 --overwrite --addReweights --logLevel INFO --interpolationOrder 2  --sample WWW #SPLIT200
#python genPostProcessing.py --targetDir v4 --overwrite --addReweights --logLevel INFO --interpolationOrder 2  --sample WWZ #SPLIT200
#python genPostProcessing.py --targetDir v4 --overwrite --addReweights --logLevel INFO --interpolationOrder 2  --sample WZZ #SPLIT200
#python genPostProcessing.py --targetDir v4 --overwrite --addReweights --logLevel INFO --interpolationOrder 2  --sample ZZZ #SPLIT200
