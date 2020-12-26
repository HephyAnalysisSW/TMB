#!/bin/sh

#use cardfile delphes_card_ATLAS or delphes_card_CMS

python genPostProcessing.py --targetDir v1 --logLevel INFO --interpolationOrder 2  --sample WW #SPLIT100 
python genPostProcessing.py --targetDir v1 --logLevel INFO --interpolationOrder 2  --sample WZ #SPLIT100 
python genPostProcessing.py --targetDir v1 --logLevel INFO --interpolationOrder 2  --sample ZZ #SPLIT100 
python genPostProcessing.py --targetDir v1 --logLevel INFO --interpolationOrder 2  --sample WWW #SPLIT200
python genPostProcessing.py --targetDir v1 --logLevel INFO --interpolationOrder 2  --sample WWZ #SPLIT200
python genPostProcessing.py --targetDir v1 --logLevel INFO --interpolationOrder 2  --sample WZZ #SPLIT200
python genPostProcessing.py --targetDir v1 --logLevel INFO --interpolationOrder 2  --sample ZZZ #SPLIT200

#######
##CMS##
#######


# ttZ Signal samples with reference point
#python genPostProcessing.py --overwrite target --targetDir v23 --logLevel DEBUG --sample fwlite_ttZ_ll_LO_order2_15weights_ref --addReweights --interpolationOrder 2  --delphes --delphesEra RunII #SPLIT200

### ttZ background samples without reference point
#python genPostProcessing.py --overwrite target --targetDir v23 --logLevel DEBUG --sample fwlite_WZ_lep_LO_order2_15weights --interpolationOrder 2  --delphes --delphesEra RunII #SPLIT200
#python genPostProcessing.py --overwrite target --targetDir v23 --logLevel DEBUG --sample fwlite_tZq_LO_order2_15weights --interpolationOrder 2  --delphes --delphesEra RunII #SPLIT200
#python genPostProcessing.py --overwrite target --targetDir v23 --logLevel DEBUG --sample fwlite_tWZ_LO_order2_15weights --interpolationOrder 2  --delphes --delphesEra RunII #SPLIT200
#python genPostProcessing.py --overwrite target --targetDir v23 --logLevel DEBUG --sample fwlite_ttgamma_bg_LO_order2_15weights --interpolationOrder 2  --delphes --delphesEra RunII #SPLIT200



# ttgamma Signal samples with reference point
#python genPostProcessing.py --overwrite target --targetDir v18 --removeDelphesFiles --logLevel DEBUG --sample fwlite_ttgamma_LO_order2_15weights_ref --addReweights --interpolationOrder 2  --delphesCard delphes_card_CMS #SPLIT200
# large sample
#python genPostProcessing.py --overwrite target --targetDir v18 --logLevel DEBUG --sample fwlite_ttgammaLarge_LO_order2_15weights_ref --addReweights --interpolationOrder 2  --delphes --delphesEra RunII #SPLIT2000

# ttgamma background samples without reference point
#python genPostProcessing.py --overwrite target --targetDir v18 --logLevel DEBUG --sample fwlite_Zgamma_LO_order2_15weights --interpolationOrder 2  --delphes --delphesEra RunII #SPLIT200
#python genPostProcessing.py --overwrite target --targetDir v20 --logLevel DEBUG --sample fwlite_tt_dilep_LO_order2_15weights --interpolationOrder 2  --delphes --delphesEra RunII #SPLIT2200
#python genPostProcessing.py --overwrite target --targetDir v20 --logLevel DEBUG --sample fwlite_tt_nonhad_LO_order2_15weights --interpolationOrder 2  --delphes --delphesEra RunII #SPLIT200
#python genPostProcessing.py --overwrite target --targetDir v20 --logLevel DEBUG --sample fwlite_tW_LO_order2_15weights --interpolationOrder 2  --delphes --delphesEra RunII #SPLIT200
# large sample
#python genPostProcessing.py --overwrite target --targetDir v18 --logLevel DEBUG --sample fwlite_tt_full_LO_order2_15weights --interpolationOrder 2  --delphes --delphesEra RunII #SPLIT2000

## backgrounds according to Lukas:
#python genPostProcessing.py --targetDir RunII_v03 --logLevel DEBUG --interpolationOrder 2  --delphes --delphesEra RunII --sample fwlite_tt_full_LO_order2_15weights #SPLIT200
#python genPostProcessing.py --targetDir RunII_v03 --logLevel DEBUG --interpolationOrder 2  --delphes --delphesEra RunII --sample fwlite_WZ_lep_LO_order2_15weights #SPLIT200
#python genPostProcessing.py --targetDir RunII_v03 --logLevel DEBUG --interpolationOrder 2  --delphes --delphesEra RunII --sample fwlite_tZq_LO_order2_15weights #SPLIT200
#python genPostProcessing.py --targetDir RunII_v03 --logLevel DEBUG --interpolationOrder 2  --delphes --delphesEra RunII --sample fwlite_tW_LO_order2_15weights #SPLIT200
#python genPostProcessing.py --targetDir RunII_v03 --logLevel DEBUG --interpolationOrder 2  --delphes --delphesEra RunII --sample fwlite_tWZ_LO_order2_15weights #SPLIT200
#python genPostProcessing.py --targetDir RunII_v03 --logLevel DEBUG --interpolationOrder 2  --delphes --delphesEra RunII --sample fwlite_Zgamma_LO_order2_15weights #SPLIT200
##python genPostProcessing.py --targetDir RunII_v03 --logLevel DEBUG --interpolationOrder 2  --delphes --delphesEra RunII --sample fwlite_ttZ_ll_LO_order3_8weights #SPLIT200
#python genPostProcessing.py --targetDir RunII_v03 --logLevel DEBUG --interpolationOrder 2  --delphes --delphesEra RunII --sample fwlite_WJetsToLNu_order2_15weights #SPLIT200
##python genPostProcessing.py --targetDir RunII_v03 --logLevel DEBUG --addReweights --interpolationOrder 2  --delphes --delphesEra RunII --sample fwlite_ttZ_ll_LO_order2_15weights_ref_ext #SPLIT1999
#python genPostProcessing.py --targetDir RunII_v03 --logLevel DEBUG --addReweights --interpolationOrder 2  --delphes --delphesEra RunII --sample fwlite_ttZ_ll_LO_order2_15weights_ref #SPLIT200
##python genPostProcessing.py --targetDir RunII_v03 --logLevel DEBUG --interpolationOrder 2  --delphes --delphesEra RunII --sample fwlite_ttW_LO_order3_8weights #SPLIT200
##python genPostProcessing.py --targetDir RunII_v03 --logLevel DEBUG --interpolationOrder 2  --delphes --delphesEra RunII --sample fwlite_ttgamma_bg_LO_order2_15weights #SPLIT200


#python genPostProcessing.py --targetDir RunII_v02 --logLevel DEBUG --sample  ttbarZ --HEPMC PP       --delphes --delphesEra RunII #SPLIT100
#python genPostProcessing.py --targetDir RunII_v02 --logLevel DEBUG --sample  ttbarZ --HEPMC 10d0_HH  --delphes --delphesEra RunII #SPLIT50
#python genPostProcessing.py --targetDir RunII_v02 --logLevel DEBUG --sample  ttbarZ --HEPMC 10d0_GH  --delphes --delphesEra RunII #SPLIT50
#python genPostProcessing.py --targetDir RunII_v02 --logLevel DEBUG --sample  ttbarZ --HEPMC 10d0_HG  --delphes --delphesEra RunII #SPLIT50
#python genPostProcessing.py --targetDir RunII_v02 --logLevel DEBUG --sample  ttbarZ --HEPMC 1d0_GH   --delphes --delphesEra RunII #SPLIT50
#python genPostProcessing.py --targetDir RunII_v02 --logLevel DEBUG --sample  ttbarZ --HEPMC 1d0_HG   --delphes --delphesEra RunII #SPLIT50
#python genPostProcessing.py --targetDir RunII_v02 --logLevel DEBUG --sample  ttbarZ --HEPMC 1d0_HH   --delphes --delphesEra RunII #SPLIT50
#python genPostProcessing.py --targetDir RunII_v02 --logLevel DEBUG --sample  ttbarZ --HEPMC 0d1_HG   --delphes --delphesEra RunII #SPLIT50
#python genPostProcessing.py --targetDir RunII_v02 --logLevel DEBUG --sample  ttbarZ --HEPMC 0d1_GH   --delphes --delphesEra RunII #SPLIT50
#python genPostProcessing.py --targetDir RunII_v02 --logLevel DEBUG --sample  ttbarZ --HEPMC 0d1_HH   --delphes --delphesEra RunII #SPLIT50

#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunII --HEPMC PP #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunII --HEPMC pdf-0.01_HG #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunII --HEPMC pdf-0.1_HG #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunII --HEPMC pdf-1_HG #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunII --HEPMC pdf-10_HG #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunII --HEPMC pdf-30_HG #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunII --HEPMC pdf-35_HG #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunII --HEPMC pdf-37_HG #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunII --HEPMC pdf-40_HG #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunII --HEPMC pdf-42_HG #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunII --HEPMC pdf-43_HG #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunII --HEPMC pdf-44_HG #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunII --HEPMC pdf-45_HG #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunII --HEPMC pdf-50_HG #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunII --HEPMC pdf-75_HG #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunII --HEPMC pdf-100_GH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunII --HEPMC pdf-0.01_GH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunII --HEPMC pdf-0.1_GH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunII --HEPMC pdf-1_GH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunII --HEPMC pdf-10_GH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunII --HEPMC pdf-30_GH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunII --HEPMC pdf-35_GH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunII --HEPMC pdf-37_GH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunII --HEPMC pdf-40_GH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunII --HEPMC pdf-42_GH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunII --HEPMC pdf-43_GH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunII --HEPMC pdf-44_GH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunII --HEPMC pdf-45_GH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunII --HEPMC pdf-50_GH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunII --HEPMC pdf-75_GH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunII --HEPMC pdf-100_GH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunII --HEPMC pdf-0.01_HH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunII --HEPMC pdf-0.1_HH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunII --HEPMC pdf-1_HH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunII --HEPMC pdf-10_HH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunII --HEPMC pdf-30_HH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunII --HEPMC pdf-35_HH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunII --HEPMC pdf-37_HH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunII --HEPMC pdf-40_HH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunII --HEPMC pdf-42_HH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunII --HEPMC pdf-43_HH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunII --HEPMC pdf-44_HH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunII --HEPMC pdf-45_HH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunII --HEPMC pdf-50_HH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunII --HEPMC pdf-75_HH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunII --HEPMC pdf-100_HH #SPLIT5
#
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC PP #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-0.01_HG #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-0.1_HG #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-1_HG #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-10_HG #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-30_HG #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-35_HG #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-37_HG #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-40_HG #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-42_HG #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-43_HG #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-44_HG #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-45_HG #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-50_HG #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-75_HG #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-100_GH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-0.01_GH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-0.1_GH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-1_GH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-10_GH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-30_GH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-35_GH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-37_GH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-40_GH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-42_GH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-43_GH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-44_GH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-45_GH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-50_GH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-75_GH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-100_GH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-0.01_HH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-0.1_HH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-1_HH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-10_HH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-30_HH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-35_HH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-37_HH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-40_HH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-42_HH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-43_HH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-44_HH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-45_HH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-50_HH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-75_HH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v02/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-100_HH #SPLIT5

#python genPostProcessing.py --targetDir RunII_v03 --logLevel DEBUG --sample  ttbarZ --HEPMC PP       --delphes --delphesEra RunII #SPLIT100
#python genPostProcessing.py --targetDir RunII_v03 --logLevel DEBUG --sample  ttbarZ --HEPMC 10d0_HH  --delphes --delphesEra RunII #SPLIT50
#python genPostProcessing.py --targetDir RunII_v03 --logLevel DEBUG --sample  ttbarZ --HEPMC 10d0_GH  --delphes --delphesEra RunII #SPLIT50
#python genPostProcessing.py --targetDir RunII_v03 --logLevel DEBUG --sample  ttbarZ --HEPMC 10d0_HG  --delphes --delphesEra RunII #SPLIT50
#python genPostProcessing.py --targetDir RunII_v03 --logLevel DEBUG --sample  ttbarZ --HEPMC 1d0_GH   --delphes --delphesEra RunII #SPLIT50
#python genPostProcessing.py --targetDir RunII_v03 --logLevel DEBUG --sample  ttbarZ --HEPMC 1d0_HG   --delphes --delphesEra RunII #SPLIT50
#python genPostProcessing.py --targetDir RunII_v03 --logLevel DEBUG --sample  ttbarZ --HEPMC 1d0_HH   --delphes --delphesEra RunII #SPLIT50
#python genPostProcessing.py --targetDir RunII_v03 --logLevel DEBUG --sample  ttbarZ --HEPMC 0d1_HG   --delphes --delphesEra RunII #SPLIT50
#python genPostProcessing.py --targetDir RunII_v03 --logLevel DEBUG --sample  ttbarZ --HEPMC 0d1_GH   --delphes --delphesEra RunII #SPLIT50
#python genPostProcessing.py --targetDir RunII_v03 --logLevel DEBUG --sample  ttbarZ --HEPMC 0d1_HH   --delphes --delphesEra RunII #SPLIT50

#python genPostProcessing.py --targetDir PhaseII/22_05 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunII --HEPMC PP --nJobs 5 --job 0 &
#python genPostProcessing.py --targetDir PhaseII/22_05 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunII --HEPMC PP --nJobs 5 --job 1 &
#python genPostProcessing.py --targetDir PhaseII/22_05 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunII --HEPMC PP --nJobs 5 --job 2 &
#python genPostProcessing.py --targetDir PhaseII/22_05 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunII --HEPMC PP --nJobs 5 --job 3 &
#python genPostProcessing.py --targetDir PhaseII/22_05 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunII --HEPMC PP --nJobs 5 --job 4 &

#python genPostProcessing.py --targetDir RunII_v03/22_05 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunII --HEPMC pdf-test_HG #SPLIT5
#python genPostProcessing.py --targetDir RunII_v03/22_05 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunII --HEPMC pdf-test_GH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v03/22_05 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunII --HEPMC pdf-test_HH #SPLIT5

#python genPostProcessing.py --targetDir RunII_v03/22_05 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunII --HEPMC pdf-PP_HG #SPLIT5
#python genPostProcessing.py --targetDir RunII_v03/22_05 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunII --HEPMC pdf-PP_GH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v03/22_05 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunII --HEPMC pdf-PP_HH #SPLIT5

#python genPostProcessing.py --targetDir RunII_v03/22_05 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunII --HEPMC pdf-PPold_HG #SPLIT20
#python genPostProcessing.py --targetDir RunII_v03/22_05 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunII --HEPMC pdf-PPold_GH #SPLIT20
#python genPostProcessing.py --targetDir RunII_v03/22_05 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunII --HEPMC pdf-PPold_HH #SPLIT20
#
#python genPostProcessing.py --targetDir RunII_v03/22_05 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunII --HEPMC pdf-PPLO_HG #SPLIT20
#python genPostProcessing.py --targetDir RunII_v03/22_05 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunII --HEPMC pdf-PPLO_GH #SPLIT20
#python genPostProcessing.py --targetDir RunII_v03/22_05 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunII --HEPMC pdf-PPLO_HH #SPLIT20

#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/11_06 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunIInoDelphesIso --HEPMC PP #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/11_06 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-0.01_HG #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/11_06 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-0.1_HG #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/11_06 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-1_HG #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/11_06 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-10_HG #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/11_06 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-30_HG #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/11_06 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-35_HG #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/11_06 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-37_HG #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/11_06 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-40_HG #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/11_06 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-42_HG #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/11_06 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-43_HG #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/11_06 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-44_HG #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/11_06 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-45_HG #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/11_06 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-50_HG #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/11_06 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-75_HG #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/11_06 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-100_HG #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/11_06 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-0.01_GH #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/11_06 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-0.1_GH #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/11_06 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-1_GH #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/11_06 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-10_GH #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/11_06 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-30_GH #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/11_06 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-35_GH #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/11_06 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-37_GH #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/11_06 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-40_GH #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/11_06 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-42_GH #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/11_06 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-43_GH #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/11_06 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-44_GH #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/11_06 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-45_GH #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/11_06 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-50_GH #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/11_06 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-75_GH #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/11_06 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-100_GH #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/11_06 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-0.01_HH #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/11_06 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-0.1_HH #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/11_06 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-1_HH #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/11_06 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-10_HH #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/11_06 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-30_HH #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/11_06 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-35_HH #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/11_06 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-37_HH #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/11_06 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-40_HH #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/11_06 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-42_HH #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/11_06 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-43_HH #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/11_06 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-44_HH #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/11_06 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-45_HH #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/11_06 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-50_HH #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/11_06 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-75_HH #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/11_06 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-100_HH #SPLIT5
#
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/24_09 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunIInoDelphesIso --HEPMC PP #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/24_09 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-0.01_HG #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/24_09 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-0.1_HG #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/24_09 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-1_HG #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/24_09 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-10_HG #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/24_09 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-30_HG #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/24_09 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-35_HG #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/24_09 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-37_HG #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/24_09 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-40_HG #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/24_09 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-42_HG #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/24_09 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-43_HG #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/24_09 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-44_HG #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/24_09 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-45_HG #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/24_09 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-50_HG #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/24_09 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-75_HG #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/24_09 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-100_HG #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/24_09 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-0.01_GH #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/24_09 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-0.1_GH #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/24_09 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-1_GH #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/24_09 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-10_GH #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/24_09 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-30_GH #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/24_09 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-35_GH #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/24_09 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-37_GH #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/24_09 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-40_GH #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/24_09 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-42_GH #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/24_09 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-43_GH #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/24_09 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-44_GH #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/24_09 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-45_GH #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/24_09 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-50_GH #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/24_09 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-75_GH #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/24_09 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-100_GH #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/24_09 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-0.01_HH #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/24_09 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-0.1_HH #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/24_09 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-1_HH #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/24_09 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-10_HH #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/24_09 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-30_HH #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/24_09 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-35_HH #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/24_09 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-37_HH #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/24_09 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-40_HH #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/24_09 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-42_HH #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/24_09 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-43_HH #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/24_09 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-44_HH #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/24_09 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-45_HH #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/24_09 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-50_HH #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/24_09 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-75_HH #SPLIT5
#python genPostProcessing.py --targetDir RunIInoDelphesIso_v03/24_09 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunIInoDelphesIso --HEPMC pdf-100_HH #SPLIT5

#
#python genPostProcessing.py --targetDir RunII_v03/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC PP #SPLIT5
#python genPostProcessing.py --targetDir RunII_v03/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-0.01_HG #SPLIT5
#python genPostProcessing.py --targetDir RunII_v03/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-0.1_HG #SPLIT5
#python genPostProcessing.py --targetDir RunII_v03/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-1_HG #SPLIT5
#python genPostProcessing.py --targetDir RunII_v03/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-10_HG #SPLIT5
#python genPostProcessing.py --targetDir RunII_v03/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-30_HG #SPLIT5
#python genPostProcessing.py --targetDir RunII_v03/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-35_HG #SPLIT5
#python genPostProcessing.py --targetDir RunII_v03/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-37_HG #SPLIT5
#python genPostProcessing.py --targetDir RunII_v03/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-40_HG #SPLIT5
#python genPostProcessing.py --targetDir RunII_v03/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-42_HG #SPLIT5
#python genPostProcessing.py --targetDir RunII_v03/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-43_HG #SPLIT5
#python genPostProcessing.py --targetDir RunII_v03/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-44_HG #SPLIT5
#python genPostProcessing.py --targetDir RunII_v03/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-45_HG #SPLIT5
#python genPostProcessing.py --targetDir RunII_v03/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-50_HG #SPLIT5
#python genPostProcessing.py --targetDir RunII_v03/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-75_HG #SPLIT5
#python genPostProcessing.py --targetDir RunII_v03/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-100_HH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v03/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-0.01_GH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v03/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-0.1_GH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v03/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-1_GH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v03/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-10_GH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v03/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-30_GH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v03/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-35_GH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v03/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-37_GH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v03/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-40_GH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v03/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-42_GH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v03/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-43_GH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v03/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-44_GH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v03/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-45_GH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v03/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-50_GH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v03/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-75_GH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v03/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-100_GH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v03/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-0.01_HH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v03/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-0.1_HH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v03/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-1_HH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v03/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-10_HH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v03/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-30_HH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v03/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-35_HH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v03/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-37_HH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v03/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-40_HH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v03/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-42_HH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v03/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-43_HH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v03/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-44_HH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v03/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-45_HH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v03/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-50_HH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v03/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-75_HH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v03/22_05 --logLevel DEBUG --sample  ttbarZ --delphes --delphesEra RunII --HEPMC pdf-100_HH #SPLIT5
#
#python genPostProcessing.py --targetDir RunII_v03_noDelphesIso/13_09 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunII --HEPMC PP #SPLIT5
#python genPostProcessing.py --targetDir RunII_v03_noDelphesIso/13_09 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunII --HEPMC pdf-x_HG #SPLIT5
#python genPostProcessing.py --targetDir RunII_v03_noDelphesIso/13_09 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunII --HEPMC pdf-x_GH #SPLIT5
#python genPostProcessing.py --targetDir RunII_v03_noDelphesIso/13_09 --logLevel DEBUG --sample  ttbar --delphes --delphesEra RunII --HEPMC pdf-x_HH #SPLIT5

#python genPostProcessing.py --targetDir RunII_v03 --logLevel DEBUG --interpolationOrder 2  --delphes --delphesEra RunII --sample fwlite_tt_full_LO_order2_15weights #SPLIT200

#python genPostProcessing.py --logLevel DEBUG --sample --interpolationOrder 4 --delphes --delphesEra RunII --sample #SPLIT2000

