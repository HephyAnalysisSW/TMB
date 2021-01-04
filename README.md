# TTXPheno
## Get code for TTXPheno studies

```
cmsrel CMSSW_10_6_0
cd CMSSW_10_6_0/src
cmsenv
git cms-init
git clone https://github.com/schoef/TMB
git clone https://github.com/HephyAnalysisSW/RootTools
scram b -j40
```

## Delphes
```
cd $CMSSW_BASE/..
git clone https://github.com/TTXPheno/delphes.git
#patch $CMSSW_BASE/../delphes/cards/delphes_card_CMS.tcl < $CMSSW_BASE/src/TTXPheno/patches/slim_delphes.diff # Reduce Delphes output
cd delphes
./configure
sed -i -e 's/c++0x/c++17/g' Makefile
make -j 4 
```
