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
export PYTHIA8=/cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/pythia8/223-mlhled2/
export LD_LIBRARY_PATH=$PYTHIA8/lib:$LD_LIBRARY_PATH
cd $CMSSW_BASE/..
git clone https://github.com/TTXPheno/delphes.git
#patch $CMSSW_BASE/../delphes/cards/delphes_card_CMS.tcl < $CMSSW_BASE/src/TTXPheno/patches/slim_delphes.diff # Reduce Delphes output
cd delphes
./configure
sed -i -e 's/c++0x/c++1y/g' Makefile
make -j 4 
```
