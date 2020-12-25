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

## Prerequisites
Add your user-specific locations to [user.py](https://github.com/TTXPheno/TTXPheno/blob/master/Tools/python/user.py).
You can do that also later and thus fix potential 'ImportError' from `user.py`  when they occur. 

## How to make a flat ntuple from a GEN file
Define the sample in e.g. [fwlite_benchmarks.py](https://github.com/TTXPheno/TTXPheno/blob/master/samples/python/fwlite_benchmarks.py) and do 
```
voms-proxy-init -voms cms -out ~/private/.proxy
cd TTXPheno/postprocessing
python genPostProcessing.py --overwrite --sample fwlite_ttZ_ll_LO_highStat_scan --addReweights --small
```
Remove `--small` to process the full sample. Add `--nJobs=10` and `--job=0` etc. options to run only on the first 10% of events. In Vienna you could also do 
```
export X509_USER_PROXY=~/private/.proxy
submitBatch.py --dpm genPostProcessing.sh # append `#SPLIT200` to a line
``` 
For other places we can add submission scripts.

## Examples
Polynomial parametrization:
```
python Tools/python/HyperPoly.py
```
