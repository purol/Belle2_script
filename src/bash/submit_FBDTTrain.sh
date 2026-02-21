#!/bin/bash
# submit_FBDTTrain.sh {bin name} {nTree} {depth} {shrinkage path} {subsample type} {binning} {version name} {dirname} {MC version}
# submit_FBDTTrain.sh ./bin/FastBDTGridSearch 500 3 0.1 0.7 7 Satori v000 MC15ri

export HOME=/cvmfs/belle.cern.ch/el7/externals/v02-00-01/Linux_x86_64/common
export PATH=$HOME/bin:$PATH
export LDFLAGs="-Wl,-rpath,$HOME/bin"
export LD_LIBRARY_PATH=$HOME/lib64
source $HOME/root/bin/thisroot.sh
export ROOTLIB=$HOME/root/lib/root
export ROOTSYS=$HOME/root/lib/root

$1 $2 $3 $4 $5 $6 $7 $8 $9
