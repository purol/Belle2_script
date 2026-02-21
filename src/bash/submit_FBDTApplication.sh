#!/bin/bash
# submit_FBDTApplication.sh {Analysis bin file} {dirname} {weightfile version name} {simple dirname}
# submit_FBDTApplication.sh ./bin/FastBDTClassificationApplication /home/belle2/junewoo/storage_b1/bsub/Analysis/Satori_LS_MC_side/CHG_analysis/test_v000 Satori v000 

export HOME=/cvmfs/belle.cern.ch/el7/externals/v02-00-01/Linux_x86_64/common
export PATH=$HOME/bin:$PATH
export LDFLAGs="-Wl,-rpath,$HOME/bin"
export LD_LIBRARY_PATH=$HOME/lib64
source $HOME/root/bin/thisroot.sh
export ROOTLIB=$HOME/root/lib/root
export ROOTSYS=$HOME/root/lib/root

binfile=$1

${binfile} ${2} ${3} ${4}

