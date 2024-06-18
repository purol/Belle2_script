#!/bin/bash
# submit_Deconvertor.sh {bin name} {input file path} {OBB cut} {Oqq cut} {output file path}
# submiter_Deconvertor.sh ./bin/Deconvertor /home/belle2/junewoo/storage_b1/bsub/Analysis/Satori_LS_data_off/SIGNAL_analysis/test_v000/final_output_data_after_MVA_Application -0.1 -0.1 /home/belle2/junewoo/storage_b1/bsub/Analysis/Satori_LS_data_off/SIGNAL_analysis/test_v000/final_output_root_after_MVA_Application

export HOME=/cvmfs/belle.cern.ch/el7/externals/v02-00-01/Linux_x86_64/common
export PATH=$HOME/bin:$PATH
export LDFLAGs="-Wl,-rpath,$HOME/bin"
export LD_LIBRARY_PATH=$HOME/lib64
source $HOME/root/bin/thisroot.sh
export ROOTLIB=$HOME/root/lib/root
export ROOTSYS=$HOME/root/lib/root

binfile=$1

${binfile} ${2} ${3} ${4} ${5}

