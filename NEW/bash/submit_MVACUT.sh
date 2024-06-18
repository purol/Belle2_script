#!/bin/bash
# submit_MVACUT.sh {Analysis bin file} {directory name} {train|test|validation} {CHG_analysis|MIX_analysis|UUBAR_analysis|DDBAR_analysis|SSBAR_analysis|CHARM_analysis|SIGNAL_analysis} {version name} {MC version}
# submit_MVACUT.sh ./bin/Analysis_Main test_v001 test CHG_analysis Kogasa MC15ri

export HOME=/cvmfs/belle.cern.ch/el7/externals/v02-00-01/Linux_x86_64/common
export PATH=$HOME/bin:$PATH
export LDFLAGs="-Wl,-rpath,$HOME/bin"
export LD_LIBRARY_PATH=$HOME/lib64
source $HOME/root/bin/thisroot.sh
export ROOTLIB=$HOME/root/lib/root
export ROOTSYS=$HOME/root/lib/root

binfile=$1
VERNAME=$5

if [ $4 = "CHG_analysis" ]; then
  sample="CHG"
fi
if [ $4 = "MIX_analysis" ]; then
  sample="MIX"
fi
if [ $4 = "UUBAR_analysis" ]; then
  sample="UUBAR"
fi
if [ $4 = "DDBAR_analysis" ]; then
  sample="DDBAR"
fi
if [ $4 = "SSBAR_analysis" ]; then
  sample="SSBAR"
fi
if [ $4 = "CHARM_analysis" ]; then
  sample="CHARM"
fi
if [ $4 = "SIGNAL_analysis" ]; then
  sample="SIGNAL"
fi


DIRNAME="/home/belle2/junewoo/storage_b1/bsub/Analysis/$VERNAME/$4/$2/final_output_root_after_MVA_Application" # path of input file
OUTPUTPATH="/home/belle2/junewoo/storage_b1/bsub/Analysis/$VERNAME/$4/$2" # path of output files

${binfile} ${DIRNAME} ${OUTPUTPATH} ${sample} $6 $3

