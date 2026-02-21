#!/bin/bsub
# ./submiter.sh ./submit_ReadRootFile.sh {root file path} {output destination} 

mkdir -p /home/belle2/junewoo/storage_b1/bsub/Analysis/GridSearch_BDTc/log/
mkdir -p /home/belle2/junewoo/storage_b1/bsub/Analysis/GridSearch_BDTc/err/
mkdir -p /home/belle2/junewoo/storage_b1/bsub/Analysis/GridSearch_BDTc/out/


# train off-resonance Ntuple
code="FastBDTGridSearch_BDTc_off_RD"
VerName="KumoiRD"
DirName="v007"
MCver="MC15rd"

for nTree in 500 1000 1500
do
  for depth in 1 2 3 4
  do
    for shrinkage in 0.01 0.05 0.1
    do
      for subsample in 0.01 0.3 0.4 0.5 0.6 0.7
      do
        for binning in 5 6 7 8 9
        do
          bsub -q s -o /home/belle2/junewoo/storage_b1/bsub/Analysis/GridSearch_BDTc/log/${nTree}_${depth}_${shrinkage}_${subsample}_${binning}.log -e /home/belle2/junewoo/storage_b1/bsub/Analysis/GridSearch_BDTc/err/${nTree}_${depth}_${shrinkage}_${subsample}_${binning}.err ./submit_FBDTTrain.sh "./bin/"${code} ${nTree} ${depth} ${shrinkage} ${subsample} ${binning} ${VerName} ${DirName} ${MCver}
        done
      done
    done
  done
done

for nTree in 2000
do
  for depth in 1 2 3 4
  do
    for shrinkage in 0.01 0.05 0.1
    do
      for subsample in 0.01 0.3 0.4 0.5 0.6 0.7
      do
        for binning in 5 6 7 8 9
        do
          bsub -q s -o /home/belle2/junewoo/storage_b1/bsub/Analysis/GridSearch_BDTc/log/${nTree}_${depth}_${shrinkage}_${subsample}_${binning}.log -e /home/belle2/junewoo/storage_b1/bsub/Analysis/GridSearch_BDTc/err/${nTree}_${depth}_${shrinkage}_${subsample}_${binning}.err ./submit_FBDTTrain.sh "./bin/"${code} ${nTree} ${depth} ${shrinkage} ${subsample} ${binning} ${VerName} ${DirName} ${MCver}
        done
      done
    done
  done
done
