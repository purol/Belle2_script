#!/bin/bsub
# ./submiter.sh ./submit_ReadRootFile.sh {root file path} {output destination} 

# main analysis MC
code="Analysis_main"
MCver="MC15ri"
VerName="Kumoi"
Types=("SIGNAL_analysis" "CHG_analysis" "MIX_analysis" "UUBAR_analysis" "DDBAR_analysis" "SSBAR_analysis" "CHARM_analysis")

DirName="test_v007"
SimpleDirName="test"
for Type in "${Types[@]}"; do # for test sample
    mkdir -p "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}"
    bsub -o "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/${Type}_${DirName}_${code}.log" -e "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/${Type}_${DirName}_${code}_err.txt" "./submit_Analysis.sh" "./bin/${code}" ${DirName} ${SimpleDirName} ${Type} ${VerName} ${MCver}
    sleep 0.1s
done

DirName="train_v007"
SimpleDirName="train"
for Type in "${Types[@]}"; do # for train sample
    mkdir -p "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}"
    bsub -o "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/${Type}_${DirName}_${code}.log" -e "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/${Type}_${DirName}_${code}_err.txt" "./submit_Analysis.sh" "./bin/${code}" ${DirName} ${SimpleDirName} ${Type} ${VerName} ${MCver}
    sleep 0.1s
done


# off-resonance data
code="Analysis_off"
MCver="data"
VerName="Kumoi_LS_data_off"
Types=("SIGNAL_analysis")

DirName="validation_v007"
SimpleDirName="validation"
for Type in "${Types[@]}"; do
    mkdir -p "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}"
    bsub -o "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/${Type}_${DirName}_${code}.log" -e "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/${Type}_${DirName}_${code}_err.txt" "./submit_Analysis.sh" "./bin/${code}" ${DirName} ${SimpleDirName} ${Type} ${VerName} ${MCver}
    sleep 0.1s
done


# off-resonance MC
code="Analysis_main"
MCver="MC15ri"
VerName="Kumoi_LS_MC_off"
Types=("UUBAR_analysis" "DDBAR_analysis" "SSBAR_analysis" "CHARM_analysis")

DirName="train_v007"
SimpleDirName="train"
for Type in "${Types[@]}"; do 
    mkdir -p "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}"
    bsub -o "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/${Type}_${DirName}_${code}.log" -e "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/${Type}_${DirName}_${code}_err.txt" "./submit_Analysis.sh" "./bin/${code}" ${DirName} ${SimpleDirName} ${Type} ${VerName} ${MCver}
    sleep 0.1s
done

DirName="test_v007"
SimpleDirName="test"
for Type in "${Types[@]}"; do 
    mkdir -p "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}"
    bsub -o "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/${Type}_${DirName}_${code}.log" -e "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/${Type}_${DirName}_${code}_err.txt" "./submit_Analysis.sh" "./bin/${code}" ${DirName} ${SimpleDirName} ${Type} ${VerName} ${MCver}
    sleep 0.1s
done


# syst
code="Analysis_main"
MCver="MC15ri"
VerName="Kumoi_syst"
Types=("SIGNAL_analysis")

DirName="validation_v007"
SimpleDirName="validation"
for Type in "${Types[@]}"; do
    mkdir -p "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}"
    bsub -o "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/${Type}_${DirName}_${code}.log" -e "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/${Type}_${DirName}_${code}_err.txt" "./submit_Analysis.sh" "./bin/${code}" ${DirName} ${SimpleDirName} ${Type} ${VerName} ${MCver}
    sleep 0.1s
done


# cside data
code="Analysis_main"
MCver="data"
VerName="Kumoi_LS_data_cside"
Types=("SIGNAL_analysis")

DirName="validation_v007"
SimpleDirName="validation"
for Type in "${Types[@]}"; do
    mkdir -p "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}"
    bsub -o "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/${Type}_${DirName}_${code}.log" -e "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/${Type}_${DirName}_${code}_err.txt" "./submit_Analysis.sh" "./bin/${code}" ${DirName} ${SimpleDirName} ${Type} ${VerName} ${MCver}
    sleep 0.1s
done


# cside MC
code="Analysis_main"
MCver="MC15ri"
VerName="Kumoi_LS_MC_cside"
Types=("CHG_analysis" "MIX_analysis" "UUBAR_analysis" "DDBAR_analysis" "SSBAR_analysis" "CHARM_analysis")

DirName="test_v007"
SimpleDirName="test"
for Type in "${Types[@]}"; do # for test sample
    mkdir -p "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}"
    bsub -o "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/${Type}_${DirName}_${code}.log" -e "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/${Type}_${DirName}_${code}_err.txt" "./submit_Analysis.sh" "./bin/${code}" ${DirName} ${SimpleDirName} ${Type} ${VerName} ${MCver}
    sleep 0.1s
done

DirName="train_v007"
SimpleDirName="train"
for Type in "${Types[@]}"; do # for train sample
    mkdir -p "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}"
    bsub -o "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/${Type}_${DirName}_${code}.log" -e "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/${Type}_${DirName}_${code}_err.txt" "./submit_Analysis.sh" "./bin/${code}" ${DirName} ${SimpleDirName} ${Type} ${VerName} ${MCver}
    sleep 0.1s
done

