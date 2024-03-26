#!/bin/bsub
# ./submiter.sh ./submit_ReadRootFile.sh {root file path} {output destination} 

# main analysis
code="Analysis_main"
MCver="MC15rd"
VerName="KumoiRD"
Types=("SIGNAL_analysis" "CHG_analysis" "MIX_analysis" "UUBAR_analysis" "DDBAR_analysis" "SSBAR_analysis" "CHARM_analysis")

DirName="validation_v000"
SimpleDirName="validation"
for Type in "${Types[@]}"; do # for validation sample
    mkdir -p "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}"
    bsub -q l -o "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/${Type}_${DirName}_${code}.log" -e "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/${Type}_${DirName}_${code}_err.txt" "./submit_Analysis.sh" "./bin/${code}" ${DirName} ${SimpleDirName} ${Type} ${VerName} ${MCver}
    sleep 0.1s
done


# off-resonance data
code="Analysis_off"
MCver="data"
VerName="KumoiRD_LS_data_off"
Types=("SIGNAL_analysis")

DirName="validation_v000"
SimpleDirName="validation"
for Type in "${Types[@]}"; do
    mkdir -p "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}"
    bsub -q s -o "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/${Type}_${DirName}_${code}.log" -e "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/${Type}_${DirName}_${code}_err.txt" "./submit_Analysis.sh" "./bin/${code}" ${DirName} ${SimpleDirName} ${Type} ${VerName} ${MCver}
    sleep 0.1s
done


# off-resonance MC
code="Analysis_off"
MCver="MC15rd"
VerName="KumoiRD_LS_MC_off"
Types=("UUBAR_analysis" "DDBAR_analysis" "SSBAR_analysis" "CHARM_analysis")

DirName="validation_v000"
SimpleDirName="validation"
for Type in "${Types[@]}"; do 
    mkdir -p "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}"
    bsub -q s -o "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/${Type}_${DirName}_${code}.log" -e "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/${Type}_${DirName}_${code}_err.txt" "./submit_Analysis.sh" "./bin/${code}" ${DirName} ${SimpleDirName} ${Type} ${VerName} ${MCver}
    sleep 0.1s
done


# sideband data
code="Analysis_side"
MCver="MC15rd"
VerName="KumoiRD_LS_data_side"
Types=("SIGNAL_analysis") 

DirName="validation_v000"
SimpleDirName="validation"
for Type in "${Types[@]}"; do 
    mkdir -p "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}"
    bsub -q s -o "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/${Type}_${DirName}_${code}.log" -e "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/${Type}_${DirName}_${code}_err.txt" "./submit_Analysis.sh" "./bin/${code}" ${DirName} ${SimpleDirName} ${Type} ${VerName} ${MCver}
    sleep 0.1s
done


## sideband MC
code="Analysis_side"
MCver="MC15rd"
VerName="KumoiRD_LS_MC_side"
Types=("CHG_analysis" "MIX_analysis" "UUBAR_analysis" "DDBAR_analysis" "SSBAR_analysis" "CHARM_analysis")

DirName="validation_v000"
SimpleDirName="validation"
for Type in "${Types[@]}"; do 
    mkdir -p "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}"
    bsub -q s -o "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/${Type}_${DirName}_${code}.log" -e "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/${Type}_${DirName}_${code}_err.txt" "./submit_Analysis.sh" "./bin/${code}" ${DirName} ${SimpleDirName} ${Type} ${VerName} ${MCver}
    sleep 0.1s
done


# EL sideband data
code="Analysis_main"
MCver="MC15rd"
VerName="KumoiRD_LS_data_EL"
Types=("SIGNAL_analysis")

DirName="validation_v000"
SimpleDirName="validation"
for Type in "${Types[@]}"; do
    mkdir -p "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}"
    bsub -q s -o "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/${Type}_${DirName}_${code}.log" -e "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/${Type}_${DirName}_${code}_err.txt" "./submit_Analysis.sh" "./bin/${code}" ${DirName} ${SimpleDirName} ${Type} ${VerName} ${MCver}
    sleep 0.1s
done


## EL sideband MC
code="Analysis_main"
MCver="MC15rd"
VerName="KumoiRD_LS_MC_EL"
Types=("CHG_analysis" "MIX_analysis" "UUBAR_analysis" "DDBAR_analysis" "SSBAR_analysis" "CHARM_analysis")

DirName="validation_v000"
SimpleDirName="validation"
for Type in "${Types[@]}"; do
    mkdir -p "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}"
    bsub -q s -o "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/${Type}_${DirName}_${code}.log" -e "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/${Type}_${DirName}_${code}_err.txt" "./submit_Analysis.sh" "./bin/${code}" ${DirName} ${SimpleDirName} ${Type} ${VerName} ${MCver}
    sleep 0.1s
done


# MU sideband data
code="Analysis_main"
MCver="MC15rd"
VerName="KumoiRD_LS_data_MU"
Types=("SIGNAL_analysis")

DirName="validation_v000"
SimpleDirName="validation"
for Type in "${Types[@]}"; do
    mkdir -p "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}"
    bsub -q s -o "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/${Type}_${DirName}_${code}.log" -e "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/${Type}_${DirName}_${code}_err.txt" "./submit_Analysis.sh" "./bin/${code}" ${DirName} ${SimpleDirName} ${Type} ${VerName} ${MCver}
    sleep 0.1s
done


## MU sideband MC
code="Analysis_main"
MCver="MC15rd"
VerName="KumoiRD_LS_MC_MU"
Types=("CHG_analysis" "MIX_analysis" "UUBAR_analysis" "DDBAR_analysis" "SSBAR_analysis" "CHARM_analysis")

DirName="validation_v000"
SimpleDirName="validation"
for Type in "${Types[@]}"; do
    mkdir -p "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}"
    bsub -q s -o "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/${Type}_${DirName}_${code}.log" -e "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/${Type}_${DirName}_${code}_err.txt" "./submit_Analysis.sh" "./bin/${code}" ${DirName} ${SimpleDirName} ${Type} ${VerName} ${MCver}
    sleep 0.1s
done


# Jpsi data
code="Analysis_main"
MCver="MC15rd"
VerName="KumoiRD_LS_data_Jpsi"
Types=("SIGNAL_analysis")

DirName="validation_v000"
SimpleDirName="validation"
for Type in "${Types[@]}"; do
    mkdir -p "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}"
    bsub -q s -o "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/${Type}_${DirName}_${code}.log" -e "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/${Type}_${DirName}_${code}_err.txt" "./submit_Analysis.sh" "./bin/${code}" ${DirName} ${SimpleDirName} ${Type} ${VerName} ${MCver}
    sleep 0.1s
done


# Jpsi MC
code="Analysis_main"
MCver="MC15rd"
VerName="KumoiRD_LS_MC_Jpsi"
Types=("SIGNAL_analysis" "CHG_analysis" "MIX_analysis" "UUBAR_analysis" "DDBAR_analysis" "SSBAR_analysis" "CHARM_analysis")

DirName="validation_v000"
SimpleDirName="validation"
for Type in "${Types[@]}"; do
    mkdir -p "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}"
    bsub -q s -o "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/${Type}_${DirName}_${code}.log" -e "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/${Type}_${DirName}_${code}_err.txt" "./submit_Analysis.sh" "./bin/${code}" ${DirName} ${SimpleDirName} ${Type} ${VerName} ${MCver}
    sleep 0.1s
done


# embedded 
code="Analysis_main"
MCver="MC15rd"
VerName="KumoiRD_embedded"
Types=("SIGNAL_analysis")

DirName="validation_v000"
SimpleDirName="validation"
for Type in "${Types[@]}"; do
    mkdir -p "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}"
    bsub -q s -o "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/${Type}_${DirName}_${code}.log" -e "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/${Type}_${DirName}_${code}_err.txt" "./submit_Analysis.sh" "./bin/${code}" ${DirName} ${SimpleDirName} ${Type} ${VerName} ${MCver}
    sleep 0.1s
done


# cside data
code="Analysis_main"
MCver="MC15rd"
VerName="KumoiRD_LS_data_cside"
Types=("SIGNAL_analysis")

DirName="validation_v000"
SimpleDirName="validation"
for Type in "${Types[@]}"; do
    mkdir -p "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}"
    bsub -q s -o "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/${Type}_${DirName}_${code}.log" -e "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/${Type}_${DirName}_${code}_err.txt" "./submit_Analysis.sh" "./bin/${code}" ${DirName} ${SimpleDirName} ${Type} ${VerName} ${MCver}
    sleep 0.1s
done


# cside MC
code="Analysis_main"
MCver="MC15rd"
VerName="KumoiRD_LS_MC_cside"
Types=("CHG_analysis" "MIX_analysis" "UUBAR_analysis" "DDBAR_analysis" "SSBAR_analysis" "CHARM_analysis")

DirName="validation_v000"
SimpleDirName="validation"
for Type in "${Types[@]}"; do # for validation sample
    mkdir -p "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}"
    bsub -q s -o "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/${Type}_${DirName}_${code}.log" -e "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/${Type}_${DirName}_${code}_err.txt" "./submit_Analysis.sh" "./bin/${code}" ${DirName} ${SimpleDirName} ${Type} ${VerName} ${MCver}
    sleep 0.1s
done

