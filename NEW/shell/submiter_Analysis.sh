#!/bin/bsub
# ./submiter.sh ./submit_ReadRootFile.sh {root file path} {output destination} 

# main analysis
code="Analysis_main"
MCver="MC15ri"
VerName="Kumoi"
Types=("SIGNAL_analysis" "CHG_analysis" "MIX_analysis" "UUBAR_analysis" "DDBAR_analysis" "SSBAR_analysis" "CHARM_analysis")

DirName="test_v000"
SimpleDirName="test"
for Type in "${Types[@]}"; do # for test sample
    mkdir -p "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}"
    bsub -o "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/${Type}_${DirName}_${code}.log" -e "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/${Type}_${DirName}_${code}_err.txt" "./submit_Analysis.sh" "./bin/${code}" ${DirName} ${SimpleDirName} ${Type} ${VerName} ${MCver}
    sleep 0.1s
done

DirName="train_v000"
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

DirName="validation_v000"
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

DirName="train_v000"
SimpleDirName="train"
for Type in "${Types[@]}"; do 
    mkdir -p "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}"
    bsub -o "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/${Type}_${DirName}_${code}.log" -e "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/${Type}_${DirName}_${code}_err.txt" "./submit_Analysis.sh" "./bin/${code}" ${DirName} ${SimpleDirName} ${Type} ${VerName} ${MCver}
    sleep 0.1s
done

DirName="test_v000"
SimpleDirName="test"
for Type in "${Types[@]}"; do 
    mkdir -p "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}"
    bsub -o "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/${Type}_${DirName}_${code}.log" -e "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/${Type}_${DirName}_${code}_err.txt" "./submit_Analysis.sh" "./bin/${code}" ${DirName} ${SimpleDirName} ${Type} ${VerName} ${MCver}
    sleep 0.1s
done


## sideband data
#code="Analysis_side_data"
#MCver="MC15ri"
#VerName="Kumoi_LS_data_side"
#Types=("SIGNAL_analysis") 
#
#DirName="validation_v000"
#SimpleDirName="validation"
#for Type in "${Types[@]}"; do 
#    mkdir -p "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}"
#    bsub -o "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/${Type}_${DirName}_${code}.log" -e "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/${Type}_${DirName}_${code}_err.txt" "./submit_Analysis.sh" "./bin/${code}" ${DirName} ${SimpleDirName} ${Type} ${VerName} ${MCver}
#    sleep 0.1s
#done


## sideband MC
#code="Analysis_side_MC"
#MCver="MC15ri"
#VerName="Kumoi_LS_MC_side"
#Types=("CHG_analysis" "MIX_analysis" "UUBAR_analysis" "DDBAR_analysis" "SSBAR_analysis" "CHARM_analysis")
#
#DirName="train_v000"
#SimpleDirName="train"
#for Type in "${Types[@]}"; do 
#    mkdir -p "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}"
#    bsub -o "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/${Type}_${DirName}_${code}.log" -e "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/${Type}_${DirName}_${code}_err.txt" "./submit_Analysis.sh" "./bin/${code}" ${DirName} ${SimpleDirName} ${Type} ${VerName} ${MCver}
#    sleep 0.1s
#done

#DirName="test_v000"
#SimpleDirName="test"
#for Type in "${Types[@]}"; do
#    mkdir -p "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}"
#    bsub -o "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/${Type}_${DirName}_${code}.log" -e "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/${Type}_${DirName}_${code}_err.txt" "./submit_Analysis.sh" "./bin/${code}" ${DirName} ${SimpleDirName} ${Type} ${VerName} ${MCver}
#    sleep 0.1s
#done


## Jpsi data
#code="Analysis_main"
#MCver="MC15ri"
#VerName="Kumoi_LS_data_Jpsi"
#Types=("SIGNAL_analysis")
#
#DirName="validation_v000"
#SimpleDirName="validation"
#for Type in "${Types[@]}"; do
#    mkdir -p "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}"
#    bsub -o "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/${Type}_${DirName}_${code}.log" -e "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/${Type}_${DirName}_${code}_err.txt" "./submit_Analysis.sh" "./bin/${code}" ${DirName} ${SimpleDirName} ${Type} ${VerName} ${MCver}
#    sleep 0.1s
#done


## Jpsi MC
#code="Analysis_main"
#MCver="MC15ri"
#VerName="Kumoi_LS_MC_Jpsi"
#Types=("SIGNAL_analysis" "CHG_analysis" "MIX_analysis" "UUBAR_analysis" "DDBAR_analysis" "SSBAR_analysis" "CHARM_analysis")
#
#DirName="validation_v000"
#SimpleDirName="validation"
#for Type in "${Types[@]}"; do
#    mkdir -p "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}"
#    bsub -o "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/${Type}_${DirName}_${code}.log" -e "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/${Type}_${DirName}_${code}_err.txt" "./submit_Analysis.sh" "./bin/${code}" ${DirName} ${SimpleDirName} ${Type} ${VerName} ${MCver}
#    sleep 0.1s
#done


# syst
code="Analysis_main"
MCver="MC15ri"
VerName="Kumoi_syst"
Types=("SIGNAL_analysis")

DirName="validation_v000"
SimpleDirName="validation"
for Type in "${Types[@]}"; do
    mkdir -p "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}"
    bsub -o "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/${Type}_${DirName}_${code}.log" -e "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/${Type}_${DirName}_${code}_err.txt" "./submit_Analysis.sh" "./bin/${code}" ${DirName} ${SimpleDirName} ${Type} ${VerName} ${MCver}
    sleep 0.1s
done


## embedded 
#code="Analysis_main"
#MCver="MC15ri"
#VerName="Kumoi_embedded"
#Types=("SIGNAL_analysis")
#
#DirName="validation_v000"
#SimpleDirName="validation"
#for Type in "${Types[@]}"; do
#    mkdir -p "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}"
#    bsub -o "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/${Type}_${DirName}_${code}.log" -e "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/${Type}_${DirName}_${code}_err.txt" "./submit_Analysis.sh" "./bin/${code}" ${DirName} ${SimpleDirName} ${Type} ${VerName} ${MCver}
#    sleep 0.1s
#done


# cside data
code="Analysis_main"
MCver="MC15ri"
VerName="Kumoi_LS_data_cside"
Types=("SIGNAL_analysis")

DirName="validation_v000"
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

DirName="test_v000"
SimpleDirName="test"
for Type in "${Types[@]}"; do # for test sample
    mkdir -p "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}"
    bsub -o "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/${Type}_${DirName}_${code}.log" -e "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/${Type}_${DirName}_${code}_err.txt" "./submit_Analysis.sh" "./bin/${code}" ${DirName} ${SimpleDirName} ${Type} ${VerName} ${MCver}
    sleep 0.1s
done

DirName="train_v000"
SimpleDirName="train"
for Type in "${Types[@]}"; do # for train sample
    mkdir -p "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}"
    bsub -o "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/${Type}_${DirName}_${code}.log" -e "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/${Type}_${DirName}_${code}_err.txt" "./submit_Analysis.sh" "./bin/${code}" ${DirName} ${SimpleDirName} ${Type} ${VerName} ${MCver}
    sleep 0.1s
done

