#!/bin/bsub
# ./submiter.sh ./submit_ReadRootFile.sh {root file path} {output destination} 

submit() {
  local VerName=$1
  local Type=$2
  local DirName=$3
  local SimpleDirName=$4
  local code=$5

  mkdir -p "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/final_output_root_after_MVA_Application_after_cut"
  bsub -q l -o "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/final_output_root_after_MVA_Application_after_cut/${Type}_${DirName}_${code}.log" -e "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/final_output_root_after_MVA_Application_after_cut/${Type}_${DirName}_${code}_err.txt" "./submit_MVACUT.sh" "./bin/${code}" ${DirName} ${SimpleDirName} ${Type} ${VerName} ${MCver}
}

# main analysis
MCver="MC15rd"
VerName="KumoiRD"
Types=("SIGNAL_analysis" "CHG_analysis" "MIX_analysis" "UUBAR_analysis" "DDBAR_analysis" "SSBAR_analysis" "CHARM_analysis")
code="MVACUT"

DirName="validation_v007"
SimpleDirName="validation"
for Type in "${Types[@]}"; do # for validation sample
    submit ${VerName} ${Type} ${DirName} ${SimpleDirName} ${code}
    sleep 0.1s
done


# off-resonance data
MCver="data"
VerName="KumoiRD_LS_data_off"
Types=("SIGNAL_analysis")
code="MVACUT"

DirName="validation_v007"
SimpleDirName="validation"
for Type in "${Types[@]}"; do
    submit ${VerName} ${Type} ${DirName} ${SimpleDirName} ${code}
    sleep 0.1s
done


# off-resonance MC
MCver="MC15rd"
VerName="KumoiRD_LS_MC_off"
Types=("UUBAR_analysis" "DDBAR_analysis" "SSBAR_analysis" "CHARM_analysis")
code="MVACUT"

DirName="validation_v007"
SimpleDirName="validation"
for Type in "${Types[@]}"; do 
    submit ${VerName} ${Type} ${DirName} ${SimpleDirName} ${code}
    sleep 0.1s
done


# sideband data
MCver="data"
VerName="KumoiRD_LS_data_side"
Types=("SIGNAL_analysis") 
code="MVACUT"

DirName="validation_v007"
SimpleDirName="validation"
for Type in "${Types[@]}"; do 
    submit ${VerName} ${Type} ${DirName} ${SimpleDirName} ${code}
    sleep 0.1s
done


## sideband MC
MCver="MC15rd"
VerName="KumoiRD_LS_MC_side"
Types=("CHG_analysis" "MIX_analysis" "UUBAR_analysis" "DDBAR_analysis" "SSBAR_analysis" "CHARM_analysis")
code="MVACUT"

DirName="validation_v007"
SimpleDirName="validation"
for Type in "${Types[@]}"; do
    submit ${VerName} ${Type} ${DirName} ${SimpleDirName} ${code}
    sleep 0.1s
done


# off-resonance sideband data
MCver="data"
VerName="KumoiRD_LS_data_offside"
Types=("SIGNAL_analysis")
code="MVACUT"

DirName="validation_v007"
SimpleDirName="validation"
for Type in "${Types[@]}"; do
    submit ${VerName} ${Type} ${DirName} ${SimpleDirName} ${code}
    sleep 0.1s
done


# off-resonance sideband MC
MCver="MC15rd"
VerName="KumoiRD_LS_MC_offside"
Types=("UUBAR_analysis" "DDBAR_analysis" "SSBAR_analysis" "CHARM_analysis")
code="MVACUT"

DirName="validation_v007"
SimpleDirName="validation"
for Type in "${Types[@]}"; do
    submit ${VerName} ${Type} ${DirName} ${SimpleDirName} ${code}
    sleep 0.1s
done


# EL sideband data
MCver="data"
VerName="KumoiRD_LS_data_EL"
Types=("SIGNAL_analysis")
code="MVACUT"

DirName="validation_v007"
SimpleDirName="validation"
for Type in "${Types[@]}"; do
    submit ${VerName} ${Type} ${DirName} ${SimpleDirName} ${code}
    sleep 0.1s
done


## EL sideband MC
MCver="MC15rd"
VerName="KumoiRD_LS_MC_EL"
Types=("CHG_analysis" "MIX_analysis" "UUBAR_analysis" "DDBAR_analysis" "SSBAR_analysis" "CHARM_analysis")
code="MVACUT"

DirName="validation_v007"
SimpleDirName="validation"
for Type in "${Types[@]}"; do
    submit ${VerName} ${Type} ${DirName} ${SimpleDirName} ${code}
    sleep 0.1s
done


# MU sideband data
MCver="data"
VerName="KumoiRD_LS_data_MU"
Types=("SIGNAL_analysis")
code="MVACUT"

DirName="validation_v007"
SimpleDirName="validation"
for Type in "${Types[@]}"; do
    submit ${VerName} ${Type} ${DirName} ${SimpleDirName} ${code}
    sleep 0.1s
done


## MU sideband MC
MCver="MC15rd"
VerName="KumoiRD_LS_MC_MU"
Types=("CHG_analysis" "MIX_analysis" "UUBAR_analysis" "DDBAR_analysis" "SSBAR_analysis" "CHARM_analysis")
code="MVACUT"

DirName="validation_v007"
SimpleDirName="validation"
for Type in "${Types[@]}"; do
    submit ${VerName} ${Type} ${DirName} ${SimpleDirName} ${code}
    sleep 0.1s
done


# Jpsi data
MCver="data"
VerName="KumoiRD_LS_data_Jpsi"
Types=("SIGNAL_analysis")
code="MVACUT"

DirName="validation_v007"
SimpleDirName="validation"
for Type in "${Types[@]}"; do
    submit ${VerName} ${Type} ${DirName} ${SimpleDirName} ${code}
    sleep 0.1s
done


# Jpsi MC
MCver="MC15rd"
VerName="KumoiRD_LS_MC_Jpsi"
Types=("SIGNAL_analysis" "CHG_analysis" "MIX_analysis" "UUBAR_analysis" "DDBAR_analysis" "SSBAR_analysis" "CHARM_analysis")
code="MVACUT"

DirName="validation_v007"
SimpleDirName="validation"
for Type in "${Types[@]}"; do
    submit ${VerName} ${Type} ${DirName} ${SimpleDirName} ${code}
    sleep 0.1s
done


# cside data
MCver="data"
VerName="KumoiRD_LS_data_cside"
Types=("SIGNAL_analysis")
code="MVACUT"

DirName="validation_v007"
SimpleDirName="validation"
for Type in "${Types[@]}"; do
    submit ${VerName} ${Type} ${DirName} ${SimpleDirName} ${code}
    sleep 0.1s
done


# cside MC
MCver="MC15rd"
VerName="KumoiRD_LS_MC_cside"
Types=("CHG_analysis" "MIX_analysis" "UUBAR_analysis" "DDBAR_analysis" "SSBAR_analysis" "CHARM_analysis")
code="MVACUT"

DirName="validation_v007"
SimpleDirName="validation"
for Type in "${Types[@]}"; do # for validation sample
    submit ${VerName} ${Type} ${DirName} ${SimpleDirName} ${code}
    sleep 0.1s
done


# FBDT side MC
MCver="MC15rd"
VerName="KumoiRD_LS_MC_FBDTside"
Types=("SIGNAL_analysis" "CHG_analysis" "MIX_analysis" "UUBAR_analysis" "DDBAR_analysis" "SSBAR_analysis" "CHARM_analysis")
code="MVACUT_side"

DirName="validation_v007"
SimpleDirName="validation"
for Type in "${Types[@]}"; do # for validation sample
    submit ${VerName} ${Type} ${DirName} ${SimpleDirName} ${code}
    sleep 0.1s
done


# FBDT side data
MCver="data"
VerName="KumoiRD_LS_data_FBDTside"
Types=("SIGNAL_analysis")
code="MVACUT_side"

DirName="validation_v007"
SimpleDirName="validation"
for Type in "${Types[@]}"; do
    submit ${VerName} ${Type} ${DirName} ${SimpleDirName} ${code}
    sleep 0.1s
done
