#!/bin/bsub
# ./submiter.sh ./submit_ReadRootFile.sh {root file path} {output destination} 

submit() {
  code="Deconvertor"

  local VerName=$1
  local Type=$2
  local DirName=$3

  mkdir -p "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/final_output_root_after_MVA_Application"
  bsub -q l -o "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/final_output_root_after_MVA_Application/${Type}_${DirName}_${code}.log" -e "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/final_output_root_after_MVA_Application/${Type}_${DirName}_${code}_err.txt" "./submit_Deconvertor.sh" "./bin/${code}" "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/final_output_data_after_MVA_Application" "-0.1" "-0.1" "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/final_output_root_after_MVA_Application"
}

# main analysis
VerName="KumoiRD"
Types=("SIGNAL_analysis" "CHG_analysis" "MIX_analysis" "UUBAR_analysis" "DDBAR_analysis" "SSBAR_analysis" "CHARM_analysis")
code="Deconvertor"
DirName="validation_v007"
for Type in "${Types[@]}"; do # for validation sample
    submit ${VerName} ${Type} ${DirName}
    sleep 0.1s
done


# off-resonance data
VerName="KumoiRD_LS_data_off"
Types=("SIGNAL_analysis")

DirName="validation_v007"
for Type in "${Types[@]}"; do # for validation sample
    submit ${VerName} ${Type} ${DirName}
    sleep 0.1s
done

# off-resonance MC
VerName="KumoiRD_LS_MC_off"
Types=("UUBAR_analysis" "DDBAR_analysis" "SSBAR_analysis" "CHARM_analysis")

DirName="validation_v007"
for Type in "${Types[@]}"; do # for validation sample
    submit ${VerName} ${Type} ${DirName}
    sleep 0.1s
done


# sideband data
VerName="KumoiRD_LS_data_side"
Types=("SIGNAL_analysis") 

DirName="validation_v007"
for Type in "${Types[@]}"; do # for validation sample
    submit ${VerName} ${Type} ${DirName}
    sleep 0.1s
done


## sideband MC
VerName="KumoiRD_LS_MC_side"
Types=("CHG_analysis" "MIX_analysis" "UUBAR_analysis" "DDBAR_analysis" "SSBAR_analysis" "CHARM_analysis")

DirName="validation_v007"
for Type in "${Types[@]}"; do # for validation sample
    submit ${VerName} ${Type} ${DirName}
    sleep 0.1s
done


# off-resonance sideband data
VerName="KumoiRD_LS_data_offside"
Types=("SIGNAL_analysis")

DirName="validation_v007"
for Type in "${Types[@]}"; do # for validation sample
    submit ${VerName} ${Type} ${DirName}
    sleep 0.1s
done

# off-resonance sideband MC
VerName="KumoiRD_LS_MC_offside"
Types=("UUBAR_analysis" "DDBAR_analysis" "SSBAR_analysis" "CHARM_analysis")

DirName="validation_v007"
for Type in "${Types[@]}"; do # for validation sample
    submit ${VerName} ${Type} ${DirName}
    sleep 0.1s
done


# EL sideband data
VerName="KumoiRD_LS_data_EL"
Types=("SIGNAL_analysis")

DirName="validation_v007"
for Type in "${Types[@]}"; do # for validation sample
    submit ${VerName} ${Type} ${DirName}
    sleep 0.1s
done


## EL sideband MC
VerName="KumoiRD_LS_MC_EL"
Types=("CHG_analysis" "MIX_analysis" "UUBAR_analysis" "DDBAR_analysis" "SSBAR_analysis" "CHARM_analysis")

DirName="validation_v007"
for Type in "${Types[@]}"; do # for validation sample
    submit ${VerName} ${Type} ${DirName}
    sleep 0.1s
done


# MU sideband data
VerName="KumoiRD_LS_data_MU"
Types=("SIGNAL_analysis")

DirName="validation_v007"
for Type in "${Types[@]}"; do # for validation sample
    submit ${VerName} ${Type} ${DirName}
    sleep 0.1s
done


## MU sideband MC
VerName="KumoiRD_LS_MC_MU"
Types=("CHG_analysis" "MIX_analysis" "UUBAR_analysis" "DDBAR_analysis" "SSBAR_analysis" "CHARM_analysis")

DirName="validation_v007"
for Type in "${Types[@]}"; do # for validation sample
    submit ${VerName} ${Type} ${DirName}
    sleep 0.1s
done


# Jpsi data
VerName="KumoiRD_LS_data_Jpsi"
Types=("SIGNAL_analysis")

DirName="validation_v007"
for Type in "${Types[@]}"; do # for validation sample
    submit ${VerName} ${Type} ${DirName}
    sleep 0.1s
done


# Jpsi MC
VerName="KumoiRD_LS_MC_Jpsi"
Types=("SIGNAL_analysis" "CHG_analysis" "MIX_analysis" "UUBAR_analysis" "DDBAR_analysis" "SSBAR_analysis" "CHARM_analysis")

DirName="validation_v007"
for Type in "${Types[@]}"; do # for validation sample
    submit ${VerName} ${Type} ${DirName}
    sleep 0.1s
done


# cside data
VerName="KumoiRD_LS_data_cside"
Types=("SIGNAL_analysis")

DirName="validation_v007"
for Type in "${Types[@]}"; do # for validation sample
    submit ${VerName} ${Type} ${DirName}
    sleep 0.1s
done


# cside MC
VerName="KumoiRD_LS_MC_cside"
Types=("CHG_analysis" "MIX_analysis" "UUBAR_analysis" "DDBAR_analysis" "SSBAR_analysis" "CHARM_analysis")

DirName="validation_v007"
for Type in "${Types[@]}"; do # for validation sample
    submit ${VerName} ${Type} ${DirName}
    sleep 0.1s
done


# FBDT side MC
VerName="KumoiRD_LS_MC_FBDTside"
Types=("SIGNAL_analysis" "CHG_analysis" "MIX_analysis" "UUBAR_analysis" "DDBAR_analysis" "SSBAR_analysis" "CHARM_analysis")

DirName="validation_v007"
for Type in "${Types[@]}"; do # for validation sample
    submit ${VerName} ${Type} ${DirName}
    sleep 0.1s
done


# FBDT side data
VerName="KumoiRD_LS_data_FBDTside"
Types=("SIGNAL_analysis")

DirName="validation_v007"
for Type in "${Types[@]}"; do # for validation sample
    submit ${VerName} ${Type} ${DirName}
    sleep 0.1s
done
