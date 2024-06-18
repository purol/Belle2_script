#!/bin/bsub
# ./submiter.sh ./submit_ReadRootFile.sh {root file path} {output destination} 

submit() {
  code="Deconvertor"

  local VerName=$1
  local Type=$2
  local DirName=$3

  mkdir -p "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/final_output_root_after_MVA_Application"
  bsub -o "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/final_output_root_after_MVA_Application/${Type}_${DirName}_${code}.log" -e "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/final_output_root_after_MVA_Application/${Type}_${DirName}_${code}_err.txt" "./submit_Deconvertor.sh" "./bin/${code}" "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/final_output_data_after_MVA_Application" "-0.1" "-0.1" "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/final_output_root_after_MVA_Application"
}

# main analysis
VerName="Kumoi"
Types=("SIGNAL_analysis" "CHG_analysis" "MIX_analysis" "UUBAR_analysis" "DDBAR_analysis" "SSBAR_analysis" "CHARM_analysis")
code="Deconvertor"
DirName="test_v007"
for Type in "${Types[@]}"; do # for test sample
    submit ${VerName} ${Type} ${DirName}
    sleep 0.1s
done

DirName="train_v007"
for Type in "${Types[@]}"; do # for train sample
    submit ${VerName} ${Type} ${DirName}
    sleep 0.1s
done


# off-resonance data
VerName="Kumoi_LS_data_off"
Types=("SIGNAL_analysis")

DirName="validation_v007"
for Type in "${Types[@]}"; do # for test sample
    submit ${VerName} ${Type} ${DirName}
    sleep 0.1s
done

# off-resonance MC
VerName="Kumoi_LS_MC_off"
Types=("UUBAR_analysis" "DDBAR_analysis" "SSBAR_analysis" "CHARM_analysis")

DirName="train_v007"
for Type in "${Types[@]}"; do # for test sample
    submit ${VerName} ${Type} ${DirName}
    sleep 0.1s
done

DirName="test_v007"
for Type in "${Types[@]}"; do # for test sample
    submit ${VerName} ${Type} ${DirName}
    sleep 0.1s
done


# syst
VerName="Kumoi_syst"
Types=("SIGNAL_analysis")

DirName="validation_v007"
for Type in "${Types[@]}"; do # for test sample
    submit ${VerName} ${Type} ${DirName}
    sleep 0.1s
done


# cside data
VerName="Kumoi_LS_data_cside"
Types=("SIGNAL_analysis")

DirName="validation_v007"
for Type in "${Types[@]}"; do # for test sample
    submit ${VerName} ${Type} ${DirName}
    sleep 0.1s
done


# cside MC
VerName="Kumoi_LS_MC_cside"
Types=("CHG_analysis" "MIX_analysis" "UUBAR_analysis" "DDBAR_analysis" "SSBAR_analysis" "CHARM_analysis")

DirName="test_v007"
for Type in "${Types[@]}"; do # for test sample
    submit ${VerName} ${Type} ${DirName}
    sleep 0.1s
done

DirName="train_v007"
for Type in "${Types[@]}"; do # for test sample
    submit ${VerName} ${Type} ${DirName}
    sleep 0.1s
done
