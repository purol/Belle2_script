#!/bin/bsub
# ./submiter.sh ./submit_ReadRootFile.sh {root file path} {output destination} 

code="MVACUT"

submit() {
  code="MVACUT"

  local VerName=$1
  local Type=$2
  local DirName=$3
  local SimpleDirName=$4

  mkdir -p "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/final_output_root_after_MVA_Application_after_cut"
  bsub -q s -o "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/final_output_root_after_MVA_Application_after_cut/${Type}_${DirName}_${code}.log" -e "/home/belle2/junewoo/storage_b1/bsub/Analysis/${VerName}/${Type}/${DirName}/final_output_root_after_MVA_Application_after_cut/${Type}_${DirName}_${code}_err.txt" "./submit_MVACUT.sh" "./bin/${code}" ${DirName} ${SimpleDirName} ${Type} ${VerName} ${MCver}
}

# main analysis
MCver="MC15ri"
VerName="Kumoi"
Types=("SIGNAL_analysis" "CHG_analysis" "MIX_analysis" "UUBAR_analysis" "DDBAR_analysis" "SSBAR_analysis" "CHARM_analysis")

DirName="test_v007"
SimpleDirName="test"
for Type in "${Types[@]}"; do # for test sample
    submit ${VerName} ${Type} ${DirName} ${SimpleDirName}
    sleep 0.1s
done

DirName="train_v007"
SimpleDirName="train"
for Type in "${Types[@]}"; do # for train sample
    submit ${VerName} ${Type} ${DirName} ${SimpleDirName}
    sleep 0.1s
done


# off-resonance data
MCver="data"
VerName="Kumoi_LS_data_off"
Types=("SIGNAL_analysis")

DirName="validation_v007"
SimpleDirName="validation"
for Type in "${Types[@]}"; do
    submit ${VerName} ${Type} ${DirName} ${SimpleDirName}
    sleep 0.1s
done


# off-resonance MC
MCver="MC15ri"
VerName="Kumoi_LS_MC_off"
Types=("UUBAR_analysis" "DDBAR_analysis" "SSBAR_analysis" "CHARM_analysis")

DirName="train_v007"
SimpleDirName="train"
for Type in "${Types[@]}"; do 
    submit ${VerName} ${Type} ${DirName} ${SimpleDirName}
    sleep 0.1s
done

DirName="test_v007"
SimpleDirName="test"
for Type in "${Types[@]}"; do 
    submit ${VerName} ${Type} ${DirName} ${SimpleDirName}
    sleep 0.1s
done



# syst
MCver="MC15ri"
VerName="Kumoi_syst"
Types=("SIGNAL_analysis")

DirName="validation_v007"
SimpleDirName="validation"
for Type in "${Types[@]}"; do
    submit ${VerName} ${Type} ${DirName} ${SimpleDirName}
    sleep 0.1s
done


# cside data
MCver="MC15ri"
VerName="Kumoi_LS_data_cside"
Types=("SIGNAL_analysis")

DirName="validation_v007"
SimpleDirName="validation"
for Type in "${Types[@]}"; do
    submit ${VerName} ${Type} ${DirName} ${SimpleDirName}
    sleep 0.1s
done


# cside MC
MCver="MC15ri"
VerName="Kumoi_LS_MC_cside"
Types=("CHG_analysis" "MIX_analysis" "UUBAR_analysis" "DDBAR_analysis" "SSBAR_analysis" "CHARM_analysis")

DirName="test_v007"
SimpleDirName="test"
for Type in "${Types[@]}"; do # for test sample
    submit ${VerName} ${Type} ${DirName} ${SimpleDirName}
    sleep 0.1s
done

DirName="train_v007"
SimpleDirName="train"
for Type in "${Types[@]}"; do # for train sample
    submit ${VerName} ${Type} ${DirName} ${SimpleDirName}
    sleep 0.1s
done

