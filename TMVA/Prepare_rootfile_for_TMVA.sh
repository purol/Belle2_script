#!/bin/bash
# ./Prepare_rootfile_fot_TMVA.sh {ver name} {test dir name} {train dir name} {validation dir name}

if [ $# -ne 4 ]; then
    echo "./Prepare_rootfile_fot_TMVA.sh {ver name} {train dir name} {test dir name} {validation dir name}"
    exit 1
fi

CHGNAME="CHG_analysis"
MIXNAME="MIX_analysis"
UUBARNAME="UUBAR_analysis"
DDBARNAME="DDBAR_analysis"
SSBARNAME="SSBAR_analysis"
CHARMNAME="CHARM_analysis"
SIGNALNAME="SIGNAL_analysis"

BKGPATH="/home/jwpark/storage/BKG_gbasf2"

mkdir root_file_for_TMVA
cd ./root_file_for_TMVA

scp jwpark@login.yhep:/home/jwpark/AUTO/ToDataFile.c ./

scp jwpark@login.yhep:$BKGPATH/$1/$CHGNAME/$2/final_output/Mxs_small/*.root ./
scp jwpark@login.yhep:$BKGPATH/$1/$CHGNAME/$2/final_output/Mxs_large/*.root ./
scp jwpark@login.yhep:$BKGPATH/$1/$MIXNAME/$2/final_output/Mxs_small/*.root ./
scp jwpark@login.yhep:$BKGPATH/$1/$MIXNAME/$2/final_output/Mxs_large/*.root ./
scp jwpark@login.yhep:$BKGPATH/$1/$UUBARNAME/$2/final_output/Mxs_small/*.root ./
scp jwpark@login.yhep:$BKGPATH/$1/$UUBARNAME/$2/final_output/Mxs_large/*.root ./
scp jwpark@login.yhep:$BKGPATH/$1/$DDBARNAME/$2/final_output/Mxs_small/*.root ./
scp jwpark@login.yhep:$BKGPATH/$1/$DDBARNAME/$2/final_output/Mxs_large/*.root ./
scp jwpark@login.yhep:$BKGPATH/$1/$SSBARNAME/$2/final_output/Mxs_small/*.root ./
scp jwpark@login.yhep:$BKGPATH/$1/$SSBARNAME/$2/final_output/Mxs_large/*.root ./
scp jwpark@login.yhep:$BKGPATH/$1/$CHARMNAME/$2/final_output/Mxs_small/*.root ./
scp jwpark@login.yhep:$BKGPATH/$1/$CHARMNAME/$2/final_output/Mxs_large/*.root ./
scp jwpark@login.yhep:$BKGPATH/$1/$SIGNALNAME/$2/final_output/Mxs_small/*.root ./
scp jwpark@login.yhep:$BKGPATH/$1/$SIGNALNAME/$2/final_output/Mxs_large/*.root ./

scp jwpark@login.yhep:$BKGPATH/$1/$CHGNAME/$3/final_output/Mxs_small/*.root ./
scp jwpark@login.yhep:$BKGPATH/$1/$CHGNAME/$3/final_output/Mxs_large/*.root ./
scp jwpark@login.yhep:$BKGPATH/$1/$MIXNAME/$3/final_output/Mxs_small/*.root ./
scp jwpark@login.yhep:$BKGPATH/$1/$MIXNAME/$3/final_output/Mxs_large/*.root ./
scp jwpark@login.yhep:$BKGPATH/$1/$UUBARNAME/$3/final_output/Mxs_small/*.root ./
scp jwpark@login.yhep:$BKGPATH/$1/$UUBARNAME/$3/final_output/Mxs_large/*.root ./
scp jwpark@login.yhep:$BKGPATH/$1/$DDBARNAME/$3/final_output/Mxs_small/*.root ./
scp jwpark@login.yhep:$BKGPATH/$1/$DDBARNAME/$3/final_output/Mxs_large/*.root ./
scp jwpark@login.yhep:$BKGPATH/$1/$SSBARNAME/$3/final_output/Mxs_small/*.root ./
scp jwpark@login.yhep:$BKGPATH/$1/$SSBARNAME/$3/final_output/Mxs_large/*.root ./
scp jwpark@login.yhep:$BKGPATH/$1/$CHARMNAME/$3/final_output/Mxs_small/*.root ./
scp jwpark@login.yhep:$BKGPATH/$1/$CHARMNAME/$3/final_output/Mxs_large/*.root ./
scp jwpark@login.yhep:$BKGPATH/$1/$SIGNALNAME/$3/final_output/Mxs_small/*.root ./
scp jwpark@login.yhep:$BKGPATH/$1/$SIGNALNAME/$3/final_output/Mxs_large/*.root ./

scp jwpark@login.yhep:$BKGPATH/$1/$CHGNAME/$4/final_output/Mxs_small/*.root ./
scp jwpark@login.yhep:$BKGPATH/$1/$CHGNAME/$4/final_output/Mxs_large/*.root ./
scp jwpark@login.yhep:$BKGPATH/$1/$MIXNAME/$4/final_output/Mxs_small/*.root ./
scp jwpark@login.yhep:$BKGPATH/$1/$MIXNAME/$4/final_output/Mxs_large/*.root ./
scp jwpark@login.yhep:$BKGPATH/$1/$UUBARNAME/$4/final_output/Mxs_small/*.root ./
scp jwpark@login.yhep:$BKGPATH/$1/$UUBARNAME/$4/final_output/Mxs_large/*.root ./
scp jwpark@login.yhep:$BKGPATH/$1/$DDBARNAME/$4/final_output/Mxs_small/*.root ./
scp jwpark@login.yhep:$BKGPATH/$1/$DDBARNAME/$4/final_output/Mxs_large/*.root ./
scp jwpark@login.yhep:$BKGPATH/$1/$SSBARNAME/$4/final_output/Mxs_small/*.root ./
scp jwpark@login.yhep:$BKGPATH/$1/$SSBARNAME/$4/final_output/Mxs_large/*.root ./
scp jwpark@login.yhep:$BKGPATH/$1/$CHARMNAME/$4/final_output/Mxs_small/*.root ./
scp jwpark@login.yhep:$BKGPATH/$1/$CHARMNAME/$4/final_output/Mxs_large/*.root ./
scp jwpark@login.yhep:$BKGPATH/$1/$SIGNALNAME/$4/final_output/Mxs_small/*.root ./
scp jwpark@login.yhep:$BKGPATH/$1/$SIGNALNAME/$4/final_output/Mxs_large/*.root ./

.q | root .x ./ToDataFile.c

cd ../

mkdir -p ./For_TMVA_BKG/test/Mxs_large
mkdir -p ./For_TMVA_BKG/test/Mxs_small
mkdir -p ./For_TMVA_BKG/train/Mxs_large
mkdir -p ./For_TMVA_BKG/train/Mxs_small
mkdir -p ./For_TMVA_BKG/validation/Mxs_large
mkdir -p ./For_TMVA_BKG/validation/Mxs_small
mkdir -p ./For_TMVA_SIGNAL/test/Mxs_large
mkdir -p ./For_TMVA_SIGNAL/test/Mxs_small
mkdir -p ./For_TMVA_SIGNAL/train/Mxs_large
mkdir -p ./For_TMVA_SIGNAL/train/Mxs_small
mkdir -p ./For_TMVA_SIGNAL/validation/Mxs_large
mkdir -p ./For_TMVA_SIGNAL/validation/Mxs_small

mv ./root_file_for_TMVA/*larger*SIGNAL*$3_data.root ./For_TMVA_SIGNAL/test/Mxs_large
mv ./root_file_for_TMVA/*smaller*SIGNAL*$3_data.root ./For_TMVA_SIGNAL/test/Mxs_small
mv ./root_file_for_TMVA/*larger*SIGNAL*$2_data.root ./For_TMVA_SIGNAL/train/Mxs_large
mv ./root_file_for_TMVA/*smaller*SIGNAL*$2_data.root ./For_TMVA_SIGNAL/train/Mxs_small
mv ./root_file_for_TMVA/*larger*SIGNAL*$4_data.root ./For_TMVA_SIGNAL/validation/Mxs_large
mv ./root_file_for_TMVA/*smaller*SIGNAL*$4_data.root ./For_TMVA_SIGNAL/validation/Mxs_small

mv ./root_file_for_TMVA/*larger*$3_data.root ./For_TMVA_BKG/test/Mxs_large
mv ./root_file_for_TMVA/*smaller*$3_data.root ./For_TMVA_BKG/test/Mxs_small
mv ./root_file_for_TMVA/*larger*$2_data.root ./For_TMVA_BKG/train/Mxs_large
mv ./root_file_for_TMVA/*smaller*$2_data.root ./For_TMVA_BKG/train/Mxs_small
mv ./root_file_for_TMVA/*larger*$4_data.root ./For_TMVA_BKG/validation/Mxs_large
mv ./root_file_for_TMVA/*smaller*$4_data.root ./For_TMVA_BKG/validation/Mxs_small
