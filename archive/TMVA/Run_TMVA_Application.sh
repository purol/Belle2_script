#!/bin/bash

OPTION=""

find ./For_TMVA_BKG -name "*larger*.root" -print0 | while read -d $'\0' file
do
    .q | root ./'TMVAClassificationApplication.C("'$file'")'
done

find ./For_TMVA_BKG -name "*smaller*.root" -print0 | while read -d $'\0' file
do
    .q | root ./'TMVAClassificationApplication.C("'$file'")'
done

find ./For_TMVA_SIGNAL -name "*larger*.root" -print0 | while read -d $'\0' file
do
    .q | root ./'TMVAClassificationApplication.C("'$file'")'
done

find ./For_TMVA_SIGNAL -name "*smaller*.root" -print0 | while read -d $'\0' file
do
    .q | root ./'TMVAClassificationApplication.C("'$file'")'
done

mkdir ./TMVA_output
mv ./*.root ./TMVA_output
