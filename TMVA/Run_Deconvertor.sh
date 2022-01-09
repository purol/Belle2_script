#!/bin/bash
# Deconvertor(const char* filename, double OBB, double OContinuum)

find ./TMVA_output -name "*larger*.root" -print0 | while read -d $'\0' file
do
    .q | root ./'Deconvertor.c("'$file'",0.6,0.84)'
done

find ./TMVA_output -name "*smaller*.root" -print0 | while read -d $'\0' file
do
    .q | root ./'Deconvertor.c("'$file'",0.6,0.94)'
done

