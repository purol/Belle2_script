#!/bin/bash

mkdir -p "./output_particles"
mkdir -p "./log_particles"
mkdir -p "./err_particles"

code="./gbasf2_particles.py"
for file in $(find "/home/belle2/junewoo/storage_ghi/20231118_official_SIGNAL_MC15rd_official_feihadronic/output" -maxdepth 1 -name "*.root")
do
  echo $file
  basename=$(basename -s .root $file)

  bsub -q "s" \
  -o "./log_particles/${basename}.log" \
  -e "./err_particles/${basename}.err" \
  ${code} \
  --inputfile "${file}" \
  --destination "./output_particles" \
    sleep 0.5
done