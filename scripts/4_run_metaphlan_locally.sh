#!/usr/bin/env bash
set -euo pipefail

# run metaphlan
sample_list="data/samples.txt"

mkdir -p output/metaphlan_out

# use a loop
while read -r sample; do 

  this_R1="output/micro_fastq/micro_${sample}_R1.fastq.gz"
  this_R2="output/micro_fastq/micro_${sample}_R2.fastq.gz"
  output="output/metaphlan_out/${sample}_metaphlan.txt"
  mapout="output/metaphlan_out/${sample}_map.bowtie2.bz2"

  echo "Working on ${sample}..."

  # check that the files exist  
  if [[ ! -f "$this_R1" || ! -f "$this_R2" ]]; then
    echo "Missing files for ${sample}, skipping"
    continue
  fi

  metaphlan "${this_R1},${this_R2}" --input_type fastq --mapout "${mapout}" --nproc 8 > "$output"

done < "$sample_list"
