#!/usr/bin/env bash
set -euo pipefail

############################################
# Quality check on reads with trimmommatic #
############################################

# set up paths
TRIMJAR="/Users/ixchelgonzalezramirez/software/Trimmomatic-0.39/trimmomatic-0.39.jar"
ADAPTERS="/Users/ixchelgonzalezramirez/software/Trimmomatic-0.39/adapters/TruSeq3-PE.fa"


# count how many pairs are there two process
total=$(ls data/raw_reads/*_R1_*.fastq.gz | wc -l)
count=0

# loop over files
for file in data/raw_reads/*_R1_*.fastq.gz; do

  # get the basename  
  fn="${file##*/}"                             
  base="${fn%_R1_001.fastq.gz}" 

  # update progress
  count=$((count + 1))
  echo "[$count / $total] Processing $base"

  # check that the two files exist
  fwd="data/raw_reads/${base}_R1_001.fastq.gz"      
  rev="data/raw_reads/${base}_R2_001.fastq.gz"

  if [[ -f "$fwd" && -f "$rev" ]]; then
    echo "Pair exists:"
    echo "  $fwd"
    echo "  $rev"

    # if the two files exist, then run trimmommatic
    # first get the file names
    fwd_out="data/clean_reads/${base}_R1_001.fastq.gz" 
    rev_out="data/clean_reads/${base}_R2_001.fastq.gz" 
    fwd_out_unpaired="data/clean_reads/unpaired_${base}_R1_001.fastq.gz"
    rev_out_unpaired="data/clean_reads/unpaired_${base}_R2_001.fastq.gz"

    # then run 
    java -jar "$TRIMJAR" PE -threads 8 -phred33 \
    "$fwd" "$rev" \
    "$fwd_out" "$fwd_out_unpaired" \
    "$rev_out" "$rev_out_unpaired" \
    "ILLUMINACLIP:${ADAPTERS}:2:30:10:2:True" \
    LEADING:3 TRAILING:3 MINLEN:36

  else
    echo "Missing pair for:"
    echo "  $fwd"
  fi
done

