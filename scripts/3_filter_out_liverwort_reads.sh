#!/usr/bin/env bash

# Map the reads on C.cal genome to keep only the reads that don't map to the genome

dir="data/renamed_reads/"
sample_list="data/samples_batch2.txt"
out_dir="output/micro_fastq_batch2"
ref="data/ref/GCA_036924285.1_cmAstCali1_genomic.fna"
threads=8

mkdir -p "$out_dir"

while read -r sample; do
    # safe if coming from excel
    sample=${sample%$'\r'}  

    # file paths
    this_R1="${dir}${sample}_R1.fastq.gz"
    this_R2="${dir}${sample}_R2.fastq.gz"

    micro_R1="${out_dir}/micro_${sample}_R1.fastq.gz"
    micro_R2="${out_dir}/micro_${sample}_R2.fastq.gz"
    micro_single="${out_dir}/micro_${sample}_single.fastq.gz"

    echo "Working on sample ${sample}"
    echo "... mapping + extracting unmapped pairs -> FASTQ.gz"

    # Keep pairs where BOTH reads are unmapped:
    # -f 12 so we keep only reads whose mates are also unmapped
    # -F 256 excludes secondary alignments (good hygiene)
    # -F 2048 excludes supplementary alignments (good hygiene)
    bwa mem -t "$threads" "$ref" "$this_R1" "$this_R2" | samtools view -b -f 12 -F 256 -F 2048 - | samtools fastq -@ "$threads" \
        -1 >(gzip -c > "$micro_R1") \
        -2 >(gzip -c > "$micro_R2") \
        -s >(gzip -c > "$micro_single") \
        -n -

done < "${sample_list}"

