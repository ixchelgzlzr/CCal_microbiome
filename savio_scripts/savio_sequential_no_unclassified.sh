#!/usr/bin/env bash
#SBATCH --job-name=metaphlan_seq_skipunclass
#SBATCH --account=fc_phylodiv
#SBATCH --partition=savio3
#SBATCH --cpus-per-task=20
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mail-type=START,END,FAIL
#SBATCH --mail-user=ixchel_gonzalezrmz@berkeley.edu
set -euo pipefail

cd /global/scratch/users/ixchelgonzalezrmz/microbiome
echo "changed dir: $PWD"

# --- env ---
module load anaconda3
echo "loaded anaconda: $(conda --version 2>/dev/null || true)"

source activate metaphlan || { echo "Failed to activate 'metaphlan'"; exit 1; }

module load bio/bowtie2/2.5.1-gcc-11.4.0

echo "conda env:"
conda env list | grep metaphlan || true
echo "python: $(which python)"

echo "metaphlan version:"
metaphlan --version || true

echo "bowtie2 version:"
bowtie2 --version 2>/dev/null | head -n 1 || echo "bowtie2 not found on PATH"


# --- paths ---
PROJ="/global/scratch/users/ixchelgonzalezrmz/microbiome"
SAMPLE_LIST="/global/scratch/users/ixchelgonzalezrmz/microbiome/data/samples_all.txt"
DB_DIR="/global/scratch/users/ixchelgonzalezrmz/microbiome/database/mpa_vJan25_CHOCOPhlAnSGB_202503"
DB_IDX="mpa_vJan25_CHOCOPhlAnSGB_202503"
IN_DIR="output/micro_fastq"
OUT_DIR="output/metaphlan_out_no_unclassified"


# --- loop ---
while read -r sample; do
  [[ -z "${sample}" ]] && continue

  R1="output/micro_fastq/micro_${sample}_R1.fastq.gz"
  R2="output/micro_fastq/micro_${sample}_R2.fastq.gz"
  PROFILE="output/metaphlan_out_no_unclassified/${sample}_metaphlan.txt"
  MAPOUT="output/metaphlan_out_no_unclassified/${sample}_map.bowtie2.bz2"
  LOG="output/metaphlan_out_no_unclassified/${sample}_metaphlan.log"

  echo "============================================================"
  echo "$(date)  Starting ${sample}"
  echo "R1: $R1"
  echo "R2: $R2"

  if [[ ! -f "$R1" || ! -f "$R2" ]]; then
    echo "$(date)  Missing FASTQ(s) for ${sample}, skipping" | tee -a "$LOG"
    continue
  fi

  metaphlan "${R1},${R2}" --input_type fastq --db_dir "$DB_DIR" -x mpa_vJan25_CHOCOPhlAnSGB_202503 --offline --verbose --skip_unclassified_estimation --mapout "$MAPOUT" --nproc 20 > "$PROFILE" 2> "$LOG"

  echo "$(date)  Finished ${sample}"
done < "$SAMPLE_LIST"

echo "ALL DONE: $(date)"
