#!/usr/bin/env bash
#SBATCH --job-name=metaphlan_install
#SBATCH --account=fc_phylodiv
#SBATCH --partition=savio3
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1


# Load environment (adjust if needed)
module load anaconda3
module load bio/bowtie2/2.5.1-gcc-11.4.0
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate metaphlan 

# Run database install
metaphlan --install --db_dir database