#!/bin/bash
#SBATCH --job-name=star_index
#SBATCH --ntasks=10
#SBATCH --mem-per-cpu=6G
#SBATCH --time=04:00:00
#SBATCH --output=%x_%j.log
# Mouse genome index needs ~30 GB RAM (the reason this runs on Euler). sjdbOverhang = read length - 1 (reads are 150 bp).
set -euo pipefail
source "${REPO:-$HOME/RNAseq_aging}/euler/common.sh"
[[ -s "$INDEX/SAindex" ]] && { echo "index exists"; exit 0; }
STAR --runMode genomeGenerate --runThreadN 10 --genomeDir "$INDEX" \
     --genomeFastaFiles "$GENOME" --sjdbGTFfile "$GTF" --sjdbOverhang 149
