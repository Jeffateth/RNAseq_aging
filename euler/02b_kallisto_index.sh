#!/bin/bash
#SBATCH --job-name=kallisto_index
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=8G
#SBATCH --time=00:30:00
#SBATCH --output=%x_%j.log
# kallisto index of the GENCODE vM25 transcripts (a few GB RAM, minutes).
set -euo pipefail
source "${REPO:-$HOME/RNAseq_aging}/euler/common.sh"
[[ -s "$KAL_INDEX" ]] && { echo "index exists"; exit 0; }
kallisto index -i "$KAL_INDEX" "$TRANSCRIPTS"
