#!/bin/bash
#SBATCH --job-name=kallisto_index
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=32G
#SBATCH --time=01:00:00
#SBATCH --output=%x_%j.log
# kallisto index of the GENCODE M39 transcripts (482k transcripts; 8 GB was too little: killed for memory after 14 min; 32 GB requested).
# The index is built under a temporary name and renamed only when kallisto succeeds: a killed job once left a truncated file under the final
# name, and the next run skipped the build because the file existed (every kallisto task then failed while reading it).
set -euo pipefail
source "${REPO:-$HOME/RNAseq_aging}/euler/common.sh"
[[ -s "$KAL_INDEX" ]] && { echo "index exists"; exit 0; }
rm -f "$KAL_INDEX.tmp"
kallisto index -i "$KAL_INDEX.tmp" "$TRANSCRIPTS"
mv "$KAL_INDEX.tmp" "$KAL_INDEX"
