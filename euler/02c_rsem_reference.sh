#!/bin/bash
#SBATCH --job-name=rsem_ref
#SBATCH --ntasks=4
#SBATCH --mem-per-cpu=4G
#SBATCH --time=02:00:00
#SBATCH --output=%x_%j.log
# RSEM reference (transcript sequences and gene/transcript map) from the same genome and annotation as the STAR index.
set -euo pipefail
source "${REPO:-$HOME/RNAseq_aging}/euler/common.sh"
[[ -s "$RSEM_REF.transcripts.fa" ]] && { echo "reference exists"; exit 0; }
rsem-prepare-reference --gtf "$GTF" -p 4 "$GENOME" "$RSEM_REF"
