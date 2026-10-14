#!/bin/bash
# Login node, after 00_setup.sh and 01_download.sh: build the three references/indices, then process all 27 runs.
#   light array: FastQC + fastp + kallisto (small memory), starts when the kallisto index is ready
#   heavy array: STAR + RSEM (~30 GB RAM),                  starts when the STAR index and the RSEM reference are ready
set -euo pipefail
export REPO="${REPO:-$HOME/RNAseq_aging}"
source "$REPO/euler/common.sh"
cd "$WORK/logs"                                  # job logs land here, not in the repo
kidx=$(sbatch --parsable "$REPO/euler/02b_kallisto_index.sh")
sidx=$(sbatch --parsable "$REPO/euler/02_star_index.sh")
ridx=$(sbatch --parsable "$REPO/euler/02c_rsem_reference.sh")
light=$(sbatch --parsable --dependency=afterok:"$kidx" --job-name=qc_kallisto --mem-per-cpu=2G --time=02:00:00 \
        --export=ALL,STAGES="qc+kallisto" "$REPO/euler/03_process_array.sh")
heavy=$(sbatch --parsable --dependency=afterok:"$sidx":"$ridx" --job-name=star_rsem --export=ALL,STAGES="star+rsem" "$REPO/euler/03_process_array.sh")
echo "kallisto index $kidx -> light array $light | STAR index $sidx + RSEM reference $ridx -> heavy array $heavy   (watch: squeue -u \$USER)"
