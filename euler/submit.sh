#!/bin/bash
# Login node, after 00_setup.sh and 01_download.sh: build the three references/indices, then process all 27 runs.
#   light array: FastQC + fastp + kallisto (4 GB per core; measured: 4-13 min per run, peak 3-23.7 GB incl. conversion and trimming; 2 GB per core is too little)
#                starts when the kallisto index is ready
#   STAR array:  fastp + STAR (~30 GB RAM each),           starts when the STAR index is ready
#   RSEM array:  RSEM on STAR's alignments (3 GB per core, 13-97 min per run), task i starts when STAR task i is done.
#                Measured peak memory over 13 runs: 6.0-14.8 GB (largest run 28 M pairs: 14.8 GB), so 24 GB leaves room.
#                RSEM is a separate array because it needs little memory: inside the STAR task it would block 48 GB for over an hour.
# RANGE selects the tasks of the STAR and RSEM arrays (default all 27), e.g. RANGE=7-27 bash euler/submit.sh to resume.
set -euo pipefail
export REPO="${REPO:-$HOME/RNAseq_aging}"
RANGE="${RANGE:-1-27}"
source "$REPO/euler/common.sh"
cd "$WORK/logs"                                  # job logs land here, not in the repo
kidx=$(sbatch --parsable "$REPO/euler/02b_kallisto_index.sh")
sidx=$(sbatch --parsable "$REPO/euler/02_star_index.sh")
ridx=$(sbatch --parsable "$REPO/euler/02c_rsem_reference.sh")
light=$(sbatch --parsable --dependency=afterok:"$kidx" --job-name=qc_kallisto --mem-per-cpu=4G --time=02:00:00 \
        --export=ALL,STAGES="qc+kallisto" "$REPO/euler/03_process_array.sh")
star=$(sbatch --parsable --dependency=afterok:"$sidx" --array="$RANGE%6" --job-name=star \
       --export=ALL,STAGES="star" "$REPO/euler/03_process_array.sh")
rsem=$(sbatch --parsable --dependency=aftercorr:"$star",afterok:"$ridx" --array="$RANGE%6" --job-name=rsem --mem-per-cpu=3G \
       --export=ALL,STAGES="rsem" "$REPO/euler/03_process_array.sh")
echo "indices $kidx $sidx $ridx -> light array $light | STAR array $star -> RSEM array $rsem   (watch: squeue -u \$USER)"
