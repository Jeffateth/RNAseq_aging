#!/bin/bash
# Login node, after 00_setup.sh and 01_download.sh: build the index, then map all 27 runs once it finishes.
set -euo pipefail
export REPO="${REPO:-$HOME/RNAseq_aging}"
source "$REPO/euler/common.sh"
cd "$WORK/logs"                                  # job logs land here, not in the repo
idx=$(sbatch --parsable "$REPO/euler/02_star_index.sh")
map=$(sbatch --parsable --dependency=afterok:"$idx" "$REPO/euler/03_map_array.sh")
echo "index job $idx -> mapping array $map   (watch with: squeue -u \$USER)"
