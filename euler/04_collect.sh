#!/bin/bash
# Pack the small results (gene counts + mapping logs + fastp reports, a few MB) for copying to the laptop.
set -euo pipefail
source "${REPO:-$HOME/RNAseq_aging}/euler/common.sh"
cd "$WORK"
ls mapping/*/ReadsPerGene.out.tab | wc -l | xargs echo "runs finished:"
tar czf star_results.tar.gz mapping/*/ReadsPerGene.out.tab mapping/*/Log.final.out qc/*.json
ls -lh star_results.tar.gz
