#!/bin/bash
# Pack the results for copying to the laptop (three tarballs in $WORK; paths inside are relative to $WORK):
#   star_results.tar.gz      STAR gene counts + mapping logs (~11 MB)
#   kallisto_results.tar.gz  kallisto estimates per run (~0.7 GB)
#   qc_results.tar.gz        FastQC zips (raw reads) + fastp reports (~tens of MB)
set -euo pipefail
source "${REPO:-$HOME/RNAseq_aging}/euler/common.sh"
cd "$WORK"
echo "finished: STAR $(ls mapping/*/Log.final.out 2>/dev/null | wc -l)/27 | kallisto $(ls kallisto/*/abundance.h5 2>/dev/null | wc -l)/27 | FastQC $(ls qc/raw/*_2_fastqc.zip 2>/dev/null | wc -l)/27"
tar czf star_results.tar.gz mapping/*/ReadsPerGene.out.tab mapping/*/Log.final.out
tar czf kallisto_results.tar.gz kallisto/*/abundance.tsv kallisto/*/abundance.h5 kallisto/*/run_info.json
tar czf qc_results.tar.gz qc/raw/*_fastqc.zip qc/*.fastp.json qc/*.fastp.html
ls -lh star_results.tar.gz kallisto_results.tar.gz qc_results.tar.gz
