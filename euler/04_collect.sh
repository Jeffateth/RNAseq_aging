#!/bin/bash
# Pack the results for copying to the laptop (four tarballs in $WORK). Paths inside are relative and carry no _M39 suffix, so the
# laptop side stays reference-independent:
#   star_results.tar.gz      mapping/<SRR>/  STAR gene counts + mapping logs (~11 MB)
#   kallisto_results.tar.gz  kallisto/<SRR>/ kallisto estimates per run (~0.6 GB)
#   rsem_results.tar.gz      rsem/<SRR>/     RSEM gene-level results per run (~70 MB packed)
#   qc_results.tar.gz        qc/             FastQC zips (raw reads) + fastp reports (~tens of MB)
#   qc_trimmed_results.tar.gz qc/trimmed/    FastQC zips of the trimmed reads (separate, so the four tarballs above stay as they were)
set -euo pipefail
source "${REPO:-$HOME/RNAseq_aging}/euler/common.sh"
cd "$WORK"
echo "finished: STAR $(ls mapping_M39/*/Log.final.out 2>/dev/null | wc -l)/27 | kallisto $(ls kallisto_M39/*/abundance.h5 2>/dev/null | wc -l)/27 | RSEM $(ls rsem_M39/*/*.genes.results 2>/dev/null | wc -l)/27 | FastQC $(ls qc/raw/*_2_fastqc.zip 2>/dev/null | wc -l)/27 | FastQC trimmed $(ls qc/trimmed/*_trim_2_fastqc.zip 2>/dev/null | wc -l)/27"
tar czf star_results.tar.gz --transform 's,^mapping_M39/,mapping/,' mapping_M39/*/ReadsPerGene.out.tab mapping_M39/*/Log.final.out
tar czf kallisto_results.tar.gz --transform 's,^kallisto_M39/,kallisto/,' kallisto_M39/*/abundance.tsv kallisto_M39/*/abundance.h5 kallisto_M39/*/run_info.json
tar czf rsem_results.tar.gz --transform 's,^rsem_M39/,rsem/,' rsem_M39/*/*.genes.results
tar czf qc_results.tar.gz qc/raw/*_fastqc.zip qc/*.fastp.json qc/*.fastp.html
[[ -d qc/trimmed ]] && tar czf qc_trimmed_results.tar.gz qc/trimmed/*_fastqc.zip
ls -lh *_results.tar.gz
