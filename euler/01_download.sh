#!/bin/bash
# Login node (run inside tmux/screen): GENCODE M39 (GRCm39) genome, annotation and transcripts, and all 27 .sra files. Resumable.
set -euo pipefail
source "${REPO:-$HOME/RNAseq_aging}/euler/common.sh"
base=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M39
[[ -s "$GENOME" ]] || { wget -nv -O "$GENOME.gz" $base/GRCm39.primary_assembly.genome.fa.gz && gzip -d "$GENOME.gz"; }
[[ -s "$GTF" ]]    || { wget -nv -O "$GTF.gz" $base/gencode.vM39.annotation.gtf.gz && gzip -d "$GTF.gz"; }
[[ -s "$TRANSCRIPTS" ]] || wget -nv -O "$TRANSCRIPTS" $base/gencode.vM39.transcripts.fa.gz      # kallisto reference
for srr in $(cat "$REPO/metadata/srr_list.txt"); do
  [[ -s "$SRA/$srr/$srr.sra" ]] && continue
  prefetch "$srr" --output-directory "$SRA" --max-size 100G
done
echo "downloaded: $(ls -d "$SRA"/SRR* | wc -l)/27 runs"
