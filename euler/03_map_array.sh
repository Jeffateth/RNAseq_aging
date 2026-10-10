#!/bin/bash
#SBATCH --job-name=star_map
#SBATCH --array=1-27%6
#SBATCH --ntasks=8
#SBATCH --mem-per-cpu=6G
#SBATCH --time=04:00:00
#SBATCH --output=%x_%A_%a.log
# One job per run: .sra -> FASTQ -> fastp (same settings as the kallisto pipeline) -> STAR -> delete FASTQ.
# --quantMode GeneCounts gives per-gene counts (-> DESeq2); TranscriptomeSAM keeps the input for RSEM (Tutorial 2-2-1).
set -euo pipefail
source "${REPO:-$HOME/RNAseq_aging}/euler/common.sh"
srr=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$REPO/metadata/srr_list.txt")
out="$MAP/$srr"
[[ -s "$out/Log.final.out" ]] && { echo "$srr done"; exit 0; }
mkdir -p "$out" "$TMP/$srr"; fq="$TMP/$srr"

fasterq-dump "$SRA/$srr/$srr.sra" --split-3 --threads 8 --temp "$fq" --outdir "$fq"
fastp -i "$fq/${srr}_1.fastq" -I "$fq/${srr}_2.fastq" -o "$fq/t_1.fastq.gz" -O "$fq/t_2.fastq.gz" \
      --detect_adapter_for_pe --length_required 20 --thread 8 --json "$QC/$srr.fastp.json" --html "$QC/$srr.fastp.html"
rm -f "$fq/${srr}_1.fastq" "$fq/${srr}_2.fastq"

STAR --genomeDir "$INDEX" --runThreadN 8 --readFilesIn "$fq/t_1.fastq.gz" "$fq/t_2.fastq.gz" --readFilesCommand zcat \
     --quantMode GeneCounts TranscriptomeSAM --outSAMtype BAM SortedByCoordinate --outFileNamePrefix "$out/"
rm -rf "$fq"
