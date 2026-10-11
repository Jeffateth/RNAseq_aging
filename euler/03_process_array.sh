#!/bin/bash
#SBATCH --job-name=process
#SBATCH --array=1-27%6
#SBATCH --ntasks=8
#SBATCH --mem-per-cpu=6G
#SBATCH --time=04:00:00
#SBATCH --output=%x_%A_%a.log
# One task per run. STAGES (environment variable, stage names joined by +, default "qc+kallisto+star") selects the work:
#   qc        FastQC on the raw reads (-> $QC/raw) and the fastp report (-> $QC)
#   kallisto  kallisto on the fastp-trimmed reads (-> $KAL/<SRR>)
#   star      STAR with gene counting on the same trimmed reads (-> $MAP/<SRR>)
# Stages whose output already exists are skipped, so a rerun does only what is missing. submit.sh runs two arrays:
# a light one ("qc+kallisto", small memory) and a heavy one ("star", ~30 GB for the genome index).
# fastp settings are identical to the laptop pipeline (scripts/helpers.py: trim).
set -euo pipefail
source "${REPO:-$HOME/RNAseq_aging}/euler/common.sh"
STAGES="${STAGES:-qc+kallisto+star}"
srr=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$REPO/metadata/srr_list.txt")
has() { [[ "+$STAGES+" == *"+$1+"* ]]; }

todo_qc=0; todo_kal=0; todo_star=0
if has qc       && [[ ! -s "$QC/raw/${srr}_2_fastqc.zip" || ! -s "$QC/$srr.fastp.json" ]]; then todo_qc=1; fi
if has kallisto && [[ ! -s "$KAL/$srr/abundance.h5" ]];       then todo_kal=1; fi
if has star     && [[ ! -s "$MAP/$srr/Log.final.out" ]];      then todo_star=1; fi     # Log.final.out is written last by STAR
if (( todo_qc + todo_kal + todo_star == 0 )); then echo "$srr: nothing to do ($STAGES)"; exit 0; fi

fq="$TMP/$srr"; mkdir -p "$fq"
fasterq-dump "$SRA/$srr/$srr.sra" --split-3 --threads 8 --temp "$fq" --outdir "$fq"

if (( todo_qc )) && [[ ! -s "$QC/raw/${srr}_2_fastqc.zip" ]]; then
  fastqc -q -t 8 --outdir "$QC/raw" "$fq/${srr}_1.fastq" "$fq/${srr}_2.fastq"
fi
if (( todo_kal + todo_star )); then
  json=/dev/null; html=/dev/null                                        # keep the fastp report from the qc stage only
  if (( todo_qc )); then json="$QC/$srr.fastp.json"; html="$QC/$srr.fastp.html"; fi
  fastp -i "$fq/${srr}_1.fastq" -I "$fq/${srr}_2.fastq" -o "$fq/${srr}_trim_1.fastq.gz" -O "$fq/${srr}_trim_2.fastq.gz" \
        --detect_adapter_for_pe --trim_poly_g --length_required 20 --thread 8 --json "$json" --html "$html"
fi
rm -f "$fq/${srr}_1.fastq" "$fq/${srr}_2.fastq"

if (( todo_kal )); then
  kallisto quant -i "$KAL_INDEX" -o "$KAL/$srr" --threads 8 "$fq/${srr}_trim_1.fastq.gz" "$fq/${srr}_trim_2.fastq.gz"   # 8 threads, as on the laptop
fi
if (( todo_star )); then
  mkdir -p "$MAP/$srr"
  STAR --genomeDir "$INDEX" --runThreadN 8 --readFilesIn "$fq/${srr}_trim_1.fastq.gz" "$fq/${srr}_trim_2.fastq.gz" --readFilesCommand zcat \
       --quantMode GeneCounts TranscriptomeSAM --outSAMtype BAM SortedByCoordinate --outFileNamePrefix "$MAP/$srr/"
fi
rm -rf "$fq"
