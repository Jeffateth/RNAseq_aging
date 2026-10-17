#!/bin/bash
#SBATCH --job-name=process
#SBATCH --array=1-27%6
#SBATCH --ntasks=8
#SBATCH --mem-per-cpu=6G
#SBATCH --time=04:00:00
#SBATCH --output=%x_%A_%a.log
# One task per run. STAGES (environment variable, stage names joined by +, default "qc+kallisto+star+rsem") selects the work:
#   qc        FastQC on the raw reads (-> $QC/raw) and the fastp report (-> $QC)
#   kallisto  kallisto on the fastp-trimmed reads (-> $KAL/<SRR>)
#   star      STAR with gene counting on the same trimmed reads (-> $MAP/<SRR>), also writes the transcriptome alignments for RSEM
#   rsem      RSEM on STAR's transcriptome alignments (-> $RSEM/<SRR>); needs the star stage done first
# Stages whose output already exists are skipped, so a rerun does only what is missing. submit.sh runs two arrays:
# a light one ("qc+kallisto", 4 GB per core), a STAR one ("star", ~30 GB for the genome index) and an RSEM one ("rsem", 3 GB per core;
# RSEM takes 40-90 min per run, so it is kept out of the STAR tasks that reserve 48 GB).
# fastp settings are identical to the laptop pipeline (scripts/helpers.py: trim).
set -euo pipefail
source "${REPO:-$HOME/RNAseq_aging}/euler/common.sh"
STAGES="${STAGES:-qc+kallisto+star+rsem}"
srr=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$REPO/metadata/srr_list.txt")
has() { [[ "+$STAGES+" == *"+$1+"* ]]; }

todo_qc=0; todo_kal=0; todo_star=0; todo_rsem=0
if has qc       && [[ ! -s "$QC/raw/${srr}_2_fastqc.zip" || ! -s "$QC/$srr.fastp.json" ]]; then todo_qc=1; fi
if has kallisto && [[ ! -s "$KAL/$srr/abundance.h5" ]];       then todo_kal=1; fi
if has star     && [[ ! -s "$MAP/$srr/Log.final.out" ]];      then todo_star=1; fi     # Log.final.out is written last by STAR
if has rsem     && [[ ! -s "$RSEM/$srr/$srr.genes.results" ]]; then todo_rsem=1; fi
if (( todo_qc + todo_kal + todo_star + todo_rsem == 0 )); then echo "$srr: nothing to do ($STAGES)"; exit 0; fi

fq="$TMP/$srr"
if (( todo_qc + todo_kal + todo_star )); then          # RSEM alone needs no reads, only STAR's alignments
  mkdir -p "$fq"
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
    rm -rf "$MAP/$srr"; mkdir -p "$MAP/$srr"                              # a killed earlier attempt leaves _STARtmp behind, and STAR refuses to start then
    STAR --genomeDir "$INDEX" --runThreadN 8 --readFilesIn "$fq/${srr}_trim_1.fastq.gz" "$fq/${srr}_trim_2.fastq.gz" --readFilesCommand zcat \
         --quantMode GeneCounts TranscriptomeSAM --outSAMtype None --outFileNamePrefix "$MAP/$srr/"      # no genome BAM: not used, saves time and disk
  fi
  rm -rf "$fq"
fi

if (( todo_rsem )); then
  mkdir -p "$RSEM/$srr"
  rsem-calculate-expression --alignments --paired-end --no-bam-output -p 8 -q "$MAP/$srr/Aligned.toTranscriptome.out.bam" "$RSEM_REF" "$RSEM/$srr/$srr"   # unstranded (default); --no-bam-output: no extra BAM copy of the alignments
fi
