#!/bin/bash
#SBATCH --job-name=trim_check
#SBATCH --ntasks=8
#SBATCH --mem-per-cpu=6G
#SBATCH --time=02:00:00
#SBATCH --output=%x_%j.log
# Sensitivity check asked for by the TA: does read trimming change STAR's uniquely mapped proportion?
# One run (SRR default below), two arms, identical STAR call and index:
#   raw      the untrimmed reads straight from fasterq-dump, no fastp at all
#   trimmed  the pipeline's fastp (adapters + poly-G + length >= 20 + default quality filter)
# The two arms are the two choices an automated pipeline can actually make; a hand-tuned middle ground
# would need manual inspection per run, which is what we want to avoid.
# Output: $WORK/trim_check_M39/<SRR>/<arm>/Log.final.out. Submit from $WORK/logs:
#   sbatch --export=ALL,SRR=SRR21355398 euler/06_trim_check.sh
set -euo pipefail
source "${REPO:-$HOME/RNAseq_aging}/euler/common.sh"
srr="${SRR:-SRR21355398}"          # smallest run (6.4 M pairs) and the highest adapter content (43 % of reads) -> worst case for trimming
out="$WORK/trim_check_M39/$srr"; fq="$TMP/trimcheck_$srr"
mkdir -p "$out" "$fq"

fasterq-dump "$SRA/$srr/$srr.sra" --split-3 --threads 8 --temp "$fq" --outdir "$fq"
fastp -i "$fq/${srr}_1.fastq" -I "$fq/${srr}_2.fastq" -o "$fq/tr_1.fastq.gz" -O "$fq/tr_2.fastq.gz" \
      --detect_adapter_for_pe --trim_poly_g --length_required 20 \
      --thread 8 --json "$out/trimmed.fastp.json" --html /dev/null   # identical to 03_process_array.sh

star_arm() {   # $1 arm name, $2/$3 reads, $4 read command
  mkdir -p "$out/$1"
  STAR --genomeDir "$INDEX" --runThreadN 8 --readFilesIn "$2" "$3" ${4:+--readFilesCommand $4} \
       --quantMode GeneCounts --outSAMtype None --outFileNamePrefix "$out/$1/"
}
star_arm raw     "$fq/${srr}_1.fastq" "$fq/${srr}_2.fastq" ""
star_arm trimmed "$fq/tr_1.fastq.gz"  "$fq/tr_2.fastq.gz"  zcat
rm -rf "$fq"
echo "done: $out"
