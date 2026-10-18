#!/bin/bash
# Login node: check that every run has every result and that the numbers look sane, before 04_collect.sh. Read-only.
# Prints OK / PROBLEM lines; the last line says whether to go on.
set -uo pipefail
source "${REPO:-$HOME/RNAseq_aging}/euler/common.sh"
bad=0; n=0
problem() { echo "PROBLEM: $*"; bad=$((bad+1)); }

echo "== 1. jobs still queued or running"
q=$(squeue -u "$USER" -h 2>/dev/null | wc -l); [[ $q -eq 0 ]] && echo "OK: none" || problem "$q jobs in the queue (squeue -u \$USER)"

echo "== 2. every run has every result (missing runs listed)"
for spec in "STAR:$MAP/%s/Log.final.out" "kallisto:$KAL/%s/abundance.h5" "RSEM:$RSEM/%s/%s.genes.results" "FastQC:$QC/raw/%s_2_fastqc.zip" "FastQC-trimmed:$QC/trimmed/%s_trim_2_fastqc.zip" "fastp:$QC/%s.fastp.json"; do
  name=${spec%%:*}; pat=${spec#*:}; miss=""; n=0
  for s in $(cat "$REPO/metadata/srr_list.txt"); do
    f=${pat//%s/$s}; [[ -s "$f" ]] && n=$((n+1)) || miss="$miss $s"
  done
  [[ $n -eq 27 ]] && echo "OK: $name 27/27" || problem "$name $n/27, missing:$miss"
done

echo "== 3. sanity of the numbers"
out=$(for s in $(cat "$REPO/metadata/srr_list.txt"); do
  if [[ -s "$KAL/$s/run_info.json" ]]; then
    p=$(grep -o '"p_pseudoaligned": *[0-9.]*' "$KAL/$s/run_info.json" | grep -o '[0-9.]*$')
    awk -v p="$p" -v s="$s" 'BEGIN{ if (p<50 || p>85) print "PROBLEM: kallisto " s " pseudoaligned " p " % (expected 55-75)" }'
  fi
  if [[ -s "$MAP/$s/Log.final.out" ]]; then
    u=$(grep "Uniquely mapped reads %" "$MAP/$s/Log.final.out" | awk -F'|' '{gsub(/[ \t%]/, "", $2); print $2}')
    awk -v u="$u" -v s="$s" 'BEGIN{ if (u<70 || u>95) print "PROBLEM: STAR " s " uniquely mapped " u " % (expected 77-90)" }'
  fi
  if [[ -s "$RSEM/$s/$s.genes.results" ]]; then
    l=$(wc -l < "$RSEM/$s/$s.genes.results"); [[ $l -ge 70000 ]] || echo "PROBLEM: RSEM $s has only $l lines (expected ~78,290)"
  fi
done)
if [[ -n "$out" ]]; then echo "$out"; bad=$((bad + $(grep -c PROBLEM <<< "$out"))); else echo "OK: all kallisto, STAR and RSEM numbers in the expected range"; fi

echo "== 4. jobs that ended badly in the last 5 days"
sacct -S "$(date -d '5 days ago' +%F 2>/dev/null || date -v-5d +%F)" --format=JobID%20,JobName,State,Elapsed -X 2>/dev/null | grep -E "FAILED|TIMEOUT|OUT_OF_ME" || echo "none listed (cancelled jobs that were replaced by reruns are not shown here)"

echo
[[ $bad -eq 0 ]] && echo "ALL OK: run bash euler/04_collect.sh next" || echo "$bad problem line(s) above: fix or rerun those runs before collecting"
