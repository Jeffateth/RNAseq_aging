#!/usr/bin/env bash
# OPTIONAL: convert downloaded .sra to gzipped FASTQ (paired-end -> _1/_2) in raw_data/fastq/.
# Usage: ./convert_fastq.sh SRR21355398 [SRR...]
# Needs ~7x the .sra size in free scratch space per run during conversion; runs are done one at a time.
# The .sra file is never modified or deleted.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SRA="$ROOT/raw_data/sra"; OUT="$ROOT/raw_data/fastq"; TMP="$ROOT/raw_data/tmp"
THREADS="${THREADS:-6}"
[[ $# -ge 1 ]] || { echo "Usage: $0 SRR..."; exit 1; }
command -v fasterq-dump >/dev/null || { echo "fasterq-dump not found"; exit 1; }
mkdir -p "$OUT" "$TMP"
for r in "$@"; do
  [[ -s "$SRA/$r/$r.sra" ]] || { echo "$r.sra missing; run download_sra.sh first"; continue; }
  [[ -s "$OUT/${r}_1.fastq.gz" ]] && { echo "$r FASTQ exists, skipping"; continue; }
  echo "[$(date +%T)] fasterq-dump $r"
  fasterq-dump "$SRA/$r/$r.sra" --split-3 --threads "$THREADS" --temp "$TMP" --outdir "$OUT" --progress
  gzip -f "$OUT/${r}"*.fastq
done
rmdir "$TMP" 2>/dev/null || true
