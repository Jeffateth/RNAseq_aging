#!/usr/bin/env bash
# Download selected SRR runs (SRA format) with prefetch into raw_data/sra/.
# Usage:
#   ./download_sra.sh SRR21355398            # one run (test the workflow first)
#   ./download_sra.sh --all --yes            # all 27 runs (~50 GB) - requires --yes
# Nothing is downloaded without explicit arguments. Existing files are left untouched.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
LIST="$ROOT/metadata/srr_list.txt"
OUT="$ROOT/raw_data/sra"
mkdir -p "$OUT"
command -v prefetch >/dev/null || { echo "prefetch not found: brew install sratoolkit"; exit 1; }

if [[ "${1:-}" == "--all" ]]; then
  [[ "${2:-}" == "--yes" ]] || { echo "--all downloads ~50 GB. Re-run with: $0 --all --yes"; exit 1; }
  RUNS=$(cat "$LIST")
elif [[ $# -ge 1 ]]; then
  RUNS="$*"
else
  echo "Usage: $0 SRR... | --all --yes"; exit 1
fi

for r in $RUNS; do
  grep -qx "$r" "$LIST" || { echo "$r is not in the selected 27-run list, skipping"; continue; }
  if [[ -s "$OUT/$r/$r.sra" ]]; then echo "$r already present, skipping"; continue; fi
  echo "[$(date +%T)] prefetch $r"
  prefetch "$r" --output-directory "$OUT" --max-size 100G --progress
  vdb-validate "$OUT/$r/$r.sra" 2>&1 | tail -1
done
