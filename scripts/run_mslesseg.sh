#!/usr/bin/env bash
# Run every longitudinal MSLesSeg subject, one at a time, logging to work/logs.
# Usage: scripts/run_mslesseg.sh [threads]
set -uo pipefail
cd "$(dirname "$0")/.."
THREADS=${1:-4}
MANIFEST=work/manifest_mslesseg.tsv
OUT=derivatives/mslesseg
mkdir -p work/logs "$OUT"
for s in $(tail -n +2 "$MANIFEST" | cut -f1 | sort -u -V); do
  if [ -f "$OUT/$s/summary.json" ]; then echo "$s done"; continue; fi
  start=$(date +%s)
  ~/.venvs/lesiontrack/bin/lesiontrack run "$MANIFEST" --out "$OUT" --subject "$s" --threads "$THREADS" \
    > "work/logs/$s.log" 2>&1
  echo "$s exit=$? seconds=$(( $(date +%s) - start ))" | tee -a work/logs/run_mslesseg.log
done
