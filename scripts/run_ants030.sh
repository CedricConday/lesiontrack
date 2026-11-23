#!/usr/bin/env bash
# ANTs backtests on P1 to P3 with the 0.3.0 defaults, so the engine table compares like with
# like. About 5.5 h each at 4 cores. Waits for run_slow.sh if it is still running.
# Usage: scripts/run_ants030.sh [threads]      log: work/logs/ants030.log, ends "ants030 done"
set -uo pipefail
cd "$(dirname "$0")/.."
export PYTHONPATH=$PWD
THREADS=${1:-4}
PY=~/.venvs/lesiontrack/bin/python
LT="$PY -m lesiontrack.cli"
MS=sourcedata/mslesseg/train
W=work/v031/backtests
LOG=work/logs/ants030.log
say() { echo "$(date '+%F %T') $*" | tee -a "$LOG"; }
while pgrep -f "[r]un_slow.sh" > /dev/null; do sleep 30; done
say "start, threads=$THREADS, $($LT --version)"
for s in P1 P2 P3; do
  out=$W/${s}_ants
  if [ -f "$out/backtest_metrics.json" ]; then say "${s}_ants: done"; continue; fi
  start=$(date +%s)
  $LT backtest --t1 $MS/$s/T1/${s}_T1_T1.nii.gz --flair $MS/$s/T1/${s}_T1_FLAIR.nii.gz \
    --mask $MS/$s/T1/${s}_T1_MASK.nii.gz --out "$out" --threads "$THREADS" --engine ants \
    > "work/logs/ants030_$s.log" 2>&1
  say "backtest ${s}_ants exit=$? seconds=$(( $(date +%s) - start ))"
done
say "ants030 done"
