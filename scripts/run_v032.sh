#!/usr/bin/env bash
# The P1-P3 greedy backtests again on 0.3.2 (lesion masks resliced without the similarity scale),
# from copies of the v031 runs so the cached registrations are reused and v031 stays untouched.
# Usage: scripts/run_v032.sh [threads]      log: work/logs/v032.log, ends "v032 done"
set -uo pipefail
cd "$(dirname "$0")/.."
export PYTHONPATH=$PWD
THREADS=${1:-4}
PY=~/.venvs/lesiontrack/bin/python
LT="$PY -m lesiontrack.cli"
MS=sourcedata/mslesseg/train
W=work/v032; R=results/v032
LOG=work/logs/v032.log
mkdir -p work/logs $W/backtests $R
say() { echo "$(date "+%F %T") $*" | tee -a "$LOG"; }
while pgrep -f "[r]un_cohorts_lt032.sh" > /dev/null; do sleep 60; done
say "start, threads=$THREADS, $($LT --version)"
for s in P1 P2 P3; do
  out=$W/backtests/${s}_greedy
  if [ -f "$out/backtest_metrics.json" ]; then say "$s: done"; continue; fi
  [ -d "$out" ] || cp -a work/v031/backtests/${s}_greedy "$out"
  rm -f "$out/backtest_metrics.json"
  start=$(date +%s)
  $LT backtest --t1 $MS/$s/T1/${s}_T1_T1.nii.gz --flair $MS/$s/T1/${s}_T1_FLAIR.nii.gz \
    --mask $MS/$s/T1/${s}_T1_MASK.nii.gz --out "$out" --threads "$THREADS" > work/logs/v032_$s.log 2>&1
  say "backtest ${s}_greedy exit=$? seconds=$(( $(date +%s) - start ))"
  cp "$out/backtest_metrics.json" $R/${s}_greedy_metrics.json 2>/dev/null
done
say "v032 done"
