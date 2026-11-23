#!/usr/bin/env bash
# Backtests at the rates slowly expanding lesions are defined by, not the 15 to 40 % of the
# default design: nominal 3, 5 and 8 %/yr, three scans, over two years and over one.
# P1 to P3, three seeds each, greedy, every candidate scored (--min-mean 0).
# Usage: scripts/run_slow.sh [threads]      log: work/logs/slow.log, ends "slow done"
set -uo pipefail
cd "$(dirname "$0")/.."
export PYTHONPATH=$PWD
THREADS=${1:-4}
PY=~/.venvs/lesiontrack/bin/python
LT="$PY -m lesiontrack.cli"
MS=sourcedata/mslesseg/train
W=work/slow
mkdir -p work/logs $W results/slow
LOG=work/logs/slow.log
say() { echo "$(date '+%F %T') $*" | tee -a "$LOG"; }
say "start, threads=$THREADS, $($LT --version)"
bt() {  # name subject extra-args...
  local name=$1 s=$2; shift 2
  local out=$W/$name
  if [ -f "$out/backtest_metrics.json" ]; then say "$name: done"; return; fi
  local start=$(date +%s)
  $LT backtest --t1 $MS/$s/T1/${s}_T1_T1.nii.gz --flair $MS/$s/T1/${s}_T1_FLAIR.nii.gz \
    --mask $MS/$s/T1/${s}_T1_MASK.nii.gz --out "$out" --threads "$THREADS" "$@" > "work/logs/slow_$name.log" 2>&1
  say "backtest $name exit=$? seconds=$(( $(date +%s) - start ))"
}
for s in P1 P2 P3; do bt ${s}_2y $s --factors 1.06,1.10,1.16 --dt-years 2 --timepoints 0.5,1.0; done
for s in P1 P2 P3; do bt ${s}_1y $s --factors 1.03,1.05,1.08 --dt-years 1 --timepoints 0.5,1.0; done
say "slow done"
