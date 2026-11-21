#!/usr/bin/env bash
# The benchmark runs behind docs/preprint: the same injected truth under a clinical 2D FLAIR,
# under a scanner change, and scored for segmentation-based detection with LST-AI masks.
# Waits for scripts/run_validation.sh (it uses every core), then runs one job at a time;
# finished steps are skipped on a re-run.
# Usage: scripts/run_bench.sh [threads]      log: work/logs/bench.log
set -uo pipefail
cd "$(dirname "$0")/.."
export PYTHONPATH=$PWD   # this checkout's lesiontrack, not the one installed in the venv
THREADS=${1:-4}
PY=~/.venvs/lesiontrack/bin/python
LT="$PY -m lesiontrack.cli"
MS=sourcedata/mslesseg/train
mkdir -p work/logs work/bench results/bench
LOG=work/logs/bench.log
say() { echo "$(date '+%F %T') $*" | tee -a "$LOG"; }

say "waiting for run_validation.sh"
while pgrep -f "[r]un_validation.sh" > /dev/null; do sleep 300; done
say "start, threads=$THREADS, $($LT --version)"

bt() {  # name subject extra-args...
  local name=$1 s=$2; shift 2
  local out=work/bench/$name
  if [ -f "$out/backtest_metrics.json" ]; then say "$name: done"; return; fi
  local start=$(date +%s)
  $LT backtest --t1 $MS/$s/T1/${s}_T1_T1.nii.gz --flair $MS/$s/T1/${s}_T1_FLAIR.nii.gz \
    --mask $MS/$s/T1/${s}_T1_MASK.nii.gz --out "$out" --threads "$THREADS" "$@" > "work/logs/bench_$name.log" 2>&1
  say "backtest $name exit=$? seconds=$(( $(date +%s) - start ))"
}

SUBJECTS="P1 P2 P3"
for s in $SUBJECTS; do bt ${s}_thick3 $s --flair-slices 3; done
for pct in 0.5 1.0; do for s in $SUBJECTS; do bt ${s}_scale$pct $s --scanner-scale-pct $pct; done; done

# The unperturbed backtests come from run_validation.sh; P1's is the re-scored v5 run.
base() { if [ "$1" = P1 ]; then echo work/backtest_P1_v5; else echo work/backtests/$1_greedy; fi; }
for s in $SUBJECTS P12; do
  d=$(base $s)
  [ -f "$d/backtest_metrics.json" ] || { say "lstai $s: no backtest at $d, skipped"; continue; }
  start=$(date +%s)
  $PY scripts/bench_lstai.py "$d" $MS/$s/T1/${s}_T1_T1.nii.gz $MS/$s/T1/${s}_T1_FLAIR.nii.gz "$THREADS" \
    > "work/logs/bench_lstai_$s.log" 2>&1
  say "lstai $s exit=$? seconds=$(( $(date +%s) - start ))"
done

existing() { for d in "$@"; do [ -f "$d/backtest_metrics.json" ] && echo "$d"; done; }
score() {  # label dirs...
  local label=$1; shift
  local dirs=$(existing "$@")
  [ -n "$dirs" ] || { say "score $label: nothing to score"; return; }
  $PY scripts/bench_methods.py results/bench/methods_$label $dirs > "work/logs/bench_methods_$label.log" 2>&1
  say "methods $label exit=$? ($(echo $dirs | wc -w) backtests)"
  $PY scripts/bench_rule.py results/bench/rule_$label $dirs > "work/logs/bench_rule_$label.log" 2>&1
  say "rule $label exit=$?"
}
score base work/backtest_P1_v5 work/backtests/P2_greedy work/backtests/P3_greedy work/backtests/P12_greedy \
  work/backtests/P13_greedy work/backtests/P19_greedy work/backtests/P22_greedy
score base3 work/backtest_P1_v5 work/backtests/P2_greedy work/backtests/P3_greedy
score ants work/backtests/P1_ants work/backtests/P2_ants work/backtests/P3_ants
score thick3 work/bench/P1_thick3 work/bench/P2_thick3 work/bench/P3_thick3
score scale0.5 work/bench/P1_scale0.5 work/bench/P2_scale0.5 work/bench/P3_scale0.5
score scale1.0 work/bench/P1_scale1.0 work/bench/P2_scale1.0 work/bench/P3_scale1.0
score lesjak work/backtests/lesjak_*_greedy
say "bench done"
