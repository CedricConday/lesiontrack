#!/usr/bin/env bash
# Every greedy result in the README, re-run with the 0.3.0 defaults (similarity linear step,
# mean-rate filter 10 %/yr) into fresh directories, so the 0.2.0 outputs stay for comparison.
# The ANTs backtests (about 5.5 h each) and the LST-AI arm are not re-run.
# Usage: scripts/run_v030.sh [threads]      log: work/logs/v030.log
set -uo pipefail
cd "$(dirname "$0")/.."
export PYTHONPATH=$PWD
THREADS=${1:-4}
PY=~/.venvs/lesiontrack/bin/python
LT="$PY -m lesiontrack.cli"
MS=sourcedata/mslesseg/train
W=work/v030
R=results/v030
mkdir -p work/logs $W $R
LOG=work/logs/v030.log
say() { echo "$(date '+%F %T') $*" | tee -a "$LOG"; }
say "start, threads=$THREADS, $($LT --version)"

bt() {  # name t1 flair mask extra-args...
  local name=$1 t1=$2 fl=$3 mk=$4; shift 4
  local out=$W/backtests/$name
  if [ -f "$out/backtest_metrics.json" ]; then say "$name: done"; return; fi
  local start=$(date +%s)
  $LT backtest --t1 "$t1" --flair "$fl" --mask "$mk" --out "$out" --threads "$THREADS" "$@" \
    > "work/logs/v030_$name.log" 2>&1
  say "backtest $name exit=$? seconds=$(( $(date +%s) - start ))"
}
ms() { local name=$1 s=$2; shift 2; bt "$name" $MS/$s/T1/${s}_T1_T1.nii.gz $MS/$s/T1/${s}_T1_FLAIR.nii.gz $MS/$s/T1/${s}_T1_MASK.nii.gz "$@"; }

# Scanner change first: it is the question this run answers.
for pct in 1.0 0.5; do for s in P1 P2 P3; do ms ${s}_scale$pct $s --scanner-scale-pct $pct; done; done
for s in P1 P2 P3 P12 P13 P19 P22; do ms ${s}_greedy $s; done
for s in P1 P2 P3; do ms ${s}_thick3 $s --flair-slices 3; done
for s in patient03 patient06 patient07; do
  d=work/lesjak/input/$s
  bt lesjak_${s}_greedy $d/study1_T1w.nii.gz $d/study1_FLAIR.nii.gz $d/study1_mask.nii.gz
done

existing() { for d in "$@"; do [ -f "$d/backtest_metrics.json" ] && echo "$d"; done; }
score() {  # label dirs...
  local label=$1; shift
  local dirs=$(existing "$@")
  [ -n "$dirs" ] || { say "score $label: nothing to score"; return; }
  $PY scripts/bench_methods.py $R/methods_$label $dirs > "work/logs/v030_methods_$label.log" 2>&1
  say "methods $label exit=$?"
  $PY scripts/bench_rule.py $R/rule_$label $dirs > "work/logs/v030_rule_$label.log" 2>&1
  say "rule $label exit=$?"
}
B=$W/backtests
score base $B/P{1,2,3,12,13,19,22}_greedy
score base3 $B/P{1,2,3}_greedy
score thick3 $B/P{1,2,3}_thick3
score scale0.5 $B/P{1,2,3}_scale0.5
score scale1.0 $B/P{1,2,3}_scale1.0
score lesjak $B/lesjak_*_greedy
$PY scripts/summarize_validation.py --backtests $B --out $R > work/logs/v030_summary.log 2>&1
say "summary exit=$?"

# Cohorts: MSLesSeg (24 patients) and Lesjak (20), same commands as 0.2.0.
D=derivatives/v030
for s in $(tail -n +2 work/manifest_mslesseg.tsv | cut -f1 | sort -u -V); do
  [ -f "$D/mslesseg/$s/summary.json" ] && continue
  start=$(date +%s)
  $LT run work/manifest_mslesseg.tsv --out $D/mslesseg --subject "$s" --threads "$THREADS" > "work/logs/v030_ms_$s.log" 2>&1
  say "mslesseg $s exit=$? seconds=$(( $(date +%s) - start ))"
done
$LT cohort --out $D/mslesseg > work/logs/v030_ms_cohort.log 2>&1; say "mslesseg cohort exit=$?"
$PY scripts/summarize_cohort.py $D/mslesseg > work/logs/v030_ms_summary.log 2>&1; say "summarize_cohort exit=$?"
$PY scripts/cohort_checks.py $D/mslesseg > work/logs/v030_ms_checks.log 2>&1; say "cohort_checks exit=$?"
for s in $(tail -n +2 work/lesjak/manifest_lesjak.tsv | cut -f1 | sort -u); do
  [ -f "$D/lesjak/$s/summary.json" ] && continue
  start=$(date +%s)
  $LT run work/lesjak/manifest_lesjak.tsv --out $D/lesjak --subject "$s" --threads "$THREADS" > "work/logs/v030_lesjak_$s.log" 2>&1
  say "lesjak $s exit=$? seconds=$(( $(date +%s) - start ))"
done
$LT cohort --out $D/lesjak > work/logs/v030_lesjak_cohort.log 2>&1; say "lesjak cohort exit=$?"
$PY scripts/eval_lesjak.py $D/lesjak sourcedata/lesjak/longitudinal $R/lesjak > work/logs/v030_eval_lesjak.log 2>&1
say "eval_lesjak exit=$?"
say "v030 done"
