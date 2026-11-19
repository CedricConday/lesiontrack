#!/usr/bin/env bash
# The validation runs behind README "Results": ANTs against greedy on the same synthetic
# series, the backtest on more than one MSLesSeg patient, and the Lesjak 2016 database.
# One job at a time (each uses every core); finished steps are skipped on a re-run.
# Usage: scripts/run_validation.sh [threads]      log: work/logs/validation.log
set -uo pipefail
cd "$(dirname "$0")/.."
THREADS=${1:-4}
LT=~/.venvs/lesiontrack/bin/lesiontrack
PY=~/.venvs/lesiontrack/bin/python
MS=sourcedata/mslesseg/train
mkdir -p work/logs work/backtests
LOG=work/logs/validation.log
say() { echo "$(date '+%F %T') $*" | tee -a "$LOG"; }

backtest() {  # name t1 flair mask engine
  local out=work/backtests/$1_$5
  if [ -f "$out/backtest_metrics.json" ]; then say "$1 $5: done"; return; fi
  local start=$(date +%s)
  $LT backtest --t1 "$2" --flair "$3" --mask "$4" --out "$out" --engine "$5" --threads "$THREADS" \
    > "work/logs/backtest_$1_$5.log" 2>&1
  say "backtest $1 $5 exit=$? seconds=$(( $(date +%s) - start ))"
}
ms() { backtest "$1" $MS/$1/T1/$1_T1_T1.nii.gz $MS/$1/T1/$1_T1_FLAIR.nii.gz $MS/$1/T1/$1_T1_MASK.nii.gz "$2"; }

say "start, threads=$THREADS, $($LT --version)"
ms P1 ants
for s in P2 P3 P12 P13 P19 P22; do ms $s greedy; done

start=$(date +%s)
$PY scripts/prep_lesjak.py sourcedata/lesjak/longitudinal work/lesjak "$THREADS" > work/logs/prep_lesjak.log 2>&1
say "prep_lesjak exit=$? seconds=$(( $(date +%s) - start ))"
for s in $(tail -n +2 work/lesjak/manifest_lesjak.tsv | cut -f1 | sort -u); do
  [ -f "derivatives/lesjak/$s/summary.json" ] && continue
  start=$(date +%s)
  $LT run work/lesjak/manifest_lesjak.tsv --out derivatives/lesjak --subject "$s" --threads "$THREADS" \
    > "work/logs/lesjak_$s.log" 2>&1
  say "lesjak $s exit=$? seconds=$(( $(date +%s) - start ))"
done
$LT cohort --out derivatives/lesjak > work/logs/lesjak_cohort.log 2>&1
$PY scripts/eval_lesjak.py derivatives/lesjak sourcedata/lesjak/longitudinal results/lesjak > work/logs/eval_lesjak.log 2>&1
say "eval_lesjak exit=$?"

# Backtest at clinical resolution: the first three Lesjak patients with enough lesions.
n=0
for s in $(tail -n +2 work/lesjak/manifest_lesjak.tsv | cut -f1 | sort -u); do
  d=work/lesjak/input/$s
  ok=$($PY - "$d/study1_mask.nii.gz" <<'P'
import sys, nibabel as nib, numpy as np
from lesiontrack.tracking import label_lesions
sz = np.bincount(label_lesions(np.asarray(nib.load(sys.argv[1]).dataobj) > 0).ravel())[1:]
print(int((sz >= 30).sum() >= 12))
P
)
  [ "$ok" = 1 ] || continue
  backtest "lesjak_$s" $d/study1_T1w.nii.gz $d/study1_FLAIR.nii.gz $d/study1_mask.nii.gz greedy
  n=$((n + 1)); [ $n -ge 3 ] && break
done

for s in P2 P3; do ms $s ants; done
say "all done"
