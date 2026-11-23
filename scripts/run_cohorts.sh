#!/usr/bin/env bash
# Both cohorts with the current defaults, in place over derivatives/v031 so the cached
# registrations are reused (the cache needs the same thread count as the first run: 4).
# Usage: scripts/run_cohorts.sh [results dir]     log: work/logs/cohorts.log, ends "cohorts done"
set -uo pipefail
cd "$(dirname "$0")/.."
export PYTHONPATH=$PWD
PY=~/.venvs/lesiontrack/bin/python
LT="$PY -m lesiontrack.cli"
D=derivatives/v031
R=${1:-results/v0.3.0}
mkdir -p work/logs $R/lesjak $R/mslesseg
LOG=work/logs/cohorts.log
say() { echo "$(date '+%F %T') $*" | tee -a "$LOG"; }
say "start, $($LT --version)"
for c in lesjak mslesseg; do
  man=work/manifest_mslesseg.tsv; [ $c = lesjak ] && man=work/lesjak/manifest_lesjak.tsv
  for s in $(tail -n +2 $man | cut -f1 | sort -u -V); do
    start=$(date +%s)
    $LT run $man --out $D/$c --subject "$s" --threads 4 > "work/logs/cohorts_${c}_$s.log" 2>&1
    say "$c $s exit=$? seconds=$(( $(date +%s) - start ))"
  done
  $LT cohort --out $D/$c > work/logs/cohorts_${c}_cohort.log 2>&1; say "$c cohort exit=$?"
done
$PY scripts/summarize_cohort.py $D/mslesseg > work/logs/cohorts_ms_summary.log 2>&1; say "summarize_cohort exit=$?"
$PY scripts/cohort_checks.py $D/mslesseg > work/logs/cohorts_ms_checks.log 2>&1; say "cohort_checks exit=$?"
cp $D/mslesseg/{cohort_summary.tsv,sel_cohort_by_subject.tsv,check_*.tsv} $R/mslesseg/ 2>>"$LOG"
$PY scripts/eval_lesjak.py $D/lesjak sourcedata/lesjak/longitudinal $R/lesjak > work/logs/cohorts_eval_lesjak.log 2>&1
say "eval_lesjak exit=$?"
say "cohorts done"
