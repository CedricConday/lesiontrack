"""Tabulate the validation backtests (run_validation.sh) into results/.

Writes results/validation_backtests.tsv (one row per backtest) and
results/validation_feature_auc.tsv (AUC of each candidate feature as a classifier of
injected against spurious candidates, one row per backtest).
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RUNS = {
    "P1_greedy": ROOT / "work/backtest_P1_v5",
    **{p.name: p for p in sorted((ROOT / "work/backtests").iterdir()) if p.is_dir()},
}
FIELDS = [
    "sensitivity", "false_positive_rate", "sensitivity_definite",
    "false_positive_rate_definite", "recovery_fraction", "recovery_median",
    "untouched_peak_median_pct_per_year", "n_expanded", "n_unexpanded",
    "n_candidates", "n_candidates_injected",
]
FEATURES = [
    "mean_pct_per_year", "constancy_slope_pct_per_year", "peak_pct_per_year",
    "concentricity", "constancy_residual", "score_S", "seed_voxels", "n_voxels",
]


def main() -> int:
    rows, aucs = [], []
    for name, d in RUNS.items():
        f = d / "backtest_metrics.json"
        if not f.exists():
            print(f"missing {f}", file=sys.stderr)
            return 1
        m = json.loads(f.read_text())
        subject, engine = name.rsplit("_", 1)
        dataset = "lesjak" if subject.startswith("lesjak") else "mslesseg"
        rows.append({"backtest": name, "dataset": dataset, "subject": subject,
                     "engine": engine, **{k: m[k] for k in FIELDS}})
        fa = m["feature_auc"]
        aucs.append({"backtest": name, **{k: fa[k]["auc"] for k in FEATURES}})
    out = ROOT / "results"
    for path, data in ((out / "validation_backtests.tsv", rows),
                       (out / "validation_feature_auc.tsv", aucs)):
        with path.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(data[0]), delimiter="\t")
            w.writeheader()
            for r in data:
                w.writerow({k: (f"{v:.3f}" if isinstance(v, float) else v) for k, v in r.items()})
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
