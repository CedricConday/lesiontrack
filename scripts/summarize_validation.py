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


def runs(backtests: Path | None) -> dict[str, Path]:
    """0.2.0 layout by default (P1 greedy lives in work/backtest_P1_v5); else one directory."""
    if backtests is not None:
        return {p.name: p for p in sorted(backtests.iterdir()) if p.is_dir()}  # a failed one is reported missing
    return {"P1_greedy": ROOT / "work/backtest_P1_v5",
            **{p.name: p for p in sorted((ROOT / "work/backtests").iterdir()) if p.is_dir()}}


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


def main(argv: list[str] | None = None) -> int:
    import argparse

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default=str(ROOT / "results"), help="directory for the two tables")
    ap.add_argument("--backtests", type=Path, default=None, help="directory holding one subdirectory per backtest")
    a = ap.parse_args(argv)
    rows, aucs = [], []
    for name, d in runs(a.backtests).items():
        f = d / "backtest_metrics.json"
        if not f.exists():
            print(f"missing {f}", file=sys.stderr)
            return 1
        m = json.loads(f.read_text())
        subject, engine = name.rsplit("_", 1)  # engine is the condition for bench runs (scale1.0, thick3)
        dataset = "lesjak" if subject.startswith("lesjak") else "mslesseg"
        scales = [x for seed in m.get("per_seed", []) for x in seed.get("residual_global_scale", {}).values()]
        rows.append({"backtest": name, "dataset": dataset, "subject": subject,
                     "engine": engine, **{k: m[k] for k in FIELDS},
                     "residual_global_scale_median": float(sorted(scales)[len(scales) // 2]) if scales else ""})
        fa = m["feature_auc"]
        aucs.append({"backtest": name, **{k: fa[k]["auc"] for k in FEATURES}})
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
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
