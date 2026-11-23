#!/usr/bin/env python
"""Sensitivity of the SEL candidate rule against the mean-rate filter, by injected rate.

Pools the lesion tables of several backtests per condition and reports, for each nominal
injected rate, the fraction of expanded lesions that keep a candidate when candidates
under a mean rate are dropped, and the false positive rate on untouched lesions at the
same thresholds. A backtest scored before `max_candidate_mean_pct_per_year` existed is
re-scored from its stored volumes.

Usage: summarize_slow.py <out dir> <label>=<backtest dir>[,<backtest dir>...] ...
"""
import json
import sys
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd

from lesiontrack.synth import _auc, score

THRESHOLDS = (0.0, 2.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.5)


def _load(path: Path) -> np.ndarray:
    return np.asarray(nib.load(path).dataobj)


def lesions(bt: Path) -> pd.DataFrame:
    """Per-lesion rows of one backtest, with the best candidate mean rate per lesion."""
    table = pd.read_csv(bt / "backtest_lesions.tsv", sep="\t")
    if "max_candidate_mean_pct_per_year" in table:
        return table
    meta = json.loads((bt / "backtest_metrics.json").read_text())
    parts = []
    for seed in meta["seeds"]:
        s = bt / f"seed{seed}"
        truth = json.loads((s / "synthetic_1" / "truth.json").read_text())
        run = s / "run"
        t, _ = score(truth, _load(run / "truth_lesion_labels_halfway.nii.gz").astype(np.int32),
                     _load(run / "expansion_synthetic_1_pct_per_year.nii.gz"),
                     _load(run / "sel_candidates.nii.gz"), meta["dt_years"])
        old = table[table["seed"] == seed].set_index("lesion_id")
        t["detected_definite"] = t["lesion_id"].map(old["detected_definite"])
        t.insert(0, "seed", seed)
        parts.append(t)
    return pd.concat(parts, ignore_index=True)


def main() -> None:
    out = Path(sys.argv[1])
    out.mkdir(parents=True, exist_ok=True)
    by_rate, by_cond = [], []
    for spec in sys.argv[2:]:
        label, dirs = spec.split("=")
        dirs = [Path(d) for d in dirs.split(",")]
        t = pd.concat([lesions(d).assign(backtest=d.name) for d in dirs], ignore_index=True)
        t.to_csv(out / f"lesions_{label}.tsv", sep="\t", index=False, float_format="%.5g")
        metas = [json.loads((d / "backtest_metrics.json").read_text()) for d in dirs]
        cands = pd.concat([pd.read_csv(d / "backtest_candidates.tsv", sep="\t") for d in dirs], ignore_index=True)
        best = t["max_candidate_mean_pct_per_year"].fillna(-np.inf)
        exp, un = t[t["expanded"]], t[~t["expanded"]]
        cond = {"condition": label, "dt_years": metas[0]["dt_years"], "backtests": len(dirs),
                "n_expanded": len(exp), "n_untouched": len(un),
                "recovery_fraction_median": float(np.median([m["recovery_fraction"] for m in metas])),
                "untouched_mean_rate_sd": float(un["measured_mean_pct_per_year"].std()),
                "untouched_peak_median": float(un["measured_peak_pct_per_year"].median()),
                "auc_candidate_mean_rate": _auc(cands["mean_pct_per_year"], cands["injected"]),
                "auc_lesion_mean_rate": _auc(t["measured_mean_pct_per_year"], t["expanded"]),
                "fpr_definite": float(un["detected_definite"].mean())}
        for thr in THRESHOLDS:
            cond[f"fpr_min_mean_{thr:g}"] = float((best[un.index] >= thr).mean())
        by_cond.append(cond)
        for rate, g in exp.groupby(exp["nominal_pct_per_year"].round(1)):
            row = {"condition": label, "nominal_pct_per_year": rate, "n": len(g),
                   "true_pct_per_year_median": float(g["true_pct_per_year"].median()),
                   "measured_pct_per_year_median": float(g["measured_mean_pct_per_year"].median()),
                   "sens_definite": float(g["detected_definite"].mean())}
            for thr in THRESHOLDS:
                row[f"sens_min_mean_{thr:g}"] = float((best[g.index] >= thr).mean())
            by_rate.append(row)
    R, C = pd.DataFrame(by_rate), pd.DataFrame(by_cond)
    R.to_csv(out / "sensitivity_by_rate.tsv", sep="\t", index=False, float_format="%.3g")
    C.to_csv(out / "conditions.tsv", sep="\t", index=False, float_format="%.3g")
    with pd.option_context("display.width", 250, "display.max_columns", 30):
        print(R.to_string(index=False, float_format=lambda v: f"{v:.2f}"))
        print(C.set_index("condition").T.to_string(float_format=lambda v: f"{v:.2f}"))


if __name__ == "__main__":
    main()
