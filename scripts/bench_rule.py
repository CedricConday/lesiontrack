#!/usr/bin/env python
"""Compare Elliott's cohort-normalised definite-SEL rule with rules that need no cohort.

Elliott 2019 calls a candidate a definite SEL when S = z(concentricity) - z(constancy) >= 0,
with both z-scores taken over the candidates of the cohort at hand, so the same lesion can
change class with the company it keeps. A cohort-free rule is a fixed threshold on a
feature of the candidate itself. This script takes the candidate tables of several
backtests (injected or spurious is known for every candidate) and, leaving one subject
out at a time, picks the threshold with the best Youden index on the other subjects and
scores it on the held-out one. S >= 0 is scored on the same held-out candidates.
Usage: bench_rule.py <out_prefix> <backtest_dir> [<backtest_dir> ...]
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from lesiontrack.synth import CANDIDATE_FEATURES, _auc

FEATURES = ["mean_pct_per_year", "constancy_slope_pct_per_year", "peak_pct_per_year", "seed_voxels"]


def best_threshold(values: np.ndarray, positive: np.ndarray) -> float:
    """Threshold t maximising sensitivity - FPR of ``values >= t``."""
    cuts = np.unique(values)
    j = [(values[positive] >= t).mean() - (values[~positive] >= t).mean() for t in cuts]
    return float(cuts[int(np.argmax(j))])


def rates(flag: np.ndarray, positive: np.ndarray) -> dict:
    return {"sensitivity": float(flag[positive].mean()) if positive.any() else float("nan"),
            "false_positive_rate": float(flag[~positive].mean()) if (~positive).any() else float("nan"),
            "n_injected": int(positive.sum()), "n_spurious": int((~positive).sum())}


def main() -> int:
    prefix = Path(sys.argv[1])
    prefix.parent.mkdir(parents=True, exist_ok=True)
    frames = []
    for bt in map(Path, sys.argv[2:]):
        f = bt / "backtest_candidates.tsv"
        if f.exists():
            d = pd.read_csv(f, sep="\t")
            d.insert(0, "subject", bt.name)
            frames.append(d)
    if len(frames) < 2:
        print("need candidate tables from at least two backtests")
        return 1
    df = pd.concat(frames, ignore_index=True)
    df["injected"] = df["injected"].astype(bool)
    df["definite_sel"] = df["definite_sel"].fillna(False).astype(bool)
    pos = df["injected"].to_numpy()

    held = {"S >= 0 (Elliott)": df["definite_sel"].to_numpy().copy()}
    thresholds: dict = {f: {} for f in FEATURES}
    for feat in FEATURES:
        flag = np.zeros(len(df), bool)
        for subj in df["subject"].unique():
            test = (df["subject"] == subj).to_numpy()
            train = ~test & df[feat].notna().to_numpy()
            t = best_threshold(df.loc[train, feat].to_numpy(), pos[train])
            thresholds[feat][subj] = t
            flag[test] = df.loc[test, feat].fillna(-np.inf).to_numpy() >= t
        held[f"{feat} >= t"] = flag

    summary = {"n_subjects": int(df["subject"].nunique()), "n_candidates": len(df), "rules": {}, "auc": {}}
    for name, flag in held.items():
        entry = rates(flag, pos)
        per = [rates(flag[(df["subject"] == s).to_numpy()], pos[(df["subject"] == s).to_numpy()])
               for s in df["subject"].unique()]
        entry["sensitivity_range"] = [float(np.nanmin([p["sensitivity"] for p in per])),
                                      float(np.nanmax([p["sensitivity"] for p in per]))]
        entry["false_positive_rate_range"] = [float(np.nanmin([p["false_positive_rate"] for p in per])),
                                              float(np.nanmax([p["false_positive_rate"] for p in per]))]
        feat = name.split(" >= ")[0]
        if feat in thresholds:
            ts = list(thresholds[feat].values())
            entry["threshold_median"] = float(np.median(ts))
            entry["threshold_range"] = [float(min(ts)), float(max(ts))]
        summary["rules"][name] = entry
    for feat, sign in CANDIDATE_FEATURES.items():
        d = df[df[feat].notna()] if feat in df else df.iloc[:0]
        if len(d):
            summary["auc"][feat] = _auc(sign * d[feat].to_numpy(), d["injected"].to_numpy())
    Path(f"{prefix}_summary.json").write_text(json.dumps(summary, indent=2))
    df.to_csv(f"{prefix}_candidates.tsv", sep="\t", index=False, float_format="%.5g")
    for name, e in summary["rules"].items():
        print(f"{name:42s} sens {e['sensitivity']:.2f}  FPR {e['false_positive_rate']:.2f}"
              + (f"  t {e['threshold_median']:.3g}" if "threshold_median" in e else ""))
    print("pooled AUC: " + ", ".join(f"{k} {v:.2f}" for k, v in summary["auc"].items()))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
