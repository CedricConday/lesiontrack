#!/usr/bin/env python
"""Score deformation-based and segmentation-based SEL detection on the same injected truth.

Every backtest lesion has a known expansion. This script asks, lesion by lesion, which
method would have called it expanding:

  deformation / SEL candidate   Elliott's Jacobian rule: the lesion holds a SEL candidate
  deformation / definite SEL    ... that also passes the cohort-normalised S >= 0
  deformation / mean rate       mean Jacobian expansion over the lesion >= the tracking threshold
  segmentation / truth masks    lesion volume change between masks warped with the truth
  segmentation / LST-AI         the same on LST-AI masks of each image (after bench_lstai.py)

and writes sensitivity, false-positive rate and AUC of each method's continuous measure.
Usage: bench_methods.py <out_prefix> <backtest_dir> [<backtest_dir> ...]
"""
import json
import sys
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd

from lesiontrack.config import TrackParams
from lesiontrack.synth import _auc
from lesiontrack.tracking import track_pair

PARAMS = TrackParams()


def segmentation_arm(mask_bl: Path, mask_fu: Path, truth_labels: np.ndarray, dt: float, voxel_mm3: float) -> pd.DataFrame:
    """Per truth lesion: the volume rate and class of the tracking group it falls in."""
    table, maps = track_pair(np.asarray(nib.load(mask_bl).dataobj) > 0, np.asarray(nib.load(mask_fu).dataobj) > 0,
                             dt, voxel_mm3, PARAMS)
    group_of = {}
    for _, g in table.iterrows():
        for b in filter(None, str(g["baseline_ids"]).split(",")):
            group_of[int(b)] = g
    rows = []
    lb = maps["baseline"]
    for lid in np.unique(truth_labels[truth_labels > 0]):
        hits = lb[(truth_labels == lid) & (lb > 0)]
        if hits.size == 0:  # the segmenter did not see this lesion at baseline
            rows.append({"lesion_id": int(lid), "seen": False, "rate": 0.0, "enlarging": False})
            continue
        g = group_of[int(np.bincount(hits).argmax())]
        rows.append({"lesion_id": int(lid), "seen": True, "rate": float(g["pct_per_year"]),
                     "enlarging": g["class"] == "enlarging"})
    return pd.DataFrame(rows)


def one_seed(bt: Path, sdir: Path, dt: float) -> pd.DataFrame | None:
    reg = sdir / "run" / "reg"
    tags = sorted(p.name[: -len("_halfway.mat")] for p in reg.glob("synthetic_*_halfway.mat"))
    lesions_tsv = sdir / "backtest_lesions.tsv"
    if not tags or not lesions_tsv.exists():
        return None
    tag = tags[-1]
    truth_img = nib.load(sdir / "run" / "truth_lesion_labels_halfway.nii.gz")
    truth = np.asarray(truth_img.dataobj).astype(np.int32)
    voxel_mm3 = float(np.prod(truth_img.header.get_zooms()[:3]))
    out = pd.read_csv(lesions_tsv, sep="\t")
    out = out.rename(columns={"detected": "def_candidate", "detected_definite": "def_definite",
                              "measured_mean_pct_per_year": "def_mean_rate"})
    out["def_mean_flag"] = out["def_mean_rate"] >= PARAMS.change_pct_per_year
    arms = {"seg_truth": (reg / f"{tag}_baseline_mask_halfway.nii.gz", reg / f"{tag}_followup_mask_halfway.nii.gz"),
            "seg_lstai": (reg / "lstai_baseline_mask_halfway.nii.gz", reg / "lstai_followup_mask_halfway.nii.gz")}
    for arm, (mb, mf) in arms.items():
        if not (mb.exists() and mf.exists()):
            continue
        seg = segmentation_arm(mb, mf, truth, dt, voxel_mm3).rename(
            columns={"seen": f"{arm}_seen", "rate": f"{arm}_rate", "enlarging": f"{arm}_flag"})
        out = out.merge(seg, on="lesion_id", how="left")
    out.insert(0, "backtest", bt.name)
    return out


def summarize(df: pd.DataFrame) -> dict:
    pos, neg = df[df["expanded"]], df[~df["expanded"]]
    methods = {"deformation / SEL candidate": ("def_candidate", "measured_peak_pct_per_year"),
               "deformation / definite SEL": ("def_definite", None),
               "deformation / mean rate": ("def_mean_flag", "def_mean_rate"),
               "segmentation / truth masks": ("seg_truth_flag", "seg_truth_rate"),
               "segmentation / LST-AI": ("seg_lstai_flag", "seg_lstai_rate")}
    res = {"n_expanded": len(pos), "n_unexpanded": len(neg), "methods": {}}
    for name, (flag, measure) in methods.items():
        if flag not in df or df[flag].isna().all():
            continue
        ok = df[flag].notna()
        p, n = pos[ok[pos.index]], neg[ok[neg.index]]
        entry = {"n_expanded": len(p), "n_unexpanded": len(n),
                 "sensitivity": float(p[flag].astype(bool).mean()) if len(p) else float("nan"),
                 "false_positive_rate": float(n[flag].astype(bool).mean()) if len(n) else float("nan")}
        if measure and measure in df:
            d = df[ok & df[measure].notna()]
            entry["auc"] = _auc(d[measure].to_numpy(), d["expanded"].to_numpy(bool))
        res["methods"][name] = entry
    return res


def main() -> int:
    prefix = Path(sys.argv[1])
    prefix.parent.mkdir(parents=True, exist_ok=True)
    tables = []
    for bt in map(Path, sys.argv[2:]):
        meta = bt / "backtest_metrics.json"
        if not meta.exists():
            print(f"{bt}: no backtest_metrics.json, skipped")
            continue
        dt = float(json.loads(meta.read_text())["dt_years"])
        for sdir in sorted(bt.glob("seed*")):
            t = one_seed(bt, sdir, dt)
            if t is not None:
                tables.append(t)
    if not tables:
        print("no finished backtests")
        return 1
    df = pd.concat(tables, ignore_index=True)
    df.to_csv(f"{prefix}_lesions.tsv", sep="\t", index=False, float_format="%.5g")
    summary = {"pooled": summarize(df),
               "by_backtest": {k: summarize(g) for k, g in df.groupby("backtest")},
               "by_size": {k: summarize(g) for k, g in df.groupby("size_stratum")},
               "tracking_thresholds": PARAMS.__dict__}
    Path(f"{prefix}_summary.json").write_text(json.dumps(summary, indent=2))
    for name, m in summary["pooled"]["methods"].items():
        print(f"{name:30s} sens {m['sensitivity']:.2f}  FPR {m['false_positive_rate']:.2f}  "
              f"AUC {m.get('auc', float('nan')):.2f}  (n {m['n_expanded']}+{m['n_unexpanded']})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
