#!/usr/bin/env python
"""Score a lesiontrack run on the Lesjak 2016 database against its expert change masks.

The expert mask marks voxels where a lesion appeared, grew, shrank or disappeared between
the two studies, without saying which. Lesion masks here come from LST-AI, one per study,
so every number below carries the segmenter's errors as well as the tracker's.

Usage: eval_lesjak.py <lesiontrack output dir> [database_dir] [results_dir]
"""
import json
import sys
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd
from scipy import ndimage as ndi
from scipy.stats import spearmanr

from lesiontrack.registration import run_greedy
from lesiontrack.synth import _auc

out = Path(sys.argv[1])
db = Path(sys.argv[2] if len(sys.argv) > 2 else "sourcedata/lesjak/longitudinal")
res_dir = Path(sys.argv[3] if len(sys.argv) > 3 else "results/lesjak")
res_dir.mkdir(parents=True, exist_ok=True)
ST = ndi.generate_binary_structure(3, 2)
CHANGED = {"new", "resolved", "enlarging", "shrinking"}
UNCHANGED = {"stable", "trend_up", "trend_down"}


def load(p):
    return np.asarray(nib.load(p).dataobj)


subjects, groups, lesions, comps = [], [], [], []
for sdir in sorted(out.glob("patient*")):
    if not (sdir / "summary.json").exists():
        continue
    name, fu = sdir.name, "study2"
    reg = sdir / "reg"
    gt_h = sdir / "expert_change_halfway.nii.gz"
    run_greedy(f"-d 3 -rf {reg / f'{fu}_baseline_T1w_halfway.nii.gz'} -ri NN -rt uchar "
               f"-rm {db / 'coregistered' / name / 'gt.nii.gz'} {gt_h} -r {reg / f'{fu}_halfway.mat'},-1")
    gt = load(gt_h) > 0
    lb = load(sdir / "baseline_lesion_labels.nii.gz").astype(np.int32)
    lf = load(sdir / f"lesion_labels_{fu}.nii.gz").astype(np.int32)
    exp = load(sdir / f"expansion_{fu}_pct_per_year.nii.gz")
    cand = load(sdir / "sel_candidates.nii.gz")
    track = pd.read_csv(sdir / "lesion_tracking.tsv", sep="\t", dtype={"baseline_ids": str, "followup_ids": str})
    xor = (lb > 0) ^ (lf > 0)
    tp = int((xor & gt).sum())
    gl, n_gt = ndi.label(gt, structure=ST)
    gt_sizes = np.bincount(gl.ravel())[1:]
    hit = np.zeros(n_gt + 1, bool)
    hit[np.unique(gl[xor])] = True
    in_mask = np.zeros(n_gt + 1, bool)
    in_mask[np.unique(gl[(lb > 0) | (lf > 0)])] = True
    for k in range(1, n_gt + 1):
        comps.append({"subject": name, "gt_component": k, "n_voxels": int(gt_sizes[k - 1]),
                      "hit_by_mask_change": bool(hit[k]), "inside_either_mask": bool(in_mask[k])})
    xl, n_x = ndi.label(xor, structure=ST)
    x_sizes = np.bincount(xl.ravel())
    big = np.flatnonzero(x_sizes >= 3)
    big = big[big != 0]
    x_true = np.isin(big, np.unique(xl[gt]))

    gt_near = ndi.binary_dilation(gt, structure=ST)
    for _, r in track.iterrows():
        b = [int(x) for x in str(r["baseline_ids"]).split(",") if x and x != "nan"]
        f = [int(x) for x in str(r["followup_ids"]).split(",") if x and x != "nan"]
        vox = np.isin(lb, b) | np.isin(lf, f)
        groups.append({"subject": name, "group_id": r["group_id"], "class": r["class"],
                       "n_voxels": int(vox.sum()), "pct_per_year": r["pct_per_year"],
                       "expert_changed": bool(gt_near[vox].any()),
                       "flair_dz": r.get("flair_dz", float("nan")),
                       "flair_supported": r.get("flair_supported", float("nan")),
                       "jacobian_mean_pct_per_year": float(exp[np.isin(lb, b)].mean()) if b else float("nan")})
    for k in range(1, int(lb.max()) + 1):
        sel = lb == k
        if sel.sum() < 10:
            continue
        v = exp[sel]
        lesions.append({"subject": name, "lesion_id": k, "n_voxels": int(sel.sum()),
                        "expert_changed": bool(gt_near[sel].any()),
                        "expert_change_fraction": float(gt[sel].mean()),
                        "jacobian_mean": float(v.mean()), "jacobian_abs_mean": float(abs(v.mean())),
                        "jacobian_peak": float(v.max()), "has_candidate": bool((cand[sel] > 0).any())})
    n_c = int(cand.max())
    c_gt = sum(bool(gt_near[cand == k].any()) for k in range(1, n_c + 1))
    summ = json.loads((sdir / "summary.json").read_text())
    subjects.append({
        "subject": name, "dt_years": summ["dt_years"][-1],
        "baseline_lesions": int(lb.max()), "followup_lesions": int(lf.max()),
        "expert_change_voxels": int(gt.sum()), "mask_change_voxels": int(xor.sum()),
        "voxel_dice": 2 * tp / (xor.sum() + gt.sum()) if (xor.sum() + gt.sum()) else float("nan"),
        "voxel_sensitivity": tp / gt.sum() if gt.sum() else float("nan"),
        "voxel_ppv": tp / xor.sum() if xor.sum() else float("nan"),
        "expert_components": n_gt, "expert_components_hit": int(hit[1:].sum()),
        "expert_components_inside_either_mask": int(in_mask[1:].sum()),
        "mask_change_components_ge3": int(big.size), "mask_change_components_on_expert": int(x_true.sum()),
        "sel_candidates": n_c, "sel_candidates_on_expert_change": int(c_gt),
    })
    print(f"{name}: dice {subjects[-1]['voxel_dice']:.2f}, expert components hit "
          f"{subjects[-1]['expert_components_hit']}/{n_gt}, candidates {n_c}", flush=True)

S, G, L, C = map(pd.DataFrame, (subjects, groups, lesions, comps))
for df, nm in ((S, "by_subject"), (G, "by_group"), (L, "by_baseline_lesion"), (C, "by_expert_component")):
    df.to_csv(res_dir / f"lesjak_{nm}.tsv", sep="\t", index=False, float_format="%.5g")
scored = G[G["class"].isin(CHANGED | UNCHANGED)]
called = scored["class"].isin(CHANGED)
ec = scored["expert_changed"]
# The same call, with new and resolved counted only where FLAIR intensity moved with them.
sup = scored["flair_supported"].astype("boolean").fillna(False).astype(bool)
called_sup = scored["class"].isin({"enlarging", "shrinking"}) | (scored["class"].isin({"new", "resolved"}) & sup)
matched = G[G["jacobian_mean_pct_per_year"].notna() & G["class"].ne("resolved") & G["pct_per_year"].notna()]
summary = {
    "subjects": len(S),
    "voxel": {"dice_median": float(S["voxel_dice"].median()), "dice_pooled": float(
        2 * (S["voxel_sensitivity"] * S["expert_change_voxels"]).sum() / (S["mask_change_voxels"].sum() + S["expert_change_voxels"].sum())),
        "sensitivity_pooled": float((S["voxel_sensitivity"] * S["expert_change_voxels"]).sum() / S["expert_change_voxels"].sum()),
        "ppv_pooled": float((S["voxel_ppv"] * S["mask_change_voxels"]).sum() / S["mask_change_voxels"].sum())},
    "expert_components": {"n": int(S["expert_components"].sum()), "hit": int(S["expert_components_hit"].sum()),
                          "inside_either_mask": int(S["expert_components_inside_either_mask"].sum()),
                          "hit_by_size": {f"{lo}-{hi}": [int(C[(C.n_voxels >= lo) & (C.n_voxels < hi)]["hit_by_mask_change"].sum()),
                                                         int(((C.n_voxels >= lo) & (C.n_voxels < hi)).sum())]
                                          for lo, hi in ((1, 10), (10, 50), (50, 10**9))}},
    "mask_change_components_ge3": {"n": int(S["mask_change_components_ge3"].sum()),
                                   "on_expert_change": int(S["mask_change_components_on_expert"].sum())},
    "tracking_groups": {
        "n_scored": len(scored), "adjacent_fragments_not_scored": int((G["class"] == "adjacent_fragment").sum()),
        "expert_changed": int(ec.sum()),
        "sensitivity_called_changed": float(called[ec].mean()) if ec.any() else float("nan"),
        "specificity_called_unchanged": float((~called[~ec]).mean()) if (~ec).any() else float("nan"),
        "ppv_called_changed": float(ec[called].mean()) if called.any() else float("nan"),
        "flair_supported": {
            "sensitivity_called_changed": float(called_sup[ec].mean()) if ec.any() else float("nan"),
            "specificity_called_unchanged": float((~called_sup[~ec]).mean()) if (~ec).any() else float("nan"),
            "ppv_called_changed": float(ec[called_sup].mean()) if called_sup.any() else float("nan"),
            "by_class": {c: {"n": int((scored["class"] == c).sum()), "supported": int(((scored["class"] == c) & sup).sum()),
                             "expert_changed": int(((scored["class"] == c) & ec).sum()),
                             "supported_and_expert_changed": int(((scored["class"] == c) & sup & ec).sum())}
                         for c in ("new", "resolved")}},
        "by_class": {c: {"n": int((G["class"] == c).sum()), "expert_changed": int(G[G["class"] == c]["expert_changed"].sum())}
                     for c in sorted(G["class"].unique())}},
    "jacobian_per_baseline_lesion_ge10_voxels": {
        "n": len(L), "expert_changed": int(L["expert_changed"].sum()),
        "auc_abs_mean": _auc(L["jacobian_abs_mean"], L["expert_changed"]),
        "auc_signed_mean": _auc(L["jacobian_mean"], L["expert_changed"]),
        "auc_peak": _auc(L["jacobian_peak"], L["expert_changed"]),
        "auc_has_sel_candidate": _auc(L["has_candidate"].astype(float), L["expert_changed"]),
        "spearman_jacobian_vs_mask_rate_matched_groups": float(spearmanr(matched["jacobian_mean_pct_per_year"], matched["pct_per_year"]).statistic) if len(matched) > 2 else float("nan"),
        "n_matched_groups": len(matched)},
    "sel_candidates": {"n": int(S["sel_candidates"].sum()), "median_per_subject": float(S["sel_candidates"].median()),
                       "on_expert_change": int(S["sel_candidates_on_expert_change"].sum())},
}
(res_dir / "lesjak_summary.json").write_text(json.dumps(summary, indent=2))
print(json.dumps(summary, indent=2))
