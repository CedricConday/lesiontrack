"""Cohort summary for a lesiontrack output tree: counts, classes, intervals, and the
agreement between mask-based enlargement and Jacobian-based SEL candidates per lesion."""

import json
import sys
from pathlib import Path

import pandas as pd

out = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("derivatives/mslesseg")
subjects = sorted(p.parent for p in out.glob("*/summary.json"))
rows, tracks, cands = [], [], []
for d in subjects:
    s = json.load(open(d / "summary.json"))
    last = s["follow_ups"][-1]
    t = pd.read_csv(d / "lesion_tracking.tsv", sep="\t")
    c = pd.read_csv(d / "sel_candidates.tsv", sep="\t")
    tl = t[t["follow_up"] == last]
    rows.append({
        "subject": s["subject"], "n_timepoints": len(s["follow_ups"]) + 1,
        "interval_years": round(s["dt_years"][-1], 2),
        "baseline_lesions": s["baseline_lesion_count"],
        "baseline_volume_ml": round(s["baseline_lesion_volume_mm3"] / 1000, 2),
        "volume_change_pct_per_year": round((tl["volume_followup_mm3"].sum() / tl["volume_baseline_mm3"].sum() - 1) * 100 / s["dt_years"][-1], 1),
        **{k: int(v) for k, v in s["tracking"][last]["classes"].items()},
        "sel_candidates": s["sel_candidates"],
        "sel_candidates_mean_ge8": int((c["mean_pct_per_year"] >= 8).sum()) if len(c) else 0,
        "sel_volume_ml": round(s["sel_candidate_volume_mm3"] / 1000, 2),
    })
    tracks.append(tl.assign(subject=s["subject"]))
    cands.append(c)
table = pd.DataFrame(rows).fillna(0)
for col in ("new", "resolved", "enlarging", "shrinking", "stable", "trend_up", "trend_down", "adjacent_fragment"):
    if col not in table:
        table[col] = 0
    table[col] = table[col].astype(int)
table.to_csv(out / "cohort_summary.tsv", sep="\t", index=False)
pd.set_option("display.width", 250)
print(table.to_string(index=False))
print()
print("subjects:", len(table), " with >= 1 SEL candidate:", int((table.sel_candidates > 0).sum()),
      " with >= 1 candidate of mean >= 8 %/yr:", int((table.sel_candidates_mean_ge8 > 0).sum()))
print("median per subject: candidates", table.sel_candidates.median(), " mean>=8", table.sel_candidates_mean_ge8.median(),
      " SEL volume ml", table.sel_volume_ml.median(), " interval y", table.interval_years.median())
print("lesion classes at last follow-up (all subjects):", {k: int(table[k].sum()) for k in ("new", "resolved", "enlarging", "shrinking", "stable", "trend_up", "trend_down")})

# Agreement: does a baseline lesion that the masks call enlarging carry a SEL candidate?
allc = pd.concat(cands, ignore_index=True)
allt = pd.concat(tracks, ignore_index=True)
allt = allt[allt["n_baseline"] == 1].copy()
allt["lesion_id"] = allt["baseline_ids"].astype(int)
sel_parents = set(zip(allc["subject"], allc["parent_lesion_id"]))
allt["has_candidate"] = [(s, l) in sel_parents for s, l in zip(allt["subject"], allt["lesion_id"])]
strong = allc[allc["mean_pct_per_year"] >= 8]
strong_parents = set(zip(strong["subject"], strong["parent_lesion_id"]))
allt["has_strong_candidate"] = [(s, l) in strong_parents for s, l in zip(allt["subject"], allt["lesion_id"])]
ag = allt.groupby("class")[["has_candidate", "has_strong_candidate"]].mean().round(2)
ag["n"] = allt.groupby("class").size()
print("\nshare of baseline lesions carrying a SEL candidate, by mask-based class at last follow-up:")
print(ag.to_string())
