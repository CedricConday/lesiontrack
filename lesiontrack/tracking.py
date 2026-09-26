"""Lesion identity across timepoints from binary masks in a common space.

Baseline and follow-up components are matched by voxel overlap. Because
lesions merge and split, matching is done on the connected components of the
bipartite overlap graph: each group holds the baseline and follow-up lesions
that share any voxel, and the group's volume change is classified.

Classes (rates are percent of baseline volume per year):
  new        follow-up lesion with no overlap with the baseline mask dilated by
             ``new_lesion_margin_voxels``
  resolved   baseline lesion with no follow-up overlap
  enlarging  >= +change_pct_per_year and at least min_change_voxels gained
  shrinking  <= -change_pct_per_year and at least min_change_voxels lost
  stable     within +/- stable_pct_per_year
  trend_up / trend_down   between the stable and change thresholds
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import ndimage as ndi

from .config import TrackParams


def label_lesions(mask: np.ndarray, connectivity: int = 2) -> np.ndarray:
    labels, _ = ndi.label(np.asarray(mask) > 0, structure=ndi.generate_binary_structure(3, connectivity))
    return labels.astype(np.int32)


def _overlap_pairs(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    both = (a > 0) & (b > 0)
    if not both.any():
        return np.zeros((0, 3), dtype=np.int64)
    pairs = np.stack([a[both], b[both]], axis=1)
    uniq, counts = np.unique(pairs, axis=0, return_counts=True)
    return np.column_stack([uniq, counts])


class _UnionFind:
    def __init__(self) -> None:
        self.parent: dict = {}

    def find(self, x):
        self.parent.setdefault(x, x)
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a, b) -> None:
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[rb] = ra


def classify_rate(pct_per_year: float, delta_voxels: int, params: TrackParams) -> str:
    if pct_per_year >= params.change_pct_per_year and delta_voxels >= params.min_change_voxels:
        return "enlarging"
    if pct_per_year <= -params.change_pct_per_year and -delta_voxels >= params.min_change_voxels:
        return "shrinking"
    if abs(pct_per_year) <= params.stable_pct_per_year:
        return "stable"
    return "trend_up" if pct_per_year > 0 else "trend_down"


def track_pair(mask_baseline: np.ndarray, mask_followup: np.ndarray, dt_years: float,
               voxel_mm3: float, params: TrackParams | None = None) -> tuple[pd.DataFrame, dict]:
    """Match lesions between two masks. Returns (groups table, label maps)."""
    params = params or TrackParams()
    if dt_years <= 0:
        raise ValueError("dt_years must be positive")
    lb = label_lesions(mask_baseline, params.connectivity)
    lf = label_lesions(mask_followup, params.connectivity)
    st = ndi.generate_binary_structure(3, params.connectivity)
    bl_dilated = ndi.binary_dilation(lb > 0, structure=st, iterations=params.new_lesion_margin_voxels)

    vol_b = np.bincount(lb.ravel())
    vol_f = np.bincount(lf.ravel())
    pairs = _overlap_pairs(lb, lf)

    uf = _UnionFind()
    for i, j, _ in pairs:
        uf.union(("b", int(i)), ("f", int(j)))
    groups: dict = {}
    for i in range(1, vol_b.size):
        groups.setdefault(uf.find(("b", i)), {"b": [], "f": []})["b"].append(i)
    for j in range(1, vol_f.size):
        groups.setdefault(uf.find(("f", j)), {"b": [], "f": []})["f"].append(j)

    rows = []
    cent_b = dict(zip(range(1, vol_b.size), ndi.center_of_mass(lb > 0, lb, range(1, vol_b.size)))) if vol_b.size > 1 else {}
    cent_f = dict(zip(range(1, vol_f.size), ndi.center_of_mass(lf > 0, lf, range(1, vol_f.size)))) if vol_f.size > 1 else {}
    for gid, g in enumerate(groups.values(), start=1):
        vb = int(sum(vol_b[i] for i in g["b"]))
        vf = int(sum(vol_f[j] for j in g["f"]))
        if not g["b"]:
            # follow-up only: new if it clears the dilated baseline, else a fragment
            touches = any(bl_dilated[lf == j].any() for j in g["f"])
            cls = "adjacent_fragment" if touches else "new"
            rate = float("nan")
        elif not g["f"]:
            cls, rate = "resolved", -100.0 / dt_years
        else:
            rate = (vf - vb) / vb * 100.0 / dt_years
            cls = classify_rate(rate, vf - vb, params)
        cen = cent_b[g["b"][0]] if g["b"] else cent_f[g["f"][0]]
        rows.append({
            "group_id": gid,
            "baseline_ids": ",".join(map(str, g["b"])),
            "followup_ids": ",".join(map(str, g["f"])),
            "n_baseline": len(g["b"]), "n_followup": len(g["f"]),
            "merged": len(g["b"]) > 1, "split": len(g["f"]) > 1,
            "volume_baseline_mm3": vb * voxel_mm3, "volume_followup_mm3": vf * voxel_mm3,
            "pct_per_year": rate, "class": cls,
            "cx": float(cen[0]), "cy": float(cen[1]), "cz": float(cen[2]),
        })
    table = pd.DataFrame(rows)
    new_or_enlarging = (lf > 0) & ~bl_dilated
    return table, {"baseline": lb, "followup": lf, "new_or_enlarging_voxels": new_or_enlarging}


def summarize(table: pd.DataFrame, voxel_mm3: float) -> dict:
    counts = table["class"].value_counts().to_dict() if len(table) else {}
    return {
        "n_groups": len(table),
        "classes": {k: int(v) for k, v in counts.items()},
        "volume_baseline_mm3": float(table["volume_baseline_mm3"].sum()) if len(table) else 0.0,
        "volume_followup_mm3": float(table["volume_followup_mm3"].sum()) if len(table) else 0.0,
    }
