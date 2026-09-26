"""Concentricity and constancy of expansion, and the cohort-level SEL score.

Elliott et al. 2019 scored each candidate for concentric (inside-out) and
constant (linear-in-time) expansion, z-scored both across the cohort and kept
candidates with a combined score >= 0 as SELs. Vanden Bulcke et al. 2025
(medRxiv 2025.12.03.25341192) describe the quantities as:

* concentricity: slope of a Huber regression of mean expansion per concentric
  band against the band's distance from the candidate edge (positive = more
  expansion at the core);
* constancy: mean normalised squared residual of a Huber linear fit of the
  candidate's mean expansion over the intermediate timepoints (lower = more
  constant);
* S = z(concentricity) - z(constancy residual); definite SEL if S >= 0.

Both z-scores are relative to the cohort of candidates being scored, so the
definite/possible split is a ranking within the cohort, not an absolute
threshold. This module keeps the raw quantities so users can apply their own.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import ndimage as ndi


def huber_fit(x: np.ndarray, y: np.ndarray, delta: float = 1.345, iters: int = 50,
              through_origin: bool = False) -> tuple[float, float]:
    """Robust line fit by iteratively re-weighted least squares. Returns (slope, intercept)."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.size < 2:
        raise ValueError("need at least two points")
    w = np.ones_like(x)
    slope, intercept = 0.0, 0.0
    for _ in range(iters):
        if through_origin:
            denom = np.sum(w * x * x)
            slope = float(np.sum(w * x * y) / denom) if denom > 0 else 0.0
            intercept = 0.0
        else:
            sw = w.sum()
            mx, my = np.sum(w * x) / sw, np.sum(w * y) / sw
            denom = np.sum(w * (x - mx) ** 2)
            slope = float(np.sum(w * (x - mx) * (y - my)) / denom) if denom > 0 else 0.0
            intercept = float(my - slope * mx)
        r = y - (slope * x + intercept)
        scale = np.median(np.abs(r - np.median(r))) / 0.6745
        if scale <= 1e-12:
            break
        a = np.abs(r) / scale
        w_new = np.where(a <= delta, 1.0, delta / np.maximum(a, 1e-12))
        if np.allclose(w_new, w):
            break
        w = w_new
    return slope, intercept


def concentricity(candidate: np.ndarray, expansion: np.ndarray) -> tuple[float, int]:
    """Slope of mean expansion vs. distance-from-edge band. Returns (slope, n_bands).

    Bands are 1-voxel shells from the Euclidean distance transform inside the
    candidate. With fewer than two bands the slope is undefined (NaN).
    """
    sel = np.asarray(candidate) > 0
    if sel.sum() == 0:
        return float("nan"), 0
    dist = ndi.distance_transform_edt(sel)
    bands = np.ceil(dist[sel]).astype(int)  # 1 = outermost shell
    vals = expansion[sel]
    ids = np.unique(bands)
    if ids.size < 2:
        return float("nan"), int(ids.size)
    means = np.array([vals[bands == b].mean() for b in ids])
    slope, _ = huber_fit(ids.astype(float), means)
    return float(slope), int(ids.size)


def constancy(times_years: np.ndarray, cumulative_expansion_pct: np.ndarray) -> tuple[float, float]:
    """Mean normalised squared residual of a through-origin linear fit of expansion vs time.

    ``times_years`` are measured from baseline (baseline itself may be omitted;
    it is added as (0, 0)). Returns (residual, slope_pct_per_year); residual is
    NaN with fewer than two follow-ups, because a line through the origin and
    one point has no residual.
    """
    t = np.asarray(times_years, dtype=float)
    e = np.asarray(cumulative_expansion_pct, dtype=float)
    if t.size < 2:
        return float("nan"), float("nan")
    t0 = np.concatenate([[0.0], t])
    e0 = np.concatenate([[0.0], e])
    slope, _ = huber_fit(t0, e0, through_origin=True)
    fit = slope * t0
    scale = max(abs(slope * t0[-1]), 1e-6)
    resid = float(np.mean((e0 - fit) ** 2) / scale**2)
    return resid, float(slope)


def cohort_score(table: pd.DataFrame) -> pd.DataFrame:
    """Add z-scores, the combined score S and the definite/possible label.

    Expects columns ``concentricity`` and ``constancy_residual`` (NaN allowed).
    Candidates without a constancy value (two-timepoint studies) get
    S = z(concentricity) and ``score_basis`` = "concentricity_only".
    """
    out = table.copy()

    def z(col: str) -> pd.Series:
        v = out[col].astype(float)
        sd = v.std(ddof=0)
        if not np.isfinite(sd) or sd == 0:
            return pd.Series(np.where(v.notna(), 0.0, np.nan), index=out.index)
        return (v - v.mean()) / sd

    out["z_concentricity"] = z("concentricity")
    out["z_constancy"] = z("constancy_residual")
    has_const = out["z_constancy"].notna()
    out["score_S"] = np.where(
        has_const, out["z_concentricity"] - out["z_constancy"], out["z_concentricity"]
    )
    out["score_basis"] = np.where(has_const, "concentricity_and_constancy", "concentricity_only")
    out["definite_sel"] = out["score_S"] >= 0
    return out
