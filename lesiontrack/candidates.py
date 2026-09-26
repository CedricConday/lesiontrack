"""SEL candidates: contiguous expanding regions inside pre-existing T2 lesions.

Elliott et al. 2019, "SEL candidates":
1. voxels inside baseline lesions with expansion >= JE1 (12.5 %/year);
2. 18-connected components of those voxels seed the candidates;
3. each seed grows into neighbouring voxels with expansion >= JE2 (4 %/year),
   competing so that distinct expansions stay distinct even when connected;
4. candidates under 10 voxels are discarded.

Step 3 is implemented as a marker-based watershed on the negated expansion
map, restricted to the JE2 region: each JE2 voxel joins the seed it is
reached from first along the steepest expansion path.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import ndimage as ndi
from skimage.segmentation import watershed

from .config import SELParams


def structure(connectivity: int) -> np.ndarray:
    return ndi.generate_binary_structure(3, connectivity)


def sel_candidates(expansion_pct_per_year: np.ndarray, lesion_mask: np.ndarray,
                   params: SELParams | None = None) -> np.ndarray:
    """Label map of SEL candidates (0 = background), consecutively numbered."""
    params = params or SELParams()
    exp = np.asarray(expansion_pct_per_year, dtype=np.float32)
    mask = np.asarray(lesion_mask) > 0
    if exp.shape != mask.shape:
        raise ValueError(f"shape mismatch {exp.shape} vs {mask.shape}")
    st = structure(params.connectivity)

    seeds = mask & (exp >= params.je1_pct_per_year)
    seed_labels, n_seeds = ndi.label(seeds, structure=st)
    if n_seeds == 0:
        return np.zeros(mask.shape, dtype=np.int32)

    grow = mask & (exp >= params.je2_pct_per_year)
    labels = watershed(-exp, markers=seed_labels, mask=grow, connectivity=st)

    # Reliability criterion: drop small candidates, then renumber consecutively.
    counts = np.bincount(labels.ravel())
    keep = np.flatnonzero(counts >= params.min_voxels)
    keep = keep[keep != 0]
    lut = np.zeros(counts.size, dtype=np.int32)
    lut[keep] = np.arange(1, keep.size + 1, dtype=np.int32)
    return lut[labels]


def candidate_table(labels: np.ndarray, expansion_pct_per_year: np.ndarray,
                    lesion_labels: np.ndarray | None, voxel_mm3: float,
                    params: SELParams | None = None) -> pd.DataFrame:
    """One row per candidate: size, expansion statistics, parent lesion, centroid."""
    params = params or SELParams()
    rows = []
    n = int(labels.max())
    if n == 0:
        return pd.DataFrame(columns=[
            "candidate_id", "n_voxels", "volume_mm3", "peak_pct_per_year",
            "mean_pct_per_year", "seed_voxels", "parent_lesion_id", "cx", "cy", "cz",
        ])
    idx = np.arange(1, n + 1)
    centroids = ndi.center_of_mass(labels > 0, labels, idx)
    for k, (cx, cy, cz) in zip(idx, centroids):
        sel = labels == k
        vals = expansion_pct_per_year[sel]
        parent = 0
        if lesion_labels is not None:
            parents = lesion_labels[sel]
            parents = parents[parents > 0]
            if parents.size:
                parent = int(np.bincount(parents).argmax())
        rows.append({
            "candidate_id": int(k),
            "n_voxels": int(sel.sum()),
            "volume_mm3": float(sel.sum() * voxel_mm3),
            "peak_pct_per_year": float(vals.max()),
            "mean_pct_per_year": float(vals.mean()),
            "seed_voxels": int((vals >= params.je1_pct_per_year).sum()),
            "parent_lesion_id": parent,
            "cx": float(cx), "cy": float(cy), "cz": float(cz),
        })
    return pd.DataFrame(rows)
