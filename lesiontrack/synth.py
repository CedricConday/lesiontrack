"""Synthetic follow-up scans with known lesion expansion: the backtest gate.

No public dataset carries ground truth for slowly expanding lesions, and the
field's thresholds were set by eye (Elliott 2019). Before a real result is
reported, the pipeline must recover expansions it is known to contain.

Given a real baseline (T1, FLAIR, lesion mask), :func:`make_followup` picks
lesions, expands each chosen one radially by a known volume factor with a
smooth deformation, applies a small rigid perturbation, rescales intensities
and adds noise, and writes the synthetic follow-up together with the truth:
which lesion expanded by how much and the expanded mask. :func:`score`
compares a pipeline run against that truth.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd
from scipy import ndimage as ndi
from scipy.spatial.transform import Rotation

from .tracking import label_lesions


@dataclass
class SynthLesion:
    lesion_id: int
    n_voxels: int
    volume_factor: float  # nominal follow-up volume / baseline volume; 1.0 = untouched
    cx: float
    cy: float
    cz: float
    true_mean_expansion_pct: float = float("nan")  # mean analytic expansion inside the lesion


def _radial_displacement(shape: tuple, centre: np.ndarray, radius_vox: float,
                         volume_factor: float, falloff_vox: float) -> np.ndarray:
    """Displacement (3, X, Y, Z) that pulls follow-up voxels from baseline positions.

    Inside ``radius_vox`` the field is a pure radial scaling by
    ``volume_factor ** (1/3)``; it fades to zero over ``falloff_vox`` beyond
    that, so the surrounding brain is untouched.
    """
    lam = volume_factor ** (1.0 / 3.0)
    grids = np.meshgrid(*[np.arange(n, dtype=np.float32) for n in shape], indexing="ij")
    d = np.stack([g - c for g, c in zip(grids, centre)])
    r = np.sqrt((d**2).sum(axis=0))
    # follow-up(x) = baseline(x + u(x)), u = (x - c) * (1/lam - 1) * w(r)
    w = np.clip(1.0 - (r - radius_vox) / falloff_vox, 0.0, 1.0)
    w = np.where(r <= radius_vox, 1.0, w)
    return d * ((1.0 / lam - 1.0) * w)[None]


def _rigid_field(shape: tuple, rot_deg: np.ndarray, trans_vox: np.ndarray) -> np.ndarray:
    grids = np.meshgrid(*[np.arange(n, dtype=np.float32) for n in shape], indexing="ij")
    centre = (np.array(shape, dtype=np.float32) - 1) / 2
    x = np.stack([g - c for g, c in zip(grids, centre)])
    rot = Rotation.from_euler("xyz", rot_deg, degrees=True).as_matrix().astype(np.float32)
    xr = np.tensordot(rot, x, axes=(1, 0))
    return xr - x + trans_vox[:, None, None, None].astype(np.float32)


def _jacobian_det(displacement: np.ndarray) -> np.ndarray:
    """det(I + grad u) of a displacement field shaped (3, X, Y, Z), central differences."""
    g = np.empty((3, 3, *displacement.shape[1:]), dtype=np.float32)
    for i in range(3):
        for j in range(3):
            g[i, j] = np.gradient(displacement[i], axis=j)
        g[i, i] += 1.0
    return np.linalg.det(np.moveaxis(g, (0, 1), (-2, -1)))


def _warp(image: np.ndarray, displacement: np.ndarray, order: int, cval: float = 0.0) -> np.ndarray:
    grids = np.meshgrid(*[np.arange(n, dtype=np.float32) for n in image.shape], indexing="ij")
    coords = np.stack(grids) + displacement
    return ndi.map_coordinates(image, coords, order=order, mode="constant", cval=cval)


def make_followup(t1: Path, flair: Path, mask: Path, out_dir: Path, *,
                  n_expand: int = 8, volume_factors: tuple = (1.15, 1.25, 1.40),
                  min_voxels: int = 30, margin_vox: float = 2.0, falloff_vox: float = 6.0,
                  rigid_rot_deg: float = 1.0, rigid_trans_vox: float = 1.0,
                  intensity_scale: tuple = (0.9, 1.1), noise_frac: float = 0.02,
                  seed: int = 0, time_fraction: float = 1.0) -> dict:
    """Write a synthetic follow-up and its ground truth. Returns the truth record.

    ``time_fraction`` scales the injected expansion to an intermediate
    timepoint: a lesion with volume factor ``f`` at fraction 1 has
    ``f ** time_fraction`` here, so expansion is linear in log-volume over
    time, while lesion choice and factors stay those of ``seed``. Rigid
    perturbation, intensity scale and noise are drawn independently per
    timepoint.
    """
    rng = np.random.default_rng(seed)
    rng_tp = np.random.default_rng([seed, int(round(time_fraction * 1000))])
    out_dir.mkdir(parents=True, exist_ok=True)
    t1_img = nib.load(t1)
    t1_arr = np.asarray(t1_img.dataobj, dtype=np.float32)
    fl_arr = np.asarray(nib.load(flair).dataobj, dtype=np.float32)
    mk_arr = np.asarray(nib.load(mask).dataobj) > 0
    labels = label_lesions(mk_arr)
    sizes = np.bincount(labels.ravel())
    eligible = [k for k in range(1, sizes.size) if sizes[k] >= min_voxels]
    if len(eligible) < n_expand:
        raise ValueError(f"only {len(eligible)} lesions >= {min_voxels} voxels; asked for {n_expand}")
    chosen = sorted(rng.choice(eligible, size=n_expand, replace=False).tolist())
    factors = [float(volume_factors[i % len(volume_factors)]) for i in range(n_expand)]
    rng.shuffle(factors)
    centroids = dict(zip(range(1, sizes.size), ndi.center_of_mass(mk_arr, labels, range(1, sizes.size))))

    disp = np.zeros((3, *t1_arr.shape), dtype=np.float32)
    truth = []
    for k, f in zip(chosen, factors):
        c = np.array(centroids[k], dtype=np.float32)
        # equivalent-sphere radius plus a margin so the whole lesion scales
        radius = (3 * sizes[k] / (4 * np.pi)) ** (1 / 3) + margin_vox
        f_here = f ** time_fraction
        disp += _radial_displacement(t1_arr.shape, c, radius, f_here, falloff_vox)
        truth.append(SynthLesion(int(k), int(sizes[k]), f_here, *map(float, c)))
    for k in range(1, sizes.size):
        if k not in chosen:
            c = centroids[k]
            truth.append(SynthLesion(int(k), int(sizes[k]), 1.0, *map(float, c)))

    # Analytic truth: the sampling map is x -> x + u(x); its Jacobian determinant is the
    # follow-up-to-baseline volume ratio, so the expansion the pipeline should measure is
    # its reciprocal. Computed before the rigid part, which has unit determinant.
    truth_expansion = (1.0 / _jacobian_det(disp) - 1.0) * 100.0

    rot = rng_tp.uniform(-rigid_rot_deg, rigid_rot_deg, size=3)
    trans = rng_tp.uniform(-rigid_trans_vox, rigid_trans_vox, size=3)
    disp_total = disp + _rigid_field(t1_arr.shape, rot, trans)

    def synth(arr: np.ndarray, order: int, cval: float = 0.0) -> np.ndarray:
        return _warp(arr, disp_total, order, cval)

    t1_fu = synth(t1_arr, 1)
    fl_fu = synth(fl_arr, 1)
    # Warp a signed distance function rather than the binary mask so that
    # sub-voxel boundary shifts of small lesions survive resampling.
    sdf = ndi.distance_transform_edt(~mk_arr) - ndi.distance_transform_edt(mk_arr)
    mask_fu = synth(sdf.astype(np.float32), 1, cval=1e6) <= 0.0  # outside the FOV is background
    brain = t1_arr > 0
    brain_fu = synth(brain.astype(np.float32), 1) >= 0.5
    for arr in (t1_fu, fl_fu):
        scale = rng_tp.uniform(*intensity_scale)
        sigma = noise_frac * arr[brain_fu].mean()
        arr *= scale
        arr += rng_tp.normal(0.0, sigma, size=arr.shape).astype(np.float32)
        arr[~brain_fu] = 0.0

    aff = t1_img.affine
    paths = {
        "t1": out_dir / "synth_T1w.nii.gz",
        "flair": out_dir / "synth_FLAIR.nii.gz",
        "mask": out_dir / "synth_mask.nii.gz",
        "brainmask": out_dir / "synth_brainmask.nii.gz",
        "baseline_lesion_labels": out_dir / "baseline_lesion_labels.nii.gz",
        "truth_expansion_pct": out_dir / "truth_expansion_pct.nii.gz",
    }
    nib.save(nib.Nifti1Image(truth_expansion.astype(np.float32), aff), paths["truth_expansion_pct"])
    nib.save(nib.Nifti1Image(t1_fu.astype(np.float32), aff), paths["t1"])
    nib.save(nib.Nifti1Image(fl_fu.astype(np.float32), aff), paths["flair"])
    nib.save(nib.Nifti1Image(mask_fu.astype(np.uint8), aff), paths["mask"])
    nib.save(nib.Nifti1Image(brain_fu.astype(np.uint8), aff), paths["brainmask"])
    nib.save(nib.Nifti1Image(labels.astype(np.int32), aff), paths["baseline_lesion_labels"])

    for les in truth:
        les.true_mean_expansion_pct = float(truth_expansion[labels == les.lesion_id].mean())

    record = {
        "source": {"t1": str(t1), "flair": str(flair), "mask": str(mask)},
        "seed": seed,
        "time_fraction": time_fraction,
        "rigid_rotation_deg": rot.tolist(),
        "rigid_translation_vox": trans.tolist(),
        "lesions": [asdict(x) for x in truth],
        "outputs": {k: str(v) for k, v in paths.items()},
    }
    with open(out_dir / "truth.json", "w") as fh:
        json.dump(record, fh, indent=2)
    return record


def score(truth: dict, baseline_lesion_labels: np.ndarray, expansion_pct_per_year: np.ndarray,
          candidate_labels: np.ndarray, dt_years: float,
          definite_labels: np.ndarray | None = None) -> tuple[pd.DataFrame, dict]:
    """Per-lesion recovery of known expansion, and cohort-level gate metrics.

    A lesion counts as detected when any SEL candidate overlaps it. The
    measured rate is the mean expansion inside the lesion, compared with the
    mean of the analytic Jacobian of the injected deformation inside that
    lesion (which is below the nominal factor wherever the lesion extends past
    the radial window), divided by ``dt``.

    Gate metrics:

    * ``recovery_fraction``: robust through-origin slope of measured against
      true rate over the expanded lesions (1.0 = the pipeline measures what
      was injected; 0.3 = it reports a third of it);
    * ``noise_peak_p95_pct_per_year``: 95th percentile of the per-voxel
      expansion inside untouched lesions, the level a seed threshold has to
      clear to keep false positives down;
    * sensitivity and false positive rate at the configured SEL thresholds,
      overall and by lesion size stratum.
    """
    from .scores import huber_fit

    rows = []
    untouched_voxels = []
    for les in truth["lesions"]:
        k = les["lesion_id"]
        sel = baseline_lesion_labels == k
        if not sel.any():
            continue
        true_rate = les["true_mean_expansion_pct"] / dt_years
        vals = expansion_pct_per_year[sel]
        measured = float(vals.mean())
        detected = bool((candidate_labels[sel] > 0).any())
        definite = bool((definite_labels[sel] > 0).any()) if definite_labels is not None else detected
        expanded = les["volume_factor"] > 1.0
        if not expanded:
            untouched_voxels.append(vals)
        rows.append({
            "lesion_id": k, "n_voxels": les["n_voxels"],
            "size_stratum": _stratum(les["n_voxels"]),
            "nominal_pct_per_year": (les["volume_factor"] - 1.0) * 100.0 / dt_years,
            "true_pct_per_year": true_rate,
            "measured_mean_pct_per_year": measured,
            "measured_peak_pct_per_year": float(vals.max()),
            "recovery": measured / true_rate if expanded and true_rate else float("nan"),
            "detected": detected, "detected_definite": definite, "expanded": expanded,
        })
    table = pd.DataFrame(rows)
    exp = table[table["expanded"]]
    unexp = table[~table["expanded"]]
    if len(exp) >= 2:
        slope, _ = huber_fit(exp["true_pct_per_year"].to_numpy(), exp["measured_mean_pct_per_year"].to_numpy(), through_origin=True)
    else:
        slope = float("nan")
    un_vox = np.concatenate(untouched_voxels) if untouched_voxels else np.array([np.nan])
    metrics = {
        "n_expanded": len(exp), "n_unexpanded": len(unexp),
        "sensitivity": float(exp["detected"].mean()) if len(exp) else float("nan"),
        "false_positive_rate": float(unexp["detected"].mean()) if len(unexp) else float("nan"),
        "sensitivity_definite": float(exp["detected_definite"].mean()) if len(exp) else float("nan"),
        "false_positive_rate_definite": float(unexp["detected_definite"].mean()) if len(unexp) else float("nan"),
        "recovery_fraction": float(slope),
        "recovery_median": float(exp["recovery"].median()) if len(exp) else float("nan"),
        "noise_mean_pct_per_year": float(np.nanmean(un_vox)),
        "noise_peak_p95_pct_per_year": float(np.nanpercentile(un_vox, 95)),
        "noise_peak_max_pct_per_year": float(np.nanmax(un_vox)),
        "by_stratum": {
            s: {
                "n_expanded": int((g["expanded"]).sum()),
                "sensitivity": float(g[g["expanded"]]["detected"].mean()) if g["expanded"].any() else float("nan"),
                "n_unexpanded": int((~g["expanded"]).sum()),
                "false_positive_rate": float(g[~g["expanded"]]["detected"].mean()) if (~g["expanded"]).any() else float("nan"),
                "recovery_median": float(g[g["expanded"]]["recovery"].median()) if g["expanded"].any() else float("nan"),
            }
            for s, g in table.groupby("size_stratum")
        },
    }
    return table, metrics


def _stratum(n_voxels: int) -> str:
    if n_voxels < 100:
        return "small_lt100"
    if n_voxels < 500:
        return "medium_100_500"
    return "large_ge500"
