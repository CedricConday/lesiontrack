"""Subject-level pipeline: register, expand, detect, score, track, write."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd

from . import __version__
from .candidates import candidate_table, sel_candidates
from .config import Params
from .registration import PairResult, Timepoint, jacobian_pct_per_year, register_pair
from .scores import concentricity, constancy
from .tracking import label_lesions, summarize, track_pair


@dataclass
class SubjectResult:
    subject: str
    out_dir: Path
    pairs: list[PairResult]
    candidates: pd.DataFrame
    tracking: pd.DataFrame


def _voxel_mm3(img: nib.Nifti1Image) -> float:
    return float(np.prod(img.header.get_zooms()[:3]))


def run_subject(subject: str, timepoints: list[Timepoint], out_dir: Path,
                params: Params | None = None) -> SubjectResult:
    """Run every follow-up against the baseline (first timepoint by time)."""
    params = params or Params()
    if len(timepoints) < 2:
        raise ValueError(f"{subject}: need at least two timepoints")
    tps = sorted(timepoints, key=lambda t: t.time_years)
    baseline, follows = tps[0], tps[1:]
    out_dir.mkdir(parents=True, exist_ok=True)
    log = out_dir / "greedy.log"
    log.write_text("")

    for fu in follows:
        if fu.time_years <= baseline.time_years:
            raise ValueError(f"{subject}: {fu.name} is not after baseline {baseline.name}")
    # The last follow-up defines the subject's halfway frame; earlier follow-ups are
    # registered into that same frame so all maps share one grid.
    last = register_pair(baseline, follows[-1], out_dir / "reg", params.reg, log)
    pairs = [register_pair(baseline, fu, out_dir / "reg", params.reg, log, frame=last.halfway_matrix)
             for fu in follows[:-1]] + [last]

    # Expansion maps (percent per year) for every pair, on the halfway baseline grid.
    expansion = {}
    ref_img = None
    for p in pairs:
        exp, img = jacobian_pct_per_year(p.jacobian, p.dt_years)
        expansion[p.follow_up] = exp
        ref_img = img
        nib.save(nib.Nifti1Image(exp, img.affine), out_dir / f"expansion_{p.follow_up}_pct_per_year.nii.gz")
    voxel_mm3 = _voxel_mm3(ref_img)

    # SEL candidates from baseline to the last follow-up (Elliott 2019).
    mask_bl = np.asarray(nib.load(last.baseline_mask).dataobj) > 0
    lesion_labels = label_lesions(mask_bl, params.sel.connectivity)
    cand = sel_candidates(expansion[last.follow_up], mask_bl, params.sel)
    nib.save(nib.Nifti1Image(cand.astype(np.int32), ref_img.affine), out_dir / "sel_candidates.nii.gz")
    nib.save(nib.Nifti1Image(lesion_labels, ref_img.affine), out_dir / "baseline_lesion_labels.nii.gz")
    table = candidate_table(cand, expansion[last.follow_up], lesion_labels, voxel_mm3, params.sel)

    conc, nb, const_res, const_slope = [], [], [], []
    times = np.array([p.dt_years for p in pairs])
    for k in table["candidate_id"]:
        sel = cand == k
        s, n = concentricity(sel, expansion[last.follow_up])
        conc.append(s)
        nb.append(n)
        cum = np.array([expansion[p.follow_up][sel].mean() * p.dt_years for p in pairs])
        r, sl = constancy(times, cum)
        const_res.append(r)
        const_slope.append(sl)
    table["concentricity"] = conc
    table["n_bands"] = nb
    table["constancy_residual"] = const_res
    table["constancy_slope_pct_per_year"] = const_slope
    table.insert(0, "subject", subject)
    table["dt_years"] = last.dt_years
    table["n_followups"] = len(pairs)
    table.to_csv(out_dir / "sel_candidates.tsv", sep="\t", index=False, float_format="%.5g")

    # Lesion tracking for every pair.
    track_tables = []
    for p in pairs:
        mask_fu = np.asarray(nib.load(p.followup_mask).dataobj) > 0
        t, maps = track_pair(mask_bl, mask_fu, p.dt_years, voxel_mm3, params.track)
        t.insert(0, "follow_up", p.follow_up)
        t.insert(0, "subject", subject)
        t["dt_years"] = p.dt_years
        track_tables.append(t)
        nib.save(nib.Nifti1Image(maps["followup"], ref_img.affine), out_dir / f"lesion_labels_{p.follow_up}.nii.gz")
        nib.save(nib.Nifti1Image(maps["new_or_enlarging_voxels"].astype(np.uint8), ref_img.affine),
                 out_dir / f"new_or_enlarging_{p.follow_up}.nii.gz")
    tracking = pd.concat(track_tables, ignore_index=True)
    tracking.to_csv(out_dir / "lesion_tracking.tsv", sep="\t", index=False, float_format="%.5g")

    summary = {
        "lesiontrack_version": __version__,
        "subject": subject,
        "baseline": baseline.name,
        "follow_ups": [p.follow_up for p in pairs],
        "dt_years": [p.dt_years for p in pairs],
        "voxel_mm3": voxel_mm3,
        "baseline_lesion_count": int(lesion_labels.max()),
        "baseline_lesion_volume_mm3": float(mask_bl.sum() * voxel_mm3),
        "sel_candidates": len(table),
        "sel_candidate_volume_mm3": float(table["volume_mm3"].sum()) if len(table) else 0.0,
        "tracking": {p.follow_up: summarize(tracking[tracking["follow_up"] == p.follow_up], voxel_mm3) for p in pairs},
        "params": {
            "sel": params.sel.__dict__, "track": params.track.__dict__, "reg": params.reg.__dict__,
        },
        "registration_records": [str(out_dir / "reg" / f"{p.follow_up}_reg_params.json") for p in pairs],
    }
    with open(out_dir / "summary.json", "w") as fh:
        json.dump(summary, fh, indent=2)
    return SubjectResult(subject, out_dir, pairs, table, tracking)
