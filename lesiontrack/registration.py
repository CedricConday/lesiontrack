"""Halfway-space rigid alignment, joint T1+FLAIR deformable registration, Jacobian.

Follows the Jacobian pipeline of Elliott et al. 2019 (after Nakamura et al.):
linear alignment of the two timepoints in an unbiased halfway space, non-linear
registration driven by T1 and T2/FLAIR together, then the determinant of the
Jacobian of the deformation field as local percent volume change.

The registration engine is greedy (Yushkevich), through the ``picsl_greedy``
wheel, so no FreeSurfer or ANTs installation is needed and the pipeline runs on
x86-64 and arm64 CPUs.

Convention: with the baseline as the fixed image and the follow-up as the
moving image, greedy's warp maps baseline-space points to follow-up-space
points. Its Jacobian determinant is greater than 1 where a baseline structure
occupies a larger region at follow-up, i.e. where tissue expanded. The
synthetic backtest (:mod:`lesiontrack.synth`) checks this sign empirically.
"""

from __future__ import annotations

import io
from dataclasses import dataclass
from pathlib import Path

import nibabel as nib
import numpy as np
from scipy.linalg import expm, logm

from .config import RegParams


def run_greedy(command: str, log: Path | None = None) -> str:
    """Run one greedy command line in-process and return its console output.

    greedy tokenises the command on whitespace, so a path containing a space
    cannot be passed; that is caught here rather than surfacing as a
    misleading "file does not exist".
    """
    from picsl_greedy import Greedy3D

    buf = io.StringIO()
    Greedy3D().execute(command, out=buf, err=buf)
    text = buf.getvalue()
    if log is not None:
        with open(log, "a") as fh:
            fh.write(f"$ greedy {command}\n{text}\n")
    return text


def read_matrix(path: Path) -> np.ndarray:
    return np.loadtxt(path).reshape(4, 4)


def write_matrix(path: Path, mat: np.ndarray) -> None:
    np.savetxt(path, mat, fmt="%.10f")


def half_transform(mat: np.ndarray) -> np.ndarray:
    """Real matrix square root of a rigid/affine 4x4, so that H @ H == mat."""
    half = expm(0.5 * logm(mat))
    if np.abs(half.imag).max() > 1e-8:
        raise ValueError("transform has no real square root; is it a reflection?")
    return np.real(half)


@dataclass
class Timepoint:
    name: str
    t1: Path
    flair: Path
    mask: Path
    time_years: float
    brainmask: Path | None = None


@dataclass
class PairResult:
    """Everything produced for one baseline/follow-up pair, all in halfway space."""

    follow_up: str
    dt_years: float
    baseline_t1: Path
    baseline_flair: Path
    baseline_mask: Path
    baseline_brainmask: Path | None
    followup_t1: Path
    followup_flair: Path
    followup_mask: Path
    followup_t1_warped: Path
    warp: Path
    jacobian: Path
    halfway_matrix: Path


def _check_paths(*paths: Path | None) -> None:
    """greedy tokenises its command line on whitespace; refuse paths with spaces early."""
    bad = [str(p) for p in paths if p is not None and " " in str(p)]
    if bad:
        raise ValueError(
            "greedy cannot take paths containing spaces; symlink or move these first: " + ", ".join(bad)
        )


def _reslice(ref: Path, moving: Path, out: Path, transform: str, interp: str, threads: int,
             log: Path | None) -> None:
    run_greedy(
        f"-d 3 -threads {threads} -rf {ref} -rm {moving} {out} -ri {interp} -r {transform}", log
    )


def register_pair(baseline: Timepoint, follow: Timepoint, out_dir: Path,
                  params: RegParams | None = None, log: Path | None = None) -> PairResult:
    """Align a follow-up to the baseline in halfway space and compute the Jacobian map."""
    params = params or RegParams()
    _check_paths(out_dir, baseline.t1, baseline.flair, baseline.mask, baseline.brainmask,
                 follow.t1, follow.flair, follow.mask)
    out_dir.mkdir(parents=True, exist_ok=True)
    tag = follow.name
    thr = params.threads
    full = out_dir / f"{tag}_to_baseline_rigid.mat"
    half = out_dir / f"{tag}_halfway.mat"

    # 1. Rigid follow-up -> baseline on T1, then the square root of that transform.
    run_greedy(
        f"-d 3 -threads {thr} -a -dof {params.affine_dof} -m {params.affine_metric} "
        f"-n {params.affine_iterations} -ia-image-centers "
        f"-i {baseline.t1} {follow.t1} -o {full}",
        log,
    )
    write_matrix(half, half_transform(read_matrix(full)))

    # 2. Both timepoints into halfway space on the baseline grid.
    #    baseline -> halfway is H^-1, follow-up -> halfway is H.
    def name(tp: str, what: str) -> Path:
        return out_dir / f"{tag}_{tp}_{what}_halfway.nii.gz"

    ref = baseline.t1
    bl = {}
    fu = {}
    for what, interp in (("T1w", "LINEAR"), ("FLAIR", "LINEAR"), ("mask", "NN")):
        src_bl = getattr(baseline, {"T1w": "t1", "FLAIR": "flair", "mask": "mask"}[what])
        src_fu = getattr(follow, {"T1w": "t1", "FLAIR": "flair", "mask": "mask"}[what])
        bl[what] = name("baseline", what)
        fu[what] = name("followup", what)
        _reslice(ref, src_bl, bl[what], f"{half},-1", interp, thr, log)
        _reslice(ref, src_fu, fu[what], f"{half}", interp, thr, log)
    bl_brain = None
    if baseline.brainmask is not None:
        bl_brain = name("baseline", "brainmask")
        _reslice(ref, baseline.brainmask, bl_brain, f"{half},-1", "NN", thr, log)

    # 3. Deformable, T1 and FLAIR jointly, baseline fixed, follow-up moving.
    warp = out_dir / f"{tag}_warp.nii.gz"
    mask_arg = f"-gm {bl_brain} " if bl_brain is not None else ""
    run_greedy(
        f"-d 3 -threads {thr} -m {params.deform_metric} -n {params.deform_iterations} "
        f"-e {params.deform_step} -s {params.deform_sigma_update} {params.deform_sigma_total} "
        f"{mask_arg}"
        f"-w {params.t1_weight} -i {bl['T1w']} {fu['T1w']} "
        f"-w {params.flair_weight} -i {bl['FLAIR']} {fu['FLAIR']} -o {warp}",
        log,
    )

    # 4. Jacobian determinant of the warp on the baseline grid, plus a warped T1 for QC.
    jac = out_dir / f"{tag}_jacobian.nii.gz"
    warped = out_dir / f"{tag}_followup_T1w_warped.nii.gz"
    run_greedy(
        f"-d 3 -threads {thr} -rf {bl['T1w']} -rm {fu['T1w']} {warped} -r {warp} -rj {jac}", log
    )

    return PairResult(
        follow_up=tag,
        dt_years=follow.time_years - baseline.time_years,
        baseline_t1=bl["T1w"], baseline_flair=bl["FLAIR"], baseline_mask=bl["mask"],
        baseline_brainmask=bl_brain,
        followup_t1=fu["T1w"], followup_flair=fu["FLAIR"], followup_mask=fu["mask"],
        followup_t1_warped=warped, warp=warp, jacobian=jac, halfway_matrix=half,
    )


def jacobian_pct_per_year(jacobian: Path, dt_years: float) -> tuple[np.ndarray, nib.Nifti1Image]:
    """Local volume change as percent per year: (det J - 1) * 100 / dt."""
    if dt_years <= 0:
        raise ValueError(f"dt_years must be positive, got {dt_years}")
    img = nib.load(jacobian)
    jac = np.asarray(img.dataobj, dtype=np.float32)
    return (jac - 1.0) * 100.0 / dt_years, img
