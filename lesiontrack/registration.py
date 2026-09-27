"""Halfway-space rigid alignment, joint T1+FLAIR deformable registration, Jacobian.

Follows the Jacobian pipeline of Elliott et al. 2019 (after Nakamura et al.):
linear alignment of the two timepoints in an unbiased halfway space, non-linear
registration driven by T1 and T2/FLAIR together, then the determinant of the
Jacobian of the deformation field as local percent volume change.

The registration engine is greedy (Yushkevich), through the ``picsl_greedy``
wheel, so no FreeSurfer or ANTs installation is needed and the pipeline runs on
x86-64 and arm64 CPUs.

Convention: with the baseline as the fixed image and the follow-up as the
moving image, the deformable warp maps baseline-space points to follow-up-space
points. Its Jacobian determinant is greater than 1 where a baseline structure
occupies a larger region at follow-up, i.e. where tissue expanded. The
determinant is computed here, by :func:`jacobian_determinant`, from the warp
file in the image's own voxel frame; it is not taken from the engine. greedy's
``-rj`` differentiates the physical-space (LPS, mm) displacement along voxel
index axes without the image direction matrix, so it is only right for images
whose voxel axes are LPS-aligned at 1 mm (found 2026-09-27 on a phantom with an
identity affine: det 0.68 reported for a true 1.95). The synthetic backtest
(:mod:`lesiontrack.synth`) checks the sign empirically.
"""

from __future__ import annotations

import io
import json
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
    # greedy applies -ri (and -rt) to the *next* -rm pair, so they must precede it.
    dtype = "-rt uchar " if interp == "NN" else ""
    run_greedy(
        f"-d 3 -threads {threads} -rf {ref} -ri {interp} {dtype}-rm {moving} {out} -r {transform}", log
    )


def _existing_pair(baseline: Timepoint, follow: Timepoint, out_dir: Path, tag: str) -> PairResult:
    def name(tp: str, what: str) -> Path:
        return out_dir / f"{tag}_{tp}_{what}_halfway.nii.gz"

    bl_brain = name("baseline", "brainmask") if baseline.brainmask is not None else None
    return PairResult(
        follow_up=tag,
        dt_years=follow.time_years - baseline.time_years,
        baseline_t1=name("baseline", "T1w"), baseline_flair=name("baseline", "FLAIR"),
        baseline_mask=name("baseline", "mask"), baseline_brainmask=bl_brain,
        followup_t1=name("followup", "T1w"), followup_flair=name("followup", "FLAIR"),
        followup_mask=name("followup", "mask"),
        followup_t1_warped=out_dir / f"{tag}_followup_T1w_warped.nii.gz",
        warp=out_dir / f"{tag}_warp.nii.gz", jacobian=out_dir / f"{tag}_jacobian.nii.gz",
        halfway_matrix=out_dir / f"{tag}_halfway.mat",
    )


class RegistrationError(RuntimeError):
    """An input or a transform failed a sanity check; results would be meaningless."""


def check_image_content(path: Path, min_nonzero_fraction: float = 0.02) -> float:
    """Refuse an image that is (nearly) empty; returns its nonzero fraction."""
    img = nib.load(path)
    data = np.asarray(img.dataobj)
    frac = float(np.count_nonzero(data) / data.size)
    if not np.isfinite(data).all():
        raise RegistrationError(f"{path}: contains NaN or inf")
    if frac < min_nonzero_fraction:
        raise RegistrationError(
            f"{path}: only {frac:.1%} of voxels are nonzero; the image is empty or broken"
        )
    return frac


def rigid_qc(mat: np.ndarray, max_rotation_deg: float = 20.0, max_translation_mm: float = 40.0) -> dict:
    """Rotation angle and translation of a rigid 4x4; raise if implausible for one subject."""
    rot = mat[:3, :3]
    angle = float(np.degrees(np.arccos(np.clip((np.trace(rot) - 1.0) / 2.0, -1.0, 1.0))))
    trans = float(np.linalg.norm(mat[:3, 3]))
    if angle > max_rotation_deg or trans > max_translation_mm:
        raise RegistrationError(
            f"rigid registration implausible for a same-subject pair: rotation {angle:.1f} deg, "
            f"translation {trans:.1f} mm (limits {max_rotation_deg} deg, {max_translation_mm} mm)"
        )
    return {"rotation_deg": angle, "translation_mm": trans}


def _params_record(params: RegParams, frame: Path | None) -> dict:
    return {"reg": dict(params.__dict__), "frame": str(frame) if frame else None, "jacobian": JACOBIAN_VERSION}


def _same_registration(previous: dict | None, current: dict) -> bool:
    """The cached pair was produced by the same registration settings and frame.

    Records written before the second engine carry no ``engine`` (they were greedy)
    and no ``ants_*`` settings; those settings do not touch a greedy run, so they
    are ignored whenever the engine is greedy.
    """
    if previous is None:
        return False

    def comparable(rec: dict) -> dict:
        reg = dict(rec.get("reg", {}))
        reg.setdefault("engine", "greedy")
        if reg["engine"] == "greedy":
            reg = {k: v for k, v in reg.items() if not k.startswith("ants_")}
        return {"reg": reg, "frame": rec.get("frame")}

    return comparable(previous) == comparable(current)


def register_pair(baseline: Timepoint, follow: Timepoint, out_dir: Path,
                  params: RegParams | None = None, log: Path | None = None,
                  force: bool = False, frame: Path | None = None) -> PairResult:
    """Align a follow-up to the baseline in halfway space and compute the Jacobian map.

    ``frame`` is an optional halfway matrix H from another pair of the same
    subject. When given, the baseline is resliced with that H^-1 and this
    follow-up with H^-1 @ (its own rigid to baseline), so every pair of the
    subject lives on one common grid. The pipeline uses the last follow-up's
    halfway matrix as the frame, which keeps that pair unbiased and lets
    constancy and tracking compare voxels across follow-ups directly.
    """
    params = params or RegParams()
    _check_paths(out_dir, baseline.t1, baseline.flair, baseline.mask, baseline.brainmask,
                 follow.t1, follow.flair, follow.mask)
    out_dir.mkdir(parents=True, exist_ok=True)
    tag = follow.name
    if follow.time_years <= baseline.time_years:
        raise ValueError(f"{tag}: time_years {follow.time_years} is not after baseline {baseline.time_years}")
    done = out_dir / f"{tag}_jacobian.nii.gz"
    record = out_dir / f"{tag}_reg_params.json"
    if (out_dir / f"{tag}_warp.nii.gz").exists() and record.exists() and not force:
        # Registration is the expensive step; reuse it so candidate rules can be re-scored,
        # but only when it was produced with the same settings and frame.
        previous = json.loads(record.read_text()) if record.exists() else None
        current = _params_record(params, frame)
        if _same_registration(previous, current):
            if previous.get("jacobian") != JACOBIAN_VERSION:
                # Same warp, older Jacobian (greedy's -rj, or an earlier version of ours): cheap to redo.
                jacobian_determinant(out_dir / f"{tag}_warp.nii.gz", done)
                record.write_text(json.dumps({**current, "rigid_qc": previous.get("rigid_qc")}, indent=2))
            return _existing_pair(baseline, follow, out_dir, tag)
    thr = params.threads
    full = out_dir / f"{tag}_to_baseline_rigid.mat"
    half = out_dir / f"{tag}_halfway.mat"
    for p in (baseline.t1, baseline.flair, follow.t1, follow.flair):
        check_image_content(p)

    # 1. Rigid follow-up -> baseline on T1, then the square root of that transform.
    run_greedy(
        f"-d 3 -threads {thr} -a -dof {params.affine_dof} -m {params.affine_metric} "
        f"-n {params.affine_iterations} -ia-image-centers "
        f"-i {baseline.t1} {follow.t1} -o {full}",
        log,
    )
    rigid = read_matrix(full)
    qc = rigid_qc(rigid)
    write_matrix(half, half_transform(rigid))

    # 2. Both timepoints into the common frame on the baseline grid.
    #    baseline -> frame is H^-1; follow-up -> frame is H^-1 @ rigid (== H for its own H).
    if frame is None:
        frame_mat = half
        fu_to_frame = half
    else:
        frame_mat = frame
        fu_to_frame = out_dir / f"{tag}_to_frame.mat"
        write_matrix(fu_to_frame, np.linalg.inv(read_matrix(frame)) @ rigid)

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
        _reslice(ref, src_bl, bl[what], f"{frame_mat},-1", interp, thr, log)
        _reslice(ref, src_fu, fu[what], f"{fu_to_frame}", interp, thr, log)
    bl_brain = None
    if baseline.brainmask is not None:
        bl_brain = name("baseline", "brainmask")
        _reslice(ref, baseline.brainmask, bl_brain, f"{frame_mat},-1", "NN", thr, log)

    # 3. Deformable, T1 and FLAIR jointly, baseline fixed, follow-up moving.
    # 4. Jacobian determinant of the warp on the baseline grid (ours, see jacobian_determinant),
    #    plus a warped T1 for QC.
    warp = out_dir / f"{tag}_warp.nii.gz"
    jac = out_dir / f"{tag}_jacobian.nii.gz"
    warped = out_dir / f"{tag}_followup_T1w_warped.nii.gz"
    if params.engine == "ants":
        ants_deform(bl["T1w"], fu["T1w"], bl["FLAIR"], fu["FLAIR"], bl_brain, warp, jac, warped, params, log)
    elif params.engine == "greedy":
        mask_arg = f"-gm {bl_brain} " if bl_brain is not None else ""
        sv_arg = "-sv " if params.stationary_velocity else ""
        run_greedy(
            f"-d 3 -threads {thr} -m {params.deform_metric} -n {params.deform_iterations} "
            f"-e {params.deform_step} -s {params.deform_sigma_update} {params.deform_sigma_total} "
            f"{sv_arg}{mask_arg}"
            f"-w {params.t1_weight} -i {bl['T1w']} {fu['T1w']} "
            f"-w {params.flair_weight} -i {bl['FLAIR']} {fu['FLAIR']} -o {warp}",
            log,
        )
        run_greedy(f"-d 3 -threads {thr} -rf {bl['T1w']} -rm {fu['T1w']} {warped} -r {warp}", log)
        jacobian_determinant(warp, jac)
    else:
        raise ValueError(f"unknown registration engine {params.engine!r}; use 'greedy' or 'ants'")
    record.write_text(json.dumps({**_params_record(params, frame), "rigid_qc": qc}, indent=2))

    return PairResult(
        follow_up=tag,
        dt_years=follow.time_years - baseline.time_years,
        baseline_t1=bl["T1w"], baseline_flair=bl["FLAIR"], baseline_mask=bl["mask"],
        baseline_brainmask=bl_brain,
        followup_t1=fu["T1w"], followup_flair=fu["FLAIR"], followup_mask=fu["mask"],
        followup_t1_warped=warped, warp=warp, jacobian=jac, halfway_matrix=half,
    )


def ants_deform(bl_t1: Path, fu_t1: Path, bl_flair: Path, fu_flair: Path, bl_brain: Path | None,
                warp: Path, jac: Path, warped: Path, params: RegParams, log: Path | None = None) -> None:
    """Deformable step with ANTs SyN (antspyx): baseline fixed, follow-up moving, T1 and FLAIR jointly.

    Writes the displacement field, its Jacobian determinant on the baseline grid, and the
    warped follow-up T1, with the same names greedy's path produces, so everything
    downstream is engine-blind. The field is the map from baseline to follow-up space
    and the determinant is computed by :func:`jacobian_determinant`, exactly as for
    greedy: det > 1 where the follow-up is larger than the baseline.
    """
    import os
    import shutil

    try:
        import ants
    except ImportError as e:  # pragma: no cover
        raise RuntimeError("engine 'ants' needs antspyx: pip install 'lesiontrack[ants]'") from e
    os.environ.setdefault("ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS", str(params.threads))
    fixed, moving = ants.image_read(str(bl_t1)), ants.image_read(str(fu_t1))
    extras = [(params.ants_metric, ants.image_read(str(bl_flair)), ants.image_read(str(fu_flair)),
               float(params.flair_weight), 4)]
    kwargs = dict(
        type_of_transform=params.ants_transform, syn_metric=params.ants_metric,
        reg_iterations=tuple(int(i) for i in params.ants_iterations),
        grad_step=params.ants_grad_step, flow_sigma=params.ants_flow_sigma, total_sigma=params.ants_total_sigma,
        multivariate_extras=extras, outprefix=str(warp.with_suffix("").with_suffix("")) + "_ants_",
        verbose=False,
    )
    if bl_brain is not None:
        kwargs["mask"] = ants.image_read(str(bl_brain))
    result = ants.registration(fixed, moving, **kwargs)
    field = next((t for t in result["fwdtransforms"] if t.endswith(".nii.gz")), None)
    if field is None:
        raise RegistrationError("ANTs returned no displacement field; check the transform type")
    shutil.copyfile(field, warp)
    ants.image_write(result["warpedmovout"], str(warped))
    jacobian_determinant(warp, jac)
    if log is not None:
        with open(log, "a") as fh:
            fh.write(f"$ ants.registration({params.ants_transform}, metric {params.ants_metric}, iterations {params.ants_iterations})\n")


JACOBIAN_VERSION = "voxel-frame-1"  # bump when the Jacobian computation changes; cached pairs are recomputed


def jacobian_determinant(warp: Path, out: Path | None = None) -> nib.Nifti1Image:
    """det J of x -> x + u(x) for an ITK displacement field, on the field's own grid.

    ITK (greedy, ANTs) stores displacements in physical LPS millimetres. The
    determinant is invariant to the frame only if displacement and derivative use
    the same one, so the vectors are first taken to voxel units: LPS -> RAS, then
    through the inverse of the affine's direction-and-spacing block. Central
    differences along voxel axes then give det(I + grad u). The result equals the
    ratio follow-up volume / baseline volume of the tissue at each baseline voxel.
    """
    img = nib.load(warp)
    u = np.asarray(img.dataobj, dtype=np.float64)
    u = u.reshape(u.shape[:3] + (-1,))
    if u.ndim != 4 or u.shape[-1] != 3:
        raise ValueError(f"{warp}: expected a 3-D displacement field with 3 components, got shape {u.shape}")
    to_vox = np.linalg.inv(img.affine[:3, :3]) @ np.diag([-1.0, -1.0, 1.0])
    uv = np.einsum("ij,xyzj->xyzi", to_vox, u)
    grad = np.empty(u.shape[:3] + (3, 3))
    for i in range(3):
        for j in range(3):
            grad[..., i, j] = np.gradient(uv[..., i], axis=j)
    det = np.linalg.det(np.eye(3) + grad).astype(np.float32)
    result = nib.Nifti1Image(det, img.affine)
    if out is not None:
        nib.save(result, out)
    return result


def jacobian_pct_per_year(jacobian: Path, dt_years: float) -> tuple[np.ndarray, nib.Nifti1Image]:
    """Local volume change as percent per year: (det J - 1) * 100 / dt."""
    if dt_years <= 0:
        raise ValueError(f"dt_years must be positive, got {dt_years}")
    img = nib.load(jacobian)
    jac = np.asarray(img.dataobj, dtype=np.float32)
    return (jac - 1.0) * 100.0 / dt_years, img
