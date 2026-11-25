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


# Lesion masks are resliced with the rotation and translation of the linear transform only.
# Since 0.3.0 the linear step is a similarity, and nearest-neighbour resampling under a scale
# a fraction of a percent from 1 does not shrink or grow a mask evenly: it duplicates (scale < 1)
# or drops (scale > 1) one whole voxel layer every 1/|1 - scale| voxels, about one plane per axis
# across a brain. A lesion cut by such a plane gained a layer at follow-up and lost one at
# baseline, and was called enlarging (mscard calibration, 2026-11-25: 14 of 593 unchanged lesions,
# changes of +9 to +39 %/yr). Without the scale the reslice keeps every lesion's voxel count,
# so lesion volumes are each timepoint's own mm^3, as before 0.3.0. The mask sits off its image by
# |1 - scale| times the distance from the grid centre: about 0.25 mm at 100 mm for a 0.25 % scale.
MASK_RESLICE_VERSION = "rigid-part-1"  # bump when mask reslicing changes; cached pairs reslice masks again


def rigid_part(mat: np.ndarray, centre: np.ndarray) -> np.ndarray:
    """The rotation and translation of a similarity or affine 4x4, without its scale.

    The rotation is the orthogonal polar factor of the linear part. The translation is
    chosen so the result agrees with ``mat`` at ``centre`` (physical coordinates, the
    frame the matrix is written in).
    """
    u, _, vt = np.linalg.svd(mat[:3, :3])
    rot = u @ vt
    if np.linalg.det(rot) < 0:
        raise ValueError("transform contains a reflection")
    out = np.eye(4)
    out[:3, :3] = rot
    out[:3, 3] = mat[:3, :3] @ centre + mat[:3, 3] - rot @ centre
    return out


def grid_centre(path: Path) -> np.ndarray:
    """Physical (RAS, mm) centre of an image grid, greedy's matrix frame."""
    img = nib.load(path)
    ijk = (np.asarray(img.shape[:3], dtype=float) - 1.0) / 2.0
    return img.affine[:3, :3] @ ijk + img.affine[:3, 3]


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
    baseline_mask_matrix: Path | None = None  # reslices a baseline label image like the baseline mask


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
        baseline_mask_matrix=out_dir / f"{tag}_baseline_mask.mat",
    )


def _reslice_masks(baseline: Timepoint, follow: Timepoint, out_dir: Path, tag: str,
                   frame: Path | None, threads: int, log: Path | None) -> tuple[Path, Path]:
    """Both lesion masks into the frame without the linear scale (see MASK_RESLICE_VERSION)."""
    half = out_dir / f"{tag}_halfway.mat"
    frame_mat = half if frame is None else frame
    fu_to_frame = half if frame is None else out_dir / f"{tag}_to_frame.mat"
    centre = grid_centre(baseline.t1)
    bl_mat, fu_mat = out_dir / f"{tag}_baseline_mask.mat", out_dir / f"{tag}_followup_mask.mat"
    write_matrix(bl_mat, rigid_part(np.linalg.inv(read_matrix(frame_mat)), centre))
    write_matrix(fu_mat, rigid_part(read_matrix(fu_to_frame), centre))
    bl_out = out_dir / f"{tag}_baseline_mask_halfway.nii.gz"
    fu_out = out_dir / f"{tag}_followup_mask_halfway.nii.gz"
    _reslice(baseline.t1, baseline.mask, bl_out, str(bl_mat), "NN", threads, log)
    _reslice(baseline.t1, follow.mask, fu_out, str(fu_mat), "NN", threads, log)
    return bl_out, fu_out


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


def rigid_qc(mat: np.ndarray, max_rotation_deg: float = 20.0, max_translation_mm: float = 40.0,
             max_scale_change_pct: float = 5.0) -> dict:
    """Rotation, translation and isotropic scale of a rigid or similarity 4x4.

    Raises if implausible for one subject. A similarity's scale is divided out of the
    rotation first; for a rigid matrix it is 1.
    """
    lin = mat[:3, :3]
    scale = float(np.cbrt(abs(np.linalg.det(lin))))
    rot = lin / scale
    angle = float(np.degrees(np.arccos(np.clip((np.trace(rot) - 1.0) / 2.0, -1.0, 1.0))))
    trans = float(np.linalg.norm(mat[:3, 3]))
    scale_pct = (scale - 1.0) * 100.0
    if angle > max_rotation_deg or trans > max_translation_mm or abs(scale_pct) > max_scale_change_pct:
        raise RegistrationError(
            f"linear registration implausible for a same-subject pair: rotation {angle:.1f} deg, "
            f"translation {trans:.1f} mm, scale {scale_pct:+.2f} % (limits {max_rotation_deg} deg, "
            f"{max_translation_mm} mm, {max_scale_change_pct} %)"
        )
    return {"rotation_deg": angle, "translation_mm": trans, "scale_pct": scale_pct}


FRAME_VERSION = "rigid-Hinv-1"  # bump when the frame composition changes; cached pairs are re-registered


def followup_to_frame(rigid: np.ndarray, frame: np.ndarray) -> np.ndarray:
    """Reslicing matrix for a follow-up into another pair's halfway frame H: rigid @ H^-1."""
    return rigid @ np.linalg.inv(frame)


def _params_record(params: RegParams, frame: Path | None) -> dict:
    return {"reg": dict(params.__dict__), "frame": str(frame) if frame else None, "jacobian": JACOBIAN_VERSION,
            "frame_version": FRAME_VERSION, "mask_reslice": MASK_RESLICE_VERSION}


def _same_registration(previous: dict | None, current: dict) -> bool:
    """The cached pair was produced by the same registration settings and frame.

    Records written before the second engine carry no ``engine`` (they were greedy)
    and no ``ants_*`` settings; those settings do not touch a greedy run, so they
    are ignored whenever the engine is greedy.
    """
    if previous is None:
        return False

    def comparable(rec: dict) -> dict:
        rec = json.loads(json.dumps(rec))  # tuples become lists, as they do in the stored record
        reg = dict(rec.get("reg", {}))
        reg.setdefault("engine", "greedy")
        if reg["engine"] == "greedy":
            reg = {k: v for k, v in reg.items() if not k.startswith("ants_")}
        else:  # greedy's saved-warp rounding never touches an ANTs warp
            reg.pop("warp_precision_vox", None)
        # Records without a frame_version predate the 0.3.0 composition fix; for a pair
        # registered without a frame the composition was never used, so they still match.
        fv = rec.get("frame_version") if rec.get("frame") else None
        return {"reg": reg, "frame": rec.get("frame"), "frame_version": fv}

    return comparable(previous) == comparable(current)


def register_pair(baseline: Timepoint, follow: Timepoint, out_dir: Path,
                  params: RegParams | None = None, log: Path | None = None,
                  force: bool = False, frame: Path | None = None) -> PairResult:
    """Align a follow-up to the baseline in halfway space and compute the Jacobian map.

    ``frame`` is an optional halfway matrix H from another pair of the same
    subject. When given, the baseline is resliced with that H^-1 and this
    follow-up with (its own rigid to baseline) @ H^-1, so every pair of the
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
            stale = False
            if previous.get("jacobian") != JACOBIAN_VERSION or not done.exists():
                # Same warp, older Jacobian (greedy's -rj, or an earlier version of ours): cheap to redo.
                jacobian_determinant(out_dir / f"{tag}_warp.nii.gz", done)
                stale = True
            if previous.get("mask_reslice") != MASK_RESLICE_VERSION:
                # Same warp, masks resliced the pre-0.3.2 way (with the similarity's scale): cheap to redo.
                _reslice_masks(baseline, follow, out_dir, tag, frame, params.threads, log)
                stale = True
            if stale:
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
    #    baseline -> frame is H^-1; follow-up -> frame is rigid @ H^-1 (== H for its own H).
    #    greedy reslices output(x) = moving(M x), M mapping fixed to moving, so a frame point x
    #    is baseline point H^-1 x and follow-up point rigid (H^-1 x). Until 0.3.0 this was
    #    H^-1 @ rigid, which equals the right matrix only when rigid and H commute: for the
    #    last pair, not for earlier follow-ups (sub-millimetre misalignment on those).
    if frame is None:
        frame_mat = half
        fu_to_frame = half
    else:
        frame_mat = frame
        fu_to_frame = out_dir / f"{tag}_to_frame.mat"
        write_matrix(fu_to_frame, followup_to_frame(rigid, read_matrix(frame)))

    def name(tp: str, what: str) -> Path:
        return out_dir / f"{tag}_{tp}_{what}_halfway.nii.gz"

    ref = baseline.t1
    bl = {}
    fu = {}
    for what in ("T1w", "FLAIR"):
        src_bl = getattr(baseline, {"T1w": "t1", "FLAIR": "flair"}[what])
        src_fu = getattr(follow, {"T1w": "t1", "FLAIR": "flair"}[what])
        bl[what] = name("baseline", what)
        fu[what] = name("followup", what)
        _reslice(ref, src_bl, bl[what], f"{frame_mat},-1", "LINEAR", thr, log)
        _reslice(ref, src_fu, fu[what], f"{fu_to_frame}", "LINEAR", thr, log)
    bl["mask"], fu["mask"] = _reslice_masks(baseline, follow, out_dir, tag, frame, thr, log)
    bl_brain = None
    if baseline.brainmask is not None:
        bl_brain = name("baseline", "brainmask")
        _reslice(ref, baseline.brainmask, bl_brain, f"{frame_mat},-1", "NN", thr, log)

    # 3. Deformable, T1 and FLAIR jointly, baseline fixed, follow-up moving.
    # 4. Jacobian determinant of the warp on the baseline grid (ours, see jacobian_determinant),
    #    plus a warped T1 for QC.
    warp = out_dir / f"{tag}_warp.nii.gz"
    jac = done
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
            f"-w {params.flair_weight} -i {bl['FLAIR']} {fu['FLAIR']} -wp {params.warp_precision_vox} -o {warp}",
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
        baseline_mask_matrix=out_dir / f"{tag}_baseline_mask.mat",
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
    # The last metric parameter is the CC neighbourhood radius (antspyx's default of 32 is a
    # 65-voxel window and about twenty times slower); the same value goes to both terms.
    extras = [(params.ants_metric, ants.image_read(str(bl_flair)), ants.image_read(str(fu_flair)),
               float(params.flair_weight), int(params.ants_radius))]
    kwargs = dict(
        type_of_transform=params.ants_transform, syn_metric=params.ants_metric,
        syn_sampling=int(params.ants_radius),
        # The images are already rigidly aligned in the halfway frame: no centre-of-mass
        # initialisation, so the stored field is the whole baseline -> follow-up map.
        initial_transform="identity",
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
    for t in result["fwdtransforms"]:
        if t.endswith(".mat"):
            # antsRegistration writes the initial transform as a stage; with the identity
            # initialisation it must be the identity, or the Jacobian would miss an affine.
            mat = ants.read_transform(t)
            if not np.allclose(mat.parameters, [1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0], atol=1e-6):
                raise RegistrationError(
                    f"ANTs transform {params.ants_transform!r} produced an affine stage ({t}); only the "
                    "displacement field reaches the Jacobian, so use a deformable-only transform such as SyNOnly")
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


def reference_mask(brain: np.ndarray, lesions: np.ndarray, lesion_margin_voxels: int = 3,
                   brain_erosion_voxels: int = 2) -> np.ndarray:
    """Brain tissue away from the brain edge and from every baseline lesion."""
    from scipy import ndimage as ndi

    brain = np.asarray(brain) > 0
    lesions = np.asarray(lesions) > 0
    if brain_erosion_voxels > 0:
        brain = ndi.binary_erosion(brain, iterations=brain_erosion_voxels)
    if lesion_margin_voxels > 0 and lesions.any():
        lesions = ndi.binary_dilation(lesions, iterations=lesion_margin_voxels)
    return brain & ~lesions


def global_scale(jacobian: np.ndarray, reference: np.ndarray, min_voxels: int = 1000,
                 trim: float = 0.10) -> float:
    """Trimmed mean Jacobian determinant over the reference tissue: the residual global volume ratio.

    Reported as a quality check, not divided out: dividing a single factor out of the map
    does not undo a global scale change, because the deformable registration puts it
    unevenly into the image (P1 backtest, 1 % scale: lesion-level false positive rate 0.23
    rigid, 0.17 divided, 0.00 with the similarity step that ``RegParams.affine_dof = 7``
    uses). Far from 1 means the linear step left global change behind. Trimmed (10 % off
    each tail) so that local change and edge error do not move it. Not the median:
    a smooth deformation leaves a large share of voxels at exactly det J = 1 (23 % of
    the brain on the MSLesSeg backtests), and the median lands on that spike whatever
    the global change.
    """
    vals = np.asarray(jacobian, dtype=np.float64)[np.asarray(reference, dtype=bool)]
    vals = vals[np.isfinite(vals) & (vals > 0)]
    if vals.size < min_voxels:
        raise ValueError(f"reference region has {vals.size} usable voxels, need {min_voxels}")
    vals.sort()
    k = int(vals.size * trim)
    return float(vals[k:vals.size - k].mean())
