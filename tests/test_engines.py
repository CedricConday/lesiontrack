"""Both deformable engines see the same injected expansion with the same sign."""
import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from lesiontrack.config import RegParams
from lesiontrack.registration import Timepoint, jacobian_pct_per_year, register_pair
from lesiontrack.synth import make_followup


def _phantom(root: Path, shape=(56, 56, 56)):
    g = np.indices(shape).astype(float)
    c = np.array(shape) / 2
    r = np.sqrt(((g - c[:, None, None, None]) ** 2).sum(0))
    rng = np.random.default_rng(0)
    t1 = np.where(r < 24, 300.0, 0.0) + np.where(r < 20, 150.0, 0.0)  # grey shell, bright white core
    t1 += (t1 > 0) * rng.normal(0, 4, shape)
    t1 += (t1 > 0) * 60 * np.sin(g[0] / 3) * np.cos(g[1] / 4)  # texture so a metric has something to lock on
    fl = np.where(r < 24, 200.0, 0.0) - np.where(r < 20, 40.0, 0.0) + (t1 > 0) * rng.normal(0, 4, shape)
    mask = np.zeros(shape, np.uint8)
    lesion_centre = c + np.array([6, -5, 4])
    rl = np.sqrt(((g - lesion_centre[:, None, None, None]) ** 2).sum(0))
    mask[rl < 4.5] = 1
    t1[mask > 0] *= 0.8
    fl[mask > 0] *= 1.6
    aff = np.eye(4)
    anat = root / "anat"
    anat.mkdir(parents=True)
    nib.save(nib.Nifti1Image(t1.astype(np.float32), aff), anat / "T1w.nii.gz")
    nib.save(nib.Nifti1Image(fl.astype(np.float32), aff), anat / "FLAIR.nii.gz")
    nib.save(nib.Nifti1Image(mask, aff), anat / "mask.nii.gz")
    nib.save(nib.Nifti1Image((t1 > 0).astype(np.uint8), aff), anat / "brainmask.nii.gz")
    return anat, tuple(int(v) for v in lesion_centre)


@pytest.mark.parametrize("engine", ["greedy", "ants"])
def test_both_engines_recover_the_sign_of_an_injected_expansion(tmp_path, engine):
    if engine == "ants":
        pytest.importorskip("ants")
    anat, centre = _phantom(tmp_path)
    fu_dir = tmp_path / "fu"
    truth = make_followup(anat / "T1w.nii.gz", anat / "FLAIR.nii.gz", anat / "mask.nii.gz", fu_dir,
                          n_expand=1, volume_factors=(1.6,), min_voxels=20, rigid_rot_deg=0.5, rigid_trans_vox=0.5,
                          noise_frac=0.01, seed=1)
    assert truth["lesions"], "the phantom's one lesion must have been chosen for expansion"
    base = Timepoint(name="baseline", t1=anat / "T1w.nii.gz", flair=anat / "FLAIR.nii.gz", mask=anat / "mask.nii.gz",
                     brainmask=anat / "brainmask.nii.gz", time_years=0.0)
    fu = Timepoint(name="followup", t1=fu_dir / "synth_T1w.nii.gz", flair=fu_dir / "synth_FLAIR.nii.gz", mask=fu_dir / "synth_mask.nii.gz",
                   brainmask=None, time_years=1.0)
    params = RegParams(threads=2, engine=engine, deform_iterations="40x20x0", ants_iterations=(40, 20, 0))
    res = register_pair(base, fu, tmp_path / f"reg_{engine}", params=params)
    exp, _ = jacobian_pct_per_year(res.jacobian, res.dt_years)
    lesion = np.asarray(nib.load(res.baseline_mask).dataobj) > 0
    inside = float(np.median(exp[lesion]))
    far = float(np.median(exp[(np.asarray(nib.load(res.baseline_t1).dataobj) > 0) & ~lesion]))
    # an expansion of 1.6 over one year is +60 %/yr in the truth; either engine must see it as
    # clearly positive inside the lesion and near zero elsewhere
    assert inside > 5.0, f"{engine}: median expansion inside the lesion {inside:.1f} %/yr, expected positive"
    assert inside > far + 5.0
    record = json.loads((tmp_path / f"reg_{engine}" / "followup_reg_params.json").read_text())
    assert record["reg"]["engine"] == engine


def test_old_record_reuses_the_warp_and_recomputes_only_the_jacobian(tmp_path, monkeypatch):
    """A pair registered before the Jacobian fix (record without engine or jacobian keys)."""
    from lesiontrack import registration as R

    anat, _ = _phantom(tmp_path)
    fu_dir = tmp_path / "fu"
    make_followup(anat / "T1w.nii.gz", anat / "FLAIR.nii.gz", anat / "mask.nii.gz", fu_dir,
                  n_expand=1, volume_factors=(1.6,), min_voxels=20, rigid_rot_deg=0.5, rigid_trans_vox=0.5,
                  noise_frac=0.01, seed=1)
    base = Timepoint(name="baseline", t1=anat / "T1w.nii.gz", flair=anat / "FLAIR.nii.gz", mask=anat / "mask.nii.gz",
                     brainmask=anat / "brainmask.nii.gz", time_years=0.0)
    fu = Timepoint(name="followup", t1=fu_dir / "synth_T1w.nii.gz", flair=fu_dir / "synth_FLAIR.nii.gz",
                   mask=fu_dir / "synth_mask.nii.gz", brainmask=None, time_years=1.0)
    params = RegParams(threads=2, deform_iterations="20x0x0")
    out = tmp_path / "reg"
    res = register_pair(base, fu, out, params=params)
    record = out / "followup_reg_params.json"
    old = json.loads(record.read_text())
    old.pop("jacobian")
    for k in [k for k in old["reg"] if k == "engine" or k.startswith("ants_")]:
        old["reg"].pop(k)
    record.write_text(json.dumps(old))
    (out / "followup_jacobian.nii.gz").unlink()
    monkeypatch.setattr(R, "run_greedy", lambda *a, **k: pytest.fail("the cached registration must not be redone"))
    again = register_pair(base, fu, out, params=params)
    assert again.jacobian.exists()
    assert json.loads(record.read_text())["jacobian"] == R.JACOBIAN_VERSION
    assert np.allclose(np.asarray(nib.load(again.jacobian).dataobj), np.asarray(nib.load(res.jacobian).dataobj))
