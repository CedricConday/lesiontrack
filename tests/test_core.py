import numpy as np
import pandas as pd
import pytest
from scipy import ndimage as ndi

from lesiontrack.candidates import candidate_table, sel_candidates
from lesiontrack.config import SELParams, TrackParams
from lesiontrack.registration import half_transform
from lesiontrack.scores import cohort_score, concentricity, constancy, huber_fit
from lesiontrack.tracking import track_pair


def ball(shape, centre, r):
    g = np.indices(shape).astype(float)
    d = np.sqrt(sum((g[i] - c) ** 2 for i, c in enumerate(centre)))
    return d <= r


def test_candidates_seed_grow_and_min_size():
    shape = (40, 40, 40)
    mask = ball(shape, (20, 20, 20), 8)
    exp = np.zeros(shape, np.float32)
    exp[mask] = 5.0                      # whole lesion above JE2
    exp[ball(shape, (20, 20, 20), 3)] = 20.0  # core above JE1
    lab = sel_candidates(exp, mask, SELParams())
    assert lab.max() == 1
    assert (lab > 0).sum() == mask.sum()  # hysteresis grows to every JE2 voxel in the lesion
    # a lesion that expands only below JE1 never seeds a candidate
    exp2 = np.where(mask, 10.0, 0.0).astype(np.float32)
    assert sel_candidates(exp2, mask).max() == 0
    # a seed of fewer than 10 voxels is discarded
    exp3 = np.zeros(shape, np.float32)
    exp3[20, 20, 20] = 50.0
    assert sel_candidates(exp3, mask).max() == 0


def test_two_seeds_stay_distinct_when_connected():
    shape = (60, 30, 30)
    mask = np.zeros(shape, bool)
    mask[10:50, 10:20, 10:20] = True
    exp = np.where(mask, 6.0, 0.0).astype(np.float32)
    exp[15:20, 12:18, 12:18] = 30.0
    exp[40:45, 12:18, 12:18] = 30.0
    lab = sel_candidates(exp, mask)
    assert lab.max() == 2
    t = candidate_table(lab, exp, mask.astype(int), 1.0)
    assert list(t["candidate_id"]) == [1, 2]
    assert (t["peak_pct_per_year"] == 30.0).all()


def test_concentricity_sign():
    shape = (30, 30, 30)
    m = ball(shape, (15, 15, 15), 7)
    d = ndi.distance_transform_edt(m)
    inside_out = d.astype(np.float32) * 3          # highest at core
    outside_in = (8 - d).astype(np.float32) * 3    # highest at rim
    s1, n1 = concentricity(m, inside_out)
    s2, _ = concentricity(m, outside_in)
    assert n1 >= 3 and s1 > 0 and s2 < 0


def test_constancy_linear_vs_jumpy():
    t = np.array([0.5, 1.0, 1.5])
    lin_res, lin_slope = constancy(t, np.array([5.0, 10.0, 15.0]))
    jump_res, _ = constancy(t, np.array([0.0, 0.0, 15.0]))
    assert lin_res < 1e-6 and abs(lin_slope - 10.0) < 1e-6
    assert jump_res > lin_res
    res, _ = constancy(np.array([1.0]), np.array([10.0]))
    assert np.isnan(res)


def test_huber_resists_outlier():
    x = np.arange(10.0)
    y = 2 * x + 1
    y[7] += 100
    slope, intercept = huber_fit(x, y)
    assert abs(slope - 2) < 0.2 and abs(intercept - 1) < 1.5


def test_cohort_score_two_timepoint_fallback():
    df = pd.DataFrame({"concentricity": [1.0, -1.0, 0.0], "constancy_residual": [np.nan] * 3})
    out = cohort_score(df)
    assert (out["score_basis"] == "concentricity_only").all()
    assert out["definite_sel"].tolist() == [True, False, True]


def test_tracking_classes():
    shape = (60, 60, 60)
    bl = np.zeros(shape, bool)
    fu = np.zeros(shape, bool)
    bl |= ball(shape, (15, 15, 15), 5)            # stable
    fu |= ball(shape, (15, 15, 15), 5)
    bl |= ball(shape, (40, 15, 15), 4)            # enlarging
    fu |= ball(shape, (40, 15, 15), 6)
    bl |= ball(shape, (15, 40, 15), 6)            # shrinking
    fu |= ball(shape, (15, 40, 15), 4)
    bl |= ball(shape, (40, 40, 15), 4)            # resolved
    fu |= ball(shape, (30, 30, 45), 4)            # new
    table, maps = track_pair(bl, fu, 1.0, 1.0, TrackParams())
    classes = sorted(table["class"])
    assert classes == ["enlarging", "new", "resolved", "shrinking", "stable"]
    assert maps["new_or_enlarging_voxels"].sum() > 0


def test_tracking_merge_and_split_flags():
    shape = (40, 40, 40)
    bl = ball(shape, (12, 20, 20), 4) | ball(shape, (24, 20, 20), 4)
    fu = ball(shape, (18, 20, 20), 9)
    table, _ = track_pair(bl, fu, 1.0, 1.0)
    assert len(table) == 1 and bool(table.loc[0, "merged"]) and not bool(table.loc[0, "split"])


def test_half_transform_squares_back():
    from scipy.spatial.transform import Rotation
    m = np.eye(4)
    m[:3, :3] = Rotation.from_euler("xyz", [3, -2, 5], degrees=True).as_matrix()
    m[:3, 3] = [1.5, -2.0, 0.7]
    h = half_transform(m)
    assert np.allclose(h @ h, m, atol=1e-9)


@pytest.mark.parametrize("bad", [0.0, -1.0])
def test_dt_must_be_positive(bad):
    with pytest.raises(ValueError):
        track_pair(np.zeros((5, 5, 5), bool), np.zeros((5, 5, 5), bool), bad, 1.0)


def test_paths_with_spaces_are_refused(tmp_path):
    from lesiontrack.registration import Timepoint, register_pair

    d = tmp_path / "has space"
    d.mkdir()
    tp = Timepoint("a", d / "t1.nii.gz", d / "fl.nii.gz", d / "m.nii.gz", 0.0)
    with pytest.raises(ValueError, match="spaces"):
        register_pair(tp, tp, tmp_path / "out")


def test_reslice_keeps_masks_binary(tmp_path):
    """Regression: greedy's -ri applies to the next -rm pair; a mask must not be interpolated."""
    import nibabel as nib

    from lesiontrack.registration import _reslice, write_matrix

    shape = (24, 24, 24)
    m = ball(shape, (12, 12, 12), 6).astype(np.uint8)
    aff = np.eye(4)
    nib.save(nib.Nifti1Image(m, aff), tmp_path / "mask.nii.gz")
    mat = np.eye(4)
    mat[:3, 3] = [0.4, -0.3, 0.2]  # sub-voxel shift forces interpolation if it were linear
    write_matrix(tmp_path / "shift.mat", mat)
    _reslice(tmp_path / "mask.nii.gz", tmp_path / "mask.nii.gz", tmp_path / "out.nii.gz",
             str(tmp_path / "shift.mat"), "NN", 1, None)
    out = np.asarray(nib.load(tmp_path / "out.nii.gz").dataobj)
    assert set(np.unique(out).tolist()) <= {0, 1}
    assert abs(int((out > 0).sum()) - int(m.sum())) <= 0.1 * m.sum()


def test_synthetic_jacobian_matches_nominal_factor():
    from lesiontrack.synth import _jacobian_det, _radial_displacement

    shape = (41, 41, 41)
    c = np.array([20.0, 20.0, 20.0])
    for f in (1.15, 1.4, 2.0):
        d = _radial_displacement(shape, c, 10, f, 6)
        expansion = 1.0 / _jacobian_det(d)
        assert abs(expansion[20, 20, 20] - f) < 1e-3
        assert abs(expansion[25, 20, 20] - f) < 1e-3   # still inside the uniform window
        assert abs(expansion[0, 0, 0] - 1.0) < 1e-6    # untouched far away
