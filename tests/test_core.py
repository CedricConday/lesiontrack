from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd
import pytest
from scipy import ndimage as ndi

from lesiontrack.candidates import candidate_table, sel_candidates
from lesiontrack.config import SELParams, TrackParams
from lesiontrack.registration import half_transform
from lesiontrack.scores import cohort_score, concentricity, constancy, huber_fit
from lesiontrack.tracking import track_pair

ELLIOTT = SELParams(min_mean_pct_per_year=0.0, min_mean_total_pct=0.0)  # Elliott 2019 as published, no mean-rate filter


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
    lab = sel_candidates(exp, mask, ELLIOTT)
    assert lab.max() == 1
    assert (lab > 0).sum() == mask.sum()  # hysteresis grows to every JE2 voxel in the lesion
    # a lesion that expands only below JE1 never seeds a candidate
    exp2 = np.where(mask, 10.0, 0.0).astype(np.float32)
    assert sel_candidates(exp2, mask, ELLIOTT).max() == 0
    # a seed of fewer than 10 voxels is discarded
    exp3 = np.zeros(shape, np.float32)
    exp3[20, 20, 20] = 50.0
    assert sel_candidates(exp3, mask, ELLIOTT).max() == 0


def test_two_seeds_stay_distinct_when_connected():
    shape = (60, 30, 30)
    mask = np.zeros(shape, bool)
    mask[10:50, 10:20, 10:20] = True
    exp = np.where(mask, 6.0, 0.0).astype(np.float32)
    exp[15:20, 12:18, 12:18] = 30.0
    exp[40:45, 12:18, 12:18] = 30.0
    lab = sel_candidates(exp, mask, ELLIOTT)
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


def test_min_mean_drops_weak_candidates():
    shape = (40, 40, 40)
    mask = ball(shape, (20, 20, 20), 8)
    exp = np.where(mask, 5.0, 0.0).astype(np.float32)
    exp[ball(shape, (20, 20, 20), 2)] = 13.0   # one small seed, candidate mean stays near 5
    assert sel_candidates(exp, mask, ELLIOTT).max() == 1
    assert sel_candidates(exp, mask, SELParams(min_mean_pct_per_year=8.0)).max() == 0
    # the default is a floor on total change: 10 % over the interval
    assert sel_candidates(exp, mask, dt_years=1.0).max() == 0   # 10 %/yr at one year
    assert sel_candidates(exp, mask, dt_years=4.0).max() == 1   # 2.5 %/yr at four years
    assert sel_candidates(exp, mask).max() == 1                 # no interval given: no floor to scale
    both = SELParams(min_mean_pct_per_year=8.0)                 # an absolute floor applies as well
    assert sel_candidates(exp, mask, both, dt_years=4.0).max() == 0
    assert SELParams().min_mean_at(2.0) == 5.0 and both.min_mean_at(2.0) == 8.0 and ELLIOTT.min_mean_at(0.5) == 0.0


def test_zero_margin_does_not_flood_the_baseline():
    shape = (50, 50, 50)
    bl = ball(shape, (15, 15, 15), 5)
    fu = ball(shape, (15, 15, 15), 5) | ball(shape, (38, 38, 38), 4)
    table, maps = track_pair(bl, fu, 1.0, 1.0, TrackParams(new_lesion_margin_voxels=0))
    assert sorted(table["class"]) == ["new", "stable"]
    assert maps["new_or_enlarging_voxels"].sum() > 0


def test_unscorable_candidates_are_never_definite():
    df = pd.DataFrame({"concentricity": [2.0, np.nan, -1.0], "constancy_residual": [0.1, 0.1, 0.5]})
    out = cohort_score(df)
    assert out.loc[1, "score_basis"] == "unscorable" and not bool(out.loc[1, "definite_sel"])
    assert bool(out.loc[0, "definite_sel"])


def test_swapped_je_thresholds_are_rejected():
    with pytest.raises(ValueError):
        SELParams(je1_pct_per_year=4.0, je2_pct_per_year=12.5)


def test_synthetic_composition_keeps_window_on_lesion():
    """With rigid motion, the injected expansion must still be centred on the lesion."""
    from lesiontrack.synth import _jacobian_det, _radial_displacement, _rigid_field

    shape = (61, 61, 61)
    c = np.array([45.0, 45.0, 45.0])  # off-centre so a rotation moves it
    rigid = _rigid_field(shape, np.array([1.0, -1.0, 1.0]), np.array([1.0, 0.5, -0.5]))
    grids = np.stack(np.meshgrid(*[np.arange(n, dtype=np.float32) for n in shape], indexing="ij"))
    y = grids + rigid
    d_at_y = _radial_displacement(shape, c, 6, 1.4, 6, positions=y)
    # the radial part evaluated at y is a function of y - c: exact scaling where |y - c| <= 6
    dist = np.sqrt(((y - c[:, None, None, None]) ** 2).sum(0))
    inside = dist <= 5
    lam = 1.4 ** (1 / 3)
    expected = (y - c[:, None, None, None]) * (1 / lam - 1)
    assert np.allclose(d_at_y[:, inside], expected[:, inside], atol=1e-4)
    assert abs(1 / _jacobian_det(_radial_displacement(shape, c, 6, 1.4, 6))[45, 45, 45] - 1.4) < 1e-3


def test_empty_image_is_refused(tmp_path):
    import nibabel as nib

    from lesiontrack.registration import RegistrationError, check_image_content

    a = np.zeros((20, 20, 20), np.float32)
    a[10, 10, 10] = 5.0
    nib.save(nib.Nifti1Image(a, np.eye(4)), tmp_path / "empty.nii.gz")
    with pytest.raises(RegistrationError, match="nonzero"):
        check_image_content(tmp_path / "empty.nii.gz")
    a[5:15, 5:15, 5:15] = 100.0
    nib.save(nib.Nifti1Image(a, np.eye(4)), tmp_path / "ok.nii.gz")
    assert check_image_content(tmp_path / "ok.nii.gz") > 0.1


def test_rigid_qc_rejects_implausible_transform():
    from scipy.spatial.transform import Rotation

    from lesiontrack.registration import RegistrationError, rigid_qc

    m = np.eye(4)
    m[:3, :3] = Rotation.from_euler("xyz", [2, 1, -3], degrees=True).as_matrix()
    m[:3, 3] = [3, -2, 5]
    qc = rigid_qc(m)
    assert 3.0 < qc["rotation_deg"] < 4.5 and abs(qc["translation_mm"] - np.sqrt(38)) < 1e-6
    bad = np.eye(4)
    bad[:3, :3] = Rotation.from_euler("xyz", [50, 0, 0], degrees=True).as_matrix()
    with pytest.raises(RegistrationError, match="implausible"):
        rigid_qc(bad)


def test_feature_auc_orients_each_feature_as_the_published_rule_reads_it():
    import pandas as pd

    from lesiontrack.synth import feature_auc

    cands = pd.DataFrame({
        "injected": [True, True, True, False, False, False],
        "mean_pct_per_year": [20.0, 18.0, 15.0, 6.0, 5.0, 7.0],       # separates perfectly
        "concentricity": [1.0, 2.0, 3.0, 1.0, 2.0, 3.0],              # no information
        "constancy_residual": [0.01, 0.02, float("nan"), 0.5, 0.4, 0.3],  # low = injected
    })
    auc = feature_auc(cands)
    assert auc["mean_pct_per_year"]["auc"] == 1.0
    assert auc["concentricity"]["auc"] == 0.5
    assert auc["constancy_residual"]["auc"] == 1.0 and auc["constancy_residual"]["direction"] == "low"
    assert auc["constancy_residual"]["n_scored"] == 5 and auc["constancy_residual"]["n_injected"] == 2
    assert feature_auc(cands.iloc[:0]) == {}


def _toy_subject(tmp_path, shape=(48, 48, 48)):
    rng = np.random.default_rng(0)
    grid = np.stack(np.meshgrid(*[np.arange(n) for n in shape], indexing="ij"))
    brain = ((grid - 23.5) ** 2).sum(axis=0) < 20**2
    mask = np.zeros(shape, bool)
    for c in ((16, 24, 24), (32, 24, 24), (24, 14, 24), (24, 34, 24)):
        mask |= ((grid - np.array(c)[:, None, None, None]) ** 2).sum(axis=0) < 3.2**2
    t1 = np.where(brain, 100.0, 0.0) - 30.0 * mask + rng.normal(0, 1, shape) * brain
    flair = np.where(brain, 80.0, 0.0) + 60.0 * mask + rng.normal(0, 1, shape) * brain
    paths = {}
    for name, arr in (("t1", t1), ("flair", flair), ("mask", mask.astype(np.uint8))):
        paths[name] = tmp_path / f"{name}.nii.gz"
        nib.save(nib.Nifti1Image(arr.astype(np.float32), np.eye(4)), paths[name])
    return paths


def test_scanner_scale_enlarges_the_followup_but_not_the_truth(tmp_path):
    from lesiontrack.synth import make_followup

    p = _toy_subject(tmp_path)
    kw = {"n_expand": 2, "min_voxels": 30, "rigid_rot_deg": 0.0, "rigid_trans_vox": 0.0, "noise_frac": 0.0, "seed": 1}
    plain = make_followup(p["t1"], p["flair"], p["mask"], tmp_path / "plain", **kw)
    scaled = make_followup(p["t1"], p["flair"], p["mask"], tmp_path / "scaled", scanner_scale_pct=5.0, **kw)
    vol = lambda rec: int(np.asarray(nib.load(rec["outputs"]["brainmask"]).dataobj).sum())
    growth = vol(scaled) / vol(plain)
    # 5 % linear is about 16 % in volume; a smaller step is lost in the binary mask's voxel count
    assert abs(growth - 1.05**3) < 0.02
    truth = lambda rec: np.asarray(nib.load(rec["outputs"]["truth_expansion_pct"]).dataobj)
    assert np.array_equal(truth(plain), truth(scaled))   # the nuisance is not part of the truth
    assert scaled["scanner_scale_pct"] == 5.0


def test_thicken_slices_blurs_through_plane_only():
    from lesiontrack.synth import thicken_slices

    flat = np.full((6, 6, 12), 7.0, np.float32)
    assert np.allclose(thicken_slices(flat, 3), 7.0)
    assert thicken_slices(flat, 1) is flat
    inplane = np.zeros((6, 6, 12), np.float32)
    inplane[2] = 1.0                                  # varies in-plane, constant through-plane
    assert np.allclose(thicken_slices(inplane, 3), inplane)
    spike = np.zeros((6, 6, 12), np.float32)
    spike[:, :, 4] = 3.0                              # one thin slice
    out = thicken_slices(spike, 3)
    assert out.shape == spike.shape
    assert out[0, 0].max() < 1.01                     # averaged over the 3-slice slab
    assert out[0, 0, 4] == 1.0 and out[0, 0, 0] == 0.0


def test_global_scale_reads_a_uniform_scale_change_away_from_lesions():
    from lesiontrack.registration import global_scale, reference_mask

    shape = (48, 48, 48)
    brain = ball(shape, (24, 24, 24), 20)
    lesion = ball(shape, (24, 24, 24), 4)
    s3 = 1.01 ** 3
    jac = np.where(brain, s3, 1.0).astype(np.float32)
    jac[lesion] = s3 * 1.25
    ref = reference_mask(brain, lesion)
    assert not (ref & ndi.binary_dilation(lesion, iterations=3)).any()
    assert global_scale(jac, ref) == pytest.approx(s3, rel=1e-6)


def test_global_scale_resists_local_change_and_refuses_a_tiny_reference():
    from lesiontrack.registration import global_scale

    jac = np.full((20, 20, 20), 1.02, np.float32)
    jac[:1] = 3.0                       # 5 % of the reference changes a lot locally
    jac[1:2] = 0.5
    jac[2:8] = 1.0                      # a spike at exactly 1 must not capture the estimate
    est = global_scale(jac, np.ones_like(jac, bool))
    assert 1.01 < est < 1.02
    with pytest.raises(ValueError, match="reference region"):
        global_scale(jac, np.zeros_like(jac, bool))


def test_defaults_are_mean_rate_filter_and_similarity_step():
    from lesiontrack.config import RegParams

    assert SELParams().min_mean_total_pct == 10.0 and SELParams().min_mean_pct_per_year == 0.0
    assert RegParams().affine_dof == 7


def test_qc_reads_scale_out_of_a_similarity_and_refuses_a_large_one():
    from lesiontrack.registration import RegistrationError, rigid_qc

    m = np.eye(4)
    m[:3, :3] *= 1.01
    q = rigid_qc(m)
    assert q["scale_pct"] == pytest.approx(1.0, abs=1e-9)
    assert q["rotation_deg"] == pytest.approx(0.0, abs=1e-6)
    m[:3, :3] *= 1.05
    with pytest.raises(RegistrationError, match="scale"):
        rigid_qc(m)


def _rot(axis, deg, t):
    from scipy.spatial.transform import Rotation

    m = np.eye(4)
    m[:3, :3] = Rotation.from_euler(axis, deg, degrees=True).as_matrix()
    m[:3, 3] = t
    return m


def test_earlier_followups_land_on_the_baseline_in_the_shared_frame():
    from lesiontrack.registration import followup_to_frame

    # greedy: output(x) = moving(M x). A pair's rigid A maps baseline points to follow-up
    # points; the frame is the last pair's halfway H, and the baseline is resliced with H^-1.
    a_last = _rot("z", 3.0, [8.0, -2.0, 1.0])
    a_mid = _rot("x", 2.0, [-3.0, 5.0, 2.0])
    h = half_transform(a_last)
    x = np.array([40.0, -25.0, 30.0, 1.0])
    baseline_point = np.linalg.inv(h) @ x
    fu = followup_to_frame(a_mid, h)
    assert np.allclose(fu @ x, a_mid @ baseline_point)
    assert np.allclose(followup_to_frame(a_last, h), h)  # the last pair keeps its own halfway
    assert not np.allclose(np.linalg.inv(h) @ a_mid, fu)  # the pre-0.3.0 order differs


def test_cached_framed_pairs_from_before_the_fix_are_not_reused():
    from lesiontrack.config import RegParams
    from lesiontrack.registration import _params_record, _same_registration

    cur = _params_record(RegParams(), Path("frame.mat"))
    old = {k: v for k, v in cur.items() if k != "frame_version"}
    assert not _same_registration(old, cur)
    unframed = _params_record(RegParams(), None)
    assert _same_registration({k: v for k, v in unframed.items() if k != "frame_version"}, unframed)


def test_flair_support_flags_new_and_resolved_by_intensity_change():
    from lesiontrack.tracking import flair_support

    rng = np.random.default_rng(0)
    shape = (40, 40, 40)
    bl = np.zeros(shape, bool)
    fu = np.zeros(shape, bool)
    bl[5:9, 5:9, 5:9] = True      # resolved, FLAIR really drops
    bl[5:9, 20:24, 5:9] = True    # resolved by the masks only
    fu[20:24, 5:9, 5:9] = True    # new, FLAIR really rises
    fu[20:24, 20:24, 5:9] = True  # new by the masks only
    bl[30:34, 30:34, 30:34] = fu[30:34, 30:34, 30:34] = True  # persistent
    f0 = 100 + rng.normal(0, 5, shape).astype(np.float32)
    f1 = 100 + rng.normal(0, 5, shape).astype(np.float32)
    f0[5:9, 5:9, 5:9] += 40
    f1[20:24, 5:9, 5:9] += 40
    f0[30:34, 30:34, 30:34] += 40
    f1[30:34, 30:34, 30:34] += 40
    table, maps = track_pair(bl, fu, 1.0, 1.0)
    out = flair_support(table, maps, f0, f1, np.ones(shape, bool))
    by = {(r["class"], round(r["cx"]), round(r["cy"])): r for _, r in out.iterrows()}
    assert by[("resolved", 6, 6)]["flair_supported"] and by[("resolved", 6, 6)]["flair_dz"] < -5
    assert not by[("resolved", 6, 22)]["flair_supported"]
    assert by[("new", 22, 6)]["flair_supported"] and by[("new", 22, 6)]["flair_dz"] > 5
    assert not by[("new", 22, 22)]["flair_supported"]
    persistent = out[out["class"] == "stable"].iloc[0]
    assert pd.isna(persistent["flair_supported"]) and abs(persistent["flair_dz"]) < 1
    # a uniform intensity rescale between scans is not a change
    out2 = flair_support(table, maps, f0, f1 * 1.3, np.ones(shape, bool))
    assert list(out2["flair_supported"].fillna(False)) == list(out["flair_supported"].fillna(False))


def test_flair_support_empty_table():
    from lesiontrack.tracking import flair_support

    z = np.zeros((6, 6, 6), bool)
    table, maps = track_pair(z, z, 1.0, 1.0)
    out = flair_support(table, maps, z.astype(float), z.astype(float), ~z)
    assert len(out) == 0 and {"flair_dz", "flair_supported"} <= set(out.columns)
