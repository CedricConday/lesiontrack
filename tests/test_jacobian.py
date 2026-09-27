"""The Jacobian determinant is computed in the voxel frame, whatever the image orientation."""
import numpy as np
import nibabel as nib
import pytest

from lesiontrack.registration import jacobian_determinant


def _radial_field(shape, centre, radius, k):
    """u(x) = k (x - c) inside ``radius`` (voxel units), fading to zero over 6 voxels."""
    g = np.indices(shape).astype(float)
    d = g - np.asarray(centre)[:, None, None, None]
    r = np.sqrt((d**2).sum(0))
    w = np.clip(1.0 - (r - radius) / 6.0, 0.0, 1.0)
    w = np.where(r <= radius, 1.0, w)
    return np.moveaxis(d * (k * w)[None], 0, -1), r


@pytest.mark.parametrize("diag", [(1, 1, 1), (-1, 1, 1), (-1, -1, 1), (1, -1, 1)])
def test_determinant_is_independent_of_orientation(tmp_path, diag):
    shape = (40, 40, 40)
    k = 0.25  # linear stretch 1.25 -> volume 1.953
    u_vox, r = _radial_field(shape, (20, 20, 20), 8, k)
    aff = np.diag(list(diag) + [1.0]).astype(float)
    aff[:3, 3] = [0 if d > 0 else shape[i] - 1 for i, d in enumerate(diag)]
    # store the field the ITK way: physical LPS millimetres
    u_ras = np.einsum("ij,xyzj->xyzi", aff[:3, :3], u_vox)
    u_lps = u_ras * np.array([-1.0, -1.0, 1.0])
    warp = tmp_path / "warp.nii.gz"
    nib.save(nib.Nifti1Image(u_lps[:, :, :, None, :].astype(np.float32), aff), warp)
    det = np.asarray(jacobian_determinant(warp, tmp_path / "jac.nii.gz").dataobj)
    inside = r < 6
    assert det[inside] == pytest.approx((1 + k) ** 3, rel=1e-3)
    assert np.median(det[r > 20]) == pytest.approx(1.0, abs=1e-6)
    assert (tmp_path / "jac.nii.gz").exists()


def test_anisotropic_spacing_is_honoured(tmp_path):
    shape = (40, 40, 20)
    k = 0.2
    u_vox, r = _radial_field(shape, (20, 20, 10), 6, k)
    aff = np.diag([1.0, 1.0, 2.0, 1.0])
    u_lps = np.einsum("ij,xyzj->xyzi", aff[:3, :3], u_vox) * np.array([-1.0, -1.0, 1.0])
    warp = tmp_path / "warp.nii.gz"
    nib.save(nib.Nifti1Image(u_lps[:, :, :, None, :].astype(np.float32), aff), warp)
    det = np.asarray(jacobian_determinant(warp).dataobj)
    assert det[r < 4] == pytest.approx((1 + k) ** 3, rel=1e-3)


def test_rejects_a_non_field(tmp_path):
    p = tmp_path / "scalar.nii.gz"
    nib.save(nib.Nifti1Image(np.zeros((8, 8, 8), np.float32), np.eye(4)), p)
    with pytest.raises(ValueError):
        jacobian_determinant(p)


def test_rotated_sheared_anisotropic_affine(tmp_path):
    """The whole point: a general direction matrix with anisotropic spacing."""
    from scipy.spatial.transform import Rotation

    shape = (40, 40, 24)
    k = 0.2
    u_vox, r = _radial_field(shape, (20, 20, 12), 6, k)
    aff = np.eye(4)
    aff[:3, :3] = Rotation.from_euler("xyz", [12, -7, 25], degrees=True).as_matrix() @ np.diag([0.8, 1.0, 2.5]) @ np.diag([-1, 1, 1])
    aff[:3, 3] = [10, -20, 5]
    u_lps = np.einsum("ij,xyzj->xyzi", aff[:3, :3], u_vox) * np.array([-1.0, -1.0, 1.0])
    warp = tmp_path / "warp.nii.gz"
    nib.save(nib.Nifti1Image(u_lps[:, :, :, None, :].astype(np.float32), aff), warp)
    det = np.asarray(jacobian_determinant(warp).dataobj)
    assert det[r < 4] == pytest.approx((1 + k) ** 3, rel=1e-3)
    assert np.median(det[r > 16]) == pytest.approx(1.0, abs=1e-6)
