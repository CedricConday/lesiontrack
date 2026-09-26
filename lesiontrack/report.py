"""One-page overview figure per subject."""

from __future__ import annotations

from pathlib import Path

import matplotlib
import nibabel as nib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def _load(p: Path) -> np.ndarray:
    return np.asarray(nib.load(p).dataobj)


def overview(out_dir: Path, follow_up: str, n_slices: int = 4) -> Path:
    """Axial slices through the largest lesions: T1, expansion in lesions, SEL outlines."""
    t1 = _load(out_dir / "reg" / f"{follow_up}_baseline_T1w_halfway.nii.gz").astype(np.float32)
    fu = _load(out_dir / "reg" / f"{follow_up}_followup_T1w_warped.nii.gz").astype(np.float32)
    exp = _load(out_dir / f"expansion_{follow_up}_pct_per_year.nii.gz").astype(np.float32)
    les = _load(out_dir / "baseline_lesion_labels.nii.gz")
    cand = _load(out_dir / "sel_candidates.nii.gz")
    new = _load(out_dir / f"new_or_enlarging_{follow_up}.nii.gz")

    # Slices with the most lesion voxels, at least min_gap apart so the page spans the brain.
    lesion_per_slice = (les > 0).sum(axis=(0, 1))
    min_gap = max(3, les.shape[2] // (4 * n_slices))
    slices: list[int] = []
    for z in np.argsort(lesion_per_slice)[::-1]:
        if lesion_per_slice[z] == 0:
            break
        if all(abs(int(z) - s) >= min_gap for s in slices):
            slices.append(int(z))
        if len(slices) == n_slices:
            break
    slices.sort()
    vmax = np.percentile(t1[t1 > 0], 99) if (t1 > 0).any() else 1.0

    fig, axes = plt.subplots(len(slices), 3, figsize=(10, 3.3 * len(slices)))
    axes = np.atleast_2d(axes)
    for row, z in zip(axes, slices):
        for ax, img, title in zip(row, (t1, fu, t1), ("baseline T1w", f"{follow_up} T1w warped to baseline", "expansion in lesions (%/yr)")):
            ax.imshow(np.rot90(img[:, :, z]), cmap="gray", vmin=0, vmax=vmax)
            ax.set_title(f"{title}  z={z}", fontsize=8)
            ax.axis("off")
        m = np.rot90(les[:, :, z] > 0)
        row[0].contour(m, levels=[0.5], colors="cyan", linewidths=0.6)
        row[1].contour(np.rot90(new[:, :, z] > 0), levels=[0.5], colors="yellow", linewidths=0.6)
        heat = np.rot90(np.where(les[:, :, z] > 0, exp[:, :, z], np.nan))
        row[2].imshow(heat, cmap="coolwarm", vmin=-25, vmax=25, alpha=0.9)
        row[2].contour(np.rot90(cand[:, :, z] > 0), levels=[0.5], colors="lime", linewidths=0.8)
    fig.suptitle("cyan: baseline lesions   yellow: new or enlarging voxels   lime: SEL candidates", fontsize=9)
    fig.tight_layout()
    path = out_dir / f"overview_{follow_up}.png"
    fig.savefig(path, dpi=110)
    plt.close(fig)
    return path
