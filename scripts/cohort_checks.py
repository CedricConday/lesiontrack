"""Two consistency checks on a lesiontrack output tree.

1. Mask-based versus Jacobian-based volume change, per subject and per 1:1 matched lesion.
2. Are the voxels the masks call new already FLAIR-bright at baseline, and are resolved
   voxels still bright at follow-up? Robust z within non-lesion brain; "bright" is the dimmer
   quartile of persistent lesion voxels.
"""

import json
import sys
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd
from scipy import ndimage as ndi
from scipy.stats import spearmanr

out = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("derivatives/mslesseg")
subj_rows, les_rows, flair_rows = [], [], []
for d in sorted(p.parent for p in out.glob("*/summary.json")):
    s = json.load(open(d / "summary.json"))
    last, dt = s["follow_ups"][-1], s["dt_years"][-1]
    R = d / "reg"
    exp = np.asarray(nib.load(d / f"expansion_{last}_pct_per_year.nii.gz").dataobj)
    labs = np.asarray(nib.load(d / "baseline_lesion_labels.nii.gz").dataobj)
    bl, fu = labs > 0, np.asarray(nib.load(d / f"lesion_labels_{last}.nii.gz").dataobj) > 0
    brain = np.asarray(nib.load(R / f"{last}_baseline_brainmask_halfway.nii.gz").dataobj) > 0
    subj_rows.append({"subject": s["subject"], "dt_years": round(dt, 2),
                      "mask_change_pct": (fu.sum() / bl.sum() - 1) * 100,
                      "jacobian_change_pct": float(exp[bl].mean() * dt),
                      "brain_mean_pct": float(exp[brain].mean() * dt)})
    t = pd.read_csv(d / "lesion_tracking.tsv", sep="\t")
    t = t[(t.follow_up == last) & (t.n_baseline == 1) & (t.n_followup == 1)]
    for _, r in t.iterrows():
        sel = labs == int(r.baseline_ids)
        les_rows.append({"subject": s["subject"], "n_voxels": int(sel.sum()),
                         "mask_change_pct": (r.volume_followup_mm3 / r.volume_baseline_mm3 - 1) * 100,
                         "jacobian_change_pct": float(exp[sel].mean() * dt)})
    bl_fl = np.asarray(nib.load(R / f"{last}_baseline_FLAIR_halfway.nii.gz").dataobj).astype(np.float32)
    fu_fl = np.asarray(nib.load(R / f"{last}_followup_FLAIR_halfway.nii.gz").dataobj).astype(np.float32)

    def z(img):
        ref = img[brain & ~bl & ~fu]
        med = np.median(ref)
        mad = np.median(np.abs(ref - med)) * 1.4826
        return (img - med) / mad

    zb, zf = z(bl_fl), z(fu_fl)
    new = fu & ~ndi.binary_dilation(bl, iterations=1)
    resolved = bl & ~ndi.binary_dilation(fu, iterations=1)
    persistent = bl & fu
    thr = np.percentile(zb[persistent], 25)
    flair_rows.append({"subject": s["subject"], "n_new_voxels": int(new.sum()), "n_resolved_voxels": int(resolved.sum()),
                       "new_z_baseline": float(np.median(zb[new])), "new_z_followup": float(np.median(zf[new])),
                       "new_already_bright_frac": float((zb[new] >= thr).mean()),
                       "resolved_z_baseline": float(np.median(zb[resolved])), "resolved_z_followup": float(np.median(zf[resolved])),
                       "resolved_still_bright_frac": float((zf[resolved] >= thr).mean()),
                       "persistent_z": float(np.median(zb[persistent]))})

subj = pd.DataFrame(subj_rows)
les = pd.DataFrame(les_rows)
fl = pd.DataFrame(flair_rows)
subj.to_csv(out / "check_volume_change_by_subject.tsv", sep="\t", index=False, float_format="%.4g")
les.to_csv(out / "check_volume_change_by_lesion.tsv", sep="\t", index=False, float_format="%.4g")
fl.to_csv(out / "check_flair_consistency.tsv", sep="\t", index=False, float_format="%.4g")
big = les[les.n_voxels >= 100]
print(f"subjects {len(subj)}: spearman(mask, jacobian) = {spearmanr(subj.mask_change_pct, subj.jacobian_change_pct)[0]:.2f}; "
      f"median |mask| {subj.mask_change_pct.abs().median():.1f} %, median |jacobian| {subj.jacobian_change_pct.abs().median():.1f} %, "
      f"brain mean {subj.brain_mean_pct.median():.1f} %")
print(f"lesions 1:1 matched {len(les)}: spearman = {spearmanr(les.mask_change_pct, les.jacobian_change_pct)[0]:.2f}; "
      f">= 100 voxels ({len(big)}): {spearmanr(big.mask_change_pct, big.jacobian_change_pct)[0]:.2f}")
print(f"FLAIR: new voxels z {fl.new_z_baseline.median():.2f} -> {fl.new_z_followup.median():.2f} (already bright {fl.new_already_bright_frac.median():.2f}); "
      f"resolved voxels z {fl.resolved_z_baseline.median():.2f} -> {fl.resolved_z_followup.median():.2f} (still bright {fl.resolved_still_bright_frac.median():.2f}); "
      f"persistent lesion z {fl.persistent_z.median():.2f}")
