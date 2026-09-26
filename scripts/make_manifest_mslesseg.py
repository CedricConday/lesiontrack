"""Build a lesiontrack manifest for the longitudinal patients of MSLesSeg.

MSLesSeg (Guarnera et al. 2025, Sci Data 12:920, CC BY) ships every timepoint
already in MNI152 1 mm space, skull-stripped, with a lesion mask per timepoint
and the patient's age at each scan (two decimals in clinical_data.csv), so the
interval between scans is the age difference.
"""

import sys
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd

root = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("sourcedata/MSLesSeg Dataset")
out = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("work/manifest_mslesseg.tsv")
clin = pd.read_csv(root / "info_dataset" / "clinical_data.csv", sep=";", decimal=",", encoding="utf-8-sig")
clin = clin.rename(columns=lambda c: c.strip())
rows = []
for (pat, tp), r in clin.groupby(["Patient", "Timepoint"]):
    d = root / "train" / pat / tp
    if not d.exists():
        continue
    t1 = d / f"{pat}_{tp}_T1.nii.gz"
    bm = d / f"{pat}_{tp}_brainmask.nii.gz"
    if not bm.exists():
        img = nib.load(t1)
        nib.save(nib.Nifti1Image((np.asarray(img.dataobj) > 0).astype(np.uint8), img.affine), bm)
    rows.append({
        "subject": pat, "session": tp, "time_years": float(r["Age"].iloc[0]),
        "t1": t1.absolute(), "flair": (d / f"{pat}_{tp}_FLAIR.nii.gz").absolute(),
        "mask": (d / f"{pat}_{tp}_MASK.nii.gz").absolute(), "brainmask": bm.absolute(),
    })
m = pd.DataFrame(rows)
counts = m.groupby("subject")["session"].count()
m = m[m["subject"].isin(counts[counts >= 2].index)]
m = m.sort_values(["subject", "time_years"])
m.to_csv(out, sep="\t", index=False)
print(f"{m['subject'].nunique()} longitudinal subjects, {len(m)} timepoints -> {out}")
print(m.groupby("subject")["time_years"].agg(lambda s: round(s.max() - s.min(), 2)).describe())
