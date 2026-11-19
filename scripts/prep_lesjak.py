#!/usr/bin/env python
"""Prepare the Lesjak et al. 2016 longitudinal database for lesiontrack.

The database ships two co-registered studies per patient (T1W, T2W, FLAIR at about
0.7 x 0.7 x 3 mm), a brain mask and an expert mask of lesion *change*, but no lesion
masks. This script skull-strips T1W and FLAIR with the shipped brain mask, segments
lesions in every study with LST-AI v2 (docker, CPU, fast mode) and writes a lesiontrack
manifest. Usage: prep_lesjak.py [database_dir] [work_dir] [threads]
"""
import csv
import os
import subprocess
import sys
import time
from pathlib import Path

import nibabel as nib
import numpy as np

IMAGE = "jqmcginnis/lst-ai:v2.0.0rc1-cpu"
SEG = "space-flair_seg-lst.nii.gz"

db = Path(sys.argv[1] if len(sys.argv) > 1 else "sourcedata/lesjak/longitudinal").resolve()
work = Path(sys.argv[2] if len(sys.argv) > 2 else "work/lesjak").resolve()
threads = sys.argv[3] if len(sys.argv) > 3 else "4"
days = {int(r["patient_id"]): float(r["days_between_studies"]) for r in csv.DictReader(open(db / "patient_info.csv"))}
rows = []
for pdir in sorted((db / "coregistered").glob("patient*")):
    pid = int(pdir.name.replace("patient", ""))
    out = work / "input" / pdir.name
    out.mkdir(parents=True, exist_ok=True)
    brain_img = nib.load(pdir / "brainmask.nii.gz")
    brain = np.asarray(brain_img.dataobj) > 0
    nib.save(nib.Nifti1Image(brain.astype(np.uint8), brain_img.affine), out / "brainmask.nii.gz")
    for study, t in (("study1", 0.0), ("study2", days[pid] / 365.25)):
        paths = {}
        for mod, name in (("T1W", "T1w"), ("FLAIR", "FLAIR")):
            src = nib.load(pdir / f"{study}_{mod}.nii.gz")
            if src.shape != brain.shape or not np.allclose(src.affine, brain_img.affine, atol=1e-3):
                raise SystemExit(f"{pdir.name} {study} {mod}: not on the brain mask's grid")
            arr = np.asarray(src.dataobj, dtype=np.float32)
            # N4 and resampling left small negative values; inside the brain every voxel stays
            # positive so that "nonzero" still means "brain" downstream.
            arr = np.where(brain, np.maximum(arr, 1.0), 0.0).astype(np.float32)
            paths[name] = out / f"{study}_{name}.nii.gz"
            if not paths[name].exists():
                nib.save(nib.Nifti1Image(arr, src.affine), paths[name])
        seg_dir = work / "lstai" / f"{pdir.name}_{study}"
        mask = out / f"{study}_mask.nii.gz"
        if not mask.exists():
            (seg_dir / "tmp").mkdir(parents=True, exist_ok=True)
            name = f"lesjak-{pdir.name}-{study}-{os.getpid()}"
            t0 = time.time()
            res = subprocess.run([
                "docker", "run", "--rm", "--name", name,
                "-v", f"{pdir}:/in:ro", "-v", f"{seg_dir}:/out", "-v", f"{seg_dir / 'tmp'}:/tmpd", IMAGE,
                "--t1", f"/in/{study}_T1W.nii.gz", "--flair", f"/in/{study}_FLAIR.nii.gz",
                "--output", "/out", "--temp", "/tmpd", "--device", "cpu", "--segment_only",
                "--threads", threads, "--fast-mode"], capture_output=True, text=True, check=False)
            subprocess.run(["docker", "run", "--rm", "-v", f"{seg_dir}:/p", "--entrypoint", "chown", IMAGE,
                            "-R", f"{os.getuid()}:{os.getgid()}", "/p"], capture_output=True, check=False)
            (seg_dir / "lst-ai.log").write_text(res.stdout + "\n" + res.stderr)
            seg = seg_dir / SEG
            if res.returncode != 0 or not seg.exists():
                print(f"{pdir.name} {study}: LST-AI FAILED (exit {res.returncode}), see {seg_dir / 'lst-ai.log'}", flush=True)
                continue
            s = nib.load(seg)
            if s.shape != brain.shape or not np.allclose(s.affine, brain_img.affine, atol=1e-3):
                print(f"{pdir.name} {study}: LST-AI mask is on another grid, skipped", flush=True)
                continue
            m = (np.asarray(s.dataobj) > 0) & brain
            nib.save(nib.Nifti1Image(m.astype(np.uint8), brain_img.affine), mask)
            print(f"{pdir.name} {study}: {int(m.sum())} lesion voxels, {time.time() - t0:.0f} s", flush=True)
        rows.append((pdir.name, study, f"{t:.4f}", paths["T1w"], paths["FLAIR"], mask, out / "brainmask.nii.gz"))

complete = {s for s in {r[0] for r in rows} if sum(1 for r in rows if r[0] == s and Path(r[5]).exists()) == 2}
with open(work / "manifest_lesjak.tsv", "w") as fh:
    fh.write("subject\tsession\ttime_years\tt1\tflair\tmask\tbrainmask\n")
    for r in rows:
        if r[0] in complete:
            fh.write("\t".join(map(str, r)) + "\n")
print(f"{len(complete)} patients with both masks -> {work / 'manifest_lesjak.tsv'}")
