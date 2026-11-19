#!/usr/bin/env python
"""Segment a finished backtest with LST-AI so segmentation-based SEL detection can be scored.

The backtest tracks lesions on masks warped with the injected truth, which is a segmenter
that never disagrees with itself. This script adds the realistic arm: LST-AI v2 (docker,
CPU, fast mode) segments the baseline and each seed's last synthetic follow-up from the
images alone, and both masks are resliced into the halfway space the backtest used.
Usage: bench_lstai.py <backtest_dir> <baseline_t1> <baseline_flair> [threads]
Writes <backtest_dir>/seedN/run/reg/lstai_{baseline,followup}_mask_halfway.nii.gz
"""
import os
import subprocess
import sys
import time
from pathlib import Path

import nibabel as nib
import numpy as np

from lesiontrack.registration import run_greedy

IMAGE = "jqmcginnis/lst-ai:v2.0.0rc1-cpu"
SEG = "space-flair_seg-lst.nii.gz"


def lst_ai(t1: Path, flair: Path, out_dir: Path, threads: str) -> Path | None:
    """Binary LST-AI mask on the T1's grid; cached. None if the container failed."""
    mask = out_dir / "lstai_mask.nii.gz"
    if mask.exists():
        return mask
    (out_dir / "tmp").mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    res = subprocess.run([
        "docker", "run", "--rm", "--name", f"bench-lstai-{os.getpid()}",
        "-v", f"{t1.parent.resolve()}:/in1:ro", "-v", f"{flair.parent.resolve()}:/in2:ro",
        "-v", f"{out_dir.resolve()}:/out", "-v", f"{(out_dir / 'tmp').resolve()}:/tmpd", IMAGE,
        "--t1", f"/in1/{t1.name}", "--flair", f"/in2/{flair.name}", "--output", "/out", "--temp", "/tmpd",
        "--device", "cpu", "--segment_only", "--threads", threads, "--fast-mode"],
        capture_output=True, text=True, check=False)
    subprocess.run(["docker", "run", "--rm", "-v", f"{out_dir.resolve()}:/p", "--entrypoint", "chown", IMAGE,
                    "-R", f"{os.getuid()}:{os.getgid()}", "/p"], capture_output=True, check=False)
    (out_dir / "lst-ai.log").write_text(res.stdout + "\n" + res.stderr)
    seg = out_dir / SEG
    if res.returncode != 0 or not seg.exists():
        print(f"{out_dir}: LST-AI FAILED (exit {res.returncode}), see {out_dir / 'lst-ai.log'}", flush=True)
        return None
    s, ref = nib.load(seg), nib.load(t1)
    if s.shape != ref.shape or not np.allclose(s.affine, ref.affine, atol=1e-3):
        print(f"{out_dir}: LST-AI mask is on another grid than the T1, skipped", flush=True)
        return None
    m = (np.asarray(s.dataobj) > 0).astype(np.uint8)
    nib.save(nib.Nifti1Image(m, ref.affine), mask)
    print(f"{out_dir}: {int(m.sum())} lesion voxels, {time.time() - t0:.0f} s", flush=True)
    return mask


def main() -> int:
    bt, t1, flair = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
    threads = sys.argv[4] if len(sys.argv) > 4 else "4"
    base = lst_ai(t1, flair, bt / "lstai_baseline", threads)
    if base is None:
        return 1
    failed = 0
    for sdir in sorted(bt.glob("seed*")):
        reg = sdir / "run" / "reg"
        tags = sorted(p.name[: -len("_halfway.mat")] for p in reg.glob("synthetic_*_halfway.mat"))
        if not tags:
            print(f"{sdir}: no finished registration, skipped", flush=True)
            continue
        tag = tags[-1]  # the last follow-up carries the full expansion and defines the frame
        synth = sdir / tag
        fu = lst_ai(synth / "synth_T1w.nii.gz", synth / "synth_FLAIR.nii.gz", sdir / f"lstai_{tag}", threads)
        if fu is None:
            failed += 1
            continue
        ref, half = reg / f"{tag}_baseline_T1w_halfway.nii.gz", reg / f"{tag}_halfway.mat"
        # Same transforms as the pipeline: baseline -> frame is H^-1, follow-up -> frame is H.
        run_greedy(f"-d 3 -rf {ref} -ri NN -rm {base} {reg / 'lstai_baseline_mask_halfway.nii.gz'} -r {half},-1")
        run_greedy(f"-d 3 -rf {ref} -ri NN -rm {fu} {reg / 'lstai_followup_mask_halfway.nii.gz'} -r {half}")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
