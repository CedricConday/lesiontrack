"""Command line: run a manifest, generate a synthetic follow-up, score a cohort."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd

from . import __version__
from .config import Params, RegParams
from .registration import Timepoint

MANIFEST_COLUMNS = ("subject", "session", "time_years", "t1", "flair", "mask")


def read_manifest(path: Path) -> dict[str, list[Timepoint]]:
    """TSV with columns subject, session, time_years, t1, flair, mask[, brainmask]."""
    df = pd.read_csv(path, sep="\t", dtype=str)
    missing = [c for c in MANIFEST_COLUMNS if c not in df.columns]
    if missing:
        raise SystemExit(f"manifest lacks columns: {missing}")
    base = path.resolve().parent

    def p(v) -> Path:
        q = Path(v)
        return q if q.is_absolute() else base / q

    subjects: dict[str, list[Timepoint]] = {}
    for _, r in df.iterrows():
        bm = r.get("brainmask")
        bm = p(bm) if isinstance(bm, str) and bm else None
        subjects.setdefault(r["subject"], []).append(
            Timepoint(r["session"], p(r["t1"]), p(r["flair"]), p(r["mask"]), float(r["time_years"]), bm)
        )
    return subjects


def cmd_run(a: argparse.Namespace) -> int:
    from .pipeline import run_subject
    from .report import overview

    subjects = read_manifest(Path(a.manifest))
    if a.subject:
        subjects = {k: v for k, v in subjects.items() if k in set(a.subject)}
        if not subjects:
            raise SystemExit("no matching subject in manifest")
    params = Params(reg=RegParams(threads=a.threads))
    out = Path(a.out)
    for name, tps in subjects.items():
        res = run_subject(name, tps, out / name, params)
        for p in res.pairs:
            overview(out / name, p.follow_up)
        print(f"{name}: {len(res.candidates)} SEL candidates, "
              f"{len(res.tracking)} lesion groups over {len(res.pairs)} follow-up(s) -> {out / name}")
    return 0


def cmd_cohort(a: argparse.Namespace) -> int:
    from .scores import cohort_score

    out = Path(a.out)
    files = sorted(out.glob("*/sel_candidates.tsv"))
    if not files:
        raise SystemExit(f"no */sel_candidates.tsv under {out}")
    table = pd.concat([pd.read_csv(f, sep="\t") for f in files], ignore_index=True)
    scored = cohort_score(table)
    scored.to_csv(out / "sel_cohort.tsv", sep="\t", index=False, float_format="%.5g")
    per_subject = scored.groupby("subject").agg(
        candidates=("candidate_id", "count"), definite=("definite_sel", "sum"),
        candidate_volume_mm3=("volume_mm3", "sum"),
    )
    per_subject.to_csv(out / "sel_cohort_by_subject.tsv", sep="\t", float_format="%.5g")
    print(per_subject.to_string())
    print(f"\n{len(scored)} candidates in {len(files)} subjects; "
          f"{int(scored['definite_sel'].sum())} definite (S >= 0) -> {out / 'sel_cohort.tsv'}")
    return 0


def cmd_synth(a: argparse.Namespace) -> int:
    from .synth import make_followup

    rec = make_followup(Path(a.t1), Path(a.flair), Path(a.mask), Path(a.out),
                        n_expand=a.n_expand, seed=a.seed)
    n = sum(1 for x in rec["lesions"] if x["volume_factor"] > 1)
    print(f"wrote synthetic follow-up with {n} expanded of {len(rec['lesions'])} lesions -> {a.out}/truth.json")
    return 0


def cmd_backtest(a: argparse.Namespace) -> int:
    """Generate a synthetic follow-up, run the pipeline on it, score recovery."""
    from .pipeline import run_subject
    from .synth import make_followup, score

    out = Path(a.out)
    synth_dir = out / "synthetic"
    rec = make_followup(Path(a.t1), Path(a.flair), Path(a.mask), synth_dir,
                        n_expand=a.n_expand, seed=a.seed)
    t1_img = nib.load(a.t1)
    bm = synth_dir / "baseline_brainmask.nii.gz"
    nib.save(nib.Nifti1Image((np.asarray(t1_img.dataobj) > 0).astype(np.uint8), t1_img.affine), bm)
    tps = [
        Timepoint("baseline", Path(a.t1), Path(a.flair), Path(a.mask), 0.0, bm),
        Timepoint("synthetic", synth_dir / "synth_T1w.nii.gz", synth_dir / "synth_FLAIR.nii.gz",
                  synth_dir / "synth_mask.nii.gz", a.dt_years),
    ]
    params = Params(reg=RegParams(threads=a.threads))
    res = run_subject("backtest", tps, out / "run", params)
    pair = res.pairs[0]
    labels = np.asarray(nib.load(out / "run" / "baseline_lesion_labels.nii.gz").dataobj)
    # lesion ids in the run are those of the halfway-resampled baseline mask; map back to
    # the truth's ids through the un-resampled labels resampled the same way.
    truth_labels = np.asarray(nib.load(synth_dir / "baseline_lesion_labels.nii.gz").dataobj)
    from .registration import run_greedy
    resampled = out / "run" / "truth_lesion_labels_halfway.nii.gz"
    run_greedy(f"-d 3 -rf {pair.baseline_t1} -ri NN -rt int "
               f"-rm {synth_dir / 'baseline_lesion_labels.nii.gz'} {resampled} -r {pair.halfway_matrix},-1")
    truth_labels = np.asarray(nib.load(resampled).dataobj).astype(np.int32)
    exp = np.asarray(nib.load(out / "run" / f"expansion_{pair.follow_up}_pct_per_year.nii.gz").dataobj)
    cand = np.asarray(nib.load(out / "run" / "sel_candidates.nii.gz").dataobj)
    table, metrics = score(rec, truth_labels, exp, cand, a.dt_years)
    table.to_csv(out / "backtest_lesions.tsv", sep="\t", index=False, float_format="%.5g")
    metrics["seed"] = a.seed
    metrics["dt_years"] = a.dt_years
    metrics["n_baseline_lesions_in_run"] = int(labels.max())
    with open(out / "backtest_metrics.json", "w") as fh:
        json.dump(metrics, fh, indent=2)
    print(json.dumps(metrics, indent=2))
    return 0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="lesiontrack", description=__doc__)
    ap.add_argument("--version", action="version", version=f"lesiontrack {__version__}")
    sub = ap.add_subparsers(dest="cmd", required=True)

    r = sub.add_parser("run", help="run every subject in a manifest")
    r.add_argument("manifest")
    r.add_argument("--out", required=True)
    r.add_argument("--subject", action="append", help="restrict to these subjects (repeatable)")
    r.add_argument("--threads", type=int, default=4)
    r.set_defaults(func=cmd_run)

    c = sub.add_parser("cohort", help="z-score candidates across all subjects in an output tree")
    c.add_argument("--out", required=True)
    c.set_defaults(func=cmd_cohort)

    s = sub.add_parser("synth", help="write a synthetic follow-up with known expansion")
    s.add_argument("--t1", required=True)
    s.add_argument("--flair", required=True)
    s.add_argument("--mask", required=True)
    s.add_argument("--out", required=True)
    s.add_argument("--n-expand", type=int, default=8)
    s.add_argument("--seed", type=int, default=0)
    s.set_defaults(func=cmd_synth)

    b = sub.add_parser("backtest", help="synthesize, run, and score recovery of known expansion")
    b.add_argument("--t1", required=True)
    b.add_argument("--flair", required=True)
    b.add_argument("--mask", required=True)
    b.add_argument("--out", required=True)
    b.add_argument("--n-expand", type=int, default=8)
    b.add_argument("--dt-years", type=float, default=1.0)
    b.add_argument("--seed", type=int, default=0)
    b.add_argument("--threads", type=int, default=4)
    b.set_defaults(func=cmd_backtest)

    a = ap.parse_args(argv)
    return a.func(a)


if __name__ == "__main__":
    sys.exit(main())
