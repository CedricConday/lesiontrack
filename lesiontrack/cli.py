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
from .config import Params, RegParams, SELParams
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
    """Generate synthetic follow-ups, run the pipeline on each, score recovery, aggregate."""
    from .pipeline import run_subject
    from .registration import run_greedy
    from .scores import cohort_score
    from .synth import make_followup, score

    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    t1_img = nib.load(a.t1)
    bm = out / "baseline_brainmask.nii.gz"
    nib.save(nib.Nifti1Image((np.asarray(t1_img.dataobj) > 0).astype(np.uint8), t1_img.affine), bm)
    factors = tuple(float(x) for x in a.factors.split(","))
    seeds = [int(x) for x in a.seeds.split(",")]
    params = Params(reg=RegParams(threads=a.threads), sel=SELParams(min_mean_pct_per_year=a.min_mean))
    tables, all_metrics = [], []
    fractions = [float(x) for x in a.timepoints.split(",")]
    if fractions[-1] != 1.0:
        raise SystemExit("--timepoints must end with 1.0 (the last follow-up carries the full expansion)")
    for seed in seeds:
        sdir = out / f"seed{seed}"
        tps = [Timepoint("baseline", Path(a.t1), Path(a.flair), Path(a.mask), 0.0, bm)]
        for frac in fractions:
            synth_dir = sdir / f"synthetic_{frac:g}"
            rec = make_followup(Path(a.t1), Path(a.flair), Path(a.mask), synth_dir, n_expand=a.n_expand,
                                volume_factors=factors, seed=seed, time_fraction=frac)
            tps.append(Timepoint(f"synthetic_{frac:g}", synth_dir / "synth_T1w.nii.gz",
                                 synth_dir / "synth_FLAIR.nii.gz", synth_dir / "synth_mask.nii.gz",
                                 frac * a.dt_years))
        res = run_subject("backtest", tps, sdir / "run", params)
        pair = res.pairs[-1]
        resampled = sdir / "run" / "truth_lesion_labels_halfway.nii.gz"
        run_greedy(f"-d 3 -rf {pair.baseline_t1} -ri NN -rt int "
                   f"-rm {synth_dir / 'baseline_lesion_labels.nii.gz'} {resampled} -r {pair.halfway_matrix},-1")
        truth_labels = np.asarray(nib.load(resampled).dataobj).astype(np.int32)
        exp = np.asarray(nib.load(sdir / "run" / f"expansion_{pair.follow_up}_pct_per_year.nii.gz").dataobj)
        cand = np.asarray(nib.load(sdir / "run" / "sel_candidates.nii.gz").dataobj)
        # Definite SELs: Elliott's cohort z-scoring applied to this run's candidates.
        scored = cohort_score(res.candidates) if len(res.candidates) else res.candidates
        definite_ids = set(scored.loc[scored.get("definite_sel", pd.Series(dtype=bool)) == True, "candidate_id"]) if len(scored) else set()
        definite = np.isin(cand, list(definite_ids)) if definite_ids else np.zeros_like(cand, dtype=bool)
        table, metrics = score(rec, truth_labels, exp, cand, a.dt_years, definite_labels=definite)
        table.insert(0, "seed", seed)
        table.to_csv(sdir / "backtest_lesions.tsv", sep="\t", index=False, float_format="%.5g")
        metrics["seed"] = seed
        with open(sdir / "backtest_metrics.json", "w") as fh:
            json.dump(metrics, fh, indent=2)
        tables.append(table)
        all_metrics.append(metrics)
        print(f"seed {seed}: candidates sens {metrics['sensitivity']:.2f} FPR {metrics['false_positive_rate']:.2f} | "
              f"definite sens {metrics['sensitivity_definite']:.2f} FPR {metrics['false_positive_rate_definite']:.2f} | "
              f"recovery {metrics['recovery_fraction']:.2f}  noise p95 {metrics['noise_peak_p95_pct_per_year']:.1f} %/yr")

    combined = pd.concat(tables, ignore_index=True)
    combined.to_csv(out / "backtest_lesions.tsv", sep="\t", index=False, float_format="%.5g")
    _, agg = score_aggregate(combined)
    agg["candidate_mean_auc_per_seed"] = [m["candidate_mean_auc"] for m in all_metrics]
    agg["noise_candidates_per_1000_untouched_voxels_per_seed"] = [
        m["noise_candidates_per_1000_untouched_voxels"] for m in all_metrics]
    thrs = list(all_metrics[0]["min_mean_tradeoff"].keys()) if all_metrics and all_metrics[0]["min_mean_tradeoff"] else []
    # Averaged over seeds (each seed has its own lesion set), not pooled at the lesion level.
    agg["min_mean_tradeoff_mean_over_seeds"] = {
        thr: {
            "sensitivity": float(np.nanmean([m["min_mean_tradeoff"][thr]["sensitivity"] for m in all_metrics])),
            "false_positive_rate": float(np.nanmean([m["min_mean_tradeoff"][thr]["false_positive_rate"] for m in all_metrics])),
        }
        for thr in thrs
    }
    agg.update({"seeds": seeds, "dt_years": a.dt_years, "factors": list(factors), "n_expand": a.n_expand,
                "per_seed": all_metrics, "registration": params.reg.__dict__, "sel": params.sel.__dict__})
    with open(out / "backtest_metrics.json", "w") as fh:
        json.dump(agg, fh, indent=2)
    print(json.dumps({k: v for k, v in agg.items() if k not in ("per_seed", "registration", "sel")}, indent=2))
    return 0


def score_aggregate(table: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """Pool per-lesion rows from several seeds into one set of gate metrics."""
    from .scores import huber_fit

    exp = table[table["expanded"]]
    unexp = table[~table["expanded"]]
    slope = float("nan")
    if len(exp) >= 2:
        slope, _ = huber_fit(exp["true_pct_per_year"].to_numpy(), exp["measured_mean_pct_per_year"].to_numpy(), through_origin=True)
    agg = {
        "n_expanded": len(exp), "n_unexpanded": len(unexp),
        "sensitivity": float(exp["detected"].mean()) if len(exp) else float("nan"),
        "false_positive_rate": float(unexp["detected"].mean()) if len(unexp) else float("nan"),
        "sensitivity_definite": float(exp["detected_definite"].mean()) if len(exp) else float("nan"),
        "false_positive_rate_definite": float(unexp["detected_definite"].mean()) if len(unexp) else float("nan"),
        "recovery_fraction": slope,
        "recovery_median": float(exp["recovery"].median()) if len(exp) else float("nan"),
        "untouched_peak_median_pct_per_year": float(unexp["measured_peak_pct_per_year"].median()) if len(unexp) else float("nan"),
        "by_stratum": {
            s: {
                "n_expanded": int(g["expanded"].sum()),
                "sensitivity": float(g[g["expanded"]]["detected"].mean()) if g["expanded"].any() else float("nan"),
                "n_unexpanded": int((~g["expanded"]).sum()),
                "false_positive_rate": float(g[~g["expanded"]]["detected"].mean()) if (~g["expanded"]).any() else float("nan"),
                "recovery_median": float(g[g["expanded"]]["recovery"].median()) if g["expanded"].any() else float("nan"),
            }
            for s, g in table.groupby("size_stratum")
        },
    }
    return table, agg


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
    b.add_argument("--seeds", default="0,1,2", help="comma-separated seeds, one synthetic series each")
    b.add_argument("--timepoints", default="0.5,1.0",
                   help="comma-separated fractions of --dt-years at which synthetic follow-ups are made; must end in 1.0")
    b.add_argument("--factors", default="1.15,1.25,1.40", help="comma-separated volume factors to inject")
    b.add_argument("--threads", type=int, default=4)
    b.add_argument("--min-mean", type=float, default=0.0,
                   help="drop candidates whose mean expansion is below this (%%/yr); 0 = Elliott's rule")
    b.set_defaults(func=cmd_backtest)

    a = ap.parse_args(argv)
    return a.func(a)


if __name__ == "__main__":
    sys.exit(main())
