# Validation record

Every measured result, by date and version, oldest first. Numbers from an earlier version
are kept as they were reported, with a note where a later version replaced them; the
README carries only the current ones. Working notes with the registration sweeps are in
`BACKTEST_NOTES.md`. Paths are relative to the repository root.

Versions: 0.1.0 took the Jacobian from greedy's `-rj` (wrong on RAS/LAS images); 0.2.0
computes it itself; 0.3.0 saves the warp unrounded, fits a global scale in the linear step
and adds the mean-rate filter.

## Backtest on MSLesSeg patient P1 (2026-09-26)

Three synthetic series (seeds 0, 1, 2), each with scans at 0, 0.5 and 1 year, eight of
eighteen lesions expanded by nominal factors 1.15 / 1.25 / 1.40 (analytic mean inside the
lesion 11 to 25 % because elongated lesions extend past the radial window), rigid
perturbation about 1 degree and 1 voxel per scan, intensity rescale 0.9 to 1.1, 2 % noise.
Default registration (greedy stationary-velocity, NCC 3x3x3, halfway space). Numbers
re-measured 2026-09-27 with the corrected Jacobian on the same registrations (the
2026-09-26 values, computed through greedy's `-rj`, are kept in `docs/BACKTEST_NOTES.md`).

| quantity | value |
|---|---|
| recovery fraction (measured / true, robust slope) | **1.16** |
| recovery by lesion size: >= 500 / 100-500 / < 100 voxels (median) | 1.00 / 1.22 / 1.21 |
| candidate sensitivity, Elliott thresholds | 1.00 |
| candidate false positive rate, Elliott thresholds | 0.47 (1.00 on untouched lesions >= 500 voxels, 0.07 under 100) |
| definite SEL (S >= 0, z-scored within the run): sensitivity / FPR | 0.67 / 0.37 |
| noise floor: median peak per-voxel expansion in untouched lesions | 10.8 %/yr |
| candidate mean expansion as a classifier, AUC | 0.99 to 1.00 |

Three conclusions the tool would not have earned without the gate:

1. **The Jacobian recovers the injected expansion in full**, slightly over on small
   lesions (the radial window's falloff is counted in the mean). The earlier "under-reports
   by 3x" was greedy's `-rj` on LAS images, see "The Jacobian is ours" below.
2. **Elliott's candidate rule is noise-limited on confluent lesions.** A one-voxel peak
   over 12.5 %/yr, grown by hysteresis to 4 %/yr, makes a 10-voxel candidate; every
   untouched lesion over 500 voxels produced one. Concentricity and constancy, the
   published filters, stayed weak on six more patients in 0.2.0; with unrounded warps
   (0.3.0) S reaches AUC 0.76. See "Validation beyond P1" and "0.3.0" below.
3. **The candidate's mean rate is the filter that works** (AUC 0.99). It is exposed as
   `min_mean_pct_per_year`; since 0.3.0 it is on by default at 10 (`--min-mean 0` gives
   Elliott's rule). Trade-off on this registration (mean over
   the three seeds):

| min mean (%/yr) | sensitivity | false positive rate |
|---|---|---|
| 0 (Elliott) | 1.00 | 0.47 |
| 8 | 1.00 | 0.20 |
| 10 | 1.00 | 0.07 |
| 12.5 | 0.92 | 0.00 |

Full record with the registration sweep: `docs/BACKTEST_NOTES.md`.

## MSLesSeg cohort (24 patients, 2026-09-26, Jacobian re-measured 2026-09-27)

All longitudinal MSLesSeg patients, baseline against their last scan (median interval
1.5 years, range 0.16 to 6.0), default settings. One patient (P49) was refused by the
input gate: its follow-up T1 in the released dataset has 0.1 % nonzero voxels.

**Lesion tracking, all subjects, last follow-up** (baseline lesions 573, matched groups):
new 216, resolved 250, enlarging 127, shrinking 115, stable 36, trend 45. That churn is
real signal, not just segmentation noise: voxels the masks call new go from FLAIR z 0.73
at baseline to 1.49 at follow-up and only 6 % were already lesion-bright; resolved voxels
go from 1.45 to 0.87. But masks are segmented per timepoint in MSLesSeg, so "new" and
"resolved" also carry every difference between two segmentations.

**SEL candidates.** 484 candidates in 24 subjects, median 16 per subject (12.5 with mean
>= 8 %/yr), 97 definite; the literature reports about 2 to 4 SELs per patient over two
years. On this registration the Elliott rule is noise-dominated on real data as it was on
the backtest, and `min_mean_pct_per_year` around 10 is where the backtest puts the floor.

**The Jacobian and the masks agree in rank, not in size.** Per subject, mask-based total
lesion volume change and the Jacobian integrated over baseline lesions correlate at
Spearman 0.41 (medians +12.5 % by mask, -6.3 % by Jacobian, brain-wide mean -1.5 %). Per
1:1 matched lesion (n = 276) the correlation is 0.22; 0.25 for lesions over 100 voxels.
Masks are segmented per timepoint in MSLesSeg, so their change carries every difference
between two segmentations; the Jacobian sees only displacement of the surrounding tissue,
and its median inside persisting lesions is negative.

Cohort tables: `derivatives/mslesseg/cohort_summary.tsv`, `sel_cohort.tsv`,
`check_volume_change_by_lesion.tsv`, `check_flair_consistency.tsv`
(`scripts/summarize_cohort.py`, `scripts/cohort_checks.py`).

## Validation beyond P1 (2026-11-20)

_These are 0.2.0 numbers, taken from warps greedy had rounded to 0.1 voxel; the 0.3.0
section below has them re-measured._

Same backtest design as for P1 (three seeds, eight lesions expanded per seed, 24 injected
lesions per backtest), lesiontrack 0.2.0 with the corrected Jacobian. The six extra
MSLesSeg patients (P2 P3 P12 P13 P19 P22) were picked before any result was seen.
`scripts/run_validation.sh` runs everything; `scripts/summarize_validation.py` writes
`results/validation_backtests.tsv` and `results/validation_feature_auc.tsv`.

**Seven MSLesSeg patients, greedy** (range over patients, median in brackets):

| quantity | value |
|---|---|
| candidate sensitivity, Elliott thresholds | 1.00 in all seven |
| candidate false positive rate | 0.25 to 0.81 (0.46) |
| definite SEL (S >= 0): sensitivity / FPR | 0.58 to 0.92 (0.71) / 0.06 to 0.42 (0.18) |
| recovery fraction | 0.98 to 1.25 |
| noise floor, median peak in untouched lesions | 10.3 to 18.8 %/yr |

How well each candidate feature separates injected from spurious candidates (AUC,
candidate level, median and range over the seven patients):

| feature | AUC |
|---|---|
| mean expansion rate | **1.00** (0.67 to 1.00) |
| constancy residual (lower = more constant) | 0.78 (0.60 to 0.81) |
| S = z(concentricity) - z(constancy) | 0.64 (0.56 to 0.74) |
| concentricity | 0.58 (0.45 to 0.73) |

The P1 finding holds on more patients: the published filters add little over chance on
this design and the mean rate does the work. Two patients (P12 0.67, P19 0.74) are the
exception; in both, the expanded lesions break into many small candidates (78 and 53
injected candidates for 24 injected lesions) whose mean rate overlaps the noise.

**ANTs SyN against greedy, same synthetic series** (P1 / P2 / P3):

| quantity | greedy | ANTs |
|---|---|---|
| candidates per backtest | 154 / 55 / 61 | 33 / 29 / 27 |
| candidate sensitivity | 1.00 / 1.00 / 1.00 | 1.00 / 1.00 / 1.00 |
| candidate false positive rate | 0.47 / 0.81 / 0.25 | **0.20 / 0.24 / 0.03** |
| definite SEL sensitivity | 0.71 / 0.71 / 0.79 | 0.58 / 0.50 / 0.54 |
| definite SEL FPR | 0.40 / 0.38 / 0.06 | 0.10 / 0.05 / 0.00 |
| recovery fraction | 1.16 / 1.21 / 1.14 | 1.06 / 1.04 / 1.02 |
| noise floor (%/yr) | 10.8 / 16.0 / 10.3 | 10.0 / 9.7 / 6.5 |
| run time, 4 cores | about 12 min | about 5.5 h |

ANTs is smoother: a third or less of the spurious candidates, recovery closer to 1, a lower noise
floor. It also loses about a quarter of the definite calls on truly expanding lesions, so S >= 0
gets more specific and less sensitive. At about 27 times the run time it is the engine for
a final analysis, not for a sweep. (P1 greedy here is the 2026-11-19 re-score, which
scores every candidate; its definite numbers differ slightly from the table above.)

**Lesjak et al. 2016, 20 patients** (two timepoints, 0.72 x 0.72 x 3 mm, expert masks of
lesion *change* only). Lesion masks come from LST-AI (`jqmcginnis/lst-ai:v2.0.0rc1-cpu`,
fast mode) on each timepoint; `scripts/prep_lesjak.py` prepares the inputs and
`scripts/eval_lesjak.py` scores against the expert change. There is no SEL ground truth
in this dataset, so it tests tracking, not SEL detection.

| quantity | value |
|---|---|
| LST-AI change voxels against expert change: Dice pooled / median per patient | 0.38 / 0.33 |
| expert change components touched: all / 1-10 / 10-50 / > 50 voxels | 453 of 846 / 52 of 185 / 219 of 434 / 182 of 227 |
| tracked groups called changed (new, resolved, enlarging, shrinking): sensitivity / specificity / PPV | 0.95 / 0.08 / 0.22 |
| Jacobian per baseline lesion (>= 10 voxels) as a classifier of expert change, AUC | 0.58 |
| Jacobian against mask change rate, matched groups, Spearman | 0.31 (n = 594) |
| SEL candidates: total / median per patient / on expert change | 318 / 9 / 121 |

Tracking on per-timepoint segmentations calls almost everything changed: 1,426 of 2,020
groups are new or resolved, against 282 expert-confirmed. The churn is the segmenter's,
not the tracker's, but the tool reports it as change, and that is what a user gets.
Backtests on three Lesjak patients (03, 06, 07, greedy) at this resolution: candidate
sensitivity 1.00 / 1.00 / 0.96, FPR 0.22 / 0.09 / 0.12, recovery 1.02 / 1.03 / 1.00,
noise floor 5.0 / 2.1 / 4.8 %/yr. Tables: `results/lesjak/`.

## 0.3.0: unrounded warps, similarity step, mean-rate filter on (2026-11-21)

Four changes, and every greedy number above re-run with them (`TAG=v031
scripts/run_v030.sh`, tables in `results/v031/`, cohorts in `derivatives/v031/`). ANTs and
the LST-AI arm were not re-run; their numbers above are 0.2.0, and the first change does
not touch ANTs.

1. **The warp is saved unrounded.** greedy rounds the warp it writes to 0.1 voxel unless
   told otherwise (`-wp`, default 0.1), and lesiontrack takes the Jacobian from that file.
   Every rate up to 0.2.0 sat on steps of about 5 %/yr at a one-year interval, coarser
   than Elliott's 4 %/yr growth threshold. This is the largest change below.
2. **The linear step fits one global scale** (`affine_dof = 7`, `--affine-dof 6` restores
   rigid). A simulated scanner change, every follow-up rescaled by 1 %, made untouched
   lesions look like they were expanding. The cost: whole-brain atrophy is absorbed as
   well, so lesion change is measured relative to the brain.
3. **`lesiontrack run` drops candidates under 10 %/yr mean rate** by default. `backtest`
   still scores every candidate.
4. **Earlier follow-ups land in the right frame.** For a third or later timepoint the
   rigid and halfway matrices were composed in the wrong order. At the small rotations
   seen here the two orders give the same alignment; it is fixed for larger ones.

`summary.json` now reports the residual global scale away from lesions (QC only) and the
halfway volume factor; tracking and candidate volumes are on the halfway grid, the
baseline lesion volume is native.

**P1 to P3, three seeds each, lesion level** (72 expanded, 150 untouched lesions per
condition), sensitivity / false positive rate, 0.2.0 -> 0.3.0:

| condition | Elliott candidate | definite S >= 0 | mean rate >= 10 |
|---|---|---|---|
| no scanner change | 1.00 / 0.37 -> 1.00 / **0.15** | 0.74 / 0.17 -> 0.57 / **0.04** | 0.99 / 0.03 -> 0.99 / 0.03 |
| 0.5 % scale change | 1.00 / 0.49 -> 1.00 / 0.12 | 0.74 / 0.19 -> 0.54 / 0.03 | 0.99 / 0.07 -> 0.99 / 0.03 |
| 1 % scale change | 1.00 / 0.67 -> 1.00 / 0.14 | 0.79 / 0.29 -> 0.54 / 0.05 | 1.00 / **0.18** -> 0.99 / **0.03** |
| FLAIR at 3 mm slices | 0.99 / 0.43 -> 0.99 / 0.23 | 0.64 / 0.18 -> 0.58 / 0.07 | 0.96 / 0.06 -> 0.96 / 0.05 |
| Lesjak 03, 06, 07 | 0.99 / 0.12 -> 0.97 / 0.08 | 0.61 / 0.05 -> 0.61 / 0.03 | 0.96 / 0.03 -> 0.96 / 0.03 |

A 1 % scale change now costs nothing. The Elliott rule and S >= 0 lose most of their
false positives; S >= 0 also loses sensitivity (0.74 -> 0.57), so it is now a specific
but conservative call. The mean rate was already near its ceiling and stays there.

**Seven MSLesSeg patients, greedy, 0.2.0 -> 0.3.0** (median, range):

| quantity | 0.2.0 | 0.3.0 |
|---|---|---|
| candidates per backtest | 154 (55 to 278) | 37 (29 to 152) |
| candidate FPR | 0.46 (0.25 to 0.81) | 0.27 (0.04 to 0.59) |
| definite SEL sensitivity / FPR | 0.71 / 0.18 | 0.71 / 0.10 |
| recovery fraction | 1.16 | 1.17 |
| noise floor, median peak in untouched lesions | 10.8 %/yr | 8.9 %/yr (5.1 to 14.4) |
| AUC mean rate | 1.00 (0.67 to 1.00) | 0.99 (0.87 to 1.00) |
| AUC constancy residual | 0.78 | 0.76 (0.61 to 0.94) |
| AUC S | 0.64 | **0.76** (0.63 to 0.78) |
| AUC concentricity | 0.58 | 0.62 (0.31 to 0.76) |

With unrounded rates S carries real information (0.76), though the mean rate still
separates best. P13 stays the hardest patient (candidate FPR 0.59).

**Greedy 0.3.0 against ANTs 0.2.0** (P1 / P2 / P3): candidates 37 / 35 / 29 against 33 /
29 / 27; candidate FPR 0.27 / 0.48 / 0.04 against 0.20 / 0.24 / 0.03; definite sensitivity
0.50 / 0.50 / 0.71 against 0.58 / 0.50 / 0.54; noise floor 8.9 / 12.3 / 5.1 against 10.0 /
9.7 / 6.5 %/yr. Most of the gap in the ANTs table above was greedy's rounding; ANTs
keeps an edge on P2 and in recovery (1.02 to 1.06 against 1.14 to 1.21).

**Cohorts with the defaults.** MSLesSeg (24 patients): 232 candidates, median 6 per
patient (P2 none), 68 definite; 0.2.0 without the filter gave 484, 16 and 97. Mask against
Jacobian volume change per subject: Spearman 0.40 (0.41). The halfway volume factor spans
0.975 to 1.010, which is how much global scale the linear step took up. Lesjak (20
patients): 208 candidates (median 6.5), 45 definite, 85 on expert change (318 / 62 / 121);
tracking is unchanged, it depends on the masks alone.


## Slow expansion and the mean-rate floor (2026-11-23, 0.3.0)

Raised in review: the backtests above inject 11 to 25 %/yr inside the lesion and the
10 %/yr floor was chosen on them, while SELs are defined by slow expansion.
`scripts/run_slow.sh` injects nominal 3, 5 and 8 %/yr, three scans, over two years and
over one, on P1 to P3; `scripts/summarize_slow.py` writes `results/slow/`. All six
backtests are in; the table and what it says are in the README ("Backtest: slow
expansion"), the numbers in `results/slow/sensitivity_by_rate.tsv` and `conditions.tsv`.
Recovery fraction at these rates is 1.41 (two years) and 1.67 (one year): small
expansions are over-read. Definite SEL (S >= 0) sensitivity at 3 / 5 / 8 %/yr: 0.07 /
0.11 / 0.61 at two years, 0.26 / 0.30 / 0.67 at one; false positive rate 0.00 and 0.07.

The floor became one on total change (10 % over the interval) on the two-year result.

**Cohorts re-run with that floor and the FLAIR flag** (`scripts/run_cohorts.sh`, cached
registrations, tables in `results/v0.3.0/`). MSLesSeg: 374 candidates, median 12.5 per
patient, 92 definite (0.3.0.dev0 with the fixed floor: 232 / 6 / 68); most intervals are
over a year, so the floor is lower than 10 %/yr for them. Lesjak: 162 candidates, median
7, 46 definite, 82 on expert change (dev0: 208 / 6.5 / 45 / 85); its median interval is
0.87 years, so the floor is higher. Tracking with FLAIR support on Lesjak: sensitivity
0.95 -> 0.81, specificity 0.08 -> 0.59, PPV 0.22 -> 0.35 (threshold chosen in-sample).
**ANTs at 0.3.0 against greedy at 0.3.0 (2026-11-24).** The same three synthetic series,
both engines on lesiontrack 0.3.0 with unrounded warps, 4 threads, three seeds each
(`scripts/run_ants030.sh 4`; metrics in `results/v0.3.0/ants/`). P1 / P2 / P3:

| quantity | greedy 0.3.0 | ANTs 0.3.0 |
|---|---|---|
| candidates per backtest | 37 / 35 / 29 | 36 / 29 / 28 |
| candidate sensitivity | 1.00 / 1.00 / 1.00 | 1.00 / 1.00 / 1.00 |
| candidate false positive rate | 0.27 / 0.48 / 0.04 | 0.23 / 0.24 / 0.03 |
| definite SEL sensitivity | 0.50 / 0.50 / 0.71 | 0.58 / 0.58 / 0.54 |
| definite SEL FPR | 0.10 / 0.14 / 0.00 | 0.07 / 0.05 / 0.00 |
| recovery fraction | 1.17 / 1.21 / 1.14 | **1.07 / 1.05 / 1.03** |
| noise floor (%/yr) | 8.9 / 12.3 / 5.2 | 10.4 / 9.9 / 6.9 |
| AUC mean rate | 1.00 / 0.99 / 0.99 | 0.96 / 1.00 / 1.00 |
| AUC S | 0.67 / 0.63 / 0.76 | 0.67 / 0.88 / 0.91 |
| run time per patient, 4 cores | about 12 min | about 5.2 h |

With greedy's rounding fixed, most of the gap closes. One difference is consistent:
recovery. Greedy over-reads by 14 to 21 % (per seed 1.08 to 1.23), ANTs by 3 to 7 % (0.95
to 1.10). The rest sits inside the seed-to-seed spread and should not be read as an
engine effect. P2's candidate FPR looks halved (0.48 against 0.24), but the per-seed
values are 0.14 / 0.71 / 0.57 against 0.43 / 0.29 / 0.00. S separates better on ANTs for P2
and P3 (0.88 and 0.91 against 0.63 and 0.76). The noise floor is not lower with ANTs:
it is higher on P1 and P3 and lower on P2. ANTs at 0.3.0 differs from ANTs at 0.2.0 by at
most 3 candidates and 0.08 in any rate, so the 0.3.0 changes did not reach the ANTs path,
as expected. Greedy stays the default. At about 26 times the run time ANTs is for a
final analysis, not for a sweep. Run time is from P2 and P3 (18,680 s and 18,784 s). P1 is
left out: its wall time includes the run being paused for most of a day. Three patients
with three seeds each is a small sample.

**Run-to-run reproducibility (2026-11-23).** Lesjak patient03 was registered three times
with identical settings: the v031 run (4 threads), once at 1 thread, once more at 4.
Neither repeat reproduced the first. Candidates shared by the two 4-thread runs: 131 /
184 / 24 / 39 voxels against 129 / 177 / 26 / 37, mean rates within 1 %/yr (10.05 against
9.93 for the weakest); the repeat also has a 10-voxel candidate at the size limit that the
first run lacks. So greedy's multi-threaded result is not bit-reproducible and candidates
near a threshold can flip. The registration cache keys on the thread count, so a cached
tree is only reused at the thread count it was made with.

## 0.3.2: lesion masks resliced without the similarity scale (2026-11-25)

Since 0.3.0 the linear step is a similarity. Nearest-neighbour reslicing of the lesion masks under its
scale (a fraction of a percent from 1) duplicated one voxel layer per axis at follow-up and dropped one at
baseline, so a lesion cut by such a plane read as enlarging. mscard's calibration on 24 MSLesSeg subjects
found it: 14 of 593 unchanged lesions called enlarging (+9 to +39 %/yr), against 1 of 594 under 0.2.0. On
120 lesions of four of those subjects, the same lesion pushed through the halfway transform and its inverse
differed by 8 % or more for 14; through the rotation and translation alone, for none (largest 0.2 %).
0.3.2 reslices masks that way. Lesion volumes are again each timepoint's own mm³, and a mask sits off its
image by at most |1 − scale| times its distance from the grid centre (about 0.25 mm at 100 mm).

After the fix, mscard's calibration: 2 of 597 unchanged lesions called enlarging or shrinking (smooth
profile 1 of 597). Four lesions of 5 to 27 mm³ that do shrink in the truth are now called resolved: their
follow-up remnant sits one voxel from the baseline lesion and no longer overlaps it.

The greedy backtests for P1 to P3, re-run on the cached registrations (`scripts/run_v032.sh`,
`results/v032/`), are unchanged except P2's definite-SEL figures:

| | 0.3.1 (P1/P2/P3) | 0.3.2 (P1/P2/P3) |
|---|---|---|
| FPR | 0.27/0.48/0.04 | 0.27/0.48/0.04 |
| Definite SEL sensitivity | 0.50/0.50/0.71 | 0.50/0.54/0.71 |
| Definite SEL FPR | 0.10/0.14/0.00 | 0.10/0.19/0.00 |
| Recovery | 1.17/1.21/1.14 | 1.17/1.21/1.14 |
| Noise floor (%/yr) | 8.9/12.3/5.2 | 8.9/12.3/5.2 |

The P2 change is inside that patient's per-seed spread (definite SEL FPR 0.14/0.71/0.57 across seeds at
0.3.0); part of P2's registrations were redone because the cached frame path named the old directory, and a
re-run of the same registration is not bit-identical.
