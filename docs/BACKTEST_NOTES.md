# Backtest notes (dated, append only)

## 2026-09-26 — first gate results on MSLesSeg P1, seed 0

Synthetic follow-up: 8 of 18 lesions expanded by nominal factors 1.15 / 1.25 / 1.40 over
1 year, rigid perturbation about 1 degree / 1 voxel, intensity scale 0.9–1.1, 2 % noise.

**Bug the gate caught first.** greedy applies `-ri` to the *next* `-rm` pair. With the
flag placed after the pair, every mask was linearly interpolated: baseline mask +45 % as a
float halo, follow-up mask −36 % as truncated uint8. Every tracked lesion read as
"shrinking 50–90 %". Fixed in cfe9d08 with a regression test.

**Analytic truth.** The radial window scales a sphere around the centroid, so an elongated
lesion's true mean expansion is below the nominal factor (21.5 % inside a 1354-voxel
periventricular lesion given a nominal 40 %). The generator now writes the analytic
Jacobian of the injected field and scores against it (659a16d).

**Registration sweep on the same pair** (mean measured expansion inside injected lesions,
analytic truth ≈ 21 % for the large ones; "untouched max" = highest per-voxel value inside
any untouched lesion):

| setting | time | large injected mean | untouched max |
|---|---|---|---|
| default: NCC 2x2x2, 100x60x30, e 0.7, s 2vox 0.5vox | 161 s | ≈ 5 | 27 |
| s 3vox 1vox | 137 s | 2.2 | 21.5 |
| 200x120x60 | 260 s | 3.3 | 31 |
| NCC 4x4x4 | 136 s | 3.3 | 26.5 |
| e 1.0 | 136 s | 3.2 | 27.6 |
| SSD | 61 s | −9.4 | 27 (intensity rescale breaks it) |
| 100x100x300 | 1005 s | 9.6 | 27.8 |
| **-sv (stationary velocity)** | **167 s** | **9.5** | 27.3 |
| s 1vox 0vox, 100x100x100 | 319 s | 10.7 | 110 (unusable noise) |
| factor 2.0, default | 151 s | 21 of ≈ 78 analytic | 27.4 |

Stationary-velocity mode matches ten times the iterations at a sixth of the time. A third
sweep around it plateaued: `-sv` with 100x100x100 iterations 9.3 (497 s), with step 1.0
9.7, with NCC 3x3x3 9.9 (170 s). Default set to `-sv`, NCC 3x3x3, 100x60x30, step 0.7.

**Recovery is a fraction, and it depends on lesion size.** With `-sv`: recovery fraction
0.28 overall (robust slope of measured vs analytic true); median 0.56 for lesions ≥ 500
voxels, 0.37 for 100–500, 0.16 for < 100. A 15 % volume increase on a 45-voxel lesion
moves its boundary by 0.1 voxel; no 1 mm registration sees that.

**Where the false positives come from.** Noise candidates have one or two seed voxels at
about 15.5 %/yr (just over JE1 = 12.5) and sit almost entirely in the 15 633-voxel
confluent lesion: more voxels, more chances for a peak. Concentricity does not separate
them from injected candidates (medians 1.57 vs 1.88). Elliott's filter for exactly this is
constancy across intermediate timepoints, which a two-timepoint synthetic cannot exercise;
the backtest now generates three timepoints (0, 0.5, 1 year).

**Open:** whether constancy filters the noise candidates; whether a small minimum
seed-voxel count would be a defensible addition (it is not in Elliott 2019 and is not
enabled).

## 2026-09-26 — three seeds, three timepoints (0, 0.5, 1 yr), `-sv` NCC 3x3x3

Pooled over 24 injected and 30 untouched lesions: candidate sensitivity 0.88, false
positive rate 0.47 (1.0 on untouched lesions ≥ 500 voxels, 0.07 under 100 voxels).
Elliott's definite/possible split (cohort z-scored S ≥ 0): sensitivity 0.29, false
positive rate 0.27, i.e. no separation. Recovery fraction 0.30 (median 0.55 for ≥ 500
voxels, 0.27 for 100–500, 0.33 under 100).

**What separates injected from noise candidates** (139 candidates, 50 injected, rank AUC):

| feature | AUC | injected median | noise median |
|---|---|---|---|
| mean expansion in candidate | **0.90** | 8.6 %/yr | 6.3 %/yr |
| fitted expansion slope (constancy fit) | 0.90 | 8.3 | 5.8 |
| expansion at the half-year scan | 0.76 | 8.7 | 4.3 |
| seed voxels | 0.57 | 2 | 2 |
| concentricity | 0.54 | 2.6 | 2.4 |
| peak expansion | 0.45 | 15.5 | 15.8 |
| constancy residual | 0.40 | 0.01 | 0.01 |
| candidate size | 0.23 | 31 | 111 |

Noise candidates are hysteresis regions grown from a one- or two-voxel peak just over
JE1; 62 of 89 sit in lesions over 2000 voxels. The half-year expansion of an injected
candidate equals its full-year rate (ratio 1.0, linear growth); for noise it is 0.71.

Decision: `min_mean_pct_per_year` added to `SELParams`, default 0 (Elliott unchanged); the
backtest prints the sensitivity / false-positive trade-off at 0, 6, 8, 10, 12.5 %/yr.

**Trade-off table from the re-score (three seeds pooled, lesion level):**

| min mean (%/yr) | sensitivity | false positive rate |
|---|---|---|
| 0 (Elliott) | 0.88 | 0.47 |
| 8 | 0.67 | 0.23 |
| 10 | 0.50 | 0.00 |
| 12.5 | 0.21 | 0.00 |

Candidate-mean AUC per seed 0.96 / 0.85 / 0.92; noise candidates per 1000 untouched
lesion voxels 1.7 / 2.0 / 1.1. On this registration a user who wants no false positives
pays half the sensitivity; the injected rates (15–40 % nominal, 11–25 % analytic over one
year) are at the low end of what SEL studies report, so real SELs over two years would sit
higher on the curve.

## 2026-09-26 — re-run after the review fixes (one frame per subject, composed synthetic motion)

Same design, `work/backtest_P1_v4`. Recovery 0.29 (0.57 / 0.27 / 0.29 by size stratum),
candidate sensitivity 0.88, false positive rate 0.47, definite 0.42 / 0.30, mean-rate AUC
0.96 / 0.91 / 0.91. Trade-off (mean over seeds): min mean 8 %/yr gives 0.67 / 0.10, 10 %/yr
0.38 / 0.00. The composed motion tightened the 8 %/yr operating point (false positives 0.23
before, 0.10 now); everything else moved within seed noise. These are the numbers in the
README.

## 2026-09-26 — MSLesSeg cohort, 24 subjects, default settings

Numbers as in the README. Additional detail: subject-level mask vs Jacobian volume change,
worst disagreements: P50 mask +298 % / Jacobian −6 %, P49 excluded (empty follow-up T1),
P3 −36 % / −0.2 %, P31 +52 % / −5 %. Agreements exist (P20 +107 % / +26 %, P22 −8 % / −12 %)
but are the minority. Brain-mean Jacobian −1.1 % per interval (plausible atrophy plus bias).
Elliott's definite/possible split kept 116 of 530 candidates (S >= 0 after z-scoring; 22 %
because unscorable candidates are never definite).

Failure found and gated: P49's follow-up T1 (0.1 % nonzero) produced a 52-degree, 190 mm
rigid transform that the pipeline accepted. Inputs and rigid transforms are now checked.

## 2026-09-27 — greedy's `-rj` Jacobian ignores the direction matrix; recovery re-measured

Found while adding ANTs SyN as a second engine: on a 56^3 phantom with an identity affine
and one lesion expanded by 1.6, greedy read det 0.80 inside the lesion for every setting
tried (NCC 3x3x3 / 5x5x5, SSD, T1 only, FLAIR only, with and without `-gm`), ANTs read
1.48, the masks said 1.59. A sphere growing from radius 8 to 10 (ratio 1.95) under three
affines settled it: `-rj` gives 0.675 (RAS), 1.183 (LAS), 2.034 (LPS-aligned); `-jac`
gives 2.033 for all. `RunReslice` calls `field_jacobian_det` on the physical-space warp;
`RunJacobian` converts it to voxel units first. For a radial stretch 1+k the LAS value is
(1-k)(1+k)^2 against the true (1+k)^3.

MSLesSeg is LAS. Every Jacobian in the cohort results and in the backtests above was
computed this way. The "Jacobian under-reports by about 3x" conclusion of 2026-09-26 was
this bug, not a property of 1 mm registration. lesiontrack 0.2.0 computes the determinant
itself (`registration.jacobian_determinant`), for both engines; cached registrations get
their Jacobian recomputed from the stored warp without re-registering.

Re-measured numbers (same warps as `work/backtest_P1_v4`, only the Jacobian recomputed):

| quantity | 2026-09-26 (`-rj`) | 2026-09-27 (voxel-frame) |
|---|---|---|
| recovery fraction (robust slope) | 0.29 | 1.16 |
| recovery median, >= 500 / 100-500 / < 100 voxels | 0.57 / 0.27 / 0.29 | 1.00 / 1.22 / 1.21 |
| candidate sensitivity / FPR | 0.88 / 0.47 | 1.00 / 0.47 |
| definite sensitivity / FPR | 0.42 / 0.30 | 0.67 / 0.37 |
| candidate mean AUC per seed | 0.91 to 0.96 | 0.99, 1.00, 1.00 |
| min mean 8 / 10 / 12.5: sensitivity, FPR | 0.67, 0.10 / 0.38, 0.00 / - | 1.00, 0.20 / 1.00, 0.07 / 0.92, 0.00 |
| untouched peak median | 10.9 %/yr | 10.8 %/yr |

The noise floor is unchanged (it is set by registration error, which the bug scaled by
about the same factor as the signal), so the candidate false-positive picture stands; the
recovery story is gone. Small-lesion recovery above 1 is the radial window's falloff
region being counted in the lesion mean, not a bias in the Jacobian.

Cohort, same warps: Spearman(mask change, Jacobian change) per subject 0.10 -> 0.41,
per matched lesion 0.09 -> 0.22 (>= 100 voxels 0.18 -> 0.25); candidates 530 -> 484,
definite 97; median Jacobian change inside baseline lesions -6.3 % per subject against
+12.5 % by mask.

Second engine: ANTs SyN (antspyx, CC radius 4, identity initialisation on the halfway
images) is in as `--engine ants`; on the test phantom both engines recover the injected
sign (greedy 1.72, ANTs 1.48 for a true 1.6). No cohort or backtest run with it yet.

## 2026-11-20 — validation chain: six more patients, ANTs, Lesjak

`scripts/run_validation.sh 4`, launched 2026-11-19 18:32, `all done` 2026-11-20 19:47, every
step exit 0 (log `work/logs/validation.log`). lesiontrack 0.2.0 plus `5a17a11` (every candidate
scored). Patients P2 P3 P12 P13 P19 P22 fixed before any result was seen.

Per-backtest numbers: `results/validation_backtests.tsv`; candidate-feature AUCs:
`results/validation_feature_auc.tsv` (both from `scripts/summarize_validation.py`). P1 greedy
is `work/backtest_P1_v5`, the 2026-11-19 re-score.

* Greedy, 7 MSLesSeg patients: candidate sensitivity 1.00 everywhere; FPR 0.25-0.81
  (median 0.46); definite S>=0 0.58-0.92 / FPR 0.06-0.42; recovery 0.98-1.25; noise floor
  10.3-18.8 %/yr.
* Feature AUC (candidate level), median over the 7: mean rate 1.00, constancy residual 0.78,
  S 0.64, concentricity 0.58. Gap (1) closed: the published filters stay near chance on the
  corrected Jacobian. P12 (0.67) and P19 (0.74) are the mean-rate outliers: expanded lesions
  split into many small candidates (78 and 53 injected candidates for 24 lesions).
* ANTs vs greedy, P1/P2/P3: candidate FPR 0.47/0.81/0.25 -> 0.20/0.24/0.03; definite
  sensitivity 0.71/0.71/0.79 -> 0.58/0.50/0.54; recovery 1.16/1.21/1.14 -> 1.06/1.04/1.02;
  noise floor 10.8/16.0/10.3 -> 10.0/9.7/6.5. ANTs about 20,000 s per backtest at 4 cores
  against about 730 s for greedy. Gap (2) closed.
* Lesjak, 20 patients, LST-AI fast-mode masks, scored against expert change masks: change
  Dice pooled 0.38 (median 0.33); expert components touched 453/846 (small ones 52/185);
  tracked groups called changed sens 0.95 / spec 0.08 / PPV 0.22 (1,426 of 2,020 groups new or
  resolved, 282 expert-confirmed); Jacobian per baseline lesion vs expert change AUC 0.58;
  Spearman Jacobian vs mask rate 0.31 (594 groups). No SEL ground truth: tracking only.
  Lesjak backtests 03/06/07 (greedy): sensitivity 1.00/1.00/0.96, FPR 0.22/0.09/0.12,
  recovery 1.02/1.03/1.00, noise floor 5.0/2.1/4.8 %/yr. Gap (3) closed.

## 2026-11-21 — 0.3.0.dev0 re-run (similarity step, mean-rate default, frame fix)

`scripts/run_v030.sh 4`, 10:55 to 16:55, 81 steps exit 0; one exit 1, MSLesSeg P49, the
known broken T1 (input gate, as in 0.2.0). Code `b30520d`. Tables `results/v030/`, cohorts
`derivatives/v030/`. ANTs and the LST-AI arm not re-run.

* Pooled P1-P3, lesion level, mean rate >= 10 sens / FPR, 0.2.0 -> 0.3.0: unscaled
  0.99/0.03 -> 0.97/0.03; 0.5 % scale 0.99/0.07 -> 0.97/0.06; 1 % scale 1.00/0.18 ->
  0.99/0.00; thick3 0.96/0.06 -> 0.96/0.05; Lesjak x3 0.96/0.03 -> 0.96/0.05.
* Definite S>=0 sensitivity drops: 0.74 -> 0.62 unscaled, 0.79 -> 0.56 at 1 %. FPR level.
* Noise floor (median untouched peak) over 7 patients 10.8 -> 15.5 %/yr (P1 10.8 -> 15.8).
  Not explained. Candidates to test: the frame fix (changes the intermediate timepoint,
  which feeds the peak) against the similarity step (`work/dof7/` has dof 7 before the frame
  fix on P1). Open.
* Feature AUC medians: mean rate 0.98, constancy residual 0.75, S 0.62, concentricity 0.57.
* Cohorts with defaults: MSLesSeg 246 candidates (median 6.5, P2 and P20 none), 75
  definite; Lesjak 231 / 40; Lesjak tracking identical to 0.2.0. Halfway volume factor
  0.975-1.010; residual global scale in backtests 0.993-0.996.

## 2026-11-21 late — greedy's saved warp was rounded; v031 re-run

The v030 noise-floor rise (10.8 -> 15.5 %/yr) was not the similarity step (P1 with
`--affine-dof 6` on the new code: 15.8) and not the frame order (intermediate follow-up
resliced with either order: correlation 0.9902 both). Peaks sat on ~5 %/yr steps: greedy
writes the warp rounded to 0.1 voxel (`-wp`, default 0.1; the old warps have 0 % of values
off that grid). Fixed `16d72c0` (`-wp 0`, `RegParams.warp_precision_vox`), `3684dc2` (ANTs
records not invalidated). Affects every greedy number through v030; not ANTs.

`TAG=v031 scripts/run_v030.sh 4`, 17:39 to 23:20, 81 exit 0, P49 exit 1 (broken T1, known).
Tables `results/v031/`; `results/v030/` kept as the rounded-warp record.

* Pooled P1-P3 lesion level, sens / FPR: Elliott candidate 1.00/0.15 (v030 0.39, 0.2.0
  0.37); definite S>=0 0.57/0.04; mean rate >= 10 0.99/0.03. At 1 % scale 1.00/0.14,
  0.54/0.05, 0.99/0.03. Thick3 0.99/0.23, 0.58/0.07, 0.96/0.05. Lesjak x3 0.97/0.08,
  0.61/0.03, 0.96/0.03.
* 7 patients: candidates median 37 (154 in 0.2.0); candidate FPR 0.27 (0.04-0.59); definite
  0.71/0.10; noise floor 8.9 (5.1-14.4); AUC mean 0.99, constancy residual 0.76, S 0.76,
  concentricity 0.62.
* Greedy now close to ANTs 0.2.0 on P1/P3; P2 candidate FPR 0.48 vs 0.24.
* Cohorts: MSLesSeg 232 candidates, median 6, P2 none, 68 definite; Spearman 0.40. Lesjak
  208 / 45, 85 on expert change; tracking identical.
