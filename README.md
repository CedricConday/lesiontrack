# lesiontrack

[![tests](https://github.com/CedricConday/lesiontrack/actions/workflows/tests.yml/badge.svg)](https://github.com/CedricConday/lesiontrack/actions/workflows/tests.yml)

Longitudinal multiple sclerosis lesion evolution from two or more MRI timepoints:
which lesions are new, enlarging, shrinking, stable or resolved, and which are
**slowly expanding lesions (SELs)** by the Jacobian definition of Elliott et al. 2019,
with a **synthetic backtest gate** that has to pass before a real number is reported.

> **Nothing in this repository is evidence about multiple sclerosis.** It is research
> software. It is not a medical device, it has not been validated against histology or
> clinical outcome, and it makes no diagnostic claim.

![MSLesSeg patient P20, baseline to last follow-up: baseline lesions (cyan), new voxels on the warped follow-up (yellow), Jacobian expansion inside lesions with SEL candidate outlines (lime)](https://raw.githubusercontent.com/CedricConday/lesiontrack/master/results/overview_P20.png)

## Why this exists

Every MS clinic that pays for follow-up MRI analysis gets lesion-level change and, in
trials, SEL counts. In the open there was one SEL implementation, bound to a GUI, needing
FreeSurfer and ANTs, taking hours per subject and carrying its authors' own warning not
to use it. There is no public SEL ground truth, and the field's thresholds were set by
visual assessment. This package is the missing open engine, and the backtest is the
missing check.

* one dependency for registration: [greedy](https://github.com/pyushkevich/greedy),
  installed as a wheel, x86-64 and arm64, CPU only
* minutes per subject pair, not hours
* every threshold in one file with its citation (`lesiontrack/config.py`)
* every output a plain NIfTI or TSV, on the baseline grid, ready for any viewer

## Method

1. **Halfway space.** Rigid registration of the follow-up T1w to the baseline T1w, then
   the square root of that transform moves both scans into an unbiased midpoint space
   (Nakamura et al.; Elliott et al. 2019).
2. **Deformable registration** driven by T1w and FLAIR jointly (normalised cross
   correlation, step 0.7, update smoothing 2 voxels), baseline fixed, follow-up moving.
3. **Jacobian determinant** of the deformation field, converted to local percent volume
   change per year: `(det J - 1) * 100 / years`.
4. **SEL candidates** (Elliott 2019): voxels inside baseline lesions expanding at
   >= 12.5 %/year seed 18-connected components; each seed grows by hysteresis into
   neighbours expanding at >= 4 %/year, competing so distinct expansions stay distinct
   (marker watershed on the expansion map); candidates under 10 voxels are dropped.
   Not Elliott's: candidates whose mean expansion over the whole interval is under 10 %
   are dropped too (`--min-mean-total`; `--min-mean 0` gives the published rule). The
   backtest is why, see Results.
5. **Scores** (Elliott 2019 as described by Vanden Bulcke et al. 2025): *concentricity*
   is the robust slope of mean expansion per 1-voxel shell against depth from the
   candidate edge; *constancy* is the normalised residual of a robust linear fit of
   expansion over the intermediate timepoints (needs three or more scans). Both are
   z-scored across the cohort; `S = z(concentricity) - z(constancy)`; `S >= 0` is a
   definite SEL. Because the z-scores are cohort-relative, the definite/possible split
   is a ranking within your cohort, not an absolute threshold; raw values are kept.
   Candidates with no interior voxel have no concentricity and are marked unscorable.
6. **Lesion tracking** from the binary masks: connected components at baseline and
   follow-up are matched on the bipartite overlap graph, so merges and splits are
   explicit. Volume change is annualised and classified as enlarging or shrinking beyond
   +/-9 %/year, stable within +/-4 %/year, and trend_up / trend_down between (bands of
   Vanden Bulcke et al. 2025). A follow-up component with no baseline overlap is *new* if
   it clears the 1-voxel-dilated baseline mask and *adjacent_fragment* if it touches it;
   a baseline component with no follow-up overlap is *resolved*. New and resolved calls
   rest on two segmentations, so each is checked against the images: `flair_supported`
   says whether FLAIR intensity there rose (new) or fell (resolved) by at least one
   robust z.

## The backtest gate

`lesiontrack backtest` takes a real baseline (T1w, FLAIR, lesion mask), picks lesions,
expands each by a known volume factor (1.15, 1.25, 1.40, or your own with `--factors`
and `--dt-years`) with a smooth radial
deformation, perturbs the whole head rigidly, rescales intensities and adds noise, then
runs the full pipeline on the pair and reports:

* sensitivity: expanded lesions that receive a SEL candidate
* false positive rate: untouched lesions that receive one
* recovery fraction: robust slope of measured against true expansion (1.0 = exact)
* the noise floor inside untouched lesions and, per candidate, how well its mean rate
  separates injected from spurious candidates

It is the only place the pipeline is checked against a known answer. Run it on your own
data, at the rates and interval you care about, before trusting a count.

## Install

```bash
pip install lesiontrack
```

Python 3.10+. Wheels for greedy exist for Linux x86-64 and arm64 and macOS. This README
describes 0.3.x (`lesiontrack --version`). Use 0.3.2: in 0.3.0 and 0.3.1 the lesion masks
were resliced with the similarity's scale, which added or dropped a voxel layer on some
lesions and called a few unchanged ones enlarging (see docs/VALIDATION.md, "0.3.2"). 0.3.1
is 0.3.0 with the complete validation text. 0.2.0 saved rounded warps and had no mean-rate
floor, so upgrade before comparing numbers.

## Use

Write a tab-separated manifest, one row per scan, paths relative to the manifest:

```
subject	session	time_years	t1	flair	mask	brainmask
P1	T1	28.09	P1/T1/T1.nii.gz	P1/T1/FLAIR.nii.gz	P1/T1/MASK.nii.gz	P1/T1/brain.nii.gz
P1	T2	28.51	P1/T2/T1.nii.gz	P1/T2/FLAIR.nii.gz	P1/T2/MASK.nii.gz
P1	T3	29.49	P1/T3/T1.nii.gz	P1/T3/FLAIR.nii.gz	P1/T3/MASK.nii.gz
```

`time_years` is any consistent clock (age, or years since baseline). Masks come from
whatever segmenter you trust (LST-AI, SAMSEG, nnU-Net, manual); the brain mask is
optional and only restricts the registration metric.

```bash
lesiontrack run manifest.tsv --out derivatives          # every subject
lesiontrack run manifest.tsv --out derivatives --subject P1
lesiontrack cohort --out derivatives                    # z-scores and S across subjects
lesiontrack backtest --t1 T1.nii.gz --flair FLAIR.nii.gz --mask MASK.nii.gz --out backtest
```

### Two deformable engines

The linear step (similarity since 0.3.0: rigid plus one global scale), the halfway frame
and the reslicing are greedy in both cases; `--engine`
picks what runs on the halfway-space images. `greedy` (default) needs nothing extra.
`ants` uses SyN through [antspyx](https://github.com/ANTsX/ANTsPy), cross-correlation on
T1 and FLAIR jointly, and is how the original SEL studies registered:

```bash
pip install 'lesiontrack[ants]'
lesiontrack run manifest.tsv --out derivatives_ants --engine ants
lesiontrack backtest --t1 T1.nii.gz --flair FLAIR.nii.gz --mask MASK.nii.gz --out backtest_ants --engine ants
```

Both engines write the same files, so the backtest compares them on identical inputs. The
Jacobian determinant is computed by lesiontrack from the warp, in the image's voxel frame,
for either engine; see "The Jacobian is ours" below for why.

Per subject you get, on the baseline grid in halfway space:

| file | content |
|---|---|
| `expansion_<tp>_pct_per_year.nii.gz` | Jacobian-derived local volume change |
| `sel_candidates.nii.gz`, `sel_candidates.tsv` | candidate labels; size, peak and mean expansion, concentricity, constancy |
| `baseline_lesion_labels.nii.gz`, `lesion_labels_<tp>.nii.gz` | lesion components |
| `lesion_tracking.tsv` | one row per matched group with class, volumes, rate, merge/split flags, `flair_dz`, `flair_supported` |
| `new_or_enlarging_<tp>.nii.gz` | follow-up lesion voxels outside the dilated baseline |
| `overview_<tp>.png` | four axial slices: baseline, warped follow-up, expansion with SEL outlines |
| `summary.json` | counts, volumes, intervals, every parameter used |

`lesiontrack cohort` adds `sel_cohort.tsv` with `z_concentricity`, `z_constancy`,
`score_S` and `definite_sel`.

## Results on public data

Current numbers, lesiontrack 0.3.0, greedy engine. Every earlier version's numbers, the
ANTs comparison and the per-patient tables are in [docs/VALIDATION.md](https://github.com/CedricConday/lesiontrack/blob/master/docs/VALIDATION.md).

### Backtest: injected expansion of 15 to 40 % over one year

MSLesSeg patients P1 to P3, three seeds each, scans at 0, 0.5 and 1 year; 72 expanded and
150 untouched lesions per condition. Lesion-level sensitivity / false positive rate:

| condition | Elliott candidate | definite S >= 0 | default (mean rate floor) |
|---|---|---|---|
| no scanner change | 1.00 / 0.15 | 0.57 / 0.04 | 0.99 / 0.03 |
| 0.5 % scale change | 1.00 / 0.12 | 0.54 / 0.03 | 0.99 / 0.03 |
| 1 % scale change | 1.00 / 0.14 | 0.54 / 0.05 | 0.99 / 0.03 |
| FLAIR at 3 mm slices | 0.99 / 0.23 | 0.58 / 0.07 | 0.96 / 0.05 |
| Lesjak 03, 06, 07 (0.72 x 0.72 x 3 mm) | 0.97 / 0.08 | 0.61 / 0.03 | 0.96 / 0.03 |

Over seven MSLesSeg patients (P1 P2 P3 P12 P13 P19 P22, picked before any result was
seen; median and range):

| quantity | value |
|---|---|
| recovery fraction (measured / injected, robust slope) | 1.17 |
| candidate false positive rate, Elliott thresholds | 0.27 (0.04 to 0.59) |
| noise floor: median peak per-voxel rate in untouched lesions | 8.9 %/yr (5.1 to 14.4) |
| candidate AUC, injected against spurious: mean rate | 0.99 (0.87 to 1.00) |
| candidate AUC: S / constancy residual / concentricity | 0.76 / 0.76 / 0.62 |

What the gate showed:

1. **The Jacobian recovers the injected expansion**, slightly over on small lesions.
2. **Elliott's candidate rule is noise-limited on large lesions.** One voxel over
   12.5 %/yr, grown by hysteresis to 4 %/yr, makes a candidate; S >= 0 removes most of
   those and about 40 % of the true ones with them.
3. **The candidate's mean rate separates best**, which is why `run` applies a floor on it.

### Backtest: slow expansion

The design above injects 11 to 25 %/yr inside the lesion, and the floor was first chosen
on it. SELs are defined by slow expansion, so the backtest was repeated at nominal 3, 5
and 8 %/yr, three scans, over two years and over one (`scripts/run_slow.sh`; P1 to P3,
three seeds each, 72 expanded and 150 untouched lesions per interval; 2026-11-23). Share
of lesions that keep a SEL candidate:

| interval | injected rate | Elliott's rule | default floor (10 % total) | fixed 10 %/yr (0.3.0.dev0) |
|---|---|---|---|---|
| 2 years | 3 %/yr (n = 27) | 0.15 | 0.15 | 0.04 |
| 2 years | 5 %/yr (n = 27) | 0.26 | 0.26 | 0.07 |
| 2 years | 8 %/yr (n = 18) | 0.89 | 0.89 | 0.50 |
| 2 years | untouched (false positive rate) | 0.01 | 0.01 | 0.00 |
| 1 year | 3 %/yr (n = 27) | 0.59 | 0.26 | 0.26 |
| 1 year | 5 %/yr (n = 27) | 0.70 | 0.37 | 0.37 |
| 1 year | 8 %/yr (n = 18) | 0.94 | 0.89 | 0.89 |
| 1 year | untouched (false positive rate) | 0.12 | 0.05 | 0.05 |

What this says:

1. **A fixed floor in %/yr is wrong for longer intervals.** Noise in %/yr halves over two
   years (mean rate of untouched lesions, SD 1.8 against 3.7 %/yr). At two years a
   10 %/yr floor removed nearly half the lesions expanding at 8 %/yr to remove one false
   positive in 150. The default is therefore a floor on total change, 10 % over the
   interval, which at two years costs nothing.
2. **At one year the floor is a real trade.** It takes the false positive rate from 0.12
   to 0.05 and detection of lesions expanding at 3 to 5 %/yr from about two thirds to
   about one third. `--min-mean 0` gives the published rule and the other side of that
   trade; `--min-mean 8` sits between (0.48 / 0.59 / 0.94, false positives 0.07).
3. **Uniform expansion at 3 to 5 %/yr is mostly missed over two years by Elliott's rule
   itself**: few such lesions hold a voxel at the 12.5 %/yr seed once the noise is lower.
   The one-year detections of the same lesions lean on that noise: the measured rate
   over-reads small expansions (median 6.4 %/yr measured for 3.0 injected at one year,
   4.6 at two).
4. The lesion's mean Jacobian rate separates expanded from untouched lesions at AUC 0.96
   (two years) and 0.89 (one year), so the signal is in the map more than in the
   candidates.

### Cohorts with the defaults

MSLesSeg, 24 patients, baseline against last scan, median interval 1.5 years (0.16 to
6.0): 374 SEL candidates, median 12.5 per patient, 92 definite. The fixed 10 %/yr floor
gave 232 and 6, Elliott's rule without a floor 484 and 16 (0.2.0); the literature reports
about 2 to 4 SELs per patient over two years, so these counts are still noise-heavy.
Mask-based against Jacobian volume change per subject: Spearman 0.40. Of 216 new and 250
resolved lesions by the masks, FLAIR supports 55 and 40. One patient (P49) is refused by
the input gate: its follow-up T1 is 0.1 % nonzero.

Lesjak et al. 2016, 20 patients, median interval 0.87 years, LST-AI masks, scored against
expert masks of lesion change (there is no SEL ground truth, so this tests tracking):
162 SEL candidates, median 7 per patient, 82 on expert change.

| tracked groups called changed | sensitivity | specificity | PPV |
|---|---|---|---|
| by the masks alone | 0.95 | 0.08 | 0.22 |
| new and resolved counted only when `flair_supported` | 0.81 | 0.59 | 0.35 |

The masks call 1,426 of 2,020 groups new or resolved, against 282 expert-confirmed; the
churn is the segmenter's. FLAIR supports 346 of the 710 new lesions (168 of the 200
expert-confirmed ones) and 207 of the 716 resolved (49 of 82). The support threshold was
chosen on these same 20 patients, so the second row is an in-sample number.

## Data used

* **MSLesSeg** (Guarnera et al. 2025, Sci Data 12:920, CC BY): 25 patients with 2 to 4
  timepoints, already in MNI152 1 mm space with a lesion mask per timepoint.
  `scripts/make_manifest_mslesseg.py` builds the manifest; intervals are age differences.
* **Lesjak et al. 2016** longitudinal database (CC BY): 20 patients, 2 timepoints, expert
  masks of lesion change. It has no lesion masks, so they come from LST-AI; see
  "Validation beyond P1".

## The Jacobian is ours

Until 0.2.0 the Jacobian came from greedy's `-rj` in the reslice command. That option
differentiates the composed warp while it is still in physical (LPS, millimetre) units,
along voxel index axes, without the image direction matrix, so it is only right for
images whose voxel axes are LPS-aligned at 1 mm. On a phantom with an identity (RAS)
affine, a sphere growing from radius 8 to 10 voxels (volume ratio 1.95) came back as det
0.68; with a LAS affine, as used by the MSLesSeg scans, as 1.18. For a radial stretch
1+k the true determinant is (1+k)^3; each flipped axis turns one factor into (1-k), so
RAS (two flips) reports (1-k)^2 (1+k), which is below 1 for any expansion, and LAS (one
flip) reports (1-k)(1+k)^2, roughly a third of the true volume change for small k. greedy's standalone `-jac` converts the warp first and is
correct; the reslice path does not (found 2026-09-27, draft report in
`docs/upstream/greedy_rj_direction.md`).

A second loss of precision sat in the warp file itself: greedy writes it rounded to 0.1
voxel by default (`-wp 0.1`), which put every rate on steps of about 5 %/yr. Since 0.3.0
the warp is written unrounded (`-wp 0`, found 2026-11-21).

lesiontrack now computes det(I + grad u) itself from the warp file, after taking the
vectors LPS -> RAS -> voxel units through the inverse affine, and `tests/test_jacobian.py`
checks it under four orientations and anisotropic spacing. Registrations cached under an
older record are not redone; their Jacobian is recomputed from the stored warp.

## Limitations, stated plainly

* The Jacobian sees any local volume change, including registration error at lesion
  edges, atrophy, and interpolation. The backtest measures the pipeline's own noise floor
  on untouched lesions; report it with your SEL counts.
* Two-timepoint studies cannot score constancy; S then rests on concentricity alone and
  is labelled so.
* Thresholds are Elliott's, set by eye in 2019 on 3D scans resampled to 1x1x3 mm. Vanden
  Bulcke et al. 2025 found only weak lesion-level agreement (Dice about 0.4) between SEL
  masks from different segmenters or resolutions. Subject-level counts were stable;
  lesion-level identity was not. Expect the same here.
* Slow uniform expansion is missed. Elliott's seed needs a voxel at 12.5 %/yr; lesions
  injected at 3 to 5 %/yr over two years got a candidate in about a fifth of cases. The injected expansion is a
  uniform scaling, not the centre-weighted pattern SELs are described by, so real SELs
  may fare better or worse; nothing public can tell.
* A re-run is not bit-identical. greedy with several threads does not reproduce a
  registration exactly: the same patient, settings and thread count gave candidates of
  129 and 177 voxels where the first run gave 131 and 184, and a 10-voxel candidate at the
  size limit appeared. Counts near a threshold can move by one between runs.
* The mean-rate floor is this package's, not Elliott's, and it was chosen on synthetic
  expansion. `--min-mean 0` reports the published rule.
* Masks are taken as given. A segmenter that hallucinates at follow-up will produce
  "new" lesions; the class is kept and `flair_supported` says whether the images agree.
  On Lesjak with LST-AI masks the unfiltered call has specificity 0.08 for "changed".
* Since 0.3.0 the linear step absorbs a global scale, including real whole-brain atrophy:
  lesion expansion is relative to the brain. The simulated scanner change is a pure linear
  scale; a real scanner change also distorts non-linearly, and that is not tested.
* Every SEL number here comes from synthetic expansion injected into real scans. No public
  dataset has SEL ground truth; the backtest measures the pipeline, not the biology.

## References

* Elliott C, et al. Slowly expanding/evolving lesions as a magnetic resonance imaging
  marker of chronic active multiple sclerosis lesions. Mult Scler 2019;25:1915-1925.
* Vanden Bulcke C, et al. Do slowly expanding lesions correspond to chronic active
  multiple sclerosis lesions? An integrated imaging analysis study. medRxiv 2025,
  doi:10.64898/2025.12.03.25341192.
* Guarnera F, et al. MSLesSeg: baseline and benchmarking of a new Multiple Sclerosis
  Lesion Segmentation dataset. Sci Data 2025;12:920.
* Yushkevich PA, et al. Fast automatic segmentation of hippocampal subfields and medial
  temporal lobe subregions in 3 Tesla and 7 Tesla T2-weighted MRI (greedy). 2016.

## License

MIT. Written by Cedric Conday with Claude (Anthropic) as coding partner.
