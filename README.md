# lesiontrack

Longitudinal multiple sclerosis lesion evolution from two or more MRI timepoints:
which lesions are new, enlarging, shrinking, stable or resolved, and which are
**slowly expanding lesions (SELs)** by the Jacobian definition of Elliott et al. 2019,
with a **synthetic backtest gate** that has to pass before a real number is reported.

> **Nothing in this repository is evidence about multiple sclerosis.** It is research
> software. It is not a medical device, it has not been validated against histology or
> clinical outcome, and it makes no diagnostic claim.

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
5. **Scores** (Elliott 2019 as described by Vanden Bulcke et al. 2025): *concentricity*
   is the robust slope of mean expansion per 1-voxel shell against depth from the
   candidate edge; *constancy* is the normalised residual of a robust linear fit of
   expansion over the intermediate timepoints (needs three or more scans). Both are
   z-scored across the cohort; `S = z(concentricity) - z(constancy)`; `S >= 0` is a
   definite SEL. Because the z-scores are cohort-relative, the definite/possible split
   is a ranking within your cohort, not an absolute threshold; raw values are kept.
6. **Lesion tracking** from the binary masks: connected components at baseline and
   follow-up are matched on the bipartite overlap graph, so merges and splits are
   explicit. Volume change is annualised and classified with the +/-9 %/year change and
   +/-4 %/year stable bands of Vanden Bulcke et al. 2025; a follow-up component that
   clears the 1-voxel-dilated baseline mask is *new*.

## The backtest gate

`lesiontrack backtest` takes a real baseline (T1w, FLAIR, lesion mask), picks lesions,
expands each by a known volume factor (1.15, 1.25, 1.40) with a smooth radial
deformation, perturbs the whole head rigidly, rescales intensities and adds noise, then
runs the full pipeline on the pair and reports:

* sensitivity: expanded lesions that receive a SEL candidate
* false positive rate: untouched lesions that receive one
* rate error: measured mean expansion minus the true rate, in %/year

It is the only place the pipeline is checked against a known answer. Run it on your own
data before trusting a count. The numbers for MSLesSeg are in the section below.

## Install

```bash
pip install lesiontrack            # once released; until then:
pip install git+https://github.com/CedricConday/lesiontrack
```

Python 3.10+. Wheels for greedy exist for Linux x86-64 and arm64 and macOS.

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

Per subject you get, on the baseline grid in halfway space:

| file | content |
|---|---|
| `expansion_<tp>_pct_per_year.nii.gz` | Jacobian-derived local volume change |
| `sel_candidates.nii.gz`, `sel_candidates.tsv` | candidate labels; size, peak and mean expansion, concentricity, constancy |
| `baseline_lesion_labels.nii.gz`, `lesion_labels_<tp>.nii.gz` | lesion components |
| `lesion_tracking.tsv` | one row per matched group with class, volumes, rate, merge/split flags |
| `new_or_enlarging_<tp>.nii.gz` | follow-up lesion voxels outside the dilated baseline |
| `overview_<tp>.png` | four axial slices: baseline, warped follow-up, expansion with SEL outlines |
| `summary.json` | counts, volumes, intervals, every parameter used |

`lesiontrack cohort` adds `sel_cohort.tsv` with `z_concentricity`, `z_constancy`,
`score_S` and `definite_sel`.

## Results on public data

### Backtest on MSLesSeg patient P1 (2026-09-26)

Three synthetic series (seeds 0, 1, 2), each with scans at 0, 0.5 and 1 year, eight of
eighteen lesions expanded by nominal factors 1.15 / 1.25 / 1.40 (analytic mean inside the
lesion 11 to 25 % because elongated lesions extend past the radial window), rigid
perturbation about 1 degree and 1 voxel per scan, intensity rescale 0.9 to 1.1, 2 % noise.
Default registration (greedy stationary-velocity, NCC 3x3x3, halfway space).

| quantity | value |
|---|---|
| recovery fraction (measured / true, robust slope) | **0.30** |
| recovery by lesion size: >= 500 / 100-500 / < 100 voxels | 0.55 / 0.27 / 0.33 |
| candidate sensitivity, Elliott thresholds | 0.88 |
| candidate false positive rate, Elliott thresholds | 0.47 (1.00 on untouched lesions >= 500 voxels) |
| definite SEL (cohort z-scored S >= 0): sensitivity / FPR | 0.29 / 0.27 |
| noise floor: 95th percentile of per-voxel expansion in untouched lesions | 5 to 10 %/yr |
| candidate mean expansion as a classifier, AUC | 0.85 to 0.96 |

Three conclusions the tool would not have earned without the gate:

1. **The Jacobian under-reports by about 3x** on sub-voxel boundary shifts at 1 mm, and
   more so for small lesions. SEL rates from any pipeline of this kind are lower bounds.
2. **Elliott's candidate rule is noise-limited on confluent lesions.** A one-voxel peak
   over 12.5 %/yr, grown by hysteresis to 4 %/yr, makes a 10-voxel candidate; large
   lesions contain enough voxels to produce one by chance. Concentricity and constancy,
   the published filters, were at chance here (AUC 0.54 and 0.40).
3. **The candidate's mean rate is the filter that works** (AUC 0.90). It is exposed as
   `min_mean_pct_per_year`, off by default. Trade-off on this registration:

| min mean (%/yr) | sensitivity | false positive rate |
|---|---|---|
| 0 (Elliott) | 0.88 | 0.47 |
| 8 | 0.67 | 0.23 |
| 10 | 0.50 | 0.00 |

Full record with the registration sweep: `docs/BACKTEST_NOTES.md`.

### MSLesSeg cohort

Pending: 25 patients running; filled in from `derivatives/` when done.

## Data used

* **MSLesSeg** (Guarnera et al. 2025, Sci Data 12:920, CC BY): 25 patients with 2 to 4
  timepoints, already in MNI152 1 mm space with a lesion mask per timepoint.
  `scripts/make_manifest_mslesseg.py` builds the manifest; intervals are age differences.
* **Lesjak et al. 2016** longitudinal database (CC BY): 20 patients, 2 timepoints, expert
  masks of lesion change. Not yet used; it lacks baseline lesion masks.

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
* Masks are taken as given. A segmenter that hallucinates at follow-up will produce
  "new" lesions; the tool does not second-guess it.

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
