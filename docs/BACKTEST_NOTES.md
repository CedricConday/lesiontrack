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

Stationary-velocity mode matches ten times the iterations at a sixth of the time.

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
