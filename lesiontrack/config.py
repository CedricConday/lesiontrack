"""Tunable parameters, each with its provenance.

The SEL thresholds follow Elliott et al. 2019 (Mult Scler 25:1915), which every
later SEL study except Yokote 2024 reused. They were set by visual assessment,
not by ground truth; the synthetic backtest in :mod:`lesiontrack.synth` is the
first check of what they actually recover.
"""

from dataclasses import dataclass, field


@dataclass(frozen=True)
class SELParams:
    je1_pct_per_year: float = 12.5  # seed threshold, Elliott 2019
    je2_pct_per_year: float = 4.0  # hysteresis growth threshold, Elliott 2019
    min_voxels: int = 10  # reliability criterion, Elliott 2019
    connectivity: int = 2  # scipy structure rank: 2 = 18-connectivity in 3D
    # Not in Elliott 2019. A candidate's mean expansion separates injected expansion from
    # registration noise where concentricity and constancy do not (candidate AUC 0.99 against
    # 0.62 and 0.76 over seven MSLesSeg backtests, docs/VALIDATION.md). The floor is on the
    # total change over the interval, not on the yearly rate: registration noise is a
    # property of one registration, so in %/yr it shrinks with the interval (mean rate of
    # untouched lesions, SD 3.7 %/yr at one year, 1.7 at two). 10 % total is 10 %/yr at one
    # year, where it keeps 0.99 of lesions injected at 15 %/yr or more and cuts the false
    # positive rate from 0.15 to 0.05, and 5 %/yr at two. 0.3.0.dev0 used a fixed 10 %/yr,
    # which at two years removed most lesions expanding at 8 %/yr and no false positive
    # (2026-11-23, "Slow expansion" in docs/VALIDATION.md). 0 turns it off.
    min_mean_total_pct: float = 10.0
    # An absolute floor in %/yr, applied as well when above 0. Off by default.
    min_mean_pct_per_year: float = 0.0

    def min_mean_at(self, dt_years: float | None) -> float:
        """The mean-rate floor in %/yr for a pair ``dt_years`` apart (the larger of the two)."""
        total = self.min_mean_total_pct / dt_years if dt_years and self.min_mean_total_pct > 0 else 0.0
        return max(total, self.min_mean_pct_per_year)

    def __post_init__(self) -> None:
        if self.je1_pct_per_year < self.je2_pct_per_year:
            raise ValueError("je1 (seed) must be >= je2 (growth); a swapped pair yields no candidates")


@dataclass(frozen=True)
class TrackParams:
    # Vanden Bulcke 2025 classified lesion volume slopes beyond +/-9 %/year as
    # expanding/shrinking and within +/-4 %/year as stable.
    change_pct_per_year: float = 9.0
    stable_pct_per_year: float = 4.0
    min_change_voxels: int = 3  # below this a percentage is noise on a tiny lesion
    new_lesion_margin_voxels: int = 1  # follow-up component must clear the dilated baseline
    # Not from the literature. A "new" or "resolved" call rests on two independent
    # segmentations; it is marked supported when FLAIR itself moved the same way: the median
    # change in robust z (units: MAD of non-lesion brain) over the group is at least this.
    # On Lesjak 2016 with LST-AI masks, against expert change, 1.0 keeps 0.84 of the
    # confirmed new lesions and drops 0.65 of the unconfirmed (AUC 0.83); for resolved,
    # 0.60 and 0.75 (AUC 0.70). Chosen and measured on the same 20 patients (2026-11-23).
    flair_support_z: float = 1.0
    connectivity: int = 2


@dataclass(frozen=True)
class RegParams:
    # Elliott 2019: "step size = 0.7; Gaussian sigma = 2", T1 and T2/FLAIR jointly.
    # 7 = rigid plus one isotropic scale (since 0.3.0; 0.2.0 used 6, rigid). A scanner or
    # gradient calibration change between timepoints scales the whole image; a rigid step
    # leaves that to the deformable registration, which puts it unevenly into textured
    # tissue such as lesion edges. On the P1 backtest a simulated 1 % scale change raised
    # the lesion-level false positive rate of the mean rate from 0.03 to 0.23 with dof 6;
    # with dof 7 it was 0.00 (greedy measured +0.95 %; -0.03 % with no change). Side effect
    # on real data: whole-brain atrophy is absorbed too, so lesion change is measured
    # relative to the brain. 6 restores the rigid step.
    affine_dof: int = 7
    affine_metric: str = "NMI"
    affine_iterations: str = "100x50x20"
    # Backtest sweep 2026-09-26 (docs/BACKTEST_NOTES.md): greedy's stationary-velocity mode
    # recovers as much injected expansion as ten times the iterations of the plain greedy
    # update at a sixth of the time; a 3-voxel NCC window edged out 2 and 4.
    deform_metric: str = "NCC 3x3x3"
    deform_step: float = 0.7
    deform_sigma_update: str = "2vox"
    deform_sigma_total: str = "0.5vox"
    deform_iterations: str = "100x60x30"
    stationary_velocity: bool = True
    # greedy rounds the saved warp to 0.1 voxel unless told otherwise (-wp, default 0.1).
    # The Jacobian is a finite difference of that field, so the rounding put every rate on
    # steps of about 5 %/yr at a one-year interval: Elliott's 12.5 %/yr seed is two to
    # three steps and the 4 %/yr growth threshold less than one. 0 saves the warp unrounded
    # (found 2026-11-21; 0.2.0 and earlier saved it rounded).
    warp_precision_vox: float = 0.0
    t1_weight: float = 1.0
    flair_weight: float = 1.0
    threads: int = 4
    # Second engine (2026-09-27): ANTs SyN through antspyx, on the same halfway-space images.
    # The rigid step and the reslicing stay greedy; only the deformable step and the
    # Jacobian change, so the two engines are compared on identical inputs.
    engine: str = "greedy"  # "greedy" or "ants"
    ants_transform: str = "SyNOnly"
    ants_metric: str = "CC"  # cross-correlation, the SyN default in the ANTs papers
    ants_radius: int = 4  # CC neighbourhood radius (for a histogram metric: bins), both terms
    ants_iterations: tuple = (100, 60, 30)
    ants_grad_step: float = 0.2
    ants_flow_sigma: float = 3.0
    ants_total_sigma: float = 0.0


@dataclass(frozen=True)
class Params:
    sel: SELParams = field(default_factory=SELParams)
    track: TrackParams = field(default_factory=TrackParams)
    reg: RegParams = field(default_factory=RegParams)


DEFAULT = Params()
