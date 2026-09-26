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


@dataclass(frozen=True)
class TrackParams:
    # Vanden Bulcke 2025 classified lesion volume slopes beyond +/-9 %/year as
    # expanding/shrinking and within +/-4 %/year as stable.
    change_pct_per_year: float = 9.0
    stable_pct_per_year: float = 4.0
    min_change_voxels: int = 3  # below this a percentage is noise on a tiny lesion
    new_lesion_margin_voxels: int = 1  # follow-up component must clear the dilated baseline
    connectivity: int = 2


@dataclass(frozen=True)
class RegParams:
    # Elliott 2019: "step size = 0.7; Gaussian sigma = 2", T1 and T2/FLAIR jointly.
    affine_dof: int = 6
    affine_metric: str = "NMI"
    affine_iterations: str = "100x50x20"
    deform_metric: str = "NCC 2x2x2"
    deform_step: float = 0.7
    deform_sigma_update: str = "2vox"
    deform_sigma_total: str = "0.5vox"
    deform_iterations: str = "100x60x30"
    t1_weight: float = 1.0
    flair_weight: float = 1.0
    threads: int = 4


@dataclass(frozen=True)
class Params:
    sel: SELParams = field(default_factory=SELParams)
    track: TrackParams = field(default_factory=TrackParams)
    reg: RegParams = field(default_factory=RegParams)


DEFAULT = Params()
