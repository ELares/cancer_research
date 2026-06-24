"""Map the Talkington & Kearsley (2025) multiscale ICB model parameters onto the
repo's spatial-immune model (#472).

Source (verified): A. M. Talkington & A. J. Kearsley, "Optimization of Immune
Checkpoint Blockade via a Multiscale Model System", Computational and Systems
Oncology 5(2):e70007, 2025. DOI 10.1002/cso2.70007, PMID 41322398, PMC12663532
(open access, CC BY 4.0; parameters in Table 1; no public code/data, values
read from the Europe PMC full-text XML).

This module is DOCUMENTATION + a reproducible encoding of the mapping. The
borrowed values are OPT-IN sanity ranges, NOT swapped into the production
defaults: our spatial agent-based model is uncalibrated, and a well-mixed ODE
rate cannot be transplanted directly (the spatial geometry alone shrinks the
SDT:RSL3 immune ratio ~104:1 -> ~4:1 in 3D, so a mass-action coefficient that
assumes every effector reaches every target would double-count the spatial
limitation). Only the DIMENSIONLESS quantities (the ICB efficiency `I` and the
exhaustion threshold `n`) port cleanly; the DIMENSIONAL rates do not without a
per-step time and a local cell density.

Run `python scripts/icb_param_map.py` to print the mapping table.
"""

from __future__ import annotations

import math
from dataclasses import dataclass


@dataclass(frozen=True)
class TalkingtonParam:
    """One published parameter with its meaning, value, units, and provenance."""

    symbol: str
    value: float
    units: str
    meaning: str
    provenance: str  # "cited" (mouse literature), "assumed", or "varied" (swept)


# Verified Table 1 of Talkington & Kearsley 2025 (the whole-tumor ODE system
# T = effector T cells, C = cancer cells, A = non-cancer APCs; exhaustion
# coefficient a = (1/n)*I*R1*R2). Values quoted exactly as the paper reports them.
TALKINGTON_PARAMS: dict[str, TalkingtonParam] = {
    "kill_rate": TalkingtonParam(
        symbol="r2",
        value=1.101e-7,
        units="day^-1 cells^-1",
        meaning="rate of cancer-cell killing by effector T cells (the r2*C*T term)",
        provenance="cited",  # Kuznetsov 1994 / Talkington 2018, mouse; MC-sampled
    ),
    "exhaustion_threshold": TalkingtonParam(
        symbol="n",
        value=1.0e4,
        units="interactions (cells)",
        meaning="min. T-cell/(A+C) interactions needed to exhaust one T cell",
        provenance="assumed",
    ),
    "icb_efficiency": TalkingtonParam(
        symbol="I",
        value=0.15,  # representative: the 80-90% blockade optimum is I ~ 0.1-0.2
        units="dimensionless",
        meaning=(
            "checkpoint-blockade control: I=0 perfect blockade (no exhausting "
            "interactions), I=1 no blockade (all interactions exhaust)"
        ),
        provenance="varied",  # swept over U(0,1); 80-90% blockade (I~0.1-0.2) optimal
    ),
    "t_stimulation": TalkingtonParam(
        symbol="r1",
        value=0.1245,
        units="day^-1",
        meaning="T-cell stimulation rate (r1*T*C/(k1+C))",
        provenance="cited",
    ),
    "t_death": TalkingtonParam(
        symbol="d_T",
        value=0.0412,
        units="day^-1",
        meaning="effector T-cell death rate",
        provenance="cited",
    ),
}


def kill_rate_to_per_step_probability(
    r2: float, effector_density: float, dt_days: float
) -> float:
    """Convert the bimolecular ODE kill rate `r2` (day^-1 cells^-1) to our model's
    per-cell, per-step kill PROBABILITY.

    The paper's `r2*C*T` is a mass-action flux, not a probability. The agent-based
    analog (which the paper's own ABM uses) is a per-contact / per-step hazard:

        p = 1 - exp(-r2 * effector_density * dt_days)

    so the raw `r2` alone is NOT transplantable: it needs a LOCAL effector-cell
    count and a per-step time `dt_days`. `effector_density` here is the number of
    effector T cells in a cell's local interaction neighborhood (NOT the global
    count `T`, which would re-introduce the well-mixed assumption our spatial
    model deliberately drops). Returns a probability in [0, 1).
    """
    if effector_density < 0 or dt_days < 0:
        raise ValueError("effector_density and dt_days must be non-negative")
    return 1.0 - math.exp(-r2 * effector_density * dt_days)


def icb_efficiency_to_checkpoint_residual(icb_efficiency: float) -> float:
    """Map the paper's single ICB efficiency `I` to our per-checkpoint
    residual-after-drug.

    Their `I` is the residual exhausting fraction after blockade: I=0 is perfect
    blockade (our residual 0), I=1 is no drug (our residual 1). Our multi-checkpoint
    brake `1 - prod(1 - residual_i*(1-drug_eff_i))` collapses to their single `I`
    when one checkpoint carries it and the rest are zeroed. The 80-90%-blockade
    optimum gives a plausible AGGREGATE residual ~0.1-0.2 (do NOT assign it
    per-checkpoint). Clamped to [0, 1].
    """
    return min(1.0, max(0.0, icb_efficiency))


def exhaustion_threshold_to_rate(n_interactions: float) -> float:
    """Map the paper's exhaustion threshold `n` (interactions to exhaust a T cell)
    to a scale for our `exhaustion_rate` in `1/(1 + exhaustion_rate*cumulative)`.

    Their exhaustion is mass-action with threshold `n`; our suppression is a
    saturating function of `cumulative_kills`. The natural correspondence is
    `exhaustion_rate ~ 1/n` (the cumulative-interaction count at which suppression
    becomes appreciable is `~n`). With n=1e4 that is 1e-4 per interaction, an
    ORDER-OF-MAGNITUDE anchor, not a calibration (the two functional forms differ).
    """
    if n_interactions <= 0:
        raise ValueError("n_interactions must be positive")
    return 1.0 / n_interactions


# Which of OUR parameters each published value informs, and whether it ports.
MAPPING: list[tuple[str, str, str]] = [
    (
        "per-cell immune kill probability",
        "r2 = 1.101e-7 day^-1 cells^-1 (cited, mouse)",
        "DIMENSIONAL: needs p=1-exp(-r2*local_density*dt); raw rate not transplantable "
        "(would double-count spatial locality). Use the paper's ABM per-contact form, "
        "calibrate magnitude to our grid.",
    ),
    (
        "exhaustion_rate (1/(1+rate*cumulative_kills))",
        "n = 1e4 interactions (assumed)",
        "DIMENSIONLESS scale: exhaustion_rate ~ 1/n ~ 1e-4. Order-of-magnitude anchor; "
        "functional forms differ (mass-action vs saturating).",
    ),
    (
        "multi-checkpoint residual-after-drug (PD-1/CTLA-4/LAG-3/TIM-3)",
        "I in [0,1], 80-90% blockade optimum (I~0.1-0.2) (varied)",
        "DIMENSIONLESS, ports cleanly to ONE AGGREGATE residual (not per-checkpoint). "
        "Their single lumped axis is coarser than our 4-checkpoint panel.",
    ),
    (
        "DAMP diffusion / suppressor field / IFN-gamma loop / 3D geometry",
        "(no counterpart in the paper)",
        "REPO-SPECIFIC, stays uncalibrated. The paper is well-mixed and has no spatial "
        "coupling, so these gain no external anchor here.",
    ),
]


def _print_table() -> None:
    print("Talkington & Kearsley 2025 ICB parameters (verified Table 1)")
    print("DOI 10.1002/cso2.70007 | PMID 41322398 | PMC12663532 (CC BY 4.0)\n")
    for key, p in TALKINGTON_PARAMS.items():
        print(f"  {key:22s} {p.symbol:5s} = {p.value:<12g} {p.units:18s} [{p.provenance}]")
        print(f"  {'':22s} {p.meaning}")
    print("\nMapping onto the repo's spatial-immune model:")
    for ours, theirs, note in MAPPING:
        print(f"  OURS:   {ours}")
        print(f"  THEIRS: {theirs}")
        print(f"  NOTE:   {note}\n")
    # A worked conversion example (illustrative numbers only).
    p = kill_rate_to_per_step_probability(
        TALKINGTON_PARAMS["kill_rate"].value, effector_density=1.0e6, dt_days=0.25
    )
    print(
        "Worked example: r2=1.101e-7, local effector density 1e6 cells, dt=0.25 day"
        f" -> per-step kill probability {p:.4f}"
    )
    print(
        "  (illustrative; the local density and dt are model-specific, so this is a "
        "form-of-the-conversion demonstration, not a calibrated value.)"
    )


if __name__ == "__main__":
    _print_table()
