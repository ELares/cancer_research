# PDT/SDT exo-ROS kill-threshold validation (Zhu 2015, #464)

This validates a load-bearing design choice in the simulation against measured
in-vivo photophysics. The model feeds SDT and PDT exogenous ROS into the **same**
bistable lipid-peroxidation death switch, with **one** `death_threshold` and
`sdt_ros == pdt_ros` (both 5.0). That encodes the assumption that the dose-to-kill
is a property of the cell/tissue, not of the particular ROS source.

- **Validation script:** `scripts/validate_pdt_threshold.py` (pure stdlib, runs in
  CI; drift-guards `sdt_ros`/`pdt_ros`/`death_threshold` against `params.rs`).
- **Committed result:** `analysis/calibration/pdt-threshold-validation.json`.
- **Anchor:** Zhu TC, Kim MM, Liang X, Finlay JC, Busch TM, "In-vivo singlet
  oxygen threshold doses for PDT," Photonics Lasers Med 2015 (PMID 25927018,
  PMC4410434).

## What the data says

Zhu 2015 measured the reacted-singlet-oxygen necrosis threshold [1O2]rx in vivo for
three clinically distinct photosensitizers:

| Photosensitizer | [1O2]rx threshold |
|-----------------|-------------------|
| Photofrin | ~0.56 mM |
| BPD (verteporfin) | ~0.72 mM |
| mTHPC | ~0.40 mM |

The threshold is **approximately photosensitizer-independent**: ~0.5 mM, within a
factor of **1.8** across three chemically very different sensitizers. The reacted-
singlet-oxygen dose that kills tissue is therefore a property of the target, not of
the photosensitizer that produced the ROS.

## What this validates in the model

The model uses **one** `death_threshold` (10.0 LP units) and sets `sdt_ros == pdt_ros`
(5.0), so SDT and PDT ROS feed the same bistable switch at the same kill threshold.
Zhu 2015's photosensitizer-independence is **real-world support for exactly that
design choice**: if the reacted-ROS kill dose were strongly source-dependent, the
model would need per-source thresholds; it is not, so a single source-independent
`death_threshold` is the right structure. This is the same logic that lets the model
compare SDT, PDT, and RSL3 on one death axis.

## Scope and honesty

- **Form + order of magnitude, not a unit calibration.** The model is dimensionless,
  so this validates that there *is* a single source-independent kill threshold and
  pins its order of magnitude (~0.5 mM reacted [1O2]); it does not set the
  dimensionless `death_threshold = 10.0` to a physical mM value (no rigorous unit
  bridge exists for the LP-units death axis).
- **The O2-dependence form is only partly addressed.** The model's
  `oxygen::o2_dependent_exo_factor` (#336) is **linear** in O2,
  `(1 - dep) + dep*o2_supply`, whereas Zhu's macroscopic singlet-oxygen model makes
  the singlet-oxygen quantum yield a **saturating (Michaelis-type)** function of [O2].
  The linear form is a first-order approximation. The precise oxygen-quenching
  constant is in the full text (not the abstract), so it is **flagged here, not
  fabricated**; adding a Michaelis-form option anchored to that constant is the
  remaining #464 follow-up and needs the full-text value.

## Status

The PDT/SDT exo-ROS layer's **single source-independent kill-threshold structure**
moves from an internal design assumption to **validated against measured in-vivo
photophysics** (form + order of magnitude). The absolute `death_threshold`
magnitude and the Michaelis O2-dependence constant remain uncalibrated/data-gated on
the full-text constant.
