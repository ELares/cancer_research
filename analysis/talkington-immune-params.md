# Borrowing immune parameters from the Talkington & Kearsley 2025 ICB model (#472)

## What this is

A documented, opt-in mapping of the published immune-checkpoint-blockade (ICB)
parameters from a peer-reviewed multiscale tumor-immune model onto this repo's
spatial-immune layer (`immune_spatial`). The goal is to give our otherwise
uncalibrated immune parameters an external literature anchor (a plausible range
and a direction), **not** to swap published values into the production defaults.
The mapping math is encoded and tested in `scripts/icb_param_map.py` /
`tests/test_icb_param_map.py`.

This complements, on the immune side, the ferroptosis-core ODE cross-validation
(`analysis/ode-cross-validation.md`, #344) and the Arbatskiy biochemical
cross-validation (`analysis/arbatskiy-cross-validation.md`, #471).

## Source (verified)

- **Title**: Optimization of Immune Checkpoint Blockade via a Multiscale Model System
- **Authors**: A. M. Talkington, A. J. Kearsley
- **Venue**: Computational and Systems Oncology 5(2):e70007, 2025 (brief report)
- **DOI**: 10.1002/cso2.70007 | **PMID**: 41322398 | **PMCID**: PMC12663532
- **Access**: open access, CC BY 4.0. Parameters in Table 1 of the main text.
  There is no supplement and no public code/data (data/code "available from the
  corresponding author on reasonable request"); the values below were read from
  the open Europe PMC full-text XML.

The model is genuinely multiscale: a 3-ODE whole-tumor system (effector T cells
`T`, cancer cells `C`, non-cancer APCs `A`) **plus** an agent-based model (ABM)
that re-expresses the checkpoint efficiency as a per-contact exhaustion
probability. It has an explicit effector-cell kill term (`r2*C*T`) and an
explicit exhaustion/deactivation term (`a*(A+C)*T`, with `a = (1/n)*I*R1*R2`).

## Published values (Table 1)

| Quantity | Symbol | Value | Units | Provenance |
| --- | --- | --- | --- | --- |
| Cancer kill rate | `r2` | 1.101e-7 | day^-1 cells^-1 | cited (Kuznetsov 1994 / Talkington 2018, mouse); MC-sampled |
| Exhaustion threshold | `n` | 1.0e4 | interactions (cells) | assumed |
| ICB efficiency | `I` | 0-1 (swept) | dimensionless | varied; **I=0 = perfect blockade, I=1 = none** |
| T-cell stimulation | `r1` | 0.1245 | day^-1 | cited |
| Stimulation half-saturation | `k1` | 2.019e7 | cells | cited |
| Cancer carrying capacity | `k2` | 2.0e9 | cells | cited |
| T birth / death | `b_T` / `d_T` | 13000 / 0.0412 | cells/day, day^-1 | cited |
| Cancer birth / death | `b_C` / `d_C` | 0.18 / 1.8e-10 | day^-1 | cited |
| APC birth / death | `b_A` / `d_A` | 10000 / 0.0412 | cells/day, day^-1 | cited |
| Receptor conc. (x2) | `R1`, `R2` | sqrt(3.422e-5) ~ 5.85e-3 each | cells^-1 | cited |

**Headline result of the paper**: an 80-90% blockade (i.e. `I` ~ 0.1-0.2) is the
transition point at which progressive disease becomes less likely.

**Provenance matters**: the kill rate and most rates are *mouse-derived
literature values* (Kuznetsov, Makalkin, Taylor & Perelson, Bull. Math. Biol. 56
(1994) 295; Talkington, Dantoin & Durrett, Bull. Math. Biol. 80 (2018)), not fit
in this paper. Only `I` is swept and `n` is assumed. So they give us a *sanity
range and a direction*, not a calibration.

## Mapping onto the repo's spatial-immune model

Our model has: a per-cell immune kill probability, an `exhaustion_rate` acting
through `1/(1 + exhaustion_rate * cumulative_kills)` (#243), and a
multi-checkpoint brake over PD-1/CTLA-4/LAG-3/TIM-3, each with a
residual-after-drug (#264).

| Ours | Theirs | Ports? | Conversion / caveat |
| --- | --- | --- | --- |
| per-cell immune **kill probability** | `r2 = 1.101e-7` (bimolecular) | **No (dimensional)** | Need `p = 1 - exp(-r2 * local_effector_density * dt)`; a local effector count and a per-step time. Use their *ABM* per-contact form, not the ODE rate (inserting `r2` directly double-counts the spatial locality our model already encodes). |
| **exhaustion_rate** | `n = 1e4` interactions | **Yes (dimensionless scale)** | `exhaustion_rate ~ 1/n ~ 1e-4` per interaction. Order-of-magnitude anchor; the functional forms differ (their mass-action vs our saturating `1/(1+rate*cumulative)`). |
| **multi-checkpoint residual-after-drug** | single `I in [0,1]` | **Yes (to one aggregate)** | `I` = residual exhausting fraction after drug maps onto our per-checkpoint residual; our combined brake `1 - prod(1 - residual_i*(1-drug_eff_i))` collapses to their single `I` when one checkpoint carries it. The 80-90% optimum gives a plausible **aggregate** residual ~0.1-0.2 (do NOT assign it per-checkpoint). |
| DAMP diffusion / suppressor field / IFN-gamma return loop / 3D geometry | (no counterpart) | **Repo-specific** | The paper is well-mixed and has no spatial coupling, so these gain no external anchor here and stay uncalibrated. |

## Honest caveats

- **Well-mixed ODE != our agent-based spatial model.** `r2` and `a` are
  mass-action coefficients assuming every effector can reach every target. Our
  spatial model already encodes contact locality, so directly inserting `r2`
  would double-count the spatial limitation. This is the same effect that shrinks
  our SDT:RSL3 immune ratio from ~104:1 (2D, near DAMP saturation) to ~4:1 in 3D:
  geometry alone changes the effective ratio ~25x, which a well-mixed rate cannot
  capture (see `CALIBRATION_STATUS.md`, "Immune ICD/DAMP cascade").
- **Magnitude is not importable, direction and order-of-magnitude are.** None of
  these values were fit to our biology; the kill rate is mouse-derived. Treat
  them as a sanity range with the conversion explicitly flagged.
- **Single vs multi-checkpoint mismatch.** They model one lumped ICB knob; our
  4-checkpoint panel is finer than anything they fit. The 80-90% optimum is an
  aggregate result, not a per-axis residual.
- **Not swapped into production.** Nothing in this repo's defaults changes. This
  is a documented anchor for future calibration, consistent with the repo's
  uncalibrated-immune-layer accounting.

## What this anchors vs leaves open

- **Anchored (as a range/direction)**: the exhaustion threshold scale (`n ~ 1e4`
  interactions -> `exhaustion_rate ~ 1e-4`), the plausible aggregate
  checkpoint-residual at the blockade optimum (~0.1-0.2), and the kill-rate
  order of magnitude (~1e-7 day^-1 cell^-1, mouse).
- **Still open**: the spatial DAMP/suppressor/IFN-gamma couplings and the 3D
  geometry have no counterpart in this well-mixed model. Calibrating the absolute
  per-cell kill magnitude still needs a spatial or co-culture dataset running SDT
  and RSL3 against a shared immune readout, which does not yet exist.
