# Calibration Status — the one honest accounting

This is the **single place** that states, for every simulation layer in the
suite, whether its parameters are calibrated to independent data, partially
anchored, only self-consistent, or uncalibrated placeholders. It exists because
the suite grew fast (the `sim-tme-3d` 3D track added ~11 "realism layers" on top
of the original 2D/single-cell engine), and the calibration caveats for those
layers were scattered across module doc-comments, two READMEs, and
`3d_validation_report.md`. A reader deserves one table, not a scavenger hunt.

**The headline, stated plainly:** the simulation suite is **broad but mostly
uncalibrated**. A small core is anchored to literature or sensitivity-tested;
the rest is *mechanistic scaffolding* — it shows the *shape* of an effect
(direction, qualitative behavior, spatial structure), not a *quantitative*
prediction. Breadth of coverage is **not** evidence of depth (see the repo's
`CLAUDE.md` guiding principles). None of the uncalibrated 3D layers feed the
manuscript's quantitative claims; the manuscript's numbers come from the 2D
engine and are themselves labeled order-of-magnitude (see
`article/drafts/v1.md` §8.4).

This document is **handwritten interpretation**, kept in sync by hand. The
machine-checked authorities remain:
- `parameter_provenance.md` — per-parameter table for the **core** engine
  (`Params`, `SpatialParams`, `ImmuneParams`, `RecoveryRates`) + RSL3 PK +
  persister + photosensitizer PK, with literature citations and a `Grounded?`
  column.
- `targets.yaml` — the 8 machine-evaluated calibration/self-consistency targets
  (read by `calibrate.py`).
- `3d_validation_report.md` — what the 3D capstone's Q1–Q4 self-consistency
  checks actually produce, and the 2D-vs-3D magnitude differences.
- Source code (`ferroptosis-core/src/*.rs`) — authoritative for default values.

## Status legend

| Tier | Meaning |
|------|---------|
| **Calibrated** | At least one key parameter fit/anchored to an independent published measurement. |
| **Partially anchored** | Structure or a subset of parameters has a literature basis, but the headline magnitude is not fit to data. |
| **Self-consistency only** | A regression-guard that the model reproduces its *own* hard-coded physics; not an independent calibration. |
| **Uncalibrated (illustrative)** | Placeholder parameters chosen for plausible behavior. The mechanism/direction is the claim; the magnitude is not. |

## Core engine and 2D spatial work

| Layer | Module(s) | Status | Notes / what would calibrate it |
|-------|-----------|--------|--------------------------------|
| Core ferroptosis biochemistry | `biochem`, `params` | **Partially anchored** | Mix of literature-grounded and estimated rate constants (see `parameter_provenance.md` `Grounded?` column). Directionally robust: ±50% sensitivity held in 22/22 conditions (manuscript §5, Chapter sensitivity analysis). Magnitudes are estimates. |
| GPX4 recovery kinetics | `params` / `RecoveryRates` | **Calibrated** | `targets.yaml: gpx4_recovery_rate` — consensus from multiple ferroptosis studies. |
| MUFA / PUFA protection | `biochem` | **Calibrated (composite)** | `targets.yaml: mufa_protection_factor` (18.6×) — model-derived composite of Dixon/Park (unpublished 2025 submission, no PMID) + Tesfay 2019 (PMID 31270077, SCD1 protects ovarian cancer from ferroptosis). Not a single measurement. |
| Persister FSP1/HDAC suppression | `biochem` | **Calibrated** | `targets.yaml: fsp1_hdac_persister` — Higuchi et al., Science Advances 2026 (PMID 41481741). |
| PDT optical attenuation | `physics`, `oxygen` | **Self-consistency only** | `targets.yaml: pdt_depth_attenuation` verifies the hard-coded `pdt_mu_eff` reproduces Jacques 2013 optics — checks the code, not biology. |
| SDT acoustic attenuation | `physics` | **Self-consistency only** | `targets.yaml: sdt_depth_attenuation` verifies `sdt_alpha` vs Cobbold 2007. |
| RSL3 pharmacokinetics | `tumor_pk` | **Uncalibrated (illustrative)** | "Order-of-magnitude estimates, not clinical measurements" (`parameter_provenance.md`). |
| Photosensitizer PK | `photosensitizer_pk` | **Partially anchored** | Distribution-phase + inter-drug ROS-yield normalization closed via #203 with literature scaling; absolute cellular PK still estimated. |
| Immune ICD/DAMP cascade (2D) | `immune` | **Uncalibrated (illustrative)** | DAMP diffusion, T-cell kill rates estimated. The 104:1 SDT:RSL3 ratio (manuscript §7.2) is explicitly "a theoretical ceiling … not a quantitative prediction." |

## 3D realism layers (`sim-tme-3d` track)

All off-by-default and byte-identical when off. Each was built to model a
*mechanism*; none has been fit to independent data. Where a parameter has a
literature anchor for its **structure** (e.g. vessel spacing), it is noted.

| Layer | Issue | Module | Status | What would calibrate it |
|-------|-------|--------|--------|-------------------------|
| Explicit vasculature | #191 | `vasculature` | **Partially anchored** | Inter-vessel spacing (~150 µm well / ~400 µm poor) is from Vaupel; the random-uniform *placement* and Krogh λ are uncalibrated. Calibrate vs measured hypoxic fractions / vessel maps (micro-CT, pimonidazole). |
| Vessel spatial index | #268 | `vasculature` | **N/A (exact)** | Pure performance; bit-identical to brute force. No parameters. |
| 3D radial spheroid biochem | #197 | `spheroid` | **Uncalibrated (illustrative)** | Radial phenotype split + MUFA/GSH/iron gradients are placeholder shapes. Calibrate vs spatially-resolved spheroid metabolomics / phenotype staining. |
| Patient-scale slab | #240 | `slab` | **Uncalibrated (illustrative)** | Krogh λ ~150 µm placeholder; the "<20 % kill at 4 mm depth" is *illustrative of the scale gap*, not a validated efficacy number. Calibrate vs depth-resolved kill in thick tissue / patient PK. |
| Slab + vasculature coupling | #272 | `slab`+`vasculature` | **Uncalibrated (illustrative)** | Inherits the slab λ and vessel placeholders; combine rule (element-wise max) is a first-order Krogh approximation. |
| Clonal heterogeneity + spatial expansion | #242, #266 | `clonal` | **Uncalibrated (illustrative)** | Per-subclone iron/GPX4/MUFA perturbations and `repopulation_rate` are placeholders. Calibrate vs single-cell resistance-marker distributions + lineage-tracing growth rates. |
| Drug-tolerant persisters | #241, #262 | `persister` | **Uncalibrated (illustrative)** | Acquisition/reversion rates are step-level placeholders (no published step-level kinetics). Target: Hangauer-style multi-cycle drug screen (fraction-surviving vs cycles) + off-drug reversion half-time — see `parameter_provenance.md`. |
| T-cell exhaustion | #243 | `immune_spatial` | **Uncalibrated (illustrative)** | `exhaustion_rate` placeholder. Calibrate vs longitudinal TIL cytotoxicity / exhaustion-marker time-courses. |
| Treg/MDSC suppressor field | #264 P2 | `immune_spatial` | **Uncalibrated (illustrative)** | Suppressor strength + diffusion length placeholders; niche sources heuristic. Calibrate vs spatial Treg/MDSC density + local kill suppression. |
| Multi-checkpoint brake | #264 P3 | `immune_spatial` | **Uncalibrated (illustrative)** | Per-axis (PD-1/CTLA-4/LAG-3/TIM-3) brake residuals + drug efficacies placeholders. Calibrate vs blockade-response dose-effect data per axis. |
| Time-varying dose schedules | #239 | `dose_schedule` | **Partially anchored** | Schedule *shapes* (bolus/infusion/multidose) are standard PK forms; the `FromPk` path inherits `tumor_pk`'s order-of-magnitude constants. |

### 3D self-consistency targets (in `targets.yaml`)

The three 3D targets (`3d_rsl3_o2_collapse_ratio`, `3d_immune_sdt_dominates`,
`3d_stromal_boundary_shielding`) are **self-consistency** — they regression-guard
the 3D model's own predictions, not measured values. `3d_validation_report.md`
records that, like-for-like, the 2D engine produces a **larger** immune ratio
than 3D (104:1 vs **4:1**) and a slightly more complete hypoxia collapse — i.e.
the 3D numbers are not "stronger," and the immune ratio in particular is
sensitive to grid size / volumetric DAMP dilution. Each target's
`source_description` names the published dataset that would upgrade it from
self-consistency to calibration.

## Calibration roadmap (priority order)

Calibration is **data-gated**: each item needs an independent measurement the
repo does not currently hold. Listed by leverage (how much it would change a
load-bearing claim) × tractability (how obtainable the data is):

1. **Immune coupling (2D + 3D).** The SDT:RSL3 ICD ratio is the manuscript's
   keystone simulation finding *and* the most volatile (104:1 → 4:1 across
   dimensionality). Highest leverage. Needs published spheroid SDT-immune-priming
   data (DAMP-density vs T-cell-infiltrate) — `3d_immune_sdt_dominates`'s
   `source_description` flags Nguyen 2019 or equivalent for a curator to pin.
2. **Persister kinetics.** Direction is well-supported; only the rates are
   guessed. A single Hangauer-style multi-cycle screen would fit
   acquisition/reversion/`max_fraction`. Tractable, well-scoped.
3. **Vasculature → hypoxic fraction.** Spacing is anchored; validating that the
   field reproduces a *measured* hypoxic fraction (pimonidazole / micro-CT) would
   move it from "partially anchored" to calibrated and underwrite the #272
   patient-scale depth-collapse story.
4. **Slab depth-collapse magnitude.** The "<20 % at depth" number drives the
   in-vitro-vs-patient narrative; depth-resolved kill in thick tissue would make
   it quantitative rather than illustrative.
5. **Spheroid radial gradients, clonal perturbations, suppressor/checkpoint/
   exhaustion rates.** Lower individual leverage (each is one knob in a
   composite); calibrate opportunistically as spatially-resolved or longitudinal
   datasets surface.

Until an item lands, treat the corresponding layer's output as **directional /
illustrative**, exactly as the module doc-comments and the manuscript's §8.4
already state. When a target is upgraded to `calibration` in `targets.yaml`
(add `source_pmid` + the measured value), update its row here in the same PR.
