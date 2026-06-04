# Calibration Status — the one honest accounting

This is the **single place** that states, for every simulation layer in the
suite, whether its parameters are calibrated to independent data, partially
anchored, only self-consistent, or uncalibrated placeholders. It exists because
the suite grew fast (the `sim-tme-3d` 3D track added more than a dozen "realism
layers" on top of the original 2D/single-cell engine), and the calibration
caveats for those
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
| Tissue drug penetration | `drug_transport` | **Uncalibrated (illustrative)** | Exponential Krogh-cylinder approximation across 3 tissue types; all transport parameters are estimated, so the cross-tissue comparison (manuscript §6.1) is directionally informative, not quantitative. 2D helper, not part of the sim-tme-3d track. #315 adds an off-by-default ECM-tortuosity factor (`λ_eff = λ/√τ`; denser ECM shortens penetration, Netti 2000 PMID 10811131 / Provenzano 2012 PMID 22439937). `τ=1` for all shipped tissues (identity/byte-identical), an uncalibrated placeholder encoding the direction; calibrate vs ECM-density-resolved interstitial diffusivity (pancreatic desmoplasia is the natural `τ>1` candidate). |
| Photosensitizer PK | `photosensitizer_pk` | **Partially anchored** | Distribution-phase + inter-drug ROS-yield normalization closed via #203 with literature scaling; absolute cellular PK still estimated. |
| Immune ICD/DAMP cascade (2D) | `immune` | **Uncalibrated (illustrative); direction literature-anchored (#288)** | DAMP diffusion, T-cell kill rates estimated. The 104:1 SDT:RSL3 ratio (manuscript §7.2) is "a theoretical ceiling … not a quantitative prediction." Per #288, the *direction* (SDT ≫ RSL3 immune priming) is supported by verified literature (Wiernicki 35760796; Wang 34669472; Luo 35568916; Foglietta 38232641); the 2D 104:1 over-extrapolates because the model's saturating Michaelis-Menten DAMP→activation (Kd=50) is driven deep into saturation by the dense 2D kill field, while 3D volumetric dilution keeps it sub-saturating (~4:1, more consistent with the literature). The exact ratio stays uncalibratable. |

## 3D realism layers (`sim-tme-3d` track)

All off-by-default and byte-identical when off. Each was built to model a
*mechanism*; none has been fit to independent data. Where a parameter has a
literature anchor for its **structure** (e.g. vessel spacing), it is noted.
Two rows below are **not** biology realism layers and are listed only for
accounting completeness: the **vessel spatial index** (#268) is a pure
performance optimization with no parameters, and **`dose_schedule`** (#239) is a
second-wave PK module the 3D track merely consumes (manuscript §11.2 groups it
with the second wave, not the realism layers). So the table carries a few more
rows than the "more than a dozen" realism-layer count above.

| Layer | Issue | Module | Status | What would calibrate it |
|-------|-------|--------|--------|-------------------------|
| Explicit vasculature | #191 | `vasculature` | **Partially anchored** | Inter-vessel spacing (~150 µm well / ~400 µm poor) is from Vaupel; the random-uniform *placement* and Krogh λ are uncalibrated. Calibrate vs measured hypoxic fractions / vessel maps (micro-CT, pimonidazole). |
| Vessel spatial index | #268 | `vasculature` | **N/A (exact)** | Pure performance; bit-identical to brute force. No parameters. |
| Fractal vessel topology | #268 | `vasculature` | **Uncalibrated (illustrative); structure literature-anchored** | `VesselTopology::Fractal` / `place_vessels_fractal_3d` build a fractal-branching vessel tree whose chaotic-perfusion *structure* follows Baish & Jain 2000 (PMID 10919633; tumor vasculature is fractal/space-filling, D≈1.89). The branching *constants* are placeholders: BASE_ANGLE≈35°, LENGTH_RATIO 0.80, MIN_LEN 3.0, MAX_DEPTH 14, DEAD_END_PROB 0.12, n_trunks = round(points/30) (clamped [2,64]). At matched point-count it leaves a higher hypoxic fraction than random, but the magnitude is uncalibrated. Calibrate vs micro-CT vessel-tree morphometry (segment-length/branch-angle distributions, tortuosity). Off-by-default (topology=Random) so the matrix stays byte-identical. |
| 3D radial spheroid biochem | #197, #270 | `spheroid` | **Zone geometry literature-grounded (#270); biochem gradients uncalibrated** | Zone *volumes* now match the Browning 2021 (eLife, DOI 10.7554/eLife.73020) limiting structure (necrotic core 0.73 of radius → 0.39 of volume; rim begins 0.90 of radius), fixing the prior radial-threshold inversion (core ~4% → ~39% of volume). The per-zone MUFA/GSH/iron gradient *strengths* remain placeholders — calibrate vs spatially-resolved spheroid metabolomics / phenotype staining. |
| Patient-scale slab | #240 | `slab` | **Uncalibrated (illustrative)** | Krogh λ ~150 µm placeholder; the "<20 % kill at 4 mm depth" is *illustrative of the scale gap*, not a validated efficacy number. Calibrate vs depth-resolved kill in thick tissue / patient PK. |
| Slab + vasculature coupling | #272 | `slab`+`vasculature` | **Uncalibrated (illustrative)** | Inherits the slab λ and vessel placeholders; combine rule (element-wise max) is a first-order Krogh approximation. |
| Depth-graded slab phenotype | #272 | `slab` | **Uncalibrated (illustrative); structure literature-motivated** | `apply_depth_graded_cells_3d` re-assigns the slab's flat bulk mix to a layered rim→core gradient (proliferating/glycolytic at the vessel-proximal +z face, persister-like in the chronically supply-starved deep (−z) layers), the depth-axis analog of the spheroid's radial zones (#197). Thresholds are on the planar supply `exp(-depth/λ)` (NOT geometric volume fractions like the spheroid, because the slab models an *absolute* depth: a `patient_deep()` 4 mm slab is uniformly persister-like, which is correct). `SlabPhenotypeConfig::literature()` cut-points (`glycolytic_supply` 0.5 / `oxphos_supply` 0.15) are placeholder magnitudes encoding the documented direction (deep, supply-starved tissue is enriched for tolerant phenotypes), not fit to data, so read the direction, not the layer counts. Scope: (1) phenotype tracks the planar depth gradient only (internal vessels (#272 coupling) raise *delivered* drug dynamically downstream but do not reshape the chronic phenotype, a future refinement); (2) only the phenotype is depth-graded (no separate static GSH/iron gradient as in the spheroid, since the supply field already deprives deep cells dynamically). Off-by-default (the matrix never enters slab mode) ⇒ byte-identical. Calibrate vs depth-resolved phenotype/marker histology (Ki-67 proliferative rim, hypoxia/quiescence markers at depth). |
| Oxygen-dependent SDT exo-ROS | #336 | `oxygen` | **Uncalibrated (illustrative); direction clinically grounded** | `oxygen::o2_dependent_exo_factor` + sim-tme-3d `Overrides.sdt_o2_dependence` make a configurable Type II fraction of the SDT/PDT exogenous ROS yield scale with local O2 (`(1−dep)+dep·o2_supply`). This addresses the manuscript's single most contested assumption (§7.1: SDT modeled as O2-independent, an optimistic upper bound) model-side: the lead clinical sonodynamic agent SONALA-001 is Type II / O2-dependent (Sanai 2025, Sci Transl Med, PMID 41296829), with first-in-human glioma results showing only modest cell death, so the real SDT hypoxia advantage is likely far smaller than the O2-independent default. The `dependence` value itself is an uncalibrated knob (the actual Type I/Type II split for a given sonosensitizer is not measured here); calibrate vs sonosensitizer-specific O2-dependence assays (kill vs measured pO2). `sdt_o2_dependence=0` default ⇒ factor 1.0 ⇒ byte-identical. |
| 3D radial pH gradient | #190 | `ph` | **Uncalibrated (illustrative)** | The `ph_on` sim-tme-3d toggle: edge/core pH (7.4 → ~6.5) plus the iron-release and drug ion-trap sensitivities are placeholders. The manuscript's "53% pH-driven RSL3 reduction" (§7.4) rests on `ion_trap_sensitivity`, flagged there as the model's most uncertain parameter (RSL3 is a chloroacetamide with an uncharacterized pKa). Calibrate vs intracellular-RSL3-vs-pH measurement. |
| Stromal / CAF shielding | #189 | `stromal` | **Self-consistency only** | `targets.yaml: 3d_stromal_boundary_shielding`. The CAF GSH/MUFA boost reduces boundary RSL3 kill; sim-tme-3d produces ~51.5% shielding (ratio ≈ 0.485), nearly identical to 2D's 50%. The boost rates have no textbook CAF-biology source (estimates), and the target is a self-consistency check, not independent data. Upgrade to **Calibrated** vs CAF-coculture spheroid data (PMID 34373744, cited in `stromal.rs`). |
| Clonal heterogeneity + spatial expansion | #242, #266 | `clonal` | **Uncalibrated (illustrative)** | Per-subclone iron/GPX4/MUFA perturbations and `repopulation_rate` are placeholders. Calibrate vs single-cell resistance-marker distributions + lineage-tracing growth rates. |
| Cell-cell contact resistance | #270 | `contact` | **Uncalibrated (illustrative); mechanism literature-anchored** | The mechanism (dense contact ⇒ E-cadherin/Merlin/NF2 ⇒ YAP inhibition ⇒ ACSL4/TFRC down ⇒ ferroptosis resistance) is established (Wu 2019, PMID 31341276). The `ContactConfig::literature()` strengths (`lipid_strength` 0.4 / `iron_strength` 0.2, applied as `1 − strength·contact_fraction` on the durable PUFA/iron axes) are placeholder magnitudes encoding the documented direction, not fit to data. The resulting kill-rate change is **threshold-proximity-sensitive** (PUFA enters LP seeding AND autocatalytic propagation, so it scales steeply — read the direction, not the number), and the effect is **almost entirely the lipid/ACSL4 axis**; the iron/TFRC axis is near-inert under RSL3 (small additive Fenton term) and would matter more in Fenton-dominated conditions. Calibrate vs density-resolved ferroptosis-sensitivity assays (sparse vs confluent; NF2/YAP knockdown). Geometric (sphere/spheroid only — mutually-exclusive with the slab geometry); off-by-default ⇒ byte-identical. |
| Radial nutrient gradient | #270 | `nutrient` | **Uncalibrated (illustrative); direction literature-anchored** | A radial nutrient-availability field (`exp(-depth/λ)`, the O2 field's form) scales the durable antioxidant setpoint `cell.nrf2` DOWN toward the nutrient-starved core, so the core has less glucose-derived NADPH for GSH/GPX4 regeneration and is more ferroptosis-sensitive (Dixon 2012 PMID 22632970, the foundational GSH-defense mechanism; glucose metabolic reprogramming regulates ferroptosis, PMID 42190602). `NutrientConfig::literature()` `antioxidant_strength` 0.3 is a placeholder magnitude encoding ONE documented direction; the NET effect is genuinely context-dependent (energy stress also activates AMPK, which INHIBITS ferroptosis; glutaminolysis is REQUIRED for some ferroptosis, PMID 30581146). Calibrate vs depth-resolved glucose/NADPH/GSH in spheroids. Geometric (no RNG); off-by-default ⇒ byte-identical. |
| Drug-tolerant persisters | #241, #262 | `persister` | **Uncalibrated (illustrative)** | Acquisition/reversion rates are step-level placeholders (no published step-level kinetics). Target: Hangauer-style multi-cycle drug screen (fraction-surviving vs cycles) + off-drug reversion half-time — see `parameter_provenance.md`. |
| T-cell exhaustion | #243 | `immune_spatial` | **Uncalibrated (illustrative)** | `exhaustion_rate` placeholder. Calibrate vs longitudinal TIL cytotoxicity / exhaustion-marker time-courses. |
| Treg/MDSC suppressor field | #264 P2 | `immune_spatial` | **Uncalibrated (illustrative)** | Suppressor strength + diffusion length placeholders; niche sources heuristic. Calibrate vs spatial Treg/MDSC density + local kill suppression. |
| Multi-checkpoint brake | #264 P3 | `immune_spatial` | **Uncalibrated (illustrative)** | Per-axis (PD-1/CTLA-4/LAG-3/TIM-3) brake residuals + drug efficacies placeholders. Calibrate vs blockade-response dose-effect data per axis. |
| DC subset mix (cDC1/cDC2) | #264 P4 | `immune_spatial` | **Uncalibrated (illustrative); direction literature-anchored** | `DcSubsetConfig` collapses to one uniform anti-tumor priming-efficiency scalar (`cdc1_fraction·cdc1_eff + (1−cdc1_fraction)·cdc2_eff`) multiplied into the immune kill probability. A cDC1-poor tumor primes CD8 killing less efficiently (cDC1/Batf3 cross-presenting DCs are the rare, critical anti-tumor APCs: Broz et al., Cancer Cell 2014, PMID 25446897). The `literature()` mix (10% cDC1, cdc2_eff 0.3 ⇒ efficiency 0.37) is an UNCALIBRATED placeholder; only the direction (cDC1-poor ⇒ weaker priming ⇒ fewer immune kills) is claimed. Calibrate vs tumor DC-subset abundance + response data. Off-by-default (balanced ⇒ efficiency 1.0) byte-identical. |
| Time-varying dose schedules | #239 | `dose_schedule` | **Partially anchored** | Schedule *shapes* (bolus/infusion/multidose) are standard PK forms; the `FromPk` path inherits `tumor_pk`'s order-of-magnitude constants. |

### 3D self-consistency targets (in `targets.yaml`)

The three 3D targets (`3d_rsl3_o2_collapse_ratio`, `3d_immune_sdt_dominates`,
`3d_stromal_boundary_shielding`) are **self-consistency** — they regression-guard
the 3D model's own predictions, not measured values. `3d_validation_report.md`
records that, like-for-like, the 2D engine produces a **larger** immune ratio
than 3D (104:1 vs **4:1**) and a slightly more complete hypoxia collapse — i.e.
the 3D numbers are not "stronger," and the immune ratio in particular is
sensitive to grid size / volumetric DAMP dilution.

**Immune target — direction now literature-anchored (#288).** A #288 literature
review (every PMID verified against PubMed) established that the SDT:RSL3 immune
ratio's *exact value* cannot be calibrated — no published study runs SDT and an
RSL3-only arm against a shared immune readout, and the previously-named
"Nguyen 2019" source was a **phantom citation** (no such applicable study exists;
removed from `targets.yaml`). But the *direction* (ratio > 1) and *single-digit
magnitude-class* ARE supported by verified primary literature: the
RSL3/GPX4-inhibitor denominator is non- to weakly-immunogenic (Wiernicki 2022
PMID 35760796, 0% vaccination protection; with the genuine caveat that early
ferroptosis can be immunogenic — Efimova 2020 PMID 33188036), SDT drives
DC maturation + CD8 infiltration over controls (Wang 2021 PMID 34669472;
Luo 2022 PMID 35568916; pancreatic-spheroid SDT-ICD Foglietta 2024 PMID 38232641;
the precise fold-changes are figure-level, not abstract-level). The 2D 104:1
over-extrapolates for a model-internal reason: the immune layer's saturating
Michaelis-Menten DAMP→activation (Kd=50) is driven deep into saturation by the
dense 2D kill field, whereas 3D volumetric DAMP dilution keeps it sub-saturating
(the more literature-consistent ~4:1). So `3d_immune_sdt_dominates`
stays self-consistency for its *value* but is now literature-anchored for its
*direction and magnitude bound*.

## Calibration roadmap (priority order)

Calibration is **data-gated**: each item needs an independent measurement the
repo does not currently hold. Listed by leverage (how much it would change a
load-bearing claim) × tractability (how obtainable the data is):

1. **Immune coupling (2D + 3D).** ~~Needs a published SDT-immune-priming dataset.~~
   **DONE as far as the literature allows (#288):** the direction + magnitude-class
   are now literature-anchored (see the immune-target note above). What remains
   genuinely data-blocked is a *precise* calibrated value, which requires a study
   running SDT and RSL3 against a **shared** immune readout — none exists. The
   honest residual is that the exact ratio (3D ~4:1) is uncalibratable, not that
   the claim is ungrounded; the keystone direction is supported.
2. **Persister kinetics.** Direction is well-supported; only the rates are
   guessed. A single Hangauer-style multi-cycle screen would fit
   acquisition/reversion/`max_fraction`. Tractable, well-scoped.
3. **Vasculature → hypoxic fraction.** Spacing is anchored; validating that the
   field reproduces a *measured* hypoxic fraction (pimonidazole / micro-CT) would
   move it from "partially anchored" to calibrated and underwrite the #272
   patient-scale depth-collapse story.
4. **Slab depth-collapse magnitude.** The "<20 % at depth" number drives the
   in-vitro-vs-patient narrative; depth-resolved kill in thick tissue would make
   it quantitative rather than illustrative. Remaining `slab` refinements tracked
   under #240/#272 (all data-gated or opportunistic, none blocking): the
   depth-collapse calibration above; the depth-graded slab phenotype's supply
   cut-points (validate vs depth-resolved proliferation/quiescence histology); a
   *multi-face* supply (the current planar gradient is single-face +z; a real
   tumor chunk is perfused from several boundaries); and folding the hard-coded
   `VIRTUAL_TUMOR_MM` / `depth_offset_mm` placement constants into a derived
   geometry once a patient-scale tumor model is available. The slab+vasculature
   coupling (#272) and depth-graded phenotype (#272) have shipped; the rest are
   left until the matching data exists rather than guessed.
5. **Spheroid radial gradients, clonal perturbations, suppressor/checkpoint/
   exhaustion rates.** Lower individual leverage (each is one knob in a
   composite); calibrate opportunistically as spatially-resolved or longitudinal
   datasets surface.

Until an item lands, treat the corresponding layer's output as **directional /
illustrative**, exactly as the module doc-comments and the manuscript's §8.4
already state. When a target is upgraded to `calibration` in `targets.yaml`
(add `source_pmid` + the measured value), update its row here in the same PR.
