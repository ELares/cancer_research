# Cross-validating ferroptosis-core against the Arbatskiy 2024 73-ODE model (#471)

## What this is

A documented, structural cross-validation of this repo's single-cell ferroptosis
core ODE against an independently-built, published 73-ODE / 93-species
ferroptosis model, plus an honest accounting of which of our placeholder rate
constants the published parameter set could constrain and which stay
repo-specific. This is a **structural / topological** cross-check and a pointer
to borrowable (CC BY) published constants; it is **not** a numeric transplant
(see "Why structural, not numeric" below).

It complements `analysis/ode-cross-validation.md` (#344, our bistable switch vs
Co 2024 / Seidel 2026 / Konstorum 2020) and `analysis/sbml-export.md` (#351, our
core exported to SBML), and on the immune side
`analysis/talkington-immune-params.md` (#472).

## Source (verified)

- **Title**: A Systems Biology Approach Towards a Comprehensive Understanding of Ferroptosis
- **Authors**: M. Arbatskiy, D. Balandin, I. Akberdin, A. Churov
- **Venue**: International Journal of Molecular Sciences 25(21):11782, 2024 (MDPI)
- **DOI**: 10.3390/ijms252111782 | **PMID**: 39519341 | **PMCID**: PMC11546516
- **Access / license**: CC BY 4.0 (open, redistributable). The model is **73
  differential equations / 93 species** built in BioUML (2023.1) as six modules:
  Fenton + Haber-Weiss, iron metabolism, lipid synthesis, lipid peroxidation,
  pentose phosphate pathway, antioxidant system.
- **Where the numbers live**: the per-module kinetic constants are in
  **Supplementary Tables S1-S7** (CC BY; landing page
  `https://www.mdpi.com/article/10.3390/ijms252111782/s1`). The runnable
  SBML/BioUML model itself is **"available on request"** (no BioModels deposit,
  no GitHub/Zenodo link in the paper). So the citable, reusable artifact is the
  CC BY supplementary kinetic tables; the executable model requires contacting
  the authors.

## Overlapping mechanisms (structure confirmed)

| Our mechanism | Arbatskiy module | Structural correspondence | Numbers |
| --- | --- | --- | --- |
| Fenton / Haber-Weiss iron-radical | Module 2.2 (Table S1) | Fe2+ + H2O2 -> Fe3+ + .OH, two-constant scheme + Haber-Weiss as a second .OH source; a fast .OH phase (~120 min) then a slow Fe2+-regeneration-limited phase | S1 |
| Iron metabolism (TfR1 / ferritin / ferritinophagy / labile iron) | Module 2.3 (Table S2) | TF, TFR1, STEAP3/4, DMT1, ferritin + mitochondrial ferritin, ceruloplasmin, MCU; ferritinophagy explicitly releases iron into the labile pool | S2 |
| Lipid peroxidation initiation + propagation | Module 2.5 (Table S4) | RO2. (peroxyl-radical) chain propagation; ALOX15 enzymatic initiation; alpha-tocopherol as the chain-breaking radical trap | S4 |
| GPX4 / GSH / GSSG turnover | Modules 2.6-2.7 (Tables S6, S7) | GPX4 reduces L-OOH using 2 GSH -> GSSG; GR regenerates GSH from GSSG using NADPH; Xc- / cystine -> cysteine -> GSH synthesis | S6/S7 |
| NADPH regeneration | Module 2.6 (Table S5, PPP) | NADPH via G6PD / pentose phosphate pathway, consumed by GR | S5 |

The two models share the load-bearing topology: a GSH/GPX4-gated peroxyl-radical
chain fed by Fenton iron and buffered by NADPH-regenerated glutathione, with
ferritinophagy as a labile-iron source. That an independently-constructed 73-ODE
model reaches the same recover-or-collapse structure our core encodes is the
structural cross-validation #471 asks for.

## Which of our constants the published set could constrain

| Our constant | Constrainable by Arbatskiy? | Table |
| --- | --- | --- |
| `lp_rate` (peroxidation initiation) | Yes (ALOX15-driven initiation) | S4 |
| `lp_propagation` (autocatalytic chain) | Yes (RO2. propagation, the heart of their chain) | S4 |
| `gpx4_rate` (GPX4 detox) | Yes (GPX4 catalytic turnover on L-OOH) | S6 |
| `gsh_scav_efficiency` (GSH scavenging) | Yes (GSH->GSSG stoichiometry/rate + GR regen) | S6/S7 |
| iron Fenton rate | Yes (explicit two-constant Fenton + Haber-Weiss) | S1 |
| `ferritinophagy_release` (#340) | Partially (ferritinophagy -> labile-iron release is explicit) | S2 |
| `gpx4_degradation_by_ros` | Weak/indirect (they inactivate GPX4 *via GSH depletion*, not a direct ROS-degradation term) | S6 |

## Which of our constants stay repo-specific

- **`dhodh_rate` (#338, DHODH/CoQ10)**: NOT in the model. No DHODH /
  mitochondrial GPX4-independent arm. Stays repo-specific.
- **`gch1_rate` (#338, GCH1/BH4)**: NOT in the model. Their only
  GPX4-independent radical trap is **alpha-tocopherol** (vitamin E), a *different*
  axis we do not currently parameterize, not a substitute for GCH1/BH4.
- **FSP1/CoQ (the FSP1 repair pathway)**: NOT modeled at all, so no FSP1-related
  constant is borrowable.
- The newer lipid-remodeling escape axes (#339 ether-lipid / MBOAT, #446 ALOX
  isoform mix, #444 ACSL4 status) and all spatial/TME layers are outside this
  single-cell biochemical model's scope.

## Why structural, not numeric (transplant caveats)

A direct 1:1 numeric transplant is blocked, so the honest deliverable is the
structural cross-check plus the pointer to the CC BY tables:

1. **Compartmental, molar state variables.** Their species are molar
   concentrations (uM/mM) with explicit cytoplasmic vs mitochondrial
   compartments; our core is a non-compartmental single-cell ODE in
   dimensionless pools. Second-order Fenton/peroxidation constants (units like
   M^-1 s^-1) cannot be dropped in without converting through our concentration
   scale and timestep.
2. **Different Fenton form.** Their lumped two-constant Fenton differs
   structurally from our single `iron Fenton rate`; you would refit, not copy.
3. **No clean analog for `gpx4_degradation_by_ros`.** They deplete GPX4
   *function* via GSH, not degrade the enzyme.
4. **Supplementary tables, not machine-fetchable here.** The CC BY S1-S7 numbers
   are reusable with citation but must be read from the rendered MDPI supplement
   (or the on-request BioUML/SBML model); MDPI returns 403 to automated fetches,
   so they are not committed as a derived artifact (consistent with the repo's
   offline contract: we do not vendor data we cannot fetch reproducibly).

## A free timescale cross-check (no units needed)

The paper reports a **~120 min fast .OH-generation phase** in its Fenton module
before the slow Fe2+-regeneration-limited phase. That is a units-free, directional
cross-check for our Fenton + ferritinophagy ramp (`ferritinophagy_tau`, #340):
our labile-iron release should ramp on a comparable order (tens of minutes to a
couple of hours), not instantaneously and not over many hours. This is a sanity
bound on the ramp shape, not a calibration of its rate.

## Bottom line

Treat Arbatskiy 2024 as (a) an independent **structural** confirmation that our
GSH/GPX4 -> L-OOH topology, Fenton + ferritinophagy iron flux, and bistable
collapse match a separately-built 73-ODE model; (b) a **pointer to borrowable CC
BY constants** (S1-S7) for our lipid-peroxidation, GSH/GPX4, Fenton, and
ferritinophagy rates, should a future calibration PR convert through the
concentration/timestep scale; and (c) a **scope marker** for what it does *not*
constrain (DHODH, GCH1/BH4, FSP1), which therefore stay repo-specific. No
production default changes; this is documentation and a calibration pointer.
