# Parameter Provenance

Every simulation parameter, its default value, source, and whether it is experimentally grounded or assumed.

## Core Biochemistry (`Params`)

| Parameter | Default | Source | Grounded? | Sensitivity |
|-----------|---------|--------|-----------|-------------|
| `fenton_rate` | 0.02 | Kakhlon & Cabantchik, Free Radic Biol Med 2002 | Estimated | Moderate — controls basal ROS from labile iron |
| `gsh_scav_efficiency` | 0.5 | Michaelis-Menten model; empirical fit | Assumed | Moderate — fraction of ROS quenched by GSH per collision |
| `gsh_km` | 2.0 mM | Michaelis-Menten saturation constant; empirical | Assumed | Low — GSH-ROS binding saturation threshold |
| `nrf2_gsh_rate` | 0.025 | Dodson et al., Free Radic Biol Med 2019 | Estimated | Moderate — NRF2-driven GSH resynthesis |
| `lp_rate` | 0.06 | Yang et al., Cell 2016 (PUFA + ROS coupling) | Estimated | High — direct lipid peroxidation from unscavenged ROS |
| `lp_propagation` | 0.10 | Porter et al., Chem Rev 2005 (lipid cascade kinetics) | Estimated | **Critical** — autocatalytic bistable switch gate |
| `gpx4_rate` | 0.30 | Ursini et al., Free Radic Biol Med 1995 | Estimated | Moderate — GPX4 repair efficiency |
| `fsp1_rate` | 0.08 | Bersuker et al., Nature 2019; Mao et al., Nature 2021 | Estimated | **Critical** — FSP1/DHODH CoQ10 pathway; persister phenotype has 0.15 mean |
| `scd_mufa_rate` | 0.0 (2D) / 0.01 (in vivo) | Dixon/Park, Cancer Res 2025; Tesfay et al., Cancer Res 2019 | Grounded | **Critical** — in-vivo MUFA accumulation; steady-state derived |
| `scd_mufa_max` | 0.0 (2D) / 0.50 (in vivo) | Dixon/Park 2025 lipidomics (40-60% range) | Grounded | **Critical** — maximum PUFA displacement fraction |
| `initial_mufa_protection` | 0.0 (2D) / 0.40 (in vivo) | Derived: M_ss = rate×max/(rate+decay×max) | Derived | **Critical** — pre-accumulated MUFA in established tumors |
| `scd_mufa_decay` | 0.0 (2D) / 0.005 (in vivo) | Membrane lipid half-life ~24-48h | Estimated | Moderate — natural phospholipid turnover |
| `gpx4_degradation_by_ros` | 0.002 | Mechanistic assumption | Assumed | Low — GPX4 protein degradation under oxidative stress |
| `gpx4_nrf2_upregulation` | 0.008 | Mechanistic assumption | Assumed | Moderate — NRF2-driven GPX4 mRNA/protein upregulation |
| `sdt_ros` | 5.0 | Literature-derived for ~1 MHz ultrasound | Estimated | **Critical** — exogenous ROS peak from SDT |
| `pdt_ros` | 5.0 | Matched to SDT for controlled comparison | Estimated | **Critical** — exogenous ROS peak from PDT |
| `rsl3_gpx4_inhib` | 0.92 | Literature IC50 data; pharmacokinetic models | Estimated | High — 92% GPX4 inhibition by RSL3 |
| `gsh_max` | 12.0 mM | Forman et al., Free Radic Biol Med 2009 | Grounded | Moderate — maximum intracellular GSH pool |
| `gpx4_nrf2_target_multiplier` | 1.0 | Scaling factor; unit default | Assumed | Low |
| `death_threshold` | 10.0 | Bistable threshold; empirical fit | Assumed | **Critical** — lipid peroxidation level triggering cell death |

## Spatial/Physics (`SpatialParams`)

| Parameter | Default | Source | Grounded? | Sensitivity |
|-----------|---------|--------|-----------|-------------|
| `cell_size_um` | 20.0 | Typical tumor cell diameter | Grounded | Low — grid resolution |
| `iron_diffusion_coeff` | 281.0 µm²/s | Estimate (tortuosity-reduced tissue value; free aqueous Fe²⁺ ≈ 700 µm²/s, scaled ~2.5× down for tissue). **Not** from Jacques 2013 — that optics reference applies to `pdt_mu_eff`, and was previously mis-attributed here | Assumed | Low — bystander iron diffusion |
| `iron_release_per_death` | 2.0 µM | Mechanistic estimate | Assumed | Low — iron released per dead cell |
| `pdt_mu_eff` | 0.31 /mm | Jacques SL, Phys Med Biol 2013 (630nm red light) | **Grounded** | **Critical** — PDT penetration depth (δ ≈ 3.2mm) |
| `pdt_i0` | 1.0 | Relative units | N/A | Low — incident fluence normalization |
| `sdt_alpha` | 0.7 dB/cm/MHz | Cobbold RSC, Foundations of Biomedical Ultrasound 2007 | **Grounded** | High — acoustic attenuation in tissue |
| `sdt_freq_mhz` | 1.0 | Typical SDT frequency | Grounded | Moderate — operating frequency |
| `sdt_i0` | 1.0 | Relative units | N/A | Low — incident intensity normalization |
| `neighbor_iron_fraction` | 0.1 | Mechanistic estimate (8-neighborhood) | Assumed | Low |
| `photosensitizer` | `Uniform(1.0)` (default) | `Photosensitizer::Porfimer` carries `t_half_h` (504 h, Bellnier 2006), `t_distribution_h` (default 0 for backwards-compat; ~24-48 h literature-reported, Bellnier 2006), and `phi_so2_relative` (default 1.0 = porfimer baseline; absolute porfimer phi_so2 ≈ 0.65 in solution, Wilson & Patterson 2008, Phys Med Biol 53(9):R61-109). `Params::pdt_ros = 5.0` is calibrated to porfimer at peak (yield = 1.0); other drug variants set `phi_so2_relative` to absolute_phi_so2 / 0.65 so the calibration carries through. | Estimated (porfimer t½, t_dist); Grounded (porfimer absolute phi_so2); N/A (`Uniform` default) | Low at default values; scales linearly via `Photosensitizer::yield_at` |
| `t_drug_light_interval_h` | 0.0 | Operational parameter. Hours from photosensitizer post-distribution peak to light, passed to `Photosensitizer::yield_at`. With `Porfimer::t_distribution_h > 0`, the model holds drug at peak for the first `t_distribution_h` hours then decays — so a clinical DLI from injection can be passed directly, and the math correctly returns peak concentration during the absorption phase. Earlier comments said "clinical DLI ≈ distribution_phase + this" — that caveat is obsolete now that distribution-phase is modeled explicitly. | N/A | High — at DLI < t_distribution_h has no effect; at DLI > t_distribution_h, drug decays from peak |

## Immune Cascade (`ImmuneParams`)

| Parameter | Default | Source | Grounded? | Sensitivity |
|-----------|---------|--------|-----------|-------------|
| `damp_per_lp` | 1.0 | Krysko et al., Nat Rev Cancer 2012 | Estimated | Moderate — DAMP signal proportional to LP at death |
| `dc_activation_kd` | 50.0 | Empirical (no direct measurement) | Assumed | Moderate — half-maximal DC activation threshold |
| `dc_maturation_rate` | 0.6 | Mechanistic estimate | Assumed | Low |
| `tcell_priming_rate` | 10.0 | Mechanistic estimate | Assumed | Low |
| `tcell_kill_rate` | 3.0 | Mechanistic estimate | Assumed | Low |
| `pd1_brake` | 0.7 | Clinical estimate (70% suppression) | Estimated | Moderate |
| `anti_pd1_efficacy` | 0.8 | Clinical estimate (80% brake removal) | Estimated | Moderate |

## Recovery Rates (`RecoveryRates`)

| Parameter | Default | Source | Grounded? | Sensitivity |
|-----------|---------|--------|-----------|-------------|
| `fsp1_half_recovery_days` | 7.0 | Epigenetic recovery kinetics; slowest pathway | Estimated | High — determines FSP1 restoration timing |
| `gpx4_half_recovery_days` | 3.0 | Transcriptional recovery kinetics | Estimated | **Critical** — controls RSL3 window closure (day 3 claim) |
| `nrf2_half_recovery_days` | 5.0 | Activation kinetics | Estimated | Moderate |
| `gsh_half_recovery_days` | 1.0 | Metabolic recovery; fastest pathway | Estimated | Low |

## Summary

- **Grounded** (value derived from specific published measurement): `gsh_max`, `pdt_mu_eff`, `sdt_alpha`, `cell_size_um`, `scd_mufa_rate`, `scd_mufa_max`
- **Estimated** (informed by literature ranges but not directly calibrated): most biochemistry rates
- **Assumed** (mechanistic placeholder with no direct data): `gsh_km`, `gpx4_degradation_by_ros`, `gpx4_nrf2_upregulation`, `death_threshold`, `iron_diffusion_coeff` (tortuosity-reduced estimate; previously mis-cited to Jacques 2013), immune cascade parameters
- **Derived** (calculated from other parameters): `initial_mufa_protection`

## RSL3 pharmacokinetics: known uncalibrated

`tumor_pk::TumorPKParams` and the Krogh penetration model in `drug_transport` use RSL3-like parameters (e.g., plasma t½ ≈ 30 min, `k_uptake_bulk`, `km_uptake`) that are **order-of-magnitude estimates from chemical-probe literature, not clinical measurements.** RSL3 has no published clinical PK profile — it is widely cited as a research probe with poor pharmacokinetics, not a development candidate (e.g., review in Yang et al., Nature 2023, on ferroptosis therapeutics). Sensitivity of manuscript claims to these values is bounded by the protection-factor range reported in Chapter 8.2 (4.8×–27×) — qualitative conclusions about tumor-PK barriers are robust, but absolute kill rates should be read as approximate. Anchoring PK parameters either to a clinically published ferroptosis inducer (e.g., IKE) or to a non-RSL3 reference compound is tracked in #316 (the re-derivation is data-gated on published clinical PK).

## Drug-tolerant persister kinetics (`PersisterConfig`): known uncalibrated

The persister model (#241) acquires/loses an epigenetic drug-tolerant fraction per cell. As of #262 the per-step update is a **competing-rate** integrator (both rates always active):

`frac += acquisition_rate · drug · (max_fraction − frac) − reversion_rate · frac`

so under sustained sub-saturating drug the fraction relaxes to the equilibrium `f* = acq·drug·max / (rev + acq·drug)` (below the cap; e.g. at the `enabled()` defaults a drug intensity of 0.3 settles at `f* = 0.30`, and even saturating drug settles at `0.53 < 0.80`), rather than ratcheting monotonically to `max_fraction`.

| Parameter | `enabled()` | Source | Grounded? | Target / regeneration |
|-----------|-------------|--------|-----------|-----------------------|
| `acquisition_rate` | 0.02 / step | Hangauer et al., Nature 2017 (drug-tolerant persisters arise over days–weeks of drug) | Assumed (qualitative direction only) | Fit to a persister-fraction-vs-time-under-drug curve |
| `reversion_rate` | 0.01 / step | Hangauer 2017 (tolerance is reversible on drug withdrawal) | Assumed | Fit jointly with `acquisition_rate` to the on-drug equilibrium + off-drug decay rate |
| `max_fraction` | 0.80 | Mechanistic cap | Assumed | The plateau of the fraction-surviving-vs-dose-cycles curve |
| `gpx4_resistance` | (enabled) | Hangauer 2017 (persisters are GPX4-dependent; the RSL3-specific resistance axis) | Assumed | Fit to the persister-vs-parental RSL3 IC50 shift |
| `mufa_boost_per_step` | (enabled) | Tsoi et al., Cancer Cell 2018 (dedifferentiated/MUFA-enriched state) | Assumed | Fit to the MUFA/PUFA lipidomic shift in the tolerant state |

**Calibration target.** A Hangauer-style multi-cycle drug screen gives fraction-surviving vs. number of dose cycles; fitting `acquisition_rate` / `reversion_rate` / `max_fraction` to that curve (plus the off-drug reversion half-time) would replace the current step-level placeholders. No step-level rates are published — the literature gives direction (persisters exist, are GPX4-dependent, are reversible), not kinetics — so these stay **Assumed** until such a fit lands. Treat the persister-fraction magnitudes as illustrative; the directional claims (sustained drug builds a tolerant sub-population that resists RSL3 and reverts off-drug) are robust to the exact rates.

The SDT/PDT resistance of persisters currently comes only through the (weak) MUFA axis; whether to add an explicit reduced-lipid-peroxide-vulnerability term to the exo-ROS path is deferred (#262 out-of-scope note).

## ALOX isoform peroxidation + MCFA sensitization (`alox::AloxConfig` → `Params.alox_propagation_boost`, `Params.mcfa_pufa_boost`): known uncalibrated

The core engine drives lipid peroxidation through one generic `lp_propagation`, implicitly assuming an average enzymatic-oxidation capacity. `AloxConfig` (#446) lets a consumer perturb that around the implicit baseline using an arachidonate-lipoxygenase isoform activity/expression mix and an MCFA level, collapsed to two **off-by-default additive boosts** written onto `Params`. Both default to `0.0` ⇒ ×1.0 ⇒ byte-identical (the production matrix never sets them; the FFI defaults them so the C ABI is unchanged).

| Parameter | Default | `literature()` | Source | Grounded? | Sensitivity |
|-----------|---------|----------------|--------|-----------|-------------|
| `alox_propagation_boost` | 0.0 | +0.355 (ALOX15-high) | ALOX15 is the canonical ferroptosis-driving isoform; ALOX12/ALOX5 contribute at isoform-specific rates (lipoxygenase-driven ferroptosis, PNAS 2016 PMID 27506793; Yang & Stockwell) | Assumed (direction only) | High — multiplies the autocatalytic propagation rate `1+boost` (clamped ≥0; `-1` = ALOX-null, no enzymatic propagation) |
| `mcfa_pufa_boost` | 0.0 | +0.25 (moderate MCFA) | MCFA → ACSL4/CD36 upregulation → more oxidizable PUFA incorporation (Sci Rep 2024 s41598-024-55050-4; MCFA ferroptosis sensitization PMC11901882) | Assumed (direction only) | Moderate — added to the oxidizable-PUFA augmentation alongside the ether-lipid pool |

**Calibration target.** Isoform Kcat/Km vary ~10-fold and are not fit here; the `literature()` activity weights/expression fractions and the MCFA saturation strength are placeholders, so only the DIRECTIONS are claimed (ALOX-high ⇒ more ferroptosis; MCFA ⇒ more ferroptosis; ALOX-null ⇒ resistant). Fitting would need ALOX-isoform-stratified expression (TCGA/RNA-seq) paired with isoform-knockdown and MCFA-exposure ferroptosis dose-response; absolute MCFA kinetics are deferred to the experimental E-series. Per-cell stochastic ALOX heterogeneity (per-phenotype sampling like iron/gsh) is a deferred refinement — this models a per-condition ALOX phenotype.

## ACSL4-status biomarker stratification (`acsl4` -> `Params.acsl4_status_boost`): known uncalibrated / validation data-gated

ACSL4 ligates the polyunsaturated fatty acids ferroptosis requires into membrane phospholipids, so a tumor's ACSL4 expression status sets its oxidizable-PUFA baseline and is the single most discriminating pro-ferroptotic lipid-metabolism gene (Doll et al., Nat Chem Biol 2017, PMID 27842070, ACSL4 dictates ferroptosis sensitivity; Yang et al., PNAS 2016, PMID 27506793). `acsl4::pufa_boost_from_status(status)` maps a relative status (`1.0` = wild-type) to an **off-by-default** additive boost on `Params.acsl4_status_boost`, folded into the oxidizable-PUFA augmentation (`biochem::ether_augmented_pufa`). Default `0.0` => x1.0 => byte-identical; FFI defaults it so the C ABI is unchanged.

| Status | `acsl4_status_boost` | Meaning | Source | Grounded? |
|--------|----------------------|---------|--------|-----------|
| ACSL4-negative (0.0) | -1.0 (null floor) | PUFA substrate collapses => ferroptosis-REFRACTORY (e.g. some HCC/AML) | Doll 2017 PMID 27842070 | Assumed (direction only) |
| ACSL4-low (0.5) | -0.5 | partially resistant | Doll 2017 | Assumed |
| ACSL4-normal (1.0) | 0.0 | model baseline (byte-identical) | - | N/A |
| ACSL4-high (1.5) | +0.5 | more oxidizable PUFA => sensitive (e.g. lung, ER+ breast, cervical) | Doll 2017 | Assumed (placeholder magnitude) |

**Calibration target (partially anchored, #462).** The status->boost MAGNITUDE mapping is still an uncalibrated linear placeholder, but leg (1) of #444 is now done: `scripts/fetch_acsl4_prevalence.py` pulls per-cancer-type ACSL4/GPX4/SLC7A11 mRNA from cBioPortal TCGA PanCancer Atlas (32 studies, login-free; committed `analysis/calibration/acsl4_prevalence_tcga.{csv,json}` + `acsl4-prevalence-calibration.md`). The within-cohort low-ACSL4 prevalence (z<-1: 10.8-18.8%, median 14.4%; z<-2: median 3.0%) is the real per-cancer-type prior, and `acsl4::status_from_zscore(z)=max(0,1+z/2)` bridges a real z-score onto the status scalar (its integer-z points reproduce the shipped constants exactly). HONEST NEGATIVE: bulk TCGA mRNA does NOT show HCC/AML as ACSL4-low (lihc ranks highest on raw RSEM; cross-study RSEM is batch-confounded), so the Doll-2017 refractory phenotype is a protein/subtype property, not a bulk-mRNA one. STILL data-gated: leg (2), the cell-line ACSL4-status-vs-ferroptosis-inducer dose-response meta-analysis (CTRPv2 / GDSC / DepMap + primary literature) that would fit the status->IC50-shift MAGNITUDE, plus protein-level (RPPA/IHC) refractory-subtype prevalence. The model-side stratification (the kill-ordering A/B + the falsifiable prediction that ACSL4-high cancers respond to ferroptosis inducers while ACSL4-negative do not) is the deliverable as far as public data allows.

## ESCRT-III membrane-repair brake (`repair` -> `Params.escrt_repair_rate`, `Params.escrt_repair_budget`): known uncalibrated

A genuinely new category (#465): every other parameter modulates the redox/lipid SUBSTRATE, but ESCRT-III membrane repair acts on the death-EXECUTION step. A cell whose lipid peroxide has crossed `death_threshold` can be resealed (Ca2+ -> CHMP5/CHMP6) for a finite per-cell budget, so more repair capacity means slower/blocked death and more resistance (Dai et al., BBRC 2020, PMID 31761326; CHMP5/CHMP6 knockdown sensitizes). Both fields are **off-by-default**; the death-check RNG roll is drawn only when `escrt_repair_rate > 0` and budget remains, so the default path is byte-identical. Not in the C ABI (FFI defaults both to 0.0).

| Parameter | Default | Source | Grounded? | Sensitivity |
|-----------|---------|--------|-----------|-------------|
| `escrt_repair_rate` | 0.0 | Per-step rescue probability when over threshold; Dai 2020 PMID 31761326 (direction only) | Assumed (direction only) | High when on (gates how often execution is delayed) |
| `escrt_repair_budget` | 0.0 | Finite per-cell repair capacity (number of rescue events); tracked on `CellState.escrt_budget_used` | Assumed (direction only) | High when on (gates the death-delay length / whether the cell survives the run) |

**Calibration target.** Only the DIRECTION (more ESCRT repair => slower execution => more resistance; CHMP5/CHMP6 loss => sensitization) is claimed; the rate and budget are placeholders. Fitting would need CHMP5/CHMP6-knockdown vs wild-type ferroptosis-sensitivity dose-response plus membrane-resealing kinetics (Ca2+ imaging / annexin-V time-courses) to set the per-step rescue probability and the finite repair capacity.

## Photosensitizer pharmacokinetics: plasma vs. cellular

`Photosensitizer::Porfimer.t_half_h` represents *plasma* terminal half-life. Cellular concentration is assumed to track plasma proportionally — a reasonable approximation for porfimer (slow-distributing, weeks-scale t½, ~100% serum-protein bound, Vd ≈ plasma volume per Bellnier 2006) but explicitly wrong for 5-ALA/PpIX, which accumulates intracellularly via ferrochelatase deficiency rather than decaying. ALA kinetics will require a different variant.

### Distribution-phase model (closed via #203)

`Porfimer.t_distribution_h` holds drug at peak for the first N hours after administration, then begins single-exponential decay. Default `0.0` reproduces the pre-#203 "light at peak" model bit-exactly. Bellnier 2006 reports porfimer redistribution over ~24-48 h; setting `t_distribution_h` to the midpoint (~36 h) lets users pass clinical DLI from injection directly.

### Inter-drug ROS-yield normalization (closed via #203)

`Porfimer.phi_so2_relative` scales `concentration_at` to give the per-photon ROS yield via `Photosensitizer::yield_at`. The calibration anchor is **porfimer at peak = 1.0**; absolute porfimer phi_so2 in solution is consensus-cited as ≈ 0.65 across PDT literature reviews (e.g., Wilson & Patterson 2008, Phys Med Biol 53(9):R61–109), with primary measurements varying by formulation and solvent (Spikes & Bommer 1991 and earlier work; the value is a community-anchored constant rather than a single primary citation). Other drug variants would set their `phi_so2_relative` to `absolute_phi_so2 / 0.65` so `Params::pdt_ros = 5.0` (calibrated to porfimer) carries through correctly. `physics::pdt_intensity_at_depth` calls `yield_at` rather than `concentration_at` so the new fields compose into the existing Beer-Lambert path automatically.

Caveat: tissue phi_so2 values can be lower than solution values due to aggregation and microenvironment effects (Wilson & Patterson 2008 §5). The relative-to-porfimer convention encodes the calibration baseline in the type system but does not eliminate the underlying empirical uncertainty in absolute values.

## External ODE/QSP models for cross-validation (#344)

The biochem ODEs were cross-validated *qualitatively* against independent
published ferroptosis dynamical-systems models (same kind of system, same kind
of behavior). The shared, structurally-required behavior is a bistable
GSH/GPX4-threshold recover-or-collapse switch driven by Fenton positive feedback;
ferroptosis-core reproduces it (a bimodal single-cell lipid-peroxide distribution
at the tipping dose plus a sharp population dose-response threshold). See
`analysis/ode-cross-validation.md` and `cross_validate_odes.py`. All PMIDs
verified via NCBI esummary.

| Model | PMID | Role in cross-validation |
|---|---|---|
| Co et al., *Nature* 2024, "Emergence of large-scale cell death through ferroptotic trigger waves" | 38987590 | Canonical monostable -> bistable bifurcation as antioxidant defense falls; the unstable-threshold separatrix our bimodal LP distribution reproduces. |
| Seidel et al., *Front Cell Dev Biol* 2026, "A feedback loop between cell proliferation and ROS regulates ferroptosis sensitivity" | 41960191 | Minimal 2-ODE with two stable states + threshold death; the most directly comparable continuous model (peer-reviewed version of bioRxiv 2025.09.15.676259, no PMID). |
| Konstorum et al., *J Theor Biol* 2020, "Systems biology of ferroptosis: A modeling approach" | 32114023 | Discrete logical model; structural anchor for GPX4-as-critical-brake on LOOH (qualitative, not a trajectory comparator). |
| Pannala et al., *Free Radic Res* 2014, "A mechanistic mathematical model for the catalytic action of glutathione peroxidase" | 24456207 | Enzyme-kinetic GPX reference for the functional form of the GPX4/GSH repair term (saturating in GSH). |

<!-- SOBOL-IDENTIFIABILITY-START (generated by scripts/sobol_sensitivity.py, #331) -->
## Practical identifiability from the kill observable (#331)

Sobol total-effect screening at the Persister+RSL3 operating point (see `analysis/sobol-sensitivity-report.md`) over the PRCC's biochemical rate constants (10 parameters, the PRCC's 11 minus the SDT-only `sdt_ros`) ranks how much the single-cell kill rate constrains each. Parameters the kill rate is insensitive to (total-effect ST < 0.05) are NOT constrainable from kill-rate calibration and are marked here so they are not read as data-fitted (scope: kill-rate observable only; a parameter flagged non-identifiable here may be constrainable from an LP-timecourse or GSH-depletion observable):

| Parameter | ST | constrainable from kill rate? |
|---|--:|:--:|
| `lp_propagation` | 0.504 | yes |
| `gpx4_rate` | 0.285 | yes |
| `lp_rate` | 0.177 | yes |
| `rsl3_gpx4_inhib` | 0.068 | yes |
| `gsh_scav_efficiency` | 0.048 | NO (kill-rate-insensitive) |
| `nrf2_gsh_rate` | 0.035 | NO (kill-rate-insensitive) |
| `fsp1_rate` | 0.018 | NO (kill-rate-insensitive) |
| `fenton_rate` | 0.011 | NO (kill-rate-insensitive) |
| `death_threshold` | 0.003 | NO (kill-rate-insensitive) |
| `gpx4_degradation_by_ros` | 0.000 | NO (kill-rate-insensitive) |
<!-- SOBOL-IDENTIFIABILITY-END -->
