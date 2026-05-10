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
| `iron_diffusion_coeff` | 281.0 µm²/s | Jacques SL, Phys Med Biol 2013 | Grounded | Low — bystander iron diffusion |
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

- **Grounded** (value derived from specific published measurement): `gsh_max`, `pdt_mu_eff`, `sdt_alpha`, `cell_size_um`, `iron_diffusion_coeff`, `scd_mufa_rate`, `scd_mufa_max`
- **Estimated** (informed by literature ranges but not directly calibrated): most biochemistry rates
- **Assumed** (mechanistic placeholder with no direct data): `gsh_km`, `gpx4_degradation_by_ros`, `gpx4_nrf2_upregulation`, `death_threshold`, immune cascade parameters
- **Derived** (calculated from other parameters): `initial_mufa_protection`

## RSL3 pharmacokinetics: known uncalibrated

`tumor_pk::TumorPKParams` and the Krogh penetration model in `drug_transport` use RSL3-like parameters (e.g., plasma t½ ≈ 30 min, `k_uptake_bulk`, `km_uptake`) that are **order-of-magnitude estimates from chemical-probe literature, not clinical measurements.** RSL3 has no published clinical PK profile — it is widely cited as a research probe with poor pharmacokinetics, not a development candidate (e.g., review in Yang et al., Nature 2023, on ferroptosis therapeutics). Sensitivity of manuscript claims to these values is bounded by the protection-factor range reported in Chapter 8.2 (4.8×–27×) — qualitative conclusions about tumor-PK barriers are robust, but absolute kill rates should be read as approximate. A future issue should anchor PK parameters either to a clinically published ferroptosis inducer (e.g., IKE) or to a non-RSL3 reference compound.

## Photosensitizer pharmacokinetics: plasma vs. cellular

`Photosensitizer::Porfimer.t_half_h` represents *plasma* terminal half-life. Cellular concentration is assumed to track plasma proportionally — a reasonable approximation for porfimer (slow-distributing, weeks-scale t½, ~100% serum-protein bound, Vd ≈ plasma volume per Bellnier 2006) but explicitly wrong for 5-ALA/PpIX, which accumulates intracellularly via ferrochelatase deficiency rather than decaying. ALA kinetics will require a different variant.

### Distribution-phase model (closed via #203)

`Porfimer.t_distribution_h` holds drug at peak for the first N hours after administration, then begins single-exponential decay. Default `0.0` reproduces the pre-#203 "light at peak" model bit-exactly. Bellnier 2006 reports porfimer redistribution over ~24-48 h; setting `t_distribution_h` to the midpoint (~36 h) lets users pass clinical DLI from injection directly.

### Inter-drug ROS-yield normalization (closed via #203)

`Porfimer.phi_so2_relative` scales `concentration_at` to give the per-photon ROS yield via `Photosensitizer::yield_at`. The calibration anchor is **porfimer at peak = 1.0**; absolute porfimer phi_so2 ≈ 0.65 in solution per Wilson & Patterson 2008 (Phys Med Biol 53(9):R61-109). Other drug variants would set their `phi_so2_relative` to `absolute_phi_so2 / 0.65` so `Params::pdt_ros = 5.0` (calibrated to porfimer) carries through correctly. `physics::pdt_intensity_at_depth` calls `yield_at` rather than `concentration_at` so the new fields compose into the existing Beer-Lambert path automatically.

Caveat: tissue phi_so2 values can be lower than solution values due to aggregation and microenvironment effects (Wilson & Patterson 2008 §5). The relative-to-porfimer convention encodes the calibration baseline in the type system but does not eliminate the underlying empirical uncertainty in absolute values.
