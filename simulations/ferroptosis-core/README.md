# ferroptosis-core

Embeddable ferroptosis biochemistry engine for cancer simulation.

This library provides a mechanistic model of the ferroptosis cell death pathway — from ROS generation through GSH depletion, GPX4 inactivation, lipid peroxidation, and cell death. Some parameters are grounded in published measurements; others are estimated from literature ranges or assumed as mechanistic placeholders (see parameter provenance below). It is designed to be embedded in multi-scale cancer simulators (PhysiCell, CompuCell3D, custom frameworks) or used standalone for parameter exploration.

## Quick start

```rust
use rand::prelude::*;
use ferroptosis_core::biochem::sim_cell;
use ferroptosis_core::cell::{gen_cell, Phenotype, Treatment};
use ferroptosis_core::params::Params;

let params = Params::default();
let mut rng = StdRng::seed_from_u64(42);
let cell = gen_cell(Phenotype::Persister, &mut rng);
let mut sim_rng = StdRng::seed_from_u64(43);

let (dead, lp, gsh, gpx4) = sim_cell(&cell, Treatment::SDT, &params, &mut sim_rng);
println!("Dead: {dead}, LP: {lp:.2}, GSH: {gsh:.2}, GPX4: {gpx4:.2}");
```

Run the included example: `cargo run -p ferroptosis-core --example basic_usage`

## Modules

| Module | Purpose |
|--------|---------|
| `cell` | Cell types, phenotypes (Glycolytic, OXPHOS, Persister, PersisterNrf2, Stromal), treatments, stochastic cell generation |
| `photosensitizer_pk` | Photosensitizer plasma PK and drug-light-interval scaling for PDT |
| `params` | All rate constants: `Params` (biochemistry), `SpatialParams` (physics), `ImmuneParams` (immune cascade), `RecoveryRates` |
| `biochem` | Core simulation engine: `sim_cell` (full 180-step loop), `sim_cell_step` (single timestep for spatial interleaving) |
| `stats` | Wilson confidence intervals, parallel Monte Carlo execution via rayon |
| `physics` | Depth-dependent energy deposition: Beer-Lambert (PDT), acoustic attenuation (SDT), uniform (RSL3). 2D row-based (`local_ros_multiplier`) and 3D radial-depth (`local_ros_multiplier_3d`) dispatchers share the same per-treatment depth functions (#186). |
| `grid` | 2D `TumorGrid` (8-Moore, circular) and 3D `TumorGrid3D` (26-Moore, spherical) with heterogeneous architecture, neighbor iteration, iron diffusion. `TumorGrid3D::radial_depth_um` provides per-cell signed depth from the spheroid surface for energy physics (#185, #186), hoisted via `RadialDepthGeom` for per-cell field sweeps (#289). The consuming binary is `sim-tme-3d` (#195; the standalone sim-spatial-3d #194 was closed as superseded). |
| `oxygen` | 3D radial O₂ gradients for spheroid tumors: `radial_o2_field` (per-cell `exp(-d/λ)` factor) and `radial_o2_zone_kill_rates` (normoxic / transition / hypoxic zone census). First-order Krogh approximation; pure functions for composable cycling (#187). Plus three off-by-default O2-coupling helpers (each identity at its `0` default ⇒ byte-identical) that let a consumer make the hypoxia leg of §7.1 runnable from both sides: `o2_dependent_exo_factor` (#336, the Type II / O2-dependent fraction of SDT/PDT exo-ROS, `(1−dep)+dep·o2`), `hypoxia_iron_factor` (#340, HIF/TfR1 raising labile iron in hypoxia, `1+sens·(1−o2)`), and `fenton_o2_factor` (#383, the O2-derived-H₂O₂ Fenton substrate, `(1−dep)+dep·o2`) — the last two compose multiplicatively on `cell.iron` so hypoxia raises the iron but starves the Fenton H₂O₂, correcting the #365 deep-core artifact (Haber-Weiss, Kehrer 2000 PMID 10963860; ferroptotic LPO is an O₂-consuming lipid-radical chain, reviewed in Stockwell 2017 PMID 28985560). |
| `ph` | 3D radial pH gradient for spheroid tumors: `radial_ph_field` (per-cell `pH(d) = ph_edge - delta·(1 - exp(-d/λ))`) plus pure-scalar helpers `iron_multiplier_from_ph` (ferritin destabilization) and `ion_trap_factor_from_ph` (linearized Henderson-Hasselbalch, weak-base drug bioavailability). Same form as sim-tme's 2D pH; pure functions for #195 sim-tme-3d (#190). |
| `stromal` | CAF-shielded boundary detection: `stromal_adjacency_mask_3d` / `stromal_adjacent_kill_rate_3d` (`TumorGrid3D`, 26-Moore neighbors) and `stromal_adjacency_mask_2d` / `stromal_adjacent_kill_rate_2d` (`TumorGrid`, 8-Moore). Surface-to-volume scaling means 3D shielding affects ~1.5× more cells than 2D at matched R (#189, lifted from sim-tme's binary-local 2D copy in #224). |
| `immune` | ICD/DAMP immune cascade (dimensionless, single-event): ferroptotic death quality drives dendritic cell activation and T cell priming |
| `immune_spatial` | 3D spatial immune coupling: `diffuse_damp_3d_step` (26-Moore DAMP diffusion + exponential clearance, scratch-buffer API) plus pure-scalar helpers `dc_activation` (Michaelis-Menten), `immune_kill_probability` (with sim-tme's 0.99 cap), and `exhaustion_factor` (T-cell exhaustion `1/(1+rate·cumulative_kills)`, identity at rate 0; #243). Treg/MDSC immunosuppression (#264 Phase 2): `suppressor_kill_multiplier` (`1/(1+strength·field)`, identity at strength 0 — mirrors exhaustion), `suppressor_source_mask_3d` (perivascular niches from vessel positions, else heuristic patches via an INDEPENDENT RNG), and `SuppressorConfig` (`disabled()`/`enabled()`); the consumer replenishes a second field at the source cells and diffuses it via `diffuse_damp_3d_step`. Multi-checkpoint immune brake (#264 Phase 3): `CheckpointPanel` (PD-1/CTLA-4/LAG-3/TIM-3, each a `Checkpoint { brake, drug_efficacy }`) with `combined_brake()` = `1 − Π(1 − brakeᵢ·(1−drug_efficacyᵢ))`, which reduces exactly to the single PD-1 `effective_brake` when only PD-1 is active (byte-identical) and models anti-PD-1 + anti-CTLA-4 combinations. DC subset mix (#264 Phase 4): `DcSubsetConfig` collapses the cDC1/cDC2 composition to one uniform anti-tumor `priming_efficiency()` scalar the consumer multiplies into the immune kill probability; a cDC1-poor tumor primes CD8 killing less efficiently (cDC1/Batf3 cross-presenting DCs are the rare critical anti-tumor APCs, Broz 2014 PMID 25446897), `balanced()` ⇒ efficiency 1.0 (identity/byte-identical), `literature()` is a cDC1-poor placeholder. Immunosuppressive ferroptosis (#337): `ferroptotic_immunosuppression` (`1/(1+strength·local_damp)`, identity at strength 0) scales immune kill DOWN as ferroptotic-death density rises (extracellular GPX4 / oxidized-lipid DC suppression), keyed on the SAME DAMP that drives `dc_activation` so the net immune effect can flip sign at high death density. Diffusing SASP field (#376): `sasp_field_kill_multiplier(field, strength)` is the PARACRINE extension of the cell-autonomous `senescence::sasp_immune_multiplier` (#341) — the consumer seeds a SASP field at the senescent cells (`SenescenceConfig::sasp_field_strength`), diffuses it via `diffuse_damp_3d_step`, and applies this SIGNED multiplier to EVERY exposed cell including non-senescent neighbors (`strength > 0` immunosuppressive `1/(1+strength·field)` lowers neighbor kill, Di Mitri 2014 PMID 25156255; `strength < 0` surveillance `1+|strength|·field` raises it, Kang 2011 PMID 22080947; identity at field 0 or strength 0). **Stability requirement**: `diffusion_fraction × 26 < 1` (use ≤ 0.038; sim-tme's 2D default 0.08 is UNSAFE — `assert!`-enforced). Composes with `immune` for downstream sim-tme-3d (#195). (#188) |
| `io` | JSON and CSV output helpers |
| `drug_transport` | Krogh cylinder drug penetration model. Includes an off-by-default ECM-tortuosity factor (#315): `TissueParams::ecm_tortuosity` scales the penetration length `λ_eff = λ/√τ` so a denser desmoplastic matrix shortens drug penetration; `τ=1` identity for all shipped tissues (byte-identical), an uncalibrated placeholder (Netti 2000, Provenzano 2012) |
| `tumor_pk` | Two-compartment vascular/interstitial pharmacokinetics |
| `dose_schedule` | Time-varying drug-administration schedules (Constant / Bolus / MultiDose / Infusion / FromPk); `factor_at(step)` per-step availability, identity-default for byte-identical steady-state (#239) |
| `persister` | Drug-tolerant persister cells (#241): pure helpers (`step` competing-rate per-step update / `acquire` / `revert` / `gpx4_inactivation_multiplier` / `mufa_boost_increment`) + `PersisterConfig` (identity-default ⇒ no-op). Cells acquire epigenetic ferroptosis tolerance under drug exposure and revert after clearance; the consumer applies `step` (#262 — acquisition + reversion act simultaneously, so sustained sub-saturating drug reaches a sub-cap equilibrium rather than ratcheting to the cap). Reversible-to-irreversible epigenetic locking via `step_with_locking` + `PersisterState` (#342). TWO entry routes: drug-driven acquisition, and a non-drug stress-niche entry `stress_entry(state, stress, cfg)` (#377) that raises the reversible pool from a hypoxic/nutrient-poor drug-sanctuary niche, decoupled from drug (stress drives entry, drug drives durability), off-by-default (`stress_entry_rate=0`) byte-identical. Consumer owns `CellState::persister_fraction` |
| `spheroid` | 3D spheroid radial cell biology (#197): `apply_radial_cells_3d` (re-assigns tumor cells radially — glycolytic rim / OXPHOS mid / persister core — via an INDEPENDENT RNG, with core-low GSH + core-high iron gradients, and a per-cell **MUFA cap** rim-high/core-low so position-dependent MUFA is durable, #270) + `radial_phenotype` / `radial_mufa_protection` / `radial_fraction_3d` + `SpheroidConfig`. Pairs with `Params::spheroid()` (partial-MUFA context). Opt-in ⇒ default random grid byte-identical |
| `vasculature` | Explicit 3D vessel network (#191): `place_vessels_3d` (random internal vessel seeds via an INDEPENDENT RNG; count from inter-vessel spacing) + `vessel_supply_field` (per-cell `exp(-dist_to_nearest_vessel/λ)`, a drop-in alternative to `oxygen::radial_o2_field` supplying both O2 and drug) + `hypoxic_fraction`. `VasculatureConfig::well_vascularized()` / `poorly_vascularized()`. Replaces the edge-distance proxy with patchy, non-radial oxygenation. `place_vessels_in_slab_3d` (#272) scatters vessels uniform-in-box (not the central sphere) so the patient-scale slab (#240) can carry internal vessels alongside the planar depth gradient. `vessel_supply_field`'s nearest-vessel lookup uses a uniform-grid spatial index (#268) — EXACT (byte-identical to the former brute force, verified bit-for-bit), roughly O(cells) instead of O(cells×vessels), so vasculature scales to patient-size grids (100³/1M cells/2370 vessels ≈ 105 ms). `place_vessels_fractal_3d` + `VesselTopology::{Random,Fractal}` (#268) replace uniform-random points with a fractal-branching tree (trunks enter from the periphery and bifurcate inward with high, tumor-like variability — Baish & Jain 2000) — a hierarchical-but-chaotic network with avascular gaps/dead ends, capped at the same *point-count* target as random (matches raw point count, NOT effective coverage — the clustered, 1-cell-spaced branch points cover far less unique volume than the same number of scattered points, which is why it leaves more hypoxic tissue; read the effect qualitatively). Slab geometry ignores `topology`. Off-by-default (topology=Random) byte-identical |
| `clonal` | Clonal heterogeneity (#242): `assign_subclones_3d` (Voronoi subclone map via an INDEPENDENT RNG, so `TumorGrid3D::generate`'s stream is untouched) + `ClonalConfig` / `SubclonePerturbation` (per-subclone `iron_mul` / `gpx4_mul` / `lipid_unsat_mul` the consumer applies as RNG-neutral setup mutations; `lipid_unsat_mul` is the MUFA-enrichment axis, scaling the static `Cell` PUFA field so it persists across steps; `gpx4_mul` is the **durable** antioxidant axis (#266) — the consumer scales both the initial `state.gpx4` and the static `cell.nrf2` setpoint, so a low-antioxidant subclone stays differentiated for the whole run rather than relaxing back). `single_identity()` (K=1) ⇒ byte-identical; `literature_4()` spans the mesenchymal⇄epithelial vulnerability axis. `repopulate_dead_sites_3d` adds spatial clonal **expansion** (#266 item 3): dead tumor sites are repopulated (two-phase, deterministic) from living Moore-neighbors, so resistant subclones grow their territory; gated on `ClonalConfig::with_repopulation(rate)` (`0` ⇒ static map, byte-identical) |
| `slab` | Patient-scale slab geometry (#240): `TumorGrid3D::generate_slab` (an all-tumor block, no sphere carve) + `slab_supply_field` (per-cell planar `exp(-depth/λ)` where the +z face is vessel-proximal at `depth_offset` and supply decays toward −z; the 1-D analog of `oxygen::radial_o2_field`, supplying both O2 and drug) + `scale_interpretation` + `SlabConfig` (`patient_deep()` = a deep, ~drug-deprived slab in a 10 mm virtual tumor; `surface()` = the shallow control). `KROGH_LAMBDA_UM` (~150 µm) is the default penetration length. Models the depth-dependent penetration collapse the in-vitro spheroid scale misses (magnitude is an uncalibrated first-order Krogh approximation). `apply_depth_graded_cells_3d` + `SlabPhenotypeConfig` + `depth_phenotype` + the `layer_supply` helper (#272) add an opt-in depth-graded phenotype: the flat bulk mix becomes a layered rim→core gradient (proliferating/glycolytic at the +z vessel face, persister-like in the chronically supply-starved deep (−z) layers), thresholded on the planar supply `exp(-depth/λ)` (NOT geometric volume fractions like the spheroid, since the slab models an *absolute* depth, so a `patient_deep()` slab is uniformly persister-like). Independent per-cell RNG; supply cut-points are uncalibrated placeholders. Opt-in ⇒ default spheroid byte-identical |
| `contact` | Cell-cell contact-mediated ferroptosis resistance (#270): `apply_contact_resistance_3d` + `contact_fraction_3d` + `ContactConfig`. Dense, highly-contacting tumor cells resist ferroptosis (E-cadherin junctions → Merlin/NF2 → YAP inhibition → ACSL4/TFRC down; Wu 2019, PMID 31341276). `contact = (tumor 26-Moore neighbours)/26`; the layer scales the durable `cell.lipid_unsat` (PUFA) and `cell.iron` down by `1 − strength·contact`, so interior cells resist while sparse/surface cells stay sensitive (stronger in 3D, up to 26 neighbours vs 8 in 2D). Effectively a per-cell ferroptosis-threshold modulation. Geometric (no RNG); `literature()` strengths are UNCALIBRATED placeholders; off-by-default identity ⇒ byte-identical |
| `phenotype_mufa` | Phenotype-specific SCD1/MUFA dynamics (#363 rate + #390 cap): `PhenotypeMufaConfig` (per-phenotype RATE multipliers on `scd_mufa_rate` + `*_cap` CAP multipliers on the effective `mufa_cap`; `identity()` default) + `apply_phenotype_mufa_3d` set each tumor cell's per-cell `Cell::mufa_rate` and scale its `Cell::mufa_cap` (`None`/identity ⇒ global ⇒ byte-identical), so the acute-vs-established MUFA build-up (#339) gets a phenotype-specific time constant AND steady state instead of one shared rate/cap. The cap is applied multiplicatively (cap-mul `1.0` leaves the cap untouched) AFTER the spheroid, so a phenotype cap composes with the spheroid's radial cap. DIRECTION is genuinely uncertain (MUFA enrichment resists ferroptosis, Magtanong 2019 PMID 30686757, and SCD1 protects, Tesfay 2019 PMID 31270077, but persisters are also GPX4-dependent/vulnerable, Hangauer 2017 PMID 29088702), so `literature()` is an UNCALIBRATED placeholder. Not in the C ABI. Off-by-default identity ⇒ byte-identical |
| `nutrient` | Radial nutrient gradient (#270 item 3b): `apply_nutrient_stress_3d` + `nutrient_availability` + `NutrientConfig`. Beyond O2 (#187) and pH (#190), glucose/glutamine are abundant at the rim and consumed toward the core; glucose metabolism feeds the NADPH (pentose-phosphate pathway) that regenerates GSH for the GPX4 defense, so a nutrient-starved core has less antioxidant capacity and is MORE ferroptosis-sensitive (Dixon 2012, PMID 22632970; glucose metabolic reprogramming regulates ferroptosis, PMID 42190602). availability = `exp(-radial_depth/λ)` (the O2 field's form); deprivation `1 − availability` scales the durable antioxidant setpoint `cell.nrf2` down. Direction caveat: ONE documented direction (energy stress also activates AMPK which INHIBITS ferroptosis, and glutaminolysis is REQUIRED for some, PMID 30581146); net is context-dependent. Geometric (reuses `RadialDepthGeom`, no RNG); `literature()` strength is an UNCALIBRATED placeholder; off-by-default identity ⇒ byte-identical |

## Key API

**Single-cell simulation (full loop):**
```rust
sim_cell(cell, treatment, params, rng) -> (dead, lp, gsh, gpx4)
```

**Single timestep (for embedding in spatial/multi-scale models):**
```rust
sim_cell_step(state, cell, params, step, extra_iron, rng) -> dead
```

**Cell generation:**
```rust
gen_cell(phenotype, rng) -> Cell
```

**Photosensitizer PK (PDT light-dose scaling):**
```rust
let ps: Photosensitizer = "porfimer=504,36,0.65".parse()?;  // FromStr
ps.concentration_at(t_h);  // drug present at time t_h post-administration
ps.yield_at(t_h);          // ROS yield = concentration × phi_so2_relative
```
Variants: `Uniform(c)` (constant fraction; default 1.0 = no PK model) and
`Porfimer { t_half_h, t_distribution_h, phi_so2_relative }` (single-
exponential plasma decay with optional saturating distribution-phase
hold and relative singlet-O₂ yield). All defaults preserve identity-
preserving physics. `physics::pdt_intensity_at_depth` calls `yield_at`
to compose drug presence + yield with depth-attenuated light.

**3D spheroid energy physics (#186):**
```rust
let g = TumorGrid3D::generate(40, 40, 40, 20.0, 42);
let depth_um = g.radial_depth_um(r, c, l);   // signed: + inside, − outside
let m = local_ros_multiplier_3d(depth_um, Treatment::PDT, &spatial_params);
```
Negative depths (cells outside the spheroid) are clipped to the surface
value. The 3D dispatcher reaches the same per-treatment depth functions
the 2D path uses, so the matched-depth invariant `local_ros_multiplier(row,
cell_size, ...) == local_ros_multiplier_3d(row × cell_size, ...)` holds
bit-exact across all `Treatment` variants — the physical *geometries*
differ (planar slab vs. spheroid + nearest-surface 1-D approximation),
but the dispatcher math does not.

**3D spatial immune coupling (#188):**
```rust
use ferroptosis_core::immune_spatial::{diffuse_damp_3d_step, dc_activation, immune_kill_probability};

let g = TumorGrid3D::generate(40, 40, 40, 20.0, 42);
let n = g.cells.len();
let mut damp_field = vec![0.0_f64; n];
let mut scratch = vec![0.0_f64; n];  // allocate ONCE, reuse per step

// Inject DAMP from a death event, then diffuse one step.
damp_field[g.flat_index(20, 20, 20)] = 10.0;
diffuse_damp_3d_step(&mut damp_field, &mut scratch, &g, 0.025, 0.03);  // 3D-safe fraction

// Per-cell immune kill probability.
let activation = dc_activation(damp_field[g.flat_index(20, 20, 21)], 50.0);
let kill_prob = immune_kill_probability(activation, 0.02, 0.21);
```
**⚠️ Stability**: `diffusion_fraction × 26 < 1` is `assert!`-enforced (sim-tme's 2D default 0.08 is UNSAFE in 3D — 0.08 × 26 = 2.08 silently destroys mass). Use ≤ 0.038; suggested 0.025 to match 2D's per-step total. Composes with [`immune`] (dimensionless single-event cascade) for downstream sim-tme-3d (#195).

**3D spheroid stromal shielding (#189):**
```rust
use ferroptosis_core::stromal::{stromal_adjacency_mask_3d, stromal_adjacent_kill_rate_3d};

let g = TumorGrid3D::generate(40, 40, 40, 20.0, 42);
let mask = stromal_adjacency_mask_3d(&g);   // Vec<bool>, true = boundary tumor cell
let rate = stromal_adjacent_kill_rate_3d(&g, &mask);  // kill rate in shielded shell
```
The 2D analog (`stromal_adjacency_mask_2d` / `stromal_adjacent_kill_rate_2d`, used by sim-tme) lives in the same module.
3D analog of sim-tme's 2D 8-Moore boundary detection, using 26-Moore
neighbors. The shielded shell is one cell deep; consumers apply CAF
GSH/MUFA boosts to flagged cells (mutation stays consumer-side per the
pure-functions pattern). Cross-geometry test confirms 3D boundary
fraction > 2D at matched R (surface-to-volume scaling: 3/R vs 2/R).

**3D spheroid pH gradient (#190):**
```rust
use ferroptosis_core::ph::{radial_ph_field, iron_multiplier_from_ph, ion_trap_factor_from_ph};

let g = TumorGrid3D::generate(40, 40, 40, 20.0, 42);
let (ph_edge, ph_core, lambda) = (7.4, 6.5, 120.0);
let ph_field = radial_ph_field(&g, ph_edge, ph_core, lambda);

for (i, &local_ph) in ph_field.iter().enumerate() {
    let iron_mult = iron_multiplier_from_ph(local_ph, ph_edge, 1.5);   // 2.35× at core
    let drug_factor = ion_trap_factor_from_ph(local_ph, ph_edge, 0.4); // 0.64 at core
    // consumer applies: cell.iron *= iron_mult; effective_drug = base × drug_factor
}
```
Stromal cells return `ph_edge` (well-perfused). Pure functions follow the
oxygen-module pattern — consumer chooses mutation strategy. Same code shape
as sim-tme's 2D `apply_ph_gradient`; cross-geometry test
(`matched_lambda_2d_vs_3d_acidic_fraction`) shows pure-geometry 3D
acidic-volume fraction is *smaller* than 2D at matched λ (same
cubic-vs-quadratic effect as O₂).

**3D spheroid oxygen gradient (#187):**
```rust
use ferroptosis_core::oxygen::{radial_o2_field, radial_o2_zone_kill_rates};

let g = TumorGrid3D::generate(40, 40, 40, 20.0, 42);
let o2 = radial_o2_field(&g, 100.0);                 // Vec<f64>, length = g.cells.len()
let (norm, trans, hyp) = radial_o2_zone_kill_rates(&g, 100.0);
```
Stromal cells (outside spheroid) get O₂ = 1.0 (well-oxygenated bulk
tissue). Pure functions return values rather than mutating
`cell.basal_ros`, so the consumer composes cycling by re-calling
`radial_o2_field` per step with alternating λ (no `original_ros`
snapshot needed). First-order Krogh approximation; the exact
Krogh-Riley spheroidal solution involves `sinh` ratios — same
approximation level as `sim-tme`'s 2D `apply_o2_gradient`. Note:
pure-geometry 3D *hypoxic* fraction at matched R and λ is *smaller*
than 2D (the cubic-vs-quadratic scaling dominates the normoxic shell,
not the hypoxic core — the Vaupel 1989 observation that real 3D
spheroids are more hypoxic reflects vasculature biology, not geometry).

**Drug-tolerant persister cells (#241):**
```rust
use ferroptosis_core::params::PersisterConfig;
use ferroptosis_core::persister::{step, gpx4_inactivation_multiplier, mufa_boost_increment};

let cfg = PersisterConfig::enabled();        // ::default() is the identity no-op
let mut frac = 0.0;
// Competing-rate update (#262): acquisition + reversion both act each step, so
// sustained sub-saturating drug settles at a sub-cap equilibrium.
frac = step(frac, drug_intensity, &cfg);
let resist = gpx4_inactivation_multiplier(frac, &cfg); // ∈ [1-gpx4_resistance, 1]
let mufa_inc = mufa_boost_increment(frac, &cfg);       // per-step MUFA protection
```
Pure functions (consumer owns `CellState::persister_fraction` and applies the
effects, like `oxygen` / `ph` / `stromal`). `PersisterConfig::default()` is the
identity element (every helper a no-op), so an off run is byte-identical to no
model at all. Refs: Hangauer 2017 (persister ⇄ GPX4 dependence), Tsoi 2018
(MUFA lipid rewiring), Viswanathan 2017 (mesenchymal ⇄ ferroptosis / reversion).

**Parameter contexts:**
- `Params::default()` — 2D culture baseline
- `Params::invivo()` — 3D/in-vivo with SCD1-driven MUFA lipid remodeling (M_ss = 0.40)

## Parameters

All ~30 simulation parameters are documented with literature sources and sensitivity ratings in [`parameter_provenance.md`](https://github.com/ELares/cancer_research/blob/main/simulations/calibration/parameter_provenance.md) (in the parent repository's `simulations/calibration/` directory).

## License

MIT
