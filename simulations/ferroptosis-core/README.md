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
| `grid` | 2D `TumorGrid` (8-Moore, circular) and 3D `TumorGrid3D` (26-Moore, spherical) with heterogeneous architecture, neighbor iteration, iron diffusion. `TumorGrid3D::radial_depth_um` provides per-cell signed depth from the spheroid surface for energy physics (#185, #186). 3D analytics (radial-depth curves, volumetric heatmaps) and the consuming binary land with #194. |
| `oxygen` | 3D radial O₂ gradients for spheroid tumors: `radial_o2_field` (per-cell `exp(-d/λ)` factor) and `radial_o2_zone_kill_rates` (normoxic / transition / hypoxic zone census). First-order Krogh approximation; pure functions for composable cycling (#187). |
| `ph` | 3D radial pH gradient for spheroid tumors: `radial_ph_field` (per-cell `pH(d) = ph_edge - delta·(1 - exp(-d/λ))`) plus pure-scalar helpers `iron_multiplier_from_ph` (ferritin destabilization) and `ion_trap_factor_from_ph` (linearized Henderson-Hasselbalch, weak-base drug bioavailability). Same form as sim-tme's 2D pH; pure functions for #195 sim-tme-3d (#190). |
| `stromal` | 3D CAF-shielded boundary detection for spheroid tumors: `stromal_adjacency_mask` (Vec<bool> flagging tumor cells with any stromal 26-Moore neighbor) and `stromal_adjacent_kill_rate` (dead-rate among masked cells). 3D analog of sim-tme's 2D 8-Moore detection; surface-to-volume scaling means 3D shielding affects ~1.5× more cells than 2D at matched R (#189). |
| `immune` | ICD/DAMP immune cascade (dimensionless, single-event): ferroptotic death quality drives dendritic cell activation and T cell priming |
| `immune_3d` | 3D spatial immune coupling: `diffuse_damp_3d_step` (26-Moore DAMP diffusion + exponential clearance, scratch-buffer API) plus pure-scalar helpers `dc_activation` (Michaelis-Menten) and `immune_kill_probability` (with sim-tme's 0.99 cap). **Stability requirement**: `diffusion_fraction × 26 < 1` (use ≤ 0.038; sim-tme's 2D default 0.08 is UNSAFE — `assert!`-enforced). Composes with `immune` for downstream sim-tme-3d (#195). (#188) |
| `io` | JSON and CSV output helpers |
| `drug_transport` | Krogh cylinder drug penetration model |
| `tumor_pk` | Two-compartment vascular/interstitial pharmacokinetics |

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
use ferroptosis_core::immune_3d::{diffuse_damp_3d_step, dc_activation, immune_kill_probability};

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
use ferroptosis_core::stromal::{stromal_adjacency_mask, stromal_adjacent_kill_rate};

let g = TumorGrid3D::generate(40, 40, 40, 20.0, 42);
let mask = stromal_adjacency_mask(&g);   // Vec<bool>, true = boundary tumor cell
let rate = stromal_adjacent_kill_rate(&g, &mask);  // kill rate in shielded shell
```
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

**Parameter contexts:**
- `Params::default()` — 2D culture baseline
- `Params::invivo()` — 3D/in-vivo with SCD1-driven MUFA lipid remodeling (M_ss = 0.40)

## Parameters

All ~30 simulation parameters are documented with literature sources and sensitivity ratings in [`parameter_provenance.md`](https://github.com/ELares/cancer_research/blob/main/simulations/calibration/parameter_provenance.md) (in the parent repository's `simulations/calibration/` directory).

## License

MIT
