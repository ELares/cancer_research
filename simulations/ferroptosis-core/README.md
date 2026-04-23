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
| `params` | All rate constants: `Params` (biochemistry), `SpatialParams` (physics), `ImmuneParams` (immune cascade), `RecoveryRates` |
| `biochem` | Core simulation engine: `sim_cell` (full 180-step loop), `sim_cell_step` (single timestep for spatial interleaving) |
| `stats` | Wilson confidence intervals, parallel Monte Carlo execution via rayon |
| `physics` | Depth-dependent energy deposition: Beer-Lambert (PDT), acoustic attenuation (SDT), uniform (RSL3) |
| `grid` | 2D tumor grid with heterogeneous architecture, neighbor iteration, iron diffusion |
| `immune` | ICD/DAMP immune cascade: ferroptotic death quality drives dendritic cell activation and T cell priming |
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

**Parameter contexts:**
- `Params::default()` — 2D culture baseline
- `Params::invivo()` — 3D/in-vivo with SCD1-driven MUFA lipid remodeling (M_ss = 0.40)

## Parameters

All ~30 simulation parameters are documented with literature sources and sensitivity ratings in [`parameter_provenance.md`](https://github.com/ELares/cancer_research/blob/main/simulations/calibration/parameter_provenance.md) (in the parent repository's `simulations/calibration/` directory).

## License

MIT
