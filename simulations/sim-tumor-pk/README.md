# sim-tumor-pk

Two-compartment PBPK (physiologically-based pharmacokinetic) model of plasma-to-tumor drug delivery, demonstrating that tumor-specific PK barriers (blood flow, vascular permeability, interstitial fluid pressure) create a large protection factor between 2D culture (constant drug) and in-vivo (time-varying, reduced concentration).

**Manuscript reference:** Chapter 8, Section 8.2

## What it models

1. **Plasma compartment** -- IV bolus with first-order elimination (t_half=30 min): C_plasma(t) = C0 * exp(-k_el * t)
2. **Tumor vascular compartment** -- blood flow delivers drug from plasma to tumor vasculature; vascular permeability and interstitial fluid pressure (IFP) govern extravasation
3. **Tumor interstitial compartment** -- drug concentration available to tumor cells; time-varying schedule fed into the ferroptosis biochemistry engine
4. **Five tumor types** with different PK parameters:
   - Breast (well-perfused)
   - Pancreatic (poorly-perfused, high IFP)
   - GBM/Glioblastoma (blood-brain barrier)
   - Melanoma (well-vascularized)
   - Sarcoma (poorly-vascularized)
5. **2D culture reference** -- constant drug concentration (C=1.0 for all 180 steps) to establish the baseline kill rate
6. **Spatial x temporal composition** -- C_i(t) from the temporal ODE is composed with Krogh spatial decay to produce C(r,t) kill rates at multiple distances from the vessel

## Quick start

```bash
cd simulations
cargo build --release -p sim-tumor-pk
cargo run --release -p sim-tumor-pk
```

Runtime: ~30-60 seconds (parallelized via rayon, 10K cells x multiple scenarios).

## Parameters

All parameters are hardcoded (no CLI arguments).

| Parameter | Value | Description |
|-----------|-------|-------------|
| N_CELLS | 10,000 | Cells per scenario |
| SEED | 42 | Random seed |
| N_STEPS | 180 | Biochemistry timesteps |
| Phenotype | Persister (FSP1-low) | Cell type for all conditions |
| Drug | RSL3-like IV bolus | t_half=30 min plasma elimination |

Tumor PK parameters are defined per tumor type in `ferroptosis_core::tumor_pk`:
- `breast_tumor()`, `pancreatic_tumor()`, `glioblastoma_tumor()`, `melanoma_tumor()`, `sarcoma_tumor()`

Spatial-temporal composition uses metabolism-only penetration length (lambda_met ~224 um for RSL3) to avoid double-counting cellular uptake.

Candidate radial distance bins for C(r,t): [0, 25, 50, 75, 100, 125] um.
Each tumor includes bins only up to its tissue half-distance (Breast/Melanoma
60 um, Pancreatic 125 um, GBM 75 um, Sarcoma 100 um), yielding 21 conditions.

## Output format

Output directory: `output/tumor-pk/`

### 1. `tumor_pk_summary.json` -- per-scenario results

JSON array with one entry per scenario (2D reference + 5 tumor types). Example
from the default 10,000-cell, seed-42 run, with floating-point values rounded:

```json
{
  "tumor_type": "Breast",
  "context": "tumor_pk",
  "n_cells": 10000,
  "n_dead": 251,
  "death_rate": 0.0251,
  "ci_low": 0.022212,
  "ci_high": 0.028353,
  "mean_lp": 0.882771,
  "mean_gsh": 1.078692,
  "mean_gpx4": 0.635620,
  "peak_c_interstitial": 0.408550,
  "auc_c_interstitial": 20.041354,
  "protection_factor": 16.310757
}
```

| Field | Description |
|-------|-------------|
| peak_c_interstitial | Maximum interstitial drug concentration (0.0-1.0 normalized) |
| auc_c_interstitial | Area under the interstitial concentration curve |
| protection_factor | 2D reference death rate / tumor death rate (higher = more PK barrier) |

### 2. `tumor_pk_timecourse.csv` -- concentration time series

Each row samples all three compartments at the same time, starting with the
zero vascular and interstitial initial conditions at minute 0. The 180 samples
cover minutes 0 through 179 and supply the concentration at the start of each
one-minute biochemistry step. Earlier versions advanced the tumor compartments
one minute before recording them, so previously generated timecourses and
PK-driven cell outcomes need to be regenerated with the corrected solver.

```csv
time_min,tumor_type,c_plasma,c_vascular,c_interstitial
0,Breast,1.000000,0.000000,0.000000
1,Breast,0.977160,0.764272,0.150001
...
```

The same sampling convention applies to the opt-in `sim-tme-3d --dose-sweep`
PK schedule. The default 3D matrix uses constant dosing and does not call this
solver. CSV plasma imports through `PlasmaModel::from_csv` reject nonfinite
times or concentrations rather than silently turning missing values into zero
exposure.

### 3. `tumor_pk_spatial_temporal.csv` -- C(r,t) kill rates

```csv
tumor_type,distance_um,peak_conc,death_rate,ci_low,ci_high,n_cells,n_dead
Breast,0,0.408550,0.025100,0.022212,0.028353,10000,251
Breast,25,0.365334,0.023100,0.020334,0.026233,10000,231
...
```

## Reproducing manuscript claims

**Chapter 8, Section 8.2 (2D-to-in-vivo gap):**
```bash
cargo run --release -p sim-tumor-pk
# stderr output includes "=== Protection Factor Summary ==="
# Expected:
#   2D culture ref: ~41% death rate (baseline)
#   Breast: protection 16.3x
#   Pancreatic: protection 21.9x
#   GBM: protection 26.4x (blood-brain barrier)
#   Melanoma: protection 17.5x
#   Sarcoma: protection 20.4x
# These estimated PK barriers reduce killing in the model; the factors do
# not establish the cause or magnitude of treatment failure in vivo.
```

**Spatial x Temporal:** compare each tumor's `death_rate` at distance 0 with
its rate at the furthest sampled distance in `tumor_pk_spatial_temporal.csv`.
Their ratio measures the additional spatial protection beyond temporal PK.
The default seed-42 run gives approximately 1.18x (Breast), 1.19x (Pancreatic),
1.05x (GBM), 1.14x (Melanoma), and 1.16x (Sarcoma), compared with 16.3–26.4x
protection from temporal PK alone. These are simulated outcomes under the
listed presets, not fixed properties of the model or measured tissue effects.
The binary prints the calculated rows; it no longer announces the obsolete
hard-coded 1.3–1.7x range before running them.

## Caveats

- All tumor PK parameters are ESTIMATED (no textbook coverage) -- the protection factors indicate relative magnitude, not precise predictions
- The RSL3 inactivation rate model (k_inact=0.015) is calibrated to match the sim-original Persister+RSL3 death rate (~42%) but uses a different mechanism than the CellState initialization model
- Plasma pharmacokinetics are simplified (IV bolus, first-order elimination) -- real PK includes distribution phases, protein binding, and metabolism
- The 180-step simulation represents a single dosing event, not repeated administration
- Tumor microenvironment factors (O2, immune, stromal, pH) are not included -- see sim-tme for those effects
- Blood-brain barrier (GBM) is modeled as reduced permeability only -- active efflux transporters (P-gp) are not explicitly modeled
- Composition of temporal C_i(t) with spatial Krogh decay uses metabolism-only lambda to avoid double-counting cellular uptake in the ODE, but this factorization is an approximation
