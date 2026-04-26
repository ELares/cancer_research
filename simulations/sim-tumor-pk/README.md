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

Radial distance bins for C(r,t): [0, 25, 50, 75, 100, 125] um.

## Output format

Output directory: `output/tumor-pk/`

### 1. `tumor_pk_summary.json` -- per-scenario results

JSON array with one entry per scenario (2D reference + 5 tumor types):

```json
{
  "tumor_type": "Breast",
  "context": "tumor_pk",
  "n_cells": 10000,
  "n_dead": 850,
  "death_rate": 0.085,
  "ci_low": 0.079,
  "ci_high": 0.091,
  "mean_lp": 4.20,
  "mean_gsh": 2.10,
  "mean_gpx4": 0.05,
  "peak_c_interstitial": 0.180,
  "auc_c_interstitial": 12.5,
  "protection_factor": 4.8
}
```

| Field | Description |
|-------|-------------|
| peak_c_interstitial | Maximum interstitial drug concentration (0.0-1.0 normalized) |
| auc_c_interstitial | Area under the interstitial concentration curve |
| protection_factor | 2D reference death rate / tumor death rate (higher = more PK barrier) |

### 2. `tumor_pk_timecourse.csv` -- concentration time series

```csv
time_min,tumor_type,c_plasma,c_vascular,c_interstitial
0,Breast,1.000000,0.120000,0.012000
1,Breast,0.977000,0.118000,0.014000
...
```

### 3. `tumor_pk_spatial_temporal.csv` -- C(r,t) kill rates

```csv
tumor_type,distance_um,peak_conc,death_rate,ci_low,ci_high,n_cells,n_dead
Breast,0,0.180,0.085,0.079,0.091,10000,850
Breast,25,0.160,0.072,0.067,0.078,10000,720
...
```

## Reproducing manuscript claims

**Chapter 8, Section 8.2 (2D-to-in-vivo gap):**
```bash
cargo run --release -p sim-tumor-pk
# stderr output includes "=== Protection Factor Summary ==="
# Expected:
#   2D culture ref: ~41% death rate (baseline)
#   Breast: protection ~4-5x
#   Pancreatic: protection ~16-27x
#   GBM: protection ~20-30x (blood-brain barrier)
#   Melanoma: protection ~3-5x
#   Sarcoma: protection ~8-15x
# Claim: 2D-to-in-vivo gap demonstrates why pharmacologic ferroptosis inducers
#        fail in vivo (protection factors of 3-30x from PK barriers alone)
```

**Spatial x Temporal (temporal PK dominates spatial decay):**
```bash
cargo run --release -p sim-tumor-pk 2>&1 | grep "Key finding"
# Expected: "temporal PK barrier (16-27x) dominates spatial decay (1.3-1.7x)"
# The PK barrier (getting drug to the tumor interstitium) is far more
# limiting than the radial diffusion gradient within tissue
```

## Caveats

- All tumor PK parameters are ESTIMATED (no textbook coverage) -- the protection factors indicate relative magnitude, not precise predictions
- The RSL3 inactivation rate model (k_inact=0.015) is calibrated to match the sim-original Persister+RSL3 death rate (~42%) but uses a different mechanism than the CellState initialization model
- Plasma pharmacokinetics are simplified (IV bolus, first-order elimination) -- real PK includes distribution phases, protein binding, and metabolism
- The 180-step simulation represents a single dosing event, not repeated administration
- Tumor microenvironment factors (O2, immune, stromal, pH) are not included -- see sim-tme for those effects
- Blood-brain barrier (GBM) is modeled as reduced permeability only -- active efflux transporters (P-gp) are not explicitly modeled
- Composition of temporal C_i(t) with spatial Krogh decay uses metabolism-only lambda to avoid double-counting cellular uptake in the ODE, but this factorization is an approximation
