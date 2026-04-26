# sim-original

Monte Carlo ferroptosis sensitivity baseline. Simulates 1,000,000 cells across 16 conditions (4 phenotypes x 4 treatments) to establish the core finding: persister cells are selectively vulnerable to ferroptosis, and physical ROS modalities (SDT/PDT) kill across all phenotypes while pharmacologic RSL3 is persister-selective.

**Manuscript reference:** Chapter 5, Section 5.2 (Figure 7)

## What it models

The ferroptosis biochemical cascade for individual cells:

1. **ROS generation** — basal mitochondrial ROS + treatment-specific exogenous ROS (SDT/PDT) or GPX4 inhibition (RSL3)
2. **GSH depletion** — Michaelis-Menten scavenging consumes glutathione
3. **Lipid peroxidation** — direct ROS damage + autocatalytic propagation (bistable switch)
4. **Repair pathways** — GPX4 (primary) + FSP1 (backup, downregulated in persisters)
5. **Death** — LP crosses threshold (10.0) → irreversible membrane damage

Each cell has stochastic parameter variation (±20%) around phenotype-specific means.

## Quick start

```bash
cd simulations
cargo build --release -p sim-original
cargo run --release -p sim-original > simulation_results.json
```

Runtime: ~20 seconds (16M cells, parallelized via rayon).

## Parameters

All parameters are hardcoded via `Params::default()`. No command-line arguments.

| Parameter | Default | Role |
|-----------|---------|------|
| fenton_rate | 0.02 | Basal ROS from iron-catalyzed Fenton reaction |
| gsh_scav_efficiency | 0.5 | GSH scavenging rate (Michaelis-Menten) |
| lp_rate | 0.06 | Direct lipid peroxidation from unscavenged ROS |
| lp_propagation | 0.10 | Autocatalytic LP chain reaction rate (bistable switch gate) |
| gpx4_rate | 0.30 | GPX4 enzymatic repair of lipid peroxides |
| fsp1_rate | 0.08 | FSP1/CoQ10 backup repair pathway |
| nrf2_gsh_rate | 0.025 | NRF2-driven GSH resynthesis |
| gpx4_degradation_by_ros | 0.002 | Oxidative degradation of GPX4 |
| sdt_ros | 5.0 | Exogenous ROS dose from sonodynamic therapy |
| pdt_ros | 5.0 | Exogenous ROS dose from photodynamic therapy |
| rsl3_gpx4_inhib | 0.92 | Fraction of GPX4 inhibited by RSL3 (92%) |
| death_threshold | 10.0 | LP level triggering cell death |
| post_death_steps | 5 | Steps of continued LP accumulation after death (for DAMP calculation) |

See `simulations/calibration/parameter_provenance.md` for literature sources and sensitivity ratings.

## Phenotypes

| Phenotype | Description | Key difference |
|-----------|-------------|---------------|
| Glycolytic | Standard cancer cell | High FSP1 (1.0), low basal ROS (0.2) |
| OXPHOS | Oxidative metabolism | Higher basal ROS (0.5), moderate FSP1 |
| Persister (FSP1 down) | Drug-tolerant persister | Low FSP1 (0.15), moderate basal ROS (0.25) |
| Persister + NRF2 | NRF2-compensated persister | Low FSP1 + high NRF2 (2.0x GSH recovery) |

## Treatments

| Treatment | Mechanism | Model effect |
|-----------|-----------|-------------|
| Control | No treatment | Only basal ROS drives LP |
| RSL3 | GPX4 covalent inhibitor | GPX4 activity reduced to 8% of normal |
| SDT | Sonodynamic therapy | Exogenous ROS burst (5.0 relative units) |
| PDT | Photodynamic therapy | Exogenous ROS burst (5.0 relative units) |

## Output format

JSON array of 16 objects (one per phenotype x treatment):

```json
{
  "phenotype": "Persister (FSP1↓)",
  "treatment": "SDT",
  "n_cells": 1000000,
  "n_dead": 999500,
  "death_rate": 0.9995,
  "ci_low": 0.9994,
  "ci_high": 0.9996,
  "mean_lipid_perox": 19.83,
  "mean_gsh_final": 0.03,
  "mean_gpx4_final": 0.12
}
```

| Field | Description |
|-------|-------------|
| death_rate | Fraction of cells killed (0.0-1.0) |
| ci_low, ci_high | Wilson 95% confidence interval |
| mean_lipid_perox | Average final LP (higher = more membrane damage) |
| mean_gsh_final | Average remaining GSH (lower = more depleted) |
| mean_gpx4_final | Average remaining GPX4 activity |

## Built-in validation

The binary checks four biological constraints after simulation:

1. **Control < 2% death** for all phenotypes (baseline viability)
2. **RSL3 selectively kills persisters** (Persister death >> Glycolytic death)
3. **SDT kills all phenotypes** (> 85% death across the board)
4. **NRF2 protects against RSL3 but not SDT** (NRF2 compensation failure mode)

## Sensitivity analysis

After the main simulation, the binary runs a ±50% perturbation on 11 rate constants (100K cells per perturbation, 22 conditions total). It tests whether "Persister > Glycolytic under SDT" holds across all perturbations. Result: 22/22 conditions hold.

The output is printed to stderr. A more comprehensive global sensitivity analysis (PRCC) is available via `scripts/run_prcc.py`.

## Reproducing manuscript claims

**Chapter 5, Figure 7:**
```bash
cargo run --release -p sim-original > simulation_results.json
# Parse: death_rate for each phenotype × treatment
# Expected: Persister×RSL3 ≈ 42%, Persister×SDT ≈ 100%, Glycolytic×RSL3 ≈ 0%
```

**Sensitivity analysis (Chapter 5, Finding 5):**
```bash
cargo run --release -p sim-original 2>&1 | grep "Result holds"
# Expected: "Result holds in 22/22 conditions (100%)"
```

## Caveats

- All parameters are estimated from literature, not calibrated to a specific cell line
- The 180-step simulation represents a single treatment window, not repeated dosing
- No spatial effects (uniform drug distribution) — see sim-spatial for depth-dependent modeling
- No in-vivo lipid remodeling (MUFA/SCD1) — see sim-invivo for 3D context
