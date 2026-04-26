# sim-icd

ICD (immunogenic cell death) cascade comparison across treatments, demonstrating that physical modalities (SDT/PDT) produce more DAMP release per dead cell and stronger downstream immune activation than pharmacologic RSL3.

**Manuscript reference:** Chapter 7, Section 7.2

## What it models

1. **Ferroptosis biochemistry** -- full ROS/GSH/LP/GPX4/FSP1 cascade for each cell, measuring LP at death (not just alive/dead)
2. **DAMP release** -- lipid peroxidation level at death determines DAMP quantity (higher LP = more membrane oxidation products = more immunogenic signals)
3. **Dendritic cell activation** -- DAMP concentration drives DC maturation via Michaelis-Menten kinetics
4. **T cell priming** -- mature DCs prime cytotoxic T cells proportional to DC activation fraction
5. **Immune kill cascade** -- primed T cells kill remaining tumor cells; two modes compared (without and with anti-PD1 checkpoint blockade)
6. **Three phenotypes** -- Persister, OXPHOS, and Glycolytic cells tested separately

## Quick start

```bash
cd simulations
cargo build --release -p sim-icd
cargo run --release -p sim-icd -- --output-dir output/icd
```

Runtime: ~30-60 seconds (parallelized via rayon, 100K cells x 12 conditions).

## Parameters / CLI

| Argument | Default | Description |
|----------|---------|-------------|
| `--n-cells` | 100,000 | Cells per condition |
| `--seed` | 42 | Random seed |
| `--output-dir` | `output/icd` | Directory for output files |

Immune cascade parameters are hardcoded via `ImmuneParams::default()`. Biochemistry parameters via `Params::default()`.

## Output format

### `icd_comparison.json` -- full results

JSON array of 12 objects (3 phenotypes x 4 treatments):

```json
{
  "phenotype": "Persister",
  "treatment": "SDT",
  "n_cells": 100000,
  "n_dead": 99950,
  "death_rate": 0.9995,
  "avg_lp_at_death": 19.83,
  "immune_no_pd1": {
    "total_damps": 1981500.0,
    "damp_per_dead_cell": 19.83,
    "dc_activation_fraction": 0.92,
    "mature_dcs": 920.0,
    "primed_tcells": 4600.0,
    "immune_kills": 3200.0
  },
  "immune_with_pd1": {
    "total_damps": 1981500.0,
    "damp_per_dead_cell": 19.83,
    "dc_activation_fraction": 0.92,
    "mature_dcs": 920.0,
    "primed_tcells": 4600.0,
    "immune_kills": 12800.0
  }
}
```

| Field | Description |
|-------|-------------|
| avg_lp_at_death | Mean LP level at the moment of death (higher = more membrane damage = more DAMPs) |
| damp_per_dead_cell | DAMP release proportional to LP at death |
| dc_activation_fraction | Fraction of dendritic cells activated (0.0-1.0) |
| immune_kills | Number of additional tumor cells killed by immune response |

## Reproducing manuscript claims

**Chapter 7, Section 7.2 (ICD comparison):**
```bash
cargo run --release -p sim-icd -- --output-dir output/icd
# Inspect: output/icd/icd_comparison.json
# Key finding (printed to stderr):
#   SDT produces ~1.5-2x more DAMP per dead cell than RSL3
#   SDT+anti-PD1 immune kills >> RSL3+anti-PD1 immune kills
```

**Core claim: SDT >> RSL3 for immune activation:**
```bash
cargo run --release -p sim-icd 2>&1 | grep "SDT produces"
# Expected: SDT produces X.Xx more DAMP per dead cell than RSL3
```

## Caveats

- Immune cascade is a simplified pipeline model (DAMP -> DC -> T cell -> kill), not a spatially resolved immune simulation
- DAMP quantity is derived from LP at death, which is model-determined (the threshold-locking effect means LP at death is closer to the threshold than biology might predict)
- Anti-PD1 is modeled as a binary modifier on the PD-1/PD-L1 suppression fraction, not a dose-response curve
- No spatial effects -- immune response is computed from population-level DAMP totals
- No temporal dynamics of immune priming (days-to-weeks DC maturation not modeled)
- Immune parameters are estimated, not calibrated to a specific experimental system
