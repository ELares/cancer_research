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
  "avg_final_lp": 19.83,
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
| avg_final_lp | Mean LP after all 180 steps, over the cells that died. **NOT LP at the moment of death**, which is ~10.25 for every arm against a `death_threshold` of 10.0 -- death IS the threshold crossing. This field was called `avg_lp_at_death` and described here as "LP level at the moment of death" until #928; the two quantities differ by 78% (RSL3) and 86% (SDT) WITHIN one run, measured per cell in `analysis/immune-2d-measurements/observations.json.gz` (death 10.2474/10.2657 against release 18.2658/19.1230). `null` when no cell died. |
| damp_per_dead_cell | `avg_final_lp` x `damp_per_lp`, which is 1.0 and uncalibrated. `null` when no cell died (it was `0.0`, which read as a measured potency of zero). |
| dc_activation_fraction | Fraction of dendritic cells activated (0.0-1.0) |
| immune_kills | Number of additional tumor cells killed by immune response. **Capped at the surviving population**, so an arm that kills everything biochemically earns almost none -- see the measured inversion below. |

## Reproducing manuscript claims

**The per-dead-cell DAMP comparison:**
```bash
cargo run --release -p sim-icd -- --output-dir output/icd
# Inspect: output/icd/icd_comparison.json
# Key finding (printed to stderr): SDT produces 1.133x the DAMP per dead cell
# of RSL3 at the persister state (20.0846 vs 17.7240), and 1.075x at OXPHOS.
```

**THE PREDICTION THIS SECTION USED TO STATE WAS REFUTED BY THE BINARY BESIDE
IT (#928).** Where the block above now gives the measured values, it used to
read:

```
#   SDT produces ~1.5-2x more DAMP per dead cell than RSL3
#   SDT+anti-PD1 immune kills >> RSL3+anti-PD1 immune kills
```

Neither reproduces. The ratio is **1.047x to 1.148x** across all three regimes
the repository can run -- this binary, `sim-modality-panel --tme-sweep` over 32
readable conditions, and the committed 2D spatial archive -- so the direction
holds everywhere and the magnitude reaches no part of the 1.5-2x band. The
second line inverts outright: SDT earns **3** immune kills against RSL3's
**57,591**, because SDT kills 99.997% of persisters and `immune_cascade` caps
kills at the surviving population. The lines are kept rather than deleted
because a stated expectation that failed is evidence, and `scripts/icd_damp_per_death.py`
pins this quote so the measurement and the claim cannot drift apart.

Nothing in the manuscript rests on the refuted band; Section 8.2 already
declines the per-cell reading. Full accounting:
`analysis/calibration/icd-damp-per-death.md`.

**Direction check (not a magnitude check):**
```bash
cargo run --release -p sim-icd 2>&1 | grep "SDT produces"
```

## Caveats

- Immune cascade is a simplified pipeline model (DAMP -> DC -> T cell -> kill), not a spatially resolved immune simulation
- DAMP quantity is derived from the post-death LP this binary reports, which is model-determined. The threshold-locking effect named in this list is not a minor caveat but the reason the documented contrast fails: LP at the moment of death is pinned near the threshold for EVERY arm, so a per-arm difference exists only after death (#928)
- Anti-PD1 is modeled as a binary modifier on the PD-1/PD-L1 suppression fraction, not a dose-response curve
- No spatial effects -- immune response is computed from population-level DAMP totals
- No temporal dynamics of immune priming (days-to-weeks DC maturation not modeled)
- Immune parameters are estimated, not calibrated to a specific experimental system
