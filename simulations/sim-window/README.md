# sim-window

Vulnerability window simulation modeling how persister cell ferroptosis sensitivity changes over 0-28 days after chemotherapy withdrawal as defense pathways (FSP1, GPX4, NRF2, GSH) progressively recover.

**Manuscript reference:** Chapter 6, Section 6.2

## What it models

1. **Time-dependent defense recovery** -- persister cells regain FSP1, GPX4, NRF2, and GSH activity via exponential half-recovery kinetics after chemotherapy withdrawal
2. **Ferroptosis sensitivity at 9 timepoints** -- 0h, 6h, 12h, 24h, 48h, 72h, 1 week, 2 weeks, 4 weeks
3. **Treatment comparison** -- Control, RSL3, SDT, and PDT applied at each timepoint to assess which modalities retain efficacy as defenses recover
4. **Recovery rate sensitivity analysis** -- +/-50% perturbation on each recovery half-life (FSP1, GPX4, NRF2, GSH) tested at 7 days with SDT

## Quick start

```bash
cd simulations
cargo build --release -p sim-window
cargo run --release -p sim-window -- --output-dir output/window
```

Runtime: ~30-60 seconds (parallelized via rayon, 100K cells x 36 conditions + sensitivity).

## Parameters / CLI

| Argument | Default | Description |
|----------|---------|-------------|
| `--n-cells` | 100,000 | Cells per condition |
| `--seed` | 42 | Random seed |
| `--output-dir` | `output/window` | Directory for output files |

Recovery half-lives are hardcoded in `RecoveryRates::default()`:
- FSP1 half-recovery: days
- GPX4 half-recovery: days
- NRF2 half-recovery: days
- GSH half-recovery: days

## Output format

### 1. `vulnerability_window.json` -- full results

JSON array of 36 objects (9 timepoints x 4 treatments):

```json
{
  "timepoint_hours": 168.0,
  "timepoint_days": 7.0,
  "treatment": "SDT",
  "n_cells": 100000,
  "n_dead": 95200,
  "death_rate": 0.952,
  "ci_low": 0.950,
  "ci_high": 0.954,
  "mean_lp": 18.5,
  "mean_gsh": 0.12,
  "mean_gpx4": 0.08
}
```

| Field | Description |
|-------|-------------|
| timepoint_hours | Hours after chemotherapy withdrawal |
| timepoint_days | Days after withdrawal |
| death_rate | Fraction of cells killed (0.0-1.0) |
| ci_low, ci_high | Wilson 95% confidence interval |
| mean_lp | Average final lipid peroxidation level |
| mean_gsh | Average final GSH level |
| mean_gpx4 | Average final GPX4 activity |

### 2. `vulnerability_window.csv` -- tabular summary

```csv
hours,treatment,death_rate,ci_low,ci_high
0.0,Control,0.003,0.002,0.004
0.0,RSL3,0.425,0.422,0.428
0.0,SDT,0.999,0.998,0.999
...
```

### 3. Sensitivity analysis (stderr)

Recovery rate +/-50% perturbation results printed to stderr for SDT at 7 days.

## Reproducing manuscript claims

**Chapter 6, Section 6.2 (vulnerability window):**
```bash
cargo run --release -p sim-window -- --output-dir output/window
# Inspect: output/window/vulnerability_window.json
# Expected: RSL3 death_rate drops substantially within days (window closes quickly)
# Expected: SDT death_rate remains high (>80%) through at least 2 weeks
# Claim: RSL3 vulnerability window closes in days, SDT persists weeks
```

**Sensitivity check:**
```bash
cargo run --release -p sim-window 2>&1 | grep "t.* SDT death"
# Expected: SDT death rate at 7 days remains >80% across all +/-50% perturbations
```

## Caveats

- Recovery kinetics use exponential half-life models -- real recovery may involve thresholds or feedback loops
- No spatial effects (uniform drug distribution at all timepoints)
- No in-vivo lipid remodeling (MUFA/SCD1 not modeled)
- Persister cells are generated independently at each timepoint (no population memory)
- Recovery rate parameters are estimated from literature, not calibrated to a specific cell line
- The model assumes chemotherapy cleanly selects for persisters; real tumor evolution is more complex
