# sim-invivo

2D-versus-in-vivo ferroptosis comparison demonstrating that SCD1-driven MUFA lipid remodeling protects tumor cells from pharmacologic ferroptosis inducers (RSL3) in vivo while physical ROS modalities (SDT/PDT) can still overwhelm the defense.

**Manuscript reference:** Chapter 7, Section 7.1

## What it models

1. **Three biological contexts** compared with identical cell phenotypes:
   - **2D baseline** -- MUFA pathway off (standard culture conditions)
   - **In-vivo** -- SCD1-driven MUFA membrane protection active (lipid remodeling)
   - **In-vivo + SCD1 inhibitor** -- SCD1 blocked (rate=0) but cells retain pre-existing membrane MUFA that decays over time
2. **MUFA protection mechanics** -- SCD1 converts saturated fatty acids to monounsaturated fatty acids, replacing oxidizable PUFAs in the membrane and reducing lipid peroxidation susceptibility
3. **Four phenotypes x four treatments** -- full factorial (Glycolytic, OXPHOS, Persister, Persister+NRF2) x (Control, RSL3, SDT, PDT) in each context
4. **Protection factor analysis** -- ratio of 2D death rate to in-vivo death rate quantifies the MUFA protection magnitude
5. **MUFA parameter sweep** -- systematic exploration of SCD1 rate x MUFA max across two modes (onset with mufa=0, and steady-state) for Persister cells with SDT and RSL3

## Quick start

```bash
cd simulations
cargo build --release -p sim-invivo
cargo run --release -p sim-invivo -- --output-dir output/invivo
```

Runtime: ~2-5 minutes (parallelized via rayon, 100K cells main + 50K cells per sweep point).

## Parameters / CLI

| Argument | Default | Description |
|----------|---------|-------------|
| `--n-cells` | 100,000 | Cells per condition (main comparison) |
| `--n-sweep` | 50,000 | Cells per condition (MUFA parameter sweep) |
| `--seed` | 42 | Random seed |
| `--output-dir` | `output/invivo` | Directory for output files |

In-vivo MUFA parameters are hardcoded via `Params::invivo()`:
- `scd_mufa_rate` -- SCD1 enzyme activity (MUFA synthesis rate per step)
- `scd_mufa_max` -- maximum MUFA membrane fraction
- `scd_mufa_decay` -- natural lipid turnover rate

MUFA sweep grid: rates [0.002, 0.005, 0.01, 0.02, 0.04] x max [0.20, 0.30, 0.40, 0.50, 0.60].

## Output format

### 1. `invivo_comparison.json` -- three-context comparison

JSON array of 48 objects (3 contexts x 4 phenotypes x 4 treatments):

```json
{
  "context": "invivo",
  "phenotype": "Persister (FSP1\u2193)",
  "treatment": "RSL3",
  "n_cells": 100000,
  "n_dead": 5200,
  "death_rate": 0.052,
  "ci_low": 0.050,
  "ci_high": 0.054,
  "mean_lp": 3.21,
  "mean_gsh": 2.85,
  "mean_gpx4": 0.04,
  "scd_mufa_rate": 0.01,
  "scd_mufa_max": 0.40,
  "scd_mufa_decay": 0.002,
  "initial_mufa_protection": 0.0
}
```

### 2. `mufa_sweep.json` and `mufa_sweep.csv` -- parameter sweep

```json
{
  "scd_mufa_rate": 0.01,
  "scd_mufa_max": 0.40,
  "scd_mufa_decay": 0.002,
  "initial_mufa_protection": 0.0,
  "phenotype": "Persister (FSP1\u2193)",
  "treatment": "RSL3",
  "n_cells": 50000,
  "death_rate": 0.052,
  "ci_low": 0.048,
  "ci_high": 0.056,
  "protection_factor": 8.2
}
```

| Field | Description |
|-------|-------------|
| protection_factor | 2D baseline death rate / sweep-point death rate (higher = more MUFA protection) |

## Reproducing manuscript claims

**Chapter 7, Section 7.1 (RSL3 loses efficacy in vivo):**
```bash
cargo run --release -p sim-invivo -- --output-dir output/invivo
# stderr output includes "Key Biological Predictions"
# Expected:
#   Dixon 2025 (RSL3 fails in vivo): 2D=~42% -> InVivo=~5% -- CONFIRMED
#   SDT survives MUFA defense: InVivo >50% -- YES
#   SCD1i resensitizes (Tesfay 2019): InVivo=~5% -> SCD1i=~30% -- CONFIRMED
```

**Protection factor table:**
```bash
cargo run --release -p sim-invivo 2>&1 | grep "Persister.*RSL3"
# Expected: RSL3 protection factor ~8x (2D ~42% vs in-vivo ~5%)
# Expected: SDT protection factor ~1.1-1.3x (minimal MUFA impact)
```

## Caveats

- MUFA protection is modeled as a continuous variable reducing LP accumulation rate -- the actual membrane remodeling involves discrete lipid species and asymmetric leaflet composition
- SCD1 inhibitor context starts with pre-existing MUFA that decays -- this models acute SCD1i administration, not long-term treatment
- In-vivo parameters (SCD1 rate, MUFA max, decay) are literature-estimated, not calibrated to a specific cell line or tumor model
- No spatial effects (uniform drug distribution)
- No immune component in this simulation
- The sweep explores a 5x5 grid of MUFA parameters; the actual parameter space is larger and includes interactions with other pathways not swept here
