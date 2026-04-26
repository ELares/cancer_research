# sim-combo

Three-phase combination therapy optimizer that sweeps second-line treatment timing (0-28 days post-chemo) and anti-PD1 addition across SDT, PDT, and RSL3 to find the schedule that minimizes tumor survival.

**Manuscript reference:** Chapter 6, Section 6.2 + Chapter 7, Section 7.2

## What it models

1. **Phase 1: Chemotherapy** -- assumed to kill 90% of proliferating cells, leaving 1,000 persister survivors (not simulated; establishes starting conditions)
2. **Phase 2: Ferroptosis induction** -- persister cells with time-dependent defense recovery are treated with SDT, PDT, or RSL3 at variable delays (0-28 days post-chemo)
3. **Phase 3: Immune cascade** -- dead cells release DAMPs proportional to LP at death, triggering DC maturation, T cell priming, and immune killing of residual tumor
4. **Anti-PD1 modulation** -- each schedule tested with and without checkpoint blockade
5. **42 conditions total** -- 7 delay timepoints x 3 treatments x 2 anti-PD1 modes

## Quick start

```bash
cd simulations
cargo build --release -p sim-combo
cargo run --release -p sim-combo -- --output-dir output/combo
```

Runtime: ~30-60 seconds (parallelized via rayon, 50K cells x 42 conditions).

## Parameters / CLI

| Argument | Default | Description |
|----------|---------|-------------|
| `--n-cells` | 50,000 | Cells per condition (for ferroptosis simulation) |
| `--seed` | 42 | Random seed |
| `--output-dir` | `output/combo` | Directory for output files |

Treatment delays swept: 0, 1, 3, 7, 14, 21, 28 days post-chemo.

Initial tumor population (post-chemo survivors): 1,000 persister cells.

Recovery rates, immune parameters, and biochemistry parameters are all hardcoded via their respective `Default` implementations.

## Output format

### `combo_results.json` -- full results

JSON array of 42 objects:

```json
{
  "sdt_delay_days": 0.0,
  "treatment": "SDT",
  "with_anti_pd1": true,
  "initial_tumor_cells": 1000,
  "ferroptosis_kill_rate": 0.999,
  "ferroptosis_killed": 999,
  "immune_kills": 0.8,
  "total_killed": 999,
  "survivors": 1,
  "survival_fraction": 0.001,
  "total_damps": 19830.0,
  "damp_per_dead_cell": 19.83,
  "primed_tcells": 46.0
}
```

| Field | Description |
|-------|-------------|
| sdt_delay_days | Days between chemo withdrawal and ferroptosis treatment |
| ferroptosis_kill_rate | Fraction killed by the ferroptosis modality (from 50K simulation) |
| ferroptosis_killed | Scaled to the biological population (1,000 cells) |
| immune_kills | Additional kills from the immune cascade |
| survivors | Remaining tumor cells after both phases |
| survival_fraction | survivors / initial_tumor_cells |

The binary also identifies and prints the **optimal schedule** (minimum survival fraction) to stderr.

## Reproducing manuscript claims

**Chapters 6.2 + 7.2 (optimal SDT timing):**
```bash
cargo run --release -p sim-combo -- --output-dir output/combo
# stderr output: "=== Optimal Schedule ==="
# Expected: SDT at day 0 + anti-PD1 achieves near-zero survivors
# Key finding: SDT efficacy declines more slowly than RSL3 over the 28-day window
```

**Schedule comparison:**
```bash
cargo run --release -p sim-combo 2>&1 | grep "Day"
# Inspect survival fraction at each delay x treatment x anti-PD1 combination
# Expected: RSL3 survival fraction increases sharply beyond day 3-7
# Expected: SDT survival fraction remains low through day 14-21
```

## Caveats

- Phase 1 (chemotherapy) is not simulated -- the 90% kill rate and clean persister selection are assumed
- Immune cascade is population-level, not spatially resolved
- Dead cell LP values are sampled from the 50K simulation and scaled to the 1,000-cell biological population
- Recovery kinetics use exponential half-life models
- Anti-PD1 is modeled as a binary modifier, not a dose-response curve
- No repeated dosing -- each schedule represents a single treatment event
- The "optimal schedule" depends on parameter choices that are literature-estimated, not calibrated
