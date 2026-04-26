# sim-combo-mech

Mechanistic drug combination modeling that runs pairwise drug combinations through the ferroptosis biochemistry engine and computes Bliss-independence synergy scores, revealing WHY combinations are synergistic by tracing which pathway nodes each drug depletes.

**Manuscript reference:** Chapter 6, Section 6.3

## What it models

1. **Four drugs targeting distinct ferroptosis pathway nodes:**
   - **RSL3** -- GPX4 covalent inhibitor (92% inhibition)
   - **SDT** -- exogenous ROS burst (5.0 relative units)
   - **FSP1i** -- FSP1/CoQ10 backup repair inhibitor (85% inhibition)
   - **HDACi** -- HDAC inhibitor that doubles basal mitochondrial ROS (2x multiplier)
2. **Single-drug baselines** -- death rate for each drug alone (plus untreated control)
3. **Six pairwise combinations** -- all C(4,2) drug pairs
4. **Bliss independence scoring** -- expected additive effect = P(A) + P(B) - P(A)*P(B); synergy score = observed / Bliss prediction (>1.0 = synergistic, <1.0 = antagonistic)
5. **Pathway tracing** -- final GPX4, FSP1, GSH, and LP levels for each condition reveal which pathway nodes are depleted by combinations

## Quick start

```bash
cd simulations
cargo build --release -p sim-combo-mech
cargo run --release -p sim-combo-mech
```

Runtime: ~1-5 seconds (1,000 cells x 11 conditions, single-threaded).

## Parameters

All parameters are hardcoded (no CLI arguments).

| Parameter | Value | Description |
|-----------|-------|-------------|
| N_CELLS | 1,000 | Cells per condition |
| N_STEPS | 180 | Biochemistry timesteps per cell |
| Seed | 42 | Random seed |
| Phenotype | Persister (FSP1-low) | Cell type used for all conditions |
| Context | 2D culture | Default params (no MUFA protection) |

Drug effect parameters:

| Drug | gpx4_inhibition | fsp1_inhibition | exo_ros_dose | basal_ros_multiplier |
|------|----------------|----------------|-------------|---------------------|
| RSL3 | 0.92 | 0.0 | 0.0 | 1.0 |
| SDT | 0.0 | 0.0 | 5.0 | 1.0 |
| FSP1i | 0.0 | 0.85 | 0.0 | 1.0 |
| HDACi | 0.0 | 0.0 | 0.0 | 2.0 |

## Output format

Output directory: `output/combo-mech/`

### 1. `combo_synergy.csv` -- pairwise combination results

```csv
drug_a,drug_b,rate_a,rate_b,rate_combo,bliss_prediction,synergy_score,ci_low,ci_high,n_dead,n_cells,mean_gpx4_final,mean_fsp1_final,mean_gsh_final,mean_lp_final
RSL3,SDT,0.425,0.999,1.000,0.999,1.00,0.996,1.000,1000,1000,0.002,0.150,0.010,21.50
RSL3,FSP1i,0.425,0.312,0.849,0.604,1.41,0.823,0.872,849,1000,0.004,0.020,0.340,15.80
...
```

### 2. `combo_summary.json` -- full structured output

```json
{
  "phenotype": "Persister (FSP1-low)",
  "context": "2D culture (default params)",
  "n_cells_per_condition": 1000,
  "singles": [
    {
      "drug": "RSL3",
      "death_rate": 0.425,
      "ci_low": 0.395,
      "ci_high": 0.455,
      "n_dead": 425,
      "n_cells": 1000,
      "mean_gpx4_final": 0.024,
      "mean_fsp1_final": 0.150,
      "mean_gsh_final": 1.200,
      "mean_lp_final": 8.500
    }
  ],
  "combinations": [
    {
      "drug_a": "RSL3",
      "drug_b": "FSP1i",
      "rate_a": 0.425,
      "rate_b": 0.312,
      "rate_combo": 0.849,
      "bliss_prediction": 0.604,
      "synergy_score": 1.41,
      "ci_low": 0.823,
      "ci_high": 0.872,
      "n_dead": 849,
      "n_cells": 1000,
      "mean_gpx4_final": 0.004,
      "mean_fsp1_final": 0.020,
      "mean_gsh_final": 0.340,
      "mean_lp_final": 15.80
    }
  ]
}
```

| Field | Description |
|-------|-------------|
| synergy_score | Observed death rate / Bliss prediction (>1.1 = synergistic, <0.9 = antagonistic) |
| bliss_prediction | Expected additive effect assuming independent targets |
| mean_gpx4_final | Average remaining GPX4 activity (reveals target engagement) |
| mean_fsp1_final | Average remaining FSP1 activity (reveals backup pathway status) |

## Reproducing manuscript claims

**Chapter 6, Section 6.3 (drug combination synergy):**
```bash
cargo run --release -p sim-combo-mech
# stderr output includes "=== Top Synergistic Pairs ==="
# Expected: RSL3 + FSP1i synergy score ~1.4-2.0x (SYNERGISTIC)
# Manuscript claim: RSL3 + FSP1i = 1.99x Bliss synergy
# Mechanism: RSL3 disables GPX4 (primary repair), FSP1i disables FSP1 (backup)
#            --> no repair pathway remains --> runaway LP cascade
```

**Pathway trace (why RSL3+FSP1i is synergistic):**
```bash
# In combo_summary.json, compare:
#   RSL3 alone: mean_fsp1_final ~0.15 (backup still active)
#   FSP1i alone: mean_gpx4_final ~0.30 (primary still active)
#   RSL3+FSP1i: both near zero --> complete repair pathway collapse
```

## Caveats

- Drug potency parameters are estimated, not calibrated to specific IC50 values
- Synergy scores depend on potency assumptions -- directional findings (synergistic vs antagonistic) are more robust than exact scores
- All conditions use 2D culture params (no MUFA protection) -- in-vivo synergy may differ
- Bliss independence assumes drugs act on different targets; violations indicate shared pathway coupling, which is the point of the mechanistic model
- 1,000 cells per condition provides moderate statistical power; wider confidence intervals than the larger simulations
- HDACi "doubles basal ROS" is an estimate -- the referenced paper shows increased ROS without quantifying fold-change
- FSP1i 85% inhibition is estimated potency, not calibrated to a specific compound's IC50
