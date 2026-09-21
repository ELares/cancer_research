# sim-tme

Tumor microenvironment simulation integrating four spatially-resolved features -- oxygen gradients, spatial immune coupling (DAMP diffusion + T cell activation), CAF-mediated stromal protection (GSH/MUFA supply), and pH gradients (iron release + drug ion trapping) -- on a 500x500 cell grid to explore how these model assumptions affect pharmacologic and physical ROS treatments.

**Manuscript reference:** Chapters 7.1-7.5

## What it models

1. **Oxygen gradients (Feature A)** -- exponential O2 decay from tumor edge inward: O2(d) = exp(-d/lambda). O2 modulates basal ROS via mitochondrial ETC. Sweep over lambda = [80, 100, 120, 150] um covers the literature range (Vaupel 1989). Includes square-wave O2 cycling (normoxic/hypoxic alternation, period=60 steps)
2. **Spatial immune coupling (Feature B)** -- ferroptotic cells release DAMPs proportional to LP after the five-step post-death grace period; DAMPs diffuse to Moore neighbors and decay exponentially; local DAMP concentration drives Michaelis-Menten DC activation; activated DCs yield immune kills with PD-1 brake. Tested with and without anti-PD1
3. **CAF-mediated stromal protection (Feature C)** -- cancer-associated fibroblasts supply GSH (cysteine via GGT1/SLC7A11) and MUFA (oleic acid via ACSL3) to stromal-adjacent tumor cells, boosting antioxidant capacity and membrane protection
4. **pH gradient and ion trapping (Feature D)** -- glycolytic tumors produce lactic acid (Warburg effect), creating pH 7.4 (edge) to 6.5 (core). Two competing effects: ferritin destabilization releases Fe2+ at low pH (pro-ferroptosis), and weak-base drug ion trapping reduces intracellular RSL3 (anti-ferroptosis, Henderson-Hasselbalch)
5. **Zone analysis** -- kill rates computed for three anatomical zones: normoxic shell, transition zone, hypoxic core (boundaries fixed at reference lambda=120 um)

## Quick start

```bash
cd simulations
cargo build --release -p sim-tme
cargo run --release -p sim-tme
```

Runtime: 10-30 minutes (500x500 grid x multiple conditions x 180 steps each, many runs).

## Parameters

The default matrix uses the values below. `--sdt-o2-dependence` and `FERRO_*` overrides support separate sensitivity runs; `--immune-measurements` accepts only the canonical settings (see below).

| Parameter | Value | Description |
|-----------|-------|-------------|
| GRID_SIZE | 500 | Grid rows = cols |
| CELL_SIZE_UM | 20.0 | Cell diameter in micrometers |
| N_STEPS | 180 | Biochemistry timesteps |
| SEED | 42 | Random seed |
| O2_LAMBDAS | [80, 100, 120, 150] | O2 penetration lengths (um) to sweep |
| ZONE_REF_LAMBDA | 120.0 | Fixed reference lambda for zone boundaries |

Key physics/biology parameters:

| Feature | Parameter | Value | Source |
|---------|-----------|-------|--------|
| O2 | penetration lambda | 80-150 um | Vaupel 1989 |
| Immune | DAMP diffusion fraction | 0.08/step | Estimated |
| Immune | DAMP clearance rate | 0.03/step | Estimated |
| Immune | DC activation Kd | 50.0 | Estimated |
| Immune | immune kill rate | 0.02/step | Estimated |
| Immune | PD-1 brake | 0.70 (70% suppression) | Estimated |
| Immune | anti-PD1 efficacy | 0.80 (80% brake removed) | Estimated |
| Stromal | GSH boost | 0.06/step (cap 18.0) | Estimated; PMID 34373744 |
| Stromal | MUFA boost | 0.003/step (cap 0.25) | Estimated; PMID 31813804 |
| pH | edge/core | 7.4/6.5 | Stubbs 2000 |
| pH | iron sensitivity | 1.5 | Estimated; Harrison & Arosio 1996 |
| pH | ion trap sensitivity | 0.4 | Estimated; Chemistry2e Sec.14.2 |

## Output format

Output directory: `output/tme/`

### 1. `tme_summary.json` -- all conditions

The `{"schema_version": 1, "conditions": [...]}` envelope contains one entry per condition. This baseline immune-on SDT row is reproduced from the [historical fixture](../../tests/fixtures/flagship_tme_rows.json); inactive optional pH fields are omitted:

```json
{
  "treatment": "SDT",
  "o2_condition": "gradient_120um",
  "o2_lambda_um": 120.0,
  "immune_mode": "immune_on",
  "total_tumor": 159033,
  "total_dead": 140161,
  "ferroptosis_kills": 139640,
  "immune_kills": 521,
  "overall_kill_rate": 0.8813328051410714,
  "normoxic_kill_rate": 0.9547834364588291,
  "transition_kill_rate": 0.9416085657370518,
  "hypoxic_kill_rate": 0.8695500315832497,
  "stromal_mode": "off",
  "stromal_adjacent_kill_rate": 0.9610244988864143,
  "stromal_adjacent_count": 1796
}
```

### 2. `depth_kill_curves.csv` -- depth-kill profiles for all conditions

### 3. Heatmap CSVs (500x500 matrices, u8-encoded):
- `death_{treatment}_uniform.csv` -- baseline death maps
- `death_{treatment}_o2gradient.csv` -- O2 gradient death maps (lambda=120)
- `death_{treatment}_immune_run.csv` -- immune-coupled death maps
- `death_{treatment}_stromal.csv` -- stromal protection death maps
- `death_{treatment}_ph.csv` -- pH gradient death maps
- `o2_field.csv` -- O2 concentration field (0-255)
- `damp_field_{treatment}.csv` -- DAMP concentration field (0-255)
- `stromal_mask.csv` -- stromal adjacency mask (0=stromal, 1=adjacent tumor, 2=non-adjacent)
- `ph_field.csv` -- pH field (0-255 mapping pH 6.5-7.4)

## Reproducing manuscript claims

**Chapter 7.1 (hypoxia collapses RSL3):**
```bash
cargo run --release -p sim-tme 2>&1 | grep "Sensitivity: SDT advantage"
# Expected: RSL3 hypoxic kill rate drops dramatically with O2 gradient
# SDT/RSL3 ratio in hypoxic zone consistently >1 across all lambda values
# Claim: hypoxia collapses RSL3 efficacy
```

**Chapter 8.2 (historical 2D total immune-kill ratio):**
```bash
cargo run --release -p sim-tme 2>&1 | grep "immune="
# Compare SDT immune_kills vs RSL3 immune_kills
# Canonical gradient-120, immune_on, no stromal/pH: SDT 521; RSL3 5
# 521/5 = 104.2 compares total immune kills, not DAMP or potency per dead cell
```

**Chapter 7.3 (stromal halves RSL3):**
```bash
cargo run --release -p sim-tme 2>&1 | grep "stromal_adj="
# Compare stromal_adjacent_kill_rate with stromal ON vs OFF
# Expected: RSL3 stromal-adjacent kill rate roughly halved
# SDT less affected by stromal protection
```

**Chapter 7.4 (pH halves RSL3):**
```bash
cargo run --release -p sim-tme 2>&1 | grep "ph_on\|pH"
# RSL3 with pH gradient: ion trapping reduces drug availability in acidic core
# Expected: RSL3 kill rate reduced ~50% in hypoxic/acidic zone
# SDT unaffected by pH (physical ROS, not a weak-base drug)
```

## Caveats

- O2 modulates basal_ros only (conservative) -- does not affect Fenton reaction rate or SDT/PDT energy deposition
- SDT/PDT modeled as O2-independent (Type I sonochemical mechanism assumption; if Type II O2-dependent mechanism dominates, SDT advantage in hypoxia would be smaller)
- DAMP release uses the modeled LP at the end of the five-step grace period; deaths with incomplete grace periods add DAMP after the last immune update. Neither LP nor DAMP has an experimentally calibrated unit here.
- The historical scope description says the immune model captures only the resident T cell phase (0-48h), not systemic lymph node priming (1-7 days); these labels do not calibrate the biochemical step duration (see [time mapping](../calibration/parameter_provenance.md)).
- DAMP clearance modeled as exponential decay (simplified)
- All stromal (CAF) and pH parameters are ESTIMATED -- no textbook coverage exists for CAF biology
- RSL3 pKa is not well-characterized (chloroacetamide, not a classic weak base) -- ion_trap_sensitivity is the most uncertain parameter
- Warburg effect not covered in reference textbooks (Biology2e describes lactic acid fermentation only in anaerobic contexts)
- No drug penetration model (RSL3 is uniform) -- see sim-tissue-pk and sim-tumor-pk for transport effects

## Passive measurement of the canonical immune comparison

```bash
# From simulations/: observe the three existing baseline immune-on arms.
cargo run --locked --release -p sim-tme -- --immune-measurements
```

This writes `output/tme/immune_measurements.json` in a separate mode. It records
death and release events, incomplete grace periods, actual immune kills, unique
eligible cells and repeated cell-step opportunities. Eligible exposure includes
the nonlinear activation sum and counts above the half-maximal and 90%-activation
DAMP thresholds. The observer neither modifies cell state nor draws random numbers.
Other arguments and all `FERRO_*` overrides are rejected in this mode.

For a pinned capture that also verifies the full default matrix, or for offline
reconstruction without running Rust, use the
[2D measurement protocol](../../docs/IMMUNE_2D_MEASUREMENT_PROTOCOL.md). Results
describe this fixed model realization; they do not isolate a causal mechanism,
calibrate an experimental immune readout, or turn cell-steps into physical time.
