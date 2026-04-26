# sim-spatial

2D heterogeneous tumor simulation with depth-dependent energy deposition from PDT (Beer-Lambert optical attenuation), SDT (acoustic attenuation), and RSL3 (uniform pharmacologic distribution), establishing that ultrasound penetrates centimeters while light penetrates millimeters.

**Manuscript reference:** Chapter 6, Section 6.1 (Figure 8)

## What it models

1. **2D tumor grid** -- heterogeneous cell composition (glycolytic, OXPHOS, persister, persister+NRF2, stromal) generated deterministically from seed
2. **Depth-dependent energy deposition** -- PDT follows Beer-Lambert exponential decay (millimeter-scale penetration); SDT follows acoustic attenuation (centimeter-scale penetration); RSL3 distributes uniformly
3. **Ferroptosis biochemistry per cell** -- full ROS/GSH/LP/GPX4/FSP1 cascade with stochastic variation, depth-adjusted exogenous ROS
4. **Iron diffusion** -- newly dead cells release labile iron to Moore neighbors, amplifying local Fenton ROS (bystander effect)
5. **Death heatmaps** -- spatial pattern of cell death across the tumor cross-section

## Quick start

```bash
cd simulations
cargo build --release -p sim-spatial
cargo run --release -p sim-spatial -- --output-dir output/spatial
```

Runtime: several minutes (500x500 grid x 4 treatments x 180 steps, single-threaded per cell).

## Parameters / CLI

| Argument | Default | Description |
|----------|---------|-------------|
| `--grid-size` | 500 | Grid dimensions (rows = cols) |
| `--cell-size` | 20.0 | Cell diameter in micrometers |
| `--seed` | 42 | Random seed (same seed = same tumor topology) |
| `--output-dir` | `output/spatial` | Directory for output files |
| `--n-steps` | 180 | Number of biochemistry timesteps per cell |

Biochemistry and physics parameters are hardcoded via `Params::default()` and `SpatialParams::default()`. See `simulations/calibration/parameter_provenance.md` for literature sources.

## Output format

Three output types are produced:

### 1. `spatial_summary.json` -- aggregate statistics per treatment

JSON array of 4 objects (Control, RSL3, SDT, PDT):

```json
{
  "treatment": "SDT",
  "total_tumor": 155724,
  "total_dead": 137500,
  "overall_death_rate": 0.883,
  "glycolytic": { "total": 62289, "dead": 54900 },
  "oxphos": { "total": 31144, "dead": 27800 },
  "persister": { "total": 31144, "dead": 28500 },
  "persister_nrf2": { "total": 31147, "dead": 26300 }
}
```

### 2. `depth_kill_curves.csv` -- death rate by depth for all treatments

```csv
depth_mm,Control,RSL3,SDT,PDT
0.00,0.003,0.128,0.998,0.997
0.50,0.002,0.131,0.985,0.912
1.00,0.001,0.125,0.971,0.234
...
```

### 3. `spatial_death_{treatment}.csv` -- 2D death heatmaps

Per-treatment CSV matrices (500x500) encoding cell state: 0 = stromal, 1 = dead tumor, 2 = alive tumor.

## Reproducing manuscript claims

**Chapter 6, Section 6.1, Figure 8 (depth-kill curves):**
```bash
cargo run --release -p sim-spatial -- --output-dir output/spatial
# Inspect: output/spatial/depth_kill_curves.csv
# Expected: PDT kill rate drops below 50% at ~1mm depth
# Expected: SDT kill rate remains >80% through ~5mm depth
# Expected: RSL3 is depth-independent (uniform distribution)
```

**Core claim: "PDT ~1mm, SDT ~10mm penetration":**
```bash
# Parse depth_kill_curves.csv
# PDT: effective kill depth (>50% death) ~1mm
# SDT: effective kill depth (>50% death) ~5-10mm
# RSL3: flat curve (no depth dependence)
```

## Caveats

- 2D grid (10mm x 10mm at default settings) -- real tumors are 3D
- O2 is uniform (no hypoxia gradient) -- see sim-tme for oxygen effects
- No in-vivo lipid remodeling (SCD1/MUFA) -- see sim-invivo for 3D context
- Acoustic and optical attenuation coefficients are literature estimates, not measured for a specific tumor type
- Iron diffusion is a nearest-neighbor simplification of the true extracellular diffusion field
- 180-step simulation represents a single treatment window, not repeated dosing
