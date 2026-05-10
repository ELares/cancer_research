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
| `--photosensitizer` | `uniform=1` | Photosensitizer PK model. Spec format (case-insensitive): `uniform[=N]` (constant fraction; default `1`) or `porfimer[=t_half[,t_dist[,phi]]]` (single-exponential plasma kinetics with optional saturating distribution phase and relative singlet-O₂ yield). Defaults: `t_half=504` h (Bellnier 2006), `t_dist=0` h (legacy "light at peak"; set to ~24-48 h to model porfimer absorption rise), `phi=1.0` (porfimer-equivalent yield baseline). |
| `--dli-h` | 0.0 | Drug-light interval in hours: time to light delivery. With `t_dist=0` (default), this is post-distribution-peak; with `t_dist>0`, this can be the **clinical DLI from injection** (the model holds drug at peak for the first `t_dist` hours then begins decay). Validated as finite and ≥ 0. |

Biochemistry and physics parameters are hardcoded via `Params::default()` and `SpatialParams::default()`. See `simulations/calibration/parameter_provenance.md` for literature sources.

### Photosensitizer / DLI examples

```bash
# Default — Photosensitizer::Uniform(1.0), DLI=0. Byte-identical to pre-CLI behavior.
cargo run --release -p sim-spatial -- --output-dir output/spatial

# Porfimer at one terminal half-life past peak (DLI ≈ 21 d), legacy spec form.
# Drug is at ~50% of peak, so the PDT kill rate drops materially.
cargo run --release -p sim-spatial \
    --grid-size 100 \
    --photosensitizer porfimer --dli-h 504 \
    --output-dir output/spatial-porfimer-504h

# Phi sweep — half ROS yield via the third positional spec field.
# `porfimer=504,0,0.5` reads: t_half=504, t_dist=0, phi=0.5.
cargo run --release -p sim-spatial \
    --grid-size 100 \
    --photosensitizer "porfimer=504,0,0.5" \
    --output-dir output/spatial-porfimer-phi05

# Distribution-phase: clinical DLI from injection.
# Setting t_dist=36 h (Bellnier 2006 midpoint) lets --dli-h be the
# clinical drug-light interval; the model holds drug at peak for the
# first 36 h, then begins exponential decay. Here the user asks for
# a 24-h clinical DLI, which falls inside the distribution phase, so
# the simulated drug stays at peak (PDT identical to baseline).
cargo run --release -p sim-spatial \
    --grid-size 100 \
    --photosensitizer "porfimer=504,36" --dli-h 24 \
    --output-dir output/spatial-porfimer-clinical-dli24

# Drug-presence sweep across DLI (built-in --dli-sweep is a follow-up).
for dli in 0 24 168 504 1008; do
    cargo run --release -p sim-spatial \
        --grid-size 100 --seed 42 \
        --photosensitizer porfimer --dli-h "$dli" \
        --output-dir "output/sweep/dli-${dli}h"
done
```

Only PDT kill rate is affected by `--photosensitizer` / `--dli-h`. Control / RSL3 / SDT remain identical to the default invocation regardless of these flags.

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
