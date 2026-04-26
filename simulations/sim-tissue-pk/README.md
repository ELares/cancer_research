# sim-tissue-pk

Tissue-specific drug penetration simulation computing ferroptosis kill rates as a function of distance from the nearest blood vessel across different tissue types, using a Krogh cylinder model with exponential concentration decay.

**Manuscript reference:** Chapter 8, Section 8.2

## What it models

1. **Krogh cylinder geometry** -- radial drug diffusion from a central blood vessel through tissue, with concentration decaying exponentially: C(r) = C0 * exp(-r / lambda)
2. **Tissue-specific vascular parameters** -- three tissue types with different inter-vessel distances and vascular permeability:
   - Epithelial (well-vascularized)
   - Epithelial (poorly-vascularized)
   - Neuroectodermal (CNS)
3. **Drug transport profiles** -- two drugs modeled: RSL3-like (small molecule GPX4 inhibitor) and doxorubicin (transport reference with published lambda)
4. **Distance-dependent ferroptosis** -- GPX4 inhibition scaled by local drug concentration at each radial bin; full ferroptosis biochemistry run per cell
5. **Kill depth analysis** -- effective kill depth defined as the distance where death rate drops below 10%

## Quick start

```bash
cd simulations
cargo build --release -p sim-tissue-pk
cargo run --release -p sim-tissue-pk
```

Runtime: ~10-30 seconds (50 bins x 1,000 cells x 6 tissue-drug combinations).

## Parameters

All parameters are hardcoded (no CLI arguments).

| Parameter | Value | Description |
|-----------|-------|-------------|
| N_RADIAL_BINS | 50 | Number of radial distance bins from vessel wall |
| N_CELLS_PER_BIN | 1,000 | Cells simulated at each distance |
| Seed | 42 | Random seed |
| Phenotype | Persister (FSP1-low) | Cell type used for all conditions |

Drug and tissue parameters are defined in `ferroptosis_core::drug_transport`:
- `rsl3_like()` -- RSL3 transport parameters
- `doxorubicin_transport_reference()` -- published transport reference (lambda ~40-80 um)
- `epithelial_well_vascularized()`, `epithelial_poorly_vascularized()`, `neuroectodermal_cns()` -- tissue types

## Output format

Output directory: `output/tissue-pk/`

### 1. `tissue_pk_results.csv` -- per-bin results

```csv
distance_um,concentration,death_rate,ci_low,ci_high,n_cells,n_dead,tissue,drug
0.0,1.000,0.425,0.395,0.455,1000,425,Epithelial (well-vasc),RSL3-like
20.0,0.914,0.398,0.368,0.428,1000,398,Epithelial (well-vasc),RSL3-like
...
```

| Field | Description |
|-------|-------------|
| distance_um | Distance from vessel wall in micrometers |
| concentration | Normalized drug concentration (0.0-1.0) |
| death_rate | Fraction of cells killed at this distance |
| tissue | Tissue type name |
| drug | Drug name |

### 2. `tissue_pk_summary.json` -- aggregate per tissue-drug pair

```json
{
  "tissue": "Epithelial (well-vasc)",
  "drug": "RSL3-like",
  "penetration_length_um": 224.0,
  "max_distance_um": 60.0,
  "vessel_wall_concentration": 1.0,
  "vessel_wall_death_rate": 0.425,
  "effective_kill_depth_um": 55.0,
  "overall_kill_fraction": 0.40
}
```

| Field | Description |
|-------|-------------|
| penetration_length_um | Exponential decay length (lambda) |
| max_distance_um | Maximum vessel-to-vessel half-distance for the tissue |
| effective_kill_depth_um | Distance where death rate drops below 10% |
| overall_kill_fraction | Average kill fraction across all bins |

## Reproducing manuscript claims

**Chapter 8, Section 8.2 (tissue-specific kill fractions):**
```bash
cargo run --release -p sim-tissue-pk
# stderr output includes "=== Summary ===" table
# Expected overall kill fractions (approximate):
#   Well-vascularized epithelial: ~40%
#   Poorly-vascularized epithelial: ~12%
#   CNS (neuroectodermal): ~2.6%
# Manuscript claim: 40% -> 12% -> 2.6% -> 1.8%
```

**Transport consistency check:**
```bash
cargo run --release -p sim-tissue-pk 2>&1 | grep "Doxorubicin-transport"
# Expected: lambda ~40-80 um (Minchinton 2006 published range)
# This validates the transport model parameterization
```

## Caveats

- All drugs use the RSL3/GPX4-inhibition pathway; differences reflect transport profiles only, not distinct pharmacology
- Krogh cylinder assumes uniform tissue between vessels -- real tumor vasculature is tortuous and heterogeneous
- Drug transport parameters are chosen to produce a lambda in the published range (self-consistency check, not independent calibration)
- No temporal dynamics -- steady-state concentration profile assumed
- No MUFA protection (2D culture params) -- see sim-invivo for in-vivo lipid remodeling effects
- The fourth tissue type mentioned in the manuscript claim (1.8%) may correspond to a different parameterization or tumor-PK composition; see sim-tumor-pk for the complementary PBPK model
