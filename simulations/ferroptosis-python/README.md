# ferroptosis-core

Mechanistic ferroptosis cell death simulation engine for cancer research.

Models the ferroptosis pathway from ROS generation through GSH depletion, GPX4/FSP1 repair, lipid peroxidation, and cell death. Built in Rust for speed, exposed to Python via PyO3. Simulates 10,000 cells in ~50ms.

## Install

```bash
pip install ferroptosis-core
```

No Rust toolchain required — pre-built wheels are provided for Linux, macOS, and Windows.

## Quick start

```python
import ferroptosis_core as fc

# Single cell
result = fc.sim_cell("Persister", "RSL3", seed=42)
print(result)  # {'dead': True, 'lp': 10.56, 'gsh': 0.02, 'gpx4': 0.40}

# Population (10,000 cells, parallel)
stats = fc.sim_batch("Persister", "RSL3", n=10000, seed=42)
print(f"Death rate: {stats['death_rate']:.1%}")  # ~40%
```

## API

### `default_params() -> dict`

Returns the 20 default biochemistry rate constants (2D culture context).

### `invivo_params() -> dict`

Returns parameters with SCD1/MUFA lipid remodeling enabled (in-vivo context).

### `sim_cell(phenotype, treatment, seed, context="2d", **kwargs) -> dict`

Simulate a single cell through 180 timesteps.

**Returns:** `{'dead': bool, 'lp': float, 'gsh': float, 'gpx4': float}`

### `sim_batch(phenotype, treatment, n, seed, context="2d", **kwargs) -> dict`

Simulate `n` cells in parallel and return population statistics.

**Returns:** `{'death_rate': float, 'ci_low': float, 'ci_high': float, 'n_dead': int, 'n_cells': int, 'mean_lp': float, 'mean_gsh': float, 'mean_gpx4': float}`

## Phenotypes

`"Glycolytic"`, `"OXPHOS"`, `"Persister"`, `"PersisterNrf2"`, `"Stromal"`

## Treatments

`"Control"`, `"RSL3"`, `"SDT"`, `"PDT"`

## Contexts

- `"2d"` (default) — standard 2D culture parameters
- `"invivo"` — SCD1/MUFA lipid remodeling enabled

## Parameter overrides

Any of the 20 biochemistry parameters can be overridden via keyword arguments:

```python
# Sweep GPX4 inhibition strength
for inhib in [0.0, 0.25, 0.5, 0.75, 1.0]:
    stats = fc.sim_batch("Persister", "RSL3", n=1000, seed=42,
                         rsl3_gpx4_inhib=inhib)
    print(f"inhib={inhib:.2f}  death_rate={stats['death_rate']:.1%}")
```

Use `fc.default_params()` to see all available parameters and their default values.

## Performance

`sim_batch` uses Rayon for parallel execution across all available CPU cores. The GIL is released during simulation, so Python threads are not blocked.

## Part of the Cancer Research Synthesis project

This package is the simulation engine from an open cancer research repository that combines cross-literature analysis of 4,830 research articles with Monte Carlo biochemical simulations.

Full repository: https://github.com/ELares/cancer_research

## License

MIT
