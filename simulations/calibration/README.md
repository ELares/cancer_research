# Simulation Calibration

This directory contains calibration infrastructure for the ferroptosis simulation suite.

## What's here

| File | Purpose |
|------|---------|
| `targets.yaml` | Calibration targets linking simulation observables to published data |
| `calibrate.py` | Evaluation script that compares simulation outputs to targets |
| `calibration_report.md` | Generated report (latest evaluation results) |
| `parameter_provenance.md` | Provenance document for all ~30 simulation parameters |

## Quick start

```bash
python simulations/calibration/calibrate.py --evaluate
```

This reads existing simulation output files (no recompilation needed) and reports how close each observable is to its calibration targets.

## Target types

Each target has a `target_type` field:

- **`calibration`**: The target value comes from independent experimental data. A PASS means the model reproduces a measured outcome.
- **`self-consistency`**: The target value follows from the model's own hard-coded physics assumptions (e.g., Beer-Lambert attenuation from `pdt_mu_eff`). A PASS only verifies the binary reproduces those assumptions correctly — useful as a regression check, but not independent validation.

## Staleness detection

The script checks whether simulation output files are older than the Rust source code. If so, results may reflect a previous build and a warning is printed. Re-run the relevant simulation to get current outputs.

## How it works

1. `targets.yaml` defines 5 calibration targets, each specifying:
   - Which simulation binary and output file to read
   - How to extract the observable (phenotype, treatment, depth, timepoint)
   - The target value and acceptable tolerance
   - The source publication (PMID where available)

2. `calibrate.py --evaluate` loads the targets, reads the existing JSON/CSV output files, extracts the relevant observable, and computes the residual.

3. The report shows PASS/FAIL/SKIP for each target, plus staleness warnings.

## Adding new targets

Add an entry to the `targets` list in `targets.yaml`:

```yaml
- id: new_target_name
  target_type: calibration  # or self-consistency
  description: "What this target measures"
  source_pmid: "12345678"
  binary: sim-original  # or sim-spatial, sim-window, sim-invivo
  output_file: "simulation_results.json"
  extraction:
    phenotype_contains: "Glycolytic"
    treatment: "SDT"
    field: "death_rate"
  target_value: 0.87
  tolerance: 0.05
  confidence: "high"
  parameters_affected:
    - sdt_ros
```

Then re-run `python calibrate.py --evaluate`.

## Future work

- **`--optimize` mode**: Automated parameter fitting using scipy.optimize. Requires adding `--params-file` CLI support to the Rust binaries so parameters can be passed at runtime without recompilation.
- **More targets**: Additional calibration points from GSH depletion kinetics, GPX4 expression time-courses, and cell-line-specific IC50 data.
- **Uncertainty quantification**: Propagating parameter confidence intervals through simulations to produce prediction intervals on manuscript claims.
