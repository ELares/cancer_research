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

This reads existing simulation output files (no recompilation needed) and reports how close each observable is to its published calibration target.

## How it works

1. `targets.yaml` defines 5 calibration targets, each specifying:
   - Which simulation binary and output file to read
   - How to extract the observable (phenotype, treatment, depth, timepoint)
   - The target value and acceptable tolerance
   - The source publication (PMID where available)

2. `calibrate.py --evaluate` loads the targets, reads the existing JSON/CSV output files from `simulations/output/`, extracts the relevant observable, and computes the residual.

3. The report shows PASS/FAIL/SKIP for each target.

## Adding new targets

Add an entry to the `targets` list in `targets.yaml`:

```yaml
- id: new_target_name
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
