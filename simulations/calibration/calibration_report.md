# Calibration Report

Comparison of current simulation outputs against published calibration targets.

**Results: 1 PASS, 4 STALE out of 5 targets.**

| Target | Status | Observed | Target | Residual | Tolerance | Confidence |
|--------|--------|----------|--------|----------|-----------|------------|
| **persister_rsl3_death_rate** | STALE | 0.4248 | 0.425 | -0.0002 | 0.1 | medium |
| **pdt_kill_at_9mm** | STALE | 0.0105 | 0.005 | +0.0055 | 0.015 | high |
| **rsl3_window_closure_72h** | STALE | 0.0143 | < 0.1 | -0.0857 | 0.05 | medium |
| **invivo_mufa_protection_factor** | PASS | 18.5516 | 18.6 | -0.0484 | 5.0 | low |
| **sdt_penetration_1cm** | STALE | 0.8412 | > 0.84 | -0.0012 | 0.1 | high |

## Details

- [~] **persister_rsl3_death_rate**: Persister (FSP1-low) cell death rate under RSL3 in 2D culture
  - Value within tolerance but output is stale: Output simulation_results.json is older than source code — results may be stale. Re-run the simulation.
- [~] **pdt_kill_at_9mm**: PDT kill rate at 9mm tissue depth
  - Value within tolerance but output is stale: Output depth_kill_curves.csv is older than source code — results may be stale. Re-run the simulation.
- [~] **rsl3_window_closure_72h**: RSL3 death rate drops below 10% by 72 hours post-treatment
  - Value within tolerance but output is stale: Output vulnerability_window.csv is older than source code — results may be stale. Re-run the simulation.
- [+] **invivo_mufa_protection_factor**: In-vivo MUFA provides substantial protection against RSL3 for persisters
- [~] **sdt_penetration_1cm**: SDT maintains >84% kill rate through 1cm tissue depth
  - Value within tolerance but output is stale: Output depth_kill_curves.csv is older than source code — results may be stale. Re-run the simulation.

## Interpretation

- PASS means the simulation output is within tolerance and the output file is current.
- STALE means the value is within tolerance but the output file is older than the source code — re-run the simulation to confirm.
- FAIL means the output is outside tolerance and the parameter may need adjustment.
- SKIP means the required simulation output file was not found locally.
- Targets with `target_type: self-consistency` verify that the binary reproduces its own hard-coded physics assumptions. They are useful as regression checks but do not constitute independent experimental validation.
- See `parameter_provenance.md` for the source and confidence of each parameter.
- See `targets.yaml` for the full target definitions including source PMIDs.
