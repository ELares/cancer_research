# Calibration Report

Comparison of current simulation outputs against published calibration targets.

**Results: 3 PASS, 5 STALE out of 8 targets.**

| Target | Status | Observed | Target | Residual | Tolerance | Confidence |
|--------|--------|----------|--------|----------|-----------|------------|
| **persister_rsl3_death_rate** | STALE | 0.4248 | 0.425 | -0.0002 | 0.1 | medium |
| **pdt_kill_at_9mm** | STALE | 0.0105 | 0.005 | +0.0055 | 0.015 | high |
| **rsl3_window_closure_72h** | STALE | 0.0143 | < 0.1 | -0.0857 | 0.05 | medium |
| **invivo_mufa_protection_factor** | STALE | 18.5516 | 18.6 | -0.0484 | 5.0 | low |
| **sdt_penetration_1cm** | STALE | 0.8412 | > 0.84 | -0.0012 | 0.1 | high |
| **3d_rsl3_o2_collapse_ratio** | PASS | 0.0160 | < 0.5 | -0.4840 | 0.5 | medium |
| **3d_immune_sdt_dominates** | PASS | 4.0345 | > 3.0 | -1.0345 | 1.0 | medium |
| **3d_stromal_boundary_shielding** | PASS | 0.4848 | < 0.7 | -0.2152 | 0.3 | medium |

## Details

- [~] **persister_rsl3_death_rate**: Persister (FSP1-low) cell death rate under RSL3 in 2D culture
  - Value within tolerance but output is stale: Output simulation_results.json is older than source code — results may be stale. Re-run the simulation.
- [~] **pdt_kill_at_9mm**: PDT kill rate at 9mm tissue depth
  - Value within tolerance but output is stale: Output depth_kill_curves.csv is older than source code — results may be stale. Re-run the simulation.
- [~] **rsl3_window_closure_72h**: RSL3 death rate drops below 10% by 72 hours post-treatment
  - Value within tolerance but output is stale: Output vulnerability_window.csv is older than source code — results may be stale. Re-run the simulation.
- [~] **invivo_mufa_protection_factor**: In-vivo MUFA provides substantial protection against RSL3 for persisters
  - Value within tolerance but output is stale: Output invivo_comparison.json is older than source code — results may be stale. Re-run the simulation.
- [~] **sdt_penetration_1cm**: SDT maintains >84% kill rate through 1cm tissue depth
  - Value within tolerance but output is stale: Output depth_kill_curves.csv is older than source code — results may be stale. Re-run the simulation.
- [+] **3d_rsl3_o2_collapse_ratio**: RSL3 at λ=120 µm: hypoxic-zone kill rate should be ≤ 50% of normoxic-zone kill rate (Q1 finding)
- [+] **3d_immune_sdt_dominates**: At λ=120 immune-on: SDT immune-kills should be > 3× RSL3 immune-kills (Q2 finding — the manuscript-keystone 104:1 → 4:1 reduction in 3D)
- [+] **3d_stromal_boundary_shielding**: At λ=120 immune-on RSL3: stromal-on boundary kill rate should be ≤ 70% of no-stromal boundary kill rate (Q3 finding — CAF shielding ≥ 30%)

## Interpretation

- PASS means the simulation output is within tolerance and the output file is current.
- STALE means the value is within tolerance but the output file is older than the source code — re-run the simulation to confirm.
- FAIL means the output is outside tolerance and the parameter may need adjustment.
- SKIP means the required simulation output file was not found locally.
- Targets with `target_type: self-consistency` verify that the binary reproduces its own hard-coded physics assumptions. They are useful as regression checks but do not constitute independent experimental validation.
- See `parameter_provenance.md` for the source and confidence of each parameter.
- See `targets.yaml` for the full target definitions including source PMIDs.
