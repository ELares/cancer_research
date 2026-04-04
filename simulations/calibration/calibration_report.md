# Calibration Report

Comparison of current simulation outputs against published calibration targets.

**Results: 5 PASS, 0 FAIL, 0 SKIP out of 5 targets.**

| Target | Status | Observed | Target | Residual | Tolerance | Confidence |
|--------|--------|----------|--------|----------|-----------|------------|
| **persister_rsl3_death_rate** | PASS | 0.4248 | 0.425 | -0.0002 | 0.1 | medium |
| **pdt_kill_at_9mm** | PASS | 0.0105 | 0.005 | +0.0055 | 0.015 | high |
| **rsl3_window_closure_72h** | PASS | 0.0143 | < 0.1 | -0.0857 | 0.05 | medium |
| **invivo_mufa_protection_factor** | PASS | 18.5516 | 18.6 | -0.0484 | 5.0 | medium |
| **sdt_penetration_1cm** | PASS | 0.8412 | > 0.84 | -0.0012 | 0.1 | high |

## Details

- [+] **persister_rsl3_death_rate**: Persister (FSP1-low) cell death rate under RSL3 in 2D culture
- [+] **pdt_kill_at_9mm**: PDT kill rate at 9mm tissue depth
- [+] **rsl3_window_closure_72h**: RSL3 death rate drops below 10% by 72 hours post-treatment
- [+] **invivo_mufa_protection_factor**: In-vivo MUFA provides ~18x protection against RSL3 for persisters
- [+] **sdt_penetration_1cm**: SDT maintains >84% kill rate through 1cm tissue depth

## Interpretation

- PASS means the current simulation output is within tolerance of the published target.
- FAIL means the output is outside tolerance and the parameter may need adjustment.
- SKIP means the required simulation output file was not found locally.
- See `parameter_provenance.md` for the source and confidence of each parameter.
- See `targets.yaml` for the full target definitions including source PMIDs.
