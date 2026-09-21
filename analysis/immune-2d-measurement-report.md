# Passive 2D immune measurements

These observations reconstruct the historical model comparison; they are not independent biological validation.

| Measurement | Control | RSL3 | SDT |
|---|---:|---:|---:|
| Tumor cells | 159,033 | 159,033 | 159,033 |
| Ferroptotic deaths | 1 | 163 | 139,640 |
| Completed releases | 1 | 160 | 139,640 |
| Horizon-censored releases | 0 | 3 | 0 |
| Mean death LP, all deaths | 10.1823 | 10.2474 | 10.2657 |
| Mean death LP, completed cohort | 10.1823 | 10.2504 | 10.2657 |
| Mean release LP, completed cohort | 17.4509 | 18.2658 | 19.123 |
| Mean horizon LP, censored cohort | undefined (n=0) | 15.7911 | undefined (n=0) |
| DAMP injected at completed releases | 17.4509 | 2922.53 | 2.67033e+06 |
| DAMP injected after last immune update | 0 | 47.3734 | 0 |
| Field DAMP before terminal additions | 15.4491 | 622.294 | 33457.9 |
| Field DAMP after terminal additions | 15.4491 | 669.668 | 33457.9 |
| Peak field DAMP after terminal additions | 1.28462 | 19.7513 | 0.227207 |
| Unique eligible cells | 34 | 14,521 | 23,510 |
| Eligible cell-step opportunities | 78 | 770,701 | 2,314,460 |
| Immune kills | 0 | 5 | 521 |
| Kills / unique eligible cells | 0 | 0.000344329 | 0.0221608 |
| Kills / eligible cell-step opportunities | 0 | 6.4876e-06 | 0.000225107 |
| Mean local DAMP per eligible opportunity | 0.411385 | 0.06584 | 2.06872 |
| Maximum local DAMP among eligible opportunities | 1.3661 | 2.58529 | 10.5005 |
| Mean activation per eligible opportunity | 0.00807886 | 0.00130919 | 0.0385056 |
| Eligible opportunities with DAMP ≥ Kd | 0 | 0 | 0 |
| Fraction of eligible opportunities with DAMP ≥ Kd | 0 | 0 | 0 |
| Eligible opportunities with DAMP ≥ 9Kd | 0 | 0 | 0 |
| Fraction of eligible opportunities with DAMP ≥ 9Kd | 0 | 0 | 0 |

## What was observed

The historical baseline uses a 500×500 grid, 20 µm cells, a 4,500 µm tumor radius, geometry seed 42, 180 steps, λ=120 µm oxygen gradient, and the existing O2-independent SDT setting. Immune coupling uses the 2D defaults; anti-PD-1, stromal shielding and pH modulation are off. Treatment-specific runtime seeds and all parameters are archived.

Living tumor cells are eligible immediately before immune killing at step ≥60 and local DAMP ≥0.01. Cells killed by that update are included. Repeated eligibility forms cell-step opportunities; unique eligibility counts each cell once. Both kill rates are descriptive, not causal treatment effects or competing-risk corrections.

Activation is DAMP / (DAMP + Kd), with Kd=50 model units. DAMP ≥Kd corresponds to activation ≥0.5; DAMP ≥9Kd corresponds to activation ≥0.9. The fractions use eligible cell-step opportunities, not all cells or terminal heatmap pixels. Empty means, maxima and fractions are undefined. Aggregates bound and reconcile these observations; they do not reconstruct every unarchived cell-step exposure.

Completed DAMP releases occur five steps after death, before diffusion and immune killing. Releases scheduled at step 180 or later are censored; their horizon LP and terminal additions are separate. Terminal additions occur after all immune updates and cannot explain earlier immune kills. Field totals differ from cumulative injections because of clearance. LP and DAMP remain uncalibrated model quantities.

## Interpretation and limits

The observed SDT:RSL3 total immune-kill ratio is 104.2:1. This reconstructs the counts behind the historical approximately 104:1 comparison. It does not measure DAMP potency, DC maturation or immunogenicity per dead cell. Release means refer to the completed-release cohort; deaths censored at the horizon have separate means and denominators.

This is one geometry seed and one runtime seed per arm. No replicate uncertainty or experimental calibration is added. The eligible-exposure measurements allow the saturation explanation to be inspected, without establishing it as a causal explanation. The older final non-ferroptotic-pool normalized statistic is a different denominator. The frozen 3D report remains a separate model comparison with different geometry, RNG streams and transport settings.

## Reproduction

```bash
python3 scripts/immune_2d_measurement_report.py
```

Reconstruction checks event accounting, eligibility and activation summaries, archive hashes, frozen configuration, historical rows and the unchanged full 33-condition production summary without building or running a simulation.

- [Protocol](../docs/IMMUNE_2D_MEASUREMENT_PROTOCOL.md).
- [Archive and provenance](immune-2d-measurements/manifest.json).
- [Unchanged full baseline](immune-2d-measurements/baseline-summary.json).
- [Separate frozen 3D report](immune-measurement-report.md).
- Source commit: `5872fd3e89b509b9e28f9e309cd058361147a098`.
- Captured: `2026-09-21T04:38:41.534977+00:00`.
- Platform: `macOS-26.5.1-arm64-arm-64bit-Mach-O`; `rustc 1.96.0 (ac68faa20 2026-05-25)`; execution: serial.
- Built binary SHA-256: `710f7ace1f17851cabccecdf40aa7e14493362faa2e95b0e86bb0270d7a450d9`.
- Full baseline SHA-256: `e04f9699ddca9b84145cdecbcd8e21abb0bc760dc3d909908358d4554a5eee98`.

- Uninstrumented baseline source commit: `e9ec92c022ede1b00258789fecdcb9a46d9481d8`.
- Uninstrumented baseline binary SHA-256: `8acc5804528c436f60f0a4b5ce8923ed4807bd6869249d6fd61ab7c1ffa00f35`.
- Uninstrumented baseline toolchain: `rustc 1.96.0 (ac68faa20 2026-05-25); aarch64-apple-darwin; serial execution.`.

- `immune_Control` runtime seed: `42`.
- `immune_RSL3` runtime seed: `10000042`.
- `immune_SDT` runtime seed: `20000042`.
