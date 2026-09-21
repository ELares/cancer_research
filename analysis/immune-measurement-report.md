# Passive 3D immune measurements

These are observations of the existing model, not independent biological validation.

| Measurement | Control | RSL3 | SDT |
|---|---:|---:|---:|
| Tumor cells | 82,519 | 82,519 | 82,519 |
| Ferroptotic death events | 6 | 964 | 77,710 |
| Completed releases | 4 | 937 | 77,710 |
| Deaths with release censored at horizon | 2 | 27 | 0 |
| Mean death-time LP, all ferroptotic deaths | 10.2098 | 10.2326 | 10.314 |
| Mean death-time LP, completed-release cohort | 10.1916 | 10.2342 | 10.314 |
| Mean release-time LP, completed-release cohort | 18.1499 | 17.9785 | 19.6915 |
| Mean horizon LP, censored-death cohort | 15.3262 | 13.8142 | undefined (n=0) |
| DAMP injected at completed releases | 72.5996 | 16845.9 | 1.53023e+06 |
| DAMP injected after last immune update | 30.6525 | 372.984 | 0 |
| DAMP field sum before terminal additions | 37.2575 | 3540.08 | 16812.6 |
| DAMP field sum after terminal additions | 67.91 | 3913.06 | 16812.6 |
| Unique eligible living tumor cells | 912 | 68,916 | 6,068 |
| Eligible cell-step opportunities | 13,763 | 5,727,847 | 577,423 |
| Immune kills | 0 | 28 | 111 |
| Kills / unique eligible cells | 0 | 0.000406292 | 0.0182927 |
| Kills / eligible cell-step opportunities | 0 | 4.8884e-06 | 0.000192233 |
| Mean local DAMP per eligible cell-step | 0.0370042 | 0.0496243 | 1.90024 |

## What this measures

The three canonical conditions share a 60³ grid, 20 µm cells, a 540 µm tumor radius, an oxygen gradient with λ=120 µm, 180 simulation steps, and the default 3D immune configuration. Geometry seed 42 is shared; treatment-specific runtime seeds are preserved. Stromal shielding, pH modulation, repopulation and additional killing modalities are off.

Eligibility is counted immediately before each immune update: a living tumor cell, step ≥60, and local DAMP ≥0.01. Unique cells count the union of those cells; cell-step opportunities count repeated eligibility. A cell killed in that update is included. Both rates are descriptive conditional summaries, not causal treatment effects or a correction for competing ferroptotic deaths. Different cells and exposure durations contribute across arms.

LP is in model units. DAMP is an uncalibrated model field with one injected unit per LP unit. Completed release occurs five steps after ferroptotic death, before diffusion and immune killing. A release scheduled at or beyond step 180 is right-censored; its horizon LP and terminal DAMP addition are reported separately. Terminal additions occur after every immune update and cannot explain kills already counted. Field sums differ from cumulative injections because of clearance.

## Interpretation and limits

The observed 3D SDT:RSL3 total immune-kill ratio is 3.96429. This does not measure DAMP potency per dead cell. The completed-release means compare a named event cohort; they exclude horizon-censored releases and should not be substituted for all-death means.

This is one existing geometry seed and one existing runtime seed per condition. There are no replicate uncertainty intervals, no calibrated conversion from cell-steps to hours, and no fitted per-cell experimental release threshold. The historical 104:1 figure belongs to the 2D model; this 3D observation neither reconstructs that experiment nor validates its mechanism. The older `immune_kills / max(total_tumor - ferroptosis_kills, 1)` statistic uses a final non-ferroptotic pool, not the observed eligible population.

All death, release, censoring, eligibility and kill counts reconcile from the archived records. Every measured legacy result equals its corresponding default-matrix row; the complete 24-condition summary retains SHA-256 `973d9443bd05851721d53d9a15e53dc828517bba4bf307eb1c4f0b1f60d27289`.

## Reproduction

```bash
python3 scripts/immune_measurement_report.py
```

This validates compressed raw event records, per-cell eligibility counts, every step, artifact hashes and frozen source hashes, then regenerates this report without running a simulation. See the [protocol](../docs/IMMUNE_MEASUREMENT_PROTOCOL.md) for capture and validation commands.

- Frozen source commit: `c08e16eecb11368e53d856f5d2e96fe885b3f092`.
- Captured: `2026-09-21T03:07:24.924185+00:00`.
- Platform: `macOS-26.5.1-arm64-arm-64bit-Mach-O`; `rustc 1.96.0 (ac68faa20 2026-05-25)`; Rayon threads: 8.
- Built binary SHA-256: `57be1a445388e01d910f497421efc2a2ad62d116e1b5e3f345e62dfcb8a1a52e`.
- [Archive and provenance](immune-measurements/manifest.json).
- [Unchanged production summary](immune-measurements/baseline-summary.json).

- `immune_Control` runtime seed: `10143494520844665724`.
- `immune_RSL3` runtime seed: `1491324243477256001`.
- `immune_SDT` runtime seed: `9942688334898329752`.
