# Controlled 2D DAMP source/recipient comparison

A deterministic experiment on the existing transport and activation rules. These artificial source and recipient masks are not biological validation or an attribution of the historical SDT–RSL3 immune-kill contrast.

The primary response sums DAMP/(DAMP+50) at eligible, available recipient sites during steps 60–179 and divides by the same 3,840 reference sites in every condition. Units are activation steps per reference recipient. Zero opportunities give zero exposure; conditional means remain undefined. No stochastic kills, source feedback or inference intervals are used.

## Every planned condition

| Condition | DAMP/source | Reference exposure | Opportunities | Conditional activation | Max eligible DAMP | ≥50 | ≥450 |
|---|---:|---:|---:|---:|---:|---:|---:|
| n64_m1280_t30_r3840 | 20 | 0.06187539 | 286948 | 0.0008280298 | 0.2175792 | 0 | 0 |
| n64_m1280_t30_r960 | 20 | 0.01546885 | 71737 | 0.0008280298 | 0.2175792 | 0 | 0 |
| n64_m1280_t60_r3840 | 20 | 0.1697066 | 390340 | 0.001669502 | 1.565658 | 0 | 0 |
| n64_m1280_t60_r960 | 20 | 0.04242666 | 97585 | 0.001669502 | 1.565658 | 0 | 0 |
| n64_m5120_t30_r3840 | 80 | 0.2650817 | 445804 | 0.002283321 | 0.8703185 | 0 | 0 |
| n64_m5120_t30_r960 | 80 | 0.06627043 | 111451 | 0.002283321 | 0.8703185 | 0 | 0 |
| n64_m5120_t60_r3840 | 80 | 0.6767939 | 451464 | 0.00575658 | 6.26263 | 0 | 0 |
| n64_m5120_t60_r960 | 80 | 0.1691985 | 112866 | 0.00575658 | 6.26263 | 0 | 0 |
| n256_m1280_t30_r3840 | 5 | 0.0646918 | 290884 | 0.0008540054 | 0.1231309 | 0 | 0 |
| n256_m1280_t30_r960 | 5 | 0.01617295 | 72721 | 0.0008540054 | 0.1231309 | 0 | 0 |
| n256_m1280_t60_r3840 | 5 | 0.1764902 | 404292 | 0.001676319 | 0.49591 | 0 | 0 |
| n256_m1280_t60_r960 | 5 | 0.04412254 | 101073 | 0.001676319 | 0.49591 | 0 | 0 |
| n256_m5120_t30_r3840 | 20 | 0.2761855 | 447876 | 0.002367959 | 0.4925235 | 0 | 0 |
| n256_m5120_t30_r960 | 20 | 0.06904637 | 111969 | 0.002367959 | 0.4925235 | 0 | 0 |
| n256_m5120_t60_r3840 | 20 | 0.7078078 | 458984 | 0.005921736 | 1.98364 | 0 | 0 |
| n256_m5120_t60_r960 | 20 | 0.176952 | 114746 | 0.005921736 | 1.98364 | 0 | 0 |
| zero_r3840 | 0 | 0 | 0 | undefined | undefined | 0 | 0 |
| zero_r960 | 0 | 0 | 0 | undefined | undefined | 0 | 0 |

Condition IDs name source count (n), total injection (m), release step (t), and available recipients (r). Both zero controls use a zero injection at step 60.

## All prespecified contrasts

Differences below are left minus right on the primary response. Equal-total-mass comparisons change source positions/count and the amount per source together. Equal-per-source-amount comparisons change source count and total injection together. Neither is a pure density effect holding both quantities fixed. Availability comparisons use an unchanged field and a fixed subset.

| Comparison | Left | Right | Difference |
|---|---|---|---:|
| equal_total_mass | n256_m1280_t30_r3840 | n64_m1280_t30_r3840 | 0.002816411 |
| equal_total_mass | n256_m1280_t30_r960 | n64_m1280_t30_r960 | 0.0007041027 |
| equal_total_mass | n256_m1280_t60_r3840 | n64_m1280_t60_r3840 | 0.006783506 |
| equal_total_mass | n256_m1280_t60_r960 | n64_m1280_t60_r960 | 0.001695877 |
| equal_total_mass | n256_m5120_t30_r3840 | n64_m5120_t30_r3840 | 0.01110374 |
| equal_total_mass | n256_m5120_t30_r960 | n64_m5120_t30_r960 | 0.002775935 |
| equal_total_mass | n256_m5120_t60_r3840 | n64_m5120_t60_r3840 | 0.03101391 |
| equal_total_mass | n256_m5120_t60_r960 | n64_m5120_t60_r960 | 0.007753478 |
| release_amount | n64_m5120_t30_r3840 | n64_m1280_t30_r3840 | 0.2032063 |
| release_amount | n64_m5120_t30_r960 | n64_m1280_t30_r960 | 0.05080159 |
| release_amount | n64_m5120_t60_r3840 | n64_m1280_t60_r3840 | 0.5070873 |
| release_amount | n64_m5120_t60_r960 | n64_m1280_t60_r960 | 0.1267718 |
| release_amount | n256_m5120_t30_r3840 | n256_m1280_t30_r3840 | 0.2114937 |
| release_amount | n256_m5120_t30_r960 | n256_m1280_t30_r960 | 0.05287342 |
| release_amount | n256_m5120_t60_r3840 | n256_m1280_t60_r3840 | 0.5313177 |
| release_amount | n256_m5120_t60_r960 | n256_m1280_t60_r960 | 0.1328294 |
| release_timing | n64_m1280_t60_r3840 | n64_m1280_t30_r3840 | 0.1078313 |
| release_timing | n64_m1280_t60_r960 | n64_m1280_t30_r960 | 0.02695782 |
| release_timing | n64_m5120_t60_r3840 | n64_m5120_t30_r3840 | 0.4117122 |
| release_timing | n64_m5120_t60_r960 | n64_m5120_t30_r960 | 0.1029281 |
| release_timing | n256_m1280_t60_r3840 | n256_m1280_t30_r3840 | 0.1117984 |
| release_timing | n256_m1280_t60_r960 | n256_m1280_t30_r960 | 0.02794959 |
| release_timing | n256_m5120_t60_r3840 | n256_m5120_t30_r3840 | 0.4316224 |
| release_timing | n256_m5120_t60_r960 | n256_m5120_t30_r960 | 0.1079056 |
| recipient_availability | n64_m1280_t30_r960 | n64_m1280_t30_r3840 | -0.04640654 |
| recipient_availability | n64_m1280_t60_r960 | n64_m1280_t60_r3840 | -0.12728 |
| recipient_availability | n64_m5120_t30_r960 | n64_m5120_t30_r3840 | -0.1988113 |
| recipient_availability | n64_m5120_t60_r960 | n64_m5120_t60_r3840 | -0.5075955 |
| recipient_availability | n256_m1280_t30_r960 | n256_m1280_t30_r3840 | -0.04851885 |
| recipient_availability | n256_m1280_t60_r960 | n256_m1280_t60_r3840 | -0.1323676 |
| recipient_availability | n256_m5120_t30_r960 | n256_m5120_t30_r3840 | -0.2071391 |
| recipient_availability | n256_m5120_t60_r960 | n256_m5120_t60_r3840 | -0.5308559 |
| equal_per_source_amount | n256_m5120_t30_r3840 | n64_m1280_t30_r3840 | 0.2143101 |
| equal_per_source_amount | n256_m5120_t30_r960 | n64_m1280_t30_r960 | 0.05357752 |
| equal_per_source_amount | n256_m5120_t60_r3840 | n64_m1280_t60_r3840 | 0.5381012 |
| equal_per_source_amount | n256_m5120_t60_r960 | n64_m1280_t60_r960 | 0.1345253 |

## Accounting and limits

All planned conditions and contrasts are published. Release occurs before transport and clearance, observation after both. Whole-field mass follows M[t]=0.97×(M[t−1]+injection[t]); the rectangular boundary loses no mass to unrepresented neighbors. Recipient availability does not remove DAMP or alter transport. Cell and step ledgers reconcile; selected recipients retain identical observations in both availability arms.

The archived observations contain per-step and per-recipient aggregates, not every field value at every step. Offline reconstruction checks those ledgers and analytic invariants; it cannot independently recompute transport or recover each activation evaluation. FNV fingerprints are noncryptographic consistency checks; SHA-256 identifies the archived payloads. Replay the frozen Rust driver for full transport reproduction. Neither check establishes biological accuracy.

DAMP is in model field units. Model steps have no new conversion to hours. The 0.001 diffusion cutoff and 0.01 eligibility floor remain active; general linear-superposition claims do not follow. The masks and timing are artificial, and the result is conditional on this finite design. It supplies no per-death immunogenicity estimate, calibrated assay, clinical conclusion or P5 replacement.

## Reproduction

```bash
python3 scripts/immune_2d_controlled_report.py
```

The default command validates the archive and regenerates both reports without running Rust. `--render-only` also validates and reconstructs the archive, then writes only Markdown. Report files are written sequentially, not as an atomic pair. A new capture requires a clean committed checkout, an absent archive destination, and unchanged historical production and observer outputs.

- [Frozen protocol](../docs/IMMUNE_2D_CONTROLLED_PROTOCOL.md).
- [Plan](../scripts/immune_2d_controlled_plan.json).
- [Archive manifest](immune-2d-controlled/manifest.json).
- [Complete derived data](immune-2d-controlled.json).
- Source freeze: `421603d2a347b8db6087287ed8cba426ba5dfe1c`.
- Binary SHA-256: `70c0c4817f7cf11ddab0328b75f21a10402182dae86414d2b281e827316817d4`.
