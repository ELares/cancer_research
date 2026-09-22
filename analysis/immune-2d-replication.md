# Twenty-block 2D immune activation replication

A prospectively frozen simulation study; these results are not independent biological validation.

The mean paired SDT − RSL3 activation contrast is **0.555834** accumulated activation steps per initial tumor cell. The 95% whole-block percentile bootstrap interval is **[0.554223, 0.557439]** (positive under the frozen design).

Activation is DAMP/(DAMP+50), summed over eligible cell-step opportunities and divided by the initial tumor count. It is a composite exposure measure, not a probability, immunogenicity per death, or an assay-calibrated biological quantity. Zero opportunities give zero exposure.

## All blocks

| Block | Root seed | Control exposure | RSL3 exposure | SDT exposure | SDT − RSL3 |
|---:|---:|---:|---:|---:|---:|
| 1 | 4294967338 | 0 | 0.0059305 | 0.567148 | 0.561218 |
| 2 | 8589934634 | 0 | 0.00601159 | 0.566393 | 0.560381 |
| 3 | 12884901930 | 4.37952e-05 | 0.0051101 | 0.557375 | 0.552264 |
| 4 | 17179869226 | 1.22606e-05 | 0.00536302 | 0.561108 | 0.555745 |
| 5 | 21474836522 | 0 | 0.00531931 | 0.559413 | 0.554094 |
| 6 | 25769803818 | 3.84558e-05 | 0.00515792 | 0.559092 | 0.553934 |
| 7 | 30064771114 | 1.37251e-05 | 0.00537155 | 0.565619 | 0.560247 |
| 8 | 34359738410 | 0 | 0.00588817 | 0.560794 | 0.554906 |
| 9 | 38654705706 | 0 | 0.0048967 | 0.562083 | 0.557186 |
| 10 | 42949673002 | 0 | 0.00556259 | 0.556492 | 0.55093 |
| 11 | 47244640298 | 4.38796e-05 | 0.00501735 | 0.56757 | 0.562553 |
| 12 | 51539607594 | 1.12598e-05 | 0.00750893 | 0.560188 | 0.552679 |
| 13 | 55834574890 | 0 | 0.00668419 | 0.562612 | 0.555928 |
| 14 | 60129542186 | 0 | 0.00666101 | 0.566099 | 0.559438 |
| 15 | 64424509482 | 2.92539e-05 | 0.00648949 | 0.554578 | 0.548089 |
| 16 | 68719476778 | 0 | 0.00551388 | 0.559038 | 0.553524 |
| 17 | 73014444074 | 0 | 0.00591733 | 0.560933 | 0.555016 |
| 18 | 77309411370 | 0 | 0.0057218 | 0.55907 | 0.553348 |
| 19 | 81604378666 | 0 | 0.00606219 | 0.563128 | 0.557065 |
| 20 | 85899345962 | 1.43803e-05 | 0.00587313 | 0.564001 | 0.558128 |

## Distributions across blocks

Each cell gives median [minimum, maximum]; n is the number of defined blocks. The JSON names every undefined block and retains every arm's full event-accounting summary. Empty conditional means, maxima and zero-denominator ratios are undefined, not zero.

| Measurement | Control | RSL3 | SDT |
|---|---|---|---|
| activation sum | 0 [0, 6.9783]; n=20 | 921.988 [778.738, 1194.17]; n=20 | 89220.8 [88196.2, 90262.4]; n=20 |
| activation sum per initial tumor cell | 0 [0, 4.38796e-05]; n=20 | 0.00579746 [0.0048967, 0.00750893]; n=20 | 0.561021 [0.554578, 0.56757]; n=20 |
| censored deaths | 0 [0, 1]; n=20 | 3 [0, 6]; n=20 | 0 [0, 0]; n=20 |
| completed injected damp | 0 [0, 56.6583]; n=20 | 2741.64 [2354.83, 3443.12]; n=20 | 2.66907e+06 [2.66225e+06, 2.67472e+06]; n=20 |
| completed releases | 0 [0, 3]; n=20 | 152 [130, 190]; n=20 | 139616 [139,301, 139,895]; n=20 |
| damp after terminal | 10.2606 [0, 38.3004]; n=20 | 650.67 [547.525, 791.146]; n=20 | 33441.9 [33358.6, 33517.2]; n=20 |
| damp before terminal | 0 [0, 38.3004]; n=20 | 623.919 [497.478, 744.153]; n=20 | 33441.9 [33358.6, 33517.2]; n=20 |
| damp ge 9kd opportunities | 0 [0, 0]; n=20 | 0 [0, 0]; n=20 | 0 [0, 0]; n=20 |
| damp ge 9kd opportunities fraction | 0 [0, 0]; n=8 | 0 [0, 0]; n=20 | 0 [0, 0]; n=20 |
| damp ge kd opportunities | 0 [0, 0]; n=20 | 0 [0, 0]; n=20 | 0 [0, 0]; n=20 |
| damp ge kd opportunities fraction | 0 [0, 0]; n=8 | 0 [0, 0]; n=20 | 0 [0, 0]; n=20 |
| death lp all deaths | 10.3199 [10.2075, 10.6139]; n=12 | 10.2299 [10.2037, 10.2454]; n=20 | 10.2666 [10.2656, 10.2678]; n=20 |
| death lp completed cohort | 10.3199 [10.2202, 10.6139]; n=8 | 10.2318 [10.2058, 10.2472]; n=20 | 10.2666 [10.2656, 10.2678]; n=20 |
| eligible cell steps | 0 [0, 4,040]; n=20 | 686530 [593,085, 922,795]; n=20 | 2.31595e+06 [2,283,637, 2,352,610]; n=20 |
| ferroptotic deaths | 1 [0, 3]; n=20 | 153.5 [133, 196]; n=20 | 139616 [139,301, 139,895]; n=20 |
| final damp peak | 0.219557 [0, 19.3231]; n=20 | 16.9902 [2.91804, 21.6969]; n=20 | 0.227911 [0.224025, 0.232447]; n=20 |
| horizon lp censored cohort | 15.9823 [10.3533, 19.3231]; n=5 | 14.5438 [11.3651, 21.2613]; n=18 | undefined [undefined, undefined]; n=0 |
| immune kills | 0 [0, 0]; n=20 | 5 [1, 8]; n=20 | 529 [485, 559]; n=20 |
| kills per eligible cell step | 0 [0, 0]; n=8 | 7.01506e-06 [1.53728e-06, 1.16206e-05]; n=20 | 0.000230737 [0.000209842, 0.00024346]; n=20 |
| kills per unique eligible cell | 0 [0, 0]; n=8 | 0.000368822 [7.64584e-05, 0.000576606]; n=20 | 0.0226972 [0.0206902, 0.023976]; n=20 |
| max eligible damp | 1.4542 [1.33591, 1.72379]; n=8 | 2.54785 [1.90124, 2.89772]; n=20 | 10.4032 [9.93068, 11.3494]; n=20 |
| mean activation | 0.00299976 [0.00115157, 0.00441139]; n=8 | 0.00132592 [0.00127635, 0.00139885]; n=20 | 0.0385403 [0.0383625, 0.0386943]; n=20 |
| mean damp per eligible cell step | 0.151362 [0.0578916, 0.223236]; n=8 | 0.0666863 [0.0641872, 0.0703614]; n=20 | 2.0709 [2.06064, 2.07968]; n=20 |
| release lp completed cohort | 18.5763 [17.0652, 19.732]; n=8 | 18.0444 [17.9531, 18.1991]; n=20 | 19.1186 [19.1067, 19.1322]; n=20 |
| terminal injected damp | 0 [0, 19.3231]; n=20 | 46.799 [0, 102.214]; n=20 | 0 [0, 0]; n=20 |
| total tumor | 159033 [159,033, 159,033]; n=20 | 159033 [159,033, 159,033]; n=20 | 159033 [159,033, 159,033]; n=20 |
| unique eligible cells | 0 [0, 226]; n=20 | 13177.5 [11,891, 16,985]; n=20 | 23465 [23,255, 23,863]; n=20 |

## Threshold observations

Counts refer to blocks containing any eligible opportunity at each activation threshold.

| Arm | Activation ≥0.5 | Activation ≥0.9 |
|---|---:|---:|
| Control | 0/20 | 0/20 |
| RSL3 | 0/20 | 0/20 |
| SDT | 0/20 | 0/20 |

The descriptive SDT/RSL3 immune-kill ratio has median 106.7 and range [65, 535] among 20 defined blocks. Undefined blocks: none.

Zero observed threshold-positive blocks would not establish impossibility. Secondary outcomes are descriptive; no confirmatory secondary significance tests were planned.

## Design and limits

Twenty fixed roots, 42 + block×2^32, retain the original three-arm model and within-block RNG structure. Complete blocks, never cells or arms, are resampled 10,000 times using Python random.Random(20260922). The interval uses linearly interpolated 2.5th and 97.5th percentiles of the mean paired difference. It is conditional on this model, configuration and independent-block interpretation; it excludes parameter, structural and experimental uncertainty.

A source audit verifies disjoint additive seed-address envelopes without u64 wrapping. This removes known cross-block additive aliases; it does not prove statistical independence of pseudorandom streams. The earlier seeds 42–61 reuse seed addresses and their historical bootstrap interval should not be interpreted as validated independent-run uncertainty.

The new binary reproduced both the full historical 33-condition summary and canonical observer bytes before any new block. All 20 planned blocks are included; there was no outcome-dependent replacement, tuning or stopping. The observation window, fixed initial tumor count and O2-independent SDT assumption limit interpretation. This does not establish a causal saturation mechanism or clinical efficacy.

## Reproduction

```bash
python3 scripts/immune_2d_replication_report.py
```

This validates frozen sources and payload hashes, canonical parity, all 20 event ledgers, activation summaries and seed coverage, then recomputes the endpoint, interval and both reports without Rust. `--render-only` refreshes only Markdown from the saved derived JSON; it does not revalidate the archive. JSON and Markdown publication is sequential, not an atomic pair.

- [Frozen protocol](../docs/IMMUNE_2D_REPLICATION_PROTOCOL.md).
- [Plan](../scripts/immune_2d_replication_plan.json).
- [Archive manifest](immune-2d-replication/manifest.json).
- [Complete derived data](immune-2d-replication.json).
- Source freeze: `595269cacc05abf32f484ec8ba2025fd9b317391`.
- Binary SHA-256: `c41a00026426c77b9468f121d2d1714b5869fc1ce7fdcfd04137ee4e5c0cd3ea`.
