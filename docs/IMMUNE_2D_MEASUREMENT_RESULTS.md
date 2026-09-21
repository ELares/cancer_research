# What the historical 2D immune comparison actually measures

The [passive event report](../analysis/immune-2d-measurement-report.md)
reproduces the historical 521 SDT versus 5 RSL3 immune kills without changing
the simulation. It does **not** support the manuscript's former explanation
that the 2D immune response operated deep in saturation. No eligible
cell-step opportunity reached the model's half-maximal activation threshold.
This corrects an interpretation of existing model results; it adds no
experimental validation.

The [protocol](IMMUNE_2D_MEASUREMENT_PROTOCOL.md), configuration, code, and
thresholds were committed before capture at
`5872fd3e89b509b9e28f9e309cd058361147a098`. The complete 33-condition production
summary remains byte-identical to the uninstrumented baseline. The separate
[3D archive](../analysis/immune-measurement-report.md) is unchanged.

## Exposure during the immune update

An eligible opportunity is a living tumor cell immediately before an immune
update, at step 60 or later, with local DAMP at least 0.01. Repeated
opportunities for a cell are counted separately. The activation function is
`DAMP / (DAMP + 50)`; means below average this function over the actual
opportunities, rather than applying it to mean DAMP.

| Measurement | RSL3 | SDT |
|---|---:|---:|
| Unique eligible cells | 14,521 | 23,510 |
| Eligible cell-step opportunities | 770,701 | 2,314,460 |
| Immune kills | 5 | 521 |
| Kills per eligible cell-step opportunity | 0.00000649 | 0.000225 |
| Mean local DAMP per opportunity | 0.06584 | 2.06872 |
| Maximum local DAMP among opportunities | 2.58529 | 10.50046 |
| Mean activation per opportunity | 0.001309 | 0.038506 |
| Opportunities with DAMP ≥50 (activation ≥0.5) | 0 | 0 |
| Opportunities with DAMP ≥450 (activation ≥0.9) | 0 | 0 |

The largest observed SDT activation is about 0.174. These measurements rule
out deep saturation **during eligible immune updates in this realization**.
They do not describe every cell at every time. The report also retains the
control arm, including its zero immune kills and 78 eligible opportunities.

SDT has more eligible opportunities and higher mean activation. These are
descriptions of populations shaped by different treatment histories,
ferroptotic deaths, and random streams. Dividing kills by either opportunities
or unique cells does not isolate a causal treatment effect or immunogenicity
per dead cell. The existing final non-ferroptotic-pool normalization uses yet
another denominator and cannot recover the earlier eligibility history.

## Death, release, and the terminal heatmap

SDT has 139,640 completed releases with mean release LP 19.12298. RSL3 has
160 completed releases with mean release LP 18.26580, plus three deaths whose
scheduled releases fall at or beyond the horizon. Their separate horizon LP mean
is 15.79113. These conditional cohort means are close compared with the total
kill ratio; they are not measurements of biological DAMP potency.

The three censored RSL3 deaths contribute 47.37338 DAMP units **after the last
immune update**. The terminal RSL3 field peak is 19.75132, while the largest
DAMP exposure at an eligible immune update is 2.58529. For SDT the terminal
peak is 0.22721, whereas the largest eligible exposure is 10.50046. Thus a
terminal heatmap cannot stand in for exposure during earlier immune updates.
The archive distinguishes cumulative injections, cleared field mass, and
terminal additions; all event populations and mass accounting reconcile.

## What remains unresolved

This capture uses the canonical geometry seed and one runtime seed per arm.
It measures the already-known headline realization; it is not a new
independent replication study. The earlier 20-seed count-only dispersion
study remains separate. Cells and repeated opportunities are not independent
experimental replicates, so no uncertainty intervals are inferred from their
large counts.

The 2D and 3D runs differ in geometry, population size, transport settings,
and random streams. Their roughly 104:1 and 4:1 kill ratios do not isolate a
geometry or dilution effect. A causal explanation needs a separately designed
comparison. LP and DAMP are uncalibrated model units; no physical clock or
experimental assay mapping is established here. The registered P5 within-1.5×
criterion remains unchanged and experimentally untested.

Next, separately freeze a replicate study with a defined independent unit
and observable, while resolving the [independent assay mapping and metadata
gap](INDEPENDENT_ASSAY_CANDIDATE.md). Preserve both measurement archives and
publish all declared outcomes, including null results.

## Reconstruct the evidence

Run `python3 scripts/immune_2d_measurement_report.py` from the repository root.
It validates the frozen archive and rebuilds the report without compiling
Rust or running a simulation. The [manifest](../analysis/immune-2d-measurements/manifest.json)
records source and binary hashes, toolchain, baseline identity, and capture
time. Raw observations record geometry and runtime seeds. The archived raw
events and aggregates support reconciliation;
individual DAMP values for every eligible cell-step are not archived.

After capture, the shared observer was moved to
`simulations/observers/immune_measurements.rs`, keeping the biological core
sources identical to those used by the older frozen sampling studies. Full
compatibility reruns preserved the complete 2D baseline and both 2D and 3D
observation JSON files byte for byte. The capture's source archive still
records the original implementation and layout.
