# Archived covariance proposal feature diagnostics

The completed benchmark remains **failed: 10/12 correlated runs passed**. Its estimates, thresholds and outcomes are unchanged.

This is a post hoc descriptive analysis of all 12 shared pilots and all 36 production arms.
The two highlighted features were selected after inspecting the completed failures.
Every declared feature is retained in the numerical report. No new pilot, production draw,
oracle draw, or random-number stream is run, and no biological sampler is changed.

## Highlighted feature errors across every paired pilot

| Fixture | Seed | Arm | Estimate | Signed error | Conditional MCSE | Largest deletion change | Attempt index |
|---|---:|---|---:|---:|---:|---:|---:|
| annular_cylinder | 2026092811 | bounded_diagonal | 0.494058 | -0.005942 | 0.017496 | -0.003946 | 942 |
| annular_cylinder | 2026092811 | unbounded_diagonal | 0.511171 | +0.011171 | 0.017122 | +0.005408 | 1991 |
| annular_cylinder | 2026092811 | correlated | 0.583643 | +0.083643 | 0.017428 | +0.004134 | 2057 |
| annular_cylinder | 2026092812 | bounded_diagonal | 0.518214 | +0.018214 | 0.015086 | +0.002252 | 301 |
| annular_cylinder | 2026092812 | unbounded_diagonal | 0.511982 | +0.011982 | 0.016620 | -0.004279 | 3362 |
| annular_cylinder | 2026092812 | correlated | 0.492858 | -0.007142 | 0.016199 | -0.002877 | 3948 |
| annular_cylinder | 2026092813 | bounded_diagonal | 0.480009 | -0.019991 | 0.014363 | +0.002206 | 1494 |
| annular_cylinder | 2026092813 | unbounded_diagonal | 0.485824 | -0.014176 | 0.017064 | +0.002637 | 606 |
| annular_cylinder | 2026092813 | correlated | 0.486804 | -0.013196 | 0.016897 | +0.005914 | 197 |
| shifted_rotated_box | 2026092831 | bounded_diagonal | 0.513109 | +0.013109 | 0.038501 | -0.023749 | 7737 |
| shifted_rotated_box | 2026092831 | unbounded_diagonal | 0.502347 | +0.002347 | 0.036024 | -0.013755 | 4601 |
| shifted_rotated_box | 2026092831 | correlated | 0.489673 | -0.010327 | 0.021668 | +0.004610 | 1580 |
| shifted_rotated_box | 2026092832 | bounded_diagonal | 0.550932 | +0.050932 | 0.035076 | +0.015112 | 1279 |
| shifted_rotated_box | 2026092832 | unbounded_diagonal | 0.544112 | +0.044112 | 0.036263 | -0.012784 | 7461 |
| shifted_rotated_box | 2026092832 | correlated | 0.441601 | -0.058399 | 0.020236 | +0.005315 | 1099 |
| shifted_rotated_box | 2026092833 | bounded_diagonal | 0.514891 | +0.014891 | 0.030758 | -0.007061 | 6159 |
| shifted_rotated_box | 2026092833 | unbounded_diagonal | 0.525334 | +0.025334 | 0.040022 | -0.016126 | 6200 |
| shifted_rotated_box | 2026092833 | correlated | 0.508174 | +0.008174 | 0.021536 | -0.004830 | 949 |

The annular feature is `(1 + sin(2 theta))/2`; the shifted-box feature is `P(T6 <= 0.5)`.
Both full-target truths are 0.5. Attempt indices are zero-based original production indices.
Deletion changes are the signed change in the self-normalized estimate after removing one
attempt, selecting the largest absolute change and the lowest index for an exact tie.
They are sensitivity summaries only; no point is removed from an original result.

## Exact decomposition of the selected errors

| Fixture | Seed | Arm | Between-region contribution | Within-region contribution | Total error |
|---|---:|---|---:|---:|---:|
| annular_cylinder | 2026092811 | bounded_diagonal | -0.010148 | +0.004206 | -0.005942 |
| annular_cylinder | 2026092811 | unbounded_diagonal | +0.011405 | -0.000233 | +0.011171 |
| annular_cylinder | 2026092811 | correlated | +0.067252 | +0.016391 | +0.083643 |
| annular_cylinder | 2026092812 | bounded_diagonal | +0.010728 | +0.007486 | +0.018214 |
| annular_cylinder | 2026092812 | unbounded_diagonal | +0.003196 | +0.008786 | +0.011982 |
| annular_cylinder | 2026092812 | correlated | -0.007842 | +0.000700 | -0.007142 |
| annular_cylinder | 2026092813 | bounded_diagonal | -0.021454 | +0.001464 | -0.019991 |
| annular_cylinder | 2026092813 | unbounded_diagonal | -0.020557 | +0.006381 | -0.014176 |
| annular_cylinder | 2026092813 | correlated | -0.003645 | -0.009550 | -0.013196 |
| shifted_rotated_box | 2026092831 | bounded_diagonal | +0.000000 | +0.013109 | +0.013109 |
| shifted_rotated_box | 2026092831 | unbounded_diagonal | +0.000000 | +0.002347 | +0.002347 |
| shifted_rotated_box | 2026092831 | correlated | +0.000000 | -0.010327 | -0.010327 |
| shifted_rotated_box | 2026092832 | bounded_diagonal | +0.000000 | +0.050932 | +0.050932 |
| shifted_rotated_box | 2026092832 | unbounded_diagonal | +0.000000 | +0.044112 | +0.044112 |
| shifted_rotated_box | 2026092832 | correlated | +0.000000 | -0.058399 | -0.058399 |
| shifted_rotated_box | 2026092833 | bounded_diagonal | +0.000000 | +0.014891 | +0.014891 |
| shifted_rotated_box | 2026092833 | unbounded_diagonal | +0.000000 | +0.025334 | +0.025334 |
| shifted_rotated_box | 2026092833 | correlated | +0.000000 | +0.008174 | +0.008174 |

For each region r, the error is `(p_hat_r - p_r) * mu_r + sum_i_in_r w_i * (f_i - mu_r)`,
summed across regions. The eight annular sector conditional means are `0.5 +/- 1/pi`,
with signs `++--++--`. The shifted feature has conditional mean 0.5 in every octant.
Empty regions have an unavailable conditional estimate and a zero direct within-region term.
Detailed per-region counts, masses and contributions remain in the JSON.

The annular failing run has a substantial combined sector-mass contribution even though
every individual region passed its original error limit. The shifted failure lies within
the declared octants; their identities involve different latent coordinates from T6.

## Descriptive pilot endpoint summaries

| Fixture | Seed | Pooled selected-feature mean | Island 0 mean | Island 1 mean | Unique endpoints / total |
|---|---:|---:|---:|---:|---:|
| annular_cylinder | 2026092811 | 0.517569 | 0.498549 | 0.536588 | 448 / 512 |
| annular_cylinder | 2026092812 | 0.491159 | 0.503039 | 0.479279 | 488 / 512 |
| annular_cylinder | 2026092813 | 0.491916 | 0.439666 | 0.544167 | 475 / 512 |
| shifted_rotated_box | 2026092831 | 0.429688 | 0.398438 | 0.460938 | 332 / 512 |
| shifted_rotated_box | 2026092832 | 0.562500 | 0.523438 | 0.601562 | 335 / 512 |
| shifted_rotated_box | 2026092833 | 0.439453 | 0.328125 | 0.550781 | 286 / 512 |

Endpoint means and region fractions retain particle multiplicities. Unique endpoint counts
describe repetition, not independent information. The JSON retains these summaries for every
fixture, island and declared feature. Pilot endpoints are dependent training particles;
no inferential ESS, MCSE, or posterior estimate is assigned to them.

## Numerical definitions and limits

Production weights use each archived draw’s actual full generating density, after independent
density validation. Rejected and outside-cube attempts retain zero weight and stay in n.
For normalized weights w, the feature variance estimate is
`n/(n-1) * sum(w_i^2 * (f_i - estimate)^2)`. Its square root is the delta-method MCSE
conditional on the frozen proposal. The JSON also records signed error divided by this MCSE
(unavailable when MCSE is zero or unavailable), plus ESS, maximum weight and top-1/5/20 weight shares.
MCSE and the largest deletion are unavailable when only one attempt has positive weight.
These are descriptive diagnostics, not new gates or calibrated significance tests after selection.
They do not measure pilot uncertainty or certify tails absent from the realized sample.

The shared pilots permit descriptive arm comparisons, but production draws use different streams.
Neither the decompositions, endpoint imbalances nor deletion sensitivities identify a causal
failure mechanism. One realized production sample per arm cannot separate pilot effects,
proposal geometry and production variation. No new calibration tier or coverage claim follows.

Any covariance regularization or mixture revision requires a separate frozen protocol before
evaluation, with fixed coefficients, budgets, seeds, controls and decision rules. All four
current geometries are seen development cases; prospective claims require new declared cases.

## Inputs and reproduction

Completed study SHA256: `df4bbab24034f9bea804ce7425a5c7a8cf6849ab1b60f498d33cf04e80a5d153`.
Historical numerical/protocol commit: `0b26cbdc9dccd30f7cc1bc7aac25b3b095dde0bc`.
Source bytes and the complete archive set are verified against the pinned completed report.
Endpoint arrays, proposal fits, full mixture densities, scores, acceptance decisions and
original assessments are recomputed without sampling. Component IDs and stream identities
are validated as archived metadata; generating randomness and intermediate pilot paths are
not replayed here. Archived runtime labels describe the original experiment only.
Analytic truths, identities, counts, flags and selected attempt indices remain exact;
derived numerical comparisons allow only the established narrow floating-point tolerance.

    python scripts/proposal_correlated_diagnostics.py
    python scripts/proposal_correlated_diagnostics.py --render-only

Both commands reread and validate the archives. Render-only validates stored provenance and
rebuilds all diagnostics rather than trusting cached derived values. Both output strings are
prepared before either report is replaced. The completed study and every original archive
remain unchanged.
