# Additional coverage challenges: results and limits

The unchanged sampler **failed this study's overall acceptance rule**: seven of
nine learned runs passed, while two rotated-box runs failed. All planned runs
and controls are retained. The [protocol](COVERAGE_CHALLENGE_PLAN.md) and
implementation were frozen in commit
[`e8e9b749118dfd51f49e8fbd966b1e1fcab2e589`](https://github.com/ELares/cancer_research/commit/e8e9b749118dfd51f49e8fbd966b1e1fcab2e589)
before evaluation. See the [generated report](../analysis/calibration/proposal-coverage-challenges.md)
and [full numerical results](../analysis/calibration/proposal-coverage-challenges.json).

Each learned run used 29,184 pilot attempts and 8,192 independent production
attempts. These analytic experiments made no biological simulator dose calls.
The proposal, budgets, seeds and thresholds were not changed after seeing results.

## What failed

| Geometry | Learned runs passing | Production ESS range |
|---|---:|---:|
| Correlated rotated box | 1 / 3 | 213.3–317.8 |
| Annular cylinder | 3 / 3 | 426.8–514.8 |
| Unequal disconnected balls | 3 / 3 | 1,162.8–1,205.3 |

The rotated box has coordinate correlation 143/145. Its failed runs were:

| Seed | Accepted | ESS | Maximum normalized weight | Maximum region-mass error | Failed checks |
|---|---:|---:|---:|---:|---|
| 2026092502 | 450 | 221.1 | 0.0222659780 | 0.0822690895 | Maximum weight; region masses |
| 2026092503 | 399 | 213.3 | 0.0271607910 | 0.0328608687 | Maximum weight |

The fixed maximum-weight limit was 0.02; the absolute region-mass error limit
was 0.075. For seed 2026092502, quadrant 3 had estimated mass 0.3322690895
against the exact 0.25. ESS exceeded its threshold of 200 in both failed runs,
but their usual diagnostic screens still failed the separate maximum-weight
check. ESS alone would have missed these declared failures.

Every prespecified region was observed in every learned run. Both islands
reached the final threshold in all nine runs; all learned runs passed the
accepted-count, relative-MCSE, total-mass-error and recorded-moment checks.
The observed rotated-box failures are weight concentration and, in one run,
inaccurate regional mass. They do **not** establish omission of a mode or region.
Passing seven runs cannot rescue the fixed requirement that all nine pass.

## What the controls establish

All nine full-target oracle controls passed. All nine support-hole controls
passed the usual weight screens and failed known-truth screens, as specified.
Both kinds had nominal importance ESS 8,192 because their sampled weights were
equal. The restricted controls deliberately omit half of the rotated or annular
target, or the smaller ball. Their densities integrate to one on the retained
support but do not cover the full target; they are invalid importance proposals
for that target despite their favorable ESS and estimated relative MCSE.

These controls demonstrate the specified checks' response to known support
defects and achievable accuracy with direct target sampling. They are separate
from the learned proposals, whose 20% uniform component covers the entire cube.
No oracle points or labels entered training. Do not attribute the controls'
intentional support holes to the learned runs.

## Interpretation and remaining uncertainty

This is a measured limitation at the fixed computational budget on the selected
correlated fixture. Diagonal kernels or local steps poorly aligned with the
narrow rotated direction are possible explanations, as are uneven pilot
coverage and finite production variation. The study did not isolate those
causes. A proposed explanation must be tested separately before being presented
as the mechanism of failure.

Finite regional masses and recorded moments do not identify the complete
within-region distribution. For example, a hypothetical unequal-ball
distribution with equal probability on the 14 signed coordinate-axis directions
within each ball, the correct radial CDF and correct A/B masses would match all
recorded moments while remaining angularly nonuniform. This is a mathematical limitation of the
checks, not an observed outcome of the sampler. Passing them does not establish
the full angular distribution, precise tails or complete target coverage.

The selected geometries were informed by earlier failures and are now observed
development cases. These results neither prove a corresponding defect in the
biological target nor change its historical scores. They do not provide an
independent assay validation or a pooled posterior. The
[model-to-assay acquisition and mapping requirements](INDEPENDENT_ASSAY_CANDIDATE.md)
remain unresolved.

## Requirements for a subsequent experiment

Retain this completed study unchanged. Before another evaluation, commit a
separate protocol naming the proposed change, fixed budgets, independent seeds,
error criteria and all intended comparisons. Treat reusing these geometries as
development or regression testing; choose any additional evaluation geometries
before observing their results and disclose how prior failures informed them.

If testing correlation-aware proposals, establish the implemented density's
normalization, boundary treatment, full target support and applicable Metropolis
ratios before using importance estimates. Compare against the frozen baseline
under declared budgets, keep pilot and production randomness separate, and
retain failures without replacing seeds or relaxing thresholds. Archive-only
diagnostics can generate hypotheses, but cannot retroactively change this study's
criteria or turn its failed overall outcome into a pass.
