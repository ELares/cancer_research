# Understanding the archived covariance-proposal errors

The completed [comparison](CORRELATED_PROPOSAL_RESULTS.md) passed 10 of 12
correlated-proposal runs and failed its fixed overall benchmark. This follow-up
explains the two failed feature estimates using existing archived attempts.
It adds no pilot, production draw, proposal revision or replacement outcome.

The [reproducible diagnostic](../analysis/calibration/proposal-correlated-diagnostics.md)
retains every declared feature across all 36 production arms and all 12 shared
pilots. The two features highlighted below were selected **after observing the
failures**. Their decompositions are descriptive follow-up evidence, not a new
prospective test.

## What the errors contain

For regions with exact masses `p_r` and conditional feature means `mu_r`, the
signed feature error splits into two terms:

```
estimate - truth
  = sum_r (estimated_region_mass_r - p_r) * mu_r
  + sum_r sum_i_in_region normalized_weight_i * (feature_i - mu_r)
```

The first term measures the contribution from region-mass imbalance. The second
measures the contribution from variation within regions. Computing the second
term directly avoids inventing a conditional sample mean for an empty region.
This is an exact algebraic decomposition of the observed estimate; it does not
identify which part of the sampling procedure caused either contribution.

For the annular target, `angular_sin_2 = (1 + sin(2 theta))/2`. Its exact mean is
0.5; the eight angular sectors have conditional means
`0.5 + [+, +, -, -, +, +, -, -]/pi` and equal target mass. The failed correlated
run's error of **+0.083643** splits into **+0.067252 between sectors** and
**+0.016391 within sectors**. Individually acceptable sector-mass errors can
therefore combine into an unacceptable error in this feature.

For the shifted box, `latent_cdf_6_50` asks whether the sixth latent coordinate
lies in its lower half. Its exact conditional mean is 0.5 in every declared
octant, because the octants partition only the first three latent coordinates.
The failed correlated run's **-0.058399** error is consequently entirely within
octants. The sixth coordinate is a broad, unrotated direction; this failure is
not evidence that a thin rotated direction was missed.

The shifted pilot's endpoint fraction was 0.5625. Its bounded and unbounded
diagonal production estimates were 0.550932 and 0.544112, while the correlated
estimate was 0.441601. This opposite sign cautions against attributing the
production error to the pilot's endpoint marginal alone. The report keeps all
three seeds and all three arms for each highlighted feature.

## Precision and influence

For independent production attempts from a fixed fitted proposal, the report
uses the plug-in self-normalized importance-sampling MCSE:

```
MCSE^2 = n / (n - 1) * sum_i normalized_weight_i^2 * (feature_i - estimate)^2
```

Here `n` includes every attempted draw, including rejected and outside-cube
draws with zero target weight. The estimate is conditional on the fitted
proposal. It excludes pilot variability, can miss unsampled contributions and
is not a calibrated confidence interval or an additional pass/fail criterion.
Error-to-MCSE ratios are descriptive, without a normal-reference test or a
multiple-comparison claim. Fewer than two positive weights make the MCSE
unavailable; a constant feature with sufficient positive weights can have zero
estimated MCSE.
Pilot endpoint summaries retain particle multiplicities and island identities;
they are descriptive and receive no independent-sample MCSE.

The report also computes the exact change from deleting each single weighted
attempt. The largest absolute changes for the two failed estimates are
0.004134 and 0.005315 respectively. Neither is large enough to bring that error
within the original 0.05 bound. This sensitivity calculation removes no samples
from any published estimate; it does not rule out collective influence from a
group of draws or diagnose the source of the error.

## Reproduce and interpret

```bash
python scripts/proposal_correlated_diagnostics.py --render-only
```

Reconstruction checks the completed report's pinned hash, frozen source
specification and the complete archive inventory. It verifies proposal fits,
densities, scores and decisions deterministically, then calculates diagnostics
from the stored points and their actual generating densities. It consumes no
randomness. Unlike the original study's replay, it does not reconstruct random
streams or intermediate pilot paths. A changed or missing input stops before
the diagnostic outputs are replaced.

Freshness checks allow bounded floating-point differences only in named derived
metrics. For error-to-MCSE ratios, a small denominator can amplify a last-bit
difference in the error. An extra allowance requires both errors and MCSEs to
pass the usual comparisons and each ratio to match its own operands. Any extra
ratio difference is bounded by the observed operand differences. This does not relax
archive hashes, analytic truths, counts, classifications or selected attempt
indices, and does not change the recorded diagnostic values.

The next experiment still needs a separately frozen proposal formula,
regularization or mixture coefficients, budgets, streams, controls and decision
rules. The four current geometries are now development cases; a prospective
claim requires newly declared challenges. One production realization per arm
cannot separate pilot variability from conditional production variability.
These diagnostics do not upgrade biological calibration or resolve the
[independent assay requirements](INDEPENDENT_ASSAY_CANDIDATE.md).
