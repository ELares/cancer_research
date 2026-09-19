# Joint sampling: fixed evaluation plan

This plan precedes the production importance-sampling runs. It follows the
[dose-support correction](CALIBRATION_DOSE_SUPPORT.md), whose ordinary rejection
run accepted 4 of 40,000 prior draws. It tests computational sampling adequacy;
it does not establish biological validity or independent assay validation.

## Target held fixed

Use the seven existing uniform priors, supported CTRPv2 cohorts and dose grids,
2,000 simulated cells, simulation seed 42, and the existing reference vector.
The final criterion remains `ML162 RMSE + erastin RMSE <= 1.10 * reference_distance`.
Recompute that distance at full precision from the current CSV; do not use the
rounded JSON threshold. ML210 never enters pilot selection or acceptance.

The target is the uniform prior conditioned on this deterministic tolerance
event. The simulator seed is fixed, and the target curves are fitted population
summaries. It is not a likelihood for raw experimental replicates.

## Proposal and production

Run three independent pilot-plus-production experiments with root seeds
`2026091901`, `2026091902`, and `2026091903`. Spawn separate pilot and production
random streams within each experiment. No historical accepted draw or hand-tuned
reference vector is inserted into a proposal or accepted sample.

Map parameters affinely to the unit cube. Four sequential pilot rounds each
attempt 1,024 points. The first proposal is uniform. Subsequent proposals mix
20% uniform with equal shares of the Gaussian components fitted in prior rounds.
Select the best 20% of finite, in-prior pilot scores (at least 32 when available).
Fit their mean and population covariance using inverse-current-proposal weights;
use twice that covariance plus `1e-6 * identity`. Keep every fitted component.
Pilot selection thresholds tune only the proposal; they never change the final
acceptance criterion. Pilot points do not enter production estimates.

Freeze the resulting mixture and attempt 8,192 independent production points.
Gaussian proposals are untruncated: points outside the prior cube count as
zero-weight attempts and are not resampled. For an in-prior accepted point,
the unnormalized importance weight is `1 / q(u)`, where `q` is the **whole**
frozen mixture density, including the uniform component. All other weights are
zero. This preserves the original target despite concentrating proposals.

For production only, evaluate erastin first and skip ML162 when erastin RMSE
alone exceeds the fixed threshold. RMSEs are nonnegative, so this is an exact
rejection. Record the lower bound and skipped computation explicitly. Compute
complete curves for every accepted point. Their GPX4 curves also supply the
ML210 comparison; this model uses the same dose map for those compounds.

## Prespecified diagnostic screens

Each run must have at least the existing 20 accepted points, importance-weight
ESS at least 200, maximum normalized weight at most 0.02, and estimated relative
Monte Carlo error of the normalizing constant at most 0.10. The last screen is
largely redundant with ESS and is reported for interpretability.

Across all three runs, require:

- Each parameter's median range at most 5% of its prior width.
- Each parameter's 2.5% and 97.5% quantile ranges at most 10% of its prior width.
- Held-out RMSE weighted-median range at most 0.02 viability units.
- Pairwise normalizer differences at most four combined Monte Carlo standard
  errors.

These are operational screens, not calibrated statistical tests or proof that
all modes were found. Weighted CDF Monte Carlo errors at the reported quantiles
are recorded; ESS 200 alone does not guarantee precise tails. If any screen
fails, report the specific shortfall and withhold a pooled inferential summary.
No budget, seed, tolerance or screen is relaxed after viewing results. Longer
or redesigned experiments require a separately documented plan.

The uniform component ensures support and bounds raw weights by five. At the
historical rejection rate of 4/40,000, however, each production run expects only
about 0.16 acceptable uniform-component draws. It is weak protection against an
undiscovered mode at this budget. Agreement among three adaptive pilots cannot
exclude a shared blind spot, and weight ESS is not a literal independent
posterior sample count; see [Elvira et al.](https://arxiv.org/abs/1809.04129).

Compare ESS per simulator dose call (including pilots and the reference) and
accepted counts with the existing uniform-rejection baseline. Different budgets
and exact rejection shortcuts preclude claiming a controlled wall-clock speedup.
The old four-draw run is not a posterior oracle. Known synthetic conditional
distributions test the weighting independently of the biological simulator.

## Reproducibility and interpretation

Archive the proposal parameters, all production attempts and densities, accepted
curves, target and code hashes, seeds, runtime environment, and per-run results.
Render summaries from this archive; no unweighted resampling is used to inflate
the sample count. Retain the original rejection result as the historical
baseline. Do not feed weighted intervals into the existing equal-weight
information-width null or change the production simulation parameters.

The method uses ordinary importance weights and a prior component to preserve
support; these choices follow the importance-sampling correction discussed by
[Beaumont et al. (2009)](https://arxiv.org/abs/0805.2256) and the mixture approach
of [Owen and Zhou (2000)](https://doi.org/10.1080/01621459.2000.10473909).
The prior-component construction is also described by
[Hesterberg (1995)](https://www.stat.cmu.edu/technometrics/90-00/vol-37-02/v3702185.pdf).
The particular pilot heuristic and diagnostic thresholds above are this
repository's experimental choices, not guarantees supplied by those papers.
