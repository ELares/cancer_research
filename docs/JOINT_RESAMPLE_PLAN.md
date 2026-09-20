# Joint resample/move sampling: fixed evaluation plan

This follow-up responds to the [first frozen-proposal experiment](JOINT_SAMPLING_PLAN.md),
whose pilot cutoffs remained far above the final tolerance. Its results are
historical observations used to design this study, not extra posterior samples.
The numerical implementation and this plan are committed before new biological
production runs. No favorable result is required.

## Scientific target and evaluation held fixed

Keep the seven uniform priors, supported CTRPv2 cohorts and dose grids,
2,000 simulated cells, simulation seed 42, reference vector, and full-precision
criterion `ML162 RMSE + erastin RMSE <= 1.10 * reference_distance`. The currently
computed epsilon is approximately 0.1744797963. ML210 never enters adaptation
or acceptance. Do not insert historical accepted points or the reference vector
into any particle population.

This remains a uniform prior conditioned on a deterministic tolerance event.
The input targets are fitted-curve medians, not raw independent experimental
replicates. Sampling quality is separate from biological or assay validation.

## Pilot islands and fixed production

Run root seeds **2026091911, 2026091912, and 2026091913**. Spawn independent
pilot and production streams from each seed; spawn two independent island
streams from its pilot stream. Every island starts with 256 uniform prior
points, fully scored at the unchanged simulation settings.

Each island executes **14 levels and four independence-Metropolis sweeps per
level**. Select the next intermediate threshold as the larger of final epsilon
and the upper empirical 50% quantile of current scores, capped by the preceding
threshold. The upper quantile (`method="higher"`) retains 129 of 256 distinct
scores; ties are all retained. Resample uniformly among qualifying particles
only if restriction removed particles. Once epsilon is reached, continue the
remaining levels at epsilon without unnecessary resampling.

Fit a bounded proposal to the resampled cloud and freeze it for all four sweeps
of that level. For proposal `y` from current state `x`, accept with probability

```text
I[distance(y) <= intermediate_threshold] * min(1, q(x) / q(y)).
```

The same **whole normalized mixture** enters both densities. A failed Metropolis
coin may skip the model evaluation; otherwise erastin RMSE alone may reject a
candidate when it strictly exceeds the level threshold. Both are exact
shortcuts. Every retained new state has a complete score. Do not reflect a
correlated Gaussian and silently assume the resulting proposal is symmetric.

Each bounded proposal mixes **20% uniform prior** with equal shares of Gaussian
product kernels centered on every particle, independently truncated to `[0,1]`
in each unit-cube coordinate. Kernel widths are 1.5 times the population standard
deviation of the nearest 32 **unique** centers (or all unique centers if fewer),
with floor 0.005 and ceiling 0.5. Neighbors use Euclidean distance after scaling
coordinates by the unique cloud's marginal standard deviations, floored at
0.005. Duplicate particles retain their mixture multiplicities but cannot
artificially collapse neighborhood variance. Each coordinate's truncation mass
is included in its density; sampling uses its inverse conditional CDF without
clipping or rejection.

The fixed pilot budget is **29,184 attempted model points per root seed**:
512 initial points plus `2 * 14 * 4 * 256` proposed moves. Some moves need no
model evaluation. Archive actual dose-call costs, thresholds, accepted moves,
unique particles and surviving initial ancestors; none is an independent
posterior sample count.

After all levels, fit the same bounded kernel mixture to the combined 512 final
particles and freeze it. Draw **8,192 independent production attempts** from the
separate production stream. Weight each accepted point by `1/q(u)` and every
rejected point by zero. Production erastin-first rejection uses only the fixed
final epsilon. No pilot point, clone, or resampled draw enters these estimates.
If a pilot does not reach epsilon, still run the fixed production budget and
report the failure; do not widen the final criterion or extend the budget.

## Known-target checks before biological production

Use the [synthetic study](../scripts/proposal_synthetic_validation.py), whose
truths and acceptance gates are defined before its learned-proposal evaluations.
The main fixture is a seven-dimensional union of two disjoint, unequal-volume
boxes, one touching the prior boundary. Its exact prior mass is `1/4096`, and
the smaller box has conditional mass `1/4`. A boundary-slab fixture checks that
six uninformed coordinates retain their uniform distributions.

Run positive seeds **2026092101–2026092103**, with the same pilot design and
8,192 independent production attempts per fixture and seed. Require the usual
per-run screens below, both expected regions represented, relative normalizer
error at most 15%, conditional region-mass error at most 0.075, and absolute
errors of the specified means, second moments and cross moments at most 0.05.
These are attainable finite-fixture checks, not guarantees for learned proposals
or the biological model. An analytic oracle is used only to contextualize their
sampling noise, never to initialize a learned proposal.

Also retain deliberately incomplete proposals centered only on the larger box,
with seeds **2026092001–2026092003**. These negative controls were explored while
designing the study and are not blind tests. The known-truth checks must reject
them. Some can pass ordinary ESS/weight/MCSE screens despite missing one quarter
of the target mass; report all outcomes, including any ancillary screen failure.
Passing the known-target study is required before biological production. A failed
fixture requires an explicitly documented design revision, not silent gate or
seed selection. Reusing a fixture in method development is not independent
validation of general mode-discovery reliability.

## Prespecified biological sampling screens

Keep the [first experiment's](JOINT_SAMPLING_PLAN.md) four per-run screens:
accepted count at least 20, importance-weight ESS at least 200, maximum normalized
weight at most 0.02, and relative normalizer MCSE at most 0.10.

Across all three root seeds, require every island to reach epsilon, all per-run
screens to pass, each parameter's median range at most 5% of its prior width,
each 2.5% and 97.5% quantile range at most 10%, held-out weighted-median RMSE
range at most 0.02, and each pair of normalizer estimates to differ by at most
four combined MCSEs. Report every check, even if another has already failed.

No criterion, seed, or budget is changed after seeing biological results.
Publish run-wise diagnostics rather than a pooled posterior, even on a pass.
Retain the original uniform and first importance-sampling experiments separately;
compare observed ESS per dose call including pilot/reference costs without
claiming a controlled wall-clock speedup.

## Reproducibility, scope and method sources

Archive initial and final pilot clouds/scores, per-level diagnostics, final
proposal, all production attempts/densities/decisions and evaluated curves, seeds,
and code/data/runtime hashes. Pilot movement trajectories are reproducible by
rerunning the seed and model but are not stored individually. Offline checks
replay production streams, reconstruct proposals from final clouds, reconcile
pilot costs, recalculate every complete production score or strict early bound,
and rebuild summaries from the archived attempts. Changing numerical
source during a run is an error. Historical sampling files stay unchanged.

The nested-region pilot follows the idea of
[Chiachio et al., ABC by Subset Simulation (2014)](https://arxiv.org/abs/1404.6225),
and adaptive population methods are discussed by
[Del Moral, Doucet and Jasra (2012)](https://www.cs.ubc.ca/~arnaud/delmoral_doucet_jasra_adaptiveSMCforABC.pdf).
This implementation uses those ideas for **proposal training only**; it does not
claim to reproduce either paper's algorithm or evidence estimator. Final inference
uses ordinary independent importance sampling with a defensive prior component,
as in the [previous method rationale](JOINT_SAMPLING_PLAN.md). The bounded density
uses the usual [truncated-normal normalization](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.truncnorm.html).
All budgets, kernel choices and operational screens are this repository's choices.

Neither high ESS, island diversity, three agreeing runs, nor success on these
fixtures excludes a shared unseen region in the real target. The uniform component
preserves support and bounds raw weights by five, but at the historical uniform
acceptance rate it yields only about 0.16 acceptable prior draws per production
run. It is weak protection against missed rare regions. Normalizer and CDF MCSEs
are conditional sampling diagnostics; zero estimated CDF error at an observed
extreme is not proof of a precisely estimated tail. Independent raw-assay
validation and the model-to-assay mapping remain pending.
