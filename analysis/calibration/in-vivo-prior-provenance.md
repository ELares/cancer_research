# What the "in-vivo priors" actually are (#PRIOR-PROVENANCE)

## The claim being examined

The only data-anchored calibration leg rejects the shipped defaults. Fitting the
single-cell kill switch to CTRPv2 GPX4-inhibitor dose-response requires
`lp_propagation` 0.10 → ~0.78 and `lp_rate` 0.06 → ~0.71, and the default
parameters score RMSE 0.5666 against the fitted 0.0504 — an eleven-fold gap
(`analysis/calibration/kill-switch-calibration.md`,
`analysis/calibration/joint-posterior.md`).

The repository's stated reason this does not invalidate the spatial and headline
results is a **regime disjunction**: the in-vitro posterior lies entirely above
the "in-vivo PRCC priors" that drive the prior-predictive intervals, so in-vitro
data cannot condition in-vivo outputs.

| parameter | "in-vivo PRCC range" | in-vitro posterior (95%) |
|---|---|---|
| `lp_propagation` | [0.05, 0.20] | 0.398 – 0.99 |
| `lp_rate` | [0.03, 0.12] | 0.386 – 0.98 |

The disjunction itself is real and correctly measured. This note is about what
the left-hand column is.

## What the left-hand column is

`scripts/run_prcc.py:38` states it directly:

```python
# Ranges: ±50% of default, capped at biologically meaningful bounds.
```

Against the defaults in `simulations/ferroptosis-core/src/params.rs:337-338`
(`lp_rate: 0.06`, `lp_propagation: 0.10`):

* `lp_rate` [0.03, 0.12] is **exactly ±50% of 0.06**;
* `lp_propagation` [0.05, 0.20] is −50% / +100% of 0.10, the upper end widened
  by the "biologically meaningful bounds" cap.

So the ranges are **not independently derived in-vivo measurements**. They are
the default parameters' own neighbourhood. And `context="2d"`, which every
calibration, ABC and uncertainty script runs, resolves to `Params::default()`
(`simulations/ferroptosis-python/src/lib.rs:58-60`) — the same defaults.

## Why that matters

Naming that column "in-vivo priors" gives it an authority it has not earned, and
using the disjunction to bound the falsification is circular. Unpacked, the
argument becomes:

> The in-vitro data rejects our default parameters. Our default parameters'
> ±50% neighbourhood does not overlap the in-vitro fit. Therefore the in-vitro
> data does not bear on results computed from our default parameters.

The second sentence is a restatement of the first, not independent evidence. A
genuine regime disjunction would need the in-vivo ranges to come from in-vivo
measurements; these come from the defaults being questioned.

## What is and is not affected

**Not affected:** the disjunction as a *numerical* statement, the ABC posterior,
the held-out ML210 validation, and the conclusion that the in-vitro posterior
cannot condition the spatial headlines. All of that stands. So does the
possibility that the true in-vivo values genuinely differ from in-vitro ones —
that is a reasonable prior belief about ferroptosis, and nothing here refutes it.

**Affected:** the standing of that belief as *evidence*. The honest position is:

* the one leg anchored to independent data rejects the defaults by 11x RMSE;
* no independent in-vivo measurement of these rate constants exists in this
  repository, so the in-vivo regime is asserted, not measured;
* therefore the spatial and headline numbers rest on parameters that are
  unfalsified only because they are untested in their own regime, not because a
  test has cleared them.

That is a weaker position than "the regimes are disjoint" implies, and it is the
one the reader should be given.

## What would settle it

Either an in-vivo ferroptosis dataset that maps onto these dimensionless
observables — which the repository has looked for and documented as not publicly
existing — or re-deriving every headline at the fitted cascade and reporting both,
so a reader can see which directions survive crossing the bistable tipping point.
The second is achievable now and is the cheaper of the two.

## It has now been done, and it settles the question differently than expected

`scripts/headline_at_fitted.py` re-derives every headline under three parameter
sets — the in-vivo defaults, the #330 CTRPv2 point fit, and the #500 joint
posterior medians — and writes `analysis/headline-at-fitted-cascade.md`.

The expected outcome was a list of directions that survive the crossing and
directions that do not. That is not what came back. **Both fitted sets are
inadmissible for these models.** Under either, the untreated
Persister/Control population dies at 99.97%, against the
model's own stated constraint of under 2% — a factor of
50. Every headline then degenerates: the Bliss
ratio is exactly 1.0 because both single arms saturate, and all three penetration
tissues kill 100%, so their ordering is trivially preserved and says nothing.

So the disjunction documented above is no longer only an argument about parameter
ranges. Substituting the in-vitro values into the in-vivo and spatial models does
not merely fail to condition them — it breaks the baseline-viability constraint
outright. The cheaper of the two routes has been taken and it does not lead where
it was expected to: it rules the substitution out rather than grading the
directions. Conditioning these headlines still requires an in-vivo dataset, and
this is now demonstrated rather than argued.

One direction is worth recording despite the degeneracy: the hypoxia
kill-collapse gap stays positive at all three sets
(0.866, 0.606,
0.498). That is weak evidence — a
positive gap in a saturated regime is nearly guaranteed — but it is the only
headline that does not become meaningless, and it points the same way the
manuscript reports.

## Provenance

Every claim above is checkable from committed files: `scripts/run_prcc.py:38-39`,
`simulations/ferroptosis-core/src/params.rs:337-338`,
`simulations/ferroptosis-python/src/lib.rs:58-60`,
`analysis/prcc-results.json` (`metadata.parameter_ranges`),
`analysis/calibration/joint-posterior.md`.
