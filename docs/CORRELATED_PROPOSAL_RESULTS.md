# Correlated final proposals: gains and remaining failures

The revision **failed the prespecified overall benchmark**: 10 of 12 correlated
runs passed. It passed all three runs on the original rotated box, where the
historical proposal passed one. One annular run and one new shifted-box run
failed a moment/CDF accuracy check. The revision remains experimental and has
not replaced the biological sampler.

The [protocol](CORRELATED_PROPOSAL_PLAN.md), implementation and unit tests were
committed at [`0b26cbdc`](https://github.com/ELares/cancer_research/commit/0b26cbdc9dccd30f7cc1bc7aac25b3b095dde0bc)
before any study pilot or production run. There was one declared configuration,
with no outcome-driven tuning. The [complete report](../analysis/calibration/proposal-correlated-study.md)
and [numerical results](../analysis/calibration/proposal-correlated-study.json)
retain all 36 production arms and 24 analytic controls.

## What changed and what was compared

Each three-arm comparison shares one unchanged 29,184-attempt pilot. Every arm
then makes 8,192 independent attempts from a frozen proposal. The new proposal
uses local full covariance and an exactly evaluable mixture of unrestricted
Gaussians plus a 20% uniform-cube component. Draws outside the cube remain
attempts with zero target weight, without clipping or redrawing.

The unbounded diagonal control uses exactly the same centers, multiplicities,
uniform weight and covariance diagonals as the correlated arm; only the
off-diagonal entries differ. Comparisons between these two arms address that
specific change conditional on the shared pilots. Comparisons with the old
bounded arm also change boundary handling and marginal scales.

| Target | Historical bounded pass | Unbounded diagonal pass | Correlated pass | Historical ESS range | Correlated ESS range |
|---|---:|---:|---:|---:|---:|
| Original rotated box | 1/3 | 1/3 | 3/3 | 185.8–269.6 | 456.0–546.2 |
| Annular cylinder | 3/3 | 3/3 | 2/3 | 457.2–540.1 | 387.9–418.9 |
| Unequal balls | 3/3 | 3/3 | 3/3 | 974.5–1270.2 | 701.6–889.8 |
| New shifted rotated box | 0/3 | 0/3 | 2/3 | 165.1–260.9 | 532.7–581.7 |

The first three geometries are previously observed regression cases with new
seeds. The shifted box is a prospective, rarer target with a different rotation
near the cube boundary. It was designed knowing earlier failures. These are
selected finite-budget challenges, not a representative or blinded sample of
possible targets. The three-seed comparisons are descriptive, without a
superiority test or claim of universal improvement.

## The two correlated-arm failures

| Fixture and seed | Feature | Estimate | Truth | Absolute error | Limit |
|---|---|---:|---:|---:|---:|
| Annulus, 2026092811 | `(1 + sin(2 theta))/2` | 0.583643 | 0.5 | 0.083643 | 0.05 |
| Shifted box, 2026092832 | `P(T6 <= 0.5)` | 0.441601 | 0.5 | 0.058399 | 0.05 |

Both failed only the moment/CDF screen. Every declared region was observed in
all 36 arms, and every shared pilot island reached the final threshold.
All 12 correlated runs passed the accepted-count, ESS, maximum-weight,
relative-MCSE, mass-accuracy and region-mass screens. A high ESS and accurate
regional masses still did not guarantee all the selected distributional checks.
The experiment does not identify whether the remaining errors arose mainly
from pilot coverage, local covariance fitting or production variation. An
[archive-only follow-up](CORRELATED_PROPOSAL_DIAGNOSTICS.md) decomposes the
observed errors and measures conditional precision and single-draw influence
without changing any study outcome.

All 12 full-target oracles passed. All 12 deliberately restricted oracles
passed the usual weight screens but failed truth checks. Their missing support
is intentional; it is not a diagnosis that the learned defensive proposals
have zero density in an unobserved region. The latter have density at least
0.2 throughout the cube and valid finite-variance mass estimates conditional
on their pilots, despite a finite sample's possible imprecision.

## What follows

Correlation-aware production improved ESS on the two selected elongated boxes,
but lowered ESS on the annulus and balls and did not pass the all-run rule.
The archived weighted feature errors have now been
[diagnosed descriptively](CORRELATED_PROPOSAL_DIAGNOSTICS.md). Any regularization
or mixture revision must still be frozen separately before evaluating
it. All four geometries are now seen development cases; another prospective
claim requires newly declared evaluation cases. Do not replace seeds, relax
thresholds or reinterpret this completed benchmark as a pass.

The existing biological sampling results, earlier failed studies and production
simulation matrix are unchanged. No biological inference was rerun and no
calibration tier was upgraded. Mapping model outputs to an independent assay
and obtaining raw experimental replicates remain [separate requirements](INDEPENDENT_ASSAY_CANDIDATE.md).

## Reproduce

    python scripts/run_correlated_study.py --render-only

This verifies source and archive hashes, rebuilds proposals from the archived
pilot endpoints, replays all production draws and recomputes every assessment.
Controls are regenerated from their declared seeds. The twelve compressed
archives retain complete learned attempts; intermediate pilot trajectories
are not independently replayed. The report and protocol state the offline commit-metadata
and floating-point limits. Historical sources and studies retain their hashes.

For a partial study, use the maintained entry point:

    python scripts/run_correlated_study.py --resume

Before generating any missing run, it validates the complete existing archive
inventory, replays each existing record, and checks that their runtime and
implementation provenance agree. An incomplete study also requires the current
Python, NumPy, SciPy, bit generator and implementation commit to match those
records. A mismatch stops before any new run is started; retain the archives
and restore their recorded environment to resume. A complete set can be replayed
under another supported runtime without generating new study runs.

The same entry point without flags starts a study only when no planned archives
already exist. The original `proposal_correlated_study.py` remains unchanged as
part of the frozen numerical implementation. Use the maintained entry point for
new or resumed work; this operational safeguard does not revise the completed
experiment or its failed overall outcome.
