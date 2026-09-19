# First bounded-pilot validation: failed prerequisite

The implementation and [evaluation plan](JOINT_RESAMPLE_PLAN.md) were frozen in
commit `49a3934e` before the six learned synthetic runs. The complete
[study](../analysis/calibration/proposal-synthetic-validation.md) retains all
results and the three previously declared negative controls. No biological
production run was performed under this design.

All three boundary-slab runs passed. All three separated-region runs failed:

| Seed | Accepted A / B | ESS | Relative mass error | Usual weight screens | Known-truth screens |
|---|---:|---:|---:|---|---|
| 2026092101 | 3 / 1133 | 35.0 | +0.1% | Fail | Pass |
| 2026092102 | 0 / 1461 | 1101.0 | −28.5% | Pass | Fail |
| 2026092103 | 1 / 1451 | 8.4 | +16.0% | Fail | Fail |

The accepted target has two separated regions with true conditional masses
25% and 75%. Every pilot island reached the final threshold, but the combined
clouds retained only 6, 2 and 5 unique points in the smaller region. Their
54, 10 and 7 total particles included many resampling copies.

The bandwidth rule uses 32 unique neighbors. For these smaller-region centers,
it necessarily included points from the distant region. This explains the broad
kernels and weak production coverage there: density at the smaller region's
center was approximately 5.85, 5.02 and 1.48, versus 1649, 2597 and 2289 at the
larger region's center. These are post-run diagnostics, not prespecified gates
or a complete causal isolation of each algorithmic choice.

Seed 2026092102 illustrates why an apparently good ESS and small reported MCSE
cannot establish coverage. It missed a region containing one quarter of the
target and underestimated total mass by 28.5%, despite passing every usual
weight screen. Independent analytic truth caught the failure. The biological
driver therefore refuses to run with this prerequisite.

A revised design must have its own recorded protocol, source identity and
fresh validation seeds. The failed numerical implementation and artifacts stay
available for offline replay; these fixtures are now development challenges,
not unseen evidence of general sampling reliability. The scientific target,
held-out compound and biological adequacy screens remain unchanged.
