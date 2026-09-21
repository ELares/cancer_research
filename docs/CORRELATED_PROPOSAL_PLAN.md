# Final-proposal covariance experiment

This protocol and its numerical implementation must be committed before any
planned pilot or production run. The driver checks that every numerical source
and this protocol match committed bytes before evaluation. This is a separate
experiment following the [completed geometry challenge](COVERAGE_CHALLENGE_RESULTS.md).
Its failed outcomes and all earlier sampler sources/archives remain unchanged.

## Question and scope

Can a final proposal fitted with local covariance meet the existing accuracy
screens at the fixed budget? The earlier rotated-box failures motivate the
revision but do not establish a cause. We change only independent production
sampling, leaving the pilot training, intermediate proposals and Metropolis
moves unchanged. No biological simulation is run, no production model parameter
is changed, and no pooled posterior or independent assay validation is claimed.

There is one proposed configuration, chosen analytically before evaluating it.
There is no pilot-outcome search or tuning phase. The three old geometries are
**seen regression cases**, even though their seeds here are new. The shifted
rotated box is a new prospective challenge selected with knowledge of earlier
failures, not a blinded or representative sample of possible targets.

## Shared pilots and three final proposals

For each fixture/seed, train the historical `resample_move_local` pilot with its
unchanged two islands, 256 particles, 14 levels and four move sweeps per level:
29,184 total pilot attempts. Both islands must reach the final threshold one.
All three arms use exactly the same final endpoint cloud. No pilot point enters
an importance estimate. No oracle point, target mass, moment or region identity
enters training or fitting.

1. **bounded_diagonal:** the historical final proposal, including its gap-limited
   32-neighbor product Gaussian kernels independently truncated to the cube.
2. **unbounded_diagonal:** the covariance revision below with every off-diagonal
   entry set to zero *after* eigenvalue clipping. It retains the full revision's
   centers, particle multiplicities, marginal variances and uniform weight.
3. **correlated:** the full local covariance revision.

The last two arms isolate the use of off-diagonal covariance conditional on a
shared pilot. Comparing either to the historical arm also changes boundary
handling and potentially marginal variances; that comparison cannot uniquely
attribute differences to correlation or boundary handling.

The covariance fit uses only unique endpoint coordinates for neighborhood
geometry. Standardize each coordinate by its population standard deviation
floored at 0.005. Select at most 32 neighbors by standardized Euclidean distance,
using lexicographically ordered unique centers and a stable distance sort.
Stop before the first positive-distance jump exceeding a factor of three,
using exactly the historical neighborhood convention. Compute the neighborhood
population covariance about its own mean and multiply it by 1.5 squared.
Clip eigenvalues to [0.005 squared, 0.5 squared] and reconstruct a symmetric
covariance. A singleton uses 0.04 squared times the identity. Each original
endpoint remains a kernel center with equal mixture weight: duplicate particles
retain their original multiplicity but do not erase neighborhood variance.

## Exact density and boundary handling

For the two new arms, on all of R^7,

    q(x) = 0.2 * 1_[0,1]^7(x) + (0.8/K) * sum_k Normal_7(x; mu_k, Sigma_k).

Each Gaussian integrates to one on R^7, and the uniform cube has volume one, so
q is normalized without a multivariate truncation constant. Sampling and density
evaluation use the same positive-definite covariance and Cholesky factor in the
existing `importance_sampling.GaussianMixture` implementation.

Every arm makes exactly 8,192 independent production attempts. An unbounded
Gaussian draw outside the cube remains an attempted draw with positive Gaussian
proposal density and **zero target weight**. Do not clip, reflect, drop or redraw
it. The uniform term is zero outside the cube; the Gaussian terms remain.
Every importance weight uses the entire mixture, never only the sampled kernel.
The target is one exactly inside both the cube and the fixture acceptance set.

Consequently q is at least 0.2 everywhere on the cube, and raw target/q weights
are at most five. The mass estimator averages weights over all attempted draws
and is conditionally unbiased for any fixed pilot. Self-normalized region masses
and moments have finite-sample bias. Correct weighting does not ensure that a
finite run observes all important regions. No new Metropolis ratio is needed:
these new proposals are never used in pilot moves.

## Fixtures and fixed random streams

The first three targets, feature functions, masses and regions are unchanged
from [the earlier protocol](COVERAGE_CHALLENGE_PLAN.md): rotated box, annular
cylinder and unequal disconnected balls. The fourth fixture is defined below.

| Fixture | Shared-pilot root seeds | Full-target oracle seeds | Support-hole oracle seeds |
|---|---|---|---|
| rotated_box | 2026092801, 2026092802, 2026092803 | 2026092901, 2026092902, 2026092903 | 2026093001, 2026093002, 2026093003 |
| annular_cylinder | 2026092811, 2026092812, 2026092813 | 2026092911, 2026092912, 2026092913 | 2026093011, 2026093012, 2026093013 |
| unequal_balls | 2026092821, 2026092822, 2026092823 | 2026092921, 2026092922, 2026092923 | 2026093021, 2026093022, 2026093023 |
| shifted_rotated_box | 2026092831, 2026092832, 2026092833 | 2026092931, 2026092932, 2026092933 | 2026093031, 2026093032, 2026093033 |

Each learned root spawns four NumPy SeedSequence children. Child (0,) supplies
the shared pilot (and its island children); (1,), (2,), (3,) respectively supply
the three production arms in the order listed above. Controls use their own
listed roots. All roots are distinct from each other and earlier study roots.
There are 12 shared pilots, 36 production arms and 24 oracle controls. Unit tests
use small artificial clouds/targets and separate seeds; they never train pilots
on a prospective fixture before the freeze.

### Shifted rotated box

Let x = c + Q y, where c=(9/50,9/50,9/50,1/2,1/2,1/2,1/2).
The first three columns of Q in the first three coordinates are

    u=(1,1,1)/sqrt(3), v=(1,-1,0)/sqrt(2), w=(1,1,-2)/sqrt(6).

Q is identity on coordinates four through seven and has no cross-block entries.
It is orthogonal with determinant one. Accept |y_j| <= a_j, with

    a=(1/4,1/48,1/32,1/4,1/4,1/4,1/4).

The score is max_j |y_j|/a_j. First/second-coordinate half-ranges are
1/(4 sqrt(3))+1/(48 sqrt(2))+1/(32 sqrt(6)), about 0.171826; the third is
1/(4 sqrt(3))+1/(16 sqrt(6)), about 0.169853. Both are less than 9/50,
so the entire box lies inside the cube, near its lower faces. Its exact mass
is product_j(2a_j) = **1/12288**. It is rarer than the previous rotated box;
this deliberate change also tests a different orientation and boundary proximity.

Define eight regions by 4*(y1>=0)+2*(y2>=0)+(y3>=0); each has mass 1/8.
The raw mean vector is c, covariance Q diag(a_j squared/3) Q transpose.
Record all seven raw means, raw second moments of x1 and x2, and cross moments
x1*x2 and x2*x3. Also record T_j=(y_j/a_j+1)/2: all seven means (truth 1/2),
all seven second moments (1/3), all 21 pair products (1/4), and indicators
T_j<=t for t in {1/4,1/2,3/4} (truth t). The radius transform
U=(max_j |y_j/a_j|)^7 is uniform: record its mean 1/2 and second moment 1/3.
All assessed features are bounded in [0,1] on the target. These finite checks
still cannot establish the entire distribution or precise tails.

The full oracle draws seven independent uniform latent coordinates and applies
Q and c. The negative oracle retains only y1>=0, samples uniformly there,
and uses its actual restricted density 2/Z. It omits four regions deliberately.
The old fixtures retain their original full/support-hole transforms. Controls
test achievable accuracy and specified support defects; no oracle information
is supplied to learned proposals.

## Fixed decision rules

For each learned arm, retain every earlier threshold:

- At least 20 accepted attempts and importance ESS at least 200.
- Maximum normalized weight at most 0.02; relative mass-estimate MCSE at most 0.10.
- Absolute relative error in exact mass at most 0.15.
- Every declared region observed; maximum absolute region-mass error at most 0.075.
- Every declared moment/CDF feature has absolute error at most 0.05.
- Both shared pilot islands reach the final threshold.

The revision passes this benchmark only if **all 12 correlated runs pass**, all
12 full-target oracles pass, and all 12 support-hole controls pass the usual
weight screens but fail the truth screens. Baseline and diagonal-ablation
outcomes are all reported, including failures; they are not additional gates
for the revision. Arm differences are descriptive across three pilots per
fixture, without a prespecified superiority test, pooling or general coverage claim.

There is no seed replacement, stopping on success, budget extension, threshold
relaxation, or tuning after inspecting outcomes. An implementation defect would
invalidate affected runs, requiring an explicitly documented separate protocol;
scientific failures stay completed failures. Passing these fixtures does not
establish biological accuracy or a new calibration tier.

## Records and replay

Each compressed archive records the source/protocol hashes, implementation
commit, runtime, shared pilot endpoints/history and all three frozen proposals,
production points, component IDs, densities, scores and decisions. Reassembly
checks the fixed specification, fits, independent RNG streams and outputs. Pilot
initial randomness/endpoints and counters are audited by the historical reader;
intermediate trajectories are not archived or independently replayed. It allows
the same narrow numerical replay tolerances as that reader; raw stored acceptance
and counts must match their recomputation. All planned archives are required,
and extras are rejected. Existing archives are never overwritten.

    python scripts/proposal_correlated_study.py
    python scripts/proposal_correlated_study.py --resume
    python scripts/proposal_correlated_study.py --render-only

The first two commands require committed numerical sources. The last uses archived learned attempts and regenerates the declared oracle
controls. It validates source/archive hashes and regenerates both reports
without pilot training. Offline replay verifies source bytes by their hashes;
the recorded commit label is informational and its Git history is not authenticated
by that reader, so replay does not require a Git object database. Both output strings are prepared before replacing any
report. Runtime labels describe archived learned runs; controls are recomputed
in the reader's environment. Existing biological and synthetic archives remain
unchanged, and the [independent-assay gap](INDEPENDENT_ASSAY_CANDIDATE.md) remains open.
