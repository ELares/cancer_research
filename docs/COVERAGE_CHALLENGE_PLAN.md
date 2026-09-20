# Frozen sampler: additional geometry challenges

This protocol and its new numerical implementation must be committed before any
of the prospective runs below. The local-move sampler is the exact implementation
used in the [previous study](JOINT_RESAMPLE_LOCAL_PLAN.md). Its seven-dimensional
separated boxes and boundary slab informed development. The three geometries here
were selected analytically, without running that sampler on them, to challenge
correlation, curvature and unequal disconnected scales.

These are deliberately chosen challenges informed by known failure modes, not a
random or blinded sample of possible targets. No biological parameters, assay
mapping, source files or historical experiment records are changed. The
[independent-assay acquisition gap](INDEPENDENT_ASSAY_CANDIDATE.md) remains open.
Results belong in [the separate report](../analysis/calibration/proposal-coverage-challenges.md).

## Fixed target definitions

The prior is uniform on the seven-dimensional unit cube. The final threshold is
one. Every analytic score below is evaluated completely, with zero simulator dose
calls; intermediate thresholds change only the pilot conditioning region.
All intervals include their endpoints, with partition ties assigned as specified.

### Correlated rotated box

Let y1 = (x1+x2-1)/sqrt(2), y2 = (x1-x2)/sqrt(2), and yj = xj-1/2 for j=3..7.
The accepted region has |yj| <= aj, with a1=3/8, a2=1/32 and a3..7=1/6.
The score is max_j |yj|/aj. Rotation has unit absolute determinant, and the
entire box is inside the cube: its first two coordinates differ from 1/2 by
at most (3/8+1/32)/sqrt(2) < 1/2.

Its exact mass is (2a1)(2a2) product_j=3..7(2aj) = **1/5184**.
The correlation of x1 and x2 is **143/145**. Four regions have equal conditional
mass 1/4, indexed by 2*(y1>=0)+(y2>=0).

All coordinate means are 1/2. E[x1²]=E[x2²]=1681/6144,
E[x1*x2]=1679/6144 and E[x2*x3]=1/4. In addition to these eleven raw-coordinate
moments, record all seven E[Tj] where Tj=(yj/aj+1)/2, both E[T1²] and E[T2²],
and E[T1*T2]. Their truths are respectively 1/2, 1/3 and 1/4. These normalized
features expose narrow-direction errors that raw-coordinate moments can hide.

### Annular cylinder

Let r be the radius around (1/2,1/2) in the first two coordinates. Accept
**1/4 <= r <= 3/8** and x3..7 in [3/8,5/8]. The score is the maximum of
|r-5/16|/(1/16) and max_j=3..7 |xj-1/2|/(1/8).
The whole annulus lies inside the unit square. Its exact mass is
pi*((3/8)²-(1/4)²)*(1/4)^5 = **5*pi/65536**.

Eight angular octants, starting on the positive x1 direction, each have mass
1/8. Angles are mapped to [0,2*pi); a floating-point wrap exactly at 2*pi is
assigned to the final octant rather than leaving a point unassigned.
All coordinate means are 1/2. E[x1²]=E[x2²]=77/256, and both recorded cross
moments E[x1*x2] and E[x2*x3] equal 1/4.
Record T=(r²-(1/4)²)/((3/8)²-(1/4)²), with E[T]=1/2 and E[T²]=1/3, plus
(1+cos(k*theta))/2 and (1+sin(k*theta))/2 for k=1,2,4, each with mean 1/2.
Angular and radial features assess a curved target without fitting a straight
coordinate box to its truth.

### Unequal disconnected balls

Ball A is centered at (1/4,...,1/4) with radius 1/5. Ball B is centered at
(3/4,...,3/4) with radius 1/4. Both lie inside the cube, and their centers are
sqrt(7)/2 apart, exceeding the sum of their radii.
The score is the minimum of the two center distances divided by their radii.
The volume of the seven-dimensional unit ball is V7=16*pi³/105. Thus
Z = V7*(5^-7+4^-7), and the two region masses are
**pA=16384/94509** and **pB=78125/94509**.

All coordinate means equal 3/4-pA/2. The two recorded second moments equal
pA*(1/16+1/225)+pB*(9/16+1/144), and both recorded cross moments equal
9/16-pA/2. Within each point's unique ball, define uj=(xj-cj)/R and
T=(||x-c||/R)^7. Record E[T]=1/2, E[T²]=1/3, all seven E[(uj+1)/2]=1/2,
and all seven E[uj²]=1/9. These are unconditional mixture features normalized
within the observed component; the separate A/B mass checks assess component loss.

## Sampler and budgets

The unchanged strategy uses two islands, 256 particles each, 14 levels and four
move sweeps per level: **29,184 total pilot attempts** per learned run.
It retains the median upper-quantile threshold rule, tie handling, 50:50 reflected
local/global moves, gap-limited 32-neighbor kernels, all bandwidth choices and
20% uniform defensive component from the earlier frozen design.

Each learned-run root seed spawns pilot stream (0,) and production stream (1,). The pilot
spawns island streams (0,0) and (0,1). Fit one final bounded proposal and draw
exactly **8,192 independent production attempts**, using its entire normalized
mixture density in every importance weight. Pilot points never enter inference.
No oracle sample, component label, target moment or exact mass enters fitting,
adaptation, or the learned proposal density.

The following three learned runs, three full-target controls and three
restricted controls are fixed for each fixture:

| Geometry | Learned seeds | Full-target oracle seeds | Support-hole oracle seeds |
|---|---|---|---|
| Rotated box | 2026092501, 2026092502, 2026092503 | 2026092601, 2026092602, 2026092603 | 2026092701, 2026092702, 2026092703 |
| Annular cylinder | 2026092511, 2026092512, 2026092513 | 2026092611, 2026092612, 2026092613 | 2026092711, 2026092712, 2026092713 |
| Unequal balls | 2026092521, 2026092522, 2026092523 | 2026092621, 2026092622, 2026092623 | 2026092721, 2026092722, 2026092723 |

All 27 root seeds are distinct and unused by the previous studies. Unit tests
use separate seeds, artificial traces and deterministic geometry integrations;
they do not train a sampler on these geometries before the protocol freeze.

## Criteria and analytic controls

Keep the numerical thresholds of the earlier synthetic study:

- Accepted count at least 20 and importance ESS at least 200.
- Maximum normalized weight at most 0.02 and relative mass-estimate MCSE at most 0.10.
- Absolute relative error in exact prior mass at most 0.15.
- Every prespecified region observed and absolute region-mass error at most 0.075.
- Every raw and geometry-specific bounded moment has absolute error at most 0.05.
- Both pilot islands reach the final threshold for a learned run.

All nine learned runs must pass individually; no averaging across seeds or
fixtures rescues a failed check. Every outcome is retained. No early stopping,
seed replacement, budget extension, threshold relaxation or sampler revision is
permitted after observing these outcomes. Any later design change requires a
separate protocol and keeps this study as a completed historical experiment.

Full-target oracle controls sample the known target directly with density 1/Z.
Rotated-box latent coordinates are independent uniforms. Annular squared radius
and angle are uniform, with independent slab coordinates. Ball component
probabilities are their exact masses; an isotropic Gaussian direction and radius
R*U^(1/7) give a uniform ball point. These controls test achievable accuracy
under ideal sampling and do not supply a proposal for learned runs.

Restricted controls use the same transforms but retain only y1>=0 for the
rotated box, the upper angular half of the annulus, or ball B. Their retained
mass fractions are respectively 1/2, 1/2 and pB; each sampled point has density
1/(Z*retained_fraction). Their missing support makes them intentionally invalid
importance proposals for the full target. Their equal weights can pass the
usual ESS/weight screens while known-truth checks detect the omitted regions.
They are distinct from the earlier study's support-preserving defensive-mixture
failure examples. The overall control checks require all nine full-target
oracles to pass and all nine restricted controls to pass usual screens but fail
truth screens. Failure of a control is reported and not tuned away.

## Records and replay

The implementation refuses changes to the frozen sampler source hashes recorded
in the previous v2 study and pins that study's exact JSON SHA256. New geometry,
driver, pilot-audit and protocol files are hashed with every learned archive.
Python, NumPy, SciPy and bit-generator labels are recorded. All learned runs
must share the same specification and runtime labels.
Report runtime labels describe those archived learned runs. Analytic controls
are recomputed in the reader's environment during reconstruction; these labels
do not assert that a reader is running the original numerical environment.

Nine compressed archives retain complete production points, component IDs,
densities, analytic scores and decisions, plus pilot initial/final states,
level histories, ancestors, costs and fitted kernels. Production randomness is
replayed and each score/decision is recomputed. Initial pilot randomness and
endpoint scores are checked; aggregate counts, thresholds and first-level
selection reconcile. Recomputed scores allow 1e-12 relative / 1e-14 absolute
roundoff; a recomputed final pilot score may exceed its archived conditioning
threshold by at most eight floating-point ULPs. Stored scores must remain
exactly at or below that threshold. Intermediate pilot trajectories are not independently
replayed because they are not archived. Oracle controls are reconstructed from
their fixed independent seeds and analytic transforms.

The CLI refuses to overwrite an existing archive. Its --resume mode validates
and retains existing planned runs. A final report requires every planned archive
and rejects extra archive files. --render-only verifies hashes and reconstructs
assessments from the archives instead of trusting cached results. Machine
roundoff may affect the final bits of derived assessments across platforms;
pass/fail decisions, counts, raw records and provenance remain exact.

Run the study once after committing this protocol and its implementation:

    python scripts/proposal_coverage_challenges.py

Resume an interrupted study without replacing completed runs:

    python scripts/proposal_coverage_challenges.py --resume

Reconstruct the report without training:

    python scripts/proposal_coverage_challenges.py --render-only

Passing these selected geometries cannot establish complete biological target
coverage, precise tails, or experimental validity. Failure identifies a sampler
limitation on the named geometry, not proof of the same defect on an unknown
biological acceptance region. Neither outcome changes the historical scores,
produces a pooled posterior, or resolves the independent-assay mapping gap.
