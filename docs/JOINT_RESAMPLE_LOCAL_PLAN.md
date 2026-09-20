# Revised joint proposal: local exploration and separated neighborhoods

This is a separate protocol after the first bounded design's
[failed prerequisite](JOINT_RESAMPLE_FAILURE.md). The old implementation and
study remain unchanged. This file and the revised numerical implementation must
be committed before the prospective synthetic runs below. No biological run has
been performed under either bounded design when this plan is written.

## Development record and selection

The failed study exposed broad minority-region kernels caused by including
distant centers among 32 unique neighbors. A declared development comparison
used seeds **2026092301, 2026092302 and 2026092303**, the separated-box fixture,
and the original budgets. Four variants compared 32 versus 8 neighbors and pure
independence moves versus a 50:50 local/global mixture. They passed respectively
0/3, 0/3, 2/3 and 0/3 runs. A second declared comparison added a gap-limited
32-neighbor rule to the mixed design, and compared it with two local sweeps
followed by a refit and two global sweeps. Both passed 3/3 development runs.

Select **mixed moves with gap-limited 32-neighbor kernels**: it preserves the
original once-per-level fit and passed the declared development checks. These
six configurations and all 18 results are retained in the
[development archive](../analysis/calibration/proposal-development-20260919/README.md).
The prospective implementation was checked against all three selected
development runs, reproducing final clouds, proposals and independent production
assessments exactly. Development results are not prospective evidence.

## Scientific target and fixed budgets

Keep the scientific target in the [original bounded protocol](JOINT_RESAMPLE_PLAN.md):
the same seven uniform priors, supported CTRPv2 cohorts and grids, reference
vector, 2,000 cells, simulation seed 42, and full-precision final tolerance
**0.17447979628221513**. ML210 remains excluded from adaptation and acceptance.
It is a same-screen compound holdout with overlapping cell lines, not an
independent raw-assay validation set.

Each root seed spawns a pilot stream `(0,)` and independent production stream
`(1,)`; the pilot stream spawns two island streams `(0, 0)` and `(0, 1)`.
Each island starts with 256 uniform draws and runs 14 levels, four sweeps per
level. Thresholds and resampling follow the original protocol: upper empirical
50% quantile, monotone restriction, no threshold below epsilon, retain all ties,
and no resampling when every particle already qualifies. There are exactly
**29,184 pilot attempts** and **8,192 independent production attempts** per run.
No stopping on a favorable sample count or extending a failed run is allowed.

## Revised proposal and moves

The final and each level's fitted proposal retain 20% uniform cube density and
80% equally weighted product truncated-normal kernels, one per particle.
The full normalized mixture density is evaluated for every importance weight
and global Hastings ratio. Particle duplicates retain their mixture weight.

For bandwidths, deduplicate centers, standardize coordinates by their global
population standard deviation with floor 0.005, and order centers by Euclidean
distance with stable tie ordering. Consider at most 32 neighbors including self.
Ignore the zero self-distance. Stop before the first adjacent positive-distance
jump exceeding a factor of **3**. Use 1.5 times the coordinate population standard
deviation of the remaining neighborhood, clipped to [0.005, 0.5]. A singleton
neighborhood uses width 0.04 before those bounds. This heuristic can miss
singleton regions and can mistake a sparse neighborhood for a gap; it is not a
certified mode detector.

Fit once after resampling at each level and freeze the proposal for four sweeps.
Compute one shared diagonal local step vector as 0.5 times the median of all
particle-kernel widths, including duplicates, clipped to [0.0025, 0.25]. Freeze
this vector for the same four sweeps. For each particle and sweep independently,
choose local versus global movement with probability 1/2, independent of state.

Local candidates add independent zero-mean Gaussian increments with the shared
step vector and repeatedly reflect each coordinate into [0, 1] using modulo-2
folding. Accept exactly when the candidate meets the current threshold. The
coordinate reflection kernel is symmetric for these fixed independent Gaussian
increments; per-particle scales, correlated reflection, clipping, or truncating
the proposal would not justify this acceptance rule. Global candidates use the
original independence-MH rule, including `min(1, q(old)/q(new))`. A failed global
MH coin skips simulation; a strict erastin lower bound may also reject a point.

The reproducible sweep draws global candidates for all particles, then local
choice uniforms, then local Gaussian increments, then MH uniforms for all
particles. Some candidates/uniforms are unused; they cost no model calls. Every
local candidate is evaluated. Archive local/global attempt and acceptance
counters alongside the common evaluation and dose-call counters.

These finite adaptive clouds train a proposal. Their particles, clones and
ancestry counts are not independent posterior observations. Inference uses only
the fresh, independent production stream from the final frozen proposal.

## Prospective prerequisite and biological evaluation

Run both unchanged analytic fixtures at **2026092401, 2026092402 and 2026092403**.
Use the original 8,192-attempt budget, analytic truths, usual weight screens,
mass/mode/moment error limits and required island threshold attainment. Reuse
and clearly label all three original negative controls. Preserve the first
study's SHA256 and all numerical dependency hashes. The revised validator writes
separate `proposal-synthetic-validation-v2.json` and `.md` artifacts.

Every prospective run must pass every declared check before biological
production. Fresh randomness on development fixtures is a reproducibility and
precision check, not unseen validation of general region discovery. If the
prerequisite fails, retain the failure and stop this design before biology.

Only on a prerequisite pass, use biological root seeds **2026091911,
2026091912 and 2026091913**, which have not yet been used. Keep every per-run and
across-run screen from the original bounded protocol, including all islands
reaching epsilon. No criterion, seed or budget may change after viewing these
results. Retain all runs and evaluated production curves and rebuild summaries
offline. Publish run-wise diagnostics, with no pooled posterior even on a pass.

## Interpretation

Passing these checks does not rule out an unseen region, establish precise
tails, validate an assay mapping, or resolve independent biological validation.
The defensive prior component preserves support but provides little protection
against extremely rare undiscovered regions at this budget. Local movement may
fail to recover a region after all its particles have disappeared.

Method background remains the references in the original protocol; the local
step is the usual symmetric random-walk Metropolis construction, described in
[Oxford's advanced simulation notes](https://www.stats.ox.ac.uk/~deligian/pdf/sc5/notes/notes6.pdf).
All numerical choices and operational screens here are repository design
choices, not guarantees inherited from those references.
