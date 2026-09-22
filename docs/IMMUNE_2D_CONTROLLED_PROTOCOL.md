# Controlled 2D DAMP sources and recipient availability

This protocol fixes a deterministic experiment on the existing 2D DAMP
transport and activation rules before evaluating its production conditions.
It asks how source count, injected amount, release timing and externally fixed
recipient availability change model activation exposure. Its interventions are
artificial source pulses and recipient masks, not treatments or experimental
assays. It does not validate biological immunogenicity, establish the cause of
the historical SDT:RSL3 immune-kill ratio, or replace P5's experimental criterion.

## Prior knowledge and the implementation freeze

The [canonical passive study](IMMUNE_2D_MEASUREMENT_PROTOCOL.md) and
[twenty-block replication](IMMUNE_2D_REPLICATION_PROTOCOL.md) are already known.
They found a larger SDT activation exposure without an eligible opportunity
reaching activation 0.5 in the declared runs. Those selected living populations
do not separate source number, release per source, release timing or recipient
survival. This is an intervention study motivated by those observations, not an
unseen test of their original interpretation.

The machine-readable contract is
[`scripts/immune_2d_controlled_plan.json`](../scripts/immune_2d_controlled_plan.json).
Freeze this protocol, that complete plan, the Rust implementation, the offline
reader and their tests in a committed checkout before evaluating any complete
production condition. Synthetic and reduced-grid software checks may precede
the freeze. Fix defects found before capture, declare any resulting protocol
changes and commit a new freeze before evaluation. A methods change after
production evaluation requires a separately identified design and freeze; never
tune masks, budgets, timing, denominators or reporting to observed outcomes.

## Fixed model and deterministic source masks

Retain a 500 by 500 grid, 20 micrometre site spacing and 180 steps numbered
0 through 179. The spacing is the existing model geometry, not a new spatial
calibration. DAMP starts at zero everywhere. Retain Moore-eight transport with
diffusion fraction 0.08 per neighbor, source-spread cutoff 0.001, clearance 0.03
per step, eligibility threshold 0.01 and activation Kd 50. Transport covers the
whole grid, including sites outside the source/recipient region. It does not
wrap at boundaries. Sources below the spread cutoff retain their mass before
clearance. Preserve the existing neighbor and accumulation order.

Use the central square `218 <= row < 282`, `218 <= col < 282`. This 64 by 64
region lies within the canonical tumor circle. Define local coordinates
`u = row - 218`, `v = col - 218`; tile coordinates `a = u // 8`, `b = v // 8`;
and within-tile coordinates `p = u % 8`, `q = v % 8`. The region contains
64 tiles, each with 64 sites. The masks use these formulas and no random draws:

- The dense source mask contains sites with `p in {2, 5}` and `q in {2, 5}`.
  There are four sites per tile and 256 sources in total.
- The sparse source mask contains one site per tile: `p = 2` when `a` is even
  and `p = 5` when `a` is odd; `q = 2` when `b` is even and `q = 5` when `b`
  is odd. There are 64 sources, all members of the dense mask.
- Reserve all 256 potential source sites in every condition. The reference
  recipient mask contains every other site in the central square: 3,840 sites.
  In particular, unused dense-source sites do not become recipients in sparse
  conditions.
- The reduced recipient mask contains reference recipients with `u % 2 == 0`
  and `v % 2 == 0`. Each tile contributes 15 such sites, giving 960 recipients.

Archive every source and recipient mask as sorted row-major indices
`row * 500 + col`. Validate these identities independently from the recorded
counts. The sparse/dense source densities are 64/4,096 and 256/4,096 sources per
site in this declared region. A density contrast also changes the declared
spatial source pattern; it does not establish an effect for arbitrary layouts.
The recipient contrast similarly applies to these particular fixed masks.

## Conditions, release phases and availability

Cross the four factors below in the listed order, with recipient count varying
fastest. The plan lists every condition and its identity explicitly.

| Factor | First level | Second level |
|---|---:|---:|
| Source count | 64 | 256 |
| Total injected DAMP | 1,280 | 5,120 |
| Release step | 30 | 60 |
| Available recipients | 3,840 | 960 |

Each source releases the same amount once. Divide the total budget by source
count: the four source-count/budget combinations give per-source amounts of
20, 80, 5 and 20 DAMP model units, respectively. These round, declared budgets
are interventions, not fitted release estimates or assay measurements. A release
step specifies actual injection; do not add a death-to-release grace period.

Append `zero_r3840` and `zero_r960` after the sixteen positive-budget conditions.
Both zero controls retain the 64-source mask and release step 60, with total
and per-source release equal to zero. The full study therefore has 18 conditions.
Retaining a source mask in a zero control does not create a positive release.

At each step, inject any scheduled source pulse, diffuse DAMP, apply clearance,
and then observe recipients. Observe immune opportunities only at steps 60
through 179 inclusive and only where local DAMP is at least 0.01. The timing
contrast compares a pulse thirty steps before the immune window with a pulse
at its onset, over the same 120-step observation window. It changes the age and
transport history of the field at observation; it is not a claim about a
physical immune delay.

Available recipient identities remain fixed through every step. Recipients
neither die nor remove DAMP; source sites remain excluded even after release.
There is no biochemical progression, ferroptotic death, immune-kill draw,
endogenous survival selection, stromal effect, pH effect or treatment arm in
this experiment. Each field is determined by sources, release amount and timing.
Changing only recipient availability must leave every field step identical.
All outcomes are deterministic under this model and implementation. Repeated
execution is a reproducibility check, not an independent replicate.

## Endpoints and the frozen contrasts

For an available recipient `i` at step `t`, let `D(i,t)` be its DAMP after
transport and clearance. Include it only if `60 <= t <= 179` and
`D(i,t) >= 0.01`. Its activation is `D(i,t) / (D(i,t) + 50)`.

The primary endpoint, `activation_sum_per_reference_recipient`, is the sum of
these activations divided by **3,840 in every condition**. This fixed denominator
includes unavailable reference recipients in reduced-mask conditions, so their
absence remains part of the availability effect. The endpoint is zero when
there are no eligible opportunities. Its units are accumulated activation steps
per reference recipient, not a probability, kill count or per-source potency.

Retain per-step and per-recipient eligible-opportunity counts, summed DAMP,
summed activation, maximum eligible DAMP, and counts at DAMP >=50 and >=450.
Report unique eligible recipients separately from repeated recipient-step
opportunities. Report the total activation sum; its primary normalization;
its normalization by the condition's available-recipient census; mean activation
and mean local DAMP per eligible opportunity; maximum eligible DAMP; and each
threshold count and fraction of eligible opportunities. The two threshold
levels correspond to activation at least 0.5 and 0.9. Every conditional mean,
maximum or fraction is null when its population is empty; counts and sums are
zero. Available-recipient-normalized exposure uses 3,840 or 960, as declared,
and remains distinct from the primary fixed-denominator endpoint.

The plan fixes 36 directed contrasts of the primary endpoint. Each is its
listed `left` endpoint minus its listed `right` endpoint:

| Contrast kind | Left minus right | Held fixed | Number |
|---|---|---|---:|
| Equal total mass | 256 sources minus 64 sources | Total budget, timing, recipients | 8 |
| Release amount | Budget 5,120 minus 1,280 | Source mask, timing, recipients | 8 |
| Release timing | Step 60 minus step 30 | Source mask, total budget, recipients | 8 |
| Recipient availability | 960 recipients minus 3,840 | Source mask, total budget, timing | 8 |
| Equal per-source amount | 256 sources at budget 5,120 minus 64 at budget 1,280 | Per-source amount 20, timing, recipients | 4 |

At fixed total mass, changing source count necessarily changes per-source
release. Those contrasts measure partitioning of the same injection budget
among the declared source positions. Conversely, the equal-per-source contrast
changes both source count and total mass. Do not label either comparison as
holding all three quantities fixed. The crossed design and complementary
contrasts make their different interventions explicit.

Publish every condition and contrast without selecting by sign or magnitude.
There is no success threshold, hypothesis-test p-value, confidence interval,
bootstrap, seed resampling or claim of independent uncertainty. Do not average
the different source or recipient sites as if they were replicate experiments.
The zero controls and ledger checks test software behavior, not biological
success. An invalid or absent condition is not a zero result.

## Reconciliation and what the archive can reconstruct

Validate finite nonnegative values, complete ordered condition coverage, exact
factor/configuration identity, exact masks and all injection totals. Integer
counts and mask identities must agree exactly. Floating-point reconciliation
uses relative tolerance `1e-10` and absolute tolerance `1e-8`, accepting
`abs(actual - expected) <= 1e-8 + 1e-10 * abs(expected)`. Apply these tolerances
to mass accounting, per-step/per-recipient sums and recomputed endpoints. This
does not relax any exact identity, condition coverage or threshold-count check.

From zero initial field, each step's retained mass must equal
`(previous_retained_mass + injected_mass) * 0.97` within that tolerance.
Diffusion transfers mass within the whole grid; recipients do not consume it.
Record injections, per-step whole-field mass and a deterministic field
fingerprint. Fingerprints check identical field trajectories across recipient
masks; they are not a cryptographic integrity check or proof that transport is
correct. Archive file integrity uses SHA-256 separately.

Per-step and per-recipient counts and sums must reconcile. The ordered step
ledger must contain every step and no eligible opportunities before step 60;
each available recipient has at most 120 opportunities. Individual recipients'
first and last eligible steps are not archived, so these aggregates cannot
reconstruct their eligibility chronology. Maxima must agree across aggregations.
Threshold counts must nest: count at DAMP >=450 <= count at DAMP >=50 <= eligible
opportunity count. Enforce activation and DAMP bounds implied by eligibility.
The reduced recipient identities must be a subset of the reference recipients;
their retained per-recipient observations must match the corresponding full-mask
observations when the source condition is identical. Both zero controls must
have zero injections, field mass and eligible observations.

Archive per-step and per-recipient aggregates, not every local-DAMP value at
every recipient step. Offline reconstruction verifies those ledgers and
recomputes every declared endpoint and contrast. It cannot reconstruct the
unarchived field trajectory or independently establish the transport history
from aggregate conservation alone. A replay of the frozen binary/source is the
transport reproduction path. This limitation must remain explicit in the report.

## Executable isolation, capture and publication

Use the standalone `sim-tme --immune-controlled` mode, with no other arguments.
Reject every `FERRO_*` environment override. It writes
`output/tme/immune_controlled.json` in its private run directory. Keep its output,
schema and reports separate from the canonical observer and replication modes.
The production kernel and controlled mode share the same extracted 2D transport
routine; the extraction must preserve the former arithmetic and iteration order.

Before evaluating the new complete conditions, build the clean committed freeze
with the pinned Rust 1.96.0 toolchain, locate the Cargo-reported executable and
copy it into private execution storage. Record the full source snapshot, source
commit, binary hash, toolchain, Python implementation/version and platform.
Require the complete historical 33-condition summary to match its frozen golden
and the canonical observer output to match the original 2D archive byte for
byte. These are compatibility checks using already-seen observations. A mismatch
stops capture; never change a golden to accept it.

Capture and reconstruct through
`scripts/immune_2d_controlled_report.py`, following its `--capture` contract and
offline default mode. Use a new archive directory, refuse an existing archive
or dirty checkout, and retain reports outside the archive. Preserve all old
measurement and replication archives. Execute all 18 conditions in the declared
order without outcome-dependent retries, replacements, dropping conditions or
early stopping. An operational rerun with unchanged methods must be identified.
A simulation or reconciliation failure prevents publication of a complete study;
retain its logs and disclose the failure.

Archive compressed observations, the baseline and canonical parity evidence,
the complete source snapshot including imported test fixtures, capture logs and
a manifest of hashes and provenance. Validate the complete archive in a sibling
staging directory before publishing it. Public JSON and Markdown are derived
outputs and must reconstruct from the archive without Rust execution. Their
writes are sequential and are not an atomic pair against disk failure.

## Limits fixed before evaluation

DAMP amounts are uncalibrated model units, with total mass meaning the sum of
the lattice field. They are not concentrations or molecular counts measured by
an assay. Source density is expressed per model lattice site in the declared
region. Activation is dimensionless; accumulated activation carries a model-step
factor. There is no new mapping of steps to hours, DAMP to a named assay, or
activation to measured dendritic-cell maturation.

This study can attribute changes in its deterministic response to its declared
interventions within these fixed rules and masks. It does not estimate the
contribution of each mechanism to the original treatment comparison, establish
an effect for other layouts, model endogenous recipient survival, or validate a
biological mechanism. Existing P5 text and its experimental falsification
threshold remain unchanged. Independent assay mapping and raw biological
replicates remain necessary for biological validation.
