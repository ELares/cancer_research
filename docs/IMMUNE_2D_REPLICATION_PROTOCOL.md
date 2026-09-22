# Replicate the 2D immune activation-window measurements

This protocol fixes a new simulation study before its production observations.
The question is whether the previously measured activation-window exposure
contrast and lack of threshold saturation persist across whole-run seed blocks
of the unchanged 2D model. It does not validate an assay, isolate a causal
mechanism, establish a physical clock, or replace P5's experimental criterion.

## Prior knowledge and seed reuse

The canonical passive measurement is already known: 521 SDT versus five RSL3
immune kills, and no eligible opportunity reaching activation 0.5. Its full
archive and the separate 3D archive remain unchanged. The earlier count-only
study used consecutive seeds 42–61 and reported a median defined kill ratio
near 107, with one undefined ratio. These outcomes informed this question;
this is a declared extension, not an unseen test of the original headline.

The 2D engine constructs each cell's RNG seed by addition. For example,
SDT initialization uses `master_seed + treatment_offset + cell_index`. Thus seed
42 at cell index `i+1` and seed 43 at index `i` initialize the same RNG state
when both positions are tumor cells. Different master seeds do not establish
independent trajectories. The initial seed-address overlap is exactly
countable from the tumor circle; later realized draws depend on survival and
eligibility. Reuse alone does not determine the magnitude or sign of outcome
correlation or interval error. Preserve the earlier numbers as descriptive
dispersion; its ordinary bootstrap intervals do not establish independent-run
uncertainty for this design.

Freeze 20 new roots, `42 + block_id * 2**32`, for block IDs 1–20. A block is
one complete three-arm simulation realization. Geometry and initial cell
properties are shared across its arms; treatment-specific runtime offsets
remain 0, 10,000,000 and 20,000,000. Existing within-block seed reuse across
arms and time is preserved. Every permitted RNG seed address in each block
must fit in its recorded inclusive namespace envelope, and the 20 envelopes
must be disjoint without unsigned-integer wrapping. This removes the known
cross-block additive aliases; it is not a proof of PRNG statistical independence.

## Frozen conditions and observations

Run Control, RSL3 and SDT, in that order, with the existing 500 × 500 grid,
20 µm spacing, circular tumor radius 225 cells, 180 steps, oxygen penetration
120 µm, immune-on baseline block, stromal and pH effects off, and the existing
oxygen-independent SDT setting. Retain all biochemical, spatial and immune
defaults from `immune_2d_measurement_config_v1.json`. Only the root seed changes.
Do not include anti-PD-1 conditions or change the 3D simulation.

The opt-in executable mode is `sim-tme --immune-replicate BLOCK_ID`. It rejects
unknown block IDs, extra arguments and every `FERRO_*` environment override.
It writes `output/tme/immune_replicate.json` in a private run directory. The
canonical `--immune-measurements` mode and default production matrix retain
their existing contracts and outputs.

Use the existing passive observer without changing its event phases, eligibility
rule, DAMP response or random draws. Preserve every block's death, release,
censoring, eligible-cell, per-step and immune-kill records. Reconcile every
ledger with the existing event/activation validator, including tumor identity,
cause-of-death disjointness, phase bounds, threshold nesting, population sums
and DAMP conservation. An unavailable or failed block is not a zero observation.

## Endpoint and analysis fixed before capture

The primary endpoint in each arm is **total activation over eligible cell-step
opportunities divided by the initial tumor count**. An opportunity's activation
is the unchanged `DAMP / (DAMP + 50)`. The numerator sums these activations;
the denominator is the initial tumor census, not the selected eligible population.
The value is zero when no opportunities occur and has units of accumulated
activation steps per initial tumor cell. It is not a probability, per-death
immunogenicity, or a biological assay measurement.

For each block, subtract RSL3's endpoint from SDT's. Report all 20 paired
differences and their arithmetic mean. Form a percentile bootstrap interval
for that mean using 10,000 resamples of the 20 **whole blocks**, with Python's
`random.Random(20260922)`. Each resample selects 20 block indices with
replacement. Use linear interpolation at the 2.5th and 97.5th percentiles of
the sorted bootstrap means. Never resample cells, events, arms, or opportunities
as independent units. The interval describes simulation variation under this
fixed design and its pseudorandom-block assumption, not parameter, structural,
or experimental uncertainty.

An interval entirely above zero supports a positive mean contrast under this
design; one entirely below zero supports a negative mean contrast; an interval
including zero is unresolved by this criterion. Publish whichever occurs.
Do not choose a different endpoint, interval method, seed set or stopping rule
after observing the results.

For every arm and block, also report immune kills; unique eligible cells and
opportunity counts; conditional mean activation; maximum eligible DAMP; counts
and shares with DAMP ≥50 and ≥450; ferroptotic deaths, completed releases and
horizon-censored deaths; and the corresponding LP/DAMP summaries. The two
thresholds correspond to activation at least 0.5 and 0.9. Report how many blocks
reach each threshold and retain every positive occurrence. Zero occurrences
in these 20 runs do not establish impossibility in other seeds or settings.

Report secondary summaries as medians and full ranges across blocks with an
explicit number defined. Keep undefined conditional means and zero-denominator
kill ratios as null; do not replace them with zero, infinity, or an imputed
denominator. Per-block SDT/RSL3 kill ratios are descriptive and exclude
undefined ratios only with their count and identities disclosed. Secondary
endpoints receive no confirmatory significance claims or replacement P5 threshold.

## Freeze, capture, and publication

1. Test the new mode on reduced grids for all-arm observer passivity and seed
   identity. Commit this protocol, code, plan and tests before producing any
   of the 20 full-grid replicate blocks.
2. Build the pinned Rust 1.96.0 executable from that clean committed checkout
   using Cargo's reported artifact path, copy it to private execution storage,
   and record its hash, toolchain and complete source snapshot.
3. Before the new blocks, require the full default 33-condition summary to
   match the existing golden SHA and the canonical observer output to match
   the existing 2D archive byte for byte. Never update a golden to admit a
   mismatch. These are compatibility checks using already-seen outcomes.
4. Run all 20 blocks in their declared order without outcome-dependent retries,
   dropping blocks, seed replacement, or tuning. Record runtime and provenance.
   A simulation or reconciliation failure stops publication of a complete study;
   preserve logs and report the failure. An unchanged rerun after an operational
   interruption must be identified; a methods change requires a new freeze.
5. Publish the complete compressed observations, source archive, baseline
   parity evidence, manifest, derived JSON and Markdown. Use a separate archive
   directory that must not already exist. Prepare and validate in a sibling
   staging directory before publication. Public JSON/Markdown writes are not
   an atomic pair against disk failure.
6. Provide offline reconstruction from the archive, with checksums, frozen
   configuration/seed coverage, event reconciliation and recomputation of every
   endpoint and bootstrap result. No Rust execution is needed for reconstruction.

The capture entry point is `python3 scripts/immune_2d_replication_report.py --capture`.
It refuses an existing archive or a dirty checkout. The same command without
`--capture` reconstructs the published reports offline. Record the Python
implementation/version as well as Rust and platform provenance. Source
validation requires the critical protocol, reader and test files and the frozen
86-file Rust/Cargo inventory independently of the manifest's declared keys.

The study preserves historical numerical artifacts. Its interpretation must
name the model, seed unit and selected populations, including null outcomes.
