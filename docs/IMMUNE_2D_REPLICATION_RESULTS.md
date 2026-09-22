# The 2D activation contrast persists across twenty seed blocks

All 20 blocks of the [frozen replication study](IMMUNE_2D_REPLICATION_PROTOCOL.md)
completed without retries or replacements. The mean paired SDT − RSL3 contrast
is **0.555834 accumulated activation steps per initial tumor cell**, with a
95% whole-block percentile bootstrap interval of **[0.554223, 0.557439]**.
The contrast is positive under the declared criterion. No eligible opportunity
reached activation 0.5 or 0.9 in any arm in any block.

The [complete report](../analysis/immune-2d-replication.md) and
[derived data](../analysis/immune-2d-replication.json) retain all 60 arm summaries,
including zero and undefined outcomes. These results extend the historical
single-run exposure finding; they do not validate biological immunogenicity,
calibrate a physical clock, or establish a causal mechanism.

## What was fixed before capture

Protocol, plan, Rust mode, numerical readers and tests were committed at
`595269cacc05abf32f484ec8ba2025fd9b317391` before any full-grid block was run.
Each block contains Control, RSL3 and SDT, with the unchanged 500×500 grid,
159,033 initial tumor cells, 180 steps, λ=120 µm oxygen gradient and existing
O2-independent SDT assumption. The opt-in wrapper changes only the root seed.
The biochemical core and shared passive observer are unchanged.

Roots are `42 + block×2^32` for blocks 1–20. An integer source audit verifies
disjoint additive RNG seed-address envelopes without wrapping. The previous
consecutive roots 42–61 reuse addresses: roots 42 and 43 share 249,999 SDT
initialization seed arguments across the grid, including 158,582 tumor-only
arguments. This does not establish the sign or magnitude of outcome correlation;
the historical bootstrap intervals remain descriptive. Disjoint new seed
addresses remove that known aliasing but do not prove PRNG independence.
The entire three-arm block is the resampling unit; within-block RNG reuse is
retained with the original model.

The primary endpoint sums `DAMP/(DAMP+50)` over eligible cell-step opportunities
and divides by the initial tumor census. It is accumulated exposure, not a
probability or a measure per dead cell. The analysis subtracts RSL3 from SDT
within each block and averages the 20 differences. The interval uses 10,000
whole-block resamples and the prespecified Python seed and linear percentiles.
It describes simulation variation conditional on this model and block design,
not parameter, structural or experimental uncertainty.

## Observed populations and null outcomes

The table gives medians and full ranges across 20 blocks. Every RSL3 and SDT
entry below is defined in all 20 blocks.

| Measurement | RSL3 | SDT |
|---|---:|---:|
| Activation steps per initial tumor cell | 0.00579746 [0.00489670, 0.00750893] | 0.561021 [0.554578, 0.567570] |
| Immune kills | 5 [1, 8] | 529 [485, 559] |
| Eligible cell-step opportunities | 686,529.5 [593,085, 922,795] | 2,315,949 [2,283,637, 2,352,610] |
| Mean activation per eligible opportunity | 0.00132592 [0.00127635, 0.00139885] | 0.0385403 [0.0383625, 0.0386943] |
| Maximum eligible DAMP | 2.54785 [1.90124, 2.89772] | 10.4032 [9.93068, 11.3494] |
| Completed releases | 152 [130, 190] | 139,616.5 [139,301, 139,895] |
| Horizon-censored deaths | 3 [0, 6] | 0 [0, 0] |
| Mean LP at completed release | 18.0444 [17.9531, 18.1991] | 19.1186 [19.1067, 19.1322] |

Control produced zero immune kills in every block. Twelve Control blocks had
no eligible opportunities: 1, 2, 5, 8, 9, 10, 13, 14, 16, 17, 18 and 19.
Their primary exposure is zero, while conditional activation means, maxima and
eligibility-normalized rates are undefined. Those populations are not dropped
or imputed. The eight other Control blocks have defined conditional summaries.

The per-block SDT/RSL3 immune-kill ratio has median **106.7** and range
**65–535**, with all 20 ratios defined. This descriptive secondary endpoint
remains sensitive to RSL3's small denominator; a narrow interval for the primary
exposure contrast does not establish a precise immune-kill ratio. No secondary
significance claim was planned.

## What this resolves and what it does not

The absence of half-maximal activation now holds across all 20 declared blocks,
not only the historical root-42 realization. The earlier deep-saturation
explanation is unsupported by these observed eligible populations. Zero
threshold-positive blocks do not establish impossibility at other seeds or
settings. The larger SDT exposure also persists, but these passive observations
do not separate death density, release amount, timing, spatial distribution and
survival. A causal comparison must fix its interventions and matched quantities
before evaluation.

The earlier 2D and 3D archives, the historical consecutive-seed numerical JSON,
and all golden outputs are unchanged. Their different geometries, transport
settings and random streams still prevent a causal 2D-versus-3D interpretation.
The registered P5 criterion remains unchanged and experimentally untested.
[Independent assay mapping and raw replicate metadata](INDEPENDENT_ASSAY_CANDIDATE.md)
remain the main barrier to biological validation.

## Reproduce the analysis

```bash
python3 scripts/immune_2d_replication_report.py
```

The offline reader verifies every compressed payload and frozen source hash,
the frozen source inventory, all event/activation ledgers, seed coverage and
configuration, then recomputes every endpoint and the bootstrap interval. It
requires neither Rust nor a new simulation. The archive contains all 20
observation files, capture logs, the complete default baseline, canonical
observer evidence and the source snapshot; its
[manifest](../analysis/immune-2d-replication/manifest.json) records binary,
Python, Rust, platform and source-commit provenance.

Before new blocks, the built executable reproduced the full historical
33-condition summary and the canonical observer JSON byte for byte. The
published archive was validated before a directory rename; the derived JSON
and Markdown writes remain sequential rather than an atomic pair.

### Reproduce the frozen source tests

Review found that the original source snapshot includes all runtime analysis
dependencies but omits `tests/test_immune_2d_measurement_report.py`, a fixture
module imported by both replication test modules. Offline reconstruction above
is unaffected. The original archive and protocol remain unchanged; future
captures include this fixture, with an isolated test-collection regression guard.

To run the tests bundled with the original snapshot, start in a Git clone of
this repository with Python 3.10 or later and pytest installed. Fetch the PR's
history (the pre-production commit is not on the squash-merged main history),
extract into a new temporary directory and recover the fixture from that exact
commit, checking its SHA-256 before running anything. The temporary Git index
is needed because the capture fixtures enumerate tracked Rust source paths
with `git ls-files`.

```bash
repo_root="$(git rev-parse --show-toplevel)"
reproduction_dir="$(mktemp -d)"
tar -xzf "$repo_root/analysis/immune-2d-replication/sources.tar.gz" -C "$reproduction_dir"
git -C "$repo_root" fetch origin refs/pull/892/head
git -C "$repo_root" show \
  595269cacc05abf32f484ec8ba2025fd9b317391:tests/test_immune_2d_measurement_report.py \
  > "$reproduction_dir/tests/test_immune_2d_measurement_report.py"
(
  set -e
  cd "$reproduction_dir"
  printf '%s  %s\n' \
    '1303aabd6e3afd818ce1d9889a7b75f47429e73dc64a1be0609ed3ba9ea4eef3' \
    'tests/test_immune_2d_measurement_report.py' | shasum -a 256 -c -
  git init -q
  git add .
  python3 -m pytest -q tests/test_immune_seed_blocks.py \
    tests/test_immune_2d_replication_report.py \
    tests/test_immune_2d_replication_capture.py
)
```

Expected for these frozen tests: **124 passed, 1 skipped**. The skipped test
requires the production observations, which are not included in the source
snapshot. The commands run synthetic fixtures only; they do not build Rust or
rerun simulations. Keep the temporary directory if you want to inspect it.
