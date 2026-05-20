# sim-tme-3d

3D spheroid tumor microenvironment simulation. Capstone binary for the spheroid-validation series (#185–#197) — the first consumer of all five library primitives landed in `ferroptosis-core` v0.7.0–v0.11.0.

## What it does

Runs a matrix of 24 conditions on a 60³ spheroid (~82.5k tumor cells, ~540 µm radius), integrating:

- **3D energy physics** (#186) via `physics::local_ros_multiplier_3d`
- **3D radial O₂ gradient** (#187) via `oxygen::radial_o2_field`
- **3D radial pH gradient** (#190) via `ph::radial_ph_field` + `iron_multiplier_from_ph` + `ion_trap_factor_from_ph`
- **3D CAF-shielded boundary detection** (#189) via `stromal::stromal_adjacency_mask`
- **3D spatial DAMP diffusion + activation** (#188) via `immune_spatial::diffuse_damp_3d_step` + `dc_activation` + `immune_kill_probability`

Emits a JSON summary that the Python comparison script pairs with `sim-tme`'s existing 2D output to answer the four manuscript-keystone questions from issue #195.

## ⚠️ Scale mismatch with sim-tme

| | sim-tme (2D) | sim-tme-3d (3D) |
|---|---|---|
| Grid | 500 × 500 | 60³ |
| Total cells | 250k | 216k |
| Tumor cells | ~159k | ~82.5k |
| Tumor radius | ~4500 µm | ~540 µm |
| Tumor diameter | ~9 mm (large in-vivo) | ~1.1 mm (upper end of spheroids) |

3D fundamentally can't match 2D's biological scale at any feasible grid size (500³ = 125M cells × ~170B ≈ 21 GB). **Compare via RATIOS, not absolute counts.** The Python comparison script does this automatically.

## Why λ sweep is `[80, 100, 120]` (skip 150)

The 3D hypoxic-zone threshold is `3λ`. At λ=150 µm, the threshold is 450 µm — but the 60³ grid's tumor radius is only 540 µm. The hypoxic zone would be ≤ 1 cell — statistically meaningless. The three smaller λs give meaningful hypoxic shells.

## Why `damp_diffusion_fraction = 0.025` (NOT sim-tme's 0.08)

Sim-tme's 2D default `0.08` is **unsafe in 3D**: with up to 26 Moore neighbors, `0.08 × 26 = 2.08 > 1` would mass-destroy the DAMP field. The library function `immune_spatial::diffuse_damp_3d_step` enforces the stability invariant with `assert!` — release-mode panic if violated. We use `0.025` (matches 2D's per-step total diffusion of ~64%).

## Usage

```bash
cd simulations
cargo run --release -p sim-tme-3d
# → output/tme-3d/summary.json (24 conditions)

# Generate the 2D vs 3D comparison. The comparison script reads BOTH
# output/tme/tme_summary.json (sim-tme, 2D) and output/tme-3d/summary.json
# (sim-tme-3d, 3D). Neither output is tracked in git — only .gitkeep is.
# On a clean checkout, run sim-tme first (10-30 min, 500x500 grid):
#   cargo run --release -p sim-tme
python3 ../scripts/generate_3d_comparison_table.py
# → output/tme-3d/comparison_2d_vs_3d.csv
# → output/tme-3d/key_questions.txt
```

The 3D 24-condition run takes ~15-30 seconds on 8 cores (rayon condition-level parallelism). The 2D prerequisite (`sim-tme`) is much heavier — 10-30 minutes on the same hardware.

## Condition matrix

| Category | Conditions | Description |
|---|---|---|
| Baseline | 3 | Uniform O₂, no toggles, {Control, RSL3, SDT} |
| O₂ sweep | 9 | λ ∈ {80, 100, 120} × {Control, RSL3, SDT} |
| Immune | 3 | Immune coupling on, λ=120, {Control, RSL3, SDT} |
| Stromal | 3 | CAF shielding on, λ=120, {Control, RSL3, SDT} |
| pH | 3 | pH gradient on, λ=120, {Control, RSL3, SDT} |
| Combined | 3 | Immune + stromal + pH on, λ=120, {Control, RSL3, SDT} |

Total: 24 conditions. Smaller than sim-tme's ~45 (no anti-PD-1, no O₂ cycling for v1 — see follow-ups below).

## Output

`output/tme-3d/summary.json` — per-condition kill rates + metadata.

Schema mirrors `sim-tme`'s `tme_summary.json` (both wrapped in an envelope object since #224) so the comparison script can read both.

### Schema versioning

Both `output/tme/tme_summary.json` (sim-tme) and `output/tme-3d/summary.json` (this binary) emit a `schema_version: u32` field at the top level. **The current schema version is `1`.**

```json
{
  "schema_version": 1,
  "grid_dim": 60,
  ...
  "conditions": [ /* ConditionResult[] */ ]
}
```

**Bump the version when the shape changes.** Adding optional fields under `conditions[]` is non-breaking and does not require a bump. Renaming or removing top-level keys, changing a field's semantics, or reshaping `conditions[]` does require a bump in both binaries together. The Python comparison script (`scripts/generate_3d_comparison_table.py`) asserts both files have the **same** `schema_version` equal to its `EXPECTED_SCHEMA_VERSION` constant — schema drift between the two binaries fails loudly there instead of silently producing `None`-filled rows.

## Tests

```bash
cargo test --release -p sim-tme-3d
```

Three smoke tests:
1. `condition_matrix_is_non_empty` — matrix sanity
2. `single_condition_runs_end_to_end` — full orchestration on baseline Control
3. `same_seed_same_output` — determinism

The library primitives (`physics`, `oxygen`, `ph`, `stromal`, `immune_spatial`) are exhaustively tested in `ferroptosis-core`'s 160+ unit tests. This binary tests orchestration, not the math.

## Manuscript-keystone questions (issue #195)

After running both `sim-tme` and `sim-tme-3d` and generating the comparison table. Each bullet states the pre-run **hypothesis** from issue #195 and the **observed** result from the canonical 60³ × 180-step run (full details in `simulations/calibration/3d_validation_report.md`).

1. **Does the hypoxia RSL3 collapse hold in 3D?**
   - Hypothesis (#195): yes, possibly stronger.
   - **Observed**: yes qualitatively (within-zone collapse 98.4% at λ=120). Like-for-like, **2D collapses more completely** on both metrics — within-zone 2D 0.0064 < 3D 0.016; overall 2D 0.028 < 3D 0.222. The "possibly stronger" hypothesis was wrong; 3D collapse is robust but smaller magnitude than 2D. See `key_questions.txt` Q1.

2. **Does the immune 104:1 ratio hold in 3D?**
   - Hypothesis (#195): unknown — DAMP density may decrease in 3D volume.
   - **Observed**: direction holds, magnitude much smaller — SDT/RSL3 immune-kills = 4.0× in 3D vs 104.2× in 2D. The ~2× tumor-cell gap (82.5 k 3D vs 159 k 2D) is too small to fully explain the ~25× ratio gap; volumetric DAMP dilution and per-cell activation density also contribute. See Q2.

3. **Does stromal shielding have MORE impact in 3D?**
   - Hypothesis (#195): yes — ~1.5× boundary fraction per #189 cross-geometry test.
   - **Observed**: no — per-cell shielding is essentially geometry-independent. Boundary shielding = 51.5% (3D) vs 50.0% (2D). The cubic-vs-quadratic scaling from #189 affects HOW MANY cells are shielded, not the per-cell magnitude. See Q3.

4. **Does pH ion trapping produce similar RSL3 reduction in 3D?**
   - Hypothesis (#195): similar — same chemistry.
   - **Observed**: yes — 46.1% kill reduction in 3D vs 54.2% in 2D, within noise. See Q4.

## Follow-ups deferred to subsequent PRs

- **Lift `PhConfig` / `StromalConfig` / `ImmuneConfig` to `ferroptosis-core::params`** — currently duplicated from sim-tme. Per the consolidated cleanup checklist on issue #195.
- **O₂ cycling** (square-wave λ alternation) — sim-tme has it, sim-tme-3d skipped for v1 scope.
- **Anti-PD-1 sweep** — included in sim-tme; skipped here for v1.
- **3D volumetric visualization** — see #193.
- **Larger grids** (80³+) — gated on #194 perf work.
- **Empirical pimonidazole validation** — see #196.
