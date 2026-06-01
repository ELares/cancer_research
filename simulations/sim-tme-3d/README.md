# sim-tme-3d

3D spheroid tumor microenvironment simulation. Capstone binary for the spheroid-validation series (#185–#197) — the first consumer of all five library primitives landed in `ferroptosis-core` v0.7.0–v0.11.0.

## What it does

Runs a matrix of 24 conditions on a 60³ spheroid (~82.5k tumor cells, ~540 µm radius), integrating:

- **3D energy physics** (#186) via `physics::local_ros_multiplier_3d`
- **3D radial O₂ gradient** (#187) via `oxygen::radial_o2_field`
- **3D radial pH gradient** (#190) via `ph::radial_ph_field` + `iron_multiplier_from_ph` + `ion_trap_factor_from_ph`
- **3D CAF-shielded boundary detection** (#189) via `stromal::stromal_adjacency_mask_3d`
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

The spheroid is also still **in-vitro** scale (~1.1 mm), where most cells are within drug/O2 penetration depth. [Patient-scale slab mode (#240)](#patient-scale-slab-240) addresses the *other* end of the scale gap — a slab at depth in a 5 mm+ virtual tumor, where penetration collapses and efficacy drops sharply — so the same engine can report numbers at both the spheroid and the patient scale.

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

## Trajectory snapshot (`--snapshot[=NAME]`, #193)

For visualization, pass `--snapshot` to run **one** condition (RSL3 + immune_on + stromal_on + ph_on at λ=120 µm — the most visually rich cell of the matrix) with per-step state capture:

**Presets** (`--snapshot=NAME`; bare `--snapshot` resolves to `combined`,
so the original UX is preserved). An unknown name prints the list and
exits 2; add a preset via a one-entry change to the `SNAPSHOTS` registry
in `main.rs`. All presets are RSL3 at λ=120 µm; the toggles vary, and the
output files use the same names regardless of preset (a rerun overwrites).

| Name | Condition |
|---|---|
| `combined` (default) | RSL3 + immune_on + stromal_on + ph_on. All TME protections active, the most visually rich cell of the matrix. |
| `bare` | RSL3 with none of the TME protections. The death front sweeps the spheroid more visibly (~3x higher kill rate than `combined`). |

```bash
cd simulations
cargo run --release -p sim-tme-3d -- --snapshot          # = --snapshot=combined
cargo run --release -p sim-tme-3d -- --snapshot=bare     # unprotected baseline (overwrites the files)
# → output/tme-3d/trajectory_dead.npy   (180 × 60 × 60 × 60, u8)
# → output/tme-3d/trajectory_damp.npy   (180 × 60 × 60 × 60, f32)
# → output/tme-3d/trajectory_lp.npy     (180 × 60 × 60 × 60, f32)
# → output/tme-3d/trajectory_meta.json  (schema + condition descriptor)

# Render an animated axial mid-slice GIF (+ MP4 if ffmpeg available)
python3 ../scripts/render_tme_3d_trajectory.py
# → output/tme-3d/trajectory_axial.gif  (~4 MB, 180 frames @ 15 fps = 12s)
# → output/tme-3d/trajectory_axial.mp4  (if ffmpeg on PATH)
```

The default 24-condition matrix path (no `--snapshot` flag) is **byte-identical** to before #193 — `summary.json` hash is unchanged. Only the snapshot path touches the new trajectory capture code.

**On-disk size**: ~333 MB total (37 MB dead + 148 MB damp + 148 MB lp). `output/` is git-ignored; the trajectory is meant to be regenerated locally.

**Schema versioning**: `trajectory_meta.json` carries its own `schema_version: u32` (currently `1`), separate from `summary.json`'s schema. The Python renderer hard-asserts the version to fail loudly on drift.

### Snapshot presets (`--snapshot=NAME`)

| name | treatment | toggles | schedule |
|---|---|---|---|
| `combined` (default) | RSL3 | immune + stromal + pH | constant |
| `bare` | RSL3 | none | constant |
| `multidose` | SDT | immune | **4-dose multi-dose (#239)** |
| `persister` | SDT | immune + **persister (#241)** | multi-dose |
| `clonal` | SDT | immune + **clonal 4-subclone (#242)** | multi-dose |
| `vasculature` | RSL3 | **explicit vessels (#191)** | constant |
| `spheroid` | SDT | immune + **radial biochem (#197)** | constant |
| `slab` | SDT | immune + **patient-scale slab (#240)** | constant |

The `slab` preset visualizes the **depth-graded supply** of a patient-scale slab: a surface slab (+z face = vessel) where the top layers are well-perfused and die while the deeper layers go drug/O2-deprived and survive. The depth axis is the layer (z) axis, which the renderer's mid-slice spans — so the death front in the existing dead/DAMP/LP panels *is* the visualization (no extra static overlay). See [Patient-scale slab](#patient-scale-slab-240) below.

The `multidose` preset shows **death waves synced to each dose**: four SDT ROS pulses at steps 10/55/100/145, each triggering a ferroptotic death wave + DAMP bloom + immune response. The renderer draws a red frame border + `💉 DOSE` marker on each dose step.

The `persister` preset adds the drug-tolerant persister model (#241) and a **4th render panel** colouring each cell by its `persister_fraction` (0..1): tolerance accumulates in survivors across the death waves and reverts between doses. It writes an extra `trajectory_persister.npy` (f32), and `summary.json`-equivalent runs report `persister_mean`. Only this preset emits the persister file, so the other presets render the original three panels unchanged.

## Time-varying dosing (`--dose-sweep` + `DoseSchedule`, #239)

By default the simulator models drug as present at **constant** strength for the whole run. `ferroptosis-core::dose_schedule::DoseSchedule` adds time-varying administration: `Constant` (the byte-identical default), `Bolus`, `MultiDose`, `Infusion`, and `FromPk` (driven by the two-compartment `tumor_pk` ODE).

- **SDT/PDT**: the schedule scales the per-step exogenous-ROS bolus. The intrinsic single-bolus decay envelope inside `sim_cell_step` is divided out on the dosed path (`biochem::exo_decay_factor`), so the schedule is the sole availability envelope — otherwise later doses would be wrongly crushed by the from-t=0 decay.
- **RSL3**: the schedule drives per-step covalent GPX4 inactivation (`gpx4 -= RSL3_INACTIVATION_RATE · availability · gpx4`, the validated `tumor_pk::sim_cell_with_pk` model) instead of the one-shot init knockdown. pH ion-trapping composes as a per-cell spatial availability multiplier.

`--dose-sweep` runs RSL3 across all five protocols at a fixed combined-TME context and writes `dose_comparison.json`:

```bash
cargo run --release -p sim-tme-3d -- --dose-sweep
# → output/tme-3d/dose_comparison.json  (one row per protocol, shared grid + RNG seed)
```

All protocols share the same tumor grid and runtime RNG seed, so kill-rate differences reflect the dosing protocol alone, not stochastic noise. On this machine the steady-state `constant` model kills ~10–60× more than any realistic time-varying protocol — a direct quantification of how much the "drug present at full strength forever" idealization overestimates efficacy. **Absolute magnitudes are uncalibrated v1** (`RSL3_INACTIVATION_RATE` was tuned for sustained `conc=1.0`); the informative signal is the cross-protocol ordering.

**Bit-identical guarantee**: when every condition uses `DoseSchedule::Constant` (the entire default 24-condition matrix), the run is byte-identical to pre-#239 — `summary.json` SHA unchanged. The dose machinery is gated behind `DoseSchedule::is_constant()`; `--dose-sweep` writes a separate file and never touches `summary.json`.

## Drug-tolerant persister cells (#241)

Cells that survive a ferroptosis inducer can enter an epigenetic **drug-tolerant persister** state (Hangauer 2017, Tsoi 2018): they resist the covalent GPX4 knockdown and enrich protective MUFA membrane lipids, then revert once the drug clears. The model lives in `ferroptosis-core::persister` (pure helpers) + `PersisterConfig`; this binary wires it into the per-step loop.

- **Off by default.** The matrix path passes no config, so the persister code path is inert and `summary.json` stays byte-identical (guarded by the #253 production SHA check + the `persister_off_is_inert_and_unreported` unit test). `persister_mean` is omitted from `summary.json` unless the model is on.
- **Two axes.** Under drug exposure each surviving cell's `persister_fraction` grows logistically toward a cap; it (a) attenuates the per-step RSL3 GPX4 inactivation (`gpx4_inactivation_multiplier`) and (b) adds MUFA protection (`mufa_boost_increment`). Between doses it reverts exponentially.
- **Observable decline.** `persister_reduces_multidose_kills` exercises a repeated-dose RSL3 run (uniform O₂, 20³ × 120): enabling the model cuts kills from 79 to 27 as survivors acquire tolerance — the Hangauer 2017 multi-cycle effect.
- **Visualization.** `--snapshot=persister` (SDT multi-dose + immune + persister) writes `trajectory_persister.npy` and the renderer adds a 4th panel colouring each cell by `persister_fraction` (0..1); a representative run reached `persister_mean ≈ 0.49` in survivors.

Parameters in `PersisterConfig::enabled()` are **plausible placeholders pending calibration** (the literature gives qualitative direction, not step-level rates).

## T-cell exhaustion (#243, Phase 1)

The spatial immune model (`immune_spatial`, #188) is a 0–48 h resident T-cell cascade with a single PD-1 brake. Phase 1 of #243 adds **T-cell exhaustion**: sustained killing in a region drives local T cells toward dysfunction, lowering their per-encounter kill probability (Wherry, Nat Immunol 2011; Snell et al., Cell 2018).

- **Model.** A per-cell `cumulative_kills` field counts immune kills accumulated in each cell's Moore-26 neighborhood; the kill probability is scaled by `ferroptosis_core::immune_spatial::exhaustion_factor` = `1 / (1 + exhaustion_rate · cumulative_kills)`.
- **Off by default.** `SpatialImmuneConfig::for_3d()` sets `exhaustion_rate = 0.0`, so `exhaustion_factor ≡ 1.0` and the `cumulative_kills` field is never allocated — `summary.json` is byte-identical (golden tests + the #253 production-SHA guard pass). The scatter that updates neighborhoods runs only when the rate is > 0.
- **Effect.** With exhaustion on, total immune kills decline as killing clusters (the "cold tumor" emergence). The `exhaustion_reduces_immune_kills` test shows a dense SDT + immune run dropping ≈20% (174 → 139) when exhaustion is enabled, all else equal. Note this also shifts a few deaths into the ferroptosis tally (a cell spared an apoptotic immune kill can die ferroptotically and release iron), so only the immune-kill count is asserted.
- **Scope.** Phase 1 only. Later phases (Treg/MDSC suppressor field, multi-checkpoint axis, DC subsets) stay separate per #243. `exhaustion_rate` is an uncalibrated placeholder.

## Clonal heterogeneity (#242)

Real tumors are genetic mosaics: 4–8+ subclones with measurably different ferroptosis vulnerabilities coexist in spatial patches, and the between-subclone variance often exceeds the between-treatment variance in real drug screens (Marusyk 2014; Conrad 2018; Heindl 2019).

- **Model.** `ferroptosis_core::clonal::assign_subclones_3d` partitions tumor cells into K Voronoi patches (seed points placed by an **independent** RNG, so `TumorGrid3D::generate`'s stream — and the cell grid — is untouched). Each subclone carries a `SubclonePerturbation` applied as RNG-neutral setup mutations, like the O2/pH gradients: `iron_mul` and `lipid_unsat_mul` (the MUFA-enrichment axis) scale **static** `Cell` fields so they persist across steps; `gpx4_mul` scales the initial `state.gpx4`, which strongly shapes the early autocatalytic window but relaxes toward the NRF2 setpoint over the run (a fully durable GPX4 axis would also scale the static `cell.nrf2`; deferred to calibration). `ClonalConfig::literature_4()` spans the mesenchymal⇄epithelial axis (subclone 1 iron-loaded/GPX4-low = vulnerable … subclone 4 GPX4-high/low-PUFA = resistant).
- **Off by default / byte-identical.** The matrix passes no clonal config; `ClonalConfig::single_identity()` (K=1, identity perturbation) is a no-op. The `clonal_k1_identity_is_byte_identical` test + the #253 production-SHA guard hold.
- **Reporting.** When enabled, `summary.json` gains a `subclone_kills` array (per-subclone `total_tumor` / `total_dead` / `kill_rate`); omitted otherwise. `clonal_subclones_differ_in_kill_rate` confirms the vulnerable subclone out-dies the resistant one under RSL3.
- **Visualization.** `--snapshot=clonal` writes a static `subclone.npy` (u8, no time axis) and the renderer adds a discrete-colored subclone-id panel alongside the death/DAMP/LP animation.
- Perturbation values are **placeholders pending calibration**. Clonal evolution (subclone selection over time) and inter-tumor heterogeneity are out of scope.

## Explicit vasculature (#191)

The 2D sims and the default 3D `oxygen::radial_o2_field` use "distance from the tumor edge" as a vasculature proxy: the spheroid surface is the only O2/drug source, so supply decays smoothly inward. Real 3D tumors carry **internal** vessels, so oxygenation is patchy — well-supplied near a vessel, hypoxic in the inter-vessel gaps (Option A of #191).

- **Model.** `ferroptosis_core::vasculature::place_vessels_3d` scatters random internal vessel seed points (count from the target inter-vessel spacing; **independent** RNG ⇒ the cell grid is untouched), and `vessel_supply_field` gives each cell `exp(-dist_to_nearest_vessel / λ)`. This is a **unified per-cell supply factor**: it **replaces** the edge-distance O2 factor (× `cell.basal_ros`) AND scales drug delivery on every path — the SDT/PDT exo dose at init (constant + dosed), the RSL3 constant one-shot knockdown (a post-init correction mirroring the pH ion-trap one), and the dosed RSL3 `rsl3_drug_avail`. So the `--snapshot=vasculature` RSL3 run exercises both axes. `VasculatureConfig::well_vascularized()` / `poorly_vascularized()` set the vessel density.
- **Off by default / byte-identical.** The matrix passes no vasculature config, so the edge-distance `radial_o2_field` path runs unchanged (golden + #253 production-SHA guard hold).
- **Comparison (AC).** The library test `vessel_field_oxygenates_the_core_unlike_the_edge_proxy` shows internal vessels reach the deep core (which the edge proxy leaves uniformly hypoxic); the sim test `vasculature_oxygenates_core_and_changes_rsl3_kills` shows the vessel field materially changes RSL3 kills (direction is config-dependent — a sparse internal vessel set covers the surface shell less uniformly than the edge proxy). `summary.json` reports `vascular_hypoxic_fraction` when enabled.
- **Visualization.** `--snapshot=vasculature` (RSL3 + vessels) writes a static `vessel_supply.npy` and the renderer adds a cividis O2-supply panel — bright vessel neighborhoods, dark hypoxic gaps.
- Vessel placement is random (Option A); a fractal-branching tree or imported micro-CT networks (Options B/C) are out of scope. Inter-vessel spacings are literature-ranged placeholders.

## 3D spheroid radial biochemistry (#197)

Beyond the *environmental* fields (O2 #187, pH #190, vessels #191), the ferroptosis cascade itself is position-dependent in a spheroid: the phenotype distribution is radially structured, and MUFA/GSH/iron vary with distance from the nutrient supply.

- **Model.** Opt-in via `Overrides.spheroid`. Runs under `Params::spheroid()` (MUFA partially active). `ferroptosis_core::spheroid::apply_radial_cells_3d` **re-generates** each tumor cell with its radial-position phenotype — **glycolytic rim / OXPHOS mid / persister-like core** — using an **independent** RNG (so `generate`'s stream is untouched), then applies the GSH (core-low, cysteine-limited) and iron (core-high, HIF-driven import) gradients to the cell. Position-dependent MUFA (peripheral high, core low) is set on the initialized `CellState`; under `Params::spheroid()` (partial SCD1) it is non-zero (unlike the #265 2D-params case where it would reset to 0 on step 1), but it is an **initial condition** that relaxes toward the uniform M_ss (≈0.20) over the run and is capped at `scd_mufa_max` (0.25) — so it shapes the early killing window rather than persisting. Of the three radial axes, only **iron** (static `cell.iron`) is fully durable; `gsh` likewise sets the initial value. Fully durable position-dependent MUFA needs per-cell `scd_mufa_max` (deferred — see #197 refinements).
- **Off by default / byte-identical.** The matrix passes no spheroid config, so it runs `Params::default()` on the random grid — unchanged (golden + #253 production-SHA guard hold). The radial re-gen and `Params::spheroid()` switch are both gated on the opt-in.
- **Comparison (AC).** `radial_spheroid_changes_kills_vs_random` shows the radial structure materially changes the kill outcome vs the random-phenotype grid (answering "does radial structure change the kill rate, or just redistribute deaths?" — it changes it). `spheroid_params_differ_from_default` pins `Params::spheroid()` between 2D (M_ss 0) and in-vivo (0.40).
- **Visualization.** `--snapshot=spheroid` writes a static `phenotype.npy` (u8) and the renderer adds a discrete phenotype panel (stroma + glycolytic/OXPHOS/persister zones).
- The GSH and iron gradients are the issue's "optional v1" items, included here. Cell-cell-contact (E-cadherin/Merlin) and glucose/glutamine/lactate gradients remain out of scope; gradient strengths are placeholders pending calibration.

## Patient-scale slab (#240)

The spheroid is ~540 µm radius — **in-vitro** scale, where drug/O2 penetration (Krogh ~150 µm) reaches most cells. Real patient tumors are **5–50 mm**, where penetration drops catastrophically with depth and the bulk of the tumor is essentially drug-deprived. The `slab` mode models a chunk of such a tumor at a configurable depth, so the efficacy numbers reflect the patient-scale penetration limit rather than the optimistic spheroid scale.

- **Geometry.** Opt-in via `Overrides.slab`. `TumorGrid3D::generate_slab` builds an **all-tumor block** (no sphere carve, no stroma — a slab of bulk tumor), and `ferroptosis_core::slab::slab_supply_field` gives each cell a **planar** depth-graded supply: the **+z face** (top layer) is vessel-proximal at `depth_offset_mm`, and supply decays `exp(-depth/λ)` going **−z** (deeper into the tumor). This is the 1-D analog of `oxygen::radial_o2_field`, and like the vasculature field it is a **unified supply factor** — it replaces the O2 factor (× `cell.basal_ros`) AND scales drug delivery on every path (exo dose at init, RSL3 constant knockdown correction, dosed `rsl3_drug_avail`). `SlabConfig::patient_5mm()` places the +z face at 4 mm (the slab spans ~4.0–5.2 mm, essentially fully deprived); `SlabConfig::surface()` places it at 0 mm (the shallow, vessel-fed control).
- **Boundary conditions.** +z = vessel (the supply maximum, immune-accessible via the existing DAMP cascade). −z = continuing tumor, modeled as a **no-flux / reflective** face — already satisfied because the iron/DAMP diffusion uses bounded (no-wrap) neighbors, so there is no flux across any grid face. The depth gradient (high at +z, ~0 at −z) supplies the directionality; the reflective −z boundary holds because the planar supply replaces the field rather than diffusing it.
- **Off by default / byte-identical.** The matrix passes no slab config, so it runs the spheroid grid + edge-distance/vessel O2 unchanged (golden + #253 production-SHA guard hold). The slab grid and planar supply are both gated on the opt-in.
- **Reporting.** A slab run adds a `scale_interpretation` string to its `ConditionResult` (e.g. *"slab spanning depth 4.0–5.2 mm of a 5 mm virtual tumor (1.2 mm thick)"*), omitted on non-slab runs — so the output says plainly what physical scale it represents.
- **Comparison (AC).** `patient_scale_slab_kills_far_less_than_spheroid` shows a 4 mm-deep slab kills **<20%** of what the in-vitro spheroid does under the same SDT treatment — the penetration collapse the spheroid scale misses. `shallow_slab_outkills_deep_slab` confirms the depth-dependence within slab mode (surface slab out-kills the deep one).
- **Visualization.** `--snapshot=slab` (SDT + immune on a **surface** slab) renders the depth-graded death front directly in the dead/DAMP/LP panels: the renderer's mid-slice fixes the row axis and spans `(col, layer)`, so the layer (depth) axis is shown — top (+z, well-perfused) dies, bottom (−z, deprived) survives. No extra static overlay file. A deep `patient_5mm()` slab would kill ~nothing (a less illustrative GIF), so the surface slab is the preset; the depth comparison lives in the tests.

## Dimensionality-dependent assumptions (#197 AC)

Which modeling choices are 2D-vs-3D-specific, surfaced for honesty:
- **Radial phenotype structure** only exists in 3D spheroid mode; the 2D/default model assigns phenotypes by a coarse core/periphery split with random rolls.
- **MUFA**: 2D culture (`Params::default`) has none; the spheroid context has partial, position-graded MUFA; full in-vivo has uniform M_ss = 0.40.
- **O2 source**: 2D/default uses the edge-distance proxy (surface = only source); 3D can use explicit internal vessels (#191). The "hypoxia collapses RSL3" magnitude differs between the two (see #191).
- These are model *contexts*, not calibrated transitions; the spheroid parameters are placeholders.

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

### Production byte-identity regression (#253)

Two layers guard the load-bearing invariant that the #239 multi-dose work relies on (default matrix = all `DoseSchedule::Constant` = byte-identical `summary.json`):

1. **Per-PR (fast):** `constant_path_golden_kill_counts` pins integer kill counts at a small 20³ × 80 config, so a structural regression fails in ordinary `cargo test`.
2. **Off-PR (full scale):** the `.github/workflows/sim-tme-3d-regression.yml` job runs the real 60³ × 180 matrix and asserts `output/tme-3d/summary.json`'s SHA-256 against the checked-in `expected_summary.sha256`. It runs weekly and on `workflow_dispatch` (not per-PR: the release run + build is too heavy to gate every PR, and a scale-dependent change could pass the small golden yet perturb the full output).

The expected hash is **toolchain-specific** (pinned 1.96.0; confirmed identical on 1.92.0). When an *intentional* output change lands, regenerate it per the instructions inside `expected_summary.sha256`.

## Performance & scalability (`--bench`, #192)

`--bench` runs one representative condition (combined-TME RSL3) at a configurable grid size and prints wall-clock + throughput. Grid/steps come from env so a sweep is scriptable:

```bash
cd simulations
BENCH_GRID_DIM=200 BENCH_N_STEPS=180 cargo run --release -p sim-tme-3d -- --bench
# peak RSS: wrap with `/usr/bin/time -l` (macOS) or `/usr/bin/time -v` (Linux)
```

**All figures below are direct `--bench` wall-clock measurements at 180 steps** (no projections), rust 1.96.0, 10-core machine, `size_of::<GridCell>() = 144 B`. Serial figures are from the pre-parallelization commit (commit 1 of this PR); parallel figures from the within-condition rayon path (commit 2). Throughput columns are each path's own `cell_steps_per_s`. Peak RSS is process-wide, measured via `/usr/bin/time -l`.

| Grid | Cells | Peak RSS | Serial 180-step | Parallel 180-step | Speedup | Parallel throughput |
|---|---|---|---|---|---|---|
| 50³  | 125 k   | 23 MB   | 3.7 s   | 1.0 s  | 3.8× | 2.3e7 cell·step/s |
| 100³ | 1.0 M   | 164 MB  | 26.4 s  | 6.1 s  | 4.3× | 3.0e7 |
| 150³ | 3.375 M | 546 MB  | 87.6 s  | 18.3 s | 4.8× | 3.3e7 |
| 200³ | 8.0 M   | **1.29 GB** | 201.3 s (3.4 min) | **40.8 s** | 4.9× | 3.5e7 |

Both **performance targets are met even serially** (100³ < 2 min, 200³ < 15 min); the within-condition rayon parallelism (#192) adds a **3.8×–4.9× speedup** on 10 cores (the ratio grows with grid size as the fixed per-step rayon-join overhead amortizes). Serial throughput is ~7e6 cell·step/s across sizes.

**Memory verdict (feeds #240 patient-scale):** dense **200³ fits the 2 GB budget at 1.29 GB** with ~35% headroom — no sparse grid needed at this scale. Sparse/adaptive only become compelling at 300³+ (≈ 6 GB) or when running many large conditions concurrently; deferred to a follow-up issue.

**Recommended grid size:** 60³ for the 24-condition matrix (throughput); up to 150³ comfortably for single high-resolution runs; 200³ feasible at ~1.29 GB / ~43 s. Do **not** run 200³ across many concurrent conditions (24 × 1.29 GB would OOM) — throttle with `RAYON_NUM_THREADS` if needed.

**Parallelism note:** the default matrix parallelizes across the 24 conditions (`par_iter`); the biochem + immune-kill loops parallelize *within* a condition (rayon, byte-identical via position-independent per-cell RNG). Iron + DAMP diffusion stay serial (cross-cell dependencies). A single large `--bench` run has no condition-level parallelism, so within-condition rayon is what makes it fast.

The within-condition rayon is **nested** inside the condition-level `par_iter` on the default 24-condition matrix. Measured before/after, the **matrix wall-clock is unchanged** (~15 s on 10 cores, serial-within-condition vs parallel-within-condition, within run-to-run noise): the 24 conditions already saturate the pool, so the inner `par_iter_mut` adds no measurable overhead and finds no idle workers to exploit until the tail. The speedup applies to **single large-grid runs** (the #240 patient-scale direction), not the everyday 60³ matrix — which is the intended target.

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

- ~~**Lift `PhConfig` / `StromalConfig` / `ImmuneConfig` to `ferroptosis-core::params`**~~ — **done** in #220/#224 (lifted as `PhConfig` / `StromalConfig` / `SpatialImmuneConfig`).
- **O₂ cycling** (square-wave λ alternation) — sim-tme has it, sim-tme-3d skipped for v1 scope.
- **Anti-PD-1 sweep** — included in sim-tme; skipped here for v1.
- **3D volumetric visualization** — delivered in #193/#238 (axial-slice GIF/MP4 renderer + `.npy` volumetric trajectory arrays + 2D-vs-3D comparison table). #193 closed as substantially delivered; a ParaView-grade VTK/HDF5 export remains optional polish if a manuscript figure ever needs it.
- ~~**Larger grids**~~ — **demonstrated feasible** in #192: up to 200³ at ~1.29 GB / ~43 s (see Performance & scalability above). A standalone `sim-spatial-3d` binary (#194) was closed as superseded: `--snapshot=bare` already provides the unprotected depth-physics baseline.
- **Empirical pimonidazole validation** — see #196.
