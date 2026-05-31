//! sim-tme-3d: 3D spheroid tumor microenvironment simulation.
//!
//! Capstone 3D binary for the spheroid-validation series (#185–#197).
//! Integrates all five library primitives landed in v0.7.0–v0.11.0:
//! - 3D energy physics (#186) via `physics::local_ros_multiplier_3d`
//! - 3D radial O₂ gradient (#187) via `oxygen::radial_o2_field`
//! - 3D radial pH gradient (#190) via `ph::radial_ph_field` + helpers
//! - 3D CAF-shielded boundary detection (#189) via `stromal::stromal_adjacency_mask_3d`
//! - 3D spatial DAMP diffusion + activation (#188) via `immune_spatial::*`
//!
//! Produces a per-condition matrix of kill rates that the Python
//! comparison script (`scripts/generate_3d_comparison_table.py`) pairs
//! with sim-tme's existing 2D output to answer the four key questions
//! in issue #195:
//! 1. Does the hypoxia RSL3 collapse hold in 3D?
//! 2. Does the immune 104:1 ratio hold in 3D?
//! 3. Does stromal shielding have MORE impact in 3D?
//! 4. Does pH ion trapping produce similar RSL3 reduction in 3D?
//!
//! ## ⚠️ Scale-mismatch with sim-tme
//!
//! sim-tme uses a 500×500 grid (tumor radius ≈ 4500 µm — large in-vivo
//! tumor). sim-tme-3d uses a **60³ grid** (tumor radius ≈ 540 µm —
//! upper end of in-vitro spheroids) as the default for the 24-condition
//! matrix. Larger single-condition grids are now feasible after the #192
//! perf work: measured 200³ × 180 ≈ 41 s at ~1.29 GB on 10 cores (see
//! `--bench` + the "Performance & scalability" section of the README).
//! Only 500³+ (≈ 18 GB dense) remains out of reach. The 60³ default is a
//! deliberate matrix-throughput choice, not a ceiling.
//!
//! **The Python comparison script reports RATIOS** (e.g., RSL3 hypoxic
//! kill / RSL3 normoxic kill) which are dimensionally meaningful at
//! different scales. Absolute kill counts are NOT directly comparable
//! between the two binaries.
//!
//! ## Why λ sweep skips 150 µm
//!
//! 3D hypoxic-zone threshold = 3λ. At λ=150 µm, threshold = 450 µm,
//! but the 60³ grid's tumor radius is only 540 µm. The hypoxic zone
//! would be ≤1 cell — statistically meaningless. We sweep [80, 100,
//! 120] which give meaningful hypoxic shells at this scale.
//!
//! ## Stability requirement (immune diffusion)
//!
//! Sim-tme's 2D `damp_diffusion_fraction = 0.08` is **unsafe in 3D**
//! (0.08 × 26 = 2.08 > 1, would mass-destroy). We use 0.025 (matches
//! 2D's per-step total diffusion of ~64% — see `immune_spatial` rustdoc).
//! `immune_spatial::diffuse_damp_3d_step` enforces the stability invariant
//! with `assert!` (release-mode panic).

use std::fs;
use std::path::Path;
use std::time::Instant;

mod npy;
mod snapshot;

use ferroptosis_core::biochem::{exo_decay_factor, sim_cell_step, CellState};
use ferroptosis_core::cell::Treatment;
use ferroptosis_core::dose_schedule::DoseSchedule;
use ferroptosis_core::grid::{TumorGrid3D, TUMOR_RADIUS_FRACTION};
use ferroptosis_core::immune_spatial::{
    dc_activation, diffuse_damp_3d_step, immune_kill_probability, DAMP_KILL_THRESHOLD,
};
use ferroptosis_core::oxygen::radial_o2_field;
use ferroptosis_core::params::{
    Params, PhConfig, SpatialImmuneConfig, SpatialParams, StromalConfig,
};
use ferroptosis_core::ph::{ion_trap_factor_from_ph, iron_multiplier_from_ph, radial_ph_field};
use ferroptosis_core::physics::local_ros_multiplier_3d;
use ferroptosis_core::stromal::stromal_adjacency_mask_3d;
use ferroptosis_core::tumor_pk::RSL3_INACTIVATION_RATE;
use rand::distributions::Distribution;
use rand::prelude::*;
use rand::rngs::StdRng;
use rand_distr::Normal;
use rayon::prelude::*;
use serde::Serialize;

// ============================================================
// Constants
// ============================================================

const GRID_DIM: usize = 60;
const CELL_SIZE_UM: f64 = 20.0;
const N_STEPS: u32 = 180;
const SEED: u64 = 42;

/// O₂ penetration sweep — 3λ must comfortably fit inside the
/// tumor_radius (60·0.45·20 = 540 µm). Skips λ=150 (3·150=450 →
/// only ~1 cell hypoxic zone, statistically meaningless).
const O2_LAMBDAS: &[f64] = &[80.0, 100.0, 120.0];

/// Zone-analysis reference λ — matches sim-tme.
const ZONE_REF_LAMBDA: f64 = 120.0;

/// DAMP diffusion fraction for 3D, retained as a const for legacy code
/// paths that emit the value into output metadata (see `RunMetadata` /
/// `summary.json`). Identical to `SpatialImmuneConfig::for_3d()
/// .damp_diffusion_fraction`. **Not** the source of truth for the
/// runtime value — that's the library const propagated through
/// `SpatialImmuneConfig`.
const DAMP_DIFFUSION_FRACTION_3D: f64 = 0.025;

/// Step at which immune activity begins (matches sim-tme).
const IMMUNE_START_STEP: u32 = 60;

/// Iron diffusion neighbor fraction for 3D — scaled from sim-tme's 2D
/// default `0.1` (which is per-Moore-8) to the per-Moore-26 equivalent:
/// `0.1 × 8/26 ≈ 0.0308`. Without scaling, 2D's `0.1` would over-spread
/// in 3D (10% × 26 = 260% per source-cell loss, vs 2D's 80%).
const IRON_DIFFUSE_FRACTION_3D: f64 = 0.1 * 8.0 / 26.0;

// `DAMP_KILL_THRESHOLD`, `SpatialImmuneConfig`, `StromalConfig`,
// `PhConfig` are imported from `ferroptosis-core` (#220). Previously
// these were binary-local copies duplicated with sim-tme.

// ============================================================
// Output records
// ============================================================

/// Per-condition kill-rate result. Mirrors `sim-tme::ConditionResult`
/// structure so the comparison script can pair fields across the two
/// JSONs.
#[derive(Clone, Debug, Serialize)]
struct ConditionResult {
    treatment: String,
    o2_condition: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    o2_lambda_um: Option<f64>,
    immune_mode: String,
    total_tumor: usize,
    total_dead: usize,
    #[serde(skip_serializing_if = "Option::is_none")]
    ferroptosis_kills: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    immune_kills: Option<usize>,
    overall_kill_rate: f64,
    normoxic_kill_rate: f64,
    transition_kill_rate: f64,
    hypoxic_kill_rate: f64,
    /// Peak per-cell DAMP value at end of simulation (max over all cells).
    /// Always populated (the binary computes the mask + DAMP regardless of
    /// immune-on toggle so cross-condition comparison is apples-to-apples).
    /// PR body's "DAMP fields" claim is vindicated by this + total_damp.
    peak_damp: f64,
    /// Sum of DAMP across all cells at end of simulation. Mass-balance
    /// indicator for the diffusion-clearance dynamics.
    total_damp: f64,
    #[serde(skip_serializing_if = "Option::is_none")]
    stromal_mode: Option<String>,
    /// Kill rate among cells flagged by `stromal_adjacency_mask_3d` (boundary
    /// cells with ≥1 stromal Moore-26 neighbor). **Always populated** —
    /// the binary computes the mask regardless of `stromal_on` toggle so
    /// Q3 in the comparison script can pair stromal-on adjacent rates
    /// with no-stromal baseline adjacent rates (reviewer-flagged fix).
    /// `None` only when grid produces zero boundary cells.
    #[serde(skip_serializing_if = "Option::is_none")]
    stromal_adjacent_kill_rate: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    stromal_adjacent_count: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    ph_mode: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    ph_edge: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    ph_core: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    ph_lambda_um: Option<f64>,
}

/// Schema version for `summary.json`. Bump when the output shape
/// changes; `scripts/generate_3d_comparison_table.py` cross-checks
/// this against `sim-tme/tme_summary.json`'s `schema_version` to catch
/// schema drift across the two binaries (#224 item 2).
const TME_3D_SCHEMA_VERSION: u32 = 1;

#[derive(Clone, Debug, Serialize)]
struct SimulationSummary {
    /// Bumped when the shape of this JSON changes. Currently `1`.
    schema_version: u32,
    grid_dim: usize,
    cell_size_um: f64,
    tumor_radius_um: f64,
    n_steps: u32,
    o2_lambdas: Vec<f64>,
    damp_diffusion_fraction_3d: f64,
    conditions: Vec<ConditionResult>,
    /// Scale-mismatch caveat embedded in the output so consumers see it.
    note: String,
}

// ============================================================
// Condition spec — every entry in the matrix is one of these.
// ============================================================

#[derive(Clone, Debug)]
struct Condition {
    name: String,
    treatment: Treatment,
    treatment_name: String,
    /// `None` = uniform O₂ (baseline); `Some(λ)` = radial O₂ gradient.
    o2_lambda: Option<f64>,
    immune_on: bool,
    stromal_on: bool,
    ph_on: bool,
    /// Drug-administration schedule over time (#239). `Constant` (the
    /// default) reproduces the historical steady-state behavior exactly;
    /// non-constant schedules drive per-step drug modulation.
    dose_schedule: DoseSchedule,
}

// ============================================================
// Zone-kill-rate analog of sim-tme's `zone_kill_rates`. Uses
// TumorGrid3D::radial_depth_um directly — same depth thresholds,
// geometry-agnostic.
// ============================================================

fn zone_kill_rates_3d(grid: &TumorGrid3D, shell_depth_um: f64) -> (f64, f64, f64) {
    let deep_threshold_um = shell_depth_um * 3.0;
    let (mut nd, mut nt, mut td, mut tt, mut hd, mut ht) = (0, 0, 0, 0, 0, 0);
    for r in 0..grid.rows {
        for c in 0..grid.cols {
            for l in 0..grid.layers {
                let gc = grid.get(r, c, l);
                if !gc.is_tumor {
                    continue;
                }
                let depth = grid.radial_depth_um(r, c, l);
                let (dead, total) = if depth < shell_depth_um {
                    (&mut nd, &mut nt)
                } else if depth < deep_threshold_um {
                    (&mut td, &mut tt)
                } else {
                    (&mut hd, &mut ht)
                };
                *total += 1;
                if gc.state.dead {
                    *dead += 1;
                }
            }
        }
    }
    let rate = |d: usize, t: usize| if t > 0 { d as f64 / t as f64 } else { 0.0 };
    (rate(nd, nt), rate(td, tt), rate(hd, ht))
}

// ============================================================
// Per-condition runner — the 3D analog of sim-tme's
// `run_spatial_with_immune`. Uses library primitives throughout.
// ============================================================

fn norm(rng: &mut StdRng, mean: f64, std: f64) -> f64 {
    if std <= 0.0 {
        return mean;
    }
    let dist = Normal::new(mean, std).expect("valid normal");
    dist.sample(rng)
}

/// Grid/step config for a run. Production uses `RunConfig::production()`;
/// tests use `RunConfig::for_test()` to avoid running the full 60³ × 180-step
/// simulation in CI debug builds (reviewer-flagged perf concern).
#[derive(Clone, Copy, Debug)]
struct RunConfig {
    grid_dim: usize,
    n_steps: u32,
}

impl RunConfig {
    fn production() -> Self {
        RunConfig {
            grid_dim: GRID_DIM,
            n_steps: N_STEPS,
        }
    }

    /// Tiny config for unit tests — exercises every code path but ~1000× cheaper.
    /// 10³ = 1000 cells, 20 steps. Runs in <100ms in debug mode.
    #[cfg(test)]
    fn for_test() -> Self {
        RunConfig {
            grid_dim: 10,
            n_steps: 20,
        }
    }
}

fn run_one_condition(condition: &Condition) -> ConditionResult {
    run_one_condition_with_config(condition, RunConfig::production(), None)
}

fn run_one_condition_with_config(
    condition: &Condition,
    run_cfg: RunConfig,
    mut snapshot: Option<&mut snapshot::SnapshotBuffers>,
) -> ConditionResult {
    let params = Params::default();
    let spatial_params = SpatialParams {
        cell_size_um: CELL_SIZE_UM,
        // SpatialParams::default()'s `neighbor_iron_fraction = 0.1` is
        // 2D-tuned (per-Moore-8). Override here with the 3D-scaled value
        // so the field is correct in the struct, not just at the inline
        // call site (reviewer-flagged future-refactor footgun: if a
        // future `local_ros_multiplier_3d` ever reads this field, the
        // 2D value would silently propagate).
        neighbor_iron_fraction: IRON_DIFFUSE_FRACTION_3D,
        ..Default::default()
    };

    // Per-condition deterministic seed for the *runtime RNG* (per-cell ROS
    // noise, immune-kill rolls, etc.). NOT used for grid generation — sim-tme
    // uses the fixed `SEED` so every condition operates on the SAME tumor
    // geometry (same cell phenotypes, same Persister cluster placements).
    // Using a per-condition seed for the grid would mean each treatment is
    // evaluated on a DIFFERENT tumor — combining treatment effects with
    // grid-geometry noise. At ~82.5k tumor cells (60³ sphere × 0.45 radius
    // fraction), with ~268 cells/cluster × 10-20 random clusters, that
    // noise is material relative to subtle effects like the immune ratio.
    let cond_seed = SEED.wrapping_add(hash_condition_name(&condition.name));

    // Grid uses fixed SEED — matches sim-tme/src/main.rs (line 921 and
    // every other TumorGrid::generate call uses bare `SEED`). Same tumor
    // geometry across all conditions ⇒ valid treatment-effect comparison.
    let mut grid = TumorGrid3D::generate(
        run_cfg.grid_dim,
        run_cfg.grid_dim,
        run_cfg.grid_dim,
        CELL_SIZE_UM,
        SEED,
    );
    let n_cells = grid.cells.len();

    // --- Apply O₂ gradient if requested (mutates cell.basal_ros) ---
    if let Some(lambda) = condition.o2_lambda {
        let o2_factors = radial_o2_field(&grid, lambda);
        for (idx, &factor) in o2_factors.iter().enumerate() {
            if grid.cells[idx].is_tumor {
                grid.cells[idx].cell.basal_ros *= factor;
            }
        }
    }

    // --- Apply pH gradient if requested (mutates cell.iron via library helper) ---
    let ph_field = if condition.ph_on {
        let cfg = PhConfig::default();
        let field = radial_ph_field(&grid, cfg.ph_edge, cfg.ph_core, cfg.lambda_ph_um);
        for (idx, &local_ph) in field.iter().enumerate() {
            if grid.cells[idx].is_tumor {
                let iron_mult =
                    iron_multiplier_from_ph(local_ph, cfg.ph_edge, cfg.iron_ph_sensitivity);
                grid.cells[idx].cell.iron *= iron_mult;
            }
        }
        Some((field, cfg))
    } else {
        None
    };

    // --- Compute stromal adjacency mask ALWAYS ---
    // The mask is grid-geometry-dependent (not stromal_on-toggle-dependent),
    // so it's identical across all conditions on the same grid. Computing it
    // for every condition (cheap: O(N·26) one-time) lets the comparison
    // script's Q3 pair stromal-on adjacent rates against the matching
    // no-stromal baseline adjacent rates — measuring CAF shielding at the
    // boundary directly, not as a diluted whole-tumor delta.
    let adjacency_mask = stromal_adjacency_mask_3d(&grid);

    // Apply CAF GSH/MUFA boosts only when stromal_on is true.
    let stromal_cfg = if condition.stromal_on {
        Some(StromalConfig::default())
    } else {
        None
    };

    // --- Treatment-specific ROS multiplier (3D version) ---
    // PDT arm is intentionally kept for sim-tme parity (sim-tme/main.rs:573
    // has the same pattern); `generate_conditions()` does not include PDT
    // in the v1 matrix (deferred to follow-up — manuscript focuses on the
    // RSL3-vs-SDT comparison, with PDT in a separate sim-spatial track).
    // The arm is dead code in this binary today but lets a future caller
    // pass `Treatment::PDT` without a match-exhaustiveness error.
    let base_ros = match condition.treatment {
        Treatment::SDT => params.sdt_ros,
        Treatment::PDT => params.pdt_ros,
        Treatment::RSL3 | Treatment::Control => 0.0,
    };

    // Time-varying dose schedule (#239). `dosed == false` for the default
    // `Constant` schedule, in which case EVERY dose-related branch below is
    // skipped and the run is byte-identical to the pre-#239 behavior.
    let schedule = &condition.dose_schedule;
    let dosed = !schedule.is_constant();
    // `base_exo` is read ONLY by the SDT/PDT per-step rescale arm. Allocate
    // (and init-write) it only when that arm can run — i.e. a non-constant
    // schedule on an exogenous-ROS treatment. Empty (never indexed) on the
    // Constant path AND on the RSL3 dosed path, which keeps the
    // "empty ⇒ never read" invariant tight (review #7).
    let dose_modulates_exo =
        dosed && matches!(condition.treatment, Treatment::SDT | Treatment::PDT);
    let mut base_exo: Vec<f64> = if dose_modulates_exo {
        vec![0.0; n_cells]
    } else {
        Vec::new()
    };

    // Initialize cell states with per-cell ROS peak
    for r in 0..grid.rows {
        for c in 0..grid.cols {
            for l in 0..grid.layers {
                let idx = grid.flat_index(r, c, l);
                let exo_ros_peak =
                    if matches!(condition.treatment, Treatment::Control | Treatment::RSL3) {
                        0.0
                    } else {
                        let depth_um = grid.radial_depth_um(r, c, l);
                        let ros_multiplier =
                            local_ros_multiplier_3d(depth_um, condition.treatment, &spatial_params);
                        let mut rng = StdRng::seed_from_u64(cond_seed.wrapping_add(idx as u64));
                        let peak = base_ros * ros_multiplier;
                        norm(&mut rng, peak, peak * 0.2).max(0.0)
                    };
                if dose_modulates_exo {
                    base_exo[idx] = exo_ros_peak;
                }
                let gc = &mut grid.cells[idx];
                // For RSL3 under a non-constant schedule, skip the one-shot
                // init GPX4 knockdown — the schedule drives per-step covalent
                // inactivation instead (the `tumor_pk::sim_cell_with_pk`
                // model). Every other case keeps the historical one-shot
                // init via `from_cell_with_ros` (= `..._opts(.., true)`),
                // preserving byte-identical output on the default path.
                gc.state = if dosed && condition.treatment == Treatment::RSL3 {
                    CellState::from_cell_with_ros_opts(
                        &gc.cell,
                        condition.treatment,
                        &params,
                        exo_ros_peak,
                        false,
                    )
                } else {
                    CellState::from_cell_with_ros(
                        &gc.cell,
                        condition.treatment,
                        &params,
                        exo_ros_peak,
                    )
                };
                gc.extra_iron = 0.0;
                gc.newly_dead = false;
                // Init to NaN so any code path that reads `lp_at_death`
                // before writing it (grace-end write or end-of-sim
                // catch-all) produces NaN downstream — calibration trips
                // instead of silently using a stale value (#225 review).
                gc.lp_at_death = f64::NAN;
            }
        }
    }

    // pH-dependent RSL3 ion trapping correction (consumer-side, same pattern
    // as sim-tme). This is the CONSTANT-schedule path: it corrects the
    // one-shot init knockdown for spatial drug availability. For a non-constant
    // schedule the init knockdown was skipped, so this correction would have
    // nothing to correct — instead the per-cell availability is folded into
    // the per-step inactivation via `rsl3_drug_avail` below.
    if !dosed {
        if let (Some((ph_map, cfg)), Treatment::RSL3) = (&ph_field, condition.treatment) {
            for (idx, &local_ph) in ph_map.iter().enumerate() {
                if !grid.cells[idx].is_tumor {
                    continue;
                }
                let drug_factor =
                    ion_trap_factor_from_ph(local_ph, cfg.ph_edge, cfg.ion_trap_sensitivity);
                // Correct GPX4: from (1-inhib) to (1-inhib*drug_factor) — matches sim-tme:614
                let correction =
                    (1.0 - params.rsl3_gpx4_inhib * drug_factor) / (1.0 - params.rsl3_gpx4_inhib);
                grid.cells[idx].state.gpx4 *= correction;
            }
        }
    }

    // Per-cell spatial drug-availability multiplier for the DOSED RSL3 path
    // (#239): pH ion-trapping reduces effective drug in acidic regions. The
    // value is composed multiplicatively with the schedule's `factor_at(step)`
    // in the per-step inactivation below. `1.0` everywhere when pH is off;
    // empty (never indexed) on the Constant path or for non-RSL3 treatments.
    let rsl3_drug_avail: Vec<f64> = if dosed && condition.treatment == Treatment::RSL3 {
        let mut avail = vec![1.0_f64; n_cells];
        if let Some((ph_map, cfg)) = &ph_field {
            for (idx, &local_ph) in ph_map.iter().enumerate() {
                if grid.cells[idx].is_tumor {
                    avail[idx] =
                        ion_trap_factor_from_ph(local_ph, cfg.ph_edge, cfg.ion_trap_sensitivity);
                }
            }
        }
        avail
    } else {
        Vec::new()
    };

    // --- Main 180-step loop ---
    let immune_cfg = SpatialImmuneConfig::for_3d();
    let mut damp_field = vec![0.0_f64; n_cells];
    let mut damp_scratch = vec![0.0_f64; n_cells];
    let mut ferroptosis_kills = 0usize;
    let mut immune_kills = 0usize;

    // Hoisted out of the per-cell hot loop (these invariants don't vary by
    // cell or step): on the dosed path the relevant availability vec must be
    // populated. Compiled out in release.
    if dosed {
        match condition.treatment {
            Treatment::RSL3 => debug_assert!(
                !rsl3_drug_avail.is_empty(),
                "rsl3_drug_avail must be populated on the dosed RSL3 path"
            ),
            Treatment::SDT | Treatment::PDT => debug_assert!(
                !base_exo.is_empty(),
                "base_exo must be populated on the dosed SDT/PDT path"
            ),
            Treatment::Control => {}
        }
    }

    for step in 0..run_cfg.n_steps {
        // Drug availability for this step (#239). `1.0` on the Constant
        // default path (and the `dosed` guard below skips all modulation
        // there anyway, so it's never actually read for Constant).
        let dose_factor = if dosed { schedule.factor_at(step) } else { 1.0 };

        // Ferroptosis biochem + stromal protection.
        //
        // Parallelized over cells with rayon (#192). This is safe and
        // BYTE-IDENTICAL to the old serial r/c/l triple loop because:
        //   - `enumerate()` over the flat `cells` Vec yields the same `idx`
        //     as `flat_index(r,c,l)` (row-major), so the per-cell RNG seed
        //     `(cond_seed, idx, step)` is unchanged and position-independent;
        //   - each cell reads/writes ONLY its own `GridCell` and its own
        //     `damp_field[idx]` slot (paired via the zipped `par_iter_mut`);
        //     `extra_iron` is read-then-zeroed per own cell (cross-cell iron
        //     spread happens later in the serial `diffuse_iron` phase);
        //   - `ferroptosis_kills` is an integer sum, associative+commutative,
        //     so the reduction is order-independent.
        // Reads of `base_exo` / `rsl3_drug_avail` / `adjacency_mask` are
        // shared-immutable by own index. Iron + DAMP diffusion stay serial
        // (cross-cell deps) and run after the rayon join (a per-step barrier).
        //
        // The `zip` below relies on `cells` and `damp_field` having equal
        // length (both allocated `n_cells`); a shorter `damp_field` would
        // silently TRUNCATE the loop and drop cells. They're allocated full
        // size unconditionally above, but assert it so a future conditional
        // allocation can't introduce a silent correctness bug (review #255).
        debug_assert_eq!(grid.cells.len(), damp_field.len());
        let died_this_step: usize = grid
            .cells
            .par_iter_mut()
            .zip(damp_field.par_iter_mut())
            .enumerate()
            .map(|(idx, (gc, damp_slot))| {
                if !gc.is_tumor {
                    return 0;
                }
                if gc.state.dead {
                    if let Some(ds) = gc.state.death_step {
                        let grace_end = ds + params.post_death_steps;
                        if step == grace_end {
                            // `lp_at_death` = "LP at end of post-death grace"
                            // (misnamed; matches sim-tme — rename deferred to
                            // the #195 cleanup checklist).
                            gc.lp_at_death = gc.state.lp;
                            // DAMP release gated on immune_on (else damp_field
                            // is never read/aggregated — PR #219 third-pass).
                            if condition.immune_on {
                                *damp_slot += gc.lp_at_death * immune_cfg.damp_per_lp;
                            }
                        }
                        if step >= grace_end {
                            return 0;
                        }
                        // else: grace period — fall through to sim_cell_step
                        // for post-death LP accumulation.
                    } else {
                        return 0;
                    }
                }

                let mut rng = StdRng::seed_from_u64(
                    cond_seed
                        .wrapping_add(500_000)
                        .wrapping_add(idx as u64)
                        .wrapping_add(step as u64 * 1_000_000),
                );

                let extra_iron = gc.extra_iron;
                gc.extra_iron = 0.0;

                // Time-varying drug modulation (#239). Skipped on the Constant
                // default path (`dosed == false`) → byte-identical there.
                // Live cells only; a dead cell's stored drug state freezes at
                // death (see the longer note removed from here — preserved in
                // git history; the SDT post-death effective-exo nuance is
                // documented on `biochem::exo_decay_factor`).
                if dosed && !gc.state.dead {
                    match condition.treatment {
                        Treatment::RSL3 => {
                            // Covalent GPX4 inactivation ∝ availability
                            // (schedule × pH ion-trap). Mirrors sim_cell_with_pk.
                            let conc = (dose_factor * rsl3_drug_avail[idx]).clamp(0.0, 1.0);
                            gc.state.gpx4 -= RSL3_INACTIVATION_RATE * conc * gc.state.gpx4;
                            gc.state.gpx4 = gc.state.gpx4.max(0.0);
                        }
                        Treatment::SDT | Treatment::PDT => {
                            // Divide out sim_cell_step's intrinsic single-bolus
                            // envelope so the schedule is the sole envelope
                            // (#239). `.max(1e-9)` guards 0/0; the ~1000×
                            // amplification at late steps is finite and never
                            // escapes sim_cell_step's immediate re-multiply.
                            let decay = exo_decay_factor(step).max(1e-9);
                            gc.state.exo_ros_peak = base_exo[idx] * dose_factor / decay;
                        }
                        Treatment::Control => {}
                    }
                }

                let died =
                    sim_cell_step(&mut gc.state, &gc.cell, &params, step, extra_iron, &mut rng);
                if died {
                    // `newly_dead` is consumed by the later serial
                    // `diffuse_iron` to spread released iron to live neighbors.
                    gc.newly_dead = true;
                }

                // Stromal CAF protection for alive cells (stromal_on only).
                if !died && !gc.state.dead {
                    if let Some(cfg) = &stromal_cfg {
                        if adjacency_mask[idx] {
                            gc.state.gsh =
                                (gc.state.gsh + cfg.gsh_boost_per_step).min(cfg.gsh_boost_cap);
                            gc.state.mufa_protection = (gc.state.mufa_protection
                                + cfg.mufa_boost_per_step)
                                .min(cfg.mufa_boost_cap);
                        }
                    }
                }

                usize::from(died)
            })
            .sum();
        ferroptosis_kills += died_this_step;

        // Iron diffusion via TumorGrid3D. Now uses the value from
        // `spatial_params.neighbor_iron_fraction` (single source of truth);
        // the inline `0.031` is gone.
        grid.diffuse_iron(
            spatial_params.iron_release_per_death,
            spatial_params.neighbor_iron_fraction,
        );

        // DAMP diffusion if immune is enabled
        if condition.immune_on {
            diffuse_damp_3d_step(
                &mut damp_field,
                &mut damp_scratch,
                &grid,
                immune_cfg.damp_diffusion_fraction,
                immune_cfg.damp_clearance_rate,
            );

            // Immune kill (after delay). Parallelized over cells with rayon
            // (#192) — byte-identical to the old serial triple loop: each cell
            // reads its own `damp_field[idx]` (immutable here; DAMP diffusion
            // already done) and writes only its own `state.dead`; the per-cell
            // RNG seed `(cond_seed, idx, step)` is position-independent; and
            // `immune_kills` is an order-independent integer sum.
            if step >= IMMUNE_START_STEP {
                let effective_brake = immune_cfg.effective_brake();
                let killed_this_step: usize = grid
                    .cells
                    .par_iter_mut()
                    .enumerate()
                    .map(|(idx, gc)| {
                        if gc.state.dead || !gc.is_tumor {
                            return 0;
                        }
                        let local_damp = damp_field[idx];
                        if local_damp < DAMP_KILL_THRESHOLD {
                            return 0;
                        }
                        let activation = dc_activation(local_damp, immune_cfg.dc_activation_kd);
                        let kill_prob = immune_kill_probability(
                            activation,
                            immune_cfg.immune_kill_rate,
                            effective_brake,
                        );
                        let mut rng = StdRng::seed_from_u64(
                            cond_seed
                                .wrapping_add(900_000_000)
                                .wrapping_add(idx as u64)
                                .wrapping_add(step as u64 * 2_000_000),
                        );
                        if rng.gen::<f64>() < kill_prob {
                            // Immune kills set `state.dead` but deliberately
                            // NOT `death_step`/`newly_dead`: immune-killed cells
                            // are apoptotic (no post-death LP grace, no DAMP
                            // burst, no iron release) — a modeling choice
                            // consistent with sim-tme, not a side-effect.
                            gc.state.dead = true;
                            1
                        } else {
                            0
                        }
                    })
                    .sum();
                immune_kills += killed_this_step;
            }
        }

        // Per-step trajectory capture for `--snapshot` runs. No-op for the
        // default 24-condition matrix path (snapshot is None there).
        // Captured *after* all per-step work (biochem + DAMP + immune kill)
        // so the snapshot reflects end-of-step state.
        if let Some(buf) = snapshot.as_mut() {
            buf.capture_step(&grid, &damp_field);
        }
    }

    // Late DAMP release for cells still in their post-death grace period at
    // the end of the simulation. Iterate paired with damp_field by index to
    // satisfy clippy::needless_range_loop while keeping the dual-Vec access.
    // Run only when immune was on — otherwise damp_field is never read or
    // aggregated and this loop is dead writes (reviewer-flagged).
    if condition.immune_on {
        for (idx, gc) in grid.cells.iter_mut().enumerate() {
            if !(gc.is_tumor && gc.state.dead) {
                continue;
            }
            if let Some(ds) = gc.state.death_step {
                let grace_end = ds + params.post_death_steps;
                if grace_end >= run_cfg.n_steps {
                    gc.lp_at_death = gc.state.lp;
                    damp_field[idx] += gc.lp_at_death * immune_cfg.damp_per_lp;
                }
            }
        }
    }

    // --- Aggregate results ---
    let census = grid.census();
    let overall = census.total_dead as f64 / census.total_tumor.max(1) as f64;
    let (norm_r, trans_r, hyp_r) = zone_kill_rates_3d(&grid, ZONE_REF_LAMBDA);

    // DAMP aggregations — vindicate the PR-body's "DAMP fields" claim and
    // give the comparison script real numbers to work with. When immune
    // was off these are both 0.0 (damp_field was zero-initialized and
    // never written by the diffusion path).
    let peak_damp = damp_field.iter().cloned().fold(0.0_f64, f64::max);
    let total_damp: f64 = damp_field.iter().sum();

    // Stromal adjacency rates — ALWAYS computed (mask is independent of
    // the stromal_on toggle), so the comparison script's Q3 can pair
    // stromal-on adjacent rates with no-stromal baseline adjacent rates.
    // `stromal_mode` still only set when stromal_on is true (matches
    // sim-tme JSON convention).
    let mut adj_dead = 0usize;
    let mut adj_total = 0usize;
    for (idx, gc) in grid.cells.iter().enumerate() {
        if gc.is_tumor && adjacency_mask[idx] {
            adj_total += 1;
            if gc.state.dead {
                adj_dead += 1;
            }
        }
    }
    let stromal_adjacent_kill_rate = if adj_total > 0 {
        Some(adj_dead as f64 / adj_total as f64)
    } else {
        None
    };
    let stromal_adjacent_count = if adj_total > 0 { Some(adj_total) } else { None };
    let stromal_mode = if condition.stromal_on {
        Some("stromal_on".to_string())
    } else {
        None
    };

    let (ph_mode, ph_edge_v, ph_core_v, ph_lambda_v) = if let Some((_, cfg)) = &ph_field {
        (
            Some("ph_on".to_string()),
            Some(cfg.ph_edge),
            Some(cfg.ph_core),
            Some(cfg.lambda_ph_um),
        )
    } else {
        (None, None, None, None)
    };

    ConditionResult {
        treatment: condition.treatment_name.clone(),
        o2_condition: match condition.o2_lambda {
            None => "uniform".to_string(),
            Some(_) => "gradient".to_string(),
        },
        o2_lambda_um: condition.o2_lambda,
        // Match sim-tme's JSON convention exactly (sim-tme uses "immune_on" /
        // "immune_anti_pd1" / "off"). My earlier "on"/"off" prevented the
        // comparison script from finding 2D immune conditions — caught in
        // adversarial review.
        immune_mode: if condition.immune_on {
            "immune_on".to_string()
        } else {
            "off".to_string()
        },
        total_tumor: census.total_tumor,
        total_dead: census.total_dead,
        ferroptosis_kills: if condition.immune_on {
            Some(ferroptosis_kills)
        } else {
            None
        },
        immune_kills: if condition.immune_on {
            Some(immune_kills)
        } else {
            None
        },
        overall_kill_rate: overall,
        normoxic_kill_rate: norm_r,
        transition_kill_rate: trans_r,
        hypoxic_kill_rate: hyp_r,
        peak_damp,
        total_damp,
        stromal_mode,
        stromal_adjacent_kill_rate,
        stromal_adjacent_count,
        ph_mode,
        ph_edge: ph_edge_v,
        ph_core: ph_core_v,
        ph_lambda_um: ph_lambda_v,
    }
}

// Cheap deterministic hash of condition name → seed offset.
fn hash_condition_name(name: &str) -> u64 {
    let mut h: u64 = 0xcbf29ce484222325; // FNV-1a offset
    for b in name.bytes() {
        h ^= b as u64;
        h = h.wrapping_mul(0x100000001b3);
    }
    h
}

// ============================================================
// Condition matrix
// ============================================================

fn generate_conditions() -> Vec<Condition> {
    let mut conditions = Vec::new();
    let treatments = [
        (Treatment::Control, "Control"),
        (Treatment::RSL3, "RSL3"),
        (Treatment::SDT, "SDT"),
    ];

    // Baseline: uniform O₂, no immune/stromal/pH
    for (tx, name) in &treatments {
        conditions.push(Condition {
            name: format!("baseline_{}", name),
            treatment: *tx,
            treatment_name: name.to_string(),
            o2_lambda: None,
            immune_on: false,
            stromal_on: false,
            ph_on: false,
            dose_schedule: DoseSchedule::Constant,
        });
    }

    // O₂ gradient sweep (no immune/stromal/pH yet)
    for &lambda in O2_LAMBDAS {
        for (tx, name) in &treatments {
            conditions.push(Condition {
                name: format!("o2_{}_{}", lambda as i32, name),
                treatment: *tx,
                treatment_name: name.to_string(),
                o2_lambda: Some(lambda),
                immune_on: false,
                stromal_on: false,
                ph_on: false,
                dose_schedule: DoseSchedule::Constant,
            });
        }
    }

    // Immune on (at ZONE_REF_LAMBDA O₂)
    for (tx, name) in &treatments {
        conditions.push(Condition {
            name: format!("immune_{}", name),
            treatment: *tx,
            treatment_name: name.to_string(),
            o2_lambda: Some(ZONE_REF_LAMBDA),
            immune_on: true,
            stromal_on: false,
            ph_on: false,
            dose_schedule: DoseSchedule::Constant,
        });
    }

    // Stromal on + immune on (matches sim-tme's coupling: 2D stromal-on rows
    // all have `immune_mode == "immune_on"`. Standalone stromal-off-immune-off
    // wouldn't exist in 2D, so the comparison script couldn't pair against
    // anything matching — reviewer-flagged confounding in PR #219 review).
    for (tx, name) in &treatments {
        conditions.push(Condition {
            name: format!("stromal_{}", name),
            treatment: *tx,
            treatment_name: name.to_string(),
            o2_lambda: Some(ZONE_REF_LAMBDA),
            immune_on: true,
            stromal_on: true,
            ph_on: false,
            dose_schedule: DoseSchedule::Constant,
        });
    }

    // pH on + immune on (same coupling as stromal — sim-tme's pH-on rows
    // all have `immune_mode == "immune_on"`).
    for (tx, name) in &treatments {
        conditions.push(Condition {
            name: format!("ph_{}", name),
            treatment: *tx,
            treatment_name: name.to_string(),
            o2_lambda: Some(ZONE_REF_LAMBDA),
            immune_on: true,
            stromal_on: false,
            ph_on: true,
            dose_schedule: DoseSchedule::Constant,
        });
    }

    // Combined: immune + stromal + pH (at ZONE_REF_LAMBDA O₂)
    for (tx, name) in &treatments {
        conditions.push(Condition {
            name: format!("combined_{}", name),
            treatment: *tx,
            treatment_name: name.to_string(),
            o2_lambda: Some(ZONE_REF_LAMBDA),
            immune_on: true,
            stromal_on: true,
            ph_on: true,
            dose_schedule: DoseSchedule::Constant,
        });
    }

    conditions
}

// ============================================================
// Main
// ============================================================

/// Parse a `--snapshot[=NAME]` argument from the CLI args iterator.
/// Returns `Some(name)` if the flag is present (defaults to
/// `"combined"` if no `=NAME` suffix), `None` otherwise. Unknown
/// preset names are validated downstream by [`resolve_snapshot`].
fn parse_snapshot_arg<I: IntoIterator<Item = String>>(args: I) -> Option<String> {
    for a in args {
        if a == "--snapshot" {
            return Some("combined".to_string());
        }
        if let Some(name) = a.strip_prefix("--snapshot=") {
            return Some(name.to_string());
        }
    }
    None
}

/// One entry in the `--snapshot=NAME` registry. Each entry pins a
/// reproducible condition for visualization. Names are stable; adding
/// a new variant doesn't change existing ones. See [`SNAPSHOTS`].
struct SnapshotPreset {
    /// CLI name (e.g. `combined`, `bare`). Used as `--snapshot=NAME`.
    name: &'static str,
    /// Short human-readable description for the eprintln banner.
    desc: &'static str,
    /// Treatment modality for this preset.
    treatment: Treatment,
    /// Display name for `treatment` (avoids a Debug-format dependency).
    treatment_name: &'static str,
    /// True if immune kills + DAMP coupling are on (immune_mode = "immune_on").
    immune_on: bool,
    /// True if stromal protection (CAF GSH/MUFA shielding) is on.
    stromal_on: bool,
    /// True if the pH gradient + iron release + ion-trap are on.
    ph_on: bool,
    /// True if this preset uses the multi-dose schedule (#239) instead of
    /// the steady-state `Constant` default. The concrete `DoseSchedule` is
    /// built at runtime in `run_snapshot` (a `const` can't hold the
    /// `MultiDose` Vec).
    multidose: bool,
}

/// Visualization presets for `--snapshot=NAME`. Keep this list small —
/// each entry costs ~333 MB on disk when its trajectory is generated.
const SNAPSHOTS: &[SnapshotPreset] = &[
    SnapshotPreset {
        name: "combined",
        desc: "RSL3 + immune + stromal + pH (all TME protections active)",
        treatment: Treatment::RSL3,
        treatment_name: "RSL3",
        immune_on: true,
        stromal_on: true,
        ph_on: true,
        multidose: false,
    },
    SnapshotPreset {
        name: "bare",
        desc: "RSL3, no immune / stromal / pH (the unprotected baseline)",
        treatment: Treatment::RSL3,
        treatment_name: "RSL3",
        immune_on: false,
        stromal_on: false,
        ph_on: false,
        multidose: false,
    },
    SnapshotPreset {
        name: "multidose",
        desc: "SDT multi-dose (4 pulses) + immune — death waves sync to each dose (#239)",
        treatment: Treatment::SDT,
        treatment_name: "SDT",
        immune_on: true,
        stromal_on: false,
        ph_on: false,
        multidose: true,
    },
];

/// The multi-dose schedule used by the `multidose` snapshot preset:
/// four ROS pulses across the 180-step run. Sharp half-life (8 steps) so
/// each pulse rises and fades distinctly, producing visible death waves.
///
/// These parameters are **illustrative** — chosen for visual clarity of
/// the death-wave dynamics, NOT calibrated to a clinical SDT protocol.
fn multidose_snapshot_schedule() -> DoseSchedule {
    DoseSchedule::MultiDose {
        dose_steps: vec![10, 55, 100, 145],
        peak: 1.0,
        half_life_steps: 8.0,
    }
}

/// Look up a snapshot preset by name, or print available choices and exit.
fn resolve_snapshot(name: &str) -> &'static SnapshotPreset {
    SNAPSHOTS
        .iter()
        .find(|s| s.name == name)
        .unwrap_or_else(|| {
            eprintln!("ERROR: unknown --snapshot name `{name}`. Available presets:");
            for s in SNAPSHOTS {
                eprintln!("  {:<10} — {}", s.name, s.desc);
            }
            std::process::exit(2);
        })
}

/// Run a single condition with per-step trajectory capture, then write
/// `trajectory_{dead,damp,lp}.npy` + `trajectory_meta.json` to
/// `output_dir`. Driven by the `--snapshot[=NAME]` CLI flag (#193, #239).
///
/// Treatment, TME toggles, and dose schedule vary by [`SnapshotPreset`].
/// Files are written under the same names regardless of preset — running
/// a second time overwrites; rename manually to keep several side by side.
fn run_snapshot(output_dir: &Path, tumor_radius_um: f64, name: &str) {
    let preset = resolve_snapshot(name);
    let dose_schedule = if preset.multidose {
        multidose_snapshot_schedule()
    } else {
        DoseSchedule::Constant
    };
    // Administration steps for metadata / renderer dose markers (empty for
    // the steady-state Constant presets).
    let dose_steps = dose_schedule.dose_steps();

    let condition = Condition {
        name: format!("snapshot_{}_{}", preset.name, preset.treatment_name),
        treatment: preset.treatment,
        treatment_name: preset.treatment_name.to_string(),
        o2_lambda: Some(ZONE_REF_LAMBDA),
        immune_on: preset.immune_on,
        stromal_on: preset.stromal_on,
        ph_on: preset.ph_on,
        dose_schedule,
    };

    eprintln!(
        "=== --snapshot={}: {} ({}³ × {} steps) ===",
        preset.name, preset.desc, GRID_DIM, N_STEPS,
    );

    let run_cfg = RunConfig::production();
    let mut buffers = snapshot::SnapshotBuffers::new(run_cfg.grid_dim, run_cfg.n_steps);
    let result = run_one_condition_with_config(&condition, run_cfg, Some(&mut buffers));
    eprintln!(
        "  done — total_kill={:.1}%, ferroptosis_kills={}, immune_kills={}, captured {} steps",
        result.overall_kill_rate * 100.0,
        result.ferroptosis_kills.unwrap_or(0),
        result.immune_kills.unwrap_or(0),
        buffers.steps_captured(),
    );

    buffers
        .write(output_dir)
        .expect("Failed to write trajectory .npy files");

    let meta = snapshot::TrajectoryMeta {
        schema_version: snapshot::TRAJECTORY_SCHEMA_VERSION,
        grid_dim: run_cfg.grid_dim,
        cell_size_um: CELL_SIZE_UM,
        tumor_radius_um,
        n_steps: run_cfg.n_steps,
        dose_steps,
        condition: snapshot::TrajectoryCondition {
            treatment: preset.treatment_name.to_string(),
            o2_condition: "gradient".to_string(),
            o2_lambda_um: Some(ZONE_REF_LAMBDA),
            immune_mode: if preset.immune_on { "immune_on" } else { "off" }.to_string(),
            stromal_mode: preset.stromal_on.then(|| "stromal_on".to_string()),
            ph_mode: preset.ph_on.then(|| "ph_on".to_string()),
        },
    };
    let meta_path = output_dir.join("trajectory_meta.json");
    fs::write(&meta_path, serde_json::to_string_pretty(&meta).unwrap())
        .expect("Failed to write trajectory_meta.json");

    eprintln!(
        "Wrote {} + trajectory_{{dead,damp,lp}}.npy",
        meta_path.display()
    );
    eprintln!("Render with: python3 scripts/render_tme_3d_trajectory.py");
}

/// Schema version for `dose_comparison.json` (#239). Independent of
/// `summary.json`'s schema.
const DOSE_SWEEP_SCHEMA_VERSION: u32 = 1;

/// One row of the dose-sweep comparison: RSL3 kill outcome under a single
/// dosing protocol, all else equal.
#[derive(Serialize)]
struct DoseSweepEntry {
    /// Short protocol label (constant / bolus / multidose / infusion / frompk).
    schedule: String,
    /// Human-readable protocol description.
    description: String,
    total_tumor: usize,
    total_dead: usize,
    overall_kill_rate: f64,
    ferroptosis_kills: usize,
    immune_kills: usize,
}

#[derive(Serialize)]
struct DoseSweepResult {
    schema_version: u32,
    treatment: String,
    /// The fixed biological context shared by every protocol.
    context: String,
    grid_dim: usize,
    n_steps: u32,
    /// One entry per dosing protocol, all sharing the same grid + RNG seed
    /// so differences reflect the protocol alone, not stochastic noise.
    schedules: Vec<DoseSweepEntry>,
}

/// Normalized RSL3 interstitial-concentration factor series from the
/// two-compartment tumor PK ODE (#239 `DoseSchedule::FromPk` bridge).
///
/// Solves `tumor_pk::solve_tumor_pk` for an IV bolus into a breast-tumor
/// compartment, then normalizes the interstitial timecourse to peak 1.0 so
/// it reads as a drug-availability factor. This is what finally wires the
/// (previously orphaned) PK ODE into the spatial grid, without coupling the
/// grid to the solver — the grid just consumes the resulting `&[f64]`.
fn rsl3_pk_factor_series(n_steps: u32) -> Vec<f64> {
    use ferroptosis_core::tumor_pk::{breast_tumor, solve_tumor_pk, PlasmaModel};
    // IV bolus, ~35-step plasma half-life (k_el = ln2 / 35 ≈ 0.0198 /min).
    let plasma = PlasmaModel::IvBolus {
        c0: 1.0,
        k_el: 0.0198,
    };
    let tumor = breast_tumor();
    let res = solve_tumor_pk(&plasma, &tumor, n_steps as usize, 10);
    let max = res
        .c_interstitial
        .iter()
        .cloned()
        .fold(0.0_f64, f64::max)
        .max(1e-9);
    res.c_interstitial
        .iter()
        .map(|c| (c / max).clamp(0.0, 1.0))
        .collect()
}

/// `--dose-sweep`: run RSL3 across all five dosing protocols at the
/// combined-TME context (immune + stromal + pH, λ=120) and write
/// `dose_comparison.json` (#239).
///
/// **Controlled comparison**: every protocol uses the SAME tumor grid
/// (fixed `SEED`) and the SAME runtime RNG seed (all conditions share one
/// name), so the only thing that varies between rows is the dosing
/// schedule. Differences in kill rate therefore reflect the protocol, not
/// stochastic noise.
fn run_dose_sweep(output_dir: &Path) {
    eprintln!(
        "=== --dose-sweep: RSL3 across dosing protocols ({}³ × {} steps) ===",
        GRID_DIM, N_STEPS
    );

    // (label, description, schedule). All run RSL3 + immune + stromal + pH.
    //
    // The half-life / peak / level numbers below are ILLUSTRATIVE v1 values,
    // not calibrated to clinical protocols (RSL3_INACTIVATION_RATE itself was
    // tuned for sustained conc=1.0). The informative output is the
    // cross-protocol ORDERING of kill rates, not the absolute magnitudes —
    // dose_comparison.json's `context` field and the README record this.
    let protocols: Vec<(&str, &str, DoseSchedule)> = vec![
        (
            "constant",
            "Steady-state full availability (one-shot GPX4 knockdown)",
            DoseSchedule::Constant,
        ),
        (
            "bolus",
            "Single bolus at step 10, 20-step half-life",
            DoseSchedule::Bolus {
                dose_step: 10,
                peak: 1.0,
                half_life_steps: 20.0,
            },
        ),
        (
            "multidose",
            "4 doses at steps 10/55/100/145, 18-step half-life",
            DoseSchedule::MultiDose {
                dose_steps: vec![10, 55, 100, 145],
                peak: 1.0,
                half_life_steps: 18.0,
            },
        ),
        (
            "infusion",
            "Continuous infusion, level 0.5, steps 10-170",
            DoseSchedule::Infusion {
                start: 10,
                end: 170,
                level: 0.5,
            },
        ),
        (
            "frompk",
            "tumor_pk two-compartment IV-bolus interstitial curve (breast tumor)",
            DoseSchedule::FromPk {
                series: rsl3_pk_factor_series(N_STEPS),
            },
        ),
    ];

    let entries: Vec<DoseSweepEntry> = protocols
        .par_iter()
        .map(|(label, desc, sched)| {
            let cond = Condition {
                // Shared name → shared cond_seed → matched RNG across protocols.
                name: "dosesweep_RSL3".to_string(),
                treatment: Treatment::RSL3,
                treatment_name: "RSL3".to_string(),
                o2_lambda: Some(ZONE_REF_LAMBDA),
                immune_on: true,
                stromal_on: true,
                ph_on: true,
                dose_schedule: sched.clone(),
            };
            let r = run_one_condition(&cond);
            eprintln!(
                "  {label:<10} → kill={:.2}% (ferro={}, immune={})",
                r.overall_kill_rate * 100.0,
                r.ferroptosis_kills.unwrap_or(0),
                r.immune_kills.unwrap_or(0),
            );
            DoseSweepEntry {
                schedule: (*label).to_string(),
                description: (*desc).to_string(),
                total_tumor: r.total_tumor,
                total_dead: r.total_dead,
                overall_kill_rate: r.overall_kill_rate,
                ferroptosis_kills: r.ferroptosis_kills.unwrap_or(0),
                immune_kills: r.immune_kills.unwrap_or(0),
            }
        })
        .collect();

    let result = DoseSweepResult {
        schema_version: DOSE_SWEEP_SCHEMA_VERSION,
        treatment: "RSL3".to_string(),
        context: "immune + stromal + pH at λ=120 µm; shared grid + RNG seed across protocols"
            .to_string(),
        grid_dim: GRID_DIM,
        n_steps: N_STEPS,
        schedules: entries,
    };

    let path = output_dir.join("dose_comparison.json");
    fs::write(&path, serde_json::to_string_pretty(&result).unwrap())
        .expect("Failed to write dose_comparison.json");
    eprintln!("Wrote {}", path.display());
}

/// Read a positive `usize` env var, falling back to `default` if unset or
/// unparseable. Used by `--bench` to override grid size without touching the
/// byte-identical-locked `GRID_DIM` constant.
fn bench_env_usize(name: &str, default: usize) -> usize {
    std::env::var(name)
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .filter(|&n| n > 0)
        .unwrap_or(default)
}

/// Read a positive `u32` env var (parsed directly into `u32` — values above
/// `u32::MAX` are rejected, not silently truncated). For `BENCH_N_STEPS`.
fn bench_env_u32(name: &str, default: u32) -> u32 {
    std::env::var(name)
        .ok()
        .and_then(|v| v.parse::<u32>().ok())
        .filter(|&n| n > 0)
        .unwrap_or(default)
}

/// `--bench`: run ONE representative condition (the combined-TME RSL3 cell —
/// immune + stromal + pH, the heaviest per-cell path) at a configurable grid
/// size and print wall-clock timing. The performance/scalability harness for
/// issue #192.
///
/// Grid size and step count come from env (`BENCH_GRID_DIM`, `BENCH_N_STEPS`;
/// defaults 60/180) so a sweep is scriptable, e.g.
/// `BENCH_GRID_DIM=200 BENCH_N_STEPS=180 cargo run --release -p sim-tme-3d -- --bench`.
/// Capture peak RSS externally: `/usr/bin/time -v <that command>`.
///
/// Runs a SINGLE condition with no condition-level rayon, so it measures the
/// within-condition cost — exactly the work the within-condition parallelism
/// targets. Writes no output files; the default matrix path is untouched.
fn run_bench() {
    let grid_dim = bench_env_usize("BENCH_GRID_DIM", GRID_DIM);
    let n_steps = bench_env_u32("BENCH_N_STEPS", N_STEPS);
    // saturating_pow: a huge BENCH_GRID_DIM saturates rather than overflow-
    // panicking in debug (the grid allocation would fail first anyway).
    let n_cells = grid_dim.saturating_pow(3);

    eprintln!(
        "=== --bench: single combined-TME RSL3 condition, {grid_dim}³ = {n_cells} cells × {n_steps} steps ==="
    );
    eprintln!(
        "(override via BENCH_GRID_DIM / BENCH_N_STEPS; wrap with `/usr/bin/time -v` for peak RSS)"
    );

    let condition = Condition {
        name: "bench_combined_RSL3".to_string(),
        treatment: Treatment::RSL3,
        treatment_name: "RSL3".to_string(),
        o2_lambda: Some(ZONE_REF_LAMBDA),
        immune_on: true,
        stromal_on: true,
        ph_on: true,
        dose_schedule: DoseSchedule::Constant,
    };
    let cfg = RunConfig { grid_dim, n_steps };

    let t0 = Instant::now();
    let r = run_one_condition_with_config(&condition, cfg, None);
    let elapsed = t0.elapsed();

    let secs = elapsed.as_secs_f64();
    let cell_steps = (n_cells as f64) * (n_steps as f64);
    eprintln!(
        "  done in {secs:.2}s — total_kill={:.2}% (ferro={}, immune={})",
        r.overall_kill_rate * 100.0,
        r.ferroptosis_kills.unwrap_or(0),
        r.immune_kills.unwrap_or(0),
    );
    // Machine-readable line for sweep collection (grep 'BENCH_RESULT').
    eprintln!(
        "BENCH_RESULT grid_dim={grid_dim} n_cells={n_cells} n_steps={n_steps} \
         wall_s={secs:.3} cell_steps_per_s={:.3e}",
        cell_steps / secs.max(1e-9)
    );
}

fn main() {
    // Guard against silent drift between this binary's metadata const and
    // the library's runtime value: if a future PR tunes
    // `SpatialImmuneConfig::for_3d().damp_diffusion_fraction` without also
    // touching `DAMP_DIFFUSION_FRACTION_3D`, `summary.json` would report a
    // stale value while the simulation actually used the new one. This
    // debug_assert catches that in tests/dev runs without affecting release
    // bit-identical output.
    debug_assert_eq!(
        DAMP_DIFFUSION_FRACTION_3D,
        SpatialImmuneConfig::for_3d().damp_diffusion_fraction,
        "DAMP_DIFFUSION_FRACTION_3D metadata const out of sync with SpatialImmuneConfig::for_3d()",
    );

    // Use the library's TUMOR_RADIUS_FRACTION constant instead of the bare
    // 0.45 literal — keeps `tumor_radius_um` in lockstep with whatever
    // value `TumorGrid3D::generate` actually uses (reviewer-flagged drift
    // hazard if the library ever tunes the fraction).
    let tumor_radius_um = (GRID_DIM as f64) * TUMOR_RADIUS_FRACTION * CELL_SIZE_UM;
    eprintln!("=== sim-tme-3d: 3D Spheroid TME Simulation ===");
    eprintln!(
        "Grid: {0}³ ({1:.1} mm × {1:.1} mm × {1:.1} mm)",
        GRID_DIM,
        GRID_DIM as f64 * CELL_SIZE_UM / 1000.0
    );
    eprintln!("Tumor radius: {:.0} µm", tumor_radius_um);
    eprintln!(
        "O₂ λ sweep: {:?} µm (λ=150 skipped — 3λ > tumor radius)",
        O2_LAMBDAS
    );
    eprintln!(
        "DAMP diffusion fraction: {} (3D-safe; 2D's 0.08 would trigger immune_spatial stability assert!)",
        DAMP_DIFFUSION_FRACTION_3D
    );
    eprintln!();
    eprintln!(
        "⚠️  Scale mismatch with sim-tme (2D): tumor radius {:.0} µm vs sim-tme's ~4500 µm.",
        tumor_radius_um
    );
    eprintln!("    Compare via RATIOS, not absolute counts.");
    eprintln!();

    let output_dir = Path::new("output/tme-3d");
    fs::create_dir_all(output_dir).expect("Failed to create output/tme-3d");

    // `--snapshot[=NAME]` runs ONE visualization-focused condition with
    // per-step state capture for the Python animation. Default path (no
    // flag) is unchanged: runs the full 24-condition matrix →
    // summary.json. Bit-identical when the flag is absent (snapshot path
    // doesn't touch the matrix). Names listed in `SNAPSHOTS`; bare
    // `--snapshot` defaults to `combined`.
    if let Some(name) = parse_snapshot_arg(std::env::args()) {
        run_snapshot(output_dir, tumor_radius_um, &name);
        return;
    }

    // `--dose-sweep` runs RSL3 across all five dosing protocols (Constant /
    // Bolus / MultiDose / Infusion / FromPk) at the combined-TME context and
    // writes `dose_comparison.json` — the quantitative "does dosing protocol
    // change efficacy?" answer (#239). Separate file; does NOT touch
    // summary.json (which stays byte-identical to the 24-condition matrix).
    if std::env::args().any(|a| a == "--dose-sweep") {
        run_dose_sweep(output_dir);
        return;
    }

    // `--bench` runs ONE representative condition at a configurable grid size
    // (env BENCH_GRID_DIM / BENCH_N_STEPS, defaults 60/180) and prints
    // wall-clock timing — the performance/scalability harness for #192. It
    // measures the WITHIN-condition cost (no condition-level rayon), which is
    // exactly what the within-condition parallelism targets. Does NOT write
    // summary.json; the default matrix path below is untouched (byte-identical).
    if std::env::args().any(|a| a == "--bench") {
        run_bench();
        return;
    }

    let conditions = generate_conditions();
    eprintln!(
        "Running {} conditions in parallel via rayon...",
        conditions.len()
    );

    // Condition-level parallelism. Each run_one_condition allocates its own
    // grid + DAMP/scratch buffers; no shared mutable state.
    let results: Vec<ConditionResult> = conditions
        .par_iter()
        .map(|cond| {
            eprintln!("  starting: {}", cond.name);
            let r = run_one_condition(cond);
            eprintln!(
                "  done:     {} — total_kill={:.1}% (norm={:.1}%, trans={:.1}%, hyp={:.1}%)",
                cond.name,
                r.overall_kill_rate * 100.0,
                r.normoxic_kill_rate * 100.0,
                r.transition_kill_rate * 100.0,
                r.hypoxic_kill_rate * 100.0
            );
            r
        })
        .collect();

    let summary = SimulationSummary {
        schema_version: TME_3D_SCHEMA_VERSION,
        grid_dim: GRID_DIM,
        cell_size_um: CELL_SIZE_UM,
        tumor_radius_um,
        n_steps: N_STEPS,
        o2_lambdas: O2_LAMBDAS.to_vec(),
        damp_diffusion_fraction_3d: DAMP_DIFFUSION_FRACTION_3D,
        conditions: results,
        note: format!(
            "sim-tme-3d uses {}³ grid (tumor radius {:.0} µm) vs sim-tme's 500×500 \
             (tumor radius ~4500 µm). Compare conditions via RATIOS, not absolute counts. \
             See generate_3d_comparison_table.py.",
            GRID_DIM, tumor_radius_um
        ),
    };

    let json_path = output_dir.join("summary.json");
    let json = serde_json::to_string_pretty(&summary).expect("serialize summary");
    fs::write(&json_path, json).expect("write summary.json");
    eprintln!("\nWrote {}", json_path.display());
    eprintln!("Done.");
}

// ============================================================
// Tests (smoke-level — primitives are already tested in
// ferroptosis-core)
// ============================================================

#[cfg(test)]
mod tests {
    use super::*;

    // Golden integer kill counts for the Constant-path regression guard
    // (`constant_path_golden_kill_counts`), captured from SDT + immune +
    // stromal + pH at grid_dim=20, n_steps=80. Deterministic on a fixed
    // platform/toolchain (seeded RNG + IEEE f64); may differ by an ULP-edge
    // count on other libm/CPU builds — see the test's doc. Update ONLY after
    // confirming a default-path change is intentional.
    const GOLDEN_TOTAL_TUMOR: usize = 3071;
    const GOLDEN_TOTAL_DEAD: usize = 2992;
    const GOLDEN_FERRO_KILLS: usize = 2990;
    const GOLDEN_IMMUNE_KILLS: usize = 2;
    // Dosed-path golden (MultiDose SDT + immune + stromal + pH, 20³×80).
    const GOLDEN_DOSED_TOTAL_DEAD: usize = 2660;
    const GOLDEN_DOSED_FERRO_KILLS: usize = 2658;
    const GOLDEN_DOSED_IMMUNE_KILLS: usize = 2;

    /// `--snapshot` (no `=NAME`) defaults to `combined`.
    #[test]
    fn parse_snapshot_arg_bare_defaults_to_combined() {
        let args = vec!["sim-tme-3d".to_string(), "--snapshot".to_string()];
        assert_eq!(parse_snapshot_arg(args), Some("combined".to_string()));
    }

    /// `--snapshot=bare` extracts the name after `=`.
    #[test]
    fn parse_snapshot_arg_equals_form_extracts_name() {
        let args = vec!["sim-tme-3d".to_string(), "--snapshot=bare".to_string()];
        assert_eq!(parse_snapshot_arg(args), Some("bare".to_string()));
    }

    /// Absent flag returns None (the default 24-condition matrix path).
    #[test]
    fn parse_snapshot_arg_absent_returns_none() {
        let args = vec!["sim-tme-3d".to_string()];
        assert_eq!(parse_snapshot_arg(args), None);
    }

    /// All `SNAPSHOTS` entries have unique names (used as the `=NAME`
    /// match key downstream). A typo or copy-paste duplicate would
    /// make `resolve_snapshot` ambiguous (first-match-wins).
    #[test]
    fn snapshot_preset_names_are_unique() {
        let mut names: Vec<&str> = SNAPSHOTS.iter().map(|s| s.name).collect();
        names.sort();
        let original = names.len();
        names.dedup();
        assert_eq!(names.len(), original, "duplicate snapshot preset name");
    }

    /// Smoke test: condition matrix is non-empty and well-formed.
    #[test]
    fn condition_matrix_is_non_empty() {
        let conditions = generate_conditions();
        assert!(!conditions.is_empty(), "expected at least some conditions");
        // Names should be unique (used as seed-hash inputs).
        let mut names: Vec<&str> = conditions.iter().map(|c| c.name.as_str()).collect();
        names.sort();
        let original_count = names.len();
        names.dedup();
        assert_eq!(
            names.len(),
            original_count,
            "condition names must be unique"
        );
    }

    /// Smoke test: one baseline condition (Control, no toggles) runs
    /// end-to-end and produces sensible output.
    ///
    /// Uses `RunConfig::for_test()` (10³ × 20 steps) instead of production
    /// (60³ × 180 steps) — reviewer flagged that `cargo test --workspace`
    /// runs in DEBUG mode, where the full 60³ × 180 sim was costing minutes
    /// per test. The smaller config exercises every code path at ~1000×
    /// less cost.
    #[test]
    fn single_condition_runs_end_to_end() {
        let cond = Condition {
            name: "test_baseline_Control".to_string(),
            treatment: Treatment::Control,
            treatment_name: "Control".to_string(),
            o2_lambda: None,
            immune_on: false,
            stromal_on: false,
            ph_on: false,
            dose_schedule: DoseSchedule::Constant,
        };
        let r = run_one_condition_with_config(&cond, RunConfig::for_test(), None);
        assert_eq!(r.treatment, "Control");
        assert_eq!(r.o2_condition, "uniform");
        assert_eq!(r.immune_mode, "off");
        assert!(r.total_tumor > 0, "expected some tumor cells in 10³ grid");
        // Control has zero exo-ROS, no immune pressure, no pH stress. Production
        // 60³ run gives ~0.15% baseline kill rate. Threshold tightened from
        // <5% (too loose — would let a 10× baseline regression sail through)
        // to <2% (still loose for safety margin given the small test grid;
        // production-rate × 5 reasonable margin).
        assert!(
            r.overall_kill_rate < 0.02,
            "Control should have <2% kill rate, got {:.1}%",
            r.overall_kill_rate * 100.0
        );
    }

    /// Determinism: same condition seed → identical output.
    #[test]
    fn same_seed_same_output() {
        let cond = Condition {
            name: "deterministic_test".to_string(),
            treatment: Treatment::Control,
            treatment_name: "Control".to_string(),
            o2_lambda: None,
            immune_on: false,
            stromal_on: false,
            ph_on: false,
            dose_schedule: DoseSchedule::Constant,
        };
        let r1 = run_one_condition_with_config(&cond, RunConfig::for_test(), None);
        let r2 = run_one_condition_with_config(&cond, RunConfig::for_test(), None);
        assert_eq!(r1.total_dead, r2.total_dead);
        assert_eq!(r1.total_tumor, r2.total_tumor);
        assert_eq!(r1.overall_kill_rate, r2.overall_kill_rate);
    }

    /// **Reviewer-flagged invariant guard**: immune-on conditions with a
    /// treatment that produces ferroptotic kills (RSL3 or SDT) must
    /// produce at least one immune kill. If the immune block silently
    /// no-ops (e.g., DAMP threshold never crossed, or a bug in the
    /// activation chain), `immune_kills` would still be 0 and the JSON
    /// would still serialize — bug only surfaces in the manuscript.
    ///
    /// Uses RSL3 + O₂ gradient at λ=120 (the canonical "immune on"
    /// condition) and asserts `immune_kills > 0` over 20 steps. SDT
    /// also kills heavily so picking RSL3 makes the assertion stricter
    /// (RSL3 produces fewer ferroptotic deaths → fewer DAMP sources →
    /// harder for any immune kills to fire spuriously).
    #[test]
    fn immune_on_actually_fires() {
        let cond = Condition {
            name: "test_immune_RSL3".to_string(),
            treatment: Treatment::RSL3,
            treatment_name: "RSL3".to_string(),
            o2_lambda: Some(ZONE_REF_LAMBDA),
            immune_on: true,
            stromal_on: false,
            ph_on: false,
            dose_schedule: DoseSchedule::Constant,
        };
        // Larger config than the smoke tests: need IMMUNE_START_STEP=60 +
        // buffer for DAMP to accumulate above the 0.01 kill threshold.
        // 25³ × 130 steps: ~1500 tumor cells, 70 steps post-immune-start
        // for DAMP to spread and reach kill threshold. Empirically reliable
        // (production 60³×180 produces 29 RSL3 immune kills; volume scales
        // 60³/25³ ≈ 14×, so expect ~2 kills at this config — comfortable
        // margin over the >0 invariant).
        let cfg = RunConfig {
            grid_dim: 25,
            n_steps: 130,
        };
        let r = run_one_condition_with_config(&cond, cfg, None);
        assert_eq!(r.immune_mode, "immune_on");
        let im_kills = r
            .immune_kills
            .expect("immune_on must populate immune_kills");
        // Strict invariant: at least one immune kill in this RSL3 + immune-on
        // run. If this fails the immune block has gone no-op.
        assert!(
            im_kills > 0,
            "immune-on RSL3 should produce ≥1 immune kill in {} steps on {grid}³; got {im_kills}. \
             Likely cause: DAMP threshold never crossed, or activation chain broken.",
            cfg.n_steps,
            grid = cfg.grid_dim
        );
    }

    /// **Reviewer-flagged invariant guard**: for immune-on conditions,
    /// `total_dead` must equal `ferroptosis_kills + immune_kills`.
    /// Catches double-counting drift (e.g., a cell counted in both
    /// branches) and missed-kill drift (e.g., a kill that mutates
    /// `state.dead` but doesn't increment either counter).
    #[test]
    fn total_dead_equals_ferro_plus_immune() {
        let cond = Condition {
            name: "test_invariant_RSL3".to_string(),
            treatment: Treatment::RSL3,
            treatment_name: "RSL3".to_string(),
            o2_lambda: Some(ZONE_REF_LAMBDA),
            immune_on: true,
            stromal_on: false,
            ph_on: false,
            dose_schedule: DoseSchedule::Constant,
        };
        let cfg = RunConfig {
            grid_dim: 15,
            n_steps: 80,
        };
        let r = run_one_condition_with_config(&cond, cfg, None);
        let ferro = r.ferroptosis_kills.expect("immune_on populates this");
        let im = r.immune_kills.expect("immune_on populates this");
        assert_eq!(
            r.total_dead,
            ferro + im,
            "total_dead {} != ferroptosis_kills {} + immune_kills {} — kill accounting drifted",
            r.total_dead,
            ferro,
            im
        );
    }

    // ============================================================
    // Time-varying dose schedule (#239)
    // ============================================================

    /// A zero-availability schedule (a bolus scheduled past the end of the
    /// run, so `factor_at ≡ 0`) must suppress SDT kills relative to the
    /// Constant full-availability default — proving the SDT exo-ROS rescale
    /// path is actually wired in.
    #[test]
    fn dosed_zero_availability_schedule_suppresses_sdt_kills() {
        let cfg = RunConfig {
            grid_dim: 20,
            n_steps: 60,
        };
        let constant = run_one_condition_with_config(
            &Condition {
                name: "test_sdt_constant".to_string(),
                treatment: Treatment::SDT,
                treatment_name: "SDT".to_string(),
                o2_lambda: Some(ZONE_REF_LAMBDA),
                immune_on: false,
                stromal_on: false,
                ph_on: false,
                dose_schedule: DoseSchedule::Constant,
            },
            cfg,
            None,
        );
        let zero = run_one_condition_with_config(
            &Condition {
                name: "test_sdt_zero_drug".to_string(),
                treatment: Treatment::SDT,
                treatment_name: "SDT".to_string(),
                o2_lambda: Some(ZONE_REF_LAMBDA),
                immune_on: false,
                stromal_on: false,
                ph_on: false,
                // Dose never arrives within the run → factor_at ≡ 0 → no SDT ROS.
                dose_schedule: DoseSchedule::Bolus {
                    dose_step: 9999,
                    peak: 1.0,
                    half_life_steps: 10.0,
                },
            },
            cfg,
            None,
        );
        assert!(
            constant.total_dead > 0,
            "Constant SDT should kill cells at 20³×60; got {}",
            constant.total_dead
        );
        assert!(
            zero.total_dead < constant.total_dead,
            "zero-availability schedule must suppress SDT kills: zero={}, constant={}",
            zero.total_dead,
            constant.total_dead
        );
    }

    /// A gentle single-bolus RSL3 schedule (continuous covalent inactivation)
    /// must kill fewer cells than the Constant one-shot 92% GPX4 knockdown —
    /// proving the RSL3 dosed path is live and the no-init-knockdown +
    /// per-step inactivation mechanism differs from the steady-state default.
    #[test]
    fn dosed_rsl3_bolus_kills_less_than_constant_oneshot() {
        let cfg = RunConfig {
            grid_dim: 20,
            n_steps: 60,
        };
        let constant = run_one_condition_with_config(
            &Condition {
                name: "test_rsl3_constant".to_string(),
                treatment: Treatment::RSL3,
                treatment_name: "RSL3".to_string(),
                o2_lambda: Some(ZONE_REF_LAMBDA),
                immune_on: false,
                stromal_on: false,
                ph_on: false,
                dose_schedule: DoseSchedule::Constant,
            },
            cfg,
            None,
        );
        let bolus = run_one_condition_with_config(
            &Condition {
                name: "test_rsl3_bolus".to_string(),
                treatment: Treatment::RSL3,
                treatment_name: "RSL3".to_string(),
                o2_lambda: Some(ZONE_REF_LAMBDA),
                immune_on: false,
                stromal_on: false,
                ph_on: false,
                dose_schedule: DoseSchedule::Bolus {
                    dose_step: 2,
                    peak: 1.0,
                    half_life_steps: 8.0,
                },
            },
            cfg,
            None,
        );
        assert!(
            constant.total_dead > 0,
            "Constant RSL3 (one-shot 92% knockdown) should kill cells; got {}",
            constant.total_dead
        );
        assert!(
            bolus.total_dead < constant.total_dead,
            "a single decaying RSL3 bolus (gentle continuous inactivation) must kill \
             fewer than the Constant one-shot: bolus={}, constant={}",
            bolus.total_dead,
            constant.total_dead
        );
    }

    /// **Byte-identity regression guard for the Constant default path.**
    /// Pins the exact integer kill counts of a representative Constant
    /// condition (SDT + immune + stromal + pH) at a fixed small config,
    /// catching ANY drift in the default (non-dosed) path — the load-bearing
    /// "summary.json byte-identical" property — at unit-test speed, without
    /// the full 60³×180 production run.
    ///
    /// **Scope: deterministic on a fixed platform + toolchain.** The kill
    /// counts are integers, but the death decision is a strict
    /// `lp > threshold` on values flowing through `powf` and the Ziggurat
    /// normal sampler, which are not bit-identical across libm/CPU
    /// implementations — so a cell within an ULP of the threshold could flip
    /// a count on a different platform. (This is a pre-existing
    /// whole-simulation property, not introduced by #239.) This test is the
    /// same-platform CI guard against accidental drift; the production
    /// SHA-256 check (done manually in #239's PR) is the cross-build
    /// authority. If these numbers change, the byte-identical claim is broken
    /// — investigate before updating.
    #[test]
    fn constant_path_golden_kill_counts() {
        let cond = Condition {
            name: "golden_constant_SDT".to_string(),
            treatment: Treatment::SDT,
            treatment_name: "SDT".to_string(),
            o2_lambda: Some(ZONE_REF_LAMBDA),
            immune_on: true,
            stromal_on: true,
            ph_on: true,
            dose_schedule: DoseSchedule::Constant,
        };
        let cfg = RunConfig {
            grid_dim: 20,
            n_steps: 80,
        };
        let r = run_one_condition_with_config(&cond, cfg, None);
        // Golden values captured from the Constant path. A change here means
        // the default-path numerics drifted — investigate before updating.
        assert_eq!(
            r.total_tumor, GOLDEN_TOTAL_TUMOR,
            "tumor-cell count drifted"
        );
        assert_eq!(
            r.total_dead, GOLDEN_TOTAL_DEAD,
            "Constant-path total_dead drifted"
        );
        assert_eq!(
            r.ferroptosis_kills,
            Some(GOLDEN_FERRO_KILLS),
            "Constant-path ferroptosis_kills drifted"
        );
        assert_eq!(
            r.immune_kills,
            Some(GOLDEN_IMMUNE_KILLS),
            "Constant-path immune_kills drifted"
        );
    }

    /// **Golden guard for the DOSED parallel path** (review #255 substantive
    /// point 2). `constant_path_golden_kill_counts` pins the Constant path;
    /// this pins a MultiDose SDT + immune + stromal + pH condition, which
    /// exercises the machinery the Constant golden does NOT: the per-step
    /// SDT exo-envelope divide-out, the grace-end DAMP write, and the
    /// immune-kill loop — all through the rayon-parallelized loops. Same
    /// fixed-platform/toolchain caveat as the Constant golden (integer
    /// counts; `powf`/Ziggurat aren't bit-identical across libm).
    #[test]
    fn dosed_path_golden_kill_counts() {
        let cond = Condition {
            name: "golden_multidose_SDT".to_string(),
            treatment: Treatment::SDT,
            treatment_name: "SDT".to_string(),
            o2_lambda: Some(ZONE_REF_LAMBDA),
            immune_on: true,
            stromal_on: true,
            ph_on: true,
            dose_schedule: DoseSchedule::MultiDose {
                dose_steps: vec![5, 30, 55],
                peak: 1.0,
                half_life_steps: 8.0,
            },
        };
        let cfg = RunConfig {
            grid_dim: 20,
            n_steps: 80,
        };
        let r = run_one_condition_with_config(&cond, cfg, None);
        assert_eq!(
            r.total_tumor, GOLDEN_TOTAL_TUMOR,
            "tumor-cell count drifted"
        );
        assert_eq!(
            r.total_dead, GOLDEN_DOSED_TOTAL_DEAD,
            "dosed-path total_dead drifted"
        );
        assert_eq!(
            r.ferroptosis_kills,
            Some(GOLDEN_DOSED_FERRO_KILLS),
            "dosed-path ferroptosis_kills drifted"
        );
        assert_eq!(
            r.immune_kills,
            Some(GOLDEN_DOSED_IMMUNE_KILLS),
            "dosed-path immune_kills drifted"
        );
    }

    /// The full dosed stack (RSL3 + immune + stromal + pH + MultiDose) must
    /// run end-to-end without panicking and be deterministic across runs.
    /// Exercises `rsl3_drug_avail` indexing (pH on), the no-knockdown init,
    /// per-step inactivation, and the immune/stromal interactions together.
    #[test]
    fn dosed_rsl3_full_stack_runs_and_is_deterministic() {
        let cond = Condition {
            name: "test_rsl3_multidose_full".to_string(),
            treatment: Treatment::RSL3,
            treatment_name: "RSL3".to_string(),
            o2_lambda: Some(ZONE_REF_LAMBDA),
            immune_on: true,
            stromal_on: true,
            ph_on: true,
            dose_schedule: DoseSchedule::MultiDose {
                dose_steps: vec![2, 10, 18],
                peak: 1.0,
                half_life_steps: 6.0,
            },
        };
        let cfg = RunConfig {
            grid_dim: 20,
            n_steps: 30,
        };
        let r1 = run_one_condition_with_config(&cond, cfg, None);
        let r2 = run_one_condition_with_config(&cond, cfg, None);
        assert_eq!(
            r1.total_dead, r2.total_dead,
            "dosed full-stack run must be deterministic"
        );
        assert!(r1.total_tumor > 0, "expected tumor cells");
        assert!(
            (0.0..=1.0).contains(&r1.overall_kill_rate),
            "kill rate must be a valid fraction, got {}",
            r1.overall_kill_rate
        );
    }

    /// End-to-end coverage for the `Infusion` variant (review #10): runs
    /// RSL3 under a continuous infusion through the full per-step dosed
    /// path and asserts it's deterministic and produces a valid fraction.
    #[test]
    fn dosed_rsl3_infusion_runs_and_is_deterministic() {
        let cond = Condition {
            name: "test_rsl3_infusion".to_string(),
            treatment: Treatment::RSL3,
            treatment_name: "RSL3".to_string(),
            o2_lambda: Some(ZONE_REF_LAMBDA),
            immune_on: true,
            stromal_on: false,
            ph_on: true,
            dose_schedule: DoseSchedule::Infusion {
                start: 2,
                end: 28,
                level: 1.0,
            },
        };
        let cfg = RunConfig {
            grid_dim: 20,
            n_steps: 30,
        };
        let r1 = run_one_condition_with_config(&cond, cfg, None);
        let r2 = run_one_condition_with_config(&cond, cfg, None);
        assert_eq!(
            r1.total_dead, r2.total_dead,
            "Infusion run must be deterministic"
        );
        assert!(r1.total_tumor > 0);
        assert!((0.0..=1.0).contains(&r1.overall_kill_rate));
    }

    /// End-to-end coverage for the `FromPk` variant (review #10): runs RSL3
    /// driven by the tumor_pk-derived availability series through the full
    /// per-step dosed path. Exercises the ODE-bridge → grid wiring beyond
    /// the `factor_at`-level unit test.
    #[test]
    fn dosed_rsl3_frompk_runs_and_is_deterministic() {
        let cfg = RunConfig {
            grid_dim: 20,
            n_steps: 30,
        };
        let cond = Condition {
            name: "test_rsl3_frompk".to_string(),
            treatment: Treatment::RSL3,
            treatment_name: "RSL3".to_string(),
            o2_lambda: Some(ZONE_REF_LAMBDA),
            immune_on: true,
            stromal_on: false,
            ph_on: true,
            dose_schedule: DoseSchedule::FromPk {
                series: rsl3_pk_factor_series(cfg.n_steps),
            },
        };
        let r1 = run_one_condition_with_config(&cond, cfg, None);
        let r2 = run_one_condition_with_config(&cond, cfg, None);
        assert_eq!(
            r1.total_dead, r2.total_dead,
            "FromPk run must be deterministic"
        );
        assert!(r1.total_tumor > 0);
        assert!((0.0..=1.0).contains(&r1.overall_kill_rate));
    }

    /// The FromPk bridge series must be a valid normalized availability
    /// factor: correct length, all values in [0, 1], and peak ≈ 1.0
    /// (normalization target). Guards the `tumor_pk` → spatial-grid wiring.
    #[test]
    fn rsl3_pk_factor_series_is_normalized() {
        let series = rsl3_pk_factor_series(180);
        assert_eq!(series.len(), 180, "series length must match n_steps");
        assert!(
            series.iter().all(|&f| (0.0..=1.0).contains(&f)),
            "all factors must be in [0, 1]"
        );
        let peak = series.iter().cloned().fold(0.0_f64, f64::max);
        assert!(
            (peak - 1.0).abs() < 1e-9,
            "series must be normalized to peak 1.0, got {peak}"
        );
    }

    /// Empirical NaN check on `norm()` per reviewer ask: confirms that a
    /// NaN `std` propagates to a panic (via `Normal::new(...).expect(...)`),
    /// not a silent NaN value. Documents the failure mode rather than
    /// making it "safe".
    #[test]
    #[should_panic(expected = "valid normal")]
    fn norm_panics_on_nan_std() {
        let mut rng = StdRng::seed_from_u64(0);
        norm(&mut rng, 1.0, f64::NAN);
    }

    /// **Reviewer-flagged invariant**: par_iter must produce identical
    /// results across runs. Per-condition closures are independent (each
    /// allocates its own grid + DAMP/scratch buffers), so determinism
    /// should hold — but an unverified invariant is one rayon refactor
    /// away from silent drift. This test exercises the actual
    /// `par_iter().map(run_one_condition)` path (smaller config for CI
    /// debug speed) and compares the resulting JSON.
    #[test]
    fn rayon_run_is_deterministic() {
        // Subset of 6 conditions (small enough for debug mode, large
        // enough to exercise inter-thread scheduling).
        let small_conditions: Vec<_> = generate_conditions().into_iter().take(6).collect();
        let cfg = RunConfig::for_test();

        let run = |conds: &[Condition]| -> Vec<ConditionResult> {
            conds
                .par_iter()
                .map(|c| run_one_condition_with_config(c, cfg, None))
                .collect()
        };

        let r1 = run(&small_conditions);
        let r2 = run(&small_conditions);

        assert_eq!(r1.len(), r2.len());
        for (a, b) in r1.iter().zip(r2.iter()) {
            assert_eq!(a.treatment, b.treatment);
            assert_eq!(a.total_dead, b.total_dead);
            assert_eq!(a.overall_kill_rate, b.overall_kill_rate);
            assert_eq!(a.peak_damp, b.peak_damp);
            assert_eq!(a.total_damp, b.total_damp);
        }
    }

    /// **Reviewer-flagged invariant**: the library functions
    /// `radial_o2_field`, `radial_ph_field`, and `stromal_adjacency_mask_3d`
    /// all return `Vec<T>` of length `grid.cells.len()`, and this binary
    /// indexes those vecs by `grid.flat_index(r, c, l)`. The implicit
    /// contract is that the library iterates in the SAME row-major order
    /// as the grid's flat-index formula. If a library refactor ever
    /// changes the iteration order (e.g., column-major instead of
    /// row-major), every cell would silently get the wrong O₂/pH/stromal
    /// value and nothing in the smoke tests would catch it.
    ///
    /// This test pins the contract with two spot-checks:
    /// 1. Vec lengths match `grid.cells.len()`.
    /// 2. A known cell has expected qualitative properties:
    ///    - The center cell has higher O₂ depth (less O₂) than a corner cell
    ///    - The center cell has lower pH than a corner cell
    ///    - The center cell is NOT in the stromal-adjacency mask (interior),
    ///      while a near-surface tumor cell IS.
    ///
    /// A column-major refactor would scramble these properties and fail.
    #[test]
    fn library_field_order_matches_flat_index() {
        use ferroptosis_core::oxygen::radial_o2_field;
        use ferroptosis_core::ph::radial_ph_field;

        let grid = TumorGrid3D::generate(20, 20, 20, CELL_SIZE_UM, SEED);
        let n = grid.cells.len();

        let o2 = radial_o2_field(&grid, 100.0);
        let ph = radial_ph_field(&grid, 7.4, 6.5, 120.0);
        let mask = stromal_adjacency_mask_3d(&grid);
        assert_eq!(o2.len(), n, "radial_o2_field length contract");
        assert_eq!(ph.len(), n, "radial_ph_field length contract");
        assert_eq!(mask.len(), n, "stromal_adjacency_mask_3d length contract");

        // Index the same cells via the binary's `flat_index` access path:
        let center_idx = grid.flat_index(10, 10, 10);
        let surface_idx = grid.flat_index(10, 10, 19); // depth=0 surface (per grid::tests_3d)

        // Center cell is interior tumor → high pH gradient depth → lower pH
        // than surface; lower O₂ than surface; mask=false (no stromal neighbors).
        assert!(
            grid.cells[center_idx].is_tumor,
            "test precondition: center is tumor"
        );
        assert!(
            grid.cells[surface_idx].is_tumor,
            "test precondition: surface is tumor"
        );
        assert!(
            ph[center_idx] < ph[surface_idx],
            "center pH {} should be lower than surface pH {} — \
             if this fails, library field order is desync'd from flat_index",
            ph[center_idx],
            ph[surface_idx]
        );
        assert!(
            o2[center_idx] < o2[surface_idx],
            "center O₂ {} should be lower than surface O₂ {} (deeper = more hypoxic)",
            o2[center_idx],
            o2[surface_idx]
        );
        // Surface cell at (10,10,19) is at the spheroid edge → has stromal
        // neighbors at (10,10,20)? No — (10,10,20) is OUTSIDE the grid bounds
        // for a 20³ grid. Surface cell's neighbors include some that are
        // outside the sphere → mask = true.
        assert!(
            mask[surface_idx],
            "surface tumor cell should be in adjacency mask (has stromal neighbors)"
        );
        assert!(
            !mask[center_idx],
            "interior center cell should NOT be in adjacency mask"
        );
    }
}
