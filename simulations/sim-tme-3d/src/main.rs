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
//! upper end of in-vitro spheroids). Larger 3D grids are infeasible
//! at this stage (a 500³ grid would need ~21 GB; even 100³ × 180 steps
//! is hours per condition without #194's perf optimizations).
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

mod npy;

use ferroptosis_core::biochem::{sim_cell_step, CellState};
use ferroptosis_core::cell::Treatment;
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
    run_one_condition_with_config(condition, RunConfig::production())
}

fn run_one_condition_with_config(condition: &Condition, run_cfg: RunConfig) -> ConditionResult {
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
                let gc = &mut grid.cells[idx];
                gc.state = CellState::from_cell_with_ros(
                    &gc.cell,
                    condition.treatment,
                    &params,
                    exo_ros_peak,
                );
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

    // pH-dependent RSL3 ion trapping correction (consumer-side, same pattern as sim-tme)
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

    // --- Main 180-step loop ---
    let immune_cfg = SpatialImmuneConfig::for_3d();
    let mut damp_field = vec![0.0_f64; n_cells];
    let mut damp_scratch = vec![0.0_f64; n_cells];
    let mut ferroptosis_kills = 0usize;
    let mut immune_kills = 0usize;

    for step in 0..run_cfg.n_steps {
        // Ferroptosis biochem + stromal protection
        for r in 0..grid.rows {
            for c in 0..grid.cols {
                for l in 0..grid.layers {
                    let idx = grid.flat_index(r, c, l);
                    if !grid.cells[idx].is_tumor {
                        continue;
                    }
                    if grid.cells[idx].state.dead {
                        if let Some(ds) = grid.cells[idx].state.death_step {
                            let grace_end = ds + params.post_death_steps;
                            if step == grace_end {
                                // `lp_at_death` is misleadingly named — it's
                                // actually "LP at the end of the post-death
                                // grace period." Matches sim-tme's naming;
                                // renaming would touch both binaries and is
                                // deferred per the consolidated #195 cleanup
                                // checklist.
                                grid.cells[idx].lp_at_death = grid.cells[idx].state.lp;
                                // DAMP release gated on `immune_on` —
                                // without immune coupling, the damp_field
                                // isn't read or aggregated, so writes here
                                // would produce nonsense partial values
                                // in serialized output (reviewer-flagged in
                                // PR #219 third-pass).
                                if condition.immune_on {
                                    damp_field[idx] +=
                                        grid.cells[idx].lp_at_death * immune_cfg.damp_per_lp;
                                }
                            }
                            if step >= grace_end {
                                continue;
                            }
                        } else {
                            continue;
                        }
                    }

                    let mut rng = StdRng::seed_from_u64(
                        cond_seed
                            .wrapping_add(500_000)
                            .wrapping_add(idx as u64)
                            .wrapping_add(step as u64 * 1_000_000),
                    );

                    let extra_iron = grid.cells[idx].extra_iron;
                    grid.cells[idx].extra_iron = 0.0;

                    let gc = &mut grid.cells[idx];
                    let died =
                        sim_cell_step(&mut gc.state, &gc.cell, &params, step, extra_iron, &mut rng);

                    if died {
                        // `newly_dead` is consumed by `TumorGrid3D::diffuse_iron`
                        // later this step to spread released iron to live
                        // 26-Moore neighbors (grid.rs:586). Not vestigial —
                        // load-bearing for iron-driven Fenton bystander effects.
                        gc.newly_dead = true;
                        ferroptosis_kills += 1;
                        // `lp_at_death` is set later, at grace-end, just before
                        // being read for the DAMP-field write. The moment-of-death
                        // write that used to live here was dead (always overwritten
                        // before read); removed in #220.
                    }

                    // Stromal CAF protection for alive cells — apply only when
                    // stromal_on. The adjacency_mask is always computed (used
                    // for kill-rate accounting later) but the boost helpers
                    // gate on the toggle.
                    if !died && !gc.state.dead {
                        if let Some(cfg) = &stromal_cfg {
                            if adjacency_mask[idx] {
                                let gc = &mut grid.cells[idx];
                                gc.state.gsh =
                                    (gc.state.gsh + cfg.gsh_boost_per_step).min(cfg.gsh_boost_cap);
                                gc.state.mufa_protection = (gc.state.mufa_protection
                                    + cfg.mufa_boost_per_step)
                                    .min(cfg.mufa_boost_cap);
                            }
                        }
                    }
                }
            }
        }

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

            // Immune kill (after delay)
            if step >= IMMUNE_START_STEP {
                let effective_brake = immune_cfg.effective_brake();
                for r in 0..grid.rows {
                    for c in 0..grid.cols {
                        for l in 0..grid.layers {
                            let idx = grid.flat_index(r, c, l);
                            if grid.cells[idx].state.dead || !grid.cells[idx].is_tumor {
                                continue;
                            }
                            let local_damp = damp_field[idx];
                            if local_damp < DAMP_KILL_THRESHOLD {
                                continue;
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
                                // do NOT set `death_step` or `newly_dead`.
                                // **Modeling asymmetry**: ferroptosis-killed
                                // cells with `Some(death_step)` fall through
                                // to `sim_cell_step` for `post_death_steps`
                                // more iterations (post-death LP overshoot
                                // dynamics); immune-killed cells with
                                // `None` death_step hit the `else { continue }`
                                // branch and exit immediately. Three
                                // intentional consequences:
                                //   (1) immune-killed cells skip grace-period
                                //       biochem — they're apoptotic, not
                                //       ferroptotic (no LP buildup beyond
                                //       moment of death);
                                //   (2) no DAMP burst — T-cell perforin
                                //       death doesn't release the DAMPs
                                //       (calreticulin, HMGB1) that ferroptosis
                                //       releases;
                                //   (3) no iron release via `diffuse_iron`
                                //       (no `newly_dead` flag).
                                // This is a real modeling choice (consistent
                                // with sim-tme/main.rs:749), not a side-effect.
                                grid.cells[idx].state.dead = true;
                                immune_kills += 1;
                            }
                        }
                    }
                }
            }
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
        });
    }

    conditions
}

// ============================================================
// Main
// ============================================================

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
        };
        let r = run_one_condition_with_config(&cond, RunConfig::for_test());
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
        };
        let r1 = run_one_condition_with_config(&cond, RunConfig::for_test());
        let r2 = run_one_condition_with_config(&cond, RunConfig::for_test());
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
        let r = run_one_condition_with_config(&cond, cfg);
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
        };
        let cfg = RunConfig {
            grid_dim: 15,
            n_steps: 80,
        };
        let r = run_one_condition_with_config(&cond, cfg);
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
                .map(|c| run_one_condition_with_config(c, cfg))
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
