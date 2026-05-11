//! sim-tme-3d: 3D spheroid tumor microenvironment simulation.
//!
//! Capstone 3D binary for the spheroid-validation series (#185–#197).
//! Integrates all five library primitives landed in v0.7.0–v0.11.0:
//! - 3D energy physics (#186) via `physics::local_ros_multiplier_3d`
//! - 3D radial O₂ gradient (#187) via `oxygen::radial_o2_field`
//! - 3D radial pH gradient (#190) via `ph::radial_ph_field` + helpers
//! - 3D CAF-shielded boundary detection (#189) via `stromal::stromal_adjacency_mask`
//! - 3D spatial DAMP diffusion + activation (#188) via `immune_3d::*`
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
//! 2D's per-step total diffusion of ~64% — see `immune_3d` rustdoc).
//! `immune_3d::diffuse_damp_3d_step` enforces the stability invariant
//! with `assert!` (release-mode panic).

use std::fs;
use std::path::Path;

use ferroptosis_core::biochem::{sim_cell_step, CellState};
use ferroptosis_core::cell::Treatment;
use ferroptosis_core::grid::{GridCensus, TumorGrid3D};
use ferroptosis_core::immune_3d::{
    dc_activation, diffuse_damp_3d_step, immune_kill_probability,
};
use ferroptosis_core::oxygen::radial_o2_field;
use ferroptosis_core::params::{Params, SpatialParams};
use ferroptosis_core::ph::{ion_trap_factor_from_ph, iron_multiplier_from_ph, radial_ph_field};
use ferroptosis_core::physics::local_ros_multiplier_3d;
use ferroptosis_core::stromal::stromal_adjacency_mask;
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

/// DAMP diffusion fraction. **3D-safe value**: sim-tme's 2D default
/// (0.08) is unsafe in 3D — `immune_3d::diffuse_damp_3d_step`'s
/// stability `assert!` would panic. Use 0.025 to match 2D's per-step
/// total diffusion (~64%).
const DAMP_DIFFUSION_FRACTION_3D: f64 = 0.025;

/// Step at which immune activity begins (matches sim-tme).
const IMMUNE_START_STEP: u32 = 60;

// ============================================================
// Config structs (duplicated from sim-tme for v1 — lift to
// ferroptosis-core::params is a follow-up cleanup PR per the
// consolidated checklist on issue #195).
// ============================================================

/// Parameters for the spatial immune model. Mirror of
/// `sim-tme::ImmuneConfig` — the duplication is documented technical
/// debt; see issue #195 follow-up tracking.
#[derive(Clone, Copy, Debug)]
struct ImmuneConfig {
    damp_per_lp: f64,
    /// 3D-safe diffusion fraction (NOT sim-tme's 2D 0.08 — that
    /// triggers the immune_3d stability `assert!`).
    damp_diffusion_fraction: f64,
    damp_clearance_rate: f64,
    dc_activation_kd: f64,
    immune_kill_rate: f64,
    pd1_brake: f64,
    anti_pd1_efficacy: f64,
}

impl ImmuneConfig {
    fn default_no_pd1() -> Self {
        ImmuneConfig {
            damp_per_lp: 1.0,
            damp_diffusion_fraction: DAMP_DIFFUSION_FRACTION_3D,
            damp_clearance_rate: 0.03,
            dc_activation_kd: 50.0,
            immune_kill_rate: 0.02,
            pd1_brake: 0.7,
            anti_pd1_efficacy: 0.0,
        }
    }

    fn effective_brake(&self) -> f64 {
        self.pd1_brake * (1.0 - self.anti_pd1_efficacy)
    }
}

/// CAF-mediated stromal protection params. Mirror of
/// `sim-tme::StromalConfig`.
#[derive(Clone, Copy, Debug)]
struct StromalConfig {
    gsh_boost_per_step: f64,
    gsh_boost_cap: f64,
    mufa_boost_per_step: f64,
    mufa_boost_cap: f64,
}

impl StromalConfig {
    fn default() -> Self {
        StromalConfig {
            gsh_boost_per_step: 0.06,
            gsh_boost_cap: 18.0,
            mufa_boost_per_step: 0.003,
            mufa_boost_cap: 0.25,
        }
    }
}

/// pH gradient + iron / ion-trap sensitivity params. Mirror of
/// `sim-tme::PhConfig`.
#[derive(Clone, Copy, Debug)]
struct PhConfig {
    ph_edge: f64,
    ph_core: f64,
    lambda_ph_um: f64,
    iron_ph_sensitivity: f64,
    ion_trap_sensitivity: f64,
}

impl PhConfig {
    fn default() -> Self {
        PhConfig {
            ph_edge: 7.4,
            ph_core: 6.5,
            lambda_ph_um: 120.0,
            iron_ph_sensitivity: 1.5,
            ion_trap_sensitivity: 0.4,
        }
    }
}

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
    #[serde(skip_serializing_if = "Option::is_none")]
    stromal_mode: Option<String>,
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

#[derive(Clone, Debug, Serialize)]
struct SimulationSummary {
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

fn run_one_condition(condition: &Condition) -> ConditionResult {
    let params = Params::default();
    let spatial_params = SpatialParams {
        cell_size_um: CELL_SIZE_UM,
        ..Default::default()
    };

    // Per-condition deterministic seed.
    let cond_seed = SEED.wrapping_add(hash_condition_name(&condition.name));

    let mut grid = TumorGrid3D::generate(GRID_DIM, GRID_DIM, GRID_DIM, CELL_SIZE_UM, cond_seed);
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

    // --- Compute stromal adjacency mask if needed ---
    let stromal_mask = if condition.stromal_on {
        Some((stromal_adjacency_mask(&grid), StromalConfig::default()))
    } else {
        None
    };

    // --- Treatment-specific ROS multiplier (3D version) ---
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
                let exo_ros_peak = if matches!(
                    condition.treatment,
                    Treatment::Control | Treatment::RSL3
                ) {
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
                gc.state =
                    CellState::from_cell_with_ros(&gc.cell, condition.treatment, &params, exo_ros_peak);
                gc.extra_iron = 0.0;
                gc.newly_dead = false;
                gc.lp_at_death = 0.0;
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
            let correction = (1.0 - params.rsl3_gpx4_inhib * drug_factor)
                / (1.0 - params.rsl3_gpx4_inhib);
            grid.cells[idx].state.gpx4 *= correction;
        }
    }

    // --- Main 180-step loop ---
    let immune_cfg = ImmuneConfig::default_no_pd1();
    let mut damp_field = vec![0.0_f64; n_cells];
    let mut damp_scratch = vec![0.0_f64; n_cells];
    let mut ferroptosis_kills = 0usize;
    let mut immune_kills = 0usize;

    for step in 0..N_STEPS {
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
                                grid.cells[idx].lp_at_death = grid.cells[idx].state.lp;
                                damp_field[idx] +=
                                    grid.cells[idx].lp_at_death * immune_cfg.damp_per_lp;
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
                        gc.newly_dead = true;
                        gc.lp_at_death = gc.state.lp;
                        ferroptosis_kills += 1;
                    }

                    // Stromal CAF protection for alive cells
                    if !died && !gc.state.dead {
                        if let Some((mask, cfg)) = &stromal_mask {
                            if mask[idx] {
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

        // Iron diffusion via TumorGrid3D
        grid.diffuse_iron(
            spatial_params.iron_release_per_death,
            // 3D-natural neighbor fraction (per #185 grid rustdoc):
            // 0.1 * 8/26 ≈ 0.031
            0.031,
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
                            if local_damp < 0.01 {
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
    for (idx, gc) in grid.cells.iter_mut().enumerate() {
        if !(gc.is_tumor && gc.state.dead) {
            continue;
        }
        if let Some(ds) = gc.state.death_step {
            let grace_end = ds + params.post_death_steps;
            if grace_end >= N_STEPS {
                gc.lp_at_death = gc.state.lp;
                damp_field[idx] += gc.lp_at_death * immune_cfg.damp_per_lp;
            }
        }
    }

    // --- Aggregate results ---
    let census = censuses_for_grid(&grid);
    let overall = census.total_dead as f64 / census.total_tumor.max(1) as f64;
    let (norm_r, trans_r, hyp_r) = zone_kill_rates_3d(&grid, ZONE_REF_LAMBDA);

    let (stromal_mode, stromal_kill, stromal_count) = if let Some((mask, _)) = &stromal_mask {
        let mut dead = 0usize;
        let mut total = 0usize;
        for (idx, gc) in grid.cells.iter().enumerate() {
            if gc.is_tumor && mask[idx] {
                total += 1;
                if gc.state.dead {
                    dead += 1;
                }
            }
        }
        let rate = if total > 0 {
            dead as f64 / total as f64
        } else {
            0.0
        };
        (Some("stromal_on".to_string()), Some(rate), Some(total))
    } else {
        (None, None, None)
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
        immune_mode: if condition.immune_on {
            "on".to_string()
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
        stromal_mode,
        stromal_adjacent_kill_rate: stromal_kill,
        stromal_adjacent_count: stromal_count,
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

fn censuses_for_grid(grid: &TumorGrid3D) -> GridCensus {
    grid.census()
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

    // Stromal on (at ZONE_REF_LAMBDA O₂)
    for (tx, name) in &treatments {
        conditions.push(Condition {
            name: format!("stromal_{}", name),
            treatment: *tx,
            treatment_name: name.to_string(),
            o2_lambda: Some(ZONE_REF_LAMBDA),
            immune_on: false,
            stromal_on: true,
            ph_on: false,
        });
    }

    // pH on (at ZONE_REF_LAMBDA O₂)
    for (tx, name) in &treatments {
        conditions.push(Condition {
            name: format!("ph_{}", name),
            treatment: *tx,
            treatment_name: name.to_string(),
            o2_lambda: Some(ZONE_REF_LAMBDA),
            immune_on: false,
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
    let tumor_radius_um = (GRID_DIM as f64) * 0.45 * CELL_SIZE_UM;
    eprintln!("=== sim-tme-3d: 3D Spheroid TME Simulation ===");
    eprintln!(
        "Grid: {0}³ ({1:.1} mm × {1:.1} mm × {1:.1} mm)",
        GRID_DIM,
        GRID_DIM as f64 * CELL_SIZE_UM / 1000.0
    );
    eprintln!("Tumor radius: {:.0} µm", tumor_radius_um);
    eprintln!("O₂ λ sweep: {:?} µm (λ=150 skipped — 3λ > tumor radius)", O2_LAMBDAS);
    eprintln!(
        "DAMP diffusion fraction: {} (3D-safe; 2D's 0.08 would trigger immune_3d stability assert!)",
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
    eprintln!("Running {} conditions in parallel via rayon...", conditions.len());

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
    /// end-to-end and produces sensible output. Validates the full
    /// orchestration without paying for the 30-condition matrix.
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
        let r = run_one_condition(&cond);
        assert_eq!(r.treatment, "Control");
        assert_eq!(r.o2_condition, "uniform");
        assert_eq!(r.immune_mode, "off");
        assert!(r.total_tumor > 0, "expected some tumor cells in 60³ grid");
        // Control has zero exo-ROS so kill rate should be very low (baseline only).
        assert!(
            r.overall_kill_rate < 0.05,
            "Control should have <5% kill rate, got {:.1}%",
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
        let r1 = run_one_condition(&cond);
        let r2 = run_one_condition(&cond);
        assert_eq!(r1.total_dead, r2.total_dead);
        assert_eq!(r1.total_tumor, r2.total_tumor);
        assert_eq!(r1.overall_kill_rate, r2.overall_kill_rate);
    }
}
