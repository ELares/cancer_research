//! Tumor Microenvironment Simulation
//!
//! Models two TME features:
//! A) Oxygen gradients — hypoxia protects tumor core from ferroptosis
//! B) Spatial immune zones — DAMP release from ferroptotic cells triggers
//!    local immune activation (ICD-to-kill coupling)
//!
//! Biology (cross-referenced against textbooks in books/ directory):
//! - O2 → basal ROS via mitochondrial ETC (Biology2e Ch.7-8)
//! - O2 penetration ~100-150μm (Anatomy & Physiology 2e, Vaupel 1989)
//! - DAMP release from ferroptotic cells triggers innate immunity
//!   (Biology2e Ch.42-43, Microbiology Ch.15-17, Krysko Nat Rev Cancer 2012)
//! - DC activation follows Michaelis-Menten kinetics (Chemistry2e Ch.12-14)
//! - Spatial immune model valid for resident T cell phase (0-48h);
//!   systemic lymph node priming (1-7 days) is NOT modeled
//!
//! Key finding: LP at death is ~10.0 for ALL treatments (threshold-locked),
//! so DAMP per dead cell is approximately equal. The immune differential
//! comes from kill DENSITY (SDT kills 88% = dense DAMP field) not DAMP
//! quality (which is similar across treatments). This is an honest finding
//! that corrects the initial hypothesis about ICD quality differences.
//!
//! Caveats:
//! - O2 modulates basal_ros only (conservative)
//! - SDT/PDT modeled as O2-independent (Type I mechanism, conservative)
//! - LP at death ~10.0 underestimates true DAMP quality differential by
//!   ~30-50% (biologically, SDT should drive LP to 15-20 post-threshold)
//! - DAMP clearance modeled as exponential decay (simplified)
//! - Immune kill is local/resident phase only (no systemic priming)
//!
//! Usage: `cargo run --release --bin sim-tme`

use std::fs;
use std::path::Path;

use rand::prelude::*;
use serde::Serialize;

use ferroptosis_core::biochem::{sim_cell_step, CellState};
use ferroptosis_core::cell::{norm, Treatment};
use ferroptosis_core::grid::{depth_kill_curve, death_heatmap, TumorGrid};
use ferroptosis_core::io::{write_depth_curves_csv, write_heatmap_csv, write_json};
use ferroptosis_core::params::{Params, SpatialParams};
use ferroptosis_core::physics::local_ros_multiplier;

const GRID_SIZE: usize = 500;
const CELL_SIZE_UM: f64 = 20.0;
const N_STEPS: u32 = 180;
const SEED: u64 = 42;

/// O2 penetration lengths to sweep (μm).
/// Literature range: 100-150μm (Vaupel, Cancer Res 1989).
const O2_LAMBDAS: &[f64] = &[80.0, 100.0, 120.0, 150.0];

/// Fixed reference λ for zone boundary definitions (μm).
/// Zone kill rates use this CONSTANT regardless of the sweep λ, so the
/// sensitivity table compares a fixed anatomical region across conditions.
const ZONE_REF_LAMBDA: f64 = 120.0;

// ============================================================
// O2 field computation
// ============================================================

/// Compute steady-state O2 concentration for each tumor cell based on
/// distance from the tumor edge. Stromal cells (outside tumor) are
/// well-oxygenated; O2 decays exponentially into the tumor interior.
///
/// O2(d) = exp(-d / λ) where d = distance from tumor edge in μm.
///
/// Modifies `cell.basal_ros` in place: `basal_ros *= o2_factor`.
/// Returns a Vec of (row, col, o2_factor) for heatmap generation.
fn apply_o2_gradient(
    grid: &mut TumorGrid,
    penetration_um: f64,
) -> Vec<(usize, usize, f64)> {
    let rows = grid.rows;
    let cols = grid.cols;
    let cell_size = grid.cell_size_um;
    let center_r = rows as f64 / 2.0;
    let center_c = cols as f64 / 2.0;
    let tumor_radius = (rows.min(cols) as f64) * 0.45;

    let mut o2_map = Vec::with_capacity(rows * cols);

    for r in 0..rows {
        for c in 0..cols {
            let gc = grid.get_mut(r, c);
            if !gc.is_tumor {
                o2_map.push((r, c, 1.0));
                continue;
            }

            let dist_from_center =
                ((r as f64 - center_r).powi(2) + (c as f64 - center_c).powi(2)).sqrt();
            let depth_from_edge_um =
                (tumor_radius - dist_from_center).max(0.0) * cell_size;

            let o2_factor = (-depth_from_edge_um / penetration_um).exp().clamp(0.0, 1.0);

            gc.cell.basal_ros *= o2_factor;
            o2_map.push((r, c, o2_factor));
        }
    }

    o2_map
}

// ============================================================
// Spatial simulation (reuses sim-spatial's pattern)
// ============================================================

fn run_spatial(
    grid: &mut TumorGrid,
    tx: Treatment,
    params: &Params,
    spatial_params: &SpatialParams,
    seed: u64,
) {
    let base_ros = match tx {
        Treatment::SDT => params.sdt_ros,
        Treatment::PDT => params.pdt_ros,
        Treatment::RSL3 | Treatment::Control => 0.0,
    };

    let rows = grid.rows;
    let cols = grid.cols;
    let cell_size = grid.cell_size_um;

    // Initialize cell states with depth-dependent ROS
    for r in 0..rows {
        let ros_multiplier = local_ros_multiplier(r, cell_size, tx, spatial_params);
        for c in 0..cols {
            let exo_ros_peak = if tx == Treatment::Control || tx == Treatment::RSL3 {
                0.0
            } else {
                let mut rng = StdRng::seed_from_u64(seed.wrapping_add((r * cols + c) as u64));
                let peak = base_ros * ros_multiplier;
                norm(&mut rng, peak, peak * 0.2).max(0.0)
            };

            let gc = grid.get_mut(r, c);
            gc.state = CellState::from_cell_with_ros(&gc.cell, tx, params, exo_ros_peak);
            gc.extra_iron = 0.0;
            gc.newly_dead = false;
            gc.lp_at_death = 0.0;
        }
    }

    // 180-step loop with interleaved iron diffusion
    for step in 0..N_STEPS {
        for r in 0..rows {
            for c in 0..cols {
                let idx = r * cols + c;
                if grid.cells[idx].state.dead || !grid.cells[idx].is_tumor {
                    continue;
                }

                let mut rng = StdRng::seed_from_u64(
                    seed.wrapping_add(500_000)
                        .wrapping_add(idx as u64)
                        .wrapping_add(step as u64 * 1_000_000),
                );

                let extra_iron = grid.cells[idx].extra_iron;
                grid.cells[idx].extra_iron = 0.0;

                let gc = &mut grid.cells[idx];
                let died = sim_cell_step(
                    &mut gc.state,
                    &gc.cell,
                    params,
                    step,
                    extra_iron,
                    &mut rng,
                );

                if died {
                    gc.newly_dead = true;
                    gc.lp_at_death = gc.state.lp;
                }
            }
        }

        grid.diffuse_iron(
            spatial_params.iron_release_per_death,
            spatial_params.neighbor_iron_fraction,
        );
    }
}

// ============================================================
// Spatial immune coupling (Feature B)
// ============================================================

/// Parameters for the spatial immune model.
struct ImmuneConfig {
    /// DAMP release per unit LP at death (from ImmuneParams default).
    damp_per_lp: f64,
    /// Fraction of DAMP shared with each Moore neighbor per step.
    damp_diffusion_fraction: f64,
    /// Exponential decay rate per step (models immune clearance of DAMPs).
    damp_clearance_rate: f64,
    /// Michaelis-Menten Kd for DC activation by DAMP concentration.
    dc_activation_kd: f64,
    /// Per-step immune kill rate (absorbs DC maturation + T cell priming + kill).
    immune_kill_rate: f64,
    /// PD-1 suppression fraction (0.0 = no brake, 1.0 = full suppression).
    pd1_brake: f64,
    /// Anti-PD-1 efficacy (fraction of brake removed).
    anti_pd1_efficacy: f64,
    /// LP overshoot multiplier for physical modalities (SDT/PDT).
    /// Estimates the post-threshold LP cascade: LP reaches ~2× threshold
    /// for high-ROS treatments (Biology2e Ch.7-8: autocatalytic propagation).
    physical_modality_overshoot: f64,
    /// LP overshoot multiplier for pharmacologic treatments (RSL3, Control).
    /// Minimal momentum past threshold for slow LP accumulation.
    pharmacologic_overshoot: f64,
}

impl ImmuneConfig {
    fn default_no_pd1() -> Self {
        ImmuneConfig {
            damp_per_lp: 1.0,
            damp_diffusion_fraction: 0.08,
            damp_clearance_rate: 0.03,
            dc_activation_kd: 50.0,
            immune_kill_rate: 0.02,
            pd1_brake: 0.7,
            anti_pd1_efficacy: 0.0,
            physical_modality_overshoot: 2.0,
            pharmacologic_overshoot: 1.05,
        }
    }

    fn with_anti_pd1(&self) -> Self {
        ImmuneConfig {
            anti_pd1_efficacy: 0.8,
            ..*self
        }
    }

    fn effective_brake(&self) -> f64 {
        self.pd1_brake * (1.0 - self.anti_pd1_efficacy)
    }
}

/// Run spatial sim WITH immune coupling: DAMP diffusion + immune kill.
/// Returns (ferroptosis_kills, immune_kills, final_damp_field).
fn run_spatial_with_immune(
    grid: &mut TumorGrid,
    tx: Treatment,
    params: &Params,
    spatial_params: &SpatialParams,
    immune: &ImmuneConfig,
    seed: u64,
) -> (usize, usize, Vec<f64>) {
    let base_ros = match tx {
        Treatment::SDT => params.sdt_ros,
        Treatment::PDT => params.pdt_ros,
        Treatment::RSL3 | Treatment::Control => 0.0,
    };

    let rows = grid.rows;
    let cols = grid.cols;
    let cell_size = grid.cell_size_um;
    let n_cells = rows * cols;

    // Initialize cell states
    for r in 0..rows {
        let ros_multiplier = local_ros_multiplier(r, cell_size, tx, spatial_params);
        for c in 0..cols {
            let exo_ros_peak = if tx == Treatment::Control || tx == Treatment::RSL3 {
                0.0
            } else {
                let mut rng = StdRng::seed_from_u64(seed.wrapping_add((r * cols + c) as u64));
                let peak = base_ros * ros_multiplier;
                norm(&mut rng, peak, peak * 0.2).max(0.0)
            };
            let gc = grid.get_mut(r, c);
            gc.state = CellState::from_cell_with_ros(&gc.cell, tx, params, exo_ros_peak);
            gc.extra_iron = 0.0;
            gc.newly_dead = false;
            gc.lp_at_death = 0.0;
        }
    }

    // External DAMP field (not in GridCell — zero ferroptosis-core changes)
    let mut damp_field = vec![0.0_f64; n_cells];
    let mut damp_delta = vec![0.0_f64; n_cells]; // reused each step to avoid allocation churn
    let mut ferroptosis_kills = 0usize;
    let mut immune_kills = 0usize;
    let immune_start_step = 60_u32; // immune activation delay

    for step in 0..N_STEPS {
        // --- Ferroptosis biochemistry ---
        for r in 0..rows {
            for c in 0..cols {
                let idx = r * cols + c;
                if grid.cells[idx].state.dead || !grid.cells[idx].is_tumor {
                    continue;
                }

                let mut rng = StdRng::seed_from_u64(
                    seed.wrapping_add(500_000)
                        .wrapping_add(idx as u64)
                        .wrapping_add(step as u64 * 1_000_000),
                );

                let extra_iron = grid.cells[idx].extra_iron;
                grid.cells[idx].extra_iron = 0.0;

                let gc = &mut grid.cells[idx];
                let died = sim_cell_step(
                    &mut gc.state, &gc.cell, params, step, extra_iron, &mut rng,
                );

                if died {
                    gc.newly_dead = true;
                    gc.lp_at_death = gc.state.lp;
                    ferroptosis_kills += 1;
                    // Release DAMPs into the field
                    // LP overshoot: biologically, the autocatalytic LP cascade
                    // continues 1-3 steps post-threshold for high-ROS treatments
                    // (Biology2e Ch.7-8: chain reaction propagation). SDT/PDT drive
                    // LP to ~20 (2× threshold) while RSL3 barely exceeds ~10.5.
                    // This is an estimated multiplier (Option C from issue #82);
                    // Option A (emergent overshoot from dynamics) is a follow-up.
                    let overshoot = match tx {
                        Treatment::SDT | Treatment::PDT => immune.physical_modality_overshoot,
                        _ => immune.pharmacologic_overshoot,
                    };
                    damp_field[idx] += gc.lp_at_death * immune.damp_per_lp * overshoot;
                }
            }
        }

        // --- Iron diffusion ---
        grid.diffuse_iron(
            spatial_params.iron_release_per_death,
            spatial_params.neighbor_iron_fraction,
        );

        // --- DAMP diffusion (neighbor spread + clearance) ---
        damp_delta.fill(0.0);
        for r in 0..rows {
            for c in 0..cols {
                let idx = r * cols + c;
                let local = damp_field[idx];
                if local < 0.001 {
                    continue;
                }
                let share = local * immune.damp_diffusion_fraction;
                let (neighbors, count) = grid.neighbors(r, c);
                for &(nr, nc) in &neighbors[..count] {
                    damp_delta[nr * cols + nc] += share;
                }
                damp_delta[idx] -= share * count as f64;
            }
        }
        for i in 0..n_cells {
            damp_field[i] = (damp_field[i] + damp_delta[i]).max(0.0);
            // Clearance decay
            damp_field[i] *= 1.0 - immune.damp_clearance_rate;
        }

        // --- Immune kill (after delay) ---
        if step >= immune_start_step {
            let effective_brake = immune.effective_brake();
            for r in 0..rows {
                for c in 0..cols {
                    let idx = r * cols + c;
                    if grid.cells[idx].state.dead || !grid.cells[idx].is_tumor {
                        continue;
                    }

                    let local_damp = damp_field[idx];
                    if local_damp < 0.01 {
                        continue;
                    }

                    let activation = local_damp / (local_damp + immune.dc_activation_kd);
                    let kill_prob = (activation * immune.immune_kill_rate * (1.0 - effective_brake))
                        .min(0.99);

                    let mut rng = StdRng::seed_from_u64(
                        seed.wrapping_add(900_000_000)
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

    (ferroptosis_kills, immune_kills, damp_field)
}

// ============================================================
// Output types
// ============================================================

#[derive(Serialize)]
struct ConditionResult {
    treatment: String,
    o2_condition: String,
    o2_lambda_um: Option<f64>,
    immune_mode: String,
    total_tumor: usize,
    total_dead: usize,
    ferroptosis_kills: Option<usize>,
    immune_kills: Option<usize>,
    overall_kill_rate: f64,
    normoxic_kill_rate: f64,
    transition_kill_rate: f64,
    hypoxic_kill_rate: f64,
}

/// Compute kill rates for three O2-defined zones:
/// - Normoxic shell: within `shell_depth_um` of tumor edge (O2 > 0.37)
/// - Transition zone: between shell and hypoxic core
/// - Hypoxic core: deeper than `3 * shell_depth_um` (O2 < 0.05)
///
/// Returns (normoxic_rate, transition_rate, hypoxic_rate).
fn zone_kill_rates(grid: &TumorGrid, shell_depth_um: f64) -> (f64, f64, f64) {
    let center_r = grid.rows as f64 / 2.0;
    let center_c = grid.cols as f64 / 2.0;
    let tumor_radius = (grid.rows.min(grid.cols) as f64) * 0.45;
    let cell_size = grid.cell_size_um;

    let deep_threshold_um = shell_depth_um * 3.0;

    let (mut norm_dead, mut norm_total) = (0usize, 0usize);
    let (mut trans_dead, mut trans_total) = (0usize, 0usize);
    let (mut hyp_dead, mut hyp_total) = (0usize, 0usize);

    for r in 0..grid.rows {
        for c in 0..grid.cols {
            let gc = grid.get(r, c);
            if !gc.is_tumor {
                continue;
            }
            let dist_from_center =
                ((r as f64 - center_r).powi(2) + (c as f64 - center_c).powi(2)).sqrt();
            let depth_from_edge_um = (tumor_radius - dist_from_center).max(0.0) * cell_size;

            let (dead_count, total_count) = if depth_from_edge_um < shell_depth_um {
                (&mut norm_dead, &mut norm_total)
            } else if depth_from_edge_um < deep_threshold_um {
                (&mut trans_dead, &mut trans_total)
            } else {
                (&mut hyp_dead, &mut hyp_total)
            };

            *total_count += 1;
            if gc.state.dead {
                *dead_count += 1;
            }
        }
    }

    let rate = |d: usize, t: usize| if t > 0 { d as f64 / t as f64 } else { 0.0 };
    (rate(norm_dead, norm_total), rate(trans_dead, trans_total), rate(hyp_dead, hyp_total))
}

// ============================================================
// Main
// ============================================================

fn main() {
    eprintln!("=== Tumor Microenvironment: Oxygen Gradients ===");
    eprintln!(
        "Grid: {}×{} ({:.1}mm × {:.1}mm)",
        GRID_SIZE, GRID_SIZE,
        GRID_SIZE as f64 * CELL_SIZE_UM / 1000.0,
        GRID_SIZE as f64 * CELL_SIZE_UM / 1000.0,
    );
    eprintln!("O2 penetration sweep: {:?} μm", O2_LAMBDAS);
    eprintln!("Caveats: O2 modulates basal_ros only (not Fenton/SDT); steady-state field\n");

    let params = Params::default();
    let spatial_params = SpatialParams {
        cell_size_um: CELL_SIZE_UM,
        ..Default::default()
    };

    let treatments = [
        (Treatment::Control, "Control"),
        (Treatment::RSL3, "RSL3"),
        (Treatment::SDT, "SDT"),
    ];

    let output_dir = Path::new("output/tme");
    fs::create_dir_all(output_dir).expect("Failed to create output directory");

    let mut all_results: Vec<ConditionResult> = Vec::new();
    let mut all_depth_curves: Vec<(String, Vec<(f64, f64, usize)>)> = Vec::new();

    // --- Baseline: uniform O2 (no gradient) ---
    eprintln!("=== Baseline (uniform O2) ===\n");
    for (tx, tx_name) in &treatments {
        let mut grid = TumorGrid::generate(GRID_SIZE, GRID_SIZE, CELL_SIZE_UM, SEED);
        run_spatial(
            &mut grid, *tx, &params, &spatial_params,
            SEED.wrapping_add((*tx as u64) * 10_000_000),
        );

        let census = grid.census();
        let overall = census.total_dead as f64 / census.total_tumor.max(1) as f64;
        let (norm_r, trans_r, hyp_r) = zone_kill_rates(&grid, ZONE_REF_LAMBDA);
        eprintln!("  {}: overall={:.1}%, normoxic={:.1}%, transition={:.1}%, hypoxic={:.1}%",
            tx_name, overall * 100.0, norm_r * 100.0, trans_r * 100.0, hyp_r * 100.0);

        all_results.push(ConditionResult {
            treatment: tx_name.to_string(),
            o2_condition: "uniform".to_string(),
            o2_lambda_um: None,
            immune_mode: "off".to_string(),
            total_tumor: census.total_tumor,
            total_dead: census.total_dead,
            ferroptosis_kills: None,
            immune_kills: None,
            overall_kill_rate: overall,
            normoxic_kill_rate: norm_r,
            transition_kill_rate: trans_r,
            hypoxic_kill_rate: hyp_r,
        });

        let label = format!("{}_uniform", tx_name);
        all_depth_curves.push((label, depth_kill_curve(&grid)));

        let heatmap = death_heatmap(&grid);
        let path = output_dir.join(format!("death_{}_uniform.csv", tx_name.to_lowercase()));
        write_heatmap_csv(&path, &heatmap).expect("Failed to write heatmap");
    }

    // --- O2 gradient conditions ---
    for &lambda in O2_LAMBDAS {
        eprintln!("\n=== O2 gradient (λ = {} μm) ===\n", lambda);

        for (tx, tx_name) in &treatments {
            let mut grid = TumorGrid::generate(GRID_SIZE, GRID_SIZE, CELL_SIZE_UM, SEED);

            // Apply O2 gradient BEFORE simulation (modifies cell.basal_ros)
            let o2_map = apply_o2_gradient(&mut grid, lambda);

            run_spatial(
                &mut grid, *tx, &params, &spatial_params,
                SEED.wrapping_add((*tx as u64) * 10_000_000),
            );

            let census = grid.census();
            let overall = census.total_dead as f64 / census.total_tumor.max(1) as f64;
            let (norm_r, trans_r, hyp_r) = zone_kill_rates(&grid, ZONE_REF_LAMBDA);
            eprintln!("  {}: overall={:.1}%, normoxic={:.1}%, transition={:.1}%, hypoxic={:.1}%",
                tx_name, overall * 100.0, norm_r * 100.0, trans_r * 100.0, hyp_r * 100.0);

            all_results.push(ConditionResult {
                treatment: tx_name.to_string(),
                o2_condition: format!("gradient_{}um", lambda as u64),
                o2_lambda_um: Some(lambda),
                immune_mode: "off".to_string(),
                total_tumor: census.total_tumor,
                total_dead: census.total_dead,
                ferroptosis_kills: None,
                immune_kills: None,
                overall_kill_rate: overall,
                normoxic_kill_rate: norm_r,
                transition_kill_rate: trans_r,
                hypoxic_kill_rate: hyp_r,
            });

            let label = format!("{}_{}", tx_name, lambda as u64);
            all_depth_curves.push((label, depth_kill_curve(&grid)));

            // Save heatmaps for the default lambda (120μm) only to avoid file bloat
            if (lambda - 120.0).abs() < 1.0 {
                let death_hm = death_heatmap(&grid);
                let path = output_dir.join(format!("death_{}_o2gradient.csv", tx_name.to_lowercase()));
                write_heatmap_csv(&path, &death_hm).expect("Failed to write heatmap");

                // O2 heatmap
                let mut o2_hm = ndarray::Array2::<u8>::zeros((grid.rows, grid.cols));
                for &(r, c, o2) in &o2_map {
                    o2_hm[[r, c]] = (o2 * 255.0).round() as u8;
                }
                let path = output_dir.join("o2_field.csv");
                write_heatmap_csv(&path, &o2_hm).expect("Failed to write O2 heatmap");
            }
        }
    }

    // --- Immune coupling (Feature B) at λ=120μm ---
    let immune_modes: Vec<(&str, ImmuneConfig)> = vec![
        ("immune_on", ImmuneConfig::default_no_pd1()),
        ("immune_anti_pd1", ImmuneConfig::default_no_pd1().with_anti_pd1()),
    ];

    eprintln!("\n=== Spatial Immune Coupling (O2 gradient λ=120μm) ===");
    eprintln!("NOTE: LP at death ≈ 10.0 for all treatments (threshold-locked).");
    eprintln!("DAMP differential comes from kill DENSITY, not per-cell DAMP quality.");
    eprintln!("Immune model: resident T cell phase only (0-48h), not systemic.\n");

    for (immune_label, immune_cfg) in &immune_modes {
        eprintln!("--- Immune mode: {} (brake={:.0}%) ---\n",
            immune_label, immune_cfg.effective_brake() * 100.0);

        for (tx, tx_name) in &treatments {
            let mut grid = TumorGrid::generate(GRID_SIZE, GRID_SIZE, CELL_SIZE_UM, SEED);
            apply_o2_gradient(&mut grid, 120.0);

            let (ferr_kills, imm_kills, final_damp) = run_spatial_with_immune(
                &mut grid, *tx, &params, &spatial_params, immune_cfg,
                SEED.wrapping_add((*tx as u64) * 10_000_000),
            );

            let census = grid.census();
            let overall = census.total_dead as f64 / census.total_tumor.max(1) as f64;
            let (norm_r, trans_r, hyp_r) = zone_kill_rates(&grid, ZONE_REF_LAMBDA);
            eprintln!("  {}: overall={:.1}% (ferr={}, immune={}), hypoxic={:.1}%",
                tx_name, overall * 100.0, ferr_kills, imm_kills, hyp_r * 100.0);

            // Export DAMP and immune-kill heatmaps for the first immune mode only
            if *immune_label == "immune_on" {
                // DAMP concentration heatmap (scaled to u8 for CSV export)
                let damp_max = final_damp.iter().cloned().fold(0.0_f64, f64::max).max(1.0);
                let mut damp_hm = ndarray::Array2::<u8>::zeros((grid.rows, grid.cols));
                for r in 0..grid.rows {
                    for c in 0..grid.cols {
                        let val = final_damp[r * grid.cols + c] / damp_max;
                        damp_hm[[r, c]] = (val * 255.0).round() as u8;
                    }
                }
                let path = output_dir.join(format!("damp_field_{}.csv", tx_name.to_lowercase()));
                write_heatmap_csv(&path, &damp_hm).expect("Failed to write DAMP heatmap");

                // Final death map for the immune-enabled run (same encoding as other
                // death heatmaps: 0=stromal, 1=dead tumor, 2=alive tumor; does NOT
                // distinguish immune kills from ferroptotic kills)
                let death_hm = death_heatmap(&grid);
                let path = output_dir.join(format!("death_{}_immune_run.csv", tx_name.to_lowercase()));
                write_heatmap_csv(&path, &death_hm).expect("Failed to write death heatmap");
            }

            all_results.push(ConditionResult {
                treatment: tx_name.to_string(),
                o2_condition: "gradient_120um".to_string(),
                o2_lambda_um: Some(120.0),
                immune_mode: immune_label.to_string(),
                total_tumor: census.total_tumor,
                total_dead: census.total_dead,
                ferroptosis_kills: Some(ferr_kills),
                immune_kills: Some(imm_kills),
                overall_kill_rate: overall,
                normoxic_kill_rate: norm_r,
                transition_kill_rate: trans_r,
                hypoxic_kill_rate: hyp_r,
            });

            let label = format!("{}_120_{}", tx_name, immune_label);
            all_depth_curves.push((label, depth_kill_curve(&grid)));
        }
        eprintln!();
    }

    // --- Write aggregated outputs ---
    let curves_path = output_dir.join("depth_kill_curves.csv");
    write_depth_curves_csv(&curves_path, &all_depth_curves).expect("Failed to write depth curves");

    let summary_path = output_dir.join("tme_summary.json");
    write_json(&summary_path, &all_results).expect("Failed to write summary");

    // --- Print comparison table ---
    eprintln!("\n=== Comparison Table ===\n");
    eprintln!(
        "{:<10} {:<20} {:<18} {:>10} {:>10} {:>10} {:>10}",
        "Treatment", "O2 Condition", "Immune", "Overall", "Normoxic", "Transit.", "Hypoxic"
    );
    eprintln!("{}", "-".repeat(100));
    for r in &all_results {
        eprintln!(
            "{:<10} {:<20} {:<18} {:>9.1}% {:>9.1}% {:>9.1}% {:>9.1}%",
            r.treatment, r.o2_condition, r.immune_mode,
            r.overall_kill_rate * 100.0,
            r.normoxic_kill_rate * 100.0,
            r.transition_kill_rate * 100.0,
            r.hypoxic_kill_rate * 100.0,
        );
    }

    // --- Sensitivity check: is the SDT-vs-RSL3 differential robust? ---
    eprintln!("\n=== Sensitivity: SDT advantage over RSL3 in hypoxic zone ===\n");
    eprintln!("{:<10} {:>15} {:>15} {:>15}", "λ (μm)", "RSL3 hypoxic", "SDT hypoxic", "SDT/RSL3 ratio");
    eprintln!("{}", "-".repeat(57));
    for &lambda in O2_LAMBDAS {
        let cond = format!("gradient_{}um", lambda as u64);
        let rsl3_hyp = all_results.iter()
            .find(|r| r.treatment == "RSL3" && r.o2_condition == cond)
            .map(|r| r.hypoxic_kill_rate).unwrap_or(0.0);
        let sdt_hyp = all_results.iter()
            .find(|r| r.treatment == "SDT" && r.o2_condition == cond)
            .map(|r| r.hypoxic_kill_rate).unwrap_or(0.0);
        let ratio = if rsl3_hyp > 0.001 { sdt_hyp / rsl3_hyp } else { f64::INFINITY };
        eprintln!("{:<10} {:>14.1}% {:>14.1}% {:>15.1}×", lambda, rsl3_hyp * 100.0, sdt_hyp * 100.0, ratio);
    }

    eprintln!("\nZone definitions (fixed at λ_ref = {ZONE_REF_LAMBDA}μm): normoxic = within {ZONE_REF_LAMBDA}μm of edge, transition = {ZONE_REF_LAMBDA}-{}μm, hypoxic = deeper than {}μm.", ZONE_REF_LAMBDA * 3.0, ZONE_REF_LAMBDA * 3.0);
    eprintln!("If SDT/RSL3 ratio is consistently >1 across all λ values,");
    eprintln!("the finding is robust: SDT maintains efficacy in hypoxia better than RSL3.");

    eprintln!("\nOutputs saved to {}/", output_dir.display());
}
