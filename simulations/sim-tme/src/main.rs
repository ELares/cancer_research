//! Tumor Microenvironment Simulation: Oxygen Gradients and Hypoxia
//!
//! Extends the spatial tumor model with oxygen gradients to quantify how
//! hypoxia protects tumor core cells from ferroptosis. Compares treatment
//! efficacy under uniform vs gradient O2 conditions.
//!
//! Key prediction: SDT (exogenous ROS, O2-independent) should maintain
//! efficacy in hypoxic cores better than RSL3 (depends on basal ROS,
//! O2-dependent). This differential should be robust to the O2 penetration
//! length parameter.
//!
//! Biology: mitochondrial ETC generates basal ROS using O2 as terminal
//! electron acceptor (Biology2e Ch.7-8, Murphy Biochem J 2009). In hypoxia,
//! ETC activity drops → less basal ROS → less lipid peroxidation cascade.
//! O2 penetration ~100-150μm from vasculature (Vaupel, Cancer Res 1989).
//!
//! Caveats:
//! - O2 modulates basal_ros only, not Fenton directly (conservative — Fenton
//!   substrate H2O2 also requires O2 via superoxide→SOD pathway)
//! - SDT/PDT modeled as O2-independent (Type I mechanism; Type II requires
//!   O2 for singlet oxygen — conservative, may overstate physical modality
//!   efficacy in hypoxia)
//! - Steady-state O2 field (no consumption dynamics)
//! - Distance from tumor edge as O2 source (no explicit vasculature)
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
// Output types
// ============================================================

#[derive(Serialize)]
struct ConditionResult {
    treatment: String,
    o2_condition: String,
    o2_lambda_um: Option<f64>,
    total_tumor: usize,
    total_dead: usize,
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
        let (norm_r, trans_r, hyp_r) = zone_kill_rates(&grid, 120.0);
        eprintln!("  {}: overall={:.1}%, normoxic={:.1}%, transition={:.1}%, hypoxic={:.1}%",
            tx_name, overall * 100.0, norm_r * 100.0, trans_r * 100.0, hyp_r * 100.0);

        all_results.push(ConditionResult {
            treatment: tx_name.to_string(),
            o2_condition: "uniform".to_string(),
            o2_lambda_um: None,
            total_tumor: census.total_tumor,
            total_dead: census.total_dead,
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
            let (norm_r, trans_r, hyp_r) = zone_kill_rates(&grid, lambda);
            eprintln!("  {}: overall={:.1}%, normoxic={:.1}%, transition={:.1}%, hypoxic={:.1}%",
                tx_name, overall * 100.0, norm_r * 100.0, trans_r * 100.0, hyp_r * 100.0);

            all_results.push(ConditionResult {
                treatment: tx_name.to_string(),
                o2_condition: format!("gradient_{}um", lambda as u64),
                o2_lambda_um: Some(lambda),
                total_tumor: census.total_tumor,
                total_dead: census.total_dead,
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

    // --- Write aggregated outputs ---
    let curves_path = output_dir.join("depth_kill_curves.csv");
    write_depth_curves_csv(&curves_path, &all_depth_curves).expect("Failed to write depth curves");

    let summary_path = output_dir.join("tme_summary.json");
    write_json(&summary_path, &all_results).expect("Failed to write summary");

    // --- Print comparison table ---
    eprintln!("\n=== Comparison Table ===\n");
    eprintln!(
        "{:<10} {:<20} {:>10} {:>10} {:>10} {:>10}",
        "Treatment", "O2 Condition", "Overall", "Normoxic", "Transit.", "Hypoxic"
    );
    eprintln!("{}", "-".repeat(82));
    for r in &all_results {
        eprintln!(
            "{:<10} {:<20} {:>9.1}% {:>9.1}% {:>9.1}% {:>9.1}%",
            r.treatment, r.o2_condition,
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

    eprintln!("\nZone definitions: normoxic = within λ of edge, transition = λ to 3λ, hypoxic = deeper than 3λ.");
    eprintln!("If SDT/RSL3 ratio is consistently >1 across all λ values,");
    eprintln!("the finding is robust: SDT maintains efficacy in hypoxia better than RSL3.");

    eprintln!("\nOutputs saved to {}/", output_dir.display());
}
