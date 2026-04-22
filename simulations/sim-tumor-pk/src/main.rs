//! Tumor PBPK compartment simulation.
//!
//! Models time-varying drug delivery from plasma to tumor interstitium,
//! then feeds the concentration schedule into the ferroptosis engine.
//! Demonstrates that tumor-specific PK barriers create a large
//! protection factor: 2D culture (constant drug) vs in-vivo (time-varying,
//! reduced by blood flow, vascular permeability, and IFP).
//!
//! Run: cargo run --release -p sim-tumor-pk

use std::fs;
use std::path::Path;

use ferroptosis_core::cell::{gen_cell, Phenotype};
use ferroptosis_core::params::Params;
use ferroptosis_core::stats::wilson_ci;
use ferroptosis_core::tumor_pk::*;
use rand::SeedableRng;
use rayon::prelude::*;
use serde::Serialize;

const N_CELLS: usize = 10_000;
const SEED: u64 = 42;
const N_STEPS: usize = 180;

#[derive(Serialize)]
struct ScenarioResult {
    tumor_type: String,
    context: String, // "tumor_pk" or "2d_reference"
    n_cells: usize,
    n_dead: usize,
    death_rate: f64,
    ci_low: f64,
    ci_high: f64,
    mean_lp: f64,
    mean_gsh: f64,
    mean_gpx4: f64,
    peak_c_interstitial: f64,
    auc_c_interstitial: f64,
    protection_factor: Option<f64>,
}

fn run_scenario(
    conc_schedule: &[f64],
    params: &Params,
    seed: u64,
) -> (usize, f64, f64, f64) {
    let results: Vec<PKCellResult> = (0..N_CELLS)
        .into_par_iter()
        .map(|i| {
            let cell_seed = seed.wrapping_add(i as u64).wrapping_add(1_000_000);
            let mut rng = rand::rngs::StdRng::seed_from_u64(cell_seed);
            let cell = gen_cell(Phenotype::Persister, &mut rng);
            sim_cell_with_pk(&cell, params, conc_schedule, params.rsl3_gpx4_inhib, cell_seed)
        })
        .collect();

    let n_dead = results.iter().filter(|r| r.dead).count();
    let sum_lp: f64 = results.iter().map(|r| r.final_lp).sum();
    let sum_gsh: f64 = results.iter().map(|r| r.final_gsh).sum();
    let sum_gpx4: f64 = results.iter().map(|r| r.final_gpx4).sum();
    let n = results.len() as f64;

    (n_dead, sum_lp / n, sum_gsh / n, sum_gpx4 / n)
}

fn main() {
    let output_dir = Path::new("output/tumor-pk");
    fs::create_dir_all(output_dir).expect("Failed to create output directory");

    let params = Params::default();

    // Define tumor scenarios
    let tumors: Vec<(&str, TumorPKParams)> = vec![
        ("Breast", breast_tumor()),
        ("Pancreatic", pancreatic_tumor()),
        ("GBM", glioblastoma_tumor()),
    ];

    let plasma = rsl3_iv_bolus();

    eprintln!("=== Tumor PBPK Compartment Model ===");
    eprintln!("Drug: RSL3-like IV bolus (t_half=30 min)");
    eprintln!("Phenotype: Persister (FSP1-low)");
    eprintln!("Cells per condition: {N_CELLS}");
    eprintln!("All tumor PK parameters ESTIMATED (no textbook coverage).\n");

    let mut all_results: Vec<ScenarioResult> = Vec::new();
    // --- 2D reference (constant drug concentration = 1.0 at all steps) ---
    // This represents ideal 2D culture where drug bathes cells directly.
    // Do NOT run through the tumor ODE — that models tumor barriers which
    // don't exist in 2D culture.
    let ref_death_rate;
    {
        let conc_schedule: Vec<f64> = vec![1.0; N_STEPS];
        let (n_dead, mean_lp, mean_gsh, mean_gpx4) =
            run_scenario(&conc_schedule, &params, SEED);
        let (ci_lo, ci_hi) = wilson_ci(N_CELLS, n_dead);
        ref_death_rate = n_dead as f64 / N_CELLS as f64;

        eprintln!(
            "  2D reference: death_rate={:.1}% [{:.1}-{:.1}], LP={:.2}, GSH={:.2}, GPX4={:.3}",
            ref_death_rate * 100.0,
            ci_lo * 100.0,
            ci_hi * 100.0,
            mean_lp,
            mean_gsh,
            mean_gpx4
        );

        all_results.push(ScenarioResult {
            tumor_type: "2D reference".to_string(),
            context: "2d_reference".to_string(),
            n_cells: N_CELLS,
            n_dead,
            death_rate: ref_death_rate,
            ci_low: ci_lo,
            ci_high: ci_hi,
            mean_lp,
            mean_gsh,
            mean_gpx4,
            peak_c_interstitial: 1.0,
            auc_c_interstitial: N_STEPS as f64,
            protection_factor: None,
        });
    }
    eprintln!();

    // --- Tumor-specific PK ---
    let mut timecourse_rows: Vec<String> = vec![
        "time_min,tumor_type,c_plasma,c_vascular,c_interstitial".to_string()
    ];

    for (tumor_name, tumor_params) in &tumors {
        let pk_result = solve_tumor_pk(&plasma, tumor_params, N_STEPS, 100);

        // Record timecourse
        for i in 0..N_STEPS {
            timecourse_rows.push(format!(
                "{},{},{:.6},{:.6},{:.6}",
                pk_result.time_min[i],
                tumor_name,
                pk_result.c_plasma[i],
                pk_result.c_vascular[i],
                pk_result.c_interstitial[i]
            ));
        }

        let peak_ci: f64 = pk_result.c_interstitial.iter().cloned().fold(0.0, f64::max);
        let auc_ci: f64 = pk_result.c_interstitial.iter().sum();

        let (n_dead, mean_lp, mean_gsh, mean_gpx4) =
            run_scenario(&pk_result.c_interstitial, &params, SEED);
        let (ci_lo, ci_hi) = wilson_ci(N_CELLS, n_dead);
        let death_rate = n_dead as f64 / N_CELLS as f64;
        let protection = if death_rate > 0.001 {
            ref_death_rate / death_rate
        } else {
            f64::INFINITY
        };

        eprintln!(
            "  {}: death_rate={:.1}% [{:.1}-{:.1}], peak_Ci={:.3}, AUC={:.1}, protection={:.1}×, LP={:.2}",
            tumor_name,
            death_rate * 100.0,
            ci_lo * 100.0,
            ci_hi * 100.0,
            peak_ci,
            auc_ci,
            protection,
            mean_lp
        );

        all_results.push(ScenarioResult {
            tumor_type: tumor_name.to_string(),
            context: "tumor_pk".to_string(),
            n_cells: N_CELLS,
            n_dead,
            death_rate,
            ci_low: ci_lo,
            ci_high: ci_hi,
            mean_lp,
            mean_gsh,
            mean_gpx4,
            peak_c_interstitial: peak_ci,
            auc_c_interstitial: auc_ci,
            protection_factor: Some(protection),
        });
    }

    // --- Write outputs ---
    let tc_path = output_dir.join("tumor_pk_timecourse.csv");
    fs::write(&tc_path, timecourse_rows.join("\n")).expect("Failed to write timecourse");

    let summary_path = output_dir.join("tumor_pk_summary.json");
    let json = serde_json::to_string_pretty(&all_results).expect("JSON serialization failed");
    fs::write(&summary_path, json).expect("Failed to write summary");

    // --- Comparison table ---
    eprintln!("\n=== Protection Factor Summary ===\n");
    eprintln!(
        "{:<25} {:>12} {:>12} {:>12} {:>12}",
        "Tumor Type", "Death Rate", "Peak C_i", "AUC C_i", "Protection"
    );
    eprintln!("{}", "-".repeat(75));
    for r in &all_results {
        let prot = r
            .protection_factor
            .map(|p| format!("{:.1}×", p))
            .unwrap_or_else(|| "reference".to_string());
        eprintln!(
            "{:<25} {:>11.1}% {:>12.3} {:>12.1} {:>12}",
            r.tumor_type,
            r.death_rate * 100.0,
            r.peak_c_interstitial,
            r.auc_c_interstitial,
            prot
        );
    }

    eprintln!("\nOutputs saved to {}/", output_dir.display());
}
