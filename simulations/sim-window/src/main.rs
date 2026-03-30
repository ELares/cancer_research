//! Ferroptosis Vulnerability Window Simulation
//!
//! Models how persister cell ferroptosis sensitivity changes over time
//! after chemotherapy withdrawal as defense pathways recover.
//!
//! Answers: "How long is the therapeutic window for SDT intervention?"

use std::path::PathBuf;

use clap::Parser;
use rand::prelude::*;
use rayon::prelude::*;

use ferroptosis_core::biochem::sim_cell;
use ferroptosis_core::cell::{gen_recovered_persister, RecoveryRates, Treatment};
use ferroptosis_core::io::{write_json, write_window_csv};
use ferroptosis_core::params::Params;
use ferroptosis_core::stats::wilson_ci;

#[derive(Parser)]
#[command(name = "sim-window", about = "Ferroptosis vulnerability window dynamics")]
struct Args {
    /// Cells per condition.
    #[arg(long, default_value_t = 100_000)]
    n_cells: usize,

    /// Random seed.
    #[arg(long, default_value_t = 42)]
    seed: u64,

    /// Output directory.
    #[arg(long, default_value = "output/window")]
    output_dir: PathBuf,
}

fn main() {
    let args = Args::parse();

    eprintln!("=== Vulnerability Window Simulation ===");
    eprintln!("Cells per condition: {}", args.n_cells);
    eprintln!("Seed: {}\n", args.seed);

    let params = Params::default();
    let recovery = RecoveryRates::default();

    // Timepoints: 0 to 28 days
    let timepoints_hours: Vec<f64> = vec![
        0.0, 6.0, 12.0, 24.0, 48.0, 72.0,
        168.0,   // 1 week
        336.0,   // 2 weeks
        672.0,   // 4 weeks
    ];

    let treatments = [
        (Treatment::Control, "Control"),
        (Treatment::RSL3, "RSL3"),
        (Treatment::SDT, "SDT"),
        (Treatment::PDT, "PDT"),
    ];

    std::fs::create_dir_all(&args.output_dir).expect("Failed to create output dir");

    let mut csv_rows = Vec::new();
    let mut json_results = Vec::new();

    for &hours in &timepoints_hours {
        let days = hours / 24.0;
        eprintln!("--- Timepoint: {:.0}h ({:.1} days) ---", hours, days);

        for (tx, tx_name) in &treatments {
            let n = args.n_cells;
            let outcomes: Vec<(bool, f64, f64, f64)> = (0..n)
                .into_par_iter()
                .map(|i| {
                    let mut cell_rng = StdRng::seed_from_u64(
                        args.seed.wrapping_add(i as u64 * 2).wrapping_add(hours as u64 * 1_000_000),
                    );
                    let mut sim_rng = StdRng::seed_from_u64(
                        args.seed.wrapping_add(i as u64 * 2 + 1).wrapping_add(hours as u64 * 1_000_000),
                    );
                    let cell = gen_recovered_persister(days, &recovery, &mut cell_rng);
                    sim_cell(&cell, *tx, &params, &mut sim_rng)
                })
                .collect();

            let dead = outcomes.iter().filter(|(d, _, _, _)| *d).count();
            let rate = dead as f64 / n as f64;
            let (ci_lo, ci_hi) = wilson_ci(n, dead);

            eprintln!(
                "  {:<8} → Death: {:7.3}% [{:.3}-{:.3}]",
                tx_name,
                rate * 100.0,
                ci_lo * 100.0,
                ci_hi * 100.0,
            );

            csv_rows.push((hours, tx_name.to_string(), rate, ci_lo, ci_hi));
            json_results.push(serde_json::json!({
                "timepoint_hours": hours,
                "timepoint_days": days,
                "treatment": tx_name,
                "n_cells": n,
                "n_dead": dead,
                "death_rate": rate,
                "ci_low": ci_lo,
                "ci_high": ci_hi,
                "mean_lp": outcomes.iter().map(|(_, l, _, _)| l).sum::<f64>() / n as f64,
                "mean_gsh": outcomes.iter().map(|(_, _, g, _)| g).sum::<f64>() / n as f64,
                "mean_gpx4": outcomes.iter().map(|(_, _, _, p)| p).sum::<f64>() / n as f64,
            }));
        }
        eprintln!();
    }

    // Save CSV
    let csv_path = args.output_dir.join("vulnerability_window.csv");
    write_window_csv(&csv_path, &csv_rows).expect("Failed to write CSV");

    // Save JSON
    let json_path = args.output_dir.join("vulnerability_window.json");
    write_json(&json_path, &json_results).expect("Failed to write JSON");

    // Recovery rate sensitivity analysis
    eprintln!("=== Sensitivity: Recovery Rate ±50% ===");
    for rate_name in &["fsp1", "gpx4", "nrf2", "gsh"] {
        for mult in [0.5, 1.5] {
            let mut r = RecoveryRates::default();
            match *rate_name {
                "fsp1" => r.fsp1_half_recovery_days *= mult,
                "gpx4" => r.gpx4_half_recovery_days *= mult,
                "nrf2" => r.nrf2_half_recovery_days *= mult,
                "gsh" => r.gsh_half_recovery_days *= mult,
                _ => {}
            }
            // Test at 7 days (1 week) with SDT
            let days = 7.0;
            let n = 50_000;
            let outcomes: Vec<(bool, f64, f64, f64)> = (0..n)
                .into_par_iter()
                .map(|i| {
                    let mut cell_rng = StdRng::seed_from_u64(args.seed.wrapping_add(i as u64 * 2 + 99999));
                    let mut sim_rng = StdRng::seed_from_u64(args.seed.wrapping_add(i as u64 * 2 + 100000));
                    let cell = gen_recovered_persister(days, &r, &mut cell_rng);
                    sim_cell(&cell, Treatment::SDT, &params, &mut sim_rng)
                })
                .collect();
            let dead = outcomes.iter().filter(|(d, _, _, _)| *d).count();
            let rate = dead as f64 / n as f64;
            eprintln!(
                "  {} t½ ×{:.1}: SDT death at 7d = {:.2}%",
                rate_name, mult, rate * 100.0
            );
        }
    }

    eprintln!("\n=== Output saved to {} ===", args.output_dir.display());
}
