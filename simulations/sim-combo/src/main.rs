//! Combination Therapy Optimizer
//!
//! Combines temporal (vulnerability window) + immune (ICD cascade) models.
//! Sweeps SDT timing × anti-PD1 timing to find optimal treatment schedule.
//!
//! 3-phase model:
//! 1. Chemotherapy kills proliferating cells, persisters survive
//! 2. SDT at variable delay (0-28 days) targets persisters
//! 3. Immune cascade from ICD ± anti-PD1

use std::path::PathBuf;

use clap::Parser;
use rand::prelude::*;
use rayon::prelude::*;

use ferroptosis_core::biochem::sim_cell;
use ferroptosis_core::cell::{gen_recovered_persister, RecoveryRates, Treatment};
use ferroptosis_core::immune::immune_cascade;
use ferroptosis_core::io::write_json;
use ferroptosis_core::params::{ImmuneParams, Params};

#[derive(Parser)]
#[command(name = "sim-combo", about = "Combination therapy optimizer")]
struct Args {
    /// Cells per condition.
    #[arg(long, default_value_t = 50_000)]
    n_cells: usize,

    /// Random seed.
    #[arg(long, default_value_t = 42)]
    seed: u64,

    /// Output directory.
    #[arg(long, default_value = "output/combo")]
    output_dir: PathBuf,
}

fn main() {
    let args = Args::parse();

    eprintln!("=== Combination Therapy Optimizer ===");
    eprintln!("Cells per condition: {}", args.n_cells);
    eprintln!("Seed: {}\n", args.seed);

    let params = Params::default();
    let recovery = RecoveryRates::default();
    let immune_params = ImmuneParams::default();

    // SDT delay: days after chemo withdrawal
    let sdt_delays_days: Vec<f64> = vec![0.0, 1.0, 3.0, 7.0, 14.0, 21.0, 28.0];

    // Anti-PD1 options
    let anti_pd1_options = [false, true];

    // Second-line treatment options
    let treatments = [
        (Treatment::SDT, "SDT"),
        (Treatment::PDT, "PDT"),
        (Treatment::RSL3, "RSL3"),
    ];

    std::fs::create_dir_all(&args.output_dir).expect("Failed to create output dir");

    let mut all_results = Vec::new();

    // Starting population: 1000 tumor cells post-chemo
    // Assume chemo killed 90% of proliferating cells, all survivors are persisters
    let initial_tumor_cells = 1000_usize;

    for &delay_days in &sdt_delays_days {
        for (tx, tx_name) in &treatments {
            for &with_pd1 in &anti_pd1_options {
                let n = args.n_cells;

                // Phase 2: simulate ferroptosis on recovered persisters
                let outcomes: Vec<(bool, f64, f64, f64)> = (0..n)
                    .into_par_iter()
                    .map(|i| {
                        let mut cell_rng = StdRng::seed_from_u64(
                            args.seed.wrapping_add(i as u64 * 2)
                                .wrapping_add(delay_days as u64 * 1_000_000),
                        );
                        let mut sim_rng = StdRng::seed_from_u64(
                            args.seed.wrapping_add(i as u64 * 2 + 1)
                                .wrapping_add(delay_days as u64 * 1_000_000),
                        );
                        let cell = gen_recovered_persister(delay_days, &recovery, &mut cell_rng);
                        sim_cell(&cell, *tx, &params, &mut sim_rng)
                    })
                    .collect();

                let dead_cell_lps: Vec<f64> = outcomes
                    .iter()
                    .filter(|(dead, _, _, _)| *dead)
                    .map(|(_, lp, _, _)| *lp)
                    .collect();

                let ferroptosis_kill_rate = dead_cell_lps.len() as f64 / n as f64;

                // Phase 3: immune cascade
                let immune = immune_cascade(
                    &dead_cell_lps,
                    initial_tumor_cells,
                    &immune_params,
                    with_pd1,
                );

                // Total tumor reduction: ferroptosis kills + immune kills
                let ferroptosis_killed = (ferroptosis_kill_rate * initial_tumor_cells as f64).round() as usize;
                let immune_killed = immune.immune_kills.round() as usize;
                let total_killed = (ferroptosis_killed + immune_killed).min(initial_tumor_cells);
                let survivors = initial_tumor_cells - total_killed;
                let survival_fraction = survivors as f64 / initial_tumor_cells as f64;

                eprintln!(
                    "  Day {:2.0} + {:<4} {}: ferro={:.1}%, immune={}, survivors={}/{} ({:.1}%)",
                    delay_days,
                    tx_name,
                    if with_pd1 { "+PD1" } else { "    " },
                    ferroptosis_kill_rate * 100.0,
                    immune_killed,
                    survivors,
                    initial_tumor_cells,
                    survival_fraction * 100.0,
                );

                all_results.push(serde_json::json!({
                    "sdt_delay_days": delay_days,
                    "treatment": tx_name,
                    "with_anti_pd1": with_pd1,
                    "initial_tumor_cells": initial_tumor_cells,
                    "ferroptosis_kill_rate": ferroptosis_kill_rate,
                    "ferroptosis_killed": ferroptosis_killed,
                    "immune_kills": immune.immune_kills,
                    "total_killed": total_killed,
                    "survivors": survivors,
                    "survival_fraction": survival_fraction,
                    "total_damps": immune.total_damps,
                    "damp_per_dead_cell": immune.damp_per_dead_cell,
                    "primed_tcells": immune.primed_tcells,
                }));
            }
        }
        eprintln!();
    }

    // Find optimal schedule
    let best = all_results
        .iter()
        .min_by(|a, b| {
            a["survival_fraction"]
                .as_f64()
                .unwrap()
                .partial_cmp(&b["survival_fraction"].as_f64().unwrap())
                .unwrap()
        });

    if let Some(best) = best {
        eprintln!("=== Optimal Schedule ===");
        eprintln!(
            "  {} at day {:.0} {}: {:.1}% survival ({} survivors / {})",
            best["treatment"],
            best["sdt_delay_days"].as_f64().unwrap(),
            if best["with_anti_pd1"].as_bool().unwrap() { "+ anti-PD1" } else { "" },
            best["survival_fraction"].as_f64().unwrap() * 100.0,
            best["survivors"],
            best["initial_tumor_cells"],
        );
    }

    // Save results
    let json_path = args.output_dir.join("combo_results.json");
    write_json(&json_path, &all_results).expect("Failed to write JSON");

    eprintln!("\n=== Output saved to {} ===", args.output_dir.display());
}
