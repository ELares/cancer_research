//! ICD-Immune Cascade Comparison
//!
//! Validates the paper's claim that physical modalities produce more
//! immunogenic cell death than pharmacologic ferroptosis inducers.
//!
//! Key insight: SDT/PDT kill via runaway LP cascade (LP >> threshold at death),
//! so they release MORE DAMPs per dead cell than RSL3 (LP ≈ threshold at death).

use std::path::PathBuf;

use clap::Parser;
use rand::prelude::*;
use rayon::prelude::*;

use ferroptosis_core::biochem::sim_cell;
use ferroptosis_core::cell::{gen_cell, Phenotype, Treatment};
use ferroptosis_core::immune::immune_cascade;
use ferroptosis_core::io::write_json;
use ferroptosis_core::params::{ImmuneParams, Params};

#[derive(Parser)]
#[command(name = "sim-icd", about = "ICD-immune cascade comparison across treatments")]
struct Args {
    /// Cells per condition.
    #[arg(long, default_value_t = 100_000)]
    n_cells: usize,

    /// Random seed.
    #[arg(long, default_value_t = 42)]
    seed: u64,

    /// Output directory.
    #[arg(long, default_value = "output/icd")]
    output_dir: PathBuf,
}

fn main() {
    let args = Args::parse();

    eprintln!("=== ICD-Immune Cascade Simulation ===");
    eprintln!("Cells per condition: {}", args.n_cells);
    eprintln!("Seed: {}\n", args.seed);

    let params = Params::default();
    let immune_params = ImmuneParams::default();

    let phenotypes = [
        (Phenotype::Persister, "Persister"),
        (Phenotype::OXPHOS, "OXPHOS"),
        (Phenotype::Glycolytic, "Glycolytic"),
    ];
    let treatments = [
        (Treatment::Control, "Control"),
        (Treatment::RSL3, "RSL3"),
        (Treatment::SDT, "SDT"),
        (Treatment::PDT, "PDT"),
    ];

    std::fs::create_dir_all(&args.output_dir).expect("Failed to create output dir");

    let mut all_results = Vec::new();

    for (pheno, pheno_name) in &phenotypes {
        eprintln!("--- Phenotype: {} ---", pheno_name);

        for (tx, tx_name) in &treatments {
            let n = args.n_cells;

            // Run biochemistry simulation, collecting LP at death for dead cells
            let outcomes: Vec<(bool, f64, f64, f64)> = (0..n)
                .into_par_iter()
                .map(|i| {
                    let mut cell_rng = StdRng::seed_from_u64(args.seed.wrapping_add(i as u64 * 2));
                    let mut sim_rng = StdRng::seed_from_u64(args.seed.wrapping_add(i as u64 * 2 + 1));
                    let cell = gen_cell(*pheno, &mut cell_rng);
                    sim_cell(&cell, *tx, &params, &mut sim_rng)
                })
                .collect();

            let dead_cell_lps: Vec<f64> = outcomes
                .iter()
                .filter(|(dead, _, _, _)| *dead)
                .map(|(_, lp, _, _)| *lp)
                .collect();

            let n_dead = dead_cell_lps.len();
            let death_rate = n_dead as f64 / n as f64;

            let avg_lp_at_death = if n_dead > 0 {
                dead_cell_lps.iter().sum::<f64>() / n_dead as f64
            } else {
                0.0
            };

            // Run immune cascade (with and without anti-PD1)
            let immune_no_pd1 = immune_cascade(&dead_cell_lps, n, &immune_params, false);
            let immune_with_pd1 = immune_cascade(&dead_cell_lps, n, &immune_params, true);

            eprintln!(
                "  {:<8} → Dead: {:6.2}%, LP/cell: {:5.1}, DAMP/cell: {:5.1}, T-kills: {:.0} / {:.0} (±PD1)",
                tx_name,
                death_rate * 100.0,
                avg_lp_at_death,
                immune_no_pd1.damp_per_dead_cell,
                immune_no_pd1.immune_kills,
                immune_with_pd1.immune_kills,
            );

            all_results.push(serde_json::json!({
                "phenotype": pheno_name,
                "treatment": tx_name,
                "n_cells": n,
                "n_dead": n_dead,
                "death_rate": death_rate,
                "avg_lp_at_death": avg_lp_at_death,
                "immune_no_pd1": {
                    "total_damps": immune_no_pd1.total_damps,
                    "damp_per_dead_cell": immune_no_pd1.damp_per_dead_cell,
                    "dc_activation_fraction": immune_no_pd1.dc_activation_fraction,
                    "mature_dcs": immune_no_pd1.mature_dcs,
                    "primed_tcells": immune_no_pd1.primed_tcells,
                    "immune_kills": immune_no_pd1.immune_kills,
                },
                "immune_with_pd1": {
                    "total_damps": immune_with_pd1.total_damps,
                    "damp_per_dead_cell": immune_with_pd1.damp_per_dead_cell,
                    "dc_activation_fraction": immune_with_pd1.dc_activation_fraction,
                    "mature_dcs": immune_with_pd1.mature_dcs,
                    "primed_tcells": immune_with_pd1.primed_tcells,
                    "immune_kills": immune_with_pd1.immune_kills,
                },
            }));
        }
        eprintln!();
    }

    // Save
    let json_path = args.output_dir.join("icd_comparison.json");
    write_json(&json_path, &all_results).expect("Failed to write JSON");

    eprintln!("=== Key Finding ===");
    // Compare SDT vs RSL3 DAMP per dead cell for persisters
    let sdt_result = all_results.iter().find(|r| {
        r["phenotype"] == "Persister" && r["treatment"] == "SDT"
    });
    let rsl3_result = all_results.iter().find(|r| {
        r["phenotype"] == "Persister" && r["treatment"] == "RSL3"
    });
    if let (Some(sdt), Some(rsl3)) = (sdt_result, rsl3_result) {
        let sdt_damp = sdt["immune_no_pd1"]["damp_per_dead_cell"].as_f64().unwrap_or(0.0);
        let rsl3_damp = rsl3["immune_no_pd1"]["damp_per_dead_cell"].as_f64().unwrap_or(0.0);
        if rsl3_damp > 0.0 {
            eprintln!(
                "SDT produces {:.1}× more DAMP per dead cell than RSL3 ({:.1} vs {:.1})",
                sdt_damp / rsl3_damp,
                sdt_damp,
                rsl3_damp,
            );
        }
        let sdt_kills = sdt["immune_with_pd1"]["immune_kills"].as_f64().unwrap_or(0.0);
        let rsl3_kills = rsl3["immune_with_pd1"]["immune_kills"].as_f64().unwrap_or(0.0);
        eprintln!(
            "SDT+anti-PD1 immune kills: {:.0} vs RSL3+anti-PD1: {:.0}",
            sdt_kills, rsl3_kills,
        );
    }

    eprintln!("\n=== Output saved to {} ===", args.output_dir.display());
}
