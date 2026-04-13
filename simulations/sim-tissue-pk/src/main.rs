//! Tissue-specific drug penetration simulation.
//!
//! Computes drug concentration profiles and ferroptosis kill rates as a
//! function of distance from the nearest blood vessel, across different
//! tissue types.
//!
//! Usage: `cargo run --release --bin sim-tissue-pk`

use std::fs;
use std::path::Path;

use rand::prelude::*;
use serde::Serialize;

use ferroptosis_core::biochem::sim_cell;
use ferroptosis_core::cell::{gen_cell, Phenotype, Treatment};
use ferroptosis_core::drug_transport::{
    self, concentration_profile, max_distance_um,
    penetration_length_um, DrugParams, TissueParams,
};
use ferroptosis_core::params::Params;
use ferroptosis_core::stats::wilson_ci;

const N_RADIAL_BINS: usize = 50;
const N_CELLS_PER_BIN: usize = 1000;

#[derive(Serialize)]
struct BinResult {
    distance_um: f64,
    concentration: f64,
    death_rate: f64,
    ci_low: f64,
    ci_high: f64,
    n_cells: usize,
    n_dead: usize,
    tissue: String,
    drug: String,
}

#[derive(Serialize)]
struct Summary {
    tissue: String,
    drug: String,
    penetration_length_um: f64,
    max_distance_um: f64,
    vessel_wall_concentration: f64,
    vessel_wall_death_rate: f64,
    effective_kill_depth_um: f64,
    overall_kill_fraction: f64,
}

fn run_tissue_drug(
    drug: &DrugParams,
    tissue: &TissueParams,
    base_params: &Params,
    phenotype: Phenotype,
    seed: u64,
) -> (Vec<BinResult>, Summary) {
    let r_max = max_distance_um(tissue);
    let lambda = penetration_length_um(drug);
    let profile = concentration_profile(drug, tissue, N_RADIAL_BINS);

    let mut results = Vec::with_capacity(N_RADIAL_BINS);
    let mut total_dead = 0usize;
    let mut total_cells = 0usize;

    for (bin_idx, &(r_um, conc)) in profile.iter().enumerate() {
        // Scale GPX4 inhibition by local drug concentration
        let mut params = base_params.clone();
        params.rsl3_gpx4_inhib *= conc;

        let mut n_dead = 0usize;
        for i in 0..N_CELLS_PER_BIN {
            let cell_seed = seed + (bin_idx as u64) * (N_CELLS_PER_BIN as u64) * 2 + (i as u64) * 2;
            let mut rng = StdRng::seed_from_u64(cell_seed);
            let cell = gen_cell(phenotype, &mut rng);
            let mut sim_rng = StdRng::seed_from_u64(cell_seed + 1);
            let (dead, _, _, _) = sim_cell(&cell, Treatment::RSL3, &params, &mut sim_rng);
            if dead {
                n_dead += 1;
            }
        }

        let death_rate = n_dead as f64 / N_CELLS_PER_BIN as f64;
        let (ci_low, ci_high) = wilson_ci(N_CELLS_PER_BIN, n_dead);

        total_dead += n_dead;
        total_cells += N_CELLS_PER_BIN;

        results.push(BinResult {
            distance_um: r_um,
            concentration: conc,
            death_rate,
            ci_low,
            ci_high,
            n_cells: N_CELLS_PER_BIN,
            n_dead,
            tissue: tissue.name.to_string(),
            drug: drug.name.to_string(),
        });
    }

    // Find effective kill depth: distance where death rate drops below 10%
    let kill_depth = results
        .iter()
        .rev()
        .find(|r| r.death_rate >= 0.10)
        .map(|r| r.distance_um)
        .unwrap_or(0.0);

    let vessel_wall_dr = results.first().map(|r| r.death_rate).unwrap_or(0.0);
    let vessel_wall_conc = results.first().map(|r| r.concentration).unwrap_or(0.0);

    let summary = Summary {
        tissue: tissue.name.to_string(),
        drug: drug.name.to_string(),
        penetration_length_um: lambda,
        max_distance_um: r_max,
        vessel_wall_concentration: vessel_wall_conc,
        vessel_wall_death_rate: vessel_wall_dr,
        effective_kill_depth_um: kill_depth,
        overall_kill_fraction: total_dead as f64 / total_cells as f64,
    };

    (results, summary)
}

fn main() {
    eprintln!("=== Tissue-Specific Drug Penetration Simulation ===");
    eprintln!("Cells per radial bin: {N_CELLS_PER_BIN}");
    eprintln!("Radial bins: {N_RADIAL_BINS}");
    eprintln!("Phenotype: Persister (FSP1-low)\n");

    let base_params = Params::default();
    let seed: u64 = 42;

    let drugs: Vec<DrugParams> = vec![drug_transport::rsl3_like(), drug_transport::doxorubicin()];

    let tissues: Vec<TissueParams> = vec![
        drug_transport::epithelial_well_vascularized(),
        drug_transport::epithelial_poorly_vascularized(),
        drug_transport::neuroectodermal_cns(),
    ];

    let output_dir = Path::new("output/tissue-pk");
    fs::create_dir_all(output_dir).expect("Failed to create output directory");

    let mut all_results: Vec<BinResult> = Vec::new();
    let mut all_summaries: Vec<Summary> = Vec::new();

    for drug in &drugs {
        eprintln!("Drug: {} (λ = {:.1} μm)", drug.name, penetration_length_um(drug));
        for tissue in &tissues {
            let (results, summary) = run_tissue_drug(
                drug,
                tissue,
                &base_params,
                Phenotype::Persister,
                seed,
            );

            eprintln!(
                "  {}: kill depth = {:.0} μm, vessel-wall death = {:.1}%, overall = {:.1}%",
                tissue.name,
                summary.effective_kill_depth_um,
                summary.vessel_wall_death_rate * 100.0,
                summary.overall_kill_fraction * 100.0,
            );

            all_summaries.push(summary);
            all_results.extend(results);
        }
        eprintln!();
    }

    // Write CSV
    let csv_path = output_dir.join("tissue_pk_results.csv");
    let mut wtr = csv::Writer::from_path(&csv_path).expect("Failed to create CSV");
    for r in &all_results {
        wtr.serialize(r).expect("Failed to write CSV row");
    }
    wtr.flush().expect("Failed to flush CSV");
    eprintln!("Written: {}", csv_path.display());

    // Write JSON summary
    let json_path = output_dir.join("tissue_pk_summary.json");
    let json = serde_json::to_string_pretty(&all_summaries).expect("Failed to serialize JSON");
    fs::write(&json_path, json).expect("Failed to write JSON");
    eprintln!("Written: {}", json_path.display());

    // Print comparison table
    eprintln!("\n=== Summary ===\n");
    eprintln!(
        "{:<15} {:<35} {:<12} {:<12} {:<12} {:<12}",
        "Drug", "Tissue", "λ (μm)", "Kill depth", "Vessel DR", "Overall"
    );
    eprintln!("{}", "-".repeat(98));
    for s in &all_summaries {
        eprintln!(
            "{:<15} {:<35} {:<12.1} {:<12.0} {:>10.1}% {:>10.1}%",
            s.drug,
            s.tissue,
            s.penetration_length_um,
            s.effective_kill_depth_um,
            s.vessel_wall_death_rate * 100.0,
            s.overall_kill_fraction * 100.0,
        );
    }

    // Calibration check
    let dox_lambda = penetration_length_um(&drug_transport::doxorubicin());
    eprintln!("\n=== Calibration Check ===");
    eprintln!(
        "Doxorubicin penetration length: {:.1} μm (literature: 40-80 μm, Minchinton 2006)",
        dox_lambda
    );
    if dox_lambda >= 30.0 && dox_lambda <= 120.0 {
        eprintln!("  PASS — within calibration range");
    } else {
        eprintln!("  WARNING — outside expected range");
    }
}
