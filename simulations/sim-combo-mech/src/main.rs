//! Mechanistic drug combination modeling.
//!
//! Runs pairwise drug combinations through the ferroptosis biochemistry
//! engine and computes Bliss-independence synergy scores. Unlike phenomenological
//! tools (SynergyFinder, CompuSyn), this reveals WHY combinations are synergistic
//! by tracing which pathway nodes each drug depletes.
//!
//! Usage: `cargo run --release --bin sim-combo-mech`

use std::fs;
use std::path::Path;

use rand::prelude::*;
use serde::Serialize;

use ferroptosis_core::biochem::{sim_cell_step, CellState};
use ferroptosis_core::cell::{gen_cell, norm, Cell, Phenotype};
use ferroptosis_core::params::Params;
use ferroptosis_core::stats::wilson_ci;

const N_CELLS: usize = 1000;
const N_STEPS: u32 = 180;

// ============================================================
// Drug effect model
// ============================================================

/// A drug's effect on the ferroptosis pathway, applied at initialization.
///
/// Each field targets a specific node. Combinations are formed by applying
/// multiple DrugEffects to the same cell — the biochemistry engine handles
/// the pathway coupling automatically.
#[derive(Clone)]
struct DrugEffect {
    name: &'static str,
    /// Fraction of GPX4 activity inhibited (0.0 = none, 0.92 = RSL3-level).
    gpx4_inhibition: f64,
    /// Fraction of FSP1 activity inhibited (0.0 = none, 0.85 = iFSP1-level).
    /// Also covers DHODH inhibitors (brequinar) since both reduce CoQ10-mediated repair.
    fsp1_inhibition: f64,
    /// Exogenous ROS peak dose (0.0 = none, 5.0 = SDT-level).
    exo_ros_dose: f64,
    /// Multiplier on basal ROS production (1.0 = none, 2.0 = HDACi doubles ROS).
    basal_ros_multiplier: f64,
}

fn rsl3() -> DrugEffect {
    DrugEffect {
        name: "RSL3",
        gpx4_inhibition: 0.92,
        fsp1_inhibition: 0.0,
        exo_ros_dose: 0.0,
        basal_ros_multiplier: 1.0,
    }
}

fn sdt() -> DrugEffect {
    DrugEffect {
        name: "SDT",
        gpx4_inhibition: 0.0,
        fsp1_inhibition: 0.0,
        exo_ros_dose: 5.0,
        basal_ros_multiplier: 1.0,
    }
}

/// FSP1 inhibitor (e.g., iFSP1/icFSP1).
/// ~80% lung tumor reduction as single agent (Wu et al., Nature 2025).
/// 85% inhibition is an estimated potency — not calibrated to a specific IC50.
fn fsp1i() -> DrugEffect {
    DrugEffect {
        name: "FSP1i",
        gpx4_inhibition: 0.0,
        fsp1_inhibition: 0.85,
        exo_ros_dose: 0.0,
        basal_ros_multiplier: 1.0,
    }
}

/// HDAC inhibitor (e.g., panobinostat).
/// Increases mitochondrial ROS in persisters by epigenetic derepression
/// of oxidative metabolism (Hangauer et al., Science Advances 2026).
/// 2x ROS multiplier is an estimate — the paper shows increased ROS without
/// quantifying the fold-change.
fn hdaci() -> DrugEffect {
    DrugEffect {
        name: "HDACi",
        gpx4_inhibition: 0.0,
        fsp1_inhibition: 0.0,
        exo_ros_dose: 0.0,
        basal_ros_multiplier: 2.0,
    }
}

// ============================================================
// Simulation engine
// ============================================================

/// Apply one or more drug effects to a cell, returning a modified Cell and
/// initialized CellState ready for the sim_cell_step loop.
fn apply_effects(
    cell: &Cell,
    effects: &[&DrugEffect],
    params: &Params,
    rng: &mut StdRng,
) -> (Cell, CellState) {
    let mut modified_cell = cell.clone();
    let mut gpx4 = cell.gpx4;
    let mut fsp1 = cell.fsp1;
    let mut exo_ros_total = 0.0_f64;

    for effect in effects {
        gpx4 *= 1.0 - effect.gpx4_inhibition;
        fsp1 *= 1.0 - effect.fsp1_inhibition;
        if effect.exo_ros_dose > 0.0 {
            exo_ros_total += norm(rng, effect.exo_ros_dose, 1.0).max(0.0);
        }
        modified_cell.basal_ros *= effect.basal_ros_multiplier;
    }

    let state = CellState {
        gsh: cell.gsh,
        gpx4,
        fsp1,
        mufa_protection: params.initial_mufa_protection,
        lp: 0.0,
        dead: false,
        death_step: None,
        exo_ros_peak: exo_ros_total,
    };

    (modified_cell, state)
}

/// Run N cells through the 180-step biochemistry with given drug effects.
/// Returns (n_dead, death_rate, ci_low, ci_high).
fn run_condition(
    effects: &[&DrugEffect],
    params: &Params,
    phenotype: Phenotype,
    n: usize,
    seed: u64,
) -> (usize, f64, f64, f64) {
    let mut n_dead = 0usize;

    for i in 0..n {
        let cell_seed = seed + (i as u64) * 2;
        let mut cell_rng = StdRng::seed_from_u64(cell_seed);
        let cell = gen_cell(phenotype, &mut cell_rng);

        let mut sim_rng = StdRng::seed_from_u64(cell_seed + 1);
        let (modified_cell, mut state) = apply_effects(&cell, effects, params, &mut sim_rng);

        for step in 0..N_STEPS {
            if sim_cell_step(&mut state, &modified_cell, params, step, 0.0, &mut sim_rng) {
                break;
            }
        }

        if state.dead {
            n_dead += 1;
        }
    }

    let rate = n_dead as f64 / n as f64;
    let (ci_low, ci_high) = wilson_ci(n, n_dead);
    (n_dead, rate, ci_low, ci_high)
}

// ============================================================
// Output types
// ============================================================

#[derive(Serialize)]
struct SingleResult {
    drug: String,
    death_rate: f64,
    ci_low: f64,
    ci_high: f64,
    n_dead: usize,
    n_cells: usize,
}

#[derive(Serialize)]
struct ComboResult {
    drug_a: String,
    drug_b: String,
    rate_a: f64,
    rate_b: f64,
    rate_combo: f64,
    bliss_prediction: f64,
    synergy_score: f64,
    ci_low: f64,
    ci_high: f64,
    n_dead: usize,
    n_cells: usize,
}

// ============================================================
// Main
// ============================================================

fn main() {
    eprintln!("=== Mechanistic Drug Combination Modeling ===");
    eprintln!("Cells per condition: {N_CELLS}");
    eprintln!("Phenotype: Persister (FSP1-low)");
    eprintln!("Context: 2D culture (default params)");
    eprintln!("NOTE: Drug potency parameters are estimated, not calibrated.\n");

    let params = Params::default();
    let seed: u64 = 42;
    let phenotype = Phenotype::Persister;

    let drugs: Vec<DrugEffect> = vec![rsl3(), sdt(), fsp1i(), hdaci()];

    // --- Single-drug baselines ---
    eprintln!("Single-drug baselines:");
    let mut single_rates: Vec<(String, f64)> = Vec::new();
    let mut single_results: Vec<SingleResult> = Vec::new();

    // Control (no drugs)
    let (n_dead, rate, ci_low, ci_high) = run_condition(&[], &params, phenotype, N_CELLS, seed);
    eprintln!("  Control: {:.1}% ({n_dead}/{N_CELLS})", rate * 100.0);
    single_rates.push(("Control".to_string(), rate));
    single_results.push(SingleResult {
        drug: "Control".to_string(), death_rate: rate, ci_low, ci_high, n_dead, n_cells: N_CELLS,
    });

    for drug in &drugs {
        let (n_dead, rate, ci_low, ci_high) =
            run_condition(&[drug], &params, phenotype, N_CELLS, seed);
        eprintln!("  {}: {:.1}% ({n_dead}/{N_CELLS})", drug.name, rate * 100.0);
        single_rates.push((drug.name.to_string(), rate));
        single_results.push(SingleResult {
            drug: drug.name.to_string(), death_rate: rate, ci_low, ci_high, n_dead, n_cells: N_CELLS,
        });
    }

    // --- Pairwise combinations ---
    eprintln!("\nPairwise combinations:");
    eprintln!(
        "{:<12} {:<12} {:>8} {:>8} {:>8} {:>8} {:>8}",
        "Drug A", "Drug B", "Rate A", "Rate B", "Combo", "Bliss", "Synergy"
    );
    eprintln!("{}", "-".repeat(76));

    let mut combo_results: Vec<ComboResult> = Vec::new();

    for i in 0..drugs.len() {
        for j in (i + 1)..drugs.len() {
            let drug_a = &drugs[i];
            let drug_b = &drugs[j];

            let rate_a = single_rates
                .iter()
                .find(|(n, _)| n == drug_a.name)
                .map(|(_, r)| *r)
                .unwrap();
            let rate_b = single_rates
                .iter()
                .find(|(n, _)| n == drug_b.name)
                .map(|(_, r)| *r)
                .unwrap();

            let (n_dead, rate_combo, ci_low, ci_high) =
                run_condition(&[drug_a, drug_b], &params, phenotype, N_CELLS, seed);

            let bliss = rate_a + rate_b - rate_a * rate_b;
            let synergy = if bliss > 0.001 {
                rate_combo / bliss
            } else {
                f64::NAN
            };

            eprintln!(
                "{:<12} {:<12} {:>7.1}% {:>7.1}% {:>7.1}% {:>7.1}% {:>8.2}",
                drug_a.name, drug_b.name,
                rate_a * 100.0, rate_b * 100.0,
                rate_combo * 100.0, bliss * 100.0,
                synergy,
            );

            combo_results.push(ComboResult {
                drug_a: drug_a.name.to_string(),
                drug_b: drug_b.name.to_string(),
                rate_a,
                rate_b,
                rate_combo,
                bliss_prediction: bliss,
                synergy_score: synergy,
                ci_low,
                ci_high,
                n_dead,
                n_cells: N_CELLS,
            });
        }
    }

    // --- Output ---
    let output_dir = Path::new("output/combo-mech");
    fs::create_dir_all(output_dir).expect("Failed to create output directory");

    let csv_path = output_dir.join("combo_synergy.csv");
    let mut wtr = csv::Writer::from_path(&csv_path).expect("Failed to create CSV");
    for r in &combo_results {
        wtr.serialize(r).expect("Failed to write CSV row");
    }
    wtr.flush().expect("Failed to flush CSV");
    eprintln!("\nWritten: {}", csv_path.display());

    let json_path = output_dir.join("combo_summary.json");
    let summary = serde_json::json!({
        "phenotype": "Persister (FSP1-low)",
        "context": "2D culture (default params)",
        "n_cells_per_condition": N_CELLS,
        "singles": single_results,
        "combinations": combo_results,
    });
    fs::write(&json_path, serde_json::to_string_pretty(&summary).unwrap())
        .expect("Failed to write JSON");
    eprintln!("Written: {}", json_path.display());

    // --- Summary ---
    eprintln!("\n=== Top Synergistic Pairs ===\n");
    let mut ranked: Vec<&ComboResult> = combo_results.iter().collect();
    ranked.sort_by(|a, b| {
        b.synergy_score
            .partial_cmp(&a.synergy_score)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    for r in &ranked {
        let label = if r.synergy_score > 1.1 {
            "SYNERGISTIC"
        } else if r.synergy_score < 0.9 {
            "ANTAGONISTIC"
        } else {
            "~ADDITIVE"
        };
        eprintln!(
            "  {} + {}: synergy={:.2} ({}) — actual {:.1}% vs Bliss {:.1}%",
            r.drug_a, r.drug_b, r.synergy_score, label,
            r.rate_combo * 100.0, r.bliss_prediction * 100.0,
        );
    }

    eprintln!("\nCaveats:");
    eprintln!("  - Drug potency parameters are estimated, not calibrated to specific IC50 values");
    eprintln!("  - Synergy scores depend on potency assumptions — directional findings are more robust than exact scores");
    eprintln!("  - All conditions use 2D culture params (no MUFA protection)");
    eprintln!("  - Bliss independence assumes drugs act on different targets; violations indicate shared pathway coupling");
}
