//! Basic usage of the ferroptosis-core library.
//!
//! Generates a persister cell and simulates ferroptosis under different
//! treatments, printing the outcome for each.
//!
//! Run with: `cargo run -p ferroptosis-core --example basic_usage`

use rand::prelude::*;

use ferroptosis_core::biochem::sim_cell;
use ferroptosis_core::cell::{gen_cell, Phenotype, Treatment};
use ferroptosis_core::params::Params;

fn main() {
    let params = Params::default();
    let seed: u64 = 42;

    let phenotypes = [
        (Phenotype::Glycolytic, "Glycolytic"),
        (Phenotype::Persister, "Persister (FSP1-low)"),
        (Phenotype::PersisterNrf2, "Persister + NRF2"),
    ];
    let treatments = [
        (Treatment::Control, "Control"),
        (Treatment::RSL3, "RSL3"),
        (Treatment::SDT, "SDT"),
    ];

    println!("ferroptosis-core v{}", env!("CARGO_PKG_VERSION"));
    println!("Parameters: default (2D culture)\n");
    println!("{:<25} {:<10} {:<8} {:<10} {:<10} {:<10}", "Phenotype", "Treatment", "Dead?", "LP", "GSH", "GPX4");
    println!("{}", "-".repeat(73));

    for (pheno, pname) in &phenotypes {
        for (tx, tname) in &treatments {
            let mut rng = StdRng::seed_from_u64(seed);
            let cell = gen_cell(*pheno, &mut rng);
            let mut sim_rng = StdRng::seed_from_u64(seed + 1);
            let (dead, lp, gsh, gpx4) = sim_cell(&cell, *tx, &params, &mut sim_rng);
            println!(
                "{:<25} {:<10} {:<8} {:<10.3} {:<10.3} {:<10.3}",
                pname, tname, if dead { "YES" } else { "no" }, lp, gsh, gpx4
            );
        }
    }

    // Demonstrate in-vivo parameters (MUFA protection) with a small population
    println!("\n--- In-vivo context (SCD1-driven MUFA protection) ---\n");
    let invivo_params = Params::invivo();
    println!("MUFA steady-state protection: {:.2}", invivo_params.initial_mufa_protection);
    let n = 500;
    let mut deaths_2d = 0;
    let mut deaths_vivo = 0;
    for i in 0..n {
        let mut rng = StdRng::seed_from_u64(i * 2);
        let cell = gen_cell(Phenotype::Persister, &mut rng);
        let mut sr = StdRng::seed_from_u64(i * 2 + 1);
        if sim_cell(&cell, Treatment::RSL3, &params, &mut sr).0 { deaths_2d += 1; }
        let mut sr = StdRng::seed_from_u64(i * 2 + 1);
        if sim_cell(&cell, Treatment::RSL3, &invivo_params, &mut sr).0 { deaths_vivo += 1; }
    }
    println!("Persister + RSL3 (2D):     {deaths_2d}/{n} dead ({:.0}%)", deaths_2d as f64 / n as f64 * 100.0);
    println!("Persister + RSL3 (in-vivo): {deaths_vivo}/{n} dead ({:.0}%)", deaths_vivo as f64 / n as f64 * 100.0);
    println!("MUFA protection factor: {:.1}x", deaths_2d as f64 / deaths_vivo.max(1) as f64);
}
