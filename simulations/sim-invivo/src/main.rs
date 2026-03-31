//! 2D vs In-Vivo Ferroptosis Comparison
//!
//! Tests the effect of SCD1/MUFA lipid remodeling on ferroptosis sensitivity.
//!
//! Runs three contexts with identical cell phenotypes:
//! 1. 2D baseline (MUFA off) — reproduces sim-original results
//! 2. In-vivo (MUFA on) — SCD1-driven MUFA protection active
//! 3. In-vivo + SCD1 inhibitor (MUFA off, same params otherwise) — resensitization control
//!
//! Then sweeps MUFA parameter space to characterize sensitivity.
//!
//! Key predictions to test:
//! - Dixon 2025: RSL3 (GPX4 inhibitor) should lose efficacy in vivo due to MUFA remodeling
//! - Exogenous ROS (SDT and PDT are modeled identically as a shared ROS burst): may still overwhelm MUFA defense
//! - SCD1 inhibition should restore ferroptosis sensitivity (Tesfay 2019)

use std::path::PathBuf;

use clap::Parser;
use rayon::prelude::*;
use rand::prelude::*;
use serde::Serialize;

use ferroptosis_core::cell::{gen_cell, Phenotype, Treatment};
use ferroptosis_core::biochem::sim_cell;
use ferroptosis_core::params::Params;
use ferroptosis_core::stats::wilson_ci;

#[derive(Parser)]
#[command(name = "sim-invivo", about = "2D vs in-vivo ferroptosis with SCD1/MUFA lipid remodeling")]
struct Args {
    /// Cells per condition for the main comparison.
    #[arg(long, default_value_t = 100_000)]
    n_cells: usize,

    /// Cells per condition for the parameter sweep.
    #[arg(long, default_value_t = 50_000)]
    n_sweep: usize,

    /// Random seed.
    #[arg(long, default_value_t = 42)]
    seed: u64,

    /// Output directory.
    #[arg(long, default_value = "output/invivo")]
    output_dir: PathBuf,
}

#[derive(Serialize, Clone, Debug)]
struct ComparisonResult {
    context: String,
    phenotype: String,
    treatment: String,
    n_cells: usize,
    n_dead: usize,
    death_rate: f64,
    ci_low: f64,
    ci_high: f64,
    mean_lp: f64,
    mean_gsh: f64,
    mean_gpx4: f64,
    scd_mufa_rate: f64,
    scd_mufa_max: f64,
    scd_mufa_decay: f64,
    initial_mufa_protection: f64,
}

#[derive(Serialize, Clone, Debug)]
struct SweepResult {
    scd_mufa_rate: f64,
    scd_mufa_max: f64,
    scd_mufa_decay: f64,
    initial_mufa_protection: f64,
    phenotype: String,
    treatment: String,
    n_cells: usize,
    death_rate: f64,
    ci_low: f64,
    ci_high: f64,
    protection_factor: f64,
}

fn run_context(
    context_name: &str,
    params: &Params,
    phenotypes: &[(Phenotype, &str)],
    treatments: &[(Treatment, &str)],
    n: usize,
    seed: u64,
) -> Vec<ComparisonResult> {
    let mut results = Vec::new();

    for (pheno, pname) in phenotypes {
        for (tx, tname) in treatments {
            let outcomes: Vec<(bool, f64, f64, f64)> = (0..n)
                .into_par_iter()
                .map(|i| {
                    let mut cell_rng = StdRng::seed_from_u64(seed.wrapping_add(i as u64 * 2));
                    let mut sim_rng = StdRng::seed_from_u64(seed.wrapping_add(i as u64 * 2 + 1));
                    let cell = gen_cell(*pheno, &mut cell_rng);
                    sim_cell(&cell, *tx, params, &mut sim_rng)
                })
                .collect();

            let dead = outcomes.iter().filter(|(d, _, _, _)| *d).count();
            let rate = dead as f64 / n as f64;
            let (ci_lo, ci_hi) = wilson_ci(n, dead);

            eprintln!(
                "  [{:<12}] {:<20} + {:<8} → Death: {:7.3}% [{:.3}-{:.3}]",
                context_name, pname, tname, rate * 100.0, ci_lo * 100.0, ci_hi * 100.0,
            );

            results.push(ComparisonResult {
                context: context_name.to_string(),
                phenotype: pname.to_string(),
                treatment: tname.to_string(),
                n_cells: n,
                n_dead: dead,
                death_rate: rate,
                ci_low: ci_lo,
                ci_high: ci_hi,
                mean_lp: outcomes.iter().map(|(_, l, _, _)| l).sum::<f64>() / n as f64,
                mean_gsh: outcomes.iter().map(|(_, _, g, _)| g).sum::<f64>() / n as f64,
                mean_gpx4: outcomes.iter().map(|(_, _, _, p)| p).sum::<f64>() / n as f64,
                scd_mufa_rate: params.scd_mufa_rate,
                scd_mufa_max: params.scd_mufa_max,
                scd_mufa_decay: params.scd_mufa_decay,
                initial_mufa_protection: params.initial_mufa_protection,
            });
        }
        eprintln!();
    }
    results
}

fn main() {
    let args = Args::parse();

    let phenotypes: Vec<(Phenotype, &str)> = vec![
        (Phenotype::Glycolytic, "Glycolytic"),
        (Phenotype::OXPHOS, "OXPHOS"),
        (Phenotype::Persister, "Persister (FSP1↓)"),
        (Phenotype::PersisterNrf2, "Persister+NRF2"),
    ];
    let treatments: Vec<(Treatment, &str)> = vec![
        (Treatment::Control, "Control"),
        (Treatment::RSL3, "RSL3"),
        (Treatment::SDT, "SDT"),
        (Treatment::PDT, "PDT"),
    ];

    std::fs::create_dir_all(&args.output_dir).expect("Failed to create output dir");

    // ============================================================
    // Part 1: Three-context comparison
    // ============================================================
    eprintln!("=== 2D vs In-Vivo Ferroptosis Comparison ===");
    eprintln!("Cells per condition: {}\n", args.n_cells);

    let params_2d = Params::default();
    let params_invivo = Params::invivo();
    // SCD1 inhibitor in vivo: cells start with pre-existing membrane MUFA
    // but SCD1 is blocked (rate=0). Natural lipid turnover (decay) gradually
    // depletes existing MUFA. This is NOT identical to 2D (no initial MUFA)
    // or in-vivo (steady-state maintenance). It models the acute phase after
    // SCD1 inhibitor administration.
    let params_scd1i = Params {
        scd_mufa_rate: 0.0,
        ..params_invivo.clone()
    };

    eprintln!("--- Context: 2D (MUFA off, scd_mufa_rate=0) ---");
    let results_2d = run_context("2d", &params_2d, &phenotypes, &treatments, args.n_cells, args.seed);

    eprintln!("--- Context: In-Vivo (MUFA on, rate={}, max={}) ---", params_invivo.scd_mufa_rate, params_invivo.scd_mufa_max);
    let results_invivo = run_context("invivo", &params_invivo, &phenotypes, &treatments, args.n_cells, args.seed);

    eprintln!("--- Context: In-Vivo + SCD1i (MUFA off, scd_mufa_rate=0) ---");
    let results_scd1i = run_context("invivo+scd1i", &params_scd1i, &phenotypes, &treatments, args.n_cells, args.seed);

    // ============================================================
    // Part 2: Protection factor analysis
    // ============================================================
    eprintln!("=== Protection Factors (2D death rate / in-vivo death rate) ===\n");
    eprintln!("  {:<20} {:<8} {:>10} {:>10} {:>10} {:>10}",
             "Phenotype", "Tx", "2D %", "InVivo %", "SCD1i %", "Prot.Fac.");

    for r2d in &results_2d {
        let riv = results_invivo.iter()
            .find(|r| r.phenotype == r2d.phenotype && r.treatment == r2d.treatment)
            .unwrap();
        let rs = results_scd1i.iter()
            .find(|r| r.phenotype == r2d.phenotype && r.treatment == r2d.treatment)
            .unwrap();

        let pf = if riv.death_rate > 0.0001 {
            r2d.death_rate / riv.death_rate
        } else if r2d.death_rate > 0.0001 {
            f64::INFINITY
        } else {
            1.0
        };

        eprintln!("  {:<20} {:<8} {:>9.2}% {:>9.2}% {:>9.2}% {:>10.2}×",
                 r2d.phenotype, r2d.treatment,
                 r2d.death_rate * 100.0, riv.death_rate * 100.0, rs.death_rate * 100.0, pf);
    }

    // ============================================================
    // Part 3: Key biological predictions
    // ============================================================
    eprintln!("\n=== Key Biological Predictions ===\n");

    let pers_rsl3_2d = results_2d.iter().find(|r| r.phenotype.contains("Persister (") && r.treatment == "RSL3").unwrap();
    let pers_rsl3_iv = results_invivo.iter().find(|r| r.phenotype.contains("Persister (") && r.treatment == "RSL3").unwrap();
    let pers_sdt_2d = results_2d.iter().find(|r| r.phenotype.contains("Persister (") && r.treatment == "SDT").unwrap();
    let pers_sdt_iv = results_invivo.iter().find(|r| r.phenotype.contains("Persister (") && r.treatment == "SDT").unwrap();
    let pers_rsl3_s = results_scd1i.iter().find(|r| r.phenotype.contains("Persister (") && r.treatment == "RSL3").unwrap();

    let dixon_pred = pers_rsl3_iv.death_rate < pers_rsl3_2d.death_rate * 0.5;
    eprintln!("  Dixon 2025 (RSL3 fails in vivo): 2D={:.1}% → InVivo={:.1}% — {}",
             pers_rsl3_2d.death_rate * 100.0, pers_rsl3_iv.death_rate * 100.0,
             if dixon_pred { "CONFIRMED" } else { "NOT CONFIRMED" });

    let sdt_survives = pers_sdt_iv.death_rate > 0.50;
    eprintln!("  SDT survives MUFA defense: InVivo={:.1}% — {}",
             pers_sdt_iv.death_rate * 100.0,
             if sdt_survives { "YES (>50% kill)" } else { "WEAKENED (<50% kill)" });

    let scd1i_resens = pers_rsl3_s.death_rate > pers_rsl3_iv.death_rate * 1.5;
    eprintln!("  SCD1i resensitizes (Tesfay 2019): InVivo={:.1}% → SCD1i={:.1}% — {}",
             pers_rsl3_iv.death_rate * 100.0, pers_rsl3_s.death_rate * 100.0,
             if scd1i_resens { "CONFIRMED" } else { "NOT CONFIRMED" });

    // ============================================================
    // Part 4: MUFA parameter sweep (Persister × SDT and RSL3)
    //
    // Two sweep modes:
    //   "onset"        — cells start at mufa=0 (freshly entering 3D context)
    //   "steady-state" — cells start at analytical steady state for each rate/max/decay
    // ============================================================
    eprintln!("\n=== MUFA Parameter Sweep (Persister) ===");
    eprintln!("Cells per point: {}\n", args.n_sweep);

    let rates = [0.002, 0.005, 0.01, 0.02, 0.04];
    let maxes = [0.20, 0.30, 0.40, 0.50, 0.60];
    let decay = params_invivo.scd_mufa_decay;

    let sweep_treatments: [(Treatment, &str); 2] = [
        (Treatment::SDT, "SDT"),
        (Treatment::RSL3, "RSL3"),
    ];

    let baseline_sdt = pers_sdt_2d.death_rate;
    let baseline_rsl3 = pers_rsl3_2d.death_rate;

    let mut sweep_results: Vec<SweepResult> = Vec::new();

    let sweep_modes: [(&str, bool); 2] = [
        ("onset (mufa=0)", false),
        ("steady-state", true),
    ];

    for (mode_label, use_steady_state) in &sweep_modes {
        eprintln!("=== Sweep mode: {} ===\n", mode_label);

        for (sweep_tx, sweep_tx_name) in &sweep_treatments {
            let baseline = match *sweep_tx {
                Treatment::SDT => baseline_sdt,
                Treatment::RSL3 => baseline_rsl3,
                _ => 1.0,
            };

            eprintln!("--- Persister + {} ({}) ---", sweep_tx_name, mode_label);
            eprintln!("  {:>8} | {}", "rate\\max", maxes.iter().map(|m| format!("{:>8.2}", m)).collect::<Vec<_>>().join(" | "));
            eprintln!("  ---------+-{}", maxes.iter().map(|_| "----------").collect::<Vec<_>>().join("-+-"));

            for &rate in &rates {
                let mut row_strs = Vec::new();
                for &max_val in &maxes {
                    let initial = if *use_steady_state {
                        rate * max_val / (rate + decay * max_val)
                    } else {
                        0.0
                    };
                    let p = Params {
                        scd_mufa_rate: rate,
                        scd_mufa_max: max_val,
                        scd_mufa_decay: decay,
                        initial_mufa_protection: initial,
                        ..Params::default()
                    };
                    let n = args.n_sweep;
                    let outcomes: Vec<(bool, f64, f64, f64)> = (0..n)
                        .into_par_iter()
                        .map(|i| {
                            let mut cell_rng = StdRng::seed_from_u64(args.seed.wrapping_add(i as u64 * 2));
                            let mut sim_rng = StdRng::seed_from_u64(args.seed.wrapping_add(i as u64 * 2 + 1));
                            let cell = gen_cell(Phenotype::Persister, &mut cell_rng);
                            sim_cell(&cell, *sweep_tx, &p, &mut sim_rng)
                        })
                        .collect();

                    let dead = outcomes.iter().filter(|(d, _, _, _)| *d).count();
                    let dr = dead as f64 / n as f64;
                    let (ci_lo, ci_hi) = wilson_ci(n, dead);
                    let pf = if dr > 0.0001 { baseline / dr } else { f64::INFINITY };

                    row_strs.push(format!("{:>7.1}%", dr * 100.0));
                    sweep_results.push(SweepResult {
                        scd_mufa_rate: rate,
                        scd_mufa_max: max_val,
                        scd_mufa_decay: decay,
                        initial_mufa_protection: initial,
                        phenotype: "Persister (FSP1↓)".to_string(),
                        treatment: sweep_tx_name.to_string(),
                        n_cells: n,
                        death_rate: dr,
                        ci_low: ci_lo,
                        ci_high: ci_hi,
                        protection_factor: pf,
                    });
                }
                eprintln!("  {:>8.3} | {}", rate, row_strs.join(" | "));
            }
            eprintln!();
        }
    }

    // ============================================================
    // Save outputs
    // ============================================================
    let mut all_comparison: Vec<&ComparisonResult> = Vec::new();
    all_comparison.extend(results_2d.iter());
    all_comparison.extend(results_invivo.iter());
    all_comparison.extend(results_scd1i.iter());

    let json_path = args.output_dir.join("invivo_comparison.json");
    let json = serde_json::to_string_pretty(&all_comparison).unwrap();
    std::fs::write(&json_path, &json).expect("Failed to write comparison JSON");

    let sweep_path = args.output_dir.join("mufa_sweep.json");
    let sweep_json = serde_json::to_string_pretty(&sweep_results).unwrap();
    std::fs::write(&sweep_path, &sweep_json).expect("Failed to write sweep JSON");

    // CSV for the sweep
    let csv_path = args.output_dir.join("mufa_sweep.csv");
    let mut wtr = csv::Writer::from_path(&csv_path).expect("Failed to create sweep CSV");
    for r in &sweep_results {
        wtr.serialize(r).unwrap();
    }
    wtr.flush().unwrap();

    eprintln!("\n=== Output saved to {} ===", args.output_dir.display());
    eprintln!("  {}", json_path.display());
    eprintln!("  {}", sweep_path.display());
    eprintln!("  {}", csv_path.display());
}
