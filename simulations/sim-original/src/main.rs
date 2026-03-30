//! Monte Carlo Ferroptosis Sensitivity Simulation v3 (workspace reproduction)
//!
//! This binary produces bitwise-identical output to the original monolithic main.rs.
//! It uses the extracted ferroptosis-core library for all types and functions.

use ferroptosis_core::cell::{Phenotype, Treatment};
use ferroptosis_core::params::Params;
use ferroptosis_core::stats::run_condition;

fn main() {
    let n: usize = 1_000_000;
    let params = Params::default();

    let phenotypes = [
        (Phenotype::Glycolytic, "Glycolytic"),
        (Phenotype::OXPHOS, "OXPHOS"),
        (Phenotype::Persister, "Persister (FSP1↓)"),
        (Phenotype::PersisterNrf2, "Persister+NRF2"),
    ];
    let treatments = [
        (Treatment::Control, "Control"),
        (Treatment::RSL3, "RSL3"),
        (Treatment::SDT, "SDT"),
        (Treatment::PDT, "PDT"),
    ];

    eprintln!("=== Monte Carlo Ferroptosis Simulation v3 ===");
    eprintln!("Features: autocatalytic LP propagation, dynamic GPX4, uncapped GSH depletion");
    eprintln!("Cells per condition: {}", n);
    eprintln!("Total: {} cells × {} conditions = {}\n", n, 16, n * 16);

    let mut results = Vec::new();

    for (pheno, pname) in &phenotypes {
        for (tx, tname) in &treatments {
            let r = run_condition(*pheno, *tx, &params, n, pname, tname);
            eprintln!("{:<20} + {:<8} → Death: {:7.3}% [{:.3}-{:.3}]  LP:{:.3} GSH:{:.2} GPX4:{:.3}",
                     pname, tname, r.death_rate*100.0, r.ci_low*100.0, r.ci_high*100.0,
                     r.mean_lipid_perox, r.mean_gsh_final, r.mean_gpx4_final);
            results.push(r);
        }
        eprintln!();
    }

    // === VALIDATIONS ===
    eprintln!("=== VALIDATION ===");

    let mut v1_pass = true;
    for r in results.iter().filter(|r| r.treatment == "Control") {
        let ok = r.death_rate < 0.02;
        eprintln!("  Baseline {}: {:.3}% — {}", r.phenotype, r.death_rate*100.0,
                 if ok { "PASS" } else { "FAIL ⚠" });
        if !ok { v1_pass = false; }
    }

    let rsl3_pers = results.iter().find(|r| r.phenotype.contains("Persister (") && r.treatment == "RSL3").unwrap();
    let v2_pass = rsl3_pers.death_rate > 0.05;
    eprintln!("  RSL3 kills persisters: {:.2}% — {}", rsl3_pers.death_rate*100.0,
             if v2_pass { "PASS (>5%)" } else { "FAIL ⚠ (should match Higuchi)" });

    let sdt_pers = results.iter().find(|r| r.phenotype.contains("Persister (") && r.treatment == "SDT").unwrap();
    let ctrl_pers = results.iter().find(|r| r.phenotype.contains("Persister (") && r.treatment == "Control").unwrap();
    let v3_pass = sdt_pers.death_rate > ctrl_pers.death_rate + 0.01;
    eprintln!("  SDT > Control for persisters: {:.2}% vs {:.2}% — {}",
             sdt_pers.death_rate*100.0, ctrl_pers.death_rate*100.0,
             if v3_pass { "PASS" } else { "FAIL ⚠" });

    let nrf2_sdt = results.iter().find(|r| r.phenotype.contains("NRF2") && r.treatment == "SDT").unwrap();
    let v4_pass = nrf2_sdt.death_rate < sdt_pers.death_rate;
    eprintln!("  NRF2 protects vs SDT: {:.2}% vs {:.2}% — {}",
             nrf2_sdt.death_rate*100.0, sdt_pers.death_rate*100.0,
             if v4_pass { "PASS" } else { "FAIL ⚠" });

    if v1_pass && v2_pass && v3_pass && v4_pass {
        eprintln!("\n  ALL VALIDATIONS PASSED ✓");
    } else {
        eprintln!("\n  SOME VALIDATIONS FAILED — parameters need tuning");
    }

    // === COMPARISONS ===
    eprintln!("\n=== Key Comparisons ===");
    for tx_name in ["RSL3", "SDT", "PDT"] {
        let g = results.iter().find(|r| r.phenotype == "Glycolytic" && r.treatment == tx_name).unwrap();
        let p = results.iter().find(|r| r.phenotype.contains("Persister (") && r.treatment == tx_name).unwrap();
        let rr = if g.death_rate > 0.0001 { p.death_rate / g.death_rate } else { f64::NAN };
        eprintln!("  {}: Persister {:.2}% vs Glycolytic {:.2}% (RR={:.1}×)",
                 tx_name, p.death_rate*100.0, g.death_rate*100.0, rr);
    }

    let sdt = results.iter().find(|r| r.phenotype.contains("Persister (") && r.treatment == "SDT").unwrap();
    let rsl3 = results.iter().find(|r| r.phenotype.contains("Persister (") && r.treatment == "RSL3").unwrap();
    eprintln!("  Persister: SDT {:.2}% vs RSL3 {:.2}%", sdt.death_rate*100.0, rsl3.death_rate*100.0);

    let pdt = results.iter().find(|r| r.phenotype.contains("Persister (") && r.treatment == "PDT").unwrap();
    eprintln!("  Persister: SDT {:.2}% vs PDT {:.2}% (expected similar)", sdt.death_rate*100.0, pdt.death_rate*100.0);

    eprintln!("\n=== NRF2 Compensation ===");
    for tx_name in ["RSL3", "SDT"] {
        let p = results.iter().find(|r| r.phenotype.contains("Persister (") && r.treatment == tx_name).unwrap();
        let n = results.iter().find(|r| r.phenotype.contains("NRF2") && r.treatment == tx_name).unwrap();
        let prot = if p.death_rate > 0.001 { (1.0 - n.death_rate / p.death_rate) * 100.0 } else { 0.0 };
        eprintln!("  {}: NRF2 protection = {:.1}% ({:.2}% → {:.2}%)",
                 tx_name, prot, p.death_rate*100.0, n.death_rate*100.0);
    }

    // === SENSITIVITY ANALYSIS ===
    eprintln!("\n=== Sensitivity Analysis ===");
    eprintln!("Perturbing each parameter ±50%, testing Persister > Glycolytic under SDT:\n");

    let param_names = [
        "fenton_rate", "gsh_scav", "lp_rate", "lp_propagation",
        "gpx4_rate", "fsp1_rate", "nrf2_gsh", "gpx4_degrad",
        "death_threshold", "sdt_ros", "rsl3_inhib",
    ];

    let n_sens = 100_000;
    let mut holds_count = 0;
    let total_tests = param_names.len() * 2;

    for pname in &param_names {
        for mult in [0.5, 1.5] {
            let mut p = Params::default();
            match *pname {
                "fenton_rate" => p.fenton_rate *= mult,
                "gsh_scav" => p.gsh_scav_efficiency *= mult,
                "lp_rate" => p.lp_rate *= mult,
                "lp_propagation" => p.lp_propagation *= mult,
                "gpx4_rate" => p.gpx4_rate *= mult,
                "fsp1_rate" => p.fsp1_rate *= mult,
                "nrf2_gsh" => p.nrf2_gsh_rate *= mult,
                "gpx4_degrad" => p.gpx4_degradation_by_ros *= mult,
                "death_threshold" => p.death_threshold *= mult,
                "sdt_ros" => p.sdt_ros *= mult,
                "rsl3_inhib" => p.rsl3_gpx4_inhib = (p.rsl3_gpx4_inhib * mult).min(0.99),
                _ => {}
            }

            let g = run_condition(Phenotype::Glycolytic, Treatment::SDT, &p, n_sens, "G", "SDT");
            let pe = run_condition(Phenotype::Persister, Treatment::SDT, &p, n_sens, "P", "SDT");
            let holds = pe.death_rate > g.death_rate;
            if holds { holds_count += 1; }

            eprintln!("  {} ×{:.1}: P={:.2}% G={:.2}% — {}",
                     pname, mult, pe.death_rate*100.0, g.death_rate*100.0,
                     if holds { "HOLDS" } else { "FAILS ⚠" });
        }
    }
    eprintln!("\n  Result holds in {}/{} conditions ({:.0}%)",
             holds_count, total_tests, holds_count as f64 / total_tests as f64 * 100.0);

    // Save
    let json = serde_json::to_string_pretty(&results).unwrap();
    std::fs::write("simulation_results.json", &json).unwrap();
    println!("{}", json);
    eprintln!("\nSaved to simulation_results.json");
}
