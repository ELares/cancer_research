//! Monte Carlo Ferroptosis Simulation: OXPHOS vs Glycolytic cells under SDT/RSL3

use rand::prelude::*;
use rand_distr::Normal;
use rayon::prelude::*;
use serde::Serialize;
use std::fs;

#[derive(Clone)]
struct Cell { iron: f64, gsh: f64, gpx4: f64, basal_ros: f64, lipid_unsat: f64, nrf2: f64 }

#[derive(Clone, Copy)] enum Treatment { Control, RSL3, SDT }
#[derive(Clone, Copy, PartialEq)] enum Phenotype { Glycolytic, OXPHOS, OXPHOSNrf2 }

#[derive(Serialize, Clone)]
struct Result { phenotype: String, treatment: String, n_cells: usize, n_dead: usize,
    death_rate: f64, ci_low: f64, ci_high: f64, mean_lp: f64, mean_gsh: f64 }

fn norm(rng: &mut StdRng, mean: f64, sd: f64) -> f64 {
    let d = Normal::new(mean, sd).unwrap();
    d.sample(rng)
}

fn gen_cell(p: Phenotype, rng: &mut StdRng) -> Cell {
    match p {
        Phenotype::Glycolytic => Cell {
            iron: norm(rng, 1.0, 0.4).max(0.2), gsh: norm(rng, 3.0, 0.8).max(0.5),
            gpx4: norm(rng, 1.0, 0.2).max(0.1), basal_ros: norm(rng, 0.3, 0.1).max(0.05),
            lipid_unsat: norm(rng, 1.0, 0.2).max(0.3), nrf2: norm(rng, 1.0, 0.2).max(0.1),
        },
        Phenotype::OXPHOS => Cell {
            iron: norm(rng, 3.5, 1.0).max(0.5), gsh: norm(rng, 2.5, 0.8).max(0.5),
            gpx4: norm(rng, 1.0, 0.2).max(0.1), basal_ros: norm(rng, 0.9, 0.3).max(0.1),
            lipid_unsat: norm(rng, 1.8, 0.3).max(0.5), nrf2: norm(rng, 1.0, 0.2).max(0.1),
        },
        Phenotype::OXPHOSNrf2 => Cell {
            iron: norm(rng, 3.5, 1.0).max(0.5), gsh: norm(rng, 8.0, 2.0).max(2.0),
            gpx4: norm(rng, 1.5, 0.3).max(0.3), basal_ros: norm(rng, 0.9, 0.3).max(0.1),
            lipid_unsat: norm(rng, 1.8, 0.3).max(0.5), nrf2: norm(rng, 3.0, 0.5).max(1.0),
        },
    }
}

fn sim_cell(cell: &Cell, tx: Treatment, rng: &mut StdRng) -> (bool, f64, f64) {
    let mut gsh = cell.gsh;
    let mut gpx4 = cell.gpx4;
    let mut lp: f64 = 0.0;

    let exo_ros: f64 = match tx {
        Treatment::Control => 0.0, Treatment::RSL3 => 0.0,
        Treatment::SDT => norm(rng, 5.0, 1.5).max(0.0),
    };
    if let Treatment::RSL3 = tx { gpx4 *= 0.15; } // 85% inhibition

    for step in 0..60_u32 {
        let fenton = cell.iron * 0.3 * norm(rng, 1.0, 0.1).max(0.0);
        let decay = 0.5_f64.powf(step as f64 / 15.0);
        let total_ros = cell.basal_ros + exo_ros * decay + fenton;
        let scav = (total_ros * 0.3 * gsh / (gsh + 1.0)).min(gsh * 0.1);
        gsh -= scav;
        gsh += cell.nrf2 * 0.02 * (10.0 - gsh).max(0.0);
        gsh = gsh.max(0.0);
        let unscav = (total_ros - scav).max(0.0);
        let repair = gpx4 * gsh / (gsh + 0.5) * 0.3;
        lp = (lp + unscav * cell.lipid_unsat * 0.1 - repair).max(0.0);
        lp += norm(rng, 0.0, 0.01);
        lp = lp.max(0.0);
    }
    (lp > 5.0, lp, gsh)
}

fn wilson(n: usize, k: usize) -> (f64, f64) {
    let (n, p, z) = (n as f64, k as f64 / n as f64, 1.96);
    let d = 1.0 + z * z / n;
    let c = (p + z * z / (2.0 * n)) / d;
    let s = z * ((p * (1.0 - p) / n + z * z / (4.0 * n * n)).sqrt()) / d;
    ((c - s).max(0.0), (c + s).min(1.0))
}

fn main() {
    let n: usize = 1_000_000;
    let phenos = [(Phenotype::Glycolytic, "Glycolytic"), (Phenotype::OXPHOS, "OXPHOS"),
                  (Phenotype::OXPHOSNrf2, "OXPHOS+NRF2")];
    let txs = [(Treatment::Control, "Control"), (Treatment::RSL3, "RSL3"),
               (Treatment::SDT, "SDT")];

    eprintln!("Monte Carlo Ferroptosis Simulation — {} cells/condition\n", n);
    let mut results = Vec::new();

    for (p, pn) in &phenos {
        for (t, tn) in &txs {
            let out: Vec<(bool, f64, f64)> = (0..n).into_par_iter().map(|i| {
                let mut rng = StdRng::seed_from_u64(42 + i as u64);
                let cell = gen_cell(*p, &mut rng);
                sim_cell(&cell, *t, &mut rng)
            }).collect();

            let dead = out.iter().filter(|(d,_,_)| *d).count();
            let dr = dead as f64 / n as f64;
            let (cl, ch) = wilson(n, dead);
            let ml = out.iter().map(|(_,l,_)| l).sum::<f64>() / n as f64;
            let mg = out.iter().map(|(_,_,g)| g).sum::<f64>() / n as f64;

            eprintln!("{:<15} + {:<8} → Death: {:6.2}% [{:.2}-{:.2}]  LP:{:.2} GSH:{:.2}",
                     pn, tn, dr*100.0, cl*100.0, ch*100.0, ml, mg);
            results.push(Result { phenotype: pn.to_string(), treatment: tn.to_string(),
                n_cells: n, n_dead: dead, death_rate: dr, ci_low: cl, ci_high: ch,
                mean_lp: ml, mean_gsh: mg });
        }
        eprintln!();
    }

    // Relative risk
    eprintln!("Relative Risk (OXPHOS / Glycolytic):");
    for (_, tn) in &txs {
        let g = results.iter().find(|r| r.phenotype == "Glycolytic" && r.treatment == *tn).unwrap();
        let o = results.iter().find(|r| r.phenotype == "OXPHOS" && r.treatment == *tn).unwrap();
        let rr = if g.death_rate > 0.0 { o.death_rate / g.death_rate } else { f64::INFINITY };
        eprintln!("  {}: RR = {:.2}x", tn, rr);
    }
    eprintln!("\nNRF2 Protection:");
    for (_, tn) in &txs {
        let o = results.iter().find(|r| r.phenotype == "OXPHOS" && r.treatment == *tn).unwrap();
        let c = results.iter().find(|r| r.phenotype == "OXPHOS+NRF2" && r.treatment == *tn).unwrap();
        let prot = if o.death_rate > 0.0 { (1.0 - c.death_rate / o.death_rate) * 100.0 } else { 0.0 };
        eprintln!("  {}: NRF2 reduces death by {:.1}%", tn, prot);
    }

    let json = serde_json::to_string_pretty(&results).unwrap();
    fs::write("simulation_results.json", &json).unwrap();
    println!("{}", json);
    eprintln!("\nSaved to simulation_results.json");
}
