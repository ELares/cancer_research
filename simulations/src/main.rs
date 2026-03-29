//! Monte Carlo Ferroptosis Sensitivity Simulation
//!
//! Models the biochemical chain: ROS → GSH depletion → GPX4/FSP1 inactivation
//! → lipid peroxidation → ferroptotic cell death
//!
//! Compares 4 cell phenotypes × 4 treatments × sensitivity analysis
//!
//! PARAMETER SOURCES:
//! - Intracellular GSH: 1-10 mM (Forman et al., Free Radic Biol Med 2009)
//! - Labile iron pool: 0.2-1.5 µM normal, 2-6 µM overloaded (Kakhlon & Cabantchik, Free Radic Biol Med 2002)
//! - Mitochondrial ROS: 1-2% O2 leak as superoxide (Murphy, Biochem J 2009)
//! - GPX4 kcat: ~40 s⁻¹ for lipid-OOH (Ursini et al., Free Radic Biol Med 1995)
//! - RSL3 GPX4 inhibition: IC50 ~10-100 nM, modeled as 85% at therapeutic dose
//! - FSP1/DHODH as GPX4-independent suppressor (Bersuker et al., Nature 2019)
//! - Persister FSP1 downregulation (Higuchi et al., Sci Adv 2026, PMID:41481741)
//!
//! KEY CONSTRAINT: All phenotypes must show <1% death under Control (no treatment).
//! Differential sensitivity emerges only under treatment, reflecting biological
//! reality that OXPHOS cells are viable but have narrower redox margins.

use rand::prelude::*;
use rand_distr::Normal;
use rayon::prelude::*;
use serde::Serialize;
use std::collections::HashMap;
use std::fs;

fn norm(rng: &mut StdRng, mean: f64, sd: f64) -> f64 {
    Normal::new(mean, sd).unwrap().sample(rng)
}

#[derive(Clone, Debug)]
struct Cell {
    /// Labile iron pool (µM). Ref: 0.2-1.5 normal, 2-6 overloaded
    iron: f64,
    /// Glutathione (mM). Ref: 1-10 mM intracellular
    gsh: f64,
    /// GPX4 activity (relative). Ref: kcat ~40/s, normalized to 1.0
    gpx4: f64,
    /// FSP1 activity (relative). GPX4-independent ferroptosis suppressor.
    /// Ref: Bersuker et al. Nature 2019. Downregulated in persisters (Higuchi 2026).
    fsp1: f64,
    /// Basal mitochondrial ROS (relative). Ref: 1-2% O2 leak, higher in OXPHOS
    basal_ros: f64,
    /// Lipid unsaturation susceptibility. OXPHOS cells have more mitochondrial membranes.
    lipid_unsat: f64,
    /// NRF2 transcriptional activity. Drives GSH synthesis, GPX4 expression.
    nrf2: f64,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
enum Treatment { Control, RSL3, SDT, PDT }

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
enum Phenotype { Glycolytic, OXPHOS, Persister, PersisterNrf2 }

#[derive(Serialize, Clone, Debug)]
struct SimResult {
    phenotype: String,
    treatment: String,
    n_cells: usize,
    n_dead: usize,
    death_rate: f64,
    ci_low: f64,
    ci_high: f64,
    mean_lipid_perox: f64,
    mean_gsh_final: f64,
    mean_iron: f64,
}

/// Rate constants — each with literature justification
struct Params {
    /// Fenton reaction rate: iron-catalyzed ROS generation
    /// Ref: labile iron catalyzes OH• at ~10³-10⁴ M⁻¹s⁻¹ (Winterbourn, Toxicol Lett 1995)
    /// Normalized so that iron=1.0 + this rate does NOT kill cells at baseline
    fenton_rate: f64,
    /// GSH scavenging efficiency for ROS
    /// Ref: GSH + ROS → GSSG, k ~10⁷ M⁻¹s⁻¹ for OH• (Buxton et al., 1988)
    /// Normalized to model units
    gsh_scav_rate: f64,
    /// Lipid peroxidation rate from unscavenged ROS
    /// Ref: PUFA peroxidation propagation ~10¹-10² M⁻¹s⁻¹ (Porter et al., 1995)
    lipid_perox_rate: f64,
    /// GPX4 repair rate for lipid peroxides
    /// Ref: kcat ~40 s⁻¹ (Ursini et al., 1995), requires GSH as cofactor
    gpx4_repair_rate: f64,
    /// FSP1 (CoQ10-dependent) lipid radical trapping
    /// Ref: Bersuker et al. Nature 2019 — GPX4-independent pathway
    fsp1_repair_rate: f64,
    /// NRF2-driven GSH resynthesis rate
    /// Ref: GCL transcription by NRF2, GSH t½ ~2-4h (Lu, Mol Aspects Med 2009)
    nrf2_gsh_synth: f64,
    /// Death threshold: lipid peroxidation level causing membrane rupture
    death_threshold: f64,
    /// SDT exogenous ROS magnitude
    sdt_ros: f64,
    /// PDT exogenous ROS magnitude (same mechanism, similar dose)
    pdt_ros: f64,
    /// RSL3 GPX4 inhibition fraction (0.85 = 85% inhibition)
    rsl3_gpx4_inhib: f64,
}

impl Default for Params {
    fn default() -> Self {
        Params {
            fenton_rate: 0.08,       // CALIBRATED: low enough that OXPHOS control < 1% death
            gsh_scav_rate: 0.4,      // GSH efficiently scavenges ROS
            lipid_perox_rate: 0.05,  // Moderate lipid peroxidation from unscavenged ROS
            gpx4_repair_rate: 0.5,   // GPX4 effectively repairs lipid-OOH (high kcat)
            fsp1_repair_rate: 0.15,  // FSP1 provides ~30% of GPX4-equivalent protection
            nrf2_gsh_synth: 0.03,    // Slow GSH resynthesis (hours timescale)
            death_threshold: 8.0,    // Higher threshold = cells tolerate more damage before death
            sdt_ros: 6.0,           // Exogenous ROS burst from SDT
            pdt_ros: 6.0,           // Same mechanism as SDT (sensitizer + ROS)
            rsl3_gpx4_inhib: 0.85,  // 85% GPX4 inhibition at therapeutic RSL3 dose
        }
    }
}

fn gen_cell(pheno: Phenotype, rng: &mut StdRng) -> Cell {
    match pheno {
        // Glycolytic (Warburg): standard cancer cell, therapy-sensitive
        Phenotype::Glycolytic => Cell {
            iron: norm(rng, 1.0, 0.3).max(0.3),       // Low labile iron
            gsh: norm(rng, 5.0, 1.0).max(1.0),         // Normal GSH
            gpx4: norm(rng, 1.0, 0.15).max(0.3),       // Normal GPX4
            fsp1: norm(rng, 1.0, 0.15).max(0.3),       // Normal FSP1
            basal_ros: norm(rng, 0.3, 0.08).max(0.05),  // Low mitochondrial ROS
            lipid_unsat: norm(rng, 1.0, 0.15).max(0.4), // Normal membranes
            nrf2: norm(rng, 1.0, 0.15).max(0.3),        // Normal NRF2
        },
        // OXPHOS-switched (resistant): survived first-line therapy via metabolic switch
        // Higher iron (ETC iron-sulfur clusters), higher basal ROS, more lipid membranes
        // BUT homeostatic — antioxidant defense is sufficient at baseline
        Phenotype::OXPHOS => Cell {
            iron: norm(rng, 2.5, 0.6).max(0.5),         // 2.5× more iron (ETC demand)
            gsh: norm(rng, 4.0, 1.0).max(1.0),           // Slightly lower (oxidative pressure)
            gpx4: norm(rng, 1.0, 0.15).max(0.3),         // Normal GPX4
            fsp1: norm(rng, 1.0, 0.15).max(0.3),         // Normal FSP1
            basal_ros: norm(rng, 0.6, 0.15).max(0.1),    // 2× more ROS (active ETC)
            lipid_unsat: norm(rng, 1.5, 0.2).max(0.6),   // 1.5× more lipid (mito membranes)
            nrf2: norm(rng, 1.2, 0.2).max(0.4),          // Slightly elevated (stress response)
        },
        // Drug-tolerant persister: OXPHOS + FSP1 downregulation (Higuchi et al. 2026)
        Phenotype::Persister => Cell {
            iron: norm(rng, 2.5, 0.6).max(0.5),
            gsh: norm(rng, 3.5, 1.0).max(0.8),           // Lower GSH (mesenchymal state)
            gpx4: norm(rng, 0.8, 0.15).max(0.2),         // Slightly reduced
            fsp1: norm(rng, 0.3, 0.1).max(0.05),         // FSP1 DOWNREGULATED (key vulnerability)
            basal_ros: norm(rng, 0.6, 0.15).max(0.1),
            lipid_unsat: norm(rng, 1.5, 0.2).max(0.6),
            nrf2: norm(rng, 0.8, 0.2).max(0.2),          // Reduced in mesenchymal state
        },
        // Persister with NRF2 compensation (the failure mode hypothesis)
        Phenotype::PersisterNrf2 => Cell {
            iron: norm(rng, 2.5, 0.6).max(0.5),
            gsh: norm(rng, 7.0, 1.5).max(2.0),           // NRF2-driven GSH elevation
            gpx4: norm(rng, 1.3, 0.2).max(0.4),          // NRF2 upregulates GPX4
            fsp1: norm(rng, 0.3, 0.1).max(0.05),         // FSP1 still down (different pathway)
            basal_ros: norm(rng, 0.6, 0.15).max(0.1),
            lipid_unsat: norm(rng, 1.5, 0.2).max(0.6),
            nrf2: norm(rng, 3.0, 0.5).max(1.0),
        },
    }
}

fn sim_cell(cell: &Cell, tx: Treatment, params: &Params, rng: &mut StdRng) -> (bool, f64, f64) {
    let mut gsh = cell.gsh;
    let mut gpx4 = cell.gpx4;
    let fsp1 = cell.fsp1;
    let mut lp: f64 = 0.0;

    // Treatment-specific effects
    let exo_ros: f64 = match tx {
        Treatment::Control | Treatment::RSL3 => 0.0,
        Treatment::SDT => norm(rng, params.sdt_ros, 1.5).max(0.0),
        Treatment::PDT => norm(rng, params.pdt_ros, 1.5).max(0.0),
    };

    if let Treatment::RSL3 = tx {
        gpx4 *= 1.0 - params.rsl3_gpx4_inhib;
    }

    // 120 time steps (representing ~2 hours post-treatment)
    for step in 0..120_u32 {
        // Fenton reaction: iron catalyzes ROS from H2O2
        let fenton = cell.iron * params.fenton_rate * norm(rng, 1.0, 0.1).max(0.0);

        // Exogenous ROS decays with half-life ~20 steps (~20 min for SDT/PDT)
        let decay = 0.5_f64.powf(step as f64 / 20.0);
        let total_ros = cell.basal_ros + exo_ros * decay + fenton;

        // GSH scavenges ROS (Michaelis-Menten-like saturation)
        let scav = (total_ros * params.gsh_scav_rate * gsh / (gsh + 2.0)).min(gsh * 0.05);
        gsh -= scav;

        // NRF2-driven GSH resynthesis (slow, represents transcriptional response)
        let gsh_max = 10.0;
        gsh += cell.nrf2 * params.nrf2_gsh_synth * ((gsh_max - gsh) / gsh_max).max(0.0);
        gsh = gsh.max(0.0);

        // Unscavenged ROS → lipid peroxidation (depends on lipid unsaturation)
        let unscav = (total_ros - scav).max(0.0);
        let lp_generation = unscav * cell.lipid_unsat * params.lipid_perox_rate;

        // Lipid peroxide repair: GPX4 (GSH-dependent) + FSP1 (GSH-independent)
        let gpx4_repair = gpx4 * (gsh / (gsh + 1.0)) * params.gpx4_repair_rate;
        let fsp1_repair = fsp1 * params.fsp1_repair_rate;
        let total_repair = gpx4_repair + fsp1_repair;

        lp = (lp + lp_generation - total_repair).max(0.0);

        // Small stochastic noise
        lp += norm(rng, 0.0, 0.005);
        lp = lp.max(0.0);
    }

    (lp > params.death_threshold, lp, gsh)
}

fn wilson_ci(n: usize, k: usize) -> (f64, f64) {
    let (nf, p, z) = (n as f64, k as f64 / n as f64, 1.96);
    let d = 1.0 + z * z / nf;
    let c = (p + z * z / (2.0 * nf)) / d;
    let s = z * ((p * (1.0 - p) / nf + z * z / (4.0 * nf * nf)).sqrt()) / d;
    ((c - s).max(0.0), (c + s).min(1.0))
}

fn run_condition(pheno: Phenotype, tx: Treatment, params: &Params, n: usize,
                 pname: &str, tname: &str) -> SimResult {
    let outcomes: Vec<(bool, f64, f64)> = (0..n).into_par_iter().map(|i| {
        let mut rng = StdRng::seed_from_u64(42 + i as u64);
        let cell = gen_cell(pheno, &mut rng);
        sim_cell(&cell, tx, params, &mut rng)
    }).collect();

    let dead = outcomes.iter().filter(|(d,_,_)| *d).count();
    let dr = dead as f64 / n as f64;
    let (cl, ch) = wilson_ci(n, dead);

    // Calculate mean iron for this phenotype
    let mean_iron: f64 = (0..1000).into_par_iter().map(|i| {
        let mut rng = StdRng::seed_from_u64(99999 + i as u64);
        gen_cell(pheno, &mut rng).iron
    }).sum::<f64>() / 1000.0;

    SimResult {
        phenotype: pname.to_string(),
        treatment: tname.to_string(),
        n_cells: n, n_dead: dead, death_rate: dr,
        ci_low: cl, ci_high: ch,
        mean_lipid_perox: outcomes.iter().map(|(_,l,_)| l).sum::<f64>() / n as f64,
        mean_gsh_final: outcomes.iter().map(|(_,_,g)| g).sum::<f64>() / n as f64,
        mean_iron,
    }
}

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

    eprintln!("=== Monte Carlo Ferroptosis Simulation ===");
    eprintln!("Cells per condition: {}", n);
    eprintln!("Phenotypes: 4 (Glycolytic, OXPHOS, Persister, Persister+NRF2)");
    eprintln!("Treatments: 4 (Control, RSL3, SDT, PDT)");
    eprintln!("Total cells simulated: {}\n", n * 16);

    let mut results = Vec::new();

    for (pheno, pname) in &phenotypes {
        for (tx, tname) in &treatments {
            let r = run_condition(*pheno, *tx, &params, n, pname, tname);
            eprintln!("{:<20} + {:<8} → Death: {:6.3}% [{:.3}-{:.3}]  LP:{:.3} GSH:{:.3} Fe:{:.2}",
                     pname, tname, r.death_rate*100.0, r.ci_low*100.0, r.ci_high*100.0,
                     r.mean_lipid_perox, r.mean_gsh_final, r.mean_iron);
            results.push(r);
        }
        eprintln!();
    }

    // Validate: Control death rate must be < 1% for all phenotypes
    eprintln!("=== VALIDATION: Baseline Viability ===");
    let mut valid = true;
    for r in results.iter().filter(|r| r.treatment == "Control") {
        let ok = r.death_rate < 0.01;
        eprintln!("  {} Control: {:.3}% — {}", r.phenotype, r.death_rate * 100.0,
                 if ok { "PASS" } else { "FAIL" });
        if !ok { valid = false; }
    }
    if !valid {
        eprintln!("  WARNING: Baseline viability check FAILED. Parameters need recalibration.");
    }

    // Relative risk analysis
    eprintln!("\n=== Relative Risk (vs Glycolytic) ===");
    for (_, tname) in &treatments {
        let glyco = results.iter().find(|r| r.phenotype == "Glycolytic" && r.treatment == *tname).unwrap();
        for (_, pname) in &phenotypes {
            if *pname == "Glycolytic" { continue; }
            let other = results.iter().find(|r| r.phenotype == *pname && r.treatment == *tname).unwrap();
            let rr = if glyco.death_rate > 0.001 { other.death_rate / glyco.death_rate } else { f64::NAN };
            eprintln!("  {} + {}: RR = {:.2}× ({:.2}% vs {:.2}%)",
                     pname, tname, rr, other.death_rate*100.0, glyco.death_rate*100.0);
        }
    }

    // SDT vs RSL3 comparison
    eprintln!("\n=== SDT vs RSL3 (head-to-head) ===");
    for (_, pname) in &phenotypes {
        let sdt = results.iter().find(|r| r.phenotype == *pname && r.treatment == "SDT").unwrap();
        let rsl3 = results.iter().find(|r| r.phenotype == *pname && r.treatment == "RSL3").unwrap();
        eprintln!("  {}: SDT {:.2}% vs RSL3 {:.2}%", pname, sdt.death_rate*100.0, rsl3.death_rate*100.0);
    }

    // SDT vs PDT comparison (should be similar — same mechanism)
    eprintln!("\n=== SDT vs PDT (depth-irrelevant biochemistry) ===");
    for (_, pname) in &phenotypes {
        let sdt = results.iter().find(|r| r.phenotype == *pname && r.treatment == "SDT").unwrap();
        let pdt = results.iter().find(|r| r.phenotype == *pname && r.treatment == "PDT").unwrap();
        eprintln!("  {}: SDT {:.2}% vs PDT {:.2}% (expected similar)", pname,
                 sdt.death_rate*100.0, pdt.death_rate*100.0);
    }

    // NRF2 protection
    eprintln!("\n=== NRF2 Compensation Effect ===");
    for (_, tname) in &treatments {
        let pers = results.iter().find(|r| r.phenotype == "Persister (FSP1↓)" && r.treatment == *tname).unwrap();
        let nrf2 = results.iter().find(|r| r.phenotype == "Persister+NRF2" && r.treatment == *tname).unwrap();
        let prot = if pers.death_rate > 0.001 {
            (1.0 - nrf2.death_rate / pers.death_rate) * 100.0
        } else { 0.0 };
        eprintln!("  {}: NRF2 reduces death by {:.1}% ({:.2}% → {:.2}%)",
                 tname, prot, pers.death_rate*100.0, nrf2.death_rate*100.0);
    }

    // === SENSITIVITY ANALYSIS ===
    eprintln!("\n=== Sensitivity Analysis (±50% on each parameter) ===");
    eprintln!("Testing whether Persister > Glycolytic under SDT is robust:\n");

    let base_params = Params::default();
    let param_names = [
        "fenton_rate", "gsh_scav_rate", "lipid_perox_rate", "gpx4_repair_rate",
        "fsp1_repair_rate", "nrf2_gsh_synth", "death_threshold", "sdt_ros",
    ];

    let n_sens = 100_000; // Fewer cells for sensitivity (speed)
    for pname in &param_names {
        for mult in [0.5, 1.5] {
            let mut p = Params::default();
            match *pname {
                "fenton_rate" => p.fenton_rate *= mult,
                "gsh_scav_rate" => p.gsh_scav_rate *= mult,
                "lipid_perox_rate" => p.lipid_perox_rate *= mult,
                "gpx4_repair_rate" => p.gpx4_repair_rate *= mult,
                "fsp1_repair_rate" => p.fsp1_repair_rate *= mult,
                "nrf2_gsh_synth" => p.nrf2_gsh_synth *= mult,
                "death_threshold" => p.death_threshold *= mult,
                "sdt_ros" => p.sdt_ros *= mult,
                _ => {}
            }

            let g = run_condition(Phenotype::Glycolytic, Treatment::SDT, &p, n_sens, "G", "SDT");
            let pe = run_condition(Phenotype::Persister, Treatment::SDT, &p, n_sens, "P", "SDT");
            let holds = pe.death_rate > g.death_rate;
            let rr = if g.death_rate > 0.001 { pe.death_rate / g.death_rate } else { f64::NAN };

            eprintln!("  {} ×{:.1}: Pers {:.2}% vs Glyco {:.2}% (RR={:.2}×) — {}",
                     pname, mult, pe.death_rate*100.0, g.death_rate*100.0, rr,
                     if holds { "HOLDS" } else { "FAILS" });
        }
    }

    // Save results
    let json = serde_json::to_string_pretty(&results).unwrap();
    fs::write("simulation_results.json", &json).unwrap();
    println!("{}", json);
    eprintln!("\nResults saved to simulation_results.json");
}
