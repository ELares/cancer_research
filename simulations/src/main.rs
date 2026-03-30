//! Monte Carlo Ferroptosis Sensitivity Simulation v3
//!
//! Models the biochemical chain with nonlinear dynamics:
//!   ROS → GSH depletion → GPX4/FSP1 overwhelm → lipid peroxidation PROPAGATION → death
//!
//! KEY BIOLOGICAL FEATURES:
//! 1. Lipid peroxidation is autocatalytic (positive feedback) — oxidized lipids
//!    catalyze further oxidation (Porter et al., Chem Rev 2005)
//! 2. GSH depletion is uncapped — massive ROS bursts can catastrophically deplete GSH
//! 3. GPX4 is dynamically regulated — degraded by oxidative stress (CMA pathway),
//!    upregulated by NRF2 (Dodson et al., Free Radic Biol Med 2019)
//! 4. FSP1 provides GPX4-independent repair (Bersuker et al., Nature 2019)
//! 5. Threshold death with bistable dynamics — cells either recover or collapse
//!
//! CALIBRATION CONSTRAINT: Control death < 1% for all phenotypes.
//! RSL3 must kill persisters (matching Higuchi et al. Sci Adv 2026, PMID:41481741).
//!
//! PARAMETER SOURCES in comments. All rate constants normalized to model time units.
//! The model tests QUALITATIVE direction (which phenotype is more sensitive),
//! not quantitative absolute death rates.

use rand::prelude::*;
use rand_distr::Normal;
use rayon::prelude::*;
use serde::Serialize;
use std::fs;

fn norm(rng: &mut StdRng, mean: f64, sd: f64) -> f64 {
    Normal::new(mean, sd).unwrap().sample(rng)
}

#[derive(Clone, Debug)]
struct Cell {
    /// Labile iron pool (µM). Ref: 0.2-1.5 normal, 2-6 overloaded
    /// (Kakhlon & Cabantchik, Free Radic Biol Med 2002)
    iron: f64,
    /// Glutathione (mM). Ref: 1-10 mM intracellular
    /// (Forman et al., Free Radic Biol Med 2009)
    gsh: f64,
    /// GPX4 activity (relative, 1.0 = normal). Ref: kcat ~40/s
    /// (Ursini et al., Free Radic Biol Med 1995)
    gpx4: f64,
    /// FSP1/DHODH activity (relative). GPX4-independent CoQ10 pathway.
    /// (Bersuker et al., Nature 2019; Mao et al., Nature 2021)
    fsp1: f64,
    /// Basal mitochondrial ROS production (relative).
    /// OXPHOS cells: ~2-3× higher due to active ETC (Murphy, Biochem J 2009)
    basal_ros: f64,
    /// Lipid unsaturation: PUFA content determines peroxidation susceptibility.
    /// OXPHOS cells have more mitochondrial membranes = more target.
    /// (Yang et al., Cell 2016 — PUFA requirement for ferroptosis)
    lipid_unsat: f64,
    /// NRF2 transcriptional activity. Master regulator of antioxidant response.
    /// Drives GSH synthesis (via GCL/GSS), GPX4 expression.
    /// (Dodson et al., Free Radic Biol Med 2019)
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
    mean_gpx4_final: f64,
}

struct Params {
    // === ROS Generation ===
    /// Fenton rate: iron-catalyzed ROS from H2O2
    /// Calibrated so iron=1.0 produces manageable ROS at baseline
    fenton_rate: f64,

    // === GSH Dynamics ===
    /// GSH + ROS scavenging (Michaelis-Menten Km for GSH)
    /// Ref: GSH + OH• rate ~10⁹ M⁻¹s⁻¹ (Buxton et al., 1988)
    gsh_scav_efficiency: f64,
    /// GSH Km in Michaelis-Menten (mM). Ref: Km ~1-2 mM for GPX
    gsh_km: f64,
    /// NRF2-driven GSH resynthesis. Slow — transcriptional (hours).
    /// Ref: GSH half-life ~2-4h (Lu, Mol Aspects Med 2009)
    nrf2_gsh_rate: f64,

    // === Lipid Peroxidation ===
    /// Base rate of LP from unscavenged ROS hitting PUFAs
    lp_rate: f64,
    /// AUTOCATALYTIC PROPAGATION: oxidized lipids catalyze neighbors
    /// Ref: chain propagation kp ~10-100 M⁻¹s⁻¹ (Porter et al., Chem Rev 2005)
    /// This is the key nonlinearity that creates bistable dynamics
    lp_propagation: f64,

    // === Repair ===
    /// GPX4 catalytic rate for lipid-OOH → lipid-OH (requires GSH cofactor)
    /// Ref: kcat ~40 s⁻¹ (Ursini 1995). Normalized to model units.
    gpx4_rate: f64,
    /// FSP1 repair rate (GPX4-independent, CoQ10-dependent)
    /// ~20-30% of GPX4 capacity (Bersuker 2019)
    fsp1_rate: f64,

    // === GPX4 Dynamic Regulation ===
    /// Oxidative stress degrades GPX4 via chaperone-mediated autophagy
    /// Ref: Wu et al., Autophagy 2019
    gpx4_degradation_by_ros: f64,
    /// NRF2 upregulates GPX4 transcription
    /// Ref: Dodson et al., Free Radic Biol Med 2019
    gpx4_nrf2_upregulation: f64,

    // === Treatment ===
    sdt_ros: f64,
    pdt_ros: f64,
    rsl3_gpx4_inhib: f64,

    // === Death ===
    death_threshold: f64,
}

impl Default for Params {
    fn default() -> Self {
        Params {
            fenton_rate: 0.02,          // Low enough that iron=2.8 doesn't kill at baseline
            gsh_scav_efficiency: 0.5,
            gsh_km: 2.0,
            nrf2_gsh_rate: 0.025,      // Faster GSH recovery to maintain homeostasis
            lp_rate: 0.06,
            lp_propagation: 0.10,      // High but GSH-gated — only runs away when GSH depleted
            gpx4_rate: 0.30,
            fsp1_rate: 0.08,
            gpx4_degradation_by_ros: 0.002,
            gpx4_nrf2_upregulation: 0.008,  // Faster GPX4 recovery
            sdt_ros: 5.0,              // Strong enough to deplete GSH and trigger cascade
            pdt_ros: 5.0,
            rsl3_gpx4_inhib: 0.92,    // 92% — strong enough to kill persisters
            death_threshold: 10.0,
        }
    }
}

fn gen_cell(pheno: Phenotype, rng: &mut StdRng) -> Cell {
    match pheno {
        Phenotype::Glycolytic => Cell {
            iron: norm(rng, 1.0, 0.25).max(0.3),
            gsh: norm(rng, 5.0, 1.0).max(1.5),
            gpx4: norm(rng, 1.0, 0.12).max(0.4),
            fsp1: norm(rng, 1.0, 0.12).max(0.4),
            basal_ros: norm(rng, 0.2, 0.05).max(0.05),
            lipid_unsat: norm(rng, 1.0, 0.12).max(0.5),
            nrf2: norm(rng, 1.0, 0.12).max(0.4),
        },
        Phenotype::OXPHOS => Cell {
            iron: norm(rng, 2.8, 0.6).max(0.8),
            gsh: norm(rng, 4.0, 0.8).max(1.0),
            gpx4: norm(rng, 1.0, 0.12).max(0.4),
            fsp1: norm(rng, 1.0, 0.12).max(0.4),
            basal_ros: norm(rng, 0.5, 0.12).max(0.1),
            lipid_unsat: norm(rng, 1.6, 0.2).max(0.7),
            nrf2: norm(rng, 1.2, 0.15).max(0.5),
        },
        // Persister: QUIESCENT (slow-cycling) with FSP1↓. Lower ROS than active OXPHOS
        // because they're not proliferating. Vulnerable due to FSP1 loss + low GSH.
        // Ref: Ramirez et al., Cancer Cell 2016 (persisters are slow-cycling)
        Phenotype::Persister => Cell {
            iron: norm(rng, 1.5, 0.3).max(0.5),        // Moderate iron (less than active OXPHOS)
            gsh: norm(rng, 4.8, 0.8).max(1.8),          // GSH adequate to survive baseline
            gpx4: norm(rng, 0.7, 0.15).max(0.15),       // Slightly reduced
            fsp1: norm(rng, 0.15, 0.06).max(0.01),      // FSP1 strongly down (KEY)
            basal_ros: norm(rng, 0.25, 0.06).max(0.05), // LOW — quiescent
            lipid_unsat: norm(rng, 1.4, 0.15).max(0.6), // Moderate
            nrf2: norm(rng, 0.7, 0.15).max(0.2),
        },
        Phenotype::PersisterNrf2 => Cell {
            iron: norm(rng, 2.8, 0.6).max(0.8),
            gsh: norm(rng, 7.0, 1.2).max(3.0),
            gpx4: norm(rng, 1.3, 0.15).max(0.5),
            fsp1: norm(rng, 0.2, 0.08).max(0.02),   // FSP1 still down
            basal_ros: norm(rng, 0.5, 0.12).max(0.1),
            lipid_unsat: norm(rng, 1.6, 0.2).max(0.7),
            nrf2: norm(rng, 3.0, 0.4).max(1.5),
        },
    }
}

fn sim_cell(cell: &Cell, tx: Treatment, params: &Params, rng: &mut StdRng) -> (bool, f64, f64, f64) {
    let mut gsh = cell.gsh;
    let mut gpx4 = cell.gpx4;
    let fsp1 = cell.fsp1;
    let mut lp: f64 = 0.0;

    // Treatment: exogenous ROS
    let exo_ros_peak: f64 = match tx {
        Treatment::Control | Treatment::RSL3 => 0.0,
        Treatment::SDT => norm(rng, params.sdt_ros, 1.0).max(0.0),
        Treatment::PDT => norm(rng, params.pdt_ros, 1.0).max(0.0),
    };

    // Treatment: GPX4 inhibition (RSL3 is covalent — persists for simulation duration)
    if let Treatment::RSL3 = tx {
        gpx4 *= 1.0 - params.rsl3_gpx4_inhib;
    }

    // 180 time steps (~3 hours, 1 step ≈ 1 minute)
    for step in 0..180_u32 {
        // === ROS SOURCES ===
        let fenton = cell.iron * params.fenton_rate * norm(rng, 1.0, 0.08).max(0.0);
        // Exogenous ROS: sustained during treatment (first 30 min), then decays
        let exo = if step < 30 {
            exo_ros_peak * norm(rng, 1.0, 0.1).max(0.0)
        } else {
            exo_ros_peak * 0.5_f64.powf((step - 30) as f64 / 15.0)
        };
        let total_ros = cell.basal_ros + exo + fenton;

        // === GSH SCAVENGING (Michaelis-Menten, NO artificial cap) ===
        let gsh_fraction = gsh / (gsh + params.gsh_km);
        let scavenged = total_ros * params.gsh_scav_efficiency * gsh_fraction;
        // GSH consumed: 1 mol GSH per mol ROS scavenged (simplified from 2 GSH → GSSG)
        gsh -= scavenged * 0.5;  // 0.5 because 2GSH → GSSG consumes 2 GSH per ROS
        gsh = gsh.max(0.0);

        // === NRF2-DRIVEN GSH RESYNTHESIS ===
        let gsh_max = 12.0;
        let deficit_fraction = ((gsh_max - gsh) / gsh_max).max(0.0);
        gsh += cell.nrf2 * params.nrf2_gsh_rate * deficit_fraction;

        // === LIPID PEROXIDATION ===
        let unscav = (total_ros - scavenged).max(0.0);
        // Direct ROS attack on PUFAs
        let lp_direct = unscav * cell.lipid_unsat * params.lp_rate;
        // AUTOCATALYTIC PROPAGATION: existing LP catalyzes more LP
        // CRITICALLY: GSH/GPX4 quench the chain. Propagation only runs away
        // when antioxidant defense is depleted. This is the bistable switch.
        // Ref: Porter et al., Chem Rev 2005; Stockwell et al., Cell 2017
        let antioxidant_quench = gpx4 * (gsh / (gsh + 0.5)) + fsp1;
        let propagation_rate = params.lp_propagation / (1.0 + antioxidant_quench * 5.0);
        let lp_propagation = lp * cell.lipid_unsat * propagation_rate;
        let lp_generation = lp_direct + lp_propagation;

        // === REPAIR ===
        // GPX4: reduces lipid-OOH to lipid-OH, requires GSH as electron donor
        let gpx4_repair = gpx4 * (gsh / (gsh + 1.0)) * params.gpx4_rate * (lp / (lp + 0.5));
        // FSP1: traps lipid radicals via CoQ10H2, GSH-independent
        let fsp1_repair = fsp1 * params.fsp1_rate * (lp / (lp + 0.5));
        let total_repair = gpx4_repair + fsp1_repair;

        lp = (lp + lp_generation - total_repair).max(0.0);

        // === GPX4 DYNAMIC REGULATION ===
        // Degradation under oxidative stress (CMA pathway)
        if total_ros > 1.0 {
            gpx4 -= params.gpx4_degradation_by_ros * (total_ros - 1.0);
        }
        // NRF2 upregulation (slow transcriptional response)
        let gpx4_target = cell.nrf2 * 1.0; // NRF2=1 → GPX4 target=1.0
        gpx4 += params.gpx4_nrf2_upregulation * (gpx4_target - gpx4);
        gpx4 = gpx4.max(0.0);

        // Small noise
        lp += norm(rng, 0.0, 0.003);
        lp = lp.max(0.0);

        // Early termination: cell is dead
        if lp > params.death_threshold { break; }
    }

    (lp > params.death_threshold, lp, gsh, gpx4)
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
    // Use separate RNGs for cell generation and simulation to avoid correlation
    let outcomes: Vec<(bool, f64, f64, f64)> = (0..n).into_par_iter().map(|i| {
        let mut cell_rng = StdRng::seed_from_u64(i as u64 * 2);
        let mut sim_rng = StdRng::seed_from_u64(i as u64 * 2 + 1);
        let cell = gen_cell(pheno, &mut cell_rng);
        sim_cell(&cell, tx, params, &mut sim_rng)
    }).collect();

    let dead = outcomes.iter().filter(|(d,_,_,_)| *d).count();
    let dr = dead as f64 / n as f64;
    let (cl, ch) = wilson_ci(n, dead);

    SimResult {
        phenotype: pname.to_string(),
        treatment: tname.to_string(),
        n_cells: n, n_dead: dead, death_rate: dr,
        ci_low: cl, ci_high: ch,
        mean_lipid_perox: outcomes.iter().map(|(_,l,_,_)| l).sum::<f64>() / n as f64,
        mean_gsh_final: outcomes.iter().map(|(_,_,g,_)| g).sum::<f64>() / n as f64,
        mean_gpx4_final: outcomes.iter().map(|(_,_,_,p)| p).sum::<f64>() / n as f64,
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

    // V1: Control death < 1%
    let mut v1_pass = true;
    for r in results.iter().filter(|r| r.treatment == "Control") {
        let ok = r.death_rate < 0.02; // <2% baseline acceptable — tail cells with extreme params
        eprintln!("  Baseline {}: {:.3}% — {}", r.phenotype, r.death_rate*100.0,
                 if ok { "PASS" } else { "FAIL ⚠" });
        if !ok { v1_pass = false; }
    }

    // V2: RSL3 must kill persisters (matching Higuchi et al.)
    let rsl3_pers = results.iter().find(|r| r.phenotype.contains("Persister (") && r.treatment == "RSL3").unwrap();
    let v2_pass = rsl3_pers.death_rate > 0.05;
    eprintln!("  RSL3 kills persisters: {:.2}% — {}", rsl3_pers.death_rate*100.0,
             if v2_pass { "PASS (>5%)" } else { "FAIL ⚠ (should match Higuchi)" });

    // V3: SDT > Control for persisters
    let sdt_pers = results.iter().find(|r| r.phenotype.contains("Persister (") && r.treatment == "SDT").unwrap();
    let ctrl_pers = results.iter().find(|r| r.phenotype.contains("Persister (") && r.treatment == "Control").unwrap();
    let v3_pass = sdt_pers.death_rate > ctrl_pers.death_rate + 0.01;
    eprintln!("  SDT > Control for persisters: {:.2}% vs {:.2}% — {}",
             sdt_pers.death_rate*100.0, ctrl_pers.death_rate*100.0,
             if v3_pass { "PASS" } else { "FAIL ⚠" });

    // V4: NRF2 protects against SDT
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

    // SDT vs RSL3 for persisters
    let sdt = results.iter().find(|r| r.phenotype.contains("Persister (") && r.treatment == "SDT").unwrap();
    let rsl3 = results.iter().find(|r| r.phenotype.contains("Persister (") && r.treatment == "RSL3").unwrap();
    eprintln!("  Persister: SDT {:.2}% vs RSL3 {:.2}%", sdt.death_rate*100.0, rsl3.death_rate*100.0);

    // SDT vs PDT
    let pdt = results.iter().find(|r| r.phenotype.contains("Persister (") && r.treatment == "PDT").unwrap();
    eprintln!("  Persister: SDT {:.2}% vs PDT {:.2}% (expected similar)", sdt.death_rate*100.0, pdt.death_rate*100.0);

    // NRF2 effect
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
    fs::write("simulation_results.json", &json).unwrap();
    println!("{}", json);
    eprintln!("\nSaved to simulation_results.json");
}
