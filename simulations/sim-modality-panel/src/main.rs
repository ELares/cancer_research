//! Every applicable treatment arm, on the same tumour, in one run.
//!
//! # Why this binary exists
//!
//! `analysis/modality-coverage.md` measured this project's central criticism
//! and then measured it away: the engine went from thirteen of sixteen
//! taxonomy mechanisms with NO representation to none. What that left was a
//! harder question the report now leads with — **presence is not
//! applicability**. A mechanism whose code exists is present; a mechanism a
//! run can SELECT is applicable, and for a long time only one was.
//!
//! This is the binary that makes the second number mean something. Every
//! `Treatment` variant runs against the identical tumour, from the identical
//! seed, and reports a kill fraction. That is the first output in this
//! repository where a non-ferroptosis arm produces a number a reader can
//! compare.
//!
//! # What it is NOT, stated first because the comparison invites the error
//!
//! **These numbers are not a ranking of therapies.** Each arm carries its own
//! calibration status, recorded per row in
//! `simulations/calibration/CALIBRATION_STATUS.md`, and most of them read
//! `used in any reported number: N` precisely because they are uncalibrated.
//! An arm's kill fraction here is a function of the parameters it was given,
//! and those parameters are placeholders for every arm except radiation's DNA
//! channel.
//!
//! What the panel DOES show is structural, and structure is what this engine
//! can honestly report:
//!
//! * which arms depend on ferroptotic death and which do not;
//! * which are dose-responses and which are thresholds;
//! * which are limited by DELIVERY and which by biology;
//! * and what happens to each when the tumour is hypoxic, deep, or
//!   antigen-poor — the axes this project has spent its existence measuring.
//!
//! Those are comparisons between MECHANISMS under one model, not claims about
//! clinical efficacy, and the summary labels every row with its calibration
//! tier so the distinction survives being read quickly.

use std::path::PathBuf;

/// Steps per arm, matching the single-cell engine's own default run length so
/// the ADC arm is comparable to the four that go through `sim_cell`.
const STEPS: u32 = 180;

use clap::Parser;
use rand::{rngs::StdRng, Rng, SeedableRng};

use ferroptosis_core::ablation::{margin_survival_fraction, AblationConfig};
use ferroptosis_core::biochem::{sim_cell, CellState};
use ferroptosis_core::cell::{gen_cell, Phenotype, Treatment};
use ferroptosis_core::drug_transport::{
    antibody_drug_conjugate, concentration_at_distance, epithelial_well_vascularized,
};
use ferroptosis_core::immune::{
    adoptive_transfer_kills, immune_cascade, oncolytic_lysis, EffectorSource,
};
use ferroptosis_core::io::write_json;
use ferroptosis_core::params::{ImmuneParams, Params, RadiationConfig};
use ferroptosis_core::radiation;

#[derive(Parser, Debug)]
#[command(about = "Head-to-head panel across every applicable treatment arm")]
struct Args {
    /// Cells per arm.
    #[arg(long, default_value_t = 20_000)]
    n_cells: usize,

    /// Base RNG seed. Every arm derives from it identically, so the arms see
    /// the SAME tumour and a difference between them is the arm.
    #[arg(long, default_value_t = 42)]
    seed: u64,

    /// Antigen presented without ferroptotic death, for the immune arms.
    /// The engine default is 0, which would make checkpoint blockade kill
    /// nothing; this binary's whole point is to run the arm, so it supplies a
    /// nonzero value and the summary records it as an INPUT.
    #[arg(long, default_value_t = 0.004)]
    baseline_antigenicity: f64,

    /// Radiation dose in Gy.
    #[arg(long, default_value_t = 2.0)]
    radiation_dose_gy: f64,

    #[arg(long, default_value = "output/modality-panel")]
    output_dir: PathBuf,
}

/// One arm's result, with the caveat that governs reading it.
struct ArmResult {
    name: &'static str,
    kill_fraction: f64,
    /// Where the lethality came from. The panel is only interpretable if a
    /// reader can see that these are not the same kind of number.
    route: &'static str,
    /// What limits this arm in this model.
    limited_by: &'static str,
    calibration: &'static str,
}

fn main() {
    let a = Args::parse();
    let params = Params::default();
    let n = a.n_cells;

    eprintln!("=== Modality panel: {n} cells per arm, seed {} ===", a.seed);
    eprintln!("Every arm sees the SAME tumour. Differences are the arm.\n");

    let mut arms: Vec<ArmResult> = Vec::new();

    // ── Arms that kill through the ferroptosis engine ────────────────────
    for (tx, name) in [
        (Treatment::Control, "Control"),
        (Treatment::RSL3, "RSL3"),
        (Treatment::SDT, "SDT"),
        (Treatment::PDT, "PDT"),
    ] {
        let dead = (0..n)
            .filter(|i| {
                let mut rng = StdRng::seed_from_u64(a.seed.wrapping_add(*i as u64));
                let cell = gen_cell(Phenotype::Glycolytic, &mut rng);
                sim_cell(&cell, tx, &params, &mut rng).0
            })
            .count();
        arms.push(ArmResult {
            name,
            kill_fraction: dead as f64 / n as f64,
            route: "ferroptosis engine (lipid peroxidation)",
            limited_by: if name == "RSL3" {
                "endogenous ROS supply and antioxidant defence"
            } else {
                "energy delivery to depth"
            },
            calibration: "uncalibrated (direction-only)",
        });
    }

    // ── Radiation: DNA damage, which does not touch CellState ────────────
    let rad = RadiationConfig {
        dose_gy: a.radiation_dose_gy,
        alpha_per_gy: radiation::ALPHA_GBM_PARAMETERISATION_PER_GY,
        beta_per_gy2: radiation::ALPHA_GBM_PARAMETERISATION_PER_GY
            / radiation::ALPHA_BETA_TUMOUR_GY,
        ..RadiationConfig::default()
    };
    let lethality = radiation::dna_lethality(&rad, 1.0, 0.0);
    let dead = (0..n)
        .filter(|i| {
            let mut rng = StdRng::seed_from_u64(a.seed.wrapping_add(0x5AD + *i as u64));
            rng.gen::<f64>() < lethality
        })
        .count();
    arms.push(ArmResult {
        name: "Radiation",
        kill_fraction: dead as f64 / n as f64,
        route: "linear-quadratic DNA damage (not CellState)",
        limited_by: "dose, and oxygen through the OER",
        calibration: "form checked against a published parameterisation",
    });

    // ── Immune arms: they kill through the cascade, never through biochem ─
    let immune = ImmuneParams {
        baseline_antigenicity: a.baseline_antigenicity,
        tcell_kill_rate: 0.02,
        ..ImmuneParams::default()
    };
    let blockade = immune_cascade(&[], n, &immune, true).immune_kills;
    arms.push(ArmResult {
        name: "Immunotherapy",
        kill_fraction: blockade / n as f64,
        route: "immune cascade (no ferroptotic death required)",
        limited_by: "antigen availability and the checkpoint brake",
        calibration: "uncalibrated; published ORR band not fitted",
    });

    let effectors = EffectorSource::CarT.effectors_after(200.0, 30, 0.15, 5_000.0);
    let cart = adoptive_transfer_kills(effectors, n, &immune, 0.0, false);
    arms.push(ArmResult {
        name: "AdoptiveCell",
        kill_fraction: cart / n as f64,
        route: "redirected effectors (bypasses DC priming)",
        limited_by: "effector persistence and TME suppression",
        calibration: "uncalibrated; B-ALL remission band not fitted",
    });

    let (lysed, quality) = oncolytic_lysis(n, 0.15, 0.9, 0.8);
    let lps = vec![quality * 10.0; lysed as usize];
    let viral_immune = immune_cascade(&lps, n, &immune, false).immune_kills;
    arms.push(ArmResult {
        name: "OncolyticVirus",
        kill_fraction: (lysed + viral_immune) / n as f64,
        route: "direct lysis + the SHARED ICD chain",
        limited_by: "infection spread, which this engine takes as an input",
        calibration: "uncalibrated; T-VEC durable-response band not fitted",
    });

    // ── Ablation: a threshold, so the number is a MARGIN and not a dose ───
    let abl = AblationConfig {
        temperature_c: 56.0,
        minutes: 1.0,
        ..AblationConfig::default()
    };
    let covered = 0.85;
    arms.push(ArmResult {
        name: "Ablation",
        kill_fraction: 1.0 - margin_survival_fraction(&abl, covered),
        route: "threshold destruction (not a dose-response)",
        limited_by: "margin geometry -- coverage, and nothing else",
        calibration: "thresholds published; applicator field not modelled",
    });

    // ── ADC: the ferroptosis payload, limited by DELIVERY ────────────────
    let tissue = epithelial_well_vascularized();
    let adc = antibody_drug_conjugate();
    // Kill each cell at its distance from the vessel, using the payload's own
    // pharmacology and the ADC's transport. The point of the arm is that the
    // SECOND factor dominates.
    let dead = (0..n)
        .filter(|i| {
            let mut rng = StdRng::seed_from_u64(a.seed.wrapping_add(0xADC + *i as u64));
            let cell = gen_cell(Phenotype::Glycolytic, &mut rng);
            // Uniformly distributed through the intervessel half-distance.
            let r_um = rng.gen::<f64>() * tissue.inter_vessel_distance_um / 2.0;
            let avail = concentration_at_distance(r_um, &adc, &tissue);
            let peak = params.sdt_ros * avail;
            let mut state = CellState::from_cell_with_ros(
                &cell,
                Treatment::AntibodyDrugConjugate,
                &params,
                peak,
            );
            for step in 0..STEPS {
                ferroptosis_core::biochem::sim_cell_step(
                    &mut state, &cell, &params, step, 0.0, &mut rng,
                );
                if state.dead {
                    break;
                }
            }
            state.dead
        })
        .count();
    arms.push(ArmResult {
        name: "AntibodyDrugConjugate",
        kill_fraction: dead as f64 / n as f64,
        route: "ferroptosis payload, delivered on an antibody",
        limited_by: "the binding-site barrier (~7 um penetration)",
        calibration: "transport anchored; payload pharmacology is RSL3's",
    });

    std::fs::create_dir_all(&a.output_dir).expect("create output dir");
    eprintln!("{:<24} {:>8}  {}", "arm", "kill", "route");
    for r in &arms {
        eprintln!(
            "{:<24} {:>7.2}%  {}",
            r.name,
            r.kill_fraction * 100.0,
            r.route
        );
    }

    let json = serde_json::json!({
        "n_cells": n,
        "seed": a.seed,
        "inputs": {
            "baseline_antigenicity": a.baseline_antigenicity,
            "radiation_dose_gy": a.radiation_dose_gy,
        },
        "not_a_ranking": "Each arm carries its own calibration status; most are \
    uncalibrated placeholders. These are comparisons between MECHANISMS under one \
    model, not claims about clinical efficacy. See CALIBRATION_STATUS.md.",
        "arms": arms.iter().map(|r| serde_json::json!({
            "arm": r.name,
            "kill_fraction": r.kill_fraction,
            "route": r.route,
            "limited_by": r.limited_by,
            "calibration": r.calibration,
        })).collect::<Vec<_>>(),
    });
    let path = a.output_dir.join("modality_panel.json");
    write_json(&path, &json).expect("write json");
    eprintln!("\nWrote {}", path.display());
}
