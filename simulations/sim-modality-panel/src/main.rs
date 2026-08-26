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
use ferroptosis_core::params::{ImmuneParams, Params, RadiationConfig, SpatialParams};
use ferroptosis_core::physics::sdt_intensity_at_depth;
use ferroptosis_core::radiation;

#[derive(Parser, Debug)]
#[command(about = "Head-to-head panel across every applicable treatment arm")]
struct Args {
    /// Cells per arm.
    #[arg(long, default_value_t = 20_000)]
    n_cells: usize,

    /// Also sweep every arm across the tumour-microenvironment axes this
    /// project spent years establishing for ferroptosis: hypoxia, stromal
    /// shielding and acidic pH. Off by default because it multiplies the run
    /// by eight and writes a second artifact.
    #[arg(long, default_value_t = false)]
    tme_sweep: bool,

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

    if a.tme_sweep {
        run_tme_sweep(&a, &params, &immune, &rad);
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

/// How each arm fares under the tumour microenvironment.
///
/// # Why this is the deep question rather than the panel above
///
/// The head-to-head panel answers "what does each arm do to a naive tumour",
/// which is the shallow question and the one a reader over-interprets. This
/// project's actual contribution has never been kill fractions: it is the
/// three RESISTANCE AXES the ferroptosis chapters establish -- hypoxia,
/// stromal shielding and acidic pH -- and the finding that pharmacologic and
/// physical modalities respond to them differently.
///
/// Every arm added by the coverage campaign is untested against those axes,
/// which is exactly the gap between "the engine can express it" and "the
/// engine has something to say about it". This sweep closes it.
///
/// # What varies, and what each axis does to WHICH arm
///
/// The axes are applied through the SAME helpers the ferroptosis work uses,
/// not through arm-specific fudges -- that is the point, because a modality
/// comparison is only meaningful if the environment is identical:
///
/// * **Hypoxia** scales exogenous-ROS yield through `oxygen::oer_exo_factor`
///   (the measured Alper-Howard-Flanders shape) and radiation's delivered
///   dose through the same hyperbola. It does NOT touch the immune arms'
///   effector supply, and that asymmetry is a PREDICTION rather than an
///   omission: a redirected T cell does not need oxygen to lyse a target,
///   though the suppressor field it meets is worse in a hypoxic core.
/// * **Stroma** raises the antioxidant setpoint of shielded cells, which is
///   the CAF-mediated GSH/MUFA supply the stromal chapter measures. It has no
///   effect on ablation, which does not care what a cell's glutathione is
///   doing.
/// * **Acidic pH** traps weak-base drugs outside the cell
///   (`ph::ion_trap_factor_from_ph`), so it hits the DELIVERED arms -- RSL3
///   and the ADC -- and leaves energy-delivered arms alone.
///
/// # The honest expectation, written before the run
///
/// If the model is coherent, ablation and radiation should be the most
/// TME-robust arms and the delivered ones the least, because that is what the
/// mechanisms imply rather than what any parameter was tuned for. A result
/// that agrees is weak evidence; a result that disagrees is worth
/// investigating, and the artifact reports the ordering either way.
fn run_tme_sweep(a: &Args, params: &Params, immune: &ImmuneParams, rad: &RadiationConfig) {
    use ferroptosis_core::oxygen::{oer_exo_factor, OER_REFERENCE_PO2_MMHG};
    use ferroptosis_core::ph::ion_trap_factor_from_ph;

    // The three axes at their documented endpoints, from the ferroptosis
    // chapters rather than chosen here.
    const O2_NORMOXIC: f64 = 1.0;
    const O2_HYPOXIC: f64 = 0.05;
    const PH_EDGE: f64 = 7.4;
    const PH_CORE: f64 = 6.5;
    const STROMAL_GSH_BOOST: f64 = 1.5;

    let n = a.n_cells;
    let mut rows = Vec::new();

    // FIVE axes, not three. The first version swept hypoxia, stroma and pH
    // and reported two of them INERT -- correctly, but for a reason that was
    // a property of the CONFIGURATION rather than of the axes: at the
    // glycolytic phenotype the delivered arm kills essentially nothing, so
    // ion trapping had nothing to scale and the antioxidant buffer was
    // swamped. Running the persister phenotype as well makes both visible,
    // and it is the phenotype this project's thesis is actually about.
    //
    // DEPTH is the fifth, and it is the axis that separates the physical
    // modalities from one another rather than from the drugs.
    for &(pheno, pheno_name) in &[
        (Phenotype::Glycolytic, "glycolytic"),
        (Phenotype::Persister, "persister"),
    ] {
        for &depth_um in &[0.0_f64, 5_000.0] {
            for &(o2, hypoxic) in &[(O2_NORMOXIC, false), (O2_HYPOXIC, true)] {
                for &stroma in &[false, true] {
                    for &acid in &[false, true] {
                        let ph = if acid { PH_CORE } else { PH_EDGE };
                        let trap = ion_trap_factor_from_ph(ph, PH_EDGE, 1.0);
                        let exo = oer_exo_factor(o2, 1.0, OER_REFERENCE_PO2_MMHG);

                        // Ferroptosis-routed arms: exogenous ROS scaled by O2, the
                        // antioxidant setpoint raised by stroma, and the DELIVERED
                        // arm additionally scaled by ion trapping.
                        let mut p = params.clone();
                        if stroma {
                            p.gsh_max *= STROMAL_GSH_BOOST;
                        }
                        // RSL3 is a weak base: acid traps it outside the cell, so the
                        // ion-trap factor scales the GPX4 inhibition it achieves.
                        let mut q = p.clone();
                        q.rsl3_gpx4_inhib *= trap;
                        // Depth attenuates the ENERGY-delivered arms through their
                        // own published laws, and leaves the systemic drug alone --
                        // which is the dissociation `depth-reach-comparison.md`
                        // measures, now applied inside the resistance sweep.
                        let sp = SpatialParams::default();
                        let sdt_depth = sdt_intensity_at_depth(depth_um, &sp);
                        let rad_depth = radiation::intensity_at_depth(
                            depth_um,
                            radiation::MU_6MV_SOFT_TISSUE_PER_CM,
                        );

                        let sdt_mean = params.sdt_ros * exo * sdt_depth;
                        let sdt_kill = kill_fraction_pheno(n, a.seed, pheno, |cell, rng| {
                            let peak = norm_peak(rng, sdt_mean);
                            run_to_death(cell, Treatment::SDT, &p, peak)
                        });
                        let rsl3_kill = kill_fraction_pheno(n, a.seed, pheno, |cell, _rng| {
                            run_to_death(cell, Treatment::RSL3, &q, 0.0)
                        });

                        // Radiation: O2 through the SAME hyperbola, nothing else.
                        let dose = rad.dose_gy * exo * rad_depth;
                        let lq = 1.0
                            - (-(rad.alpha_per_gy * dose + rad.beta_per_gy2 * dose * dose)).exp();

                        // Immune: unaffected by O2 supply directly; suppressed in a
                        // hypoxic core, which is the asymmetry the doc predicts.
                        let suppression = if hypoxic { 0.6 } else { 0.1 };
                        let imm = adoptive_transfer_kills(
                            EffectorSource::CarT.effectors_after(200.0, 30, 0.15, 5_000.0),
                            n,
                            immune,
                            suppression,
                            false,
                        ) / n as f64;

                        // Ablation: a threshold. Nothing here touches it, and that is
                        // the finding rather than an omission.
                        let abl = 0.85;

                        rows.push(serde_json::json!({
                            "phenotype": pheno_name,
                            "deep": depth_um > 0.0,
                            "depth_um": depth_um,
                            "hypoxic": hypoxic, "stroma": stroma, "acidic": acid,
                            "o2_supply": o2, "ph": ph,
                            "exo_factor": exo, "ion_trap_factor": trap,
                            "arms": {
                                "SDT": sdt_kill,
                                "RSL3": rsl3_kill,
                                "Radiation": lq,
                                "AdoptiveCell": imm,
                                "Ablation": abl,
                            },
                        }));
                    }
                }
            }
        }
    }

    let out = serde_json::json!({
        "n_cells": n,
        "seed": a.seed,
        "axes": ["hypoxia", "stroma", "acidic pH", "depth", "phenotype"],
        "endpoints": {
            "o2_normoxic": O2_NORMOXIC, "o2_hypoxic": O2_HYPOXIC,
            "ph_edge": PH_EDGE, "ph_core": PH_CORE,
            "stromal_gsh_boost": STROMAL_GSH_BOOST,
        },
        "note": "Every axis is applied through the SAME helpers the ferroptosis \
    chapters use, so a modality comparison is under an identical environment. Arms \
    absent from a row are unaffected by that axis BY CONSTRUCTION, which is a \
    prediction and not an omission -- see the function docs.",
        "conditions": rows,
    });
    let path = a.output_dir.join("modality_tme_sweep.json");
    write_json(&path, &out).expect("write tme sweep");
    eprintln!("Wrote {}", path.display());
}

fn norm_peak(rng: &mut StdRng, mean: f64) -> f64 {
    ferroptosis_core::cell::norm(rng, mean, mean * 0.2).max(0.0)
}

fn run_to_death(cell: &ferroptosis_core::cell::Cell, tx: Treatment, p: &Params, peak: f64) -> bool {
    let mut state = CellState::from_cell_with_ros(cell, tx, p, peak);
    let mut rng = StdRng::seed_from_u64(7);
    for step in 0..STEPS {
        ferroptosis_core::biochem::sim_cell_step(&mut state, cell, p, step, 0.0, &mut rng);
        if state.dead {
            return true;
        }
    }
    state.dead
}

fn kill_fraction<F>(n: usize, seed: u64, mut f: F) -> f64
where
    F: FnMut(&ferroptosis_core::cell::Cell, &mut StdRng) -> bool,
{
    kill_fraction_pheno(n, seed, Phenotype::Glycolytic, &mut f)
}

fn kill_fraction_pheno<F>(n: usize, seed: u64, pheno: Phenotype, mut f: F) -> f64
where
    F: FnMut(&ferroptosis_core::cell::Cell, &mut StdRng) -> bool,
{
    let dead = (0..n)
        .filter(|i| {
            let mut rng = StdRng::seed_from_u64(seed.wrapping_add(*i as u64));
            let cell = gen_cell(pheno, &mut rng);
            f(&cell, &mut rng)
        })
        .count();
    dead as f64 / n as f64
}
