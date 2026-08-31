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
use ferroptosis_core::adoptive::{barrier_limited_kills, effective_effectors, AdoptiveBarriers};
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
use ferroptosis_core::physics::{pdt_intensity_at_depth, sdt_intensity_at_depth};
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

    // The barriers are what make this arm different from every other one: the
    // SAME construct cures a blood cancer and does very little in a solid
    // tumour, so the arm is run twice and the pair is the result. The panel
    // row keeps the leukaemia setting (`Default`, every barrier open), which
    // is exactly what it computed before this module existed -- the arm's
    // number is unmoved and the counterfactual is new.
    let infused = EffectorSource::CarT.effectors_after(200.0, 30, 0.15, 5_000.0);
    let leukaemia = AdoptiveBarriers::default();
    let solid = AdoptiveBarriers::solid_tumour();
    let cart_kills = |b: &AdoptiveBarriers| {
        let arrived = effective_effectors(infused, b, STEPS as u32);
        barrier_limited_kills(
            adoptive_transfer_kills(arrived, n, &immune, 0.0, false),
            n as f64,
            b,
        )
    };
    let effectors = infused;
    let cart = cart_kills(&leukaemia);
    let cart_solid = cart_kills(&solid);
    arms.push(ArmResult {
        name: "AdoptiveCell",
        kill_fraction: cart / n as f64,
        route: "redirected effectors (bypasses DC priming)",
        // NOT persistence and not suppression, both of which are exactly
        // inert in this row: the panel passes suppression 0.0 and the default
        // barriers make `persistence_factor` bit-identical to 1.0. Naming
        // them made the row contradict the paragraph below it, which says
        // this is the every-barrier-open case. What actually bounds it is the
        // per-cell kill rate and the PD-1 brake; the barriers that DO bite
        // are in the solid-tumour counterfactual beside it.
        limited_by: "per-effector kill rate and the PD-1 brake (barriers open here; \
        the solid-tumour case is reported separately)",
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
        "adoptive_barriers": {
            "why": "The same infusion, run through the three barriers PMID 31848460 \
    names and then through the antigen ceiling. The panel row above is the leukaemia \
    setting; this is what the identical construct does against a solid tumour. Every \
    barrier VALUE is an uncalibrated placeholder -- the corpus establishes that the \
    barriers are general rather than antigen-specific, not which of them dominates. The \
    antigen ceiling contributes EXACTLY 1x here: delivery and persistence take the kill \
    far below it, so the cap never binds and none of the collapse is attributable to it.",
            "leukaemia_kill_fraction": cart / n as f64,
            "solid_tumour_kill_fraction": cart_solid / n as f64,
            "delivery_efficiency_solid": ferroptosis_core::adoptive::delivery_efficiency(&solid),
            "persistence_at_run_end_solid":
                ferroptosis_core::adoptive::persistence_factor(&solid, STEPS as u32),
            "antigen_ceiling_solid": ferroptosis_core::adoptive::antigen_ceiling(&solid),
        },
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
    // CLONAL HETEROGENEITY as a sixth axis. `clonal.rs` perturbs per-subclone
    // parameters; here the same effect is applied as a spread on the whole
    // population's antioxidant setpoint, because what the axis tests is
    // whether an arm's kill survives VARIANCE and not whether the variance is
    // spatially arranged. An arm that only works on the mean cell is fragile
    // in a way the mean cannot show.
    for &(clonal_spread, heterogeneous) in &[(0.0_f64, false), (0.35, true)] {
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
                            // Heterogeneity widens the antioxidant setpoint. The MEAN is
                            // held fixed so the axis is VARIANCE and not a dose change --
                            // otherwise it would be a second stroma axis under another
                            // name. What it tests is whether an arm's kill survives
                            // spread: an arm that only works on the mean cell is fragile
                            // in a way the mean cannot show.
                            let het = clonal_spread;
                            let mut sdt_lps: Vec<f64> = Vec::new();
                            let mut rsl3_lps: Vec<f64> = Vec::new();
                            let sdt_kill = kill_fraction_pheno(n, a.seed, pheno, |cell, rng| {
                                let mut c = cell.clone();
                                if het > 0.0 {
                                    c.gsh = (c.gsh * (1.0 + norm_unit(rng) * het)).max(0.1);
                                }
                                let peak = norm_peak(rng, sdt_mean);
                                let (dead, lp) = run_to_death_lp(&c, Treatment::SDT, &p, peak);
                                if dead {
                                    sdt_lps.push(lp);
                                }
                                dead
                            });
                            let rsl3_kill = kill_fraction_pheno(n, a.seed, pheno, |cell, rng| {
                                let mut c = cell.clone();
                                if het > 0.0 {
                                    c.gsh = (c.gsh * (1.0 + norm_unit(rng) * het)).max(0.1);
                                }
                                let (dead, lp) = run_to_death_lp(&c, Treatment::RSL3, &q, 0.0);
                                if dead {
                                    rsl3_lps.push(lp);
                                }
                                dead
                            });
                            // PDT: the same exogenous-ROS chemistry as SDT with its OWN
                            // depth law, which is the whole difference between them and
                            // the reason the depth axis has to be in this sweep rather
                            // than only in the panel.
                            let pdt_depth = pdt_intensity_at_depth(depth_um, &sp);
                            let pdt_mean = params.pdt_ros * exo * pdt_depth;
                            let mut pdt_lps: Vec<f64> = Vec::new();
                            let pdt_kill = kill_fraction_pheno(n, a.seed, pheno, |cell, rng| {
                                let mut c = cell.clone();
                                if het > 0.0 {
                                    c.gsh = (c.gsh * (1.0 + norm_unit(rng) * het)).max(0.1);
                                }
                                let peak = norm_peak(rng, pdt_mean);
                                let (dead, lp) = run_to_death_lp(&c, Treatment::PDT, &p, peak);
                                if dead {
                                    pdt_lps.push(lp);
                                }
                                dead
                            });

                            // ADC: the ferroptosis payload again, but reaching the cell
                            // through antibody transport. It is the arm that should be
                            // hit by BOTH depth and pH -- the first because a 150 kDa
                            // carrier barely penetrates, the second because the payload
                            // is a weak base like the free drug.
                            let adc_tissue = epithelial_well_vascularized();
                            let adc_profile = antibody_drug_conjugate();
                            let adc_avail = concentration_at_distance(
                                (depth_um).min(adc_tissue.inter_vessel_distance_um / 2.0),
                                &adc_profile,
                                &adc_tissue,
                            ) * trap;
                            let adc_mean = params.sdt_ros * adc_avail;
                            let mut adc_lps: Vec<f64> = Vec::new();
                            let adc_kill = kill_fraction_pheno(n, a.seed, pheno, |cell, rng| {
                                let mut c = cell.clone();
                                if het > 0.0 {
                                    c.gsh = (c.gsh * (1.0 + norm_unit(rng) * het)).max(0.1);
                                }
                                let peak = norm_peak(rng, adc_mean);
                                let (dead, lp) =
                                    run_to_death_lp(&c, Treatment::AntibodyDrugConjugate, &p, peak);
                                if dead {
                                    adc_lps.push(lp);
                                }
                                dead
                            });

                            // Hypoxia does not enter any immune arm's KILL term. What
                            // it changes is the suppressor field they meet in a
                            // hypoxic core -- a coupling the model asserts rather
                            // than fits, and the report labels it a prediction.
                            let suppression = if hypoxic { 0.6 } else { 0.1 };

                            // Checkpoint blockade: no oxygen term in its kill at all.
                            // What hypoxia does to it is raise the suppressor field,
                            // which is a PREDICTION of the model rather than a fitted
                            // coupling, and the report says so.
                            let blockade_immune = ImmuneParams {
                                baseline_antigenicity: a.baseline_antigenicity,
                                tcell_kill_rate: 0.02,
                                ..ImmuneParams::default()
                            };
                            let blockade = immune_cascade(&[], n, &blockade_immune, true)
                                .immune_kills
                                / n as f64
                                * (1.0 - suppression);

                            // Oncolytic virus: direct lysis plus the SHARED ICD chain.
                            // Neither term touches oxygen, so its only exposure here is
                            // the same suppressor field the other immune arms meet.
                            let (lysed, vq) = oncolytic_lysis(n, 0.15, 0.9, 0.8);
                            let viral_lps = vec![vq * 10.0; lysed as usize];
                            let viral = (lysed
                                + immune_cascade(&viral_lps, n, &blockade_immune, false)
                                    .immune_kills
                                    * (1.0 - suppression))
                                / n as f64;

                            let sdt_immune = immune_amplification(&sdt_lps, n, immune);
                            let rsl3_immune = immune_amplification(&rsl3_lps, n, immune);
                            let sdt_quality = damp_per_death(&sdt_lps, immune);
                            let rsl3_quality = damp_per_death(&rsl3_lps, immune);

                            // Radiation: O2 through the SAME hyperbola, nothing else.
                            let dose = rad.dose_gy * exo * rad_depth;
                            let lq = 1.0
                                - (-(rad.alpha_per_gy * dose + rad.beta_per_gy2 * dose * dose))
                                    .exp();

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
                                "heterogeneous": heterogeneous,
                                "clonal_spread": clonal_spread,
                                "immune_amplification": {
                                    "SDT": sdt_immune,
                                    "RSL3": rsl3_immune,
                                },
                                "damp_per_death": {
                                    "SDT": sdt_quality,
                                    "PDT": damp_per_death(&pdt_lps, immune),
                                    "RSL3": rsl3_quality,
                                    "AntibodyDrugConjugate": damp_per_death(&adc_lps, immune),
                                },
                                "deep": depth_um > 0.0,
                                "depth_um": depth_um,
                                "hypoxic": hypoxic, "stroma": stroma, "acidic": acid,
                                "o2_supply": o2, "ph": ph,
                                "exo_factor": exo, "ion_trap_factor": trap,
                                "arms": {
                                    "SDT": sdt_kill,
                                    "PDT": pdt_kill,
                                    "RSL3": rsl3_kill,
                                    "AntibodyDrugConjugate": adc_kill,
                                    "Radiation": lq,
                                    "AdoptiveCell": imm,
                                    "Immunotherapy": blockade,
                                    "OncolyticVirus": viral,
                                    "Ablation": abl,
                                },
                            }));
                        }
                    }
                }
            }
        }
    }

    let out = serde_json::json!({
        "n_cells": n,
        "seed": a.seed,
        "axes": ["hypoxia", "stroma", "acidic pH", "depth", "clonal heterogeneity"],
        "strata": ["phenotype"],
        "amplification": "immune coupling, reported per arm as the immune \
    kills its OWN death mode earns -- the LP each arm's dead cells carried, run \
    through the same cascade. An arm that kills quietly earns less than one that \
    kills loudly, which is the amplification the ferroptosis chapters measure.",
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
    run_to_death_lp(cell, tx, p, peak).0
}

/// Death, AND the lipid peroxidation the cell carried when it died.
///
/// The second value is what the immune-coupling axis needs, and it is the
/// quantity the ferroptosis chapters already argue about: a cell killed by a
/// runaway ROS cascade dies with far more LP than one killed by slow GPX4
/// inhibition, so it releases more DAMPs per death. "Quality of death" is not
/// a metaphor here -- it is this number.
fn run_to_death_lp(
    cell: &ferroptosis_core::cell::Cell,
    tx: Treatment,
    p: &Params,
    peak: f64,
) -> (bool, f64) {
    let mut state = CellState::from_cell_with_ros(cell, tx, p, peak);
    let mut rng = StdRng::seed_from_u64(7);
    let mut death_step: Option<u32> = None;
    for step in 0..STEPS {
        ferroptosis_core::biochem::sim_cell_step(&mut state, cell, p, step, 0.0, &mut rng);
        if state.dead && death_step.is_none() {
            death_step = Some(step);
        }
        // KEEP STEPPING THROUGH THE GRACE PERIOD, which is where the quality
        // difference actually lives and where the first version of this
        // measurement went wrong.
        //
        // Stopping at the threshold crossing returns LP ~= the death
        // threshold FOR EVERY ARM, by construction -- death IS the crossing.
        // Measured that way, a runaway ROS cascade and a slow GPX4 inhibition
        // both reported ~10.2, and the "quality of death" difference the
        // immune module documents vanished. It had not vanished; it was being
        // measured before it happens.
        //
        // The spatial binaries read `lp_at_grace_end` rather than `lp` for
        // exactly this reason, and the field name says so. `post_death_steps`
        // is how long the cascade keeps running after the cell is counted
        // dead, and an arm delivering exogenous ROS keeps climbing through it
        // while one that merely disabled a repair enzyme does not.
        if let Some(ds) = death_step {
            if step >= ds + p.post_death_steps {
                return (true, state.lp);
            }
        }
    }
    (state.dead, state.lp)
}

/// Immune kills the arm's OWN death mode earns, per cell treated.
///
/// The amplification axis, and the one that separates arms which kill the
/// same number of cells: it runs the LP values the arm actually produced
/// through `immune_cascade`, so an arm that kills quietly earns less than one
/// that kills loudly. `Ablation` and the DNA channel earn nothing here by
/// construction -- they do not produce lipid peroxidation at all -- and that
/// zero is a property of the death mode rather than a gap.
fn immune_amplification(lps: &[f64], n: usize, immune: &ImmuneParams) -> f64 {
    if lps.is_empty() {
        return 0.0;
    }
    immune_cascade(lps, n, immune, false).immune_kills / n as f64
}

/// DAMP release PER DEATH -- the amplification measure that is not confounded
/// by how many cells died.
///
/// The per-treated-cell figure above is confounded in BOTH directions and
/// measurably so: an arm that kills nothing releases no DAMPs and scores zero,
/// while an arm that kills everything leaves no survivors for the immune
/// system to kill and ALSO scores zero, because `immune_cascade` caps its
/// kills at the remaining population. In this panel RSL3 hits the first case
/// at the glycolytic state and SDT hits the second at the persister state, so
/// the two arms are never both readable in the same condition.
///
/// Per-death release has neither problem. It is the "quality of death"
/// quantity the ferroptosis chapters actually argue about -- a cell killed by
/// a runaway ROS cascade dies carrying far more lipid peroxidation than one
/// killed by slow GPX4 inhibition -- and it is comparable across arms whatever
/// their kill fraction.
fn damp_per_death(lps: &[f64], immune: &ImmuneParams) -> f64 {
    if lps.is_empty() {
        return 0.0;
    }
    ferroptosis_core::immune::calculate_damp_release(lps, immune).1
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

/// A standard normal draw, for the heterogeneity spread.
fn norm_unit(rng: &mut StdRng) -> f64 {
    ferroptosis_core::cell::norm(rng, 0.0, 1.0)
}
