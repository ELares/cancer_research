//! Oncolytic virus SPREAD, as opposed to a spread that is assumed.
//!
//! # Why this module exists, and it is a gap this project named itself
//!
//! `immune::oncolytic_lysis` models what happens once a virus has infected a
//! given fraction of a tumour, and its own documentation records the
//! limitation plainly: *"`infected_fraction` is an INPUT, not a result, which
//! means this cannot answer 'will the virus spread?', only 'if it spreads this
//! far, what follows?'"*. That is the difference between an arm the engine can
//! express and one it can be asked a question about, and closing it is what
//! this module is for.
//!
//! The three mechanisms that same comment names as absent are the three
//! modelled here, and each is anchored in this project's frozen corpus:
//!
//! * **Replication and spread.** A virus does not arrive everywhere at once;
//!   it infects, replicates, and infects neighbours. Adenoviral replication is
//!   gated on genes the tumour supplies and normal cells do not (PMID
//!   22509396: E1B-55K "is required for viral replication in normal cells but
//!   is dispensable in cancer cells"), which is the selectivity the whole
//!   modality rests on.
//! * **Interferon permissiveness.** Type I IFN is the reason some tumours are
//!   permissive and others are not (PMID 27119111: "inhibition of different
//!   components of the IFN response has previously been shown to increase
//!   virus replication as well as virus yield"; PMID 31694943 on VSV-M
//!   inhibiting IFN "to allow viral replication").
//! * **Antiviral clearance.** The host removes the agent, which is the
//!   ceiling every oncolytic trial runs into (PMID 27119111 again, on HDAC
//!   inhibitors undermining antiviral immunity "resulting in enhanced OV
//!   replication").
//!
//! # The shape, and why it is a logistic race rather than a growth curve
//!
//! Spread and clearance act on the same population at the same time, so the
//! outcome is a RACE and not a trajectory. That matters because the two
//! regimes are qualitatively different and a model that only produced growth
//! would misrepresent the modality: below a threshold the infection dies out
//! whatever the dose, and above it the infection saturates whatever the dose.
//! [`spread_threshold_ratio`] is where that boundary sits, and it is the
//! quantity a reader should take from this module -- not a percentage.
//!
//! # What is still NOT modelled
//!
//! Spatial structure. This is a well-mixed model, so it cannot represent a
//! virus that spreads through a tumour's rim and never reaches its core --
//! which is a real failure mode and the one `analysis/depth-reach-comparison.md`
//! would predict for an agent that has to travel cell-to-cell. A spatial
//! oncolytic model is a separate piece of work and this module does not
//! pretend to it.

use serde::{Deserialize, Serialize};

/// Oncolytic spread parameters.
///
/// [`Default`] is the identity: zero replication, so an unconfigured run has
/// no infection at all and every helper returns zero.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct OncolyticConfig {
    /// Fraction of the tumour infected at delivery.
    #[serde(default)]
    pub initial_infected: f64,
    /// Per-step replication rate -- new infections per infected cell.
    #[serde(default)]
    pub replication_rate: f64,
    /// Per-step antiviral clearance of infected cells.
    #[serde(default)]
    pub clearance_rate: f64,
    /// Tumour interferon competence in `[0, 1]`. `0.0` is a fully permissive
    /// tumour (IFN response lost, which is the state that makes a tumour a
    /// good oncolytic target); `1.0` is fully competent and suppresses
    /// replication entirely.
    #[serde(default)]
    pub interferon_competence: f64,
    /// Fraction of infected cells that lyse per step.
    #[serde(default)]
    pub lysis_rate: f64,
}

impl Default for OncolyticConfig {
    fn default() -> Self {
        Self {
            initial_infected: 0.0,
            replication_rate: 0.0,
            clearance_rate: 0.0,
            interferon_competence: 0.0,
            lysis_rate: 0.0,
        }
    }
}

/// Replication rate after the tumour's interferon response.
///
/// `rate · (1 − competence)`. Linear rather than saturating, deliberately:
/// the corpus records that inhibiting IFN components INCREASES viral yield,
/// which constrains the direction and the endpoints and says nothing about the
/// shape. A hyperbola would look more sophisticated and assert a curvature
/// nobody measured -- the same error `analysis/oxygen-form-check.md` retracted
/// for the oxygen term, in the opposite direction.
#[must_use = "the effective rate is the function's only output"]
pub fn effective_replication(cfg: &OncolyticConfig) -> f64 {
    cfg.replication_rate.max(0.0) * (1.0 - cfg.interferon_competence.clamp(0.0, 1.0))
}

/// The ratio that decides whether an infection spreads at all.
///
/// `effective_replication / clearance`. Above `1.0` the infection grows, below
/// it the infection dies out, and **the initial dose does not change which
/// side of the line a tumour is on** -- it only changes how fast the outcome
/// arrives. That is the modality's defining property and the reason
/// oncolytic trials care so much about tumour permissiveness rather than
/// about titre.
///
/// Returns `f64::INFINITY` when nothing clears the virus, which is the honest
/// answer rather than a clamp: with no clearance any replication wins.
#[must_use = "the ratio is the function's only output"]
pub fn spread_threshold_ratio(cfg: &OncolyticConfig) -> f64 {
    let clear = cfg.clearance_rate.max(0.0);
    if clear == 0.0 {
        return if effective_replication(cfg) > 0.0 {
            f64::INFINITY
        } else {
            0.0
        };
    }
    effective_replication(cfg) / clear
}

/// Infected fraction after `steps`, and the cumulative fraction lysed.
///
/// A logistic race: infection spreads into the uninfected remainder while
/// clearance removes infected cells and lysis converts them to dead ones.
/// Returns `(infected, cumulative_lysed)`, both as fractions of the tumour.
///
/// The cumulative lysed fraction is what `immune::oncolytic_lysis` previously
/// took as an input, so a consumer now DERIVES it instead of assuming it.
#[must_use = "both outputs carry the modality's state"]
pub fn simulate_spread(cfg: &OncolyticConfig, steps: u32) -> (f64, f64) {
    let r = effective_replication(cfg);
    let clear = cfg.clearance_rate.max(0.0);
    let lysis = cfg.lysis_rate.clamp(0.0, 1.0);
    let mut infected = cfg.initial_infected.clamp(0.0, 1.0);
    let mut lysed = 0.0_f64;
    for _ in 0..steps {
        let susceptible = (1.0 - infected - lysed).max(0.0);
        let new = r * infected * susceptible;
        let cleared = clear * infected;
        let died = lysis * infected;
        infected = (infected + new - cleared - died).clamp(0.0, 1.0);
        lysed = (lysed + died).clamp(0.0, 1.0);
        if infected <= 0.0 {
            break;
        }
    }
    (infected, lysed)
}

// ── The double bind, and the spread it is a bind about ───────────────────
//
// Everything above treats an infection as a race between replication and
// clearance, resolved into a threshold ratio. That is the right first
// question and it leaves out the one that makes this modality strange.
//
// An oncolytic virus has TWO mechanisms and they want opposite things from the
// immune system. The direct one is lysis, and the immune response clears the
// virus before it can spread -- so for lysis, immunity is the enemy. The
// durable one is the anti-TUMOUR response the lysis primes, and that needs an
// immune system to exist at all. A model with one immune term can express
// neither tension nor the consequence: that there is an interior optimum, and
// that suppressing immunity to help the virus spread can make the outcome
// worse.
//
// The spread itself is a travelling front, and the engine already has the
// machinery for one: `trigger_wave` solves a bistable reaction-diffusion front
// for ferroptosis. An infection is the classic Fisher case rather than the
// bistable one -- there is no ignition threshold, any inoculum grows if R0
// exceeds one -- so the speed is the Fisher-KPP result rather than the Nagumo
// one, and the two live side by side deliberately.

/// Entry-receptor density below which a cell cannot be infected.
///
/// The oncolytic analogue of the CAR-T antigen threshold in
/// [`crate::adoptive`], and it behaves the same way: a cell below it is not
/// infected slowly, it is not infected, so it bounds what any titre can reach.
/// PLACEHOLDER; entry-receptor requirements are virus-specific and nothing
/// here fits one.
pub const ENTRY_RECEPTOR_THRESHOLD: f64 = 100.0;

/// Basic reproduction number: secondary infections per infected cell.
///
/// Below 1 an inoculum dies out whatever its size, which is the same
/// statement [`spread_threshold_ratio`] makes and is repeated here in the
/// units the front-speed result needs.
#[must_use]
pub fn basic_reproduction_number(cfg: &OncolyticConfig) -> f64 {
    let cleared = cfg.clearance_rate.max(f64::MIN_POSITIVE);
    (effective_replication(cfg) / cleared).max(0.0)
}

/// Speed of the infection front, in cell diameters per step.
///
/// The Fisher-KPP result `c = 2 sqrt(D r)`, where `r` is the net growth rate
/// of the infected population and `D` its spatial spread. Below `R0 = 1` the
/// net rate is negative and there is no front: this returns 0 rather than an
/// imaginary speed, which is the honest form of "the infection does not
/// establish".
///
/// **Deliberately a different equation from [`crate::trigger_wave`].** That
/// front is bistable and has an ignition threshold; this one does not, because
/// any inoculum grows when R0 exceeds one. Using one solver for both would be
/// tidier and would misrepresent one of them.
#[must_use]
pub fn front_speed(cfg: &OncolyticConfig, diffusion: f64) -> f64 {
    let r = effective_replication(cfg) - cfg.clearance_rate;
    if r <= 0.0 || diffusion <= 0.0 {
        return 0.0;
    }
    2.0 * (diffusion * r).sqrt()
}

/// Fraction of a tumour a virus can enter at all, given receptor density.
///
/// Same shape as [`crate::adoptive::density_engagement`] and for the same
/// reason: a threshold is not a barrier, and a cell without the receptor is
/// not reachable by a larger dose.
#[must_use]
pub fn permissive_fraction(receptor_density: f64, threshold: f64, steepness: f64) -> f64 {
    let d = receptor_density.max(0.0);
    let t = threshold.max(f64::MIN_POSITIVE);
    let x = (d / t).powf(steepness.max(f64::MIN_POSITIVE));
    (x / (1.0 + x)).clamp(0.0, 1.0)
}

/// How much of the virus survives an immune response of a given competence.
///
/// Monotonically DOWN: a more competent immune system clears the virus faster,
/// so less of it is left to spread. This is the arm of the bind that wants
/// immunosuppression.
#[must_use]
pub fn virus_survival(immune_competence: f64, clearance_sensitivity: f64) -> f64 {
    let c = immune_competence.max(0.0);
    (1.0 / (1.0 + clearance_sensitivity.max(0.0) * c)).clamp(0.0, 1.0)
}

/// Anti-tumour immunity primed by lysis, given the same immune competence.
///
/// Monotonically UP, and multiplied by the lysis that did the priming: an
/// immune system that cannot respond cannot be taught, and a virus that
/// cannot lyse has nothing to teach with. This is the arm of the bind that
/// wants immunocompetence.
#[must_use]
pub fn antitumour_priming(
    immune_competence: f64,
    lysed_fraction: f64,
    priming_efficiency: f64,
) -> f64 {
    let c = immune_competence.max(0.0);
    let l = lysed_fraction.clamp(0.0, 1.0);
    (priming_efficiency.max(0.0) * l * c / (1.0 + c)).clamp(0.0, 1.0)
}

/// The durable outcome: lysis that happened, plus what the immunity it primed
/// kills afterwards.
///
/// The second term acts on the tumour the virus did NOT reach, which is the
/// only way it can add anything: a primed response that killed only the cells
/// already lysed would be arithmetic rather than biology, and the sum would be
/// a multiple of the lysis and therefore monotone in immune competence.
///
/// **Whether the sum has an interior optimum is NOT built in.** It depends on
/// whether the primed response repays the lysis it costs, and
/// [`priming_efficiency_for_interior_optimum`] finds the efficiency at which
/// it starts to. That threshold is the layer's actual finding: below it,
/// suppressing immunity to help the virus spread is the right move and the
/// model says so; above it, the same move destroys the durable arm. A model
/// that produced an optimum at every parameter would be asserting the
/// conclusion rather than deriving it.
#[must_use]
pub fn durable_outcome(
    cfg: &OncolyticConfig,
    steps: u32,
    immune_competence: f64,
    clearance_sensitivity: f64,
    priming_efficiency: f64,
    permissive: f64,
) -> f64 {
    let survival = virus_survival(immune_competence, clearance_sensitivity);
    // The SECOND element is the cumulative lysed fraction; the first is the
    // population still infected and therefore still alive. `sim-modality-panel`
    // carries a comment about having made exactly this mistake once, and this
    // function made it again -- which is what a comment on the other end of a
    // tuple cannot prevent.
    let (_still_infected, lysed_frac) = simulate_spread(cfg, steps);
    let lysed = (lysed_frac * survival * permissive.clamp(0.0, 1.0)).clamp(0.0, 1.0);
    let primed = antitumour_priming(immune_competence, lysed, priming_efficiency);
    // Acting on what the virus did not reach.
    (lysed + primed * (1.0 - lysed)).clamp(0.0, 1.0)
}

/// The priming efficiency at which the optimum stops being at zero immune
/// competence, or `None` if it never does within the range scanned.
///
/// **The layer's finding, and it is a CONDITION rather than a number.** The
/// double bind is real -- immunity clears the virus and immunity is the
/// durable mechanism -- but which side wins is not settled by naming the
/// tension. It is settled by whether the primed response kills more than the
/// lysis it cost, and this returns the efficiency where that crossover
/// happens. Below it, the clinical intuition (suppress immunity, let the virus
/// work) is what this model recommends. Above it, that move throws away the
/// only durable arm.
#[must_use]
pub fn priming_efficiency_for_interior_optimum(
    cfg: &OncolyticConfig,
    steps: u32,
    clearance_sensitivity: f64,
    permissive: f64,
    max_competence: f64,
    max_efficiency: f64,
) -> Option<f64> {
    let n = 100;
    for i in 0..=n {
        let e = max_efficiency * f64::from(i) / f64::from(n);
        let (best_c, _) = optimal_immune_competence(
            cfg,
            steps,
            clearance_sensitivity,
            e,
            permissive,
            max_competence,
        );
        // "Interior" means the optimum is not at the suppressed end. A
        // tolerance rather than a strict inequality, because the scan is on a
        // grid and the first grid point is not zero competence in any
        // meaningful sense.
        if best_c > max_competence * 0.05 {
            return Some(e);
        }
    }
    None
}

/// The immune competence that maximises the durable outcome, by scan.
///
/// Returns `(competence, outcome)`. A scan rather than a derivative because
/// `simulate_spread` is a stepper, not a closed form -- and because the shape
/// is what is being reported, so producing it the same way it is drawn keeps
/// the figure and the claim on one computation.
#[must_use]
pub fn optimal_immune_competence(
    cfg: &OncolyticConfig,
    steps: u32,
    clearance_sensitivity: f64,
    priming_efficiency: f64,
    permissive: f64,
    max_competence: f64,
) -> (f64, f64) {
    let mut best = (0.0, f64::NEG_INFINITY);
    let n = 200;
    for i in 0..=n {
        let c = max_competence * f64::from(i) / f64::from(n);
        let v = durable_outcome(
            cfg,
            steps,
            c,
            clearance_sensitivity,
            priming_efficiency,
            permissive,
        );
        if v > best.1 {
            best = (c, v);
        }
    }
    best
}

#[cfg(test)]
mod tests {
    use super::*;

    fn permissive() -> OncolyticConfig {
        OncolyticConfig {
            initial_infected: 0.01,
            replication_rate: 0.35,
            clearance_rate: 0.05,
            interferon_competence: 0.0,
            lysis_rate: 0.05,
        }
    }

    #[test]
    fn the_default_config_never_infects_anything() {
        let d = OncolyticConfig::default();
        assert_eq!(effective_replication(&d).to_bits(), 0.0_f64.to_bits());
        let (inf, lysed) = simulate_spread(&d, 200);
        assert_eq!(inf.to_bits(), 0.0_f64.to_bits());
        assert_eq!(lysed.to_bits(), 0.0_f64.to_bits());
        assert_eq!(spread_threshold_ratio(&d).to_bits(), 0.0_f64.to_bits());
    }

    /// THE defining property: the threshold decides the outcome, and the DOSE
    /// does not move it.
    #[test]
    fn the_dose_does_not_decide_whether_the_infection_spreads() {
        // Sub-threshold: clearance outruns replication.
        let dying = OncolyticConfig {
            replication_rate: 0.05,
            clearance_rate: 0.30,
            ..permissive()
        };
        assert!(spread_threshold_ratio(&dying) < 1.0);
        // A HUNDREDFOLD larger dose does not rescue it.
        for &dose in &[0.001_f64, 0.01, 0.1] {
            let c = OncolyticConfig {
                initial_infected: dose,
                ..dying
            };
            let (inf, lysed) = simulate_spread(&c, 300);
            assert!(
                inf < dose,
                "at dose {dose} the sub-threshold infection did not shrink \
                 ({inf})"
            );
            assert!(
                lysed < 0.5,
                "a sub-threshold infection lysed {lysed} of the tumour; the \
                 threshold is supposed to decide the outcome"
            );
        }
        // Super-threshold: a TINY dose still takes the tumour.
        let growing = permissive();
        assert!(spread_threshold_ratio(&growing) > 1.0);
        let (_, lysed_small) = simulate_spread(
            &OncolyticConfig {
                initial_infected: 1e-4,
                ..growing
            },
            600,
        );
        assert!(
            lysed_small > 0.5,
            "a super-threshold infection from a small dose only reached \
             {lysed_small}"
        );
    }

    /// Interferon competence is what moves a tumour across the line, which is
    /// the clinical selection criterion the modality actually uses.
    #[test]
    fn interferon_competence_moves_the_tumour_across_the_threshold() {
        let permissive_cfg = permissive();
        assert!(spread_threshold_ratio(&permissive_cfg) > 1.0);

        let competent = OncolyticConfig {
            interferon_competence: 0.95,
            ..permissive_cfg
        };
        assert!(
            spread_threshold_ratio(&competent) < 1.0,
            "a fully IFN-competent tumour should be non-permissive, got {}",
            spread_threshold_ratio(&competent)
        );
        let (_, lysed) = simulate_spread(&competent, 600);
        assert!(
            lysed < 0.2,
            "an IFN-competent tumour was still cleared: {lysed}"
        );

        // Monotone: more competence, less spread, and it must be a real
        // gradient rather than a switch at the endpoints.
        let mut prev = 1.0_f64;
        for &ifn in &[0.0_f64, 0.25, 0.5, 0.75, 1.0] {
            let c = OncolyticConfig {
                interferon_competence: ifn,
                ..permissive_cfg
            };
            let (_, l) = simulate_spread(&c, 400);
            assert!(l <= prev + 1e-12, "not monotone in IFN at {ifn}");
            prev = l;
        }
        // Fully competent means no replication at all, exactly.
        let dead = OncolyticConfig {
            interferon_competence: 1.0,
            ..permissive_cfg
        };
        assert_eq!(effective_replication(&dead).to_bits(), 0.0_f64.to_bits());
    }

    #[test]
    fn clearance_and_replication_are_a_race_not_a_trajectory() {
        // The same replication rate gives opposite outcomes depending only on
        // clearance, which is what "race" means and what a growth-curve model
        // could not express.
        let fast_clear = OncolyticConfig {
            clearance_rate: 0.5,
            ..permissive()
        };
        let slow_clear = OncolyticConfig {
            clearance_rate: 0.02,
            ..permissive()
        };
        let (_, lysed_fast) = simulate_spread(&fast_clear, 400);
        let (_, lysed_slow) = simulate_spread(&slow_clear, 400);
        assert!(
            lysed_slow > lysed_fast * 2.0,
            "clearance barely changed the outcome: {lysed_slow} vs \
             {lysed_fast}"
        );
        assert!(spread_threshold_ratio(&fast_clear) < spread_threshold_ratio(&slow_clear));
        // No clearance at all is INFINITY rather than a clamped large number,
        // because "nothing removes it" is a different statement from "a lot
        // of replication".
        let unopposed = OncolyticConfig {
            clearance_rate: 0.0,
            ..permissive()
        };
        assert!(spread_threshold_ratio(&unopposed).is_infinite());
    }

    #[test]
    fn the_population_is_conserved_and_bounded() {
        // Infected + lysed can never exceed the tumour, at any parameters --
        // including ones chosen to break it.
        for &(r, c, l) in &[
            (5.0_f64, 0.0_f64, 0.9_f64),
            (0.9, 0.9, 0.9),
            (10.0, 0.001, 0.5),
        ] {
            let cfg = OncolyticConfig {
                initial_infected: 0.5,
                replication_rate: r,
                clearance_rate: c,
                interferon_competence: 0.0,
                lysis_rate: l,
            };
            let (inf, lysed) = simulate_spread(&cfg, 500);
            assert!((0.0..=1.0).contains(&inf), "infected {inf} at r={r}");
            assert!((0.0..=1.0).contains(&lysed), "lysed {lysed} at r={r}");
            assert!(
                inf + lysed <= 1.0 + 1e-9,
                "infected {inf} + lysed {lysed} exceeds the tumour at r={r}"
            );
        }
    }

    /// The lysed fraction this module DERIVES is the quantity
    /// `immune::oncolytic_lysis` used to take as an input.
    #[test]
    fn the_derived_lysed_fraction_feeds_the_icd_chain() {
        use crate::immune::oncolytic_lysis;
        let (_, lysed) = simulate_spread(&permissive(), 400);
        assert!(lysed > 0.0);
        const N: usize = 10_000;
        let (cells, quality) = oncolytic_lysis(N, lysed, 0.9, 0.8);
        assert!(cells > 0.0 && quality > 0.0);
        // And a tumour the virus cannot take produces NO ICD, which is the
        // link the gap in `immune.rs` could not express: whether the systemic
        // effect happens at all now depends on whether the virus spreads.
        let competent = OncolyticConfig {
            interferon_competence: 1.0,
            ..permissive()
        };
        let (_, none) = simulate_spread(&competent, 400);
        let (no_cells, _) = oncolytic_lysis(N, none, 0.9, 0.8);
        assert!(
            no_cells < cells * 0.1,
            "a non-permissive tumour still produced {no_cells} lysed cells \
             against {cells}"
        );
    }

    /// P12's registered claim, asserted against the ENGINE rather than a
    /// re-implementation of it.
    ///
    /// `scripts/modality_predictions.py` mirrors this loop to derive the
    /// numbers `PREREGISTRATION.md` publishes, and its first version
    /// accumulated `lysed` from the POST-update infected count where this
    /// function uses the pre-update one -- so every published figure was wrong
    /// against the model it claimed to describe.
    ///
    /// WHAT THIS TEST DOES AND DOES NOT ESTABLISH. It asserts the PROPERTY the
    /// prediction registers, and the Python side asserts the same property, so
    /// a change that breaks titre-independence fails on both sides. It does
    /// NOT compare the two implementations' NUMBERS -- nothing in the Python
    /// suite invokes cargo -- so a change that shifts both by the same amount
    /// is invisible here. An earlier version of this comment claimed the two
    /// "cannot silently part", which overstated it.
    #[test]
    fn establishment_is_titre_independent_across_five_orders() {
        let cfg = OncolyticConfig {
            initial_infected: 0.0,
            replication_rate: 0.9,
            clearance_rate: 0.2,
            interferon_competence: 0.3,
            lysis_rate: 0.15,
        };
        let titres = [1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1];
        let lysed: Vec<f64> = titres
            .iter()
            .map(|t| {
                simulate_spread(
                    &OncolyticConfig {
                        initial_infected: *t,
                        ..cfg
                    },
                    180,
                )
                .1
            })
            .collect();
        // The verdict is the same at every titre...
        assert!(
            lysed.iter().all(|l| *l > 0.01),
            "the infection fails to establish at some titre: {lysed:?}"
        );
        // ...and the OUTCOME barely moves across five orders of magnitude,
        // which is the content of the prediction. A titre-dependent
        // establishment would show up here first.
        let spread = lysed.iter().cloned().fold(f64::MIN, f64::max)
            - lysed.iter().cloned().fold(f64::MAX, f64::min);
        assert!(
            spread < 0.01,
            "cumulative lysis spans {spread:.4} across five orders of \
             magnitude of initial titre, so the outcome is titre-DEPENDENT \
             and P12 must be rewritten rather than this bound relaxed"
        );
        // And the criterion that decides it contains no titre term: replication
        // must beat clearance PLUS lysis, which `spread_threshold_ratio`
        // (clearance alone) does not capture.
        let eff = effective_replication(&cfg);
        assert!(eff > cfg.clearance_rate + cfg.lysis_rate);
        assert!(spread_threshold_ratio(&cfg) > eff / (cfg.clearance_rate + cfg.lysis_rate));
    }
    // ── The double bind ──────────────────────────────────────────────────

    #[test]
    fn a_front_needs_r0_above_one_and_reports_zero_otherwise() {
        let spreading = OncolyticConfig {
            initial_infected: 0.01,
            replication_rate: 0.9,
            clearance_rate: 0.2,
            interferon_competence: 0.3,
            lysis_rate: 0.15,
        };
        assert!(basic_reproduction_number(&spreading) > 1.0);
        assert!(front_speed(&spreading, 0.5) > 0.0);
        let dying = OncolyticConfig {
            replication_rate: 0.05,
            clearance_rate: 0.5,
            ..spreading
        };
        assert!(basic_reproduction_number(&dying) < 1.0);
        assert!(
            front_speed(&dying, 0.5).abs() < 1e-12,
            "a sub-threshold infection reported a front speed"
        );
        // And the speed rises with both diffusion and net growth, which is the
        // Fisher form rather than a fitted curve.
        assert!(front_speed(&spreading, 2.0) > front_speed(&spreading, 0.5));
    }

    #[test]
    fn the_front_speed_scales_as_the_square_root_of_diffusion() {
        // A MUTATION SURVIVOR. Asserting the speed rises with diffusion passes
        // for a LINEAR form too, and the Fisher result is `2 sqrt(D r)` -- the
        // square root is the whole physical content, because it is what makes
        // a front slower than a doubling would suggest.
        let cfg = OncolyticConfig {
            initial_infected: 0.01,
            replication_rate: 0.9,
            clearance_rate: 0.2,
            interferon_competence: 0.3,
            lysis_rate: 0.15,
        };
        let one = front_speed(&cfg, 0.5);
        let four = front_speed(&cfg, 2.0);
        assert!((four / one - 2.0).abs() < 1e-9,
                "quadrupling the diffusion should double the speed, not scale it                  by {}", four / one);
    }

    #[test]
    fn the_outcome_reads_the_lysed_fraction_and_not_the_infected_one() {
        // A MUTATION SURVIVOR, and the SAME bug this file already fixed once:
        // `simulate_spread` returns (still infected, cumulative lysed), and
        // reading the first is reading the population that has NOT died.
        // Comparing outcomes to each other cannot catch it, because both
        // elements move together. This pins the identity.
        let cfg = OncolyticConfig {
            initial_infected: 0.01,
            replication_rate: 0.9,
            clearance_rate: 0.2,
            interferon_competence: 0.3,
            lysis_rate: 0.15,
        };
        let (still_infected, lysed) = simulate_spread(&cfg, 60);
        assert!(lysed > still_infected,
                "after sixty steps more should have lysed ({lysed}) than remain                  infected ({still_infected}) -- if not, this test's premise is                  wrong rather than the code");
        // With no immunity, no priming and full permissiveness the outcome IS
        // the cumulative lysed fraction, exactly.
        let outcome = durable_outcome(&cfg, 60, 0.0, 0.5, 0.0, 1.0);
        assert!(
            (outcome - lysed).abs() < 1e-12,
            "the outcome is {outcome} and the lysed fraction is {lysed}"
        );
    }

    #[test]
    fn the_two_arms_of_the_bind_pull_in_opposite_directions() {
        let low = virus_survival(0.5, 0.5);
        let high = virus_survival(8.0, 0.5);
        assert!(
            high < low,
            "a more competent immune system did not clear more virus"
        );
        let primed_low = antitumour_priming(0.5, 0.4, 1.0);
        let primed_high = antitumour_priming(8.0, 0.4, 1.0);
        assert!(
            primed_high > primed_low,
            "priming did not rise with competence"
        );
        // Neither term is non-monotonic on its own; that is the point.
        assert!(
            antitumour_priming(4.0, 0.0, 1.0).abs() < 1e-12,
            "priming happened with no lysis to prime it"
        );
    }

    #[test]
    fn the_interior_optimum_is_a_condition_and_not_a_built_in() {
        // THE LAYER'S FINDING. A model that produced an optimum at every
        // parameter would be asserting the conclusion. Below a priming
        // efficiency the model can name, suppressing immunity IS the right
        // move here; above it, that move throws away the durable arm.
        let cfg = OncolyticConfig {
            initial_infected: 0.01,
            replication_rate: 0.9,
            clearance_rate: 0.2,
            interferon_competence: 0.3,
            lysis_rate: 0.15,
        };
        let (weak_c, _) = optimal_immune_competence(&cfg, 60, 0.5, 0.5, 1.0, 16.0);
        assert!(
            weak_c.abs() < 1e-9,
            "with weak priming the optimum should be full suppression, got {weak_c}"
        );
        let (strong_c, _) = optimal_immune_competence(&cfg, 60, 0.5, 4.0, 1.0, 16.0);
        assert!(
            strong_c > 0.8,
            "with strong priming the optimum should be interior, got {strong_c}"
        );
        let threshold = priming_efficiency_for_interior_optimum(&cfg, 60, 0.5, 1.0, 16.0, 20.0)
            .expect("an interior optimum should appear somewhere in the range");
        assert!((0.5..=10.0).contains(&threshold),
                "the crossover efficiency is {threshold}, outside anything this                  layer could describe");
    }

    #[test]
    fn the_primed_response_must_act_on_what_the_virus_did_not_reach() {
        // If it acted only on cells already lysed, the sum would be a multiple
        // of the lysis and could not be non-monotonic at all -- the finding
        // above would be impossible by construction rather than false.
        let cfg = OncolyticConfig {
            initial_infected: 0.01,
            replication_rate: 0.9,
            clearance_rate: 0.2,
            interferon_competence: 0.3,
            lysis_rate: 0.15,
        };
        let with_priming = durable_outcome(&cfg, 60, 2.0, 0.5, 4.0, 1.0);
        let without = durable_outcome(&cfg, 60, 2.0, 0.5, 0.0, 1.0);
        assert!(
            with_priming > without,
            "priming added nothing: {with_priming} vs {without}"
        );
        assert!(with_priming <= 1.0, "the outcome passed 1");
    }

    #[test]
    fn a_tumour_without_the_entry_receptor_is_not_reachable() {
        let t = ENTRY_RECEPTOR_THRESHOLD;
        assert!((permissive_fraction(t, t, 3.0) - 0.5).abs() < 1e-12);
        assert!(permissive_fraction(t * 0.2, t, 3.0) < 0.02);
        assert!(permissive_fraction(t * 5.0, t, 3.0) > 0.98);
        let cfg = OncolyticConfig {
            initial_infected: 0.01,
            replication_rate: 0.9,
            clearance_rate: 0.2,
            interferon_competence: 0.3,
            lysis_rate: 0.15,
        };
        let reachable = durable_outcome(&cfg, 60, 1.0, 0.5, 1.0, 1.0);
        let not = durable_outcome(&cfg, 60, 1.0, 0.5, 1.0, 0.0);
        assert!(
            not.abs() < 1e-12,
            "a non-permissive tumour was still killed: {not}"
        );
        assert!(reachable > not);
    }
}
