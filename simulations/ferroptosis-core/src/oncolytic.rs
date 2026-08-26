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
}
