//! Cytotoxic chemotherapy: the modality most patients actually receive.
//!
//! # Why this module exists
//!
//! Until it did, `Treatment` had nine variants and none of them was
//! chemotherapy. This engine could express a laboratory ferroptosis inducer
//! with no clinical use in several thousand lines, and could not express
//! cisplatin at all. That is the sharpest form of the criticism this project
//! has been answering -- not that ferroptosis is over-served, but that the
//! comparator was missing.
//!
//! # The cell cycle, and the scope choice that comes with it
//!
//! Chemotherapy classes are distinguished by WHEN in the cell cycle they act,
//! and the engine had no cell cycle. This module adds one at POPULATION level:
//! a [`PhaseDistribution`] over G1, S, G2/M and quiescent G0, advanced between
//! doses by redistribution. It does NOT add a per-cell cycle position, which
//! would touch every module in the crate and change results the byte-identity
//! gates protect.
//!
//! The consequence is stated rather than hidden. A population-level cycle can
//! express phase specificity, the saturation it produces, and redistribution
//! between doses. It cannot express a cell's own progress through the cycle,
//! so nothing here models a synchronised population or the phase-specific
//! scheduling that would follow from one.
//!
//! # What the layer predicts rather than assumes
//!
//! Two things, and both are free to come out the other way.
//!
//! A phase-specific agent can only kill the cells that are in its phase, so
//! over an achievable dose range it leaves a far larger residue than a
//! phase-nonspecific agent of the same potency. That gap is not written into
//! any constant here; it falls out of multiplying a dose-response by the
//! fraction of the population that is sensitive. See
//! [`phase_specific_residue_ratio`] -- and note what it is NOT: with no phase
//! given exactly zero sensitivity, unlimited dose still kills everything, so
//! the residue is a statement about doses a patient could receive rather than
//! about an asymptote.
//!
//! And the Norton-Simon hypothesis (Norton & Simon 1986, PMID 3510732): if a
//! tumour regrows between cycles, shortening the interval at the same total
//! dose kills more. The model produces the sign, and -- the part that makes it
//! a test -- the advantage VANISHES as regrowth slows, which is the
//! comparison [`dose_density_advantage`] returns.
//!
//! # What is not calibrated, and why
//!
//! No dose-response here is fitted. The repository's CTRPv2 route reaches only
//! the five ferroptosis compounds it fetched before the DepMap download
//! catalogue moved behind a verification page, which
//! `analysis/calibration/calibration-feasibility.md` already records as an
//! ACCESS block rather than a method one. So every absolute kill fraction this
//! module produces is a placeholder, and the DIRECTIONS above are the results.

use serde::{Deserialize, Serialize};

/// Where a cell is in the division cycle.
///
/// `G0` is not a phase of the cycle but the state of being out of it, and it
/// is the one that matters most here: a quiescent cell is invisible to a
/// phase-specific agent, which is the population-level reason chemotherapy
/// leaves a residue.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum CellCyclePhase {
    /// Growth before DNA synthesis.
    G1,
    /// DNA synthesis: the target of the antimetabolites.
    S,
    /// Growth and mitosis: the target of the spindle poisons.
    G2M,
    /// Out of cycle. Quiescent, and hard to kill with anything that needs a
    /// dividing cell.
    G0,
}

/// The share of a population in each phase.
///
/// Population level, deliberately: see the module docs.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct PhaseDistribution {
    /// Fraction in G1.
    pub g1: f64,
    /// Fraction in S.
    pub s: f64,
    /// Fraction in G2/M.
    pub g2m: f64,
    /// Fraction out of cycle.
    pub g0: f64,
}

impl PhaseDistribution {
    /// A proliferating culture: most cells cycling, a minority quiescent.
    ///
    /// The shape of a flow-cytometry profile from an exponentially growing
    /// line, and a PLACEHOLDER in the sense the layer-freeze policy means --
    /// the fractions are conventional rather than fitted, and every absolute
    /// kill that rests on them inherits that.
    #[must_use]
    pub fn proliferating() -> Self {
        Self {
            g1: 0.55,
            s: 0.25,
            g2m: 0.15,
            g0: 0.05,
        }
    }

    /// A tumour whose growth fraction is low: most cells out of cycle.
    ///
    /// The state that matters clinically, because it is what a solid tumour
    /// looks like and it is why the curves measured in a dividing culture
    /// overstate what a drug does in a patient.
    #[must_use]
    pub fn quiescent_rich() -> Self {
        Self {
            g1: 0.25,
            s: 0.08,
            g2m: 0.05,
            g0: 0.62,
        }
    }

    /// The share of cells IN cycle.
    #[must_use]
    pub fn growth_fraction(&self) -> f64 {
        (self.g1 + self.s + self.g2m).clamp(0.0, 1.0)
    }

    /// The share in one phase.
    #[must_use]
    pub fn share(&self, phase: CellCyclePhase) -> f64 {
        match phase {
            CellCyclePhase::G1 => self.g1,
            CellCyclePhase::S => self.s,
            CellCyclePhase::G2M => self.g2m,
            CellCyclePhase::G0 => self.g0,
        }
    }

    /// Sum of the four shares, which a caller may check is 1.
    #[must_use]
    pub fn total(&self) -> f64 {
        self.g1 + self.s + self.g2m + self.g0
    }

    /// Renormalise after cells have been removed from some phases and not
    /// others, which is exactly what a phase-specific agent does.
    #[must_use]
    pub fn normalised(&self) -> Self {
        let t = self.total();
        if t <= 0.0 {
            return *self;
        }
        Self {
            g1: self.g1 / t,
            s: self.s / t,
            g2m: self.g2m / t,
            g0: self.g0 / t,
        }
    }
}

/// How a drug class relates to the cycle.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum ChemoClass {
    /// Alkylators and platinums: damage DNA whatever the cell is doing.
    PhaseNonspecific,
    /// Antimetabolites: act on cells synthesising DNA.
    SPhaseSpecific,
    /// Taxanes and vinca alkaloids: act on the mitotic spindle.
    MPhaseSpecific,
}

/// How sensitive a phase is to a class, as a multiplier on the kill rate.
///
/// The whole content of "phase specificity" is in this function, and the
/// numbers are conventional rather than measured: what is being asserted is
/// the ORDERING and the zeros, not the sizes. A quiescent cell is not immune
/// to an alkylator -- it is less sensitive, because it is not replicating the
/// damage into a lethal lesion -- and the value here says so rather than
/// setting it to zero.
#[must_use]
pub fn phase_sensitivity(class: ChemoClass, phase: CellCyclePhase) -> f64 {
    use CellCyclePhase::{G0, G1, G2M, S};
    use ChemoClass::{MPhaseSpecific, PhaseNonspecific, SPhaseSpecific};
    match (class, phase) {
        (PhaseNonspecific, G0) => 0.4,
        (PhaseNonspecific, _) => 1.0,
        (SPhaseSpecific, S) => 1.0,
        (SPhaseSpecific, G1) => 0.15,
        (SPhaseSpecific, G2M) => 0.1,
        (SPhaseSpecific, G0) => 0.02,
        (MPhaseSpecific, G2M) => 1.0,
        (MPhaseSpecific, G1) => 0.05,
        (MPhaseSpecific, S) => 0.1,
        (MPhaseSpecific, G0) => 0.01,
    }
}

/// Surviving fraction after one exposure, phase by phase.
///
/// Log kill within each phase -- `exp(-potency * dose * sensitivity)`, the
/// constant-fraction form Skipper's work established -- weighted by how much
/// of the population is in that phase. The weighting is what makes a
/// phase-specific agent behave differently from a phase-nonspecific one at the
/// same potency.
#[must_use]
pub fn surviving_fraction(
    dose: f64,
    potency: f64,
    class: ChemoClass,
    dist: &PhaseDistribution,
) -> f64 {
    use CellCyclePhase::{G0, G1, G2M, S};
    if dose <= 0.0 || potency <= 0.0 {
        return 1.0;
    }
    [G1, S, G2M, G0]
        .iter()
        .map(|&p| {
            let killed_rate = potency * dose * phase_sensitivity(class, p);
            dist.share(p) * (-killed_rate).exp()
        })
        .sum::<f64>()
        .clamp(0.0, 1.0)
}

/// The distribution left behind after an exposure.
///
/// A phase-specific agent does not merely kill fewer cells: it changes WHICH
/// cells are left, enriching the population for the phases it cannot reach.
/// That is the mechanism behind a resistant residue that is not genetically
/// resistant at all, and it is the same shape as the drug-tolerant persister
/// state [`crate::persister`] already models by a different route.
#[must_use]
pub fn surviving_distribution(
    dose: f64,
    potency: f64,
    class: ChemoClass,
    dist: &PhaseDistribution,
) -> PhaseDistribution {
    use CellCyclePhase::{G0, G1, G2M, S};
    let survive =
        |p: CellCyclePhase| dist.share(p) * (-(potency * dose * phase_sensitivity(class, p))).exp();
    PhaseDistribution {
        g1: survive(G1),
        s: survive(S),
        g2m: survive(G2M),
        g0: survive(G0),
    }
    .normalised()
}

/// How much more a phase-specific agent leaves behind than a
/// phase-nonspecific one, at the same dose and potency.
///
/// **The layer's first prediction, and it is emergent.** Nothing in
/// [`phase_sensitivity`] says "leave a residue"; the gap appears because a
/// phase-specific agent can only reach the cells in its phase, so raising the
/// dose kills those and leaves the others. Returns the ratio of surviving
/// fractions, which is above 1 whenever the phase-specific agent is the
/// weaker one at that dose.
///
/// **It is a property of the ACHIEVABLE dose range, not an asymptote, and the
/// first version of this function got that wrong.** It evaluated the survival
/// at an enormous dose and reported the limit as the plateau height -- which
/// is zero for every class, because no phase here has exactly zero
/// sensitivity, so given unlimited dose everything dies. That is not a defect
/// in the sensitivity table: a phase with a hard zero would produce a true
/// plateau BY ASSUMPTION, and the interesting statement is the one that
/// survives without it. So the comparison is made at a dose the caller names.
#[must_use]
pub fn phase_specific_residue_ratio(
    dose: f64,
    potency: f64,
    specific: ChemoClass,
    dist: &PhaseDistribution,
) -> f64 {
    let nonspecific = surviving_fraction(dose, potency, ChemoClass::PhaseNonspecific, dist);
    if nonspecific <= 0.0 {
        return f64::INFINITY;
    }
    surviving_fraction(dose, potency, specific, dist) / nonspecific
}

/// Gompertzian regrowth: the burden after `days`, from a starting burden.
///
/// The growth law tumours are conventionally described by, and the reason
/// Norton and Simon expected dose density to matter: growth RATE falls as the
/// tumour approaches its plateau, so a small residue regrows faster in
/// relative terms than a large mass does. `b` is the rate constant and
/// `carrying` the plateau burden.
#[must_use]
pub fn gompertz_regrowth(burden: f64, carrying: f64, rate_per_day: f64, days: f64) -> f64 {
    if burden <= 0.0 || carrying <= burden || rate_per_day <= 0.0 || days <= 0.0 {
        return burden.max(0.0);
    }
    let ln_ratio = (carrying / burden).ln();
    // The `min` is defensive and a mutation sweep confirms no test can kill
    // it: the Gompertz form approaches `carrying` from below and never
    // crosses it, so the clamp only ever catches floating-point overshoot.
    // Left in, and said out loud, rather than removed on the strength of an
    // argument about arithmetic.
    (burden * (ln_ratio * (1.0 - (-rate_per_day * days).exp())).exp()).min(carrying)
}

/// A course of chemotherapy: equal doses at a fixed interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Regimen {
    /// Number of cycles.
    pub cycles: u32,
    /// Dose per cycle, in the same arbitrary units as `potency`.
    pub dose: f64,
    /// Days between cycles.
    pub interval_days: f64,
}

impl Regimen {
    /// Total dose delivered.
    #[must_use]
    pub fn total_dose(&self) -> f64 {
        f64::from(self.cycles) * self.dose
    }

    /// Days from the first dose to the last.
    #[must_use]
    pub fn span_days(&self) -> f64 {
        f64::from(self.cycles.saturating_sub(1)) * self.interval_days
    }
}

/// Burden remaining after a regimen, with regrowth between cycles.
///
/// Each cycle kills a fraction and reshapes the phase distribution; between
/// cycles the survivors regrow and the distribution relaxes back toward its
/// starting shape, which is redistribution -- the third of radiotherapy's four
/// Rs, arriving here because chemotherapy is where the cycle had to be built.
#[must_use]
pub fn regimen_burden(
    regimen: &Regimen,
    potency: f64,
    class: ChemoClass,
    start: &PhaseDistribution,
    initial_burden: f64,
    carrying: f64,
    regrowth_per_day: f64,
    redistribution_per_day: f64,
) -> f64 {
    let mut burden = initial_burden;
    let mut dist = *start;
    for cycle in 0..regimen.cycles {
        burden *= surviving_fraction(regimen.dose, potency, class, &dist);
        dist = surviving_distribution(regimen.dose, potency, class, &dist);
        if cycle + 1 < regimen.cycles {
            burden = gompertz_regrowth(burden, carrying, regrowth_per_day, regimen.interval_days);
            dist = redistribute(&dist, start, redistribution_per_day, regimen.interval_days);
        }
    }
    burden.max(0.0)
}

/// Relaxation of a depleted phase distribution back toward its steady shape.
///
/// Cells that survived an S-phase agent are enriched in G1 and G0; left alone
/// they move on, and the population's profile returns toward what it was. The
/// relaxation is exponential with a half-life set by `rate_per_day`, which is
/// UNCALIBRATED -- the direction (a gap between doses restores sensitivity) is
/// the result.
#[must_use]
pub fn redistribute(
    current: &PhaseDistribution,
    steady: &PhaseDistribution,
    rate_per_day: f64,
    days: f64,
) -> PhaseDistribution {
    if rate_per_day <= 0.0 || days <= 0.0 {
        return *current;
    }
    let w = 1.0 - (-rate_per_day * days).exp();
    let mix = |a: f64, b: f64| a + (b - a) * w;
    PhaseDistribution {
        g1: mix(current.g1, steady.g1),
        s: mix(current.s, steady.s),
        g2m: mix(current.g2m, steady.g2m),
        g0: mix(current.g0, steady.g0),
    }
    .normalised()
}

/// How much better a dose-dense schedule does than a conventional one at the
/// SAME total dose, as a ratio of final burdens.
///
/// **The layer's second prediction.** Norton and Simon argued that if a
/// tumour regrows between cycles, giving the same total dose in less time
/// kills more; CALGB 9741 (PMID 12668651) is the randomised trial that found
/// the shorter interval better in early breast cancer.
///
/// The number this returns is not the trial's effect size and cannot be: the
/// potency is a placeholder and no dose here is in milligrams. What it is is a
/// SIGN and a dependence.
///
/// **The dependence has two ends, and only one of them was expected.** The
/// advantage vanishes when the tumour does not regrow -- there is nothing to
/// outrun, which is the obvious half. It ALSO vanishes when regrowth is fast
/// relative to the interval, because then both schedules' tumours return to
/// the Gompertz plateau between cycles and the extra week costs nothing that
/// was not already lost. So the advantage is largest at INTERMEDIATE regrowth,
/// where the interval and the doubling are comparable, and a model that
/// reported a dose-density benefit for a tumour at either extreme would be
/// reporting an artifact of its own arithmetic.
#[must_use]
#[allow(clippy::too_many_arguments)]
pub fn dose_density_advantage(
    cycles: u32,
    dose: f64,
    conventional_interval_days: f64,
    dense_interval_days: f64,
    potency: f64,
    class: ChemoClass,
    start: &PhaseDistribution,
    initial_burden: f64,
    carrying: f64,
    regrowth_per_day: f64,
    redistribution_per_day: f64,
) -> f64 {
    let run = |interval| {
        regimen_burden(
            &Regimen {
                cycles,
                dose,
                interval_days: interval,
            },
            potency,
            class,
            start,
            initial_burden,
            carrying,
            regrowth_per_day,
            redistribution_per_day,
        )
    };
    let conventional = run(conventional_interval_days);
    let dense = run(dense_interval_days);
    if dense <= 0.0 {
        return f64::INFINITY;
    }
    conventional / dense
}

/// Neutrophil count as a share of baseline, `days` after a dose.
///
/// The constraint that makes "more, sooner" the wrong answer. Myelosuppression
/// is the dose-limiting toxicity of most cytotoxic regimens, and a model
/// without it will always prefer a shorter interval -- which is how a
/// simulation ends up recommending a schedule nobody could receive.
///
/// A fall to a nadir followed by recovery, both exponential. The timings are
/// conventional (a nadir around the second week, recovery by the third) and
/// are NOT fitted here.
#[must_use]
pub fn marrow_recovery(
    days_since_dose: f64,
    nadir_day: f64,
    recovery_rate_per_day: f64,
    nadir_depth: f64,
) -> f64 {
    let depth = nadir_depth.clamp(0.0, 1.0);
    if days_since_dose <= 0.0 {
        return 1.0;
    }
    if days_since_dose <= nadir_day {
        return 1.0 - depth * (days_since_dose / nadir_day.max(1e-9));
    }
    let since = days_since_dose - nadir_day;
    (1.0 - depth * (-recovery_rate_per_day * since).exp()).clamp(0.0, 1.0)
}

/// The shortest interval at which the marrow has recovered to `threshold`.
///
/// What a dose-density argument has to clear before it means anything. Returns
/// `None` when no interval up to `max_days` reaches the threshold, which is
/// the case a growth-factor-supported regimen exists to change.
#[must_use]
pub fn min_safe_interval_days(
    threshold: f64,
    nadir_day: f64,
    recovery_rate_per_day: f64,
    nadir_depth: f64,
    max_days: f64,
) -> Option<f64> {
    let mut d = nadir_day;
    while d <= max_days {
        if marrow_recovery(d, nadir_day, recovery_rate_per_day, nadir_depth) >= threshold {
            return Some(d);
        }
        d += 0.25;
    }
    None
}

/// Multiplier on the effective dose from drug efflux.
///
/// P-glycoprotein and its relatives pump a substrate back out, so the cell
/// sees less than was given. A separate axis from the phase specificity above
/// and from the ferroptosis-side defences the rest of this crate models: a
/// cell can be sensitive in phase and still survive because the drug never
/// reached a lethal concentration inside it.
#[must_use]
pub fn efflux_dose_factor(efflux_activity: f64) -> f64 {
    1.0 / (1.0 + efflux_activity.max(0.0))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn both_phase_distributions_are_distributions() {
        for d in [
            PhaseDistribution::proliferating(),
            PhaseDistribution::quiescent_rich(),
        ] {
            assert!((d.total() - 1.0).abs() < 1e-12, "{d:?} does not sum to 1");
            assert!((d.growth_fraction() + d.g0 - 1.0).abs() < 1e-12);
        }
        // And they differ in the way their names claim, which is the only
        // reason there are two of them.
        assert!(
            PhaseDistribution::proliferating().growth_fraction()
                > 2.0 * PhaseDistribution::quiescent_rich().growth_fraction()
        );
    }

    #[test]
    fn each_class_is_most_sensitive_where_its_name_says() {
        use CellCyclePhase::{G0, G1, G2M, S};
        let s_agent = ChemoClass::SPhaseSpecific;
        assert!(phase_sensitivity(s_agent, S) > phase_sensitivity(s_agent, G1));
        assert!(phase_sensitivity(s_agent, G1) > phase_sensitivity(s_agent, G0));
        let m_agent = ChemoClass::MPhaseSpecific;
        assert!(phase_sensitivity(m_agent, G2M) > phase_sensitivity(m_agent, S));
        assert!(phase_sensitivity(m_agent, S) > phase_sensitivity(m_agent, G0));
        // The nonspecific agent is the control: equal across the cycling
        // phases, and reduced but NOT zero out of cycle.
        let flat = ChemoClass::PhaseNonspecific;
        assert!((phase_sensitivity(flat, G1) - phase_sensitivity(flat, S)).abs() < 1e-12);
        assert!((phase_sensitivity(flat, S) - phase_sensitivity(flat, G2M)).abs() < 1e-12);
        assert!(phase_sensitivity(flat, G0) < phase_sensitivity(flat, S));
        assert!(
            phase_sensitivity(flat, G0) > 0.0,
            "a quiescent cell is less sensitive to an alkylator, not immune"
        );
    }

    #[test]
    fn no_dose_kills_nothing() {
        let d = PhaseDistribution::proliferating();
        for c in [
            ChemoClass::PhaseNonspecific,
            ChemoClass::SPhaseSpecific,
            ChemoClass::MPhaseSpecific,
        ] {
            assert!((surviving_fraction(0.0, 1.0, c, &d) - 1.0).abs() < 1e-15);
            assert!((surviving_fraction(5.0, 0.0, c, &d) - 1.0).abs() < 1e-15);
        }
    }

    #[test]
    fn a_phase_specific_agent_leaves_more_behind_than_a_flat_one() {
        // The layer's first prediction, and it is emergent: nothing in the
        // sensitivity table says "leave a residue".
        let prolif = PhaseDistribution::proliferating();
        for dose in [1.0, 4.0, 8.0, 16.0] {
            for c in [ChemoClass::SPhaseSpecific, ChemoClass::MPhaseSpecific] {
                let ratio = phase_specific_residue_ratio(dose, 0.5, c, &prolif);
                assert!(ratio > 1.0, "{c:?} at dose {dose} left ratio {ratio}");
            }
        }
        // And it GROWS with dose, which is the shape the prediction is about:
        // raising the dose helps the flat agent and stops helping the
        // phase-specific one.
        let low = phase_specific_residue_ratio(1.0, 0.5, ChemoClass::SPhaseSpecific, &prolif);
        let high = phase_specific_residue_ratio(16.0, 0.5, ChemoClass::SPhaseSpecific, &prolif);
        assert!(
            high > 3.0 * low,
            "the gap did not widen with dose: {low} -> {high}"
        );
    }

    #[test]
    fn the_residue_gap_narrows_when_most_cells_are_out_of_cycle() {
        // The comparison that was NOT designed in, and the more interesting
        // half: in a quiescent tumour the flat agent is struggling too, so the
        // phase-specific agent's relative disadvantage shrinks. A model that
        // reported the same gap in both populations would be reporting its own
        // sensitivity table rather than a property of the tumour.
        let c = ChemoClass::SPhaseSpecific;
        let in_culture =
            phase_specific_residue_ratio(8.0, 0.5, c, &PhaseDistribution::proliferating());
        let in_tumour =
            phase_specific_residue_ratio(8.0, 0.5, c, &PhaseDistribution::quiescent_rich());
        assert!(in_tumour < in_culture, "{in_tumour} vs {in_culture}");
        assert!(in_tumour > 1.0);
    }

    #[test]
    fn unlimited_dose_kills_everything_in_every_class() {
        // A RETRACTION, pinned. The first version of this module reported the
        // survival at an enormous dose as a "plateau height" and claimed it
        // was substantially above zero for a phase-specific agent. It is zero
        // for every class, because no phase here has exactly zero sensitivity
        // -- and a phase given a hard zero would produce a plateau BY
        // ASSUMPTION. This test exists so nobody writes that claim again.
        let d = PhaseDistribution::quiescent_rich();
        for c in [
            ChemoClass::PhaseNonspecific,
            ChemoClass::SPhaseSpecific,
            ChemoClass::MPhaseSpecific,
        ] {
            assert!(
                surviving_fraction(1.0e6, 1.0, c, &d) < 1e-12,
                "{c:?} left a residue at unlimited dose"
            );
        }
    }

    #[test]
    fn the_survivors_are_enriched_for_the_phases_the_drug_cannot_reach() {
        let before = PhaseDistribution::proliferating();
        let after = surviving_distribution(8.0, 0.5, ChemoClass::SPhaseSpecific, &before);
        assert!(
            after.s < before.s,
            "S was not depleted: {} -> {}",
            before.s,
            after.s
        );
        assert!(
            after.g0 > before.g0,
            "G0 was not enriched: {} -> {}",
            before.g0,
            after.g0
        );
        assert!(
            (after.total() - 1.0).abs() < 1e-12,
            "the result is not a distribution"
        );
        // The flat agent barely reshapes anything, which is the control.
        let flat = surviving_distribution(8.0, 0.5, ChemoClass::PhaseNonspecific, &before);
        assert!((flat.s - before.s).abs() < (before.s - after.s).abs());
    }

    #[test]
    fn gompertz_regrowth_saturates_and_is_inert_without_a_rate() {
        let b = 1.0e9;
        let cap = 1.0e12;
        assert!((gompertz_regrowth(b, cap, 0.0, 30.0) - b).abs() < 1e-6);
        assert!((gompertz_regrowth(b, cap, 0.05, 0.0) - b).abs() < 1e-6);
        let short = gompertz_regrowth(b, cap, 0.05, 7.0);
        let long = gompertz_regrowth(b, cap, 0.05, 60.0);
        assert!(short > b && long > short, "{b} {short} {long}");
        assert!(
            long <= cap * (1.0 + 1e-12),
            "regrowth passed the carrying capacity"
        );
    }

    #[test]
    fn dose_density_helps_only_when_there_is_something_to_outrun() {
        // THE LAYER'S SECOND PREDICTION, and the shape has two ends. The
        // advantage vanishes with no regrowth -- nothing to outrun -- AND with
        // very fast regrowth, because both schedules return to the plateau
        // between cycles. Only the first end was expected; the second is what
        // the model produced.
        let d = PhaseDistribution::proliferating();
        let advantage = |regrowth| {
            dose_density_advantage(
                6,
                2.0,
                21.0,
                14.0,
                0.5,
                ChemoClass::PhaseNonspecific,
                &d,
                1.0e9,
                1.0e12,
                regrowth,
                0.2,
            )
        };
        let none = advantage(0.0);
        let middle = advantage(0.01);
        let fast = advantage(0.4);
        assert!(
            middle > 2.0,
            "no dose-density advantage at intermediate regrowth: {middle}"
        );
        assert!(
            none < 1.05,
            "an advantage appeared with nothing to outrun: {none}"
        );
        assert!(
            fast < 1.2,
            "the advantage did not vanish at fast regrowth: {fast}"
        );
        assert!(
            middle > none && middle > fast,
            "the advantage is not largest in the middle: {none} {middle} {fast}"
        );
    }

    #[test]
    fn redistribution_moves_the_population_back_toward_its_steady_shape() {
        let steady = PhaseDistribution::proliferating();
        let depleted = surviving_distribution(8.0, 0.5, ChemoClass::SPhaseSpecific, &steady);
        let relaxed = redistribute(&depleted, &steady, 0.2, 14.0);
        assert!(
            relaxed.s > depleted.s,
            "S did not recover: {} -> {}",
            depleted.s,
            relaxed.s
        );
        assert!(relaxed.s <= steady.s + 1e-12, "S overshot its steady share");
        // Inert without a rate or without time, which is what lets a regimen
        // with no gap be the same model as one with a gap.
        for (rate, days) in [(0.0, 14.0), (0.2, 0.0)] {
            let same = redistribute(&depleted, &steady, rate, days);
            assert!((same.s - depleted.s).abs() < 1e-12);
        }
    }

    #[test]
    fn the_marrow_falls_to_a_nadir_and_comes_back() {
        let (nadir_day, rate, depth) = (12.0, 0.35, 0.8);
        assert!((marrow_recovery(0.0, nadir_day, rate, depth) - 1.0).abs() < 1e-12);
        let at_nadir = marrow_recovery(nadir_day, nadir_day, rate, depth);
        assert!(
            (at_nadir - (1.0 - depth)).abs() < 1e-9,
            "nadir was {at_nadir}"
        );
        assert!(marrow_recovery(6.0, nadir_day, rate, depth) > at_nadir);
        assert!(marrow_recovery(21.0, nadir_day, rate, depth) > at_nadir);
        assert!(marrow_recovery(60.0, nadir_day, rate, depth) > 0.99);
    }

    #[test]
    fn a_schedule_has_to_clear_the_marrow_before_it_means_anything() {
        // The constraint that stops the model preferring an interval nobody
        // could receive. Without it, `dose_density_advantage` is an argument
        // for giving everything at once.
        let (nadir_day, rate, depth) = (12.0, 0.35, 0.8);
        let safe = min_safe_interval_days(0.8, nadir_day, rate, depth, 60.0)
            .expect("recovery to 80% should be reachable within sixty days");
        assert!(safe > nadir_day, "the safe interval is before the nadir");
        assert!(marrow_recovery(safe, nadir_day, rate, depth) >= 0.8);
        assert!(
            marrow_recovery(safe - 1.0, nadir_day, rate, depth) < 0.8,
            "an earlier day would also have cleared the threshold"
        );
        // And a threshold that cannot be reached returns None rather than a
        // plausible-looking number.
        assert!(min_safe_interval_days(1.0, nadir_day, rate, depth, 30.0).is_none());
    }

    #[test]
    fn efflux_reduces_the_dose_the_cell_sees_and_never_raises_it() {
        assert!((efflux_dose_factor(0.0) - 1.0).abs() < 1e-12);
        let mut last = 1.0;
        for a in [0.5, 1.0, 2.0, 8.0] {
            let f = efflux_dose_factor(a);
            assert!(f < last && f > 0.0, "activity {a} gave {f}");
            last = f;
        }
    }

    /// The two points the Python mirror pins on its side, to four decimals.
    ///
    /// `scripts/validate_chemo.py` implements this model AGAIN, in stdlib
    /// Python, and its committed artifact carries these same numbers. Neither
    /// side reads the other: the Python parses the sensitivity table out of
    /// this file and re-derives the curve, and a guard there asserts these
    /// literals are still here. So a change in either implementation breaks a
    /// test in the other, which is the only version of "two implementations
    /// agree" that is worth anything.
    #[test]
    fn the_dose_response_matches_the_independent_python_mirror() {
        let d = PhaseDistribution::proliferating();
        let s_phase = surviving_fraction(8.0, 0.5, ChemoClass::SPhaseSpecific, &d);
        let flat = surviving_fraction(8.0, 0.5, ChemoClass::PhaseNonspecific, &d);
        assert!((s_phase - 0.4531).abs() < 5e-5, "S-phase agent: {s_phase}");
        assert!(
            (flat - 0.0275).abs() < 5e-5,
            "phase-nonspecific agent: {flat}"
        );
    }

    #[test]
    fn a_regimen_kills_more_than_one_cycle_of_it() {
        let d = PhaseDistribution::proliferating();
        let one = regimen_burden(
            &Regimen {
                cycles: 1,
                dose: 2.0,
                interval_days: 21.0,
            },
            0.5,
            ChemoClass::PhaseNonspecific,
            &d,
            1.0e9,
            1.0e12,
            0.0,
            0.2,
        );
        let six = regimen_burden(
            &Regimen {
                cycles: 6,
                dose: 2.0,
                interval_days: 21.0,
            },
            0.5,
            ChemoClass::PhaseNonspecific,
            &d,
            1.0e9,
            1.0e12,
            0.0,
            0.2,
        );
        assert!(six < one, "six cycles left more than one: {six} vs {one}");
        assert!(six > 0.0);
        let r = Regimen {
            cycles: 6,
            dose: 2.0,
            interval_days: 21.0,
        };
        assert!((r.total_dose() - 12.0).abs() < 1e-12);
        assert!((r.span_days() - 105.0).abs() < 1e-12);
    }
}
