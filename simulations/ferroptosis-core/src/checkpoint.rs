//! Checkpoint blockade: the arm that owned no module, and the one whose
//! calibration failed for a reason worth fixing rather than restating.
//!
//! # Why this module exists
//!
//! Until it did, checkpoint blockade reached the panel entirely through the
//! shared immune machinery -- `immune` and `immune_spatial` -- which the depth
//! report credits to no arm, correctly, because several arms reach through
//! them. So the arm with the largest immunotherapy literature in the taxonomy
//! had zero lines of its own.
//!
//! # The calibration problem, and the way out
//!
//! `analysis/modality-calibration.md` records why fitting this arm failed. With
//! no ferroptotic death the cascade collapses to `kill = C * antigenicity`, so
//! a published response rate constrains the PRODUCT `antigenicity * kill_rate`
//! and neither factor can be recovered from it. On top of that the mapping is
//! weak on its own terms: an objective response is a 30% reduction in diameter,
//! not a dead tumour.
//!
//! **A RATIO between two strata of one trial removes the second problem and
//! makes progress on the first.** If the observed response rate is any linear
//! function of the model's kill, `ORR = m * kill`, then a ratio of two response
//! rates from the same trial, drug and endpoint cancels `m` exactly:
//!
//! ```text
//! ORR_high / ORR_low = kill_high / kill_low
//! ```
//!
//! The unknown mapping constant divides out, and what is left is a statement
//! the model can fail. KEYNOTE-158 (PMID 32919526) supplies exactly such a
//! pair, stratified by tumour mutational burden -- which moves ANTIGENICITY and
//! not, to first order, the brake. So the ratio constrains the shape of the
//! antigenicity response, which the absolute band never could.
//!
//! **What it does not do** is make both factors identifiable. One ratio is one
//! equation. It constrains the antigenicity axis and leaves the brake where it
//! was, and this module says so rather than claiming the row is now fitted.
//!
//! # The three resistance modes
//!
//! Primary, adaptive and acquired resistance are the taxonomy the clinic uses,
//! and they fail in different places in this model: a primary-resistant tumour
//! never presents to a T cell, an adaptively resistant one raises its own brake
//! in response to interferon it provoked, and an acquired-resistant one has
//! lost the machinery to be seen at all (B2M) or to respond to interferon
//! (JAK1/2) -- Zaretsky 2016, PMID 27433843.

use serde::{Deserialize, Serialize};

/// The tumour, as far as checkpoint blockade is concerned.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct CheckpointTumour {
    /// PD-L1 tumour proportion score, as a fraction in `[0, 1]`.
    ///
    /// Raises the brake. It is a modifier of the BRAKE and not of the
    /// antigenicity, which is the separation the whole module is built around.
    pub pd_l1_score: f64,
    /// Tumour mutational burden, mutations per megabase.
    ///
    /// Raises antigenicity, and not the brake.
    pub tmb_per_mb: f64,
    /// Fraction of the tumour still able to present antigen.
    ///
    /// `0.0` is a B2M-null tumour: invisible whatever its mutational burden,
    /// which is the structural claim this module can be checked on.
    pub antigen_presentation: f64,
    /// Whether the tumour can respond to interferon gamma (JAK1/2 intact).
    ///
    /// A tumour that cannot is protected from the interferon-driven arm of the
    /// response -- and, in the same stroke, cannot mount adaptive PD-L1
    /// induction either. Both directions follow from the one lesion, which is
    /// why it is a flag rather than two knobs.
    pub interferon_responsive: bool,
}

impl Default for CheckpointTumour {
    /// A tumour with no particular feature: mid PD-L1, low TMB, presentation
    /// intact.
    fn default() -> Self {
        Self {
            pd_l1_score: 0.1,
            tmb_per_mb: 3.0,
            antigen_presentation: 1.0,
            interferon_responsive: true,
        }
    }
}

/// The tumour mutational burden above which KEYNOTE-158 called a tumour
/// tTMB-high, in mutations per megabase.
pub const TMB_HIGH_THRESHOLD_PER_MB: f64 = 10.0;

/// Mutations per megabase at which the antigenicity response is half of its
/// maximum.
///
/// The shape parameter the ratio target constrains. Antigenicity cannot grow
/// without bound -- a T-cell repertoire and an antigen-presentation machinery
/// are finite -- so the response saturates, and where it saturates decides how
/// much a TMB-high tumour gains over a TMB-low one.
pub const TMB_HALF_MAX_PER_MB: f64 = 24.0;

/// Antigenicity from mutational burden, saturating.
///
/// `a(t) = floor + (1 - floor) * t / (t + k)`. The floor is what a tumour
/// presents without a high mutational burden -- viral antigens, cancer-testis
/// antigens, overexpressed self -- and it is why a TMB-low tumour responds at
/// all rather than never.
#[must_use]
pub fn antigenicity_from_tmb(tmb_per_mb: f64, half_max: f64, floor: f64) -> f64 {
    let t = tmb_per_mb.max(0.0);
    let k = half_max.max(f64::MIN_POSITIVE);
    let f = floor.clamp(0.0, 1.0);
    (f + (1.0 - f) * t / (t + k)).clamp(0.0, 1.0)
}

/// Fraction of receptors occupied by drug, `d / (d + Kd)`.
///
/// The pharmacology the arm previously replaced with a bare efficacy scalar.
/// Saturating, so the difference between a full dose and twice a full dose is
/// small -- which is the reason checkpoint dosing is flat rather than titrated,
/// and a property the arm could not express before.
#[must_use]
pub fn receptor_occupancy(dose: f64, kd: f64) -> f64 {
    let d = dose.max(0.0);
    let k = kd.max(f64::MIN_POSITIVE);
    (d / (d + k)).clamp(0.0, 1.0)
}

/// The brake on T-cell killing, after PD-L1 expression raises it and drug
/// occupancy releases it.
///
/// `brake = base * (1 + pdl1_gain * pd_l1) * (1 - occupancy)`, clamped into
/// `[0, 1]`. Returns the RESIDUAL brake: 0 is a fully released T cell.
#[must_use]
pub fn residual_brake(base_brake: f64, pd_l1_score: f64, pdl1_gain: f64, occupancy: f64) -> f64 {
    let raised = base_brake * (1.0 + pdl1_gain * pd_l1_score.clamp(0.0, 1.0));
    (raised * (1.0 - occupancy.clamp(0.0, 1.0))).clamp(0.0, 1.0)
}

/// PD-L1 induced by the interferon the response itself produces.
///
/// Adaptive resistance: a tumour that is being attacked raises the very brake
/// the drug is there to release. It is the mechanism that makes PD-L1 a
/// consequence of immune pressure as well as a cause of escape, and it is why
/// PD-L1 measured before treatment is a weaker predictor than its biology
/// suggests. A tumour that cannot respond to interferon cannot do this either.
#[must_use]
pub fn adaptive_pd_l1(
    baseline_pd_l1: f64,
    local_interferon: f64,
    induction: f64,
    interferon_responsive: bool,
) -> f64 {
    if !interferon_responsive {
        return baseline_pd_l1.clamp(0.0, 1.0);
    }
    (baseline_pd_l1 + induction * local_interferon.max(0.0)).clamp(0.0, 1.0)
}

/// The model's response index for a tumour under a given blockade.
///
/// NOT a response rate and not a kill fraction: a dimensionless quantity
/// proportional to the model's kill, which is what makes the ratio below
/// checkable while the absolute value is not. Anything that scales it -- the
/// dendritic-cell rate, the priming rate, the per-cell kill rate, the mapping
/// from kill to a radiological response -- cancels in a ratio and is therefore
/// deliberately absent here.
#[must_use]
pub fn response_index(
    tumour: &CheckpointTumour,
    occupancy: f64,
    base_brake: f64,
    pdl1_gain: f64,
    tmb_half_max: f64,
    antigenicity_floor: f64,
) -> f64 {
    let antigenicity = antigenicity_from_tmb(tumour.tmb_per_mb, tmb_half_max, antigenicity_floor);
    let presented = antigenicity * tumour.antigen_presentation.clamp(0.0, 1.0);
    let brake = residual_brake(base_brake, tumour.pd_l1_score, pdl1_gain, occupancy);
    (presented * (1.0 - brake)).max(0.0)
}

/// The ratio of response indices between two tumours.
///
/// **This is the quantity a trial can refute.** If the observed response rate
/// is any linear function of the model's kill, the unknown constant cancels
/// here, so a ratio computed from two strata of one trial is comparable to a
/// ratio computed from this model -- while neither absolute value is.
///
/// Returns `None` when the reference tumour has no response at all, where a
/// ratio is not defined rather than infinite. A B2M-null reference is exactly
/// that case.
#[must_use]
#[allow(clippy::too_many_arguments)]
pub fn response_ratio(
    high: &CheckpointTumour,
    low: &CheckpointTumour,
    occupancy: f64,
    base_brake: f64,
    pdl1_gain: f64,
    tmb_half_max: f64,
    antigenicity_floor: f64,
) -> Option<f64> {
    let denom = response_index(
        low,
        occupancy,
        base_brake,
        pdl1_gain,
        tmb_half_max,
        antigenicity_floor,
    );
    if denom <= 0.0 {
        return None;
    }
    let numer = response_index(
        high,
        occupancy,
        base_brake,
        pdl1_gain,
        tmb_half_max,
        antigenicity_floor,
    );
    Some(numer / denom)
}

/// Which resistance mode a tumour is in, if any.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum ResistanceMode {
    /// Nothing structural in the way.
    None,
    /// Never presented enough to be seen: the cold tumour.
    Primary,
    /// Raises its own brake under interferon pressure.
    Adaptive,
    /// Has lost the machinery -- presentation or interferon response.
    Acquired,
}

/// Classify a tumour's resistance mode.
///
/// The order matters and is a judgement: a tumour that has lost presentation is
/// reported as ACQUIRED even if its mutational burden is also low, because the
/// lesion is the more specific fact about it. A page reporting the milder mode
/// for a tumour that has both would understate what is wrong with it.
#[must_use]
pub fn resistance_mode(
    tumour: &CheckpointTumour,
    antigenicity_floor: f64,
    primary_threshold: f64,
    adaptive_threshold: f64,
) -> ResistanceMode {
    if tumour.antigen_presentation < 1.0 || !tumour.interferon_responsive {
        return ResistanceMode::Acquired;
    }
    let a = antigenicity_from_tmb(tumour.tmb_per_mb, TMB_HALF_MAX_PER_MB, antigenicity_floor);
    if a * tumour.antigen_presentation < primary_threshold {
        return ResistanceMode::Primary;
    }
    if tumour.pd_l1_score >= adaptive_threshold {
        return ResistanceMode::Adaptive;
    }
    ResistanceMode::None
}

#[cfg(test)]
mod tests {
    use super::*;

    const FLOOR: f64 = 0.05;
    const OCC: f64 = 0.909;
    const BRAKE: f64 = 0.6;
    const GAIN: f64 = 2.0;

    #[test]
    fn antigenicity_saturates_rather_than_growing_without_bound() {
        // The shape the trial ratio constrains. A linear antigenicity would
        // make a hypermutated tumour arbitrarily more responsive than a
        // moderately mutated one, which no series reports.
        let a = |t: f64| antigenicity_from_tmb(t, TMB_HALF_MAX_PER_MB, FLOOR);
        assert!(
            (a(0.0) - FLOOR).abs() < 1e-12,
            "no burden should give the floor"
        );
        assert!(a(1000.0) < 1.0 && a(1000.0) > 0.9);
        // Concave: equal steps in burden buy less and less.
        let first = a(10.0) - a(0.0);
        let second = a(20.0) - a(10.0);
        assert!(
            second < first,
            "the response is not saturating: {first} {second}"
        );
    }

    #[test]
    fn occupancy_saturates_which_is_why_the_dosing_is_flat() {
        assert!((receptor_occupancy(0.0, 1.0)).abs() < 1e-12);
        let single = receptor_occupancy(1.0, 0.1);
        let double = receptor_occupancy(2.0, 0.1);
        assert!(
            single > 0.9,
            "a standard dose should nearly saturate: {single}"
        );
        assert!(
            double - single < 0.05,
            "doubling the dose moved occupancy by {}",
            double - single
        );
    }

    #[test]
    fn occupancy_is_half_at_the_dissociation_constant() {
        // A MUTATION SURVIVOR. The saturation test above uses a dose ten times
        // Kd, where a LINEAR occupancy also reads 1.0 because it is clamped --
        // so replacing the Michaelis form with `d / Kd` left both assertions
        // green. At d = Kd the two forms differ by a factor of two and the
        // clamp cannot hide it.
        let at_kd = receptor_occupancy(0.5, 0.5);
        assert!(
            (at_kd - 0.5).abs() < 1e-12,
            "occupancy at the dissociation constant is {at_kd}, not one half"
        );
        let below = receptor_occupancy(0.25, 0.5);
        assert!((below - (1.0 / 3.0)).abs() < 1e-12, "{below}");
    }

    #[test]
    fn a_stronger_brake_lowers_the_response_index() {
        // ALSO A MUTATION SURVIVOR, and a worse one: every other test here
        // either cancels the brake (the ratio tests, deliberately) or reads
        // zero whatever it does (the B2M test), so multiplying by the brake
        // instead of by its complement -- inverting the drug's entire effect --
        // was invisible. This pins the DIRECTION at the one place it matters.
        let t = CheckpointTumour {
            tmb_per_mb: 20.0,
            ..CheckpointTumour::default()
        };
        let released = response_index(&t, 1.0, BRAKE, GAIN, TMB_HALF_MAX_PER_MB, FLOOR);
        let braked = response_index(&t, 0.0, BRAKE, GAIN, TMB_HALF_MAX_PER_MB, FLOOR);
        assert!(
            released > braked,
            "full drug occupancy did not raise the response: {released} vs {braked}"
        );
        let expressing = CheckpointTumour {
            pd_l1_score: 1.0,
            ..t
        };
        let bare = CheckpointTumour {
            pd_l1_score: 0.0,
            ..t
        };
        let a = response_index(&expressing, 0.5, BRAKE, GAIN, TMB_HALF_MAX_PER_MB, FLOOR);
        let b = response_index(&bare, 0.5, BRAKE, GAIN, TMB_HALF_MAX_PER_MB, FLOOR);
        assert!(
            a < b,
            "PD-L1 expression did not lower the response: {a} vs {b}"
        );
    }

    #[test]
    fn pd_l1_raises_the_brake_and_drug_releases_it() {
        let low = residual_brake(BRAKE, 0.0, GAIN, 0.0);
        let high = residual_brake(BRAKE, 1.0, GAIN, 0.0);
        assert!(high > low, "PD-L1 did not raise the brake: {low} {high}");
        let released = residual_brake(BRAKE, 1.0, GAIN, 1.0);
        assert!(
            released.abs() < 1e-12,
            "full occupancy left a brake: {released}"
        );
        // And the brake is bounded: a heavily expressing tumour cannot have a
        // residual brake above 1, which would make the kill negative.
        assert!(residual_brake(0.9, 1.0, 5.0, 0.0) <= 1.0);
    }

    #[test]
    fn a_tumour_that_cannot_present_does_not_respond_at_any_burden() {
        // THE STRUCTURAL CHECK, and it is not a fit. B2M loss is one of the
        // acquired-resistance lesions Zaretsky 2016 reports, and a model in
        // which a hypermutated tumour responds despite it would be wrong about
        // the mechanism rather than imprecise about a number.
        let b2m = CheckpointTumour {
            tmb_per_mb: 40.0,
            antigen_presentation: 0.0,
            ..CheckpointTumour::default()
        };
        let idx = response_index(&b2m, OCC, BRAKE, GAIN, TMB_HALF_MAX_PER_MB, FLOOR);
        assert!(idx.abs() < 1e-12, "a B2M-null tumour responded: {idx}");
        // A ratio against it is UNDEFINED rather than infinite, which is the
        // difference between reporting a limit and reporting a number.
        let normal = CheckpointTumour::default();
        assert!(
            response_ratio(&normal, &b2m, OCC, BRAKE, GAIN, TMB_HALF_MAX_PER_MB, FLOOR).is_none()
        );
    }

    #[test]
    fn the_ratio_cancels_everything_the_two_strata_share() {
        // THE ARGUMENT THE WHOLE CALIBRATION RESTS ON, asserted rather than
        // explained: two tumours differing only in burden must give the same
        // ratio at ANY brake, any PD-L1 and any occupancy, because those
        // factors multiply both sides. If that ever stops being true the
        // trial comparison stops being valid, and it should fail here first.
        let high = CheckpointTumour {
            tmb_per_mb: 20.0,
            ..CheckpointTumour::default()
        };
        let low = CheckpointTumour {
            tmb_per_mb: 3.0,
            ..CheckpointTumour::default()
        };
        let reference =
            response_ratio(&high, &low, OCC, BRAKE, GAIN, TMB_HALF_MAX_PER_MB, FLOOR).unwrap();
        for (occ, brake, gain) in [(0.5, 0.2, 0.5), (0.99, 0.9, 4.0), (0.1, 0.5, 1.0)] {
            let r =
                response_ratio(&high, &low, occ, brake, gain, TMB_HALF_MAX_PER_MB, FLOOR).unwrap();
            assert!(
                (r - reference).abs() < 1e-9,
                "the ratio moved with the brake ({occ} {brake} {gain}): \
                     {r} vs {reference}"
            );
        }
        // And it does NOT cancel a difference in PD-L1, which is the limit of
        // the argument: strata that differ in the brake are not comparable
        // this way, and the calibration page says so.
        let high_pdl1 = CheckpointTumour {
            pd_l1_score: 0.9,
            ..high
        };
        let r = response_ratio(
            &high_pdl1,
            &low,
            OCC,
            BRAKE,
            GAIN,
            TMB_HALF_MAX_PER_MB,
            FLOOR,
        )
        .unwrap();
        assert!(
            (r - reference).abs() > 0.05,
            "a PD-L1 difference cancelled too, which would make the \
                 stratification argument vacuous"
        );
    }

    #[test]
    fn the_measured_ratio_is_reproduced_within_the_published_band() {
        // KEYNOTE-158 (PMID 32919526): 29% against 6%, so 4.83x, and
        // conservatively 2.63x to 7.80x from the interval endpoints. The model
        // is checked against the BAND, and the point estimate is above what it
        // returns -- which the calibration page reports rather than hides.
        let high = CheckpointTumour {
            tmb_per_mb: 20.0,
            ..CheckpointTumour::default()
        };
        let low = CheckpointTumour {
            tmb_per_mb: 3.0,
            ..CheckpointTumour::default()
        };
        let r = response_ratio(&high, &low, OCC, BRAKE, GAIN, TMB_HALF_MAX_PER_MB, FLOOR).unwrap();
        assert!(
            (2.63..=7.80).contains(&r),
            "the model's TMB response ratio is {r:.2}x, outside the band the \
                 trial's own confidence interval endpoints give"
        );
        assert!(r > 1.5, "a ratio near 1 would mean burden does not matter");
    }

    #[test]
    fn adaptive_induction_needs_an_interferon_responsive_tumour() {
        let base = 0.1;
        let induced = adaptive_pd_l1(base, 0.5, 0.8, true);
        assert!(induced > base, "interferon did not raise PD-L1: {induced}");
        let deaf = adaptive_pd_l1(base, 0.5, 0.8, false);
        assert!(
            (deaf - base).abs() < 1e-12,
            "a JAK-null tumour induced PD-L1 anyway: {deaf}"
        );
        assert!(
            adaptive_pd_l1(0.9, 10.0, 1.0, true) <= 1.0,
            "PD-L1 passed 1"
        );
    }

    #[test]
    fn each_resistance_mode_is_reachable_and_the_order_is_deliberate() {
        let cold = CheckpointTumour {
            tmb_per_mb: 0.2,
            ..CheckpointTumour::default()
        };
        assert_eq!(
            resistance_mode(&cold, 0.01, 0.15, 0.5),
            ResistanceMode::Primary
        );
        let hot_braked = CheckpointTumour {
            tmb_per_mb: 30.0,
            pd_l1_score: 0.8,
            ..CheckpointTumour::default()
        };
        assert_eq!(
            resistance_mode(&hot_braked, FLOOR, 0.15, 0.5),
            ResistanceMode::Adaptive
        );
        let lost = CheckpointTumour {
            antigen_presentation: 0.0,
            ..CheckpointTumour::default()
        };
        assert_eq!(
            resistance_mode(&lost, FLOOR, 0.15, 0.5),
            ResistanceMode::Acquired
        );
        let deaf = CheckpointTumour {
            tmb_per_mb: 30.0,
            interferon_responsive: false,
            ..CheckpointTumour::default()
        };
        assert_eq!(
            resistance_mode(&deaf, FLOOR, 0.15, 0.5),
            ResistanceMode::Acquired
        );
        let ordinary = CheckpointTumour {
            tmb_per_mb: 30.0,
            ..CheckpointTumour::default()
        };
        assert_eq!(
            resistance_mode(&ordinary, FLOOR, 0.15, 0.5),
            ResistanceMode::None
        );
        // THE ORDER IS THE JUDGEMENT: a tumour with BOTH a lost lesion and a
        // low burden is reported as Acquired, because the lesion is the more
        // specific fact. Pinned so the precedence cannot be reversed quietly.
        let both = CheckpointTumour {
            tmb_per_mb: 0.2,
            antigen_presentation: 0.0,
            ..CheckpointTumour::default()
        };
        assert_eq!(
            resistance_mode(&both, 0.01, 0.15, 0.5),
            ResistanceMode::Acquired
        );
    }
}
