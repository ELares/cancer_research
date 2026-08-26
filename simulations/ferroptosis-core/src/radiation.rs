//! Ionizing radiation as a treatment arm (#726, the half the title asks for).
//!
//! The oxygen half of #726 landed first: [`crate::oxygen::oxygen_enhancement_ratio`]
//! gave the engine's exogenous-ROS O₂ dependence the measured
//! Alper–Howard-Flanders shape. This module is the modality that consumes it.
//!
//! ## Radiation is TWO channels, and collapsing them would be a lie
//!
//! The obvious implementation — feed dose in as exogenous ROS and let the
//! existing lipid-peroxidation engine do the rest — models radiation as a
//! photodynamic agent. Radiation's dominant lethal lesion is the DNA
//! double-strand break, not membrane lipid peroxidation, and a model that
//! routes all of it through the ferroptosis machinery would predict that a
//! GPX4-competent cell is radioresistant. It is not.
//!
//! So this module supplies two channels and the consumer composes them:
//!
//! 1. **DNA-damage lethality** — [`lq_survival`], the linear-quadratic model
//!    `SF = exp(−αD − βD²)`. It does not touch [`crate::biochem::CellState`] at
//!    all: it is a per-cell survival probability the consumer rolls against.
//!    This is the channel with eighty years of clonogenic-survival data behind
//!    it, and it is the calibration target this layer arrives with.
//! 2. **Radiation-induced ferroptosis** — [`exo_ros_from_dose`], which enters
//!    through the same `exo_ros_peak` path SDT and PDT already use. Radiation
//!    does induce lipid peroxidation and ferroptosis, and lipid-metabolism
//!    manipulation does radiosensitize. The direction is anchored inside this
//!    project's own frozen corpus rather than by a citation carried in from
//!    outside -- PMID 35859734 is titled on exactly this coupling, SCD
//!    targeting enhancing radiation-induced ferroptosis and immunogenic cell
//!    death. No COUNT is quoted here: a number in a doc comment is a number
//!    nothing can check, so the size of that literature belongs in a derived
//!    analysis artifact, not in Rust prose. What the literature does NOT give
//!    is the ROS-per-gray conversion, so `ros_per_gy` is an explicitly
//!    uncalibrated knob and its CALIBRATION_STATUS row says so.
//!
//! ## Oxygen enters both channels through the same hyperbola
//!
//! Radiobiology's oxygen effect is a DOSE-MODIFYING factor: a hypoxic cell
//! needs OER-times more dose for the same survival. [`oer_scaled_dose`]
//! expresses exactly that, reusing [`crate::oxygen::oer_relative_efficacy`], so
//! the O₂ shape in this module is the same measured hyperbola the SDT path
//! already validated rather than a second, independent assumption.
//!
//! That is also what makes this arm testable rather than merely present. Fit
//! `ros_per_gy` at full O₂ and the OER is then PREDICTED, not fitted: the
//! iso-survival dose ratio between anoxia and full supply should land in the
//! published 2.5–3.0 band. If it does not, this layer is wrong and the report
//! must say so.
//!
//! ## Radiation is the physical modality with no depth limit
//!
//! [`crate::physics`] gives SDT and PDT exponential attenuation with depth, and
//! that attenuation IS this project's tissue-access argument: light dies in
//! millimetres, ultrasound reaches centimetres. Megavoltage photons attenuate
//! at roughly 3 % per centimetre of soft tissue, so over a whole tumour the
//! delivered dose is nearly flat. [`intensity_at_depth`] supplies that, which
//! gives the manuscript's depth comparison a control it has never had — the
//! one physical modality whose reach is not the limiting factor.
//!
//! ## Identity default
//!
//! [`RadiationConfig::default`] has `dose_gy = 0.0`, at which [`lq_survival`]
//! returns exactly `1.0`, [`dna_lethality`] exactly `0.0` and
//! [`exo_ros_from_dose`] exactly `0.0`. A run carrying the default config is
//! byte-identical to one with no radiation model at all.

use crate::oxygen::oer_relative_efficacy;
use crate::params::RadiationConfig;

/// Linear attenuation coefficient of 6 MV photons in soft tissue, per cm.
///
/// Roughly 3 % per centimetre. Quoted as a CONVENTIONAL value for
/// megavoltage photons in water-equivalent tissue, not as a fit performed
/// here; the point it carries is the ORDER, which is what separates radiation
/// from light and ultrasound rather than any third digit.
///
/// For contrast, `crate::physics::pdt_intensity_at_depth` falls by orders of
/// magnitude over the same distance. That contrast is the reason this constant
/// exists.
pub const MU_6MV_SOFT_TISSUE_PER_CM: f64 = 0.03;

/// α/β ratio for early-responding tumour tissue, in Gy, from this project's
/// own frozen corpus.
///
/// PMID 32307022 (`corpus/by-pmid/32307022.md`) states its radiobiology
/// explicitly for glioblastoma: `SF2 = 0.5445`, `α/β = 10 Gy`, and
/// `SF(2.67 Gy) = 0.4244`. Three numbers, which is one more than a fit needs
/// -- and that is what makes this a calibration target rather than a
/// convention. Fitting α from SF2 with `β = α/(α/β)` gives
/// [`ALPHA_GLIOBLASTOMA_PER_GY`], and PREDICTS the 2.67 Gy point without
/// using it.
///
/// The textbook split (~10 Gy for tumours and early-responding normal tissue,
/// ~3 Gy for late-responding tissue) agrees, but the number here is taken from
/// a paper in the corpus rather than from the textbook, so it can be checked.
pub const ALPHA_BETA_TUMOUR_GY: f64 = 10.0;

/// α fitted to PMID 32307022's `SF2 = 0.5445` at `α/β = 10`, per Gy.
///
/// `α = −ln(SF2) / (D + D²/(α/β))` at `D = 2`. The ONE fitted parameter of the
/// DNA channel; β follows from the ratio.
pub const ALPHA_GLIOBLASTOMA_PER_GY: f64 = 0.25329;

/// The published point this module PREDICTS rather than fits: `SF = 0.4244` at
/// `2.67 Gy` (PMID 32307022). Guarded, not decorative.
pub const SF_AT_2_67_GY_PUBLISHED: f64 = 0.4244;

/// Linear-quadratic surviving fraction, `exp(−αD − βD²)`.
///
/// The standard dose-response of radiobiology, and the reason radiation has an
/// external calibration literature that nothing else in this engine has:
/// clonogenic survival curves are measured per cell line and published with
/// their α and β.
///
/// Returns exactly `1.0` at `dose_gy == 0.0`, which is what makes the whole
/// layer inert by default.
///
/// ## Why not fold this into the ROS path
///
/// Because the two channels answer different questions and the interesting
/// result is how they separate. See the module docs.
#[must_use = "the surviving fraction is the function's only output"]
pub fn lq_survival(dose_gy: f64, alpha_per_gy: f64, beta_per_gy2: f64) -> f64 {
    debug_assert!(
        dose_gy >= 0.0 && dose_gy.is_finite(),
        "lq_survival: dose_gy must be finite and non-negative, got {dose_gy}"
    );
    debug_assert!(
        alpha_per_gy >= 0.0 && beta_per_gy2 >= 0.0,
        "lq_survival: alpha and beta must be non-negative, got alpha={alpha_per_gy}, beta={beta_per_gy2}"
    );
    let d = dose_gy.max(0.0);
    (-(alpha_per_gy.max(0.0) * d + beta_per_gy2.max(0.0) * d * d))
        .exp()
        .clamp(0.0, 1.0)
}

/// Dose scaled by the oxygen effect: what `dose_gy` is WORTH at `o2_supply`.
///
/// The oxygen effect is dose-modifying, so a hypoxic cell behaves as though it
/// received less dose. This multiplies by
/// [`crate::oxygen::oer_relative_efficacy`], the same normalised hyperbola the
/// SDT path uses, so the two modalities cannot drift apart in their O₂ shape.
///
/// `dependence ∈ [0, 1]` scales how much of the effect applies, exactly as
/// `crate::oxygen::oer_exo_factor`'s does: at `0.0` this returns `dose_gy`
/// unchanged, which is the O₂-independent limit and is the identity path.
#[must_use = "the scaled dose is the function's only output"]
pub fn oer_scaled_dose(dose_gy: f64, o2_supply: f64, dependence: f64, p_full_mmhg: f64) -> f64 {
    let dep = dependence.clamp(0.0, 1.0);
    let eff = oer_relative_efficacy(o2_supply, p_full_mmhg);
    dose_gy.max(0.0) * (1.0 - dep + dep * eff)
}

/// Relative delivered dose at depth `z_um`, `exp(−μ z)`.
///
/// Near-flat over a tumour: at `MU_6MV_SOFT_TISSUE_PER_CM` a 10 cm depth still
/// receives ~74 % of the surface dose. Compare `crate::physics`'s PDT and SDT
/// attenuation over the same distance.
///
/// Negative depths (outside the tumour) clamp to `0`, matching the sign
/// convention of `crate::grid::TumorGrid3D::radial_depth_um`.
#[must_use = "the intensity fraction is the function's only output"]
pub fn intensity_at_depth(z_um: f64, mu_per_cm: f64) -> f64 {
    debug_assert!(
        mu_per_cm >= 0.0 && mu_per_cm.is_finite(),
        "intensity_at_depth: mu_per_cm must be finite and non-negative, got {mu_per_cm}"
    );
    let z_cm = (z_um.max(0.0)) / 10_000.0;
    (-mu_per_cm.max(0.0) * z_cm).exp().clamp(0.0, 1.0)
}

/// Probability this cell is killed by DNA damage alone, `1 − SF`.
///
/// The channel that does NOT pass through the ferroptosis state. A consumer
/// rolls against it once, at the moment of delivery.
#[must_use = "the lethality is the function's only output"]
pub fn dna_lethality(cfg: &RadiationConfig, o2_supply: f64, z_um: f64) -> f64 {
    let d = delivered_dose(cfg, o2_supply, z_um);
    1.0 - lq_survival(d, cfg.alpha_per_gy, cfg.beta_per_gy2)
}

/// Exogenous-ROS peak this cell receives, for the ferroptosis channel.
///
/// `ros_per_gy` is UNCALIBRATED and the CALIBRATION_STATUS row says so: the
/// radiation-ferroptosis literature establishes the direction (radiation
/// induces lipid peroxidation; GPX4 inhibition radiosensitizes) and reports no
/// conversion from gray to an intracellular ROS peak. At the default
/// `ros_per_gy = 0.0` this returns `0.0`, so the ferroptosis channel is off
/// unless a caller opts in.
#[must_use = "the ROS peak is the function's only output"]
pub fn exo_ros_from_dose(cfg: &RadiationConfig, o2_supply: f64, z_um: f64) -> f64 {
    cfg.ros_per_gy.max(0.0) * delivered_dose(cfg, o2_supply, z_um)
}

/// Dose actually delivered to a cell: prescribed, attenuated, O₂-scaled.
#[must_use = "the delivered dose is the function's only output"]
pub fn delivered_dose(cfg: &RadiationConfig, o2_supply: f64, z_um: f64) -> f64 {
    let attenuated = cfg.dose_gy.max(0.0) * intensity_at_depth(z_um, cfg.mu_per_cm);
    oer_scaled_dose(attenuated, o2_supply, cfg.o2_dependence, cfg.p_full_mmhg)
}

/// The dose-modifying factor the O₂ scaling applies — a RESTATEMENT of the
/// OER, not a prediction of it.
///
/// Named and documented this way deliberately. [`oer_scaled_dose`] multiplies
/// dose by the normalised hyperbola, so the ratio between two oxygenations
/// comes back out as the hyperbola again: at full dependence and anoxia it is
/// `m(p_full)/m(0)`, about 2.86 for `p_full = 40` mmHg. Reporting that number
/// against the published 2.5–3.0 band would be a guard computing its own
/// expectation, and this module refuses to do it.
///
/// **The layer's actual independent prediction lives one level up.** The
/// ferroptosis channel carries its own O₂ dependence through
/// `crate::oxygen::oer_exo_factor`, and the DNA channel carries this one, so
/// the OER a full run EXHIBITS — the dose ratio at equal simulated survival
/// with both channels active — is an emergent quantity that neither input
/// fixes. Fit `ros_per_gy` at full supply, then measure that ratio: it is free
/// to land outside 2.5–3.0, and if it does, this layer is wrong.
///
/// This function exists so the analysis layer can report the DNA channel's
/// contribution separately from the combined result, which is the comparison
/// that makes the emergent number interpretable.
#[must_use = "the factor is the function's only output"]
pub fn dna_channel_dose_modifying_factor(cfg: &RadiationConfig, o2_supply: f64) -> f64 {
    let hypoxic = oer_scaled_dose(1.0, o2_supply, cfg.o2_dependence, cfg.p_full_mmhg);
    let full = oer_scaled_dose(1.0, 1.0, cfg.o2_dependence, cfg.p_full_mmhg);
    if hypoxic <= 0.0 {
        return f64::INFINITY;
    }
    full / hypoxic
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::oxygen::oxygen_enhancement_ratio;
    use crate::params::SpatialParams;
    use crate::physics::{pdt_intensity_at_depth, sdt_intensity_at_depth};

    fn cfg() -> RadiationConfig {
        RadiationConfig {
            dose_gy: 2.0,
            alpha_per_gy: 0.3,
            beta_per_gy2: 0.03,
            ros_per_gy: 0.0,
            o2_dependence: 1.0,
            p_full_mmhg: crate::oxygen::OER_REFERENCE_PO2_MMHG,
            mu_per_cm: MU_6MV_SOFT_TISSUE_PER_CM,
        }
    }

    #[test]
    fn the_default_config_is_the_identity_element() {
        // The byte-identity contract: a run carrying the default must be
        // indistinguishable from a run with no radiation model. Asserted on
        // BIT PATTERNS, not with a tolerance -- `1.0 - 1e-17` would pass an
        // approximate check and still move a committed hash.
        let d = RadiationConfig::default();
        assert_eq!(d.dose_gy.to_bits(), 0.0_f64.to_bits());
        assert_eq!(
            lq_survival(d.dose_gy, d.alpha_per_gy, d.beta_per_gy2).to_bits(),
            1.0_f64.to_bits()
        );
        assert_eq!(dna_lethality(&d, 1.0, 0.0).to_bits(), 0.0_f64.to_bits());
        assert_eq!(exo_ros_from_dose(&d, 1.0, 0.0).to_bits(), 0.0_f64.to_bits());
        // And at every oxygenation and depth, not just the convenient one.
        for &o2 in &[0.0, 0.05, 0.5, 1.0] {
            for &z in &[0.0, 500.0, 10_000.0] {
                assert_eq!(dna_lethality(&d, o2, z).to_bits(), 0.0_f64.to_bits());
                assert_eq!(exo_ros_from_dose(&d, o2, z).to_bits(), 0.0_f64.to_bits());
            }
        }
    }

    #[test]
    fn the_fitted_alpha_predicts_a_published_point_it_was_not_fitted_to() {
        // The whole reason this layer satisfies the layer-freeze policy.
        // PMID 32307022 states THREE numbers; a fit needs two. So alpha comes
        // from SF2 = 0.5445 at alpha/beta = 10, and the 2.67 Gy point is a
        // PREDICTION that the model is free to get wrong.
        let beta = ALPHA_GLIOBLASTOMA_PER_GY / ALPHA_BETA_TUMOUR_GY;

        // The fit reproduces the point it was fitted to.
        let sf2 = lq_survival(2.0, ALPHA_GLIOBLASTOMA_PER_GY, beta);
        assert!(
            (sf2 - 0.5445).abs() < 1e-4,
            "the fitted alpha does not reproduce SF2 = 0.5445: {sf2}"
        );

        // And PREDICTS the one it was not.
        let sf267 = lq_survival(2.67, ALPHA_GLIOBLASTOMA_PER_GY, beta);
        assert!(
            (sf267 - SF_AT_2_67_GY_PUBLISHED).abs() < 1e-3,
            "predicted SF(2.67) = {sf267}, published {SF_AT_2_67_GY_PUBLISHED}"
        );

        // The prediction must not be trivially satisfied by any alpha: a
        // 20%-different alpha has to MISS it, or the test measures nothing.
        let wrong = lq_survival(2.67, ALPHA_GLIOBLASTOMA_PER_GY * 1.2, beta * 1.2);
        assert!(
            (wrong - SF_AT_2_67_GY_PUBLISHED).abs() > 1e-2,
            "a 20% wrong alpha still hits the published point ({wrong}), so \
             this test does not measure the fit"
        );
    }

    #[test]
    fn survival_falls_monotonically_with_dose() {
        let c = cfg();
        let mut prev = 1.0;
        for step in 0..40 {
            let d = step as f64 * 0.5;
            let sf = lq_survival(d, c.alpha_per_gy, c.beta_per_gy2);
            assert!(sf <= prev, "survival rose at {d} Gy: {sf} > {prev}");
            prev = sf;
        }
        assert!(prev < 1e-6, "20 Gy should be near-sterilising, got {prev}");
    }

    #[test]
    fn the_quadratic_term_makes_one_dose_deadlier_than_two_halves() {
        // The signature of the LQ model and the whole reason fractionation
        // works. Pinned in BOTH directions so the test cannot pass with beta
        // ignored: with beta = 0 the model is exponential and splitting the
        // dose changes nothing AT ALL (bit-equal), and with beta > 0 the
        // single dose must strictly win.
        let (a, b) = (0.3, 0.03);
        let d = 6.0;
        let one = lq_survival(d, a, b);
        let two = lq_survival(d / 2.0, a, b).powi(2);
        assert!(
            one < two,
            "beta>0 must make one dose deadlier: {one} vs {two}"
        );

        let one_lin = lq_survival(d, a, 0.0);
        let two_lin = lq_survival(d / 2.0, a, 0.0).powi(2);
        assert!(
            (one_lin - two_lin).abs() < 1e-12,
            "beta=0 must be exponential and fractionation-free: {one_lin} vs {two_lin}"
        );
    }

    #[test]
    fn the_oxygen_scaling_is_off_at_zero_dependence_and_is_the_hyperbola_at_one() {
        // Off must be BIT-identical, because this is the flag that keeps the
        // production matrix stable.
        for &o2 in &[0.0, 0.1, 0.5, 1.0] {
            assert_eq!(
                oer_scaled_dose(2.0, o2, 0.0, 40.0).to_bits(),
                2.0_f64.to_bits(),
                "dependence 0 moved the dose at o2={o2}"
            );
        }
        // On, it must be the SAME hyperbola the SDT path uses -- recomputed
        // here from `oxygen_enhancement_ratio` rather than restated, so a
        // change to that function fails this test instead of silently letting
        // the two modalities drift apart in their O2 shape.
        let m0 = oxygen_enhancement_ratio(0.0);
        let m40 = oxygen_enhancement_ratio(40.0);
        assert!(
            (oer_scaled_dose(2.0, 0.0, 1.0, 40.0) - 2.0 * m0 / m40).abs() < 1e-12,
            "the anoxic scaling is not m(0)/m(40)"
        );
    }

    #[test]
    fn the_dna_channel_factor_is_the_hyperbola_and_says_so() {
        // Asserted for what it IS -- a restatement of the OER -- so nobody
        // later reads it as the layer's prediction. The doc comment on
        // `dna_channel_dose_modifying_factor` makes the same point in prose;
        // this pins it.
        let c = cfg();
        let m0 = oxygen_enhancement_ratio(0.0);
        let m40 = oxygen_enhancement_ratio(c.p_full_mmhg);
        let f = dna_channel_dose_modifying_factor(&c, 0.0);
        assert!(
            (f - m40 / m0).abs() < 1e-12,
            "the factor should be exactly m(p_full)/m(0), got {f}"
        );
        // At full supply it is exactly 1: no modification.
        assert!((dna_channel_dose_modifying_factor(&c, 1.0) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn radiation_reaches_depth_where_light_and_ultrasound_do_not() {
        // The claim the module docs make about the manuscript's tissue-access
        // argument, measured against the two modalities it is compared with
        // rather than asserted. 1 cm is the scale that argument uses.
        let p = SpatialParams::default();
        let z = 10_000.0_f64; // 1 cm
        let rad = intensity_at_depth(z, MU_6MV_SOFT_TISSUE_PER_CM);
        let pdt = pdt_intensity_at_depth(z, &p);
        let sdt = sdt_intensity_at_depth(z, &p);
        assert!(
            rad > 0.95,
            "6 MV should be nearly flat over 1 cm, got {rad}"
        );
        assert!(
            rad > sdt,
            "radiation must outreach SDT at 1 cm: {rad} vs {sdt}"
        );
        assert!(
            rad > pdt,
            "radiation must outreach PDT at 1 cm: {rad} vs {pdt}"
        );
        // And the gap must be a GAP, not a rounding difference. Stated as
        // absolute bounds rather than a ratio: `pdt` can legitimately be zero
        // (the default drug-light interval yields no sensitizer), and a
        // division would then report `inf` as a pass.
        assert!(
            pdt < 0.1,
            "PDT should be nearly extinguished at 1 cm, got {pdt}"
        );
        assert!(sdt < 0.9, "SDT should be attenuated at 1 cm, got {sdt}");
    }

    #[test]
    fn attenuation_is_monotone_bounded_and_clamps_outside_the_tumour() {
        let mut prev = 1.0;
        for cm in 0..30 {
            let v = intensity_at_depth(cm as f64 * 10_000.0, MU_6MV_SOFT_TISSUE_PER_CM);
            assert!(v <= prev && v >= 0.0 && v <= 1.0, "at {cm} cm: {v}");
            prev = v;
        }
        // Negative depth is outside the tumour; the grid's sign convention
        // makes that reachable, and it must read as the surface rather than
        // amplifying the dose.
        assert_eq!(
            intensity_at_depth(-5_000.0, MU_6MV_SOFT_TISSUE_PER_CM).to_bits(),
            1.0_f64.to_bits()
        );
        // mu = 0 is the no-attenuation limit, exactly.
        assert_eq!(
            intensity_at_depth(10_000.0, 0.0).to_bits(),
            1.0_f64.to_bits()
        );
    }

    #[test]
    fn the_two_channels_are_independent() {
        // The design claim: DNA lethality does not depend on `ros_per_gy`, and
        // the ferroptosis channel does not depend on alpha or beta. A single
        // implementation that folded them together would fail this.
        let base = cfg();
        let mut ros_heavy = base;
        ros_heavy.ros_per_gy = 5.0;
        assert_eq!(
            dna_lethality(&base, 0.5, 0.0).to_bits(),
            dna_lethality(&ros_heavy, 0.5, 0.0).to_bits(),
            "ros_per_gy moved the DNA channel"
        );

        let mut lq_heavy = base;
        lq_heavy.alpha_per_gy = 0.9;
        lq_heavy.beta_per_gy2 = 0.5;
        let mut a = base;
        a.ros_per_gy = 2.0;
        let mut b = lq_heavy;
        b.ros_per_gy = 2.0;
        assert_eq!(
            exo_ros_from_dose(&a, 0.5, 0.0).to_bits(),
            exo_ros_from_dose(&b, 0.5, 0.0).to_bits(),
            "alpha/beta moved the ferroptosis channel"
        );
    }

    #[test]
    fn hypoxia_spares_a_cell_in_both_channels() {
        // Both channels must respond to O2 in the same DIRECTION, since both
        // read the same hyperbola. A sign error in either would show here.
        let mut c = cfg();
        c.ros_per_gy = 1.0;
        let (well, poor) = (1.0, 0.02);
        assert!(dna_lethality(&c, well, 0.0) > dna_lethality(&c, poor, 0.0));
        assert!(exo_ros_from_dose(&c, well, 0.0) > exo_ros_from_dose(&c, poor, 0.0));
    }

    #[test]
    fn radiation_through_from_cell_has_no_exogenous_ros() {
        // The CONTRACT `biochem.rs` documents, asserted rather than left as a
        // silent trap: `CellState::from_cell` cannot see dose, depth or local
        // O2, so a Radiation cell built through it carries zero exogenous ROS
        // and the ferroptosis channel is off. A caller wanting that channel
        // goes through `from_cell_with_ros` fed by `exo_ros_from_dose`.
        use crate::biochem::CellState;
        use crate::cell::{gen_cell, Phenotype, Treatment};
        use crate::params::Params;
        use rand::{rngs::StdRng, SeedableRng};

        let params = Params::default();
        let mut rng = StdRng::seed_from_u64(7);
        let cell = gen_cell(Phenotype::Glycolytic, &mut rng);
        let state = CellState::from_cell(&cell, Treatment::Radiation, &params, &mut rng);
        assert_eq!(state.exo_ros_peak.to_bits(), 0.0_f64.to_bits());

        // And the explicit path DOES deliver it, so the contract is a routing
        // rule rather than the channel being broken.
        let mut cfg = cfg();
        cfg.ros_per_gy = 2.0;
        let ros = exo_ros_from_dose(&cfg, 1.0, 0.0);
        assert!(ros > 0.0, "the dose-driven path yields nothing: {ros}");
        let state = CellState::from_cell_with_ros(&cell, Treatment::Radiation, &params, ros);
        assert_eq!(state.exo_ros_peak, ros);
    }

    #[test]
    fn the_depth_arm_is_reachable_through_the_engine_not_only_this_module() {
        // `physics::local_ros_multiplier` matches on the treatment, and the
        // depth claim is about THAT function rather than about
        // `intensity_at_depth` in isolation. Measured through the engine's
        // own entry point, at the 1 cm scale the tissue-access argument uses.
        use crate::cell::Treatment;
        use crate::params::SpatialParams;
        use crate::physics::local_ros_multiplier_3d;

        let p = SpatialParams::default();
        let z = 10_000.0_f64;
        let rad = local_ros_multiplier_3d(z, Treatment::Radiation, &p);
        let sdt = local_ros_multiplier_3d(z, Treatment::SDT, &p);
        let pdt = local_ros_multiplier_3d(z, Treatment::PDT, &p);
        assert!(rad > 0.95, "radiation should be near-flat at 1 cm: {rad}");
        assert!(rad > sdt && rad > pdt, "{rad} vs sdt {sdt} pdt {pdt}");
    }

    #[test]
    #[should_panic(expected = "dose_gy must be finite and non-negative")]
    fn a_negative_dose_is_rejected() {
        lq_survival(-1.0, 0.3, 0.03);
    }
}
