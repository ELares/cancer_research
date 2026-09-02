//! Sonodynamic therapy: what makes it a different modality from PDT (#831).
//!
//! WHY THIS MODULE DID NOT EXIST FOR SO LONG, WHICH IS THE INTERESTING PART
//! ------------------------------------------------------------------------
//! Sonodynamic therapy is one of the two modalities this project was built
//! around, and until this module it owned **zero lines** of the engine. That
//! is not an oversight. `analysis/calibration/pdt-threshold-validation.md`
//! records a deliberate, externally-validated design choice: the engine sets
//! `sdt_ros == pdt_ros` and feeds both into **one** `death_threshold`, because
//! Zhu 2015 measured the in-vivo reacted-singlet-oxygen necrosis threshold to
//! be approximately photosensitizer-independent (~0.5 mM within a factor of
//! 1.8 across three chemically unrelated sensitizers). The dose-to-kill is a
//! property of the target, not of the source that made the ROS.
//!
//! `scripts/engine_selectivity.py` then measured the consequence and reported
//! it as a defect: SDT and PDT are **bit-identical** across every phenotype,
//! so a single-cell contrast between them compares a modality with itself.
//!
//! Both statements are correct, and together they say where the difference
//! has to live. If the kill threshold is genuinely source-independent, then
//! every difference between these two arms is **upstream of the cell** — in
//! how much reacted ROS arrives, not in what it does on arrival. The engine
//! already carried one such difference (depth attenuation, in
//! [`crate::physics`]). This module carries the other, and it is the one that
//! makes SDT structurally unlike anything else in the chapter.
//!
//! THE MECHANICAL INDEX IS A THRESHOLD, NOT A GRADIENT
//! ---------------------------------------------------
//! Light delivers energy continuously: half the photons make half the singlet
//! oxygen. Ultrasound at therapeutic amplitude does not, because the dominant
//! sonochemical ROS route is **inertial cavitation** — a bubble grows over a
//! rarefaction half-cycle and collapses violently, and it either happens or
//! it does not. Apfel & Holland (1991, PMID 2053214) showed the likelihood is
//! governed by the peak rarefactional pressure divided by the square root of
//! the frequency, and that ratio is the mechanical index every clinical
//! scanner displays.
//!
//! So SDT joins the small set of arms in this chapter whose limit is a
//! THRESHOLD rather than a multiplier — CAR-T's antigen density, the ADC's
//! receptor, the oncolytic entry receptor. Below the cavitation threshold no
//! insonation time produces sonochemical ROS at all, which is a different
//! kind of failure from "not enough dose", and it is why the SDT depth limit
//! is a hard edge while the PDT one is an exponential fade.
//!
//! WHAT THIS MODULE PREDICTS AND WHERE IT IS ALREADY CONTRADICTED
//! -------------------------------------------------------------
//! Three frequency dependences act at once and two of them oppose:
//!
//!   * focusing  — the diffraction-limited focal width scales with the
//!     wavelength, so at fixed acoustic power a higher frequency concentrates
//!     the beam into a smaller spot and RAISES focal pressure
//!   * attenuation — soft tissue absorbs in proportion to frequency, so a
//!     higher frequency LOSES more on the way in, and the loss compounds with
//!     depth
//!   * the index itself — dividing by the square root of frequency penalises
//!     high frequency a third time
//!
//! The product has an interior maximum: [`optimal_frequency_mhz`]. Its closed
//! form is `f* = 10 / (alpha * ln(10) * z)`, so the model predicts the best
//! frequency falls as **1/depth** and as **1/attenuation**.
//!
//! Ellens & Hynynen (2015, Med Phys 42(8):4896, PMID 26233216) simulated the
//! same three mechanisms independently, by full-wave Rayleigh-Sommerfeld
//! propagation into a bioheat solve, and named them in the same terms —
//! "focal size (decreasing with frequency), peak pressure (generally
//! increasing with frequency), and attenuation (also increasing with
//! frequency)". They report an interior optimum, and it moves DOWN when
//! attenuation doubles (750 kHz to 500 kHz at 100 and 150 mm). Both match.
//!
//! **The depth scaling does not.** They report the SAME 750 kHz optimum at
//! 50, 100 and 150 mm, and this model says it should fall threefold across
//! that range. `scripts/validate_sonodynamic.py` treats that as a refutation
//! rather than a tuning problem, and the missing term is nameable: their
//! efficiency is limited by NEAR-FIELD heating, which they state is the rate
//! limiter for large-volume ablation, and this model has no term for the
//! tissue in front of the focus at all. Chasing the focal optimum downward
//! also lengthens the heated path, and a model that cannot see that path is
//! free to slide its optimum in a way a real applicator is not.

/// Peak rarefactional pressure divided by the square root of frequency:
/// the mechanical index (Apfel & Holland 1991, PMID 2053214).
///
/// `p_neg_mpa` in megapascals, `freq_mhz` in megahertz, result dimensionless
/// by convention (the unit MPa·MHz^-1/2 is dropped, as clinically displayed).
///
/// Zero or negative frequency returns 0.0 rather than diverging: a "wave"
/// with no frequency delivers no rarefaction half-cycle, so no bubble grows.
pub fn mechanical_index(p_neg_mpa: f64, freq_mhz: f64) -> f64 {
    if !(freq_mhz > 0.0) || !p_neg_mpa.is_finite() {
        return 0.0;
    }
    p_neg_mpa.max(0.0) / freq_mhz.sqrt()
}

/// The FDA diagnostic-imaging cap on the mechanical index.
///
/// A REGULATORY LIMIT, NOT A MEASUREMENT, and it is here as the ceiling a
/// clinical device is allowed to reach rather than as the pressure at which
/// tissue cavitates. Therapeutic ultrasound is exempt from it; a diagnostic
/// scanner is not. Nothing in this module treats it as a threshold.
pub const MI_DIAGNOSTIC_CAP: f64 = 1.9;

/// Whether inertial cavitation is expected, given an index and a threshold.
///
/// The threshold is a PARAMETER and deliberately not a constant. Published
/// in-vivo inertial-cavitation thresholds vary by more than an order of
/// magnitude with nucleation: seeded with a microbubble contrast agent or a
/// cavitation nucleus they fall well below the diagnostic cap, and in
/// nucleus-free tissue they sit far above it. A single number here would
/// assert a nucleation state the model does not represent.
pub fn cavitates(index: f64, threshold: f64) -> bool {
    index >= threshold && threshold > 0.0
}

/// Pressure amplitude surviving to `depth_cm`, as a fraction of the surface
/// amplitude.
///
/// DERIVED FROM THE ENGINE'S OWN INTENSITY LAW rather than restated: the
/// square root of [`crate::physics::sdt_intensity_at_depth`]'s
/// `10^(-alpha*f*z/10)`, because pressure amplitude goes as the square root
/// of intensity. Writing the exponent out again with a hand-chosen divisor
/// is exactly how two attenuation conventions end up in one crate, and the
/// factor-of-two between an amplitude dB and an intensity dB is the easiest
/// unit error in acoustics to make.
pub fn pressure_fraction_at_depth(depth_cm: f64, freq_mhz: f64, alpha_db_cm_mhz: f64) -> f64 {
    let db = alpha_db_cm_mhz * freq_mhz * depth_cm.max(0.0);
    10.0_f64.powf(-db / 20.0)
}

/// Focal pressure amplitude at fixed acoustic power, relative to a 1 MHz
/// reference.
///
/// A diffraction-limited focus has width proportional to the wavelength, so
/// its area goes as `1/f^2` and the intensity there as `f^2` at fixed power;
/// amplitude is the square root, so it is LINEAR in frequency. This is the
/// only one of the three frequency terms that favours going higher, and it
/// is why an optimum exists at all rather than the answer being "as low as
/// the transducer allows".
pub fn focal_pressure_gain(freq_mhz: f64) -> f64 {
    freq_mhz.max(0.0)
}

/// Mechanical index delivered to a focus at `depth_cm`, at fixed acoustic
/// power, relative to a 1 MHz surface reference.
///
/// The product of the three frequency terms, which is the whole model:
/// `focal gain * attenuation / sqrt(f)`.
pub fn delivered_index(depth_cm: f64, freq_mhz: f64, alpha_db_cm_mhz: f64) -> f64 {
    if !(freq_mhz > 0.0) {
        return 0.0;
    }
    let p = focal_pressure_gain(freq_mhz)
        * pressure_fraction_at_depth(depth_cm, freq_mhz, alpha_db_cm_mhz);
    mechanical_index(p, freq_mhz)
}

/// The frequency maximising [`delivered_index`], in closed form.
///
/// Setting `d/df [ 0.5*ln f - (alpha*z/20)*ln 10 * f ] = 0` gives
/// `f* = 10 / (alpha * ln(10) * z)`.
///
/// THE SCALE RIDES ENTIRELY ON `alpha`, so this predicts a SHAPE and not a
/// number: `f* * alpha * z` is the invariant, and quoting a megahertz value
/// without its attenuation coefficient states an assumption as a result.
/// A zero depth has no interior optimum (nothing attenuates), so it returns
/// infinity rather than dividing by zero — the honest answer, and one the
/// scanning routine has to handle rather than pretend away.
pub fn optimal_frequency_mhz(depth_cm: f64, alpha_db_cm_mhz: f64) -> f64 {
    let denom = alpha_db_cm_mhz * std::f64::consts::LN_10 * depth_cm;
    if !(denom > 0.0) {
        return f64::INFINITY;
    }
    10.0 / denom
}

/// Deepest focus at which a given applicator still clears the cavitation
/// threshold, scanning frequency freely at each depth.
///
/// This is the number the threshold makes meaningful and a fraction cannot:
/// a coverage fraction says how much of a margin was treated, while this
/// says there is a depth beyond which NO choice of frequency helps. Returns
/// 0.0 when even a surface focus fails.
///
/// `max_index_at_surface` is the index the applicator delivers at 1 MHz with
/// no attenuation, i.e. the calibration point everything else is relative to.
pub fn cavitation_depth_limit_cm(
    max_index_at_surface: f64,
    threshold: f64,
    alpha_db_cm_mhz: f64,
    max_depth_cm: f64,
    step_cm: f64,
) -> f64 {
    if !(step_cm > 0.0) || !(threshold > 0.0) {
        return 0.0;
    }
    let mut deepest = 0.0;
    let mut z = 0.0;
    while z <= max_depth_cm {
        let f = optimal_frequency_mhz(z, alpha_db_cm_mhz);
        // At z = 0 the optimum is unbounded; the delivered index there is
        // just the surface value, since nothing has attenuated yet.
        let idx = if f.is_finite() {
            max_index_at_surface * delivered_index(z, f, alpha_db_cm_mhz)
        } else {
            max_index_at_surface
        };
        if cavitates(idx, threshold) {
            deepest = z;
        } else if z > 0.0 {
            break;
        }
        z += step_cm;
    }
    deepest
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_index_is_pressure_over_root_frequency() {
        // Apfel & Holland's form, at a point easy to check by hand.
        assert!((mechanical_index(2.0, 4.0) - 1.0).abs() < 1e-12);
        // Doubling the frequency at fixed pressure divides the index by
        // sqrt(2). A mutation using f rather than sqrt(f) fails here.
        let a = mechanical_index(1.0, 1.0);
        let b = mechanical_index(1.0, 2.0);
        assert!((a / b - std::f64::consts::SQRT_2).abs() < 1e-12);
    }

    #[test]
    fn a_zero_frequency_wave_delivers_no_index() {
        assert_eq!(mechanical_index(5.0, 0.0), 0.0);
        assert_eq!(mechanical_index(5.0, -1.0), 0.0);
    }

    #[test]
    fn the_threshold_is_a_threshold_and_not_a_multiplier() {
        // The defining property, and the one that separates this arm from
        // every dose-response arm in the chapter: below the line the answer
        // is the same whatever the index is, so more of it buys nothing.
        assert!(!cavitates(0.69, 0.7));
        assert!(!cavitates(0.01, 0.7));
        assert!(cavitates(0.7, 0.7));
        // A zero threshold would make everything cavitate including silence.
        assert!(!cavitates(0.0, 0.0));
    }

    #[test]
    fn pressure_attenuates_as_the_root_of_the_engines_own_intensity_law() {
        // The guard against a factor-of-two unit error: this module's
        // amplitude law must be the square root of physics.rs's intensity
        // law at the same alpha, frequency and depth, with nothing else
        // adjusted. A mutation changing the /20 to /10 fails.
        use crate::params::SpatialParams;
        let mut sp = SpatialParams::default();
        sp.sdt_i0 = 1.0;
        sp.sdt_freq_mhz = 1.3;
        let z_cm = 4.0;
        let intensity = crate::physics::sdt_intensity_at_depth(z_cm * 10_000.0, &sp);
        let amplitude = pressure_fraction_at_depth(z_cm, sp.sdt_freq_mhz, sp.sdt_alpha);
        assert!(
            (amplitude - intensity.sqrt()).abs() < 1e-12,
            "amplitude {amplitude} is not sqrt(intensity {intensity})"
        );
    }

    #[test]
    fn the_optimum_is_interior_and_the_closed_form_finds_it() {
        // The closed form must agree with a brute scan, which is what makes
        // it a derivation rather than a claim.
        for &(z, a) in &[(3.0, 0.7), (10.0, 0.7), (5.0, 0.3), (12.0, 1.2)] {
            let closed = optimal_frequency_mhz(z, a);
            let mut best = (0.0, 0.0);
            let mut f = 0.001;
            while f < 30.0 {
                let v = delivered_index(z, f, a);
                if v > best.0 {
                    best = (v, f);
                }
                f += 0.001;
            }
            assert!(
                (closed - best.1).abs() < 0.005,
                "closed form {closed} MHz disagrees with scan {} MHz at z={z} alpha={a}",
                best.1
            );
            // Interior means strictly better than both ends, which is the
            // property a monotonic model would not have.
            assert!(delivered_index(z, closed, a) > delivered_index(z, 0.05, a));
            assert!(delivered_index(z, closed, a) > delivered_index(z, 20.0, a));
        }
    }

    #[test]
    fn the_optimum_falls_with_both_depth_and_attenuation() {
        // The two directions Ellens & Hynynen 2015 (PMID 26233216) report.
        // The attenuation one they confirm; the depth one they contradict,
        // and validate_sonodynamic.py is where that is scored. This test
        // pins what the MODEL says, so the disagreement cannot be quietly
        // tuned away.
        assert!(optimal_frequency_mhz(10.0, 0.7) < optimal_frequency_mhz(5.0, 0.7));
        assert!(optimal_frequency_mhz(5.0, 1.4) < optimal_frequency_mhz(5.0, 0.7));
        // Doubling attenuation halves the optimum exactly, which is a
        // sharper claim than "it goes down" and can fail on its own.
        let f1 = optimal_frequency_mhz(8.0, 0.5);
        let f2 = optimal_frequency_mhz(8.0, 1.0);
        assert!((f1 / f2 - 2.0).abs() < 1e-9);
    }

    #[test]
    fn a_zero_depth_has_no_interior_optimum() {
        // Nothing attenuates, so higher is always better and the answer is
        // unbounded. Returning a finite frequency here would invent a limit.
        assert!(optimal_frequency_mhz(0.0, 0.7).is_infinite());
        assert!(optimal_frequency_mhz(5.0, 0.0).is_infinite());
    }

    #[test]
    fn the_depth_limit_is_a_hard_edge_and_zero_when_the_surface_fails() {
        // An applicator too weak to cavitate at the surface fails
        // EVERYWHERE, and must report 0.0 rather than the scan limit --
        // the degenerate-row defect the oncolytic and ablation sections
        // both had to fix.
        assert_eq!(cavitation_depth_limit_cm(0.1, 0.7, 0.7, 20.0, 0.1), 0.0);
        // A strong applicator reaches some finite depth and not beyond.
        let d = cavitation_depth_limit_cm(3.0, 0.7, 0.7, 40.0, 0.1);
        assert!(d > 0.0 && d < 40.0, "expected an interior limit, got {d}");
        // More attenuation cannot reach deeper.
        let shallow = cavitation_depth_limit_cm(3.0, 0.7, 1.4, 40.0, 0.1);
        assert!(shallow <= d);
        // A higher threshold cannot reach deeper either.
        assert!(cavitation_depth_limit_cm(3.0, 1.4, 0.7, 40.0, 0.1) <= d);
    }
}
