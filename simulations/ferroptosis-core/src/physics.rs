//! Energy deposition models for PDT, SDT, and RSL3.
//!
//! Physics of how treatment energy attenuates with tissue depth,
//! converting to local ROS dose for each cell in the spatial model.

use crate::cell::Treatment;
use crate::params::SpatialParams;

/// PDT: Modified Beer-Lambert law for light in tissue, scaled by the
/// per-photon ROS yield of the photosensitizer at the drug-light interval.
///
/// I_eff(z, t_DLI) = I₀ × exp(-µ_eff × z) × Y_drug(t_DLI)
///
/// where `Y_drug(t)` is `Photosensitizer::yield_at(t)` —
/// `concentration_at(t) × phi_so2_relative` for `Porfimer`, or just
/// `concentration_at(t)` for `Uniform`. The `phi_so2_relative` factor
/// closes the inter-drug yield-comparison gap; `Params::pdt_ros` is
/// calibrated to porfimer at peak (yield = 1.0), and other drug variants
/// would set their `phi_so2_relative` to absolute_phi_so2 / porfimer's
/// (~0.65) so the calibration carries through.
///
/// µ_eff = sqrt(3 × µ_a × (µ_a + µ_s'))
///
/// At 630nm red light: δ = 1/µ_eff ≈ 3.2mm
/// At 660nm red light: δ ≈ 4-10mm
///
/// The default `Photosensitizer::Uniform(1.0)` with DLI=0 returns exactly
/// 1.0, preserving pre-PK behavior bit-for-bit.
///
/// Ref: Jacques SL, "Optical properties of biological tissues: a review",
///      Phys Med Biol 58(11):R37-61, 2013
/// Ref: Bellnier DA et al., Lasers Surg Med 38(5):439-444, 2006 — clinical
///      photosensitizer PK (porfimer, Photochlor, 5-ALA-PpIX).
/// Ref: Wilson BC, Patterson MS, "The physics, biophysics and technology
///      of photodynamic therapy", Phys Med Biol 53(9):R61-109, 2008 —
///      porfimer absolute phi_so2 ≈ 0.65 in solution.
///
/// Returns a non-negative intensity multiplier relative to surface.
///
/// For valid `Photosensitizer::Uniform(c)` with `c ≤ 1.0` or any valid
/// `Photosensitizer::Porfimer` with `phi_so2_relative ≤ 1`, the value
/// stays in `[0, 1]`. `Uniform(c)` with `c > 1.0` and
/// `phi_so2_relative > 1` are intentionally permitted (forward-compat
/// hooks for enrichment / sensitizer-engineered variants) and can
/// produce multipliers above 1.0.
///
/// Invalid configurations (NaN, negative `phi_so2_relative` or
/// `t_distribution_h`, non-positive `t_half_h`) trigger `debug_assert!`
/// failures in test/debug builds. In release builds those asserts are
/// compiled out and outputs are not bounded — call
/// [`Photosensitizer::validate`] explicitly when loading from untrusted
/// sources.
///
/// [`Photosensitizer::validate`]: crate::photosensitizer_pk::Photosensitizer::validate
pub fn pdt_intensity_at_depth(z_um: f64, params: &SpatialParams) -> f64 {
    let z_mm = z_um / 1000.0;
    let drug_yield = params
        .photosensitizer
        .yield_at(params.t_drug_light_interval_h);
    params.pdt_i0 * (-params.pdt_mu_eff * z_mm).exp() * drug_yield
}

/// SDT: Ultrasound attenuation in soft tissue.
///
/// I(z) = I₀ × 10^(-α × f × z / 10)
///
/// α = attenuation coefficient (dB/cm/MHz)
///   Soft tissue average: 0.7
///   Muscle: 1.3, Fat: 0.6, Liver: 0.45, Blood: 0.18
/// f = frequency (MHz)
/// z = depth (cm)
///
/// At 1 MHz in soft tissue:
///   1 cm: 85% intensity
///   5 cm: 45% intensity
///   10 cm: 20% intensity
///
/// Ref: Cobbold RSC, "Foundations of Biomedical Ultrasound", 2007
///      Christensen DA, "Ultrasonic Bioinstrumentation"
///
/// Returns intensity multiplier in [0, 1] relative to surface.
pub fn sdt_intensity_at_depth(z_um: f64, params: &SpatialParams) -> f64 {
    let z_cm = z_um / 10_000.0;
    let attenuation_db = params.sdt_alpha * params.sdt_freq_mhz * z_cm;
    params.sdt_i0 * 10.0_f64.powf(-attenuation_db / 10.0)
}

/// RSL3: Systemic drug — uniform concentration throughout tissue.
/// No depth dependence.
pub fn rsl3_concentration(_z_um: f64) -> f64 {
    1.0
}

/// Compute the ROS dose multiplier for a cell at a given row in the grid.
/// Energy is applied from the top (row 0 = tissue surface).
///
/// Returns a multiplier in [0, 1] that scales the base exogenous ROS.
/// For Control, always returns 0.
pub fn local_ros_multiplier(
    row: usize,
    cell_size_um: f64,
    tx: Treatment,
    params: &SpatialParams,
) -> f64 {
    let z_um = row as f64 * cell_size_um;
    match tx {
        Treatment::SDT => sdt_intensity_at_depth(z_um, params),
        Treatment::PDT => pdt_intensity_at_depth(z_um, params),
        Treatment::RSL3 => rsl3_concentration(z_um),
        Treatment::Control => 0.0,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::SpatialParams;

    #[test]
    fn sdt_1cm_about_85_percent() {
        let params = SpatialParams::default();
        let intensity = sdt_intensity_at_depth(10_000.0, &params);
        // At 1 MHz, 0.7 dB/cm: attenuation = 0.7 dB → 10^(-0.07) ≈ 0.851
        assert!(
            (intensity - 0.851).abs() < 0.01,
            "SDT at 1cm: expected ~0.851, got {intensity}"
        );
    }

    #[test]
    fn sdt_5cm_about_45_percent() {
        let params = SpatialParams::default();
        let intensity = sdt_intensity_at_depth(50_000.0, &params);
        // At 1 MHz, 0.7 dB/cm: attenuation = 3.5 dB → 10^(-0.35) ≈ 0.447
        assert!(
            (intensity - 0.447).abs() < 0.02,
            "SDT at 5cm: expected ~0.447, got {intensity}"
        );
    }

    #[test]
    fn pdt_at_delta_is_1_over_e() {
        let params = SpatialParams {
            pdt_mu_eff: 0.31,
            pdt_i0: 1.0,
            ..Default::default()
        };
        // δ = 1/0.31 ≈ 3.23 mm = 3226 µm
        let delta_um = 1000.0 / params.pdt_mu_eff;
        let intensity = pdt_intensity_at_depth(delta_um, &params);
        let expected = 1.0 / std::f64::consts::E;
        assert!(
            (intensity - expected).abs() < 0.01,
            "PDT at δ: expected {expected:.4}, got {intensity:.4}"
        );
    }

    #[test]
    fn pdt_much_deeper_is_negligible() {
        let params = SpatialParams::default();
        // At 20mm (2cm), PDT should be essentially zero
        let intensity = pdt_intensity_at_depth(20_000.0, &params);
        assert!(
            intensity < 0.01,
            "PDT at 20mm should be negligible, got {intensity}"
        );
    }

    #[test]
    fn rsl3_is_uniform() {
        assert_eq!(rsl3_concentration(0.0), rsl3_concentration(100_000.0));
    }

    #[test]
    fn pdt_with_porfimer_phi_half_is_half_of_default() {
        // Locks down the phi-scaling path in pdt_intensity_at_depth.
        // If anyone refactors the function to bypass yield_at, the
        // unit-tested phi behavior on Photosensitizer is still right
        // but the wired-through physics path would silently regress.
        // This test catches that.
        use crate::photosensitizer_pk::Photosensitizer;

        let z_um = 1000.0; // arbitrary depth
        let baseline = pdt_intensity_at_depth(z_um, &SpatialParams::default());

        let phi_half_params = SpatialParams {
            photosensitizer: Photosensitizer::Porfimer {
                t_half_h: 504.0,
                t_distribution_h: 0.0,
                phi_so2_relative: 0.5,
            },
            t_drug_light_interval_h: 0.0, // peak — only phi affects the result
            ..Default::default()
        };
        let phi_half = pdt_intensity_at_depth(z_um, &phi_half_params);

        // baseline * 0.5 is IEEE-exact (multiplying any finite f64 by 0.5
        // is exact since 0.5 has an exact binary representation), so
        // strict equality is safe here — distinct from the libm-dependent
        // exp(-ln(2)) used in the half-life test below.
        assert_eq!(phi_half, baseline * 0.5);
    }

    #[test]
    fn pdt_with_porfimer_at_one_halflife_is_half_of_default() {
        use crate::photosensitizer_pk::Photosensitizer;

        let z_um = 1000.0; // arbitrary depth
        let baseline_params = SpatialParams::default();
        let baseline = pdt_intensity_at_depth(z_um, &baseline_params);

        let pk_params = SpatialParams {
            photosensitizer: Photosensitizer::Porfimer { t_half_h: 504.0, t_distribution_h: 0.0, phi_so2_relative: 1.0 },
            t_drug_light_interval_h: 504.0,
            ..Default::default()
        };
        let with_pk = pdt_intensity_at_depth(z_um, &pk_params);

        // `exp(-ln(2))` lands on 0.5 on most libms but is not guaranteed
        // bit-exact across platforms; use a tight relative tolerance.
        let expected = baseline * 0.5;
        assert!(
            (with_pk - expected).abs() < 1e-12,
            "with_pk = {with_pk}, expected ~{expected}"
        );
    }

    #[test]
    fn pdt_with_default_photosensitizer_unchanged() {
        // The default Uniform(1.0) + DLI=0 must not perturb existing physics.
        let params = SpatialParams::default();
        let z_um = 3226.0; // 1/µ_eff at default 0.31 /mm
        let intensity = pdt_intensity_at_depth(z_um, &params);
        let expected = 1.0 / std::f64::consts::E;
        assert!(
            (intensity - expected).abs() < 0.01,
            "expected ~{expected:.4}, got {intensity:.4}"
        );
    }

    #[test]
    fn sdt_much_deeper_still_significant() {
        let params = SpatialParams::default();
        // At 10cm, SDT should still have meaningful intensity
        let intensity = sdt_intensity_at_depth(100_000.0, &params);
        // 0.7 * 1.0 * 10 = 7 dB → 10^(-0.7) ≈ 0.200
        assert!(
            intensity > 0.15,
            "SDT at 10cm should still be significant, got {intensity}"
        );
    }
}
