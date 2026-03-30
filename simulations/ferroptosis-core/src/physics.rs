//! Energy deposition models for PDT, SDT, and RSL3.
//!
//! Physics of how treatment energy attenuates with tissue depth,
//! converting to local ROS dose for each cell in the spatial model.

use crate::cell::Treatment;
use crate::params::SpatialParams;

/// PDT: Modified Beer-Lambert law for light in tissue.
///
/// I(z) = I₀ × exp(-µ_eff × z)
///
/// µ_eff = sqrt(3 × µ_a × (µ_a + µ_s'))
///
/// At 630nm red light: δ = 1/µ_eff ≈ 3.2mm
/// At 660nm red light: δ ≈ 4-10mm
///
/// Ref: Jacques SL, "Optical properties of biological tissues: a review",
///      Phys Med Biol 58(11):R37-61, 2013
///
/// Returns intensity multiplier in [0, 1] relative to surface.
pub fn pdt_intensity_at_depth(z_um: f64, params: &SpatialParams) -> f64 {
    let z_mm = z_um / 1000.0;
    params.pdt_i0 * (-params.pdt_mu_eff * z_mm).exp()
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
