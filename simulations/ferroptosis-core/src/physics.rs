//! Energy deposition models for PDT, SDT, and RSL3.
//!
//! Physics of how treatment energy attenuates with tissue depth,
//! converting to local ROS dose for each cell in the spatial model.

use crate::cell::Treatment;
use crate::params::SpatialParams;

/// PDT: Modified Beer-Lambert law for light in tissue, scaled by the
/// per-photon ROS yield of the photosensitizer at the drug-light interval.
///
/// **Validity precondition:** `params.photosensitizer` must pass
/// [`Photosensitizer::validate`] and `params.t_drug_light_interval_h`
/// must pass [`validate_dli_h`] (DLI is a free-function validator since
/// it's a plain f64, not a method on `Photosensitizer`). Invalid configs
/// trigger `debug_assert!` in tests but are not bounded in release — see
/// the bottom of this docstring. Untrusted-source callers should
/// validate both before calling.
///
/// [`validate_dli_h`]: crate::photosensitizer_pk::validate_dli_h
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

/// 3D analog of [`local_ros_multiplier`] for spheroid energy deposition.
///
/// Takes a **raw signed radial depth** in micrometers (positive inside the
/// spheroid, negative outside) — typically produced by
/// [`crate::grid::TumorGrid3D::radial_depth_um`]. Negative depths are
/// clipped to zero (cells outside the spheroid are treated as if at the
/// surface), then dispatched to the same depth-attenuation functions
/// used by the 2D path.
///
/// **Why the signature differs from the 2D version:** the 2D model lays
/// energy onto a planar surface (row 0), so `local_ros_multiplier` takes
/// `(row, cell_size)` and computes `z = row × cell_size` internally. The
/// 3D model treats energy as entering isotropically through the spheroid
/// surface, so per-cell depth is **not** a simple coordinate × size — it
/// depends on geometry. Computing depth lives in `TumorGrid3D`; this
/// function takes the already-derived depth so physics stays decoupled
/// from grid representation.
///
/// **2D ≡ 3D at matched depth:** for any treatment and parameter set, if
/// the 3D caller passes `depth = row × cell_size_um` (the same value the
/// 2D path would compute internally), this function returns the same
/// multiplier as `local_ros_multiplier(row, cell_size_um, ...)`. Locked
/// down by `pdt_2d_3d_match_at_same_depth_*` tests below.
///
/// Returns a multiplier in `[0, ∞)`. For `Control`, always returns 0.
pub fn local_ros_multiplier_3d(radial_depth_um: f64, tx: Treatment, params: &SpatialParams) -> f64 {
    // Outside-spheroid cells (negative depth) get the surface value.
    // Physically: energy enters at the surface; we don't model
    // attenuation through stromal tissue outside the tumor (yet).
    let z_um = radial_depth_um.max(0.0);
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

    /// Matched-depth invariant: when the 3D caller passes the same depth
    /// the 2D path computes internally (`row × cell_size_um`), the two
    /// dispatchers must return bit-identical values for **every** treatment
    /// variant. Locks down the contract that 3D physics is not a separate
    /// implementation — it's the same depth-functions reached through a
    /// different per-cell depth derivation.
    ///
    /// Iterates all four `Treatment` variants (including `Control` for the
    /// trivial 0.0 case and `RSL3` for the depth-independent 1.0 case) so a
    /// future refactor that introduced 3D-only physics quirks would fail
    /// here.
    #[test]
    fn local_ros_multiplier_2d_3d_match_at_same_depth() {
        let params = SpatialParams::default();
        let cell_size_um = 20.0_f64;
        // Probe several depths including row 0 (surface) and deeper.
        for row in [0usize, 1, 5, 25, 100] {
            // Compute the depth the 2D path uses internally, then pass the
            // SAME f64 to the 3D path. Any other phrasing (e.g.
            // `(row * cell_size_um) as f64` via mixed types) could
            // introduce 1-ULP drift and weaken the invariant.
            let z_um = row as f64 * cell_size_um;
            for tx in [
                Treatment::Control,
                Treatment::RSL3,
                Treatment::SDT,
                Treatment::PDT,
            ] {
                let v2d = local_ros_multiplier(row, cell_size_um, tx, &params);
                let v3d = local_ros_multiplier_3d(z_um, tx, &params);
                assert_eq!(
                    v2d, v3d,
                    "2D and 3D dispatchers diverged at row={row}, tx={tx:?}, z={z_um}: \
                     2D = {v2d}, 3D = {v3d}"
                );
            }
        }
    }

    /// Negative radial depth (cells outside the spheroid) must be clipped
    /// to zero — i.e. behave like surface intensity, not produce
    /// `exp(positive)` blowup or `10^positive` amplification.
    ///
    /// IEEE-exact: `(-z).max(0.0) == 0.0` exactly, then `exp(0) == 1.0`
    /// and `10^0 == 1.0` exactly, so the multiplier equals `I₀`. The
    /// default `SpatialParams` has `pdt_i0 = sdt_i0 = 1.0` so the
    /// surface multiplier is exactly 1.0 for PDT/SDT.
    #[test]
    fn local_ros_multiplier_3d_negative_depth_clips_to_surface() {
        let params = SpatialParams::default();
        for &z_negative in &[-1.0, -100.0, -10_000.0] {
            for tx in [Treatment::SDT, Treatment::PDT] {
                let at_surface = local_ros_multiplier_3d(0.0, tx, &params);
                let at_negative = local_ros_multiplier_3d(z_negative, tx, &params);
                assert_eq!(
                    at_negative, at_surface,
                    "negative depth z={z_negative} for {tx:?} should clip to surface value"
                );
            }
            // RSL3 is depth-independent; Control is identically 0.
            assert_eq!(
                local_ros_multiplier_3d(z_negative, Treatment::RSL3, &params),
                1.0
            );
            assert_eq!(
                local_ros_multiplier_3d(z_negative, Treatment::Control, &params),
                0.0
            );
        }
    }

    /// The photosensitizer PK path composes correctly when reached via the
    /// 3D dispatcher. Sanity check that the 0.5-yield phi factor produces
    /// exactly half the PDT multiplier the default Uniform(1.0) does, at
    /// the same depth — same property `pdt_with_porfimer_phi_half_is_half`
    /// asserts via the 2D-style direct call, but routed through
    /// `local_ros_multiplier_3d` so a future split in the PDT path would
    /// be caught.
    #[test]
    fn local_ros_multiplier_3d_composes_with_photosensitizer() {
        use crate::photosensitizer_pk::Photosensitizer;

        let baseline_params = SpatialParams::default();
        let half_params = SpatialParams {
            photosensitizer: Photosensitizer::Porfimer {
                t_half_h: 504.0,
                t_distribution_h: 0.0,
                phi_so2_relative: 0.5,
            },
            t_drug_light_interval_h: 0.0,
            ..Default::default()
        };

        let depth_um = 1000.0_f64; // 1 mm — well inside spheroid range
        let baseline = local_ros_multiplier_3d(depth_um, Treatment::PDT, &baseline_params);
        let half = local_ros_multiplier_3d(depth_um, Treatment::PDT, &half_params);

        // 0.5 is IEEE-exact, so strict equality is safe.
        assert_eq!(half, baseline * 0.5);
    }
}
