//! Drug-tolerant persister cells (#241).
//!
//! When a ferroptosis inducer (RSL3, SDT, PDT) is applied to a tumor, a
//! subpopulation of cells that survive the initial insult enters a
//! **drug-tolerant persister** state: slow-cycling, transiently
//! ferroptosis-resistant via *epigenetic* (not genetic) reprogramming. The
//! state is acquired under drug exposure and reverts over time once the drug
//! clears, which distinguishes it from a fixed genetic mutation.
//!
//! Two mechanistic axes are modeled, both well documented in the persister
//! literature:
//!   1. **Reduced lipid-peroxide vulnerability** via partial resistance to the
//!      drug's GPX4 inactivation (NANOG/SOX2-driven dedifferentiation raises
//!      the antioxidant set-point; Hangauer et al. Nature 2017 showed the
//!      persister state is GPX4-dependent and exquisitely ferroptosis-sensitive
//!      *only when GPX4 is also inhibited*).
//!   2. **Increased MUFA membrane remodeling** above the `Params::invivo`
//!      baseline (SCD1-driven monounsaturated-fatty-acid enrichment dilutes
//!      oxidizable PUFAs; Tsoi et al. Cancer Cell 2018 linked the melanoma
//!      persister/de-differentiated state to lipid-metabolism rewiring).
//!
//! Reversion is exponential: the high-mesenchymal, ZEB1-associated vulnerable
//! state re-emerges as the epigenetic mark decays after drug clearance
//! (Viswanathan et al. Nature 2017 established the mesenchymal-state ⇄
//! ferroptosis-vulnerability axis these cells move along).
//!
//! ## Design: pure functions, consumer mutates
//!
//! Like [`crate::oxygen`], [`crate::ph`], and [`crate::stromal`], this module
//! is pure: it returns the per-step increment / multiplier and the caller owns
//! the `persister_fraction` state (carried on
//! [`crate::biochem::CellState::persister_fraction`]) and decides how to apply
//! the effect. [`crate::biochem::sim_cell_step`] does **not** read the
//! persister fraction, so the core engine stays byte-identical; a consumer
//! (e.g. sim-tme-3d) composes the helpers around the step call.
//!
//! ## Per-step update
//!
//! The consumer applies [`step`] each step — a **competing-rate** integrator
//! where acquisition (drug-driven) and reversion both act simultaneously, so
//! sustained sub-saturating drug reaches a sub-cap equilibrium rather than
//! ratcheting monotonically to the cap (#262). [`acquire`] and [`revert`] remain
//! the individual rate terms (the acquisition increment and the reversion
//! decay) for use in isolation.
//!
//! ## Identity default
//!
//! [`PersisterConfig::default`] is the identity element (every rate is zero),
//! so all helpers are no-ops (`step`/`acquire`/`revert` return the input,
//! multipliers return `1.0`, the increment returns `0.0`). A run with the
//! identity config is therefore byte-identical to one with no persister model
//! at all. [`PersisterConfig::enabled`] supplies plausible (placeholder,
//! pending calibration) values that produce an observable effect.

use crate::params::PersisterConfig;

/// Increment `fraction` under drug exposure this step.
///
/// `drug_intensity` ∈ [0, 1] is the local effective drug level (e.g. a dose
/// schedule's `factor_at(step)` times spatial availability). Growth is
/// logistic toward `cfg.max_fraction`, so acquisition slows as the cell
/// saturates. Returns the input unchanged when the config is identity
/// (`acquisition_rate == 0`).
pub fn acquire(fraction: f64, drug_intensity: f64, cfg: &PersisterConfig) -> f64 {
    debug_assert!(cfg.acquisition_rate >= 0.0, "acquisition_rate must be >= 0");
    // No-op short-circuit: with the identity config `max_fraction == 0.0`, so
    // the final clamp would otherwise pull any positive `fraction` down to 0.
    // (Correctness here depends on the config invariant that acquisition_rate
    // and max_fraction are zero together — see PersisterConfig::default.)
    if cfg.acquisition_rate == 0.0 {
        return fraction;
    }
    let headroom = (cfg.max_fraction - fraction).max(0.0);
    let next = fraction + cfg.acquisition_rate * drug_intensity.clamp(0.0, 1.0) * headroom;
    next.clamp(0.0, cfg.max_fraction)
}

/// Exponentially relax `fraction` toward zero when the drug is absent
/// (epigenetic reversion). Returns the input unchanged when identity
/// (`reversion_rate == 0`).
pub fn revert(fraction: f64, cfg: &PersisterConfig) -> f64 {
    debug_assert!(
        (0.0..=1.0).contains(&cfg.reversion_rate),
        "reversion_rate must be in [0, 1]; a value > 1 inverts the decay sign"
    );
    (fraction * (1.0 - cfg.reversion_rate)).max(0.0)
}

/// One **competing-rate** per-step update of `fraction` (#262): acquisition and
/// reversion act **simultaneously**, not as an either-or keyed on whether the
/// drug is exactly zero. Per step:
///
/// `frac += acquisition_rate · drug · (max_fraction − frac)  −  reversion_rate · frac`
///
/// so under sustained sub-saturating drug the fraction relaxes to a **sub-cap
/// equilibrium** (where acquisition balances reversion) rather than ratcheting
/// monotonically to the cap — the biologically faithful behavior (Hangauer 2017;
/// Tsoi 2018). With `drug = 0` it reduces to pure reversion; with the identity
/// config (`acquisition_rate == 0`) it is a no-op, so a consumer that does not
/// opt in stays byte-identical. `drug_intensity` is clamped to `[0, 1]`; the
/// result to `[0, max_fraction]`.
///
/// This is the per-step integrator the consumer applies; [`acquire`] / [`revert`]
/// remain the individual rate terms (useful in isolation / for inspection).
pub fn step(fraction: f64, drug_intensity: f64, cfg: &PersisterConfig) -> f64 {
    debug_assert!(cfg.acquisition_rate >= 0.0, "acquisition_rate must be >= 0");
    debug_assert!(
        (0.0..=1.0).contains(&cfg.reversion_rate),
        "reversion_rate must be in [0, 1]"
    );
    // Identity config (acquisition_rate == 0, which by the config invariant
    // pairs with max_fraction == 0) ⇒ no-op ⇒ byte-identical.
    if cfg.acquisition_rate == 0.0 {
        return fraction;
    }
    let drug = drug_intensity.clamp(0.0, 1.0);
    let acquisition = cfg.acquisition_rate * drug * (cfg.max_fraction - fraction).max(0.0);
    let reversion = cfg.reversion_rate * fraction;
    (fraction + acquisition - reversion).clamp(0.0, cfg.max_fraction)
}

/// Multiplier applied to the per-step GPX4-inactivation drug effect: a
/// persister resists covalent GPX4 knockdown. Returns `1.0` at `fraction == 0`
/// or identity config (no attenuation) and `1 - gpx4_resistance` at full
/// persistence. Clamped to [0, 1].
pub fn gpx4_inactivation_multiplier(fraction: f64, cfg: &PersisterConfig) -> f64 {
    (1.0 - cfg.gpx4_resistance * fraction).clamp(0.0, 1.0)
}

/// Per-step additive MUFA-protection boost (same shape as
/// [`crate::stromal`]'s CAF MUFA supply), scaled by the persister fraction.
/// Returns `0.0` at `fraction == 0` or identity config. The caller adds this
/// to `CellState::mufa_protection` and clamps with `cfg.mufa_boost_cap`.
pub fn mufa_boost_increment(fraction: f64, cfg: &PersisterConfig) -> f64 {
    (cfg.mufa_boost_per_step * fraction).max(0.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn enabled() -> PersisterConfig {
        PersisterConfig::enabled()
    }

    #[test]
    fn identity_config_is_a_noop() {
        let id = PersisterConfig::default();
        assert!(id.is_identity());
        assert_eq!(acquire(0.3, 1.0, &id), 0.3);
        assert_eq!(revert(0.3, &id), 0.3);
        assert_eq!(gpx4_inactivation_multiplier(0.3, &id), 1.0);
        assert_eq!(mufa_boost_increment(0.3, &id), 0.0);
    }

    #[test]
    fn enabled_config_is_not_identity() {
        assert!(!enabled().is_identity());
    }

    #[test]
    fn acquisition_increases_and_saturates_at_max() {
        let cfg = enabled();
        let mut f = 0.0;
        for _ in 0..10_000 {
            f = acquire(f, 1.0, &cfg);
        }
        assert!(f > 0.0);
        assert!(f <= cfg.max_fraction + 1e-12);
        // One more step cannot exceed the cap.
        assert!(acquire(cfg.max_fraction, 1.0, &cfg) <= cfg.max_fraction + 1e-12);
    }

    #[test]
    fn acquisition_is_monotonic_nondecreasing() {
        let cfg = enabled();
        let f0 = 0.2;
        let f1 = acquire(f0, 1.0, &cfg);
        assert!(f1 >= f0);
    }

    #[test]
    fn zero_drug_intensity_does_not_acquire() {
        let cfg = enabled();
        assert_eq!(acquire(0.25, 0.0, &cfg), 0.25);
    }

    #[test]
    fn reversion_decays_toward_zero() {
        let cfg = enabled();
        let mut f = 0.8;
        for _ in 0..10_000 {
            f = revert(f, &cfg);
        }
        assert!(f < 1e-6);
    }

    #[test]
    fn gpx4_multiplier_bounds() {
        let cfg = enabled();
        assert_eq!(gpx4_inactivation_multiplier(0.0, &cfg), 1.0);
        let full = gpx4_inactivation_multiplier(1.0, &cfg);
        assert!((0.0..=1.0).contains(&full));
        assert!(full < 1.0); // some resistance when enabled
    }

    #[test]
    fn mufa_increment_scales_with_fraction() {
        let cfg = enabled();
        assert_eq!(mufa_boost_increment(0.0, &cfg), 0.0);
        assert!(mufa_boost_increment(0.5, &cfg) > 0.0);
        assert!(mufa_boost_increment(1.0, &cfg) > mufa_boost_increment(0.5, &cfg));
    }

    // ===== Competing-rate step (#262) =====

    #[test]
    fn step_identity_config_is_a_noop() {
        let id = PersisterConfig::default();
        assert_eq!(step(0.3, 1.0, &id), 0.3);
        assert_eq!(step(0.3, 0.0, &id), 0.3);
    }

    #[test]
    fn step_zero_drug_is_pure_reversion() {
        let cfg = enabled();
        let mut f = 0.8;
        for _ in 0..10_000 {
            f = step(f, 0.0, &cfg);
        }
        assert!(f < 1e-6, "no drug ⇒ reverts to zero; got {f}");
    }

    /// AC #4: under sustained **sub-saturating** drug, acquisition and reversion
    /// balance at a sub-cap equilibrium, instead of ratcheting monotonically to
    /// the cap the way the old either-or `acquire` did.
    #[test]
    fn step_reaches_sub_cap_equilibrium_under_sustained_subsaturating_drug() {
        let cfg = enabled();
        let drug = 0.3; // sub-saturating
        let mut f = 0.0;
        for _ in 0..10_000 {
            f = step(f, drug, &cfg);
        }
        // Converged strictly between 0 and the cap.
        assert!(f > 0.0, "acquisition fires under drug");
        assert!(
            f < cfg.max_fraction - 1e-6,
            "competing reversion holds it below the cap: f={f}, max={}",
            cfg.max_fraction
        );
        // Fixed point: another step barely moves it.
        assert!(
            (step(f, drug, &cfg) - f).abs() < 1e-9,
            "equilibrium is a fixed point"
        );
        // Matches the analytic equilibrium acq·drug·max / (rev + acq·drug).
        let a = cfg.acquisition_rate * drug;
        let f_star = a * cfg.max_fraction / (cfg.reversion_rate + a);
        assert!((f - f_star).abs() < 1e-3, "f={f}, analytic f*={f_star}");
    }

    /// The competing-rate signature vs the old `acquire`: reversion holds the
    /// fraction **below the cap even at saturating drug** (the old monotonic
    /// `acquire` reached the cap).
    #[test]
    fn step_equilibrium_stays_below_cap_even_at_full_drug() {
        let cfg = enabled();
        let mut f = 0.0;
        for _ in 0..10_000 {
            f = step(f, 1.0, &cfg);
        }
        let f_star =
            cfg.acquisition_rate * cfg.max_fraction / (cfg.reversion_rate + cfg.acquisition_rate);
        assert!(
            f < cfg.max_fraction,
            "reversion keeps it below the cap even at full drug: f={f}, max={}",
            cfg.max_fraction
        );
        assert!((f - f_star).abs() < 1e-3, "f={f}, analytic f*={f_star}");
    }
}
