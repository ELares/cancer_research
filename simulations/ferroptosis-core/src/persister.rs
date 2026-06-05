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
/// config (all rates zero) it is a no-op, so a consumer that does not opt in
/// stays byte-identical. `drug_intensity` is clamped to `[0, 1]`; the result to
/// `[0, max_fraction]`.
///
/// This is the per-step integrator the consumer applies; [`acquire`] / [`revert`]
/// remain the individual rate terms (useful in isolation / for inspection).
pub fn step(fraction: f64, drug_intensity: f64, cfg: &PersisterConfig) -> f64 {
    debug_assert!(cfg.acquisition_rate >= 0.0, "acquisition_rate must be >= 0");
    debug_assert!(
        (0.0..=1.0).contains(&cfg.reversion_rate),
        "reversion_rate must be in [0, 1]"
    );
    // Short-circuit only the *fully identity* config (all rates zero) ⇒ no-op ⇒
    // byte-identical. Guarding on `is_identity()` rather than just
    // `acquisition_rate == 0` keeps reversion live for a hypothetical
    // acquisition-off-but-reversion-on config (drug permanently withdrawn),
    // which should still decay rather than freeze (#262 review).
    if cfg.is_identity() {
        return fraction;
    }
    let drug = drug_intensity.clamp(0.0, 1.0);
    let acquisition = cfg.acquisition_rate * drug * (cfg.max_fraction - fraction).max(0.0);
    let reversion = cfg.reversion_rate * fraction;
    (fraction + acquisition - reversion).clamp(0.0, cfg.max_fraction)
}

/// Per-cell persister state with the irreversible (epigenetically locked)
/// sub-population split out (#342).
///
/// The base model treats all persistence as freely reversible. In reality,
/// beyond a threshold of SUSTAINED drug exposure persistence becomes
/// epigenetically locked and effectively irreversible (FSP1/HDAC-mediated
/// suppression of alternative defenses, Sci Adv 2026 PMID 41481741). This
/// splits the pool: `reversible` relaxes toward 0 after drug clearance (as
/// before), while `locked` does not revert once acquired.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PersisterState {
    /// Reversible persister fraction (relaxes toward 0 after drug clearance).
    pub reversible: f64,
    /// Locked persister fraction: acquired under sustained exposure, does NOT
    /// revert.
    pub locked: f64,
    /// Sustained-exposure tracker (an EMA of drug, see
    /// [`PersisterConfig::exposure_decay`]) that gates locking.
    pub cumulative_exposure: f64,
}

impl PersisterState {
    /// The empty state (no persistence, no exposure history).
    pub const ZERO: Self = PersisterState {
        reversible: 0.0,
        locked: 0.0,
        cumulative_exposure: 0.0,
    };

    /// Total persister fraction the consumer applies to the biochem effects
    /// (reversible + locked), clamped to the config ceiling `max_fraction`.
    /// With locking off, `locked == 0`, so this is just the reversible fraction.
    ///
    /// NOTE the `max_fraction` cap is enforced here, at the point of use, NOT
    /// jointly on the two pools: `reversible` and `locked` are each clamped to
    /// `max_fraction` individually, so under sustained acquisition + active
    /// locking the internal sum can transiently exceed `max_fraction` (up to
    /// 2×). This `total()` clamp is what bounds the value the biochem effects
    /// actually see, so the consumer should always read the pool through
    /// `total()` rather than summing the fields directly.
    pub fn total(&self, cfg: &PersisterConfig) -> f64 {
        (self.reversible + self.locked).clamp(0.0, cfg.max_fraction.max(0.0))
    }
}

/// One per-step update of the [`PersisterState`] with reversible-to-irreversible
/// locking (#342).
///
/// The `reversible` pool evolves by the existing competing-rate [`step`]. The
/// `cumulative_exposure` tracker is updated as an exponential moving average,
/// `cumulative·(1 - exposure_decay) + drug`, whose steady state is
/// `avg_drug / exposure_decay`; once it crosses `cfg.lock_threshold`, a fraction
/// `cfg.lock_rate` of the reversible pool is moved into the non-reverting
/// `locked` pool each step. Because the tracker decays during drug-off windows,
/// only SUSTAINED (continuous) exposure crosses the threshold; intermittent
/// dosing settles below it and never locks.
///
/// With `cfg.lock_rate == 0.0` (the default) this reduces EXACTLY to [`step`] on
/// the reversible pool with `locked` and `cumulative_exposure` left untouched,
/// so a consumer that has not opted into locking stays byte-identical.
pub fn step_with_locking(
    state: PersisterState,
    drug_intensity: f64,
    cfg: &PersisterConfig,
) -> PersisterState {
    // Reversible pool: the existing competing-rate update.
    let reversible = step(state.reversible, drug_intensity, cfg);
    // Locking off ⇒ behave exactly like `step` (locked + tracker untouched).
    if cfg.lock_rate == 0.0 {
        return PersisterState {
            reversible,
            locked: state.locked,
            cumulative_exposure: state.cumulative_exposure,
        };
    }
    let drug = drug_intensity.clamp(0.0, 1.0);
    let decay = cfg.exposure_decay.clamp(0.0, 1.0);
    // Sustained-exposure EMA: rises under drug, decays in drug-off windows.
    let cumulative_exposure = state.cumulative_exposure * (1.0 - decay) + drug;
    // Above threshold, ratchet reversible -> locked (irreversible).
    let (reversible, locked) = if cumulative_exposure >= cfg.lock_threshold {
        let lock_amount = (cfg.lock_rate.max(0.0) * reversible).min(reversible);
        (
            (reversible - lock_amount).max(0.0),
            (state.locked + lock_amount).clamp(0.0, cfg.max_fraction.max(0.0)),
        )
    } else {
        (reversible, state.locked)
    };
    PersisterState {
        reversible,
        locked,
        cumulative_exposure,
    }
}

/// Non-drug stress-niche persister entry (#377): raise the persister fraction in
/// a hypoxic / nutrient-poor drug-sanctuary niche, DECOUPLED from drug exposure.
///
/// The classic drug-tolerant-persister biology has a second, non-drug entry
/// route: hypoxic / nutrient-poor microenvironments drive cells into a
/// slow-cycling, drug-tolerant persister state independent of drug. `stress` ∈
/// [0, 1] is the local stress signal (e.g. `1 - o2_supply`, the hypoxia deficit).
///
/// The increment goes into the REVERSIBLE pool only — a stress-niche persister
/// reverts when the niche resolves (via the next [`step_with_locking`]'s
/// reversion) — so this does NOT feed the locking EMA (#342) or the drug-driven
/// resistance: **stress drives ENTRY, drug drives durability**. Logistic, so it
/// saturates at `cfg.max_fraction`:
///
/// `reversible += stress_entry_rate · stress · (max_fraction − reversible)`
///
/// `cfg.stress_entry_rate == 0.0` (the default) ⇒ no-op ⇒ byte-identical. Apply
/// it AFTER [`step_with_locking`] each step (the drug update first, then the
/// stress entry).
pub fn stress_entry(
    mut state: PersisterState,
    stress: f64,
    cfg: &PersisterConfig,
) -> PersisterState {
    if cfg.stress_entry_rate == 0.0 {
        return state;
    }
    let s = stress.clamp(0.0, 1.0);
    let headroom = (cfg.max_fraction - state.reversible).max(0.0);
    state.reversible =
        (state.reversible + cfg.stress_entry_rate * s * headroom).clamp(0.0, cfg.max_fraction);
    state
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
    fn stress_entry_is_noop_by_default_and_raises_reversible_under_stress() {
        // Off-by-default (`enabled()` keeps stress_entry_rate = 0): a no-op, so a
        // consumer that has not opted in is byte-identical.
        let off = enabled();
        assert_eq!(off.stress_entry_rate, 0.0);
        let s0 = PersisterState::ZERO;
        assert_eq!(stress_entry(s0, 0.9, &off), s0);

        // A stress-entry-only config (all drug rates zero, stress_entry on) is NOT
        // identity, and a HIGH-stress niche drives more reversible persister entry
        // than a LOW-stress one at the same (zero) drug.
        let cfg = PersisterConfig {
            stress_entry_rate: 0.05,
            max_fraction: 0.8,
            ..PersisterConfig::default()
        };
        assert!(
            !cfg.is_identity(),
            "a stress-entry-only config has an effect"
        );
        let run = |stress: f64| {
            let mut st = PersisterState::ZERO;
            for _ in 0..50 {
                // Drug-driven step first (zero drug here ⇒ no drug acquisition),
                // then the non-drug stress entry — the documented ordering.
                st = step_with_locking(st, 0.0, &cfg);
                st = stress_entry(st, stress, &cfg);
            }
            st
        };
        let hypoxic = run(0.9); // deep niche: 1 - o2_supply ≈ 0.9
        let normoxic = run(0.1); // rim: 1 - o2_supply ≈ 0.1
        assert!(
            hypoxic.reversible > normoxic.reversible,
            "stress entry must raise persister fraction more in the hypoxic niche: \
             hypoxic={}, normoxic={}",
            hypoxic.reversible,
            normoxic.reversible
        );
        // Stress entry feeds the REVERSIBLE pool only — it never touches the
        // locking EMA or the locked pool (stress drives entry, drug drives durability).
        assert_eq!(hypoxic.locked, 0.0);
        assert_eq!(hypoxic.cumulative_exposure, 0.0);
        // Saturates at the cap, never exceeds it.
        let saturated =
            (0..100_000).fold(PersisterState::ZERO, |st, _| stress_entry(st, 1.0, &cfg));
        assert!(saturated.reversible <= cfg.max_fraction + 1e-12);
        assert!(saturated.reversible > 0.79);
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

    /// #262 review: acquisition off but reversion on (NOT the identity config —
    /// e.g. drug permanently withdrawn) must still revert, not freeze. The
    /// short-circuit guards on `is_identity()`, so this config takes the full
    /// update and decays.
    #[test]
    fn step_reverts_when_acquisition_off_but_reversion_on() {
        let cfg = PersisterConfig {
            acquisition_rate: 0.0,
            reversion_rate: 0.1,
            max_fraction: 0.8,
            ..PersisterConfig::enabled()
        };
        assert!(!cfg.is_identity(), "rev>0 ⇒ not the identity config");
        // pure reversion: 0.5·(1−0.1) = 0.45, regardless of drug level.
        let next = step(0.5, 1.0, &cfg);
        assert!(
            (next - 0.45).abs() < 1e-12,
            "should revert toward 0: got {next}"
        );
    }

    // ===== Reversible-to-irreversible locking (#342) =====

    /// With `lock_rate == 0` (the default), `step_with_locking` reduces EXACTLY
    /// to `step` on the reversible pool, with `locked` and the exposure tracker
    /// left untouched. This is the byte-identity invariant a consumer relies on.
    #[test]
    fn step_with_locking_matches_step_when_locking_off() {
        let cfg = PersisterConfig::enabled(); // lock_rate defaults to 0.0
        let mut s = PersisterState::ZERO;
        let mut f = 0.0_f64;
        for t in 0..200 {
            let drug = if t < 100 { 0.7 } else { 0.0 };
            s = step_with_locking(s, drug, &cfg);
            f = step(f, drug, &cfg);
            assert_eq!(s.reversible, f, "reversible must track plain step exactly");
            assert_eq!(s.locked, 0.0, "no locking when lock_rate=0");
            assert_eq!(
                s.cumulative_exposure, 0.0,
                "tracker untouched when locking off"
            );
        }
        assert_eq!(s.total(&cfg), f);
    }

    /// #342: sustained CONTINUOUS dosing crosses the lock threshold and leaves
    /// an irreversible (locked) persister fraction that survives drug
    /// withdrawal, whereas INTERMITTENT dosing keeps the sustained-exposure EMA
    /// below the threshold, never locks, and reverts toward zero after washout.
    #[test]
    fn continuous_dosing_locks_persisters_but_intermittent_does_not() {
        // EMA steady state = avg_drug / exposure_decay: continuous ⇒ 1/0.1 = 10
        // (> 5), intermittent at 25% duty ⇒ 0.25/0.1 = 2.5 (< 5). So 5.0 sits
        // between them and only continuous dosing locks.
        let cfg = PersisterConfig {
            lock_rate: 0.1,
            lock_threshold: 5.0,
            exposure_decay: 0.1,
            ..PersisterConfig::enabled()
        };
        let run = |schedule: &dyn Fn(usize) -> f64| -> PersisterState {
            let mut s = PersisterState::ZERO;
            for t in 0..160 {
                s = step_with_locking(s, schedule(t), &cfg);
            }
            // Long washout (drug fully withdrawn).
            for _ in 0..160 {
                s = step_with_locking(s, 0.0, &cfg);
            }
            s
        };
        let continuous = run(&|_| 1.0);
        let intermittent = run(&|t| if t % 4 == 0 { 1.0 } else { 0.0 }); // 25% duty

        // Continuous: a locked fraction persists through washout.
        assert!(
            continuous.locked > 0.01,
            "continuous dosing should lock persisters: {continuous:?}"
        );
        assert!(
            continuous.total(&cfg) > 0.01,
            "the locked fraction survives drug withdrawal: {continuous:?}"
        );
        // Intermittent: never crosses the threshold ⇒ no locking; reverts away.
        assert!(
            intermittent.locked < 1e-9,
            "intermittent dosing must not lock: {intermittent:?}"
        );
        assert!(
            intermittent.total(&cfg) < continuous.total(&cfg),
            "intermittent reverts further than continuous: int={:?} cont={:?}",
            intermittent,
            continuous
        );
        // Determinism.
        assert_eq!(continuous, run(&|_| 1.0));
    }
}
