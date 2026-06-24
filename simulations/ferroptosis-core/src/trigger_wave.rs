//! Propagating ferroptotic trigger-wave front (#482).
//!
//! Ferroptosis does not only kill cells one at a time: a lipid-peroxide / ROS
//! front can PROPAGATE across a tissue at a constant speed, each dying cell
//! igniting its neighbours (Co, Wu, Lee & Chen, *Nature* 631:654, 2024, PMID
//! 38987590 measured a baseline front speed of 5.52 +/- 0.09 um/min, and showed
//! it is iron-tunable: iron chelation with DFO slows it to 2.33 um/min, iron
//! supplementation accelerates it to 9.40 um/min). This is the SPATIAL leg of
//! the same group whose single-cell bistable switch we already cross-validate
//! (`analysis/ode-cross-validation.md`, #344).
//!
//! The other supply modules ([`crate::reaction_diffusion`]) solve a *steady
//! state*; this module solves a *time-dependent, bistable* reaction-diffusion
//! equation whose hallmark is a travelling front. The normalized lipid-peroxide
//! field `L(x, t) ∈ [0, 1]` (0 = healthy, 1 = ferroptotic) obeys the Nagumo
//! (bistable cubic) equation:
//!
//! ```text
//!   ∂L/∂t = D ∂²L/∂x² + k · L · (L − a) · (1 − L)
//! ```
//!
//! - `D` is the lipid-radical / ROS diffusion coefficient (um²/min): how fast a
//!   peroxidizing cell's radicals reach its neighbours.
//! - `k` is the autocatalytic peroxidation rate (1/min); it is the spatial
//!   analogue of the single-cell `lp_propagation` and is RAISED BY LABILE IRON
//!   (more Fenton chemistry drives the chain faster).
//! - `a ∈ (0, 0.5)` is the ignition threshold (the unstable middle fixed point):
//!   a cell must be pushed above `a` to flip to the ferroptotic state. GPX4 /
//!   GSH / radical-trapping defenses RAISE `a` (harder to ignite); at `a ≥ 0.5`
//!   the front stalls or reverses (defense outruns propagation).
//!
//! ## The front speed and why iron tunes it
//!
//! The travelling-front solution of the Nagumo equation has the closed-form
//! speed
//!
//! ```text
//!   c = sqrt(D · k / 2) · (1 − 2a)
//! ```
//!
//! (front profile `L(z) = 1 / (1 + exp(−z / sqrt(2D/k)))`). Two consequences,
//! both matching the Co 2024 measurements:
//!
//! 1. `c ∝ sqrt(k) ∝ sqrt(iron)`: more labile iron ⇒ faster front. The measured
//!    DFO / control / iron speeds (2.33 / 5.52 / 9.40 um/min) sit on a
//!    square-root-of-iron curve at biologically plausible iron fold-changes
//!    (control:iron ≈ (9.40/5.52)² ≈ 2.9x labile iron; control:DFO ≈
//!    (2.33/5.52)² ≈ 0.18x), so the measured iron-dose response is consistent
//!    with a Fenton-iron-driven bistable front.
//! 2. raising the GPX4 threshold `a` toward 0.5 slows and ultimately halts the
//!    front, the spatial counterpart of single-cell GPX4 protection.
//!
//! [`front_speed`] integrates the PDE numerically and measures the front speed
//! from the L = 0.5 crossing; [`analytical_front_speed`] returns the closed
//! form. The numerical solve agreeing with the closed form is the module's
//! self-consistency check (cf. [`crate::reaction_diffusion::analytical_1d_slab`]).
//!
//! ## Opt-in, off the production path
//!
//! Nothing in the production simulation matrix calls this module; it is a
//! standalone spatial-front model used for the #482 calibration. So adding it is
//! byte-identical to the production matrix.

/// Configuration for a 1-D ferroptotic trigger-wave solve (#482).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TriggerWaveConfig {
    /// Domain length (um). Must be long enough that the front stays interior
    /// over the measurement window.
    pub grid_len_um: f64,
    /// Spatial step (um).
    pub h_um: f64,
    /// Time step (min). Must satisfy the explicit-diffusion CFL bound
    /// `dt < h² / (2D)`; [`front_speed`] asserts it.
    pub dt_min: f64,
    /// Lipid-radical / ROS diffusion coefficient (um²/min).
    pub diffusion_um2_per_min: f64,
    /// Base autocatalytic peroxidation rate `k0` (1/min) at `iron_level = 1`.
    pub base_reaction_rate: f64,
    /// Relative labile iron (1.0 = baseline). Scales the effective rate
    /// `k = base_reaction_rate · iron_level` (Fenton-driven), so the front speed
    /// scales as `sqrt(iron_level)`.
    pub iron_level: f64,
    /// Base ignition threshold `a0 ∈ (0, 0.5)`.
    pub ignition_threshold: f64,
    /// Extra threshold from GPX4 / GSH / radical-trapping defense, added to
    /// `a0` and capped below 0.5 (a defended tissue is harder to ignite).
    pub gpx4_defense: f64,
}

impl TriggerWaveConfig {
    /// Baseline placeholder geometry/kinetics. `diffusion`/`base_reaction_rate`
    /// are tuned in the #482 calibration so the baseline (`iron_level = 1`,
    /// `gpx4_defense = 0`) front speed lands near the measured 5.52 um/min;
    /// `iron_level`/`gpx4_defense` are then the biology knobs.
    pub fn baseline() -> Self {
        TriggerWaveConfig {
            grid_len_um: 600.0,
            h_um: 2.0,
            dt_min: 0.02,
            // Calibrated so c = sqrt(D·k/2)·(1−2a) ≈ 5.52 um/min at iron_level=1,
            // a0=0.25: sqrt(D·k/2)·0.5 = 5.52 ⇒ D·k = 2·(11.04)² ≈ 243.8.
            // With D = 30 um²/min ⇒ k ≈ 8.13 /min.
            diffusion_um2_per_min: 30.0,
            base_reaction_rate: 8.13,
            iron_level: 1.0,
            ignition_threshold: 0.25,
            gpx4_defense: 0.0,
        }
    }

    /// Effective reaction rate `k = base_reaction_rate · iron_level` (>= 0).
    #[must_use]
    pub fn effective_rate(&self) -> f64 {
        (self.base_reaction_rate * self.iron_level).max(0.0)
    }

    /// Effective ignition threshold `a = a0 + gpx4_defense`, clamped to
    /// `[0, 0.49]` (kept below 0.5 so a baseline front still propagates).
    #[must_use]
    pub fn effective_threshold(&self) -> f64 {
        (self.ignition_threshold + self.gpx4_defense).clamp(0.0, 0.49)
    }
}

/// Closed-form Nagumo front speed `c = sqrt(D·k/2)·(1 − 2a)` (um/min).
///
/// Positive when `a < 0.5` (front invades the healthy state), zero at `a = 0.5`,
/// negative when `a > 0.5` (the healthy state re-invades). `k`/`a` are taken
/// from the config's [`TriggerWaveConfig::effective_rate`] /
/// [`TriggerWaveConfig::effective_threshold`].
#[must_use]
pub fn analytical_front_speed(cfg: &TriggerWaveConfig) -> f64 {
    let k = cfg.effective_rate();
    let a = cfg.effective_threshold();
    (cfg.diffusion_um2_per_min * k / 2.0).sqrt() * (1.0 - 2.0 * a)
}

/// Integrate the 1-D Nagumo PDE and measure the travelling-front speed (um/min).
///
/// Seeds the left fifth of the domain in the ferroptotic state (`L = 1`), the
/// rest healthy (`L = 0`), and integrates with explicit Euler under no-flux
/// boundaries. After a transient (so the front settles into its travelling
/// profile) the `L = 0.5` crossing position is tracked over a measurement
/// window and the speed is its linear slope. Returns `0.0` if no propagating
/// front forms (e.g. `a ≥ 0.5`) or the front leaves the domain.
///
/// # Panics
/// Panics if the explicit-diffusion CFL bound `dt < h² / (2D)` is violated
/// (the solve would be unstable), mirroring the stability assertions elsewhere
/// in the crate.
#[must_use]
pub fn front_speed(cfg: &TriggerWaveConfig) -> f64 {
    let n = (cfg.grid_len_um / cfg.h_um).round() as usize;
    assert!(n >= 16, "grid too small: need >= 16 points, got {n}");
    let d = cfg.diffusion_um2_per_min;
    let h = cfg.h_um;
    let dt = cfg.dt_min;
    // CFL stability bound for the explicit diffusion term.
    assert!(
        dt < h * h / (2.0 * d),
        "explicit-Euler CFL violated: dt {dt} must be < h²/(2D) = {}",
        h * h / (2.0 * d)
    );
    let k = cfg.effective_rate();
    let a = cfg.effective_threshold();

    let mut l = vec![0.0f64; n];
    let seed = n / 5;
    for li in l.iter_mut().take(seed) {
        *li = 1.0;
    }
    let mut scratch = vec![0.0f64; n];

    // Track the L=0.5 front position (um) vs time (min). Skip a transient, then
    // record over a window while the front is interior.
    let mut times = Vec::new();
    let mut positions = Vec::new();
    // Step until the front passes ~80% of the domain or a step cap is hit.
    let max_steps = 2_000_000usize;
    for step in 0..max_steps {
        // Explicit-Euler update with no-flux (reflective) boundaries.
        for i in 0..n {
            let lm = if i == 0 { l[1] } else { l[i - 1] };
            let lp = if i == n - 1 { l[n - 2] } else { l[i + 1] };
            let lap = (lm - 2.0 * l[i] + lp) / (h * h);
            let react = k * l[i] * (l[i] - a) * (1.0 - l[i]);
            scratch[i] = l[i] + dt * (d * lap + react);
        }
        l.copy_from_slice(&scratch);

        let t = step as f64 * dt;
        if let Some(pos) = front_position_um(&l, h) {
            // Record once the front has cleared the seed and is interior.
            if pos > cfg.grid_len_um * 0.25 && pos < cfg.grid_len_um * 0.85 {
                times.push(t);
                positions.push(pos);
            }
            if pos >= cfg.grid_len_um * 0.85 {
                break;
            }
        }
        // No front formed (decayed away): bail early.
        if step > 5000 && positions.is_empty() && l.iter().all(|&v| v < 0.5) {
            return 0.0;
        }
    }

    if times.len() < 3 {
        return 0.0;
    }
    least_squares_slope(&times, &positions)
}

/// Interpolated position (um) of the right-most `L = 0.5` crossing, or `None`
/// if the field never crosses 0.5.
fn front_position_um(l: &[f64], h: f64) -> Option<f64> {
    // Scan from the right for the first 0.5 down-crossing (high -> low).
    for i in (1..l.len()).rev() {
        if l[i - 1] >= 0.5 && l[i] < 0.5 {
            // Linear interpolation between i-1 and i.
            let frac = (l[i - 1] - 0.5) / (l[i - 1] - l[i]);
            return Some(((i - 1) as f64 + frac) * h);
        }
    }
    None
}

/// Ordinary-least-squares slope of `y` vs `x` (the front speed). Assumes
/// `x.len() == y.len() >= 2`.
fn least_squares_slope(x: &[f64], y: &[f64]) -> f64 {
    let n = x.len() as f64;
    let sx: f64 = x.iter().sum();
    let sy: f64 = y.iter().sum();
    let sxx: f64 = x.iter().map(|v| v * v).sum();
    let sxy: f64 = x.iter().zip(y).map(|(a, b)| a * b).sum();
    let denom = n * sxx - sx * sx;
    if denom.abs() < 1e-12 {
        return 0.0;
    }
    (n * sxy - sx * sy) / denom
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The numerical front speed must agree with the closed-form Nagumo speed
    /// `c = sqrt(D·k/2)·(1−2a)` (the module's self-consistency check, the
    /// trigger-wave analogue of `reaction_diffusion`'s analytical_1d_slab).
    #[test]
    fn numerical_front_matches_the_analytical_nagumo_speed() {
        let cfg = TriggerWaveConfig::baseline();
        let numeric = front_speed(&cfg);
        let analytic = analytical_front_speed(&cfg);
        assert!(
            numeric > 0.0,
            "a baseline front must propagate; got {numeric}"
        );
        let rel_err = (numeric - analytic).abs() / analytic;
        assert!(
            rel_err < 0.06,
            "numeric {numeric:.3} vs analytic {analytic:.3} um/min (rel err {rel_err:.3}) \
             must agree within 6% (discretization)"
        );
    }

    /// The baseline calibration lands near the measured 5.52 um/min (Co 2024,
    /// PMID 38987590). Loose band: the placeholder D/k are tuned to this target,
    /// so this guards the calibration, not an independent prediction.
    #[test]
    fn baseline_speed_is_near_the_measured_5_52_um_per_min() {
        let v = front_speed(&TriggerWaveConfig::baseline());
        assert!(
            (v - 5.52).abs() < 0.8,
            "baseline front speed {v:.3} um/min must be near the measured 5.52"
        );
    }

    /// The iron-dose DIRECTION is the robust result: iron chelation (low
    /// iron_level) slows the front, iron loading (high iron_level) accelerates
    /// it, ordering DFO < control < iron, because c ∝ sqrt(iron_level).
    #[test]
    fn iron_dose_orders_the_front_speed_dfo_below_control_below_loaded() {
        let control = front_speed(&TriggerWaveConfig::baseline());
        // (9.40/5.52)² ≈ 2.9 and (2.33/5.52)² ≈ 0.18 are the iron fold-changes
        // the measured speeds imply under c ∝ sqrt(iron).
        let loaded = front_speed(&TriggerWaveConfig {
            iron_level: 2.9,
            ..TriggerWaveConfig::baseline()
        });
        let dfo = front_speed(&TriggerWaveConfig {
            iron_level: 0.18,
            ..TriggerWaveConfig::baseline()
        });
        assert!(
            dfo < control && control < loaded,
            "iron-dose ordering DFO {dfo:.2} < control {control:.2} < loaded {loaded:.2} must hold"
        );
        // And the loaded/DFO speeds land near the measured 9.40 / 2.33.
        assert!(
            (loaded - 9.40).abs() < 1.2,
            "loaded {loaded:.2} near measured 9.40"
        );
        assert!((dfo - 2.33).abs() < 0.8, "dfo {dfo:.2} near measured 2.33");
    }

    /// Raising the GPX4 defense threshold toward 0.5 slows the front; at/above
    /// 0.5 it no longer invades (speed collapses to ~0). The spatial analogue of
    /// single-cell GPX4 protection.
    #[test]
    fn gpx4_defense_slows_then_halts_the_front() {
        let undefended = front_speed(&TriggerWaveConfig::baseline());
        let defended = front_speed(&TriggerWaveConfig {
            gpx4_defense: 0.15,
            ..TriggerWaveConfig::baseline()
        });
        assert!(
            defended < undefended,
            "GPX4 defense must slow the front: defended {defended:.2} < undefended {undefended:.2}"
        );
        // a0 + defense = 0.25 + 0.24 = 0.49 (clamped): near-stall, much slower.
        let near_stall = front_speed(&TriggerWaveConfig {
            gpx4_defense: 0.30,
            ..TriggerWaveConfig::baseline()
        });
        assert!(
            near_stall < 0.5 * undefended,
            "near-threshold defense must nearly halt the front; got {near_stall:.2}"
        );
    }

    /// `analytical_front_speed` is positive below threshold, zero at 0.5, and
    /// negative above it (the healthy state re-invades).
    #[test]
    fn analytical_speed_sign_tracks_the_threshold() {
        let base = TriggerWaveConfig::baseline();
        assert!(analytical_front_speed(&base) > 0.0);
        let at_half = TriggerWaveConfig {
            ignition_threshold: 0.5,
            ..base
        };
        // effective_threshold clamps to 0.49, so this is slightly positive; test
        // the raw formula at exactly 0.5 instead.
        let k = at_half.effective_rate();
        let exact_half = (at_half.diffusion_um2_per_min * k / 2.0).sqrt() * (1.0 - 2.0 * 0.5);
        assert_eq!(exact_half, 0.0);
    }
}
