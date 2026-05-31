//! Time-varying drug-administration schedules for spatial simulations (#239).
//!
//! The spatial binaries (sim-tme, sim-tme-3d) historically modeled drug as
//! present at constant strength for the entire run. Real treatment is
//! time-varying: a bolus rises and decays, multi-dose protocols pulse,
//! infusions hold a level for a window, and full PK models produce an
//! arbitrary concentration timecourse. `DoseSchedule` captures all four
//! shapes behind a single per-step factor.
//!
//! ## Design (mirrors `photosensitizer_pk::Photosensitizer`)
//!
//! - **Identity default.** [`DoseSchedule::Constant`] returns exactly `1.0`
//!   from [`factor_at`](DoseSchedule::factor_at) for every step. When all
//!   simulation conditions use `Constant`, the run is byte-identical to the
//!   pre-#239 steady-state behavior (`x * 1.0 == x` is exact in IEEE-754).
//! - **Per-step evaluation** (like sim-tme's `run_spatial_cycling` O₂-λ
//!   recompute): the consumer calls `factor_at(step)` once per step and
//!   composes the result multiplicatively with the modality's existing
//!   drug mechanism.
//!
//! ## How consumers apply the factor
//!
//! - **SDT / PDT** (exogenous ROS): scale the per-cell exogenous-ROS bolus
//!   magnitude — `effective_peak = base_peak * factor_at(step)`.
//! - **RSL3** (covalent GPX4 inhibitor): drive per-step GPX4 inactivation —
//!   `gpx4 -= k_inact * factor_at(step) * gpx4`, mirroring the validated
//!   `tumor_pk::sim_cell_with_pk` mechanism. (The `Constant` steady-state
//!   keeps RSL3's original one-shot init knockdown instead; see
//!   `biochem::CellState::from_cell_with_ros_opts`.)
//!
//! The factor is dimensionless drug *availability* in `[0, ∞)`. For the
//! analytic shapes it is normalized so a single full-strength dose peaks at
//! `peak` (typically `1.0`).
//!
//! ## Scope / future evolution (TODO #240)
//!
//! `factor_at(step)` is a single **global, time-only** scalar: at each step
//! every cell sees the same availability. That is correct for #239's
//! uniform-distribution premise, but it does **not** model spatial drug
//! heterogeneity — radial penetration gradients, vessel-distance occlusion,
//! interstitial-pressure exclusion. The patient-scale / vasculature work
//! (#240, #191) will need per-cell availability that varies in space as well
//! as time. When that lands, the likely shape is a split: keep this type as
//! the **temporal** schedule (rename toward `DoseTiming`) and introduce a
//! separate **spatial** `DrugAvailability(cell) -> f64` that composes with
//! it multiplicatively. The cell-type-efficacy (#241) and immune
//! dose-response (#243) work will probably ride on that same split. Capturing
//! the boundary here so the eventual refactor is planned, not a surprise.

use serde::{Deserialize, Serialize};

/// A drug-administration schedule producing a per-step availability factor.
///
/// `Constant` is the default and the identity (always `1.0`). The other
/// variants are time-varying. See the module docs for how each modality
/// applies the factor.
#[derive(Clone, Debug, PartialEq, Default, Serialize, Deserialize)]
pub enum DoseSchedule {
    /// Constant full availability for the whole run: `factor_at ≡ 1.0`.
    /// The identity case — preserves byte-identical pre-#239 output.
    #[default]
    Constant,

    /// Single bolus administered at `dose_step`: availability is `0.0`
    /// before the dose, jumps to `peak` at `dose_step`, then decays
    /// exponentially with `half_life_steps`.
    Bolus {
        /// Step at which the dose is administered.
        dose_step: u32,
        /// Peak availability at the moment of administration.
        peak: f64,
        /// Exponential decay half-life, in steps. Must be > 0.
        half_life_steps: f64,
    },

    /// Repeated boluses, one per entry in `dose_steps`. The factor at any
    /// step is the **sum** of every dose's exponentially-decaying
    /// contribution that has been administered by then — so overlapping
    /// doses accumulate (clinically realistic for short inter-dose
    /// intervals). Models a multi-dose protocol.
    ///
    /// Because contributions sum, `factor_at` can **exceed `peak`** when
    /// doses overlap. How a consumer treats the overshoot is modality-
    /// dependent and asymmetric by design: the RSL3 path clamps the factor
    /// into `[0, 1]` (it's a drug *concentration*, and the validated
    /// `sim_cell_with_pk` model is calibrated for `conc ≤ 1`), whereas the
    /// SDT/PDT path does **not** clamp (exogenous ROS is additive across
    /// overlapping pulses, so >1 is physically meaningful). Pick `peak` and
    /// spacing with the consuming modality's clamp in mind.
    MultiDose {
        /// Steps at which doses are administered (need not be sorted).
        dose_steps: Vec<u32>,
        /// Peak availability contributed by each individual dose.
        peak: f64,
        /// Per-dose exponential decay half-life, in steps. Must be > 0.
        half_life_steps: f64,
    },

    /// Continuous infusion: availability is `level` for steps in
    /// `[start, end)`, `0.0` otherwise. Models a constant IV drip.
    Infusion {
        /// First step of the infusion window (inclusive).
        start: u32,
        /// One past the last step of the infusion window (exclusive).
        end: u32,
        /// Availability held during the window.
        level: f64,
    },

    /// Explicit per-step availability series — e.g. the normalized
    /// interstitial-concentration timecourse from
    /// `tumor_pk::solve_tumor_pk`. `factor_at(step)` reads `series[step]`,
    /// clamping to the last value past the end (and `0.0` if empty). The
    /// bridge that lets the validated two-compartment PK ODE drive the
    /// spatial grid without coupling the grid to the solver.
    FromPk {
        /// Per-step availability factors (index = step).
        series: Vec<f64>,
    },
}

impl DoseSchedule {
    /// Per-step drug-availability factor in `[0, ∞)`.
    ///
    /// `Constant ⇒ 1.0` exactly for all steps — the identity that keeps
    /// the default simulation path byte-identical.
    #[must_use]
    pub fn factor_at(&self, step: u32) -> f64 {
        match self {
            DoseSchedule::Constant => 1.0,
            DoseSchedule::Bolus {
                dose_step,
                peak,
                half_life_steps,
            } => bolus_contribution(step, *dose_step, *peak, *half_life_steps),
            DoseSchedule::MultiDose {
                dose_steps,
                peak,
                half_life_steps,
            } => dose_steps
                .iter()
                .map(|&d| bolus_contribution(step, d, *peak, *half_life_steps))
                .sum(),
            DoseSchedule::Infusion { start, end, level } => {
                if step >= *start && step < *end {
                    *level
                } else {
                    0.0
                }
            }
            DoseSchedule::FromPk { series } => {
                if series.is_empty() {
                    0.0
                } else {
                    series[(step as usize).min(series.len() - 1)]
                }
            }
        }
    }

    /// `true` only for [`DoseSchedule::Constant`]. Consumers use this to
    /// skip **all** per-step dose modulation on the default path, which is
    /// what guarantees byte-identical output when no schedule is set.
    #[must_use]
    pub fn is_constant(&self) -> bool {
        matches!(self, DoseSchedule::Constant)
    }

    /// The discrete administration steps for the analytic shapes (for
    /// metadata / visualization annotation). `Constant` and `FromPk` have no
    /// discrete dose events and return an empty Vec. `Bolus` returns its
    /// single `dose_step`; `MultiDose` returns its `dose_steps` (sorted);
    /// `Infusion` returns its `start` as a single "begins here" marker.
    #[must_use]
    pub fn dose_steps(&self) -> Vec<u32> {
        match self {
            DoseSchedule::Constant | DoseSchedule::FromPk { .. } => Vec::new(),
            DoseSchedule::Bolus { dose_step, .. } => vec![*dose_step],
            DoseSchedule::MultiDose { dose_steps, .. } => {
                let mut v = dose_steps.clone();
                v.sort_unstable();
                v
            }
            DoseSchedule::Infusion { start, .. } => vec![*start],
        }
    }
}

/// One bolus's exponentially-decaying contribution at `step`: `0.0` before
/// `dose_step`, `peak` at it, halving every `half_life_steps`.
#[inline]
fn bolus_contribution(step: u32, dose_step: u32, peak: f64, half_life_steps: f64) -> f64 {
    debug_assert!(
        half_life_steps > 0.0,
        "DoseSchedule half_life_steps must be > 0, got {half_life_steps}"
    );
    if step < dose_step {
        return 0.0;
    }
    let elapsed = (step - dose_step) as f64;
    // Guard against a non-positive half-life in release (debug_assert catches
    // it in dev): treat as an instantaneous spike that's gone by the next step.
    if half_life_steps <= 0.0 {
        return if step == dose_step { peak } else { 0.0 };
    }
    peak * 0.5_f64.powf(elapsed / half_life_steps)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn constant_is_identity_for_all_steps() {
        let s = DoseSchedule::Constant;
        for step in [0u32, 1, 30, 60, 179, 1000] {
            assert_eq!(s.factor_at(step), 1.0, "Constant must be exactly 1.0");
        }
        assert!(s.is_constant());
    }

    #[test]
    fn default_is_constant() {
        assert_eq!(DoseSchedule::default(), DoseSchedule::Constant);
        assert!(DoseSchedule::default().is_constant());
    }

    #[test]
    fn bolus_zero_before_dose_peak_at_dose_half_at_halflife() {
        let s = DoseSchedule::Bolus {
            dose_step: 10,
            peak: 1.0,
            half_life_steps: 12.0,
        };
        assert_eq!(s.factor_at(9), 0.0, "no drug before the dose");
        assert_eq!(s.factor_at(10), 1.0, "peak at dose step");
        assert!(
            (s.factor_at(22) - 0.5).abs() < 1e-12,
            "half at one half-life later, got {}",
            s.factor_at(22)
        );
        assert!(!s.is_constant());
    }

    #[test]
    fn multidose_sums_overlapping_contributions() {
        let s = DoseSchedule::MultiDose {
            dose_steps: vec![10, 20],
            peak: 1.0,
            half_life_steps: 12.0,
        };
        // At step 9: nothing yet.
        assert_eq!(s.factor_at(9), 0.0);
        // At step 10: first dose peaks (1.0), second not yet → 1.0.
        assert!((s.factor_at(10) - 1.0).abs() < 1e-12);
        // At step 20: first dose decayed 10 steps + second dose peaks.
        let first_at_20 = 0.5_f64.powf(10.0 / 12.0);
        let expected = first_at_20 + 1.0;
        assert!(
            (s.factor_at(20) - expected).abs() < 1e-12,
            "doses must sum: got {}, expected {}",
            s.factor_at(20),
            expected
        );
    }

    #[test]
    fn infusion_holds_level_inside_window_only() {
        let s = DoseSchedule::Infusion {
            start: 30,
            end: 90,
            level: 0.7,
        };
        assert_eq!(s.factor_at(29), 0.0);
        assert_eq!(s.factor_at(30), 0.7, "level at window start (inclusive)");
        assert_eq!(s.factor_at(89), 0.7);
        assert_eq!(s.factor_at(90), 0.0, "window end is exclusive");
    }

    #[test]
    fn frompk_reads_series_and_clamps_past_end() {
        let s = DoseSchedule::FromPk {
            series: vec![0.0, 0.5, 1.0, 0.8],
        };
        assert_eq!(s.factor_at(0), 0.0);
        assert_eq!(s.factor_at(2), 1.0);
        assert_eq!(s.factor_at(3), 0.8);
        assert_eq!(s.factor_at(99), 0.8, "clamps to last value past end");
    }

    #[test]
    fn frompk_empty_series_is_zero() {
        let s = DoseSchedule::FromPk { series: vec![] };
        assert_eq!(s.factor_at(0), 0.0);
        assert_eq!(s.factor_at(50), 0.0);
    }

    #[test]
    fn non_constant_variants_report_not_constant() {
        for s in [
            DoseSchedule::Bolus {
                dose_step: 0,
                peak: 1.0,
                half_life_steps: 1.0,
            },
            DoseSchedule::MultiDose {
                dose_steps: vec![0],
                peak: 1.0,
                half_life_steps: 1.0,
            },
            DoseSchedule::Infusion {
                start: 0,
                end: 1,
                level: 1.0,
            },
            DoseSchedule::FromPk { series: vec![1.0] },
        ] {
            assert!(!s.is_constant(), "{s:?} must not report is_constant");
        }
    }

    #[test]
    fn dose_steps_reports_administration_events() {
        assert_eq!(DoseSchedule::Constant.dose_steps(), Vec::<u32>::new());
        assert_eq!(
            DoseSchedule::Bolus {
                dose_step: 15,
                peak: 1.0,
                half_life_steps: 5.0
            }
            .dose_steps(),
            vec![15]
        );
        assert_eq!(
            DoseSchedule::MultiDose {
                dose_steps: vec![30, 10, 20],
                peak: 1.0,
                half_life_steps: 5.0
            }
            .dose_steps(),
            vec![10, 20, 30],
            "dose_steps are sorted"
        );
        assert_eq!(
            DoseSchedule::Infusion {
                start: 40,
                end: 80,
                level: 1.0
            }
            .dose_steps(),
            vec![40]
        );
        assert_eq!(
            DoseSchedule::FromPk {
                series: vec![1.0, 2.0]
            }
            .dose_steps(),
            Vec::<u32>::new()
        );
    }
}
