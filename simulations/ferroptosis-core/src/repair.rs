//! ESCRT-III membrane-repair brake on ferroptosis death execution (#465).
//!
//! Every other layer in this engine modulates the redox / lipid-substrate axis
//! (iron, GSH/GPX4, FSP1, MUFA, PUFA, ALOX, ACSL4): they set WHETHER a cell
//! peroxidizes. This layer is different. It acts on the death-EXECUTION step: a
//! cell whose lipid peroxide has already crossed the death threshold can still be
//! rescued, for a time, by ESCRT-III-dependent plasma-membrane repair.
//!
//! The biology: ferroptotic lipid-peroxidation pores let Ca2+ in, and the Ca2+
//! influx recruits ESCRT-III (CHMP5/CHMP6) to reseal the damaged membrane,
//! counterbalancing death kinetics. Knocking down CHMP5/CHMP6 SENSITIZES cells to
//! lipid-peroxidation death, so more ESCRT-III repair capacity means slower or
//! blocked execution and therefore more ferroptosis RESISTANCE (Dai et al.,
//! "ESCRT-III-dependent membrane repair blocks ferroptosis," Biochem Biophys Res
//! Commun 2020, PMID 31761326). This is the membrane-repair axis the manuscript
//! §8.4 flagged as missing.
//!
//! ## What this module models
//!
//! The brake is governed by two off-by-default [`crate::params::Params`] fields:
//!   - `escrt_repair_rate`: the per-step probability that, when a cell's LP has
//!     crossed `death_threshold`, ESCRT reseals the membrane and the cell survives
//!     THIS step (instead of dying). `0.0` (default) ⇒ the brake never fires.
//!   - `escrt_repair_budget`: the finite per-cell repair CAPACITY (number of
//!     rescue events available before the machinery is exhausted). `0.0` (default)
//!     ⇒ no rescues possible.
//!
//! The consumed budget is carried per cell on `CellState::escrt_budget_used`. The
//! decision each death-threshold crossing is [`escrt_rescue`]: draw a uniform roll
//! and rescue if it falls below the rate. The engine draws that roll ONLY when
//! `escrt_repair_rate > 0` and budget remains, so the default path makes no extra
//! RNG draw and stays byte-identical. A rescued cell does not die that step; its
//! defenses still run next step, so ESCRT buys time for GSH/GPX4 to recover (a
//! genuine rescue under transient stress) while sustained GPX4 inhibition still
//! kills once the budget is spent (a finite delay).
//!
//! ## Honesty / calibration
//!
//! The DIRECTION (more ESCRT repair ⇒ slower execution ⇒ more resistance;
//! CHMP5/CHMP6 loss ⇒ sensitization) is literature-anchored (Dai 2020). The rate
//! and budget are UNCALIBRATED placeholders; the result is the direction and the
//! finite-delay structure, not the numbers. Off-by-default
//! (`escrt_repair_rate == 0.0`) is byte-identical.

/// ESCRT membrane-repair rescue decision for one death-threshold crossing:
/// `roll < repair_rate`. `repair_rate` is the per-step rescue probability
/// (clamped to `[0, 1]`); `roll` is a uniform draw in `[0, 1)`. Returns `true`
/// when the cell is resealed (survives this step) and `false` when death proceeds.
/// `repair_rate == 0.0` ⇒ always `false` (the brake never fires).
pub fn escrt_rescue(repair_rate: f64, roll: f64) -> bool {
    roll < repair_rate.clamp(0.0, 1.0)
}

/// True when an ESCRT rescue attempt is even possible: the brake is enabled
/// (`repair_rate > 0`) and the cell has repair budget left
/// (`budget_used < budget`). The engine checks this BEFORE drawing the RNG roll,
/// so the default (`repair_rate == 0`) path draws nothing and is byte-identical.
pub fn escrt_can_attempt(repair_rate: f64, budget_used: f64, budget: f64) -> bool {
    repair_rate > 0.0 && budget_used < budget
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rescue_fires_below_rate_and_off_when_rate_zero() {
        // roll below the rate ⇒ rescued; at/above ⇒ death proceeds.
        assert!(escrt_rescue(0.5, 0.4));
        assert!(!escrt_rescue(0.5, 0.5));
        assert!(!escrt_rescue(0.5, 0.9));
        // rate 0 ⇒ never rescues, regardless of roll (the byte-identical default).
        assert!(!escrt_rescue(0.0, 0.0));
        assert!(!escrt_rescue(0.0, 0.999));
        // rate clamps to [0,1]; rate >= 1 always rescues for a valid roll < 1.
        assert!(escrt_rescue(1.0, 0.999));
        assert!(escrt_rescue(2.0, 0.5));
    }

    #[test]
    fn can_attempt_requires_enabled_and_budget() {
        // Disabled (rate 0) ⇒ never attempt (no RNG draw on the default path).
        assert!(!escrt_can_attempt(0.0, 0.0, 5.0));
        // Enabled + budget remaining ⇒ attempt.
        assert!(escrt_can_attempt(0.5, 0.0, 5.0));
        assert!(escrt_can_attempt(0.5, 4.0, 5.0));
        // Enabled but budget exhausted ⇒ no attempt (death proceeds).
        assert!(!escrt_can_attempt(0.5, 5.0, 5.0));
        assert!(!escrt_can_attempt(0.5, 6.0, 5.0));
        // Zero budget ⇒ never attempt even when enabled.
        assert!(!escrt_can_attempt(0.5, 0.0, 0.0));
    }
}
