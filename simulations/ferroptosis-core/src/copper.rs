//! Copper-ionophore / ferroptosis-cuproptosis crosstalk (#485).
//!
//! Copper ionophores (elesclomol, disulfiram) raise mitochondrial copper, and
//! the resulting copper overload simultaneously enables ferroptosis AND
//! cuproptosis: Cu(II) depletes glutathione (GSH binds copper) and drives the
//! autophagic / proteasomal degradation of GPX4, while also promoting SLC7A11
//! degradation. Elesclomol specifically degrades the ATP7A copper-efflux
//! transporter, trapping copper (Gao et al., Mol Oncol 2021, PMID 34390123); FIN + ionophore
//! synergy, PMID 37214358; review Signal Transduct Target Ther 2025). Sensitivity
//! is gated by the copper-efflux transporter ATP7A/B: a tumor that exports copper
//! efficiently keeps its intracellular copper low and RESISTS the ionophore.
//!
//! This is a metal-2 axis that bridges the repo's iron-driven ferroptosis engine
//! to a copper death mode, and it lets the suite model elesclomol + RSL3
//! combinations (the cuproptosis already scaffolded in the taxonomy but never
//! simulated).
//!
//! ## Design: pure helpers, consumer mutates
//!
//! Like [`crate::ifngamma`], this module is pure: it returns per-step
//! retention multipliers and the consumer (e.g. sim-tme-3d) multiplies them into
//! the cell GSH / GPX4 pools each step. The core single-cell engine
//! ([`crate::biochem::sim_cell_step`]) does NOT read a `CopperConfig`, so the
//! core stays byte-identical; the consumer composes the helpers around the step.
//!
//! ## Identity default
//!
//! [`CopperConfig::disabled`] applies no effect (both depletion rates `0.0`), so
//! the retention multipliers are exactly `1.0` and a consumer stays
//! byte-identical. [`CopperConfig::literature`] is an uncalibrated placeholder
//! (an elesclomol-like ionophore in an efflux-incompetent tumor).

/// Copper-ionophore configuration (#485). [`disabled`](CopperConfig::disabled) is
/// the identity (no effect, byte-identical).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CopperConfig {
    /// Per-step fractional GSH depletion at full intracellular copper (copper
    /// binds and consumes glutathione). In `[0, 1]`; `0.0` disables this arm.
    pub cu_gsh_depletion: f64,
    /// Per-step fractional GPX4 degradation at full intracellular copper
    /// (autophagic/proteasomal GPX4 loss under copper overload). In `[0, 1]`;
    /// `0.0` disables this arm.
    pub cu_gpx4_degradation: f64,
    /// ATP7A/B copper-efflux competence, in `[0, 1]`. `0.0` = no efflux (full
    /// intracellular copper, full effect); `1.0` = complete efflux (copper
    /// exported, the cell RESISTS the ionophore). The effective depletion each
    /// arm applies is scaled by `(1 - atp7b_efflux)`.
    pub atp7b_efflux: f64,
}

impl CopperConfig {
    /// Identity: no copper effect (both depletion rates `0.0`), so the retention
    /// multipliers are exactly `1.0` and a consumer is byte-identical.
    pub fn disabled() -> Self {
        CopperConfig {
            cu_gsh_depletion: 0.0,
            cu_gpx4_degradation: 0.0,
            atp7b_efflux: 0.0,
        }
    }

    /// Uncalibrated placeholder: an elesclomol-like ionophore in an
    /// efflux-incompetent (ATP7B-low) tumor, so the copper overload depletes both
    /// GSH and GPX4 each step. Magnitudes are placeholders; the direction
    /// (copper ionophore ⇒ GSH/GPX4 loss ⇒ more ferroptosis; ATP7B efflux ⇒
    /// resistance) is the result.
    pub fn literature() -> Self {
        CopperConfig {
            cu_gsh_depletion: 0.03,
            cu_gpx4_degradation: 0.03,
            atp7b_efflux: 0.0,
        }
    }

    /// True when the config applies no effect (both depletion rates `0.0`), so a
    /// consumer can skip the copper path and stay byte-identical. ATP7B efflux
    /// alone (with zero depletion rates) is still identity.
    pub fn is_disabled(&self) -> bool {
        self.cu_gsh_depletion == 0.0 && self.cu_gpx4_degradation == 0.0
    }
}

/// Per-step GSH-retention multiplier under copper overload (#485):
/// `1 - cu_gsh_depletion * (1 - atp7b_efflux)`, clamped to `[0, 1]`. A consumer
/// multiplies the cell GSH pool by this each step. Exactly `1.0` when
/// `cu_gsh_depletion == 0` (identity / byte-identical); a higher ATP7B efflux
/// RAISES it back toward `1.0` (copper exported ⇒ less GSH consumed).
#[must_use]
pub fn gsh_retention(cfg: &CopperConfig) -> f64 {
    let efflux = cfg.atp7b_efflux.clamp(0.0, 1.0);
    (1.0 - cfg.cu_gsh_depletion * (1.0 - efflux)).clamp(0.0, 1.0)
}

/// Per-step GPX4-retention multiplier under copper overload (#485):
/// `1 - cu_gpx4_degradation * (1 - atp7b_efflux)`, clamped to `[0, 1]`. A
/// consumer multiplies the cell GPX4 pool by this each step. Exactly `1.0` when
/// `cu_gpx4_degradation == 0` (identity); ATP7B efflux raises it back toward
/// `1.0`.
#[must_use]
pub fn gpx4_retention(cfg: &CopperConfig) -> f64 {
    let efflux = cfg.atp7b_efflux.clamp(0.0, 1.0);
    (1.0 - cfg.cu_gpx4_degradation * (1.0 - efflux)).clamp(0.0, 1.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn disabled_is_identity() {
        let c = CopperConfig::disabled();
        assert!(c.is_disabled());
        assert_eq!(gsh_retention(&c), 1.0);
        assert_eq!(gpx4_retention(&c), 1.0);
        assert!(!CopperConfig::literature().is_disabled());
        // ATP7B efflux alone (no depletion) is still identity.
        let efflux_only = CopperConfig {
            atp7b_efflux: 0.8,
            ..CopperConfig::disabled()
        };
        assert!(efflux_only.is_disabled());
        assert_eq!(gsh_retention(&efflux_only), 1.0);
    }

    #[test]
    fn copper_depletes_gsh_and_gpx4_and_atp7b_efflux_protects() {
        // An efflux-incompetent tumor: both pools are depleted each step.
        let c = CopperConfig {
            cu_gsh_depletion: 0.5,
            cu_gpx4_degradation: 0.4,
            atp7b_efflux: 0.0,
        };
        assert!((gsh_retention(&c) - 0.5).abs() < 1e-12);
        assert!((gpx4_retention(&c) - 0.6).abs() < 1e-12);
        // ATP7B efflux RAISES retention back toward 1.0 (copper exported).
        let with_efflux = CopperConfig {
            atp7b_efflux: 0.5,
            ..c
        };
        assert!(gsh_retention(&with_efflux) > gsh_retention(&c));
        assert!(gpx4_retention(&with_efflux) > gpx4_retention(&c));
        // Full efflux ⇒ exactly 1.0 (the tumor fully resists the ionophore).
        let full_efflux = CopperConfig {
            atp7b_efflux: 1.0,
            ..c
        };
        assert_eq!(gsh_retention(&full_efflux), 1.0);
        assert_eq!(gpx4_retention(&full_efflux), 1.0);
    }

    #[test]
    fn retention_is_bounded_and_clamps() {
        // efflux clamps to [0,1]; depletion > 1 floors retention at 0.
        let strong = CopperConfig {
            cu_gsh_depletion: 2.0,
            cu_gpx4_degradation: 2.0,
            atp7b_efflux: -1.0,
        };
        assert_eq!(gsh_retention(&strong), 0.0);
        assert_eq!(gpx4_retention(&strong), 0.0);
    }
}
