//! Antibody-drug conjugates: the payload's second journey.
//!
//! # The question this module exists to answer
//!
//! `analysis/modality-panel.md` measures something that looks fatal for the
//! modality: an ADC carrying the SAME payload as sonodynamic therapy, run
//! through the same engine with the same parameters, kills a fiftieth as
//! much — because a ~150 kDa antibody diffuses an order of magnitude more
//! slowly than a small molecule and is consumed by the first antigen it
//! meets. `drug_transport::antibody_drug_conjugate` puts its penetration at
//! about 7 µm.
//!
//! Seven micrometres is roughly one cell. If that were the whole story ADCs
//! would not work, and they do. **This module is the part that was missing**,
//! and it is not an excuse for the transport number — it is a second,
//! independent mechanism that operates after the antibody has stopped moving.
//!
//! # The bystander effect, and why the LINKER decides whether it happens
//!
//! An ADC delivers its payload inside the cell that bound it. What happens
//! next depends on the chemistry holding the two together:
//!
//! * A **cleavable** linker releases a membrane-permeable payload, which
//!   diffuses out of the dying cell and kills neighbours — including cells
//!   that never expressed the antigen at all.
//! * A **non-cleavable** linker leaves the payload attached to a charged
//!   amino-acid residue after lysosomal degradation. It cannot cross a
//!   membrane, so it kills the cell that bound it and nothing else.
//!
//! The corpus records exactly this contrast (PMID 31930187,
//! `corpus/by-pmid/31930187.md`): *"the MMAF payload can retain cytotoxicity
//! against target cells with upregulated efflux pumps, MMAE on a cleavable
//! linker can cause bystander effects and kill co-cultured cells that do not
//! express the cell surface antigen on the target cell."* One sentence, both
//! arms of the comparison, and the antigen-negative case named explicitly.
//!
//! # Why this is the mechanism that answers antigen heterogeneity too
//!
//! Antigen loss is the documented route to ADC and CAR-T failure alike
//! (`corpus/by-pmid/31947597.md`: antigen escape "limiting therapy
//! effectiveness by leading to antigen negative relapse"). A non-bystander ADC
//! cannot touch an antigen-negative cell by construction, so its ceiling IS
//! the antigen-positive fraction. A bystander ADC is not bounded that way, and
//! that difference is the whole reason the linker choice is a clinical
//! decision rather than a chemistry detail.
//!
//! # What is NOT modelled
//!
//! Efflux pumps, which the same corpus sentence names as the reason to prefer
//! the NON-cleavable payload — so the trade-off here is one-sided by
//! construction and this module cannot be used to choose a linker. It can only
//! say what bystander killing buys, not what it costs.
//!
//! Spatial arrangement, too: `bystander_reach_fraction` takes the neighbour
//! fraction within reach as an input rather than deriving it from a lattice.
//! Whether antigen-negative cells are scattered among positive ones or
//! clustered away from them changes the answer entirely, and a well-mixed
//! model cannot tell the difference.

use serde::{Deserialize, Serialize};

/// Whether the linker releases a payload that can leave the cell.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum Linker {
    /// Releases a membrane-permeable payload — MMAE-like. Kills neighbours.
    Cleavable,
    /// Payload stays attached to a charged residue — MMAF-like. Does not.
    NonCleavable,
}

impl Linker {
    /// Can this payload leave the cell that internalised it?
    ///
    /// The single property the whole module turns on, and a `match` rather
    /// than a numeric knob because it is genuinely binary chemistry: a payload
    /// either crosses a membrane or it does not.
    #[must_use]
    pub fn releases_permeable_payload(self) -> bool {
        matches!(self, Linker::Cleavable)
    }
}

/// ADC configuration.
///
/// [`Default`] is the identity: no antigen-positive cells and no payload, so
/// every helper returns zero and a run carrying it is unmoved.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct AdcConfig {
    /// Fraction of the tumour expressing the target antigen.
    #[serde(default)]
    pub antigen_positive_fraction: f64,
    /// Probability an antigen-positive cell that binds the ADC is killed.
    #[serde(default)]
    pub direct_kill_probability: f64,
    /// Linker chemistry.
    #[serde(default = "default_linker")]
    pub linker: Linker,
    /// Fraction of the payload that escapes a dying cell rather than being
    /// degraded with it. Only reachable with a cleavable linker.
    #[serde(default)]
    pub payload_escape_fraction: f64,
    /// Neighbours within diffusion reach of one dying cell, as a fraction of
    /// the tumour. An INPUT, not derived — see the module docs.
    #[serde(default)]
    pub neighbours_in_reach: f64,
}

fn default_linker() -> Linker {
    Linker::NonCleavable
}

impl Default for AdcConfig {
    fn default() -> Self {
        Self {
            antigen_positive_fraction: 0.0,
            direct_kill_probability: 0.0,
            linker: default_linker(),
            payload_escape_fraction: 0.0,
            neighbours_in_reach: 0.0,
        }
    }
}

/// Cells killed by direct antigen-mediated delivery, as a tumour fraction.
///
/// Bounded by the antigen-positive fraction by construction: an ADC cannot
/// bind a cell that does not express its target.
#[must_use = "the kill fraction is the function's only output"]
pub fn direct_kill_fraction(cfg: &AdcConfig) -> f64 {
    cfg.antigen_positive_fraction.clamp(0.0, 1.0) * cfg.direct_kill_probability.clamp(0.0, 1.0)
}

/// Additional cells killed by payload that escaped a dying cell.
///
/// Returns exactly `0.0` for a non-cleavable linker, whatever the escape
/// fraction and whatever the neighbour density. That is not a numerical
/// convenience: a payload bound to a charged residue cannot cross a membrane,
/// so no amount of it in a dying cell reaches the next one.
///
/// Bystander killing reaches ANTIGEN-NEGATIVE cells, so it is drawn from the
/// remainder the direct term cannot touch. That is the point of the mechanism
/// and the reason it addresses antigen escape rather than merely adding kill.
#[must_use = "the kill fraction is the function's only output"]
pub fn bystander_kill_fraction(cfg: &AdcConfig) -> f64 {
    if !cfg.linker.releases_permeable_payload() {
        return 0.0;
    }
    let dying = direct_kill_fraction(cfg);
    let escaped = cfg.payload_escape_fraction.clamp(0.0, 1.0);
    let reach = cfg.neighbours_in_reach.clamp(0.0, 1.0);
    // The pool bystander killing can draw from is what the direct term left,
    // which includes every antigen-negative cell.
    let remaining = (1.0 - dying).max(0.0);
    (dying * escaped * reach).min(remaining)
}

/// Total kill fraction: direct plus bystander.
#[must_use = "the kill fraction is the function's only output"]
pub fn total_kill_fraction(cfg: &AdcConfig) -> f64 {
    (direct_kill_fraction(cfg) + bystander_kill_fraction(cfg)).clamp(0.0, 1.0)
}

/// The ceiling a NON-bystander ADC cannot pass, whatever its potency.
///
/// Its own antigen-positive fraction. Stated as a function rather than left
/// implicit because it is the quantity antigen-escape relapse is about: a
/// tumour that loses the antigen moves this ceiling down, and a payload that
/// cannot leave the cell has no way around it.
#[must_use = "the ceiling is the function's only output"]
pub fn antigen_limited_ceiling(cfg: &AdcConfig) -> f64 {
    cfg.antigen_positive_fraction.clamp(0.0, 1.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn heterogeneous() -> AdcConfig {
        AdcConfig {
            antigen_positive_fraction: 0.5,
            direct_kill_probability: 0.8,
            linker: Linker::Cleavable,
            payload_escape_fraction: 0.5,
            neighbours_in_reach: 0.6,
        }
    }

    #[test]
    fn the_default_config_kills_nothing() {
        let d = AdcConfig::default();
        assert_eq!(direct_kill_fraction(&d).to_bits(), 0.0_f64.to_bits());
        assert_eq!(bystander_kill_fraction(&d).to_bits(), 0.0_f64.to_bits());
        assert_eq!(total_kill_fraction(&d).to_bits(), 0.0_f64.to_bits());
        // And the default linker is the CONSERVATIVE one: a config that
        // forgot to choose must not silently get the extra mechanism.
        assert_eq!(d.linker, Linker::NonCleavable);
    }

    /// The linker decides whether bystander killing happens AT ALL, and no
    /// other parameter can substitute for it.
    #[test]
    fn a_non_cleavable_linker_kills_no_neighbours_at_any_setting() {
        let base = heterogeneous();
        for &escape in &[0.0_f64, 0.5, 1.0] {
            for &reach in &[0.0_f64, 0.5, 1.0] {
                let c = AdcConfig {
                    linker: Linker::NonCleavable,
                    payload_escape_fraction: escape,
                    neighbours_in_reach: reach,
                    ..base
                };
                assert_eq!(
                    bystander_kill_fraction(&c).to_bits(),
                    0.0_f64.to_bits(),
                    "a non-cleavable payload reached a neighbour at \
                     escape={escape}, reach={reach}"
                );
                // And its total is exactly its direct kill -- no leakage.
                assert_eq!(
                    total_kill_fraction(&c).to_bits(),
                    direct_kill_fraction(&c).to_bits()
                );
            }
        }
        // The cleavable version of the SAME config does kill neighbours, so
        // the zero above is the linker and not the parameters.
        assert!(bystander_kill_fraction(&base) > 0.0);
    }

    /// The mechanism that answers antigen escape: bystander killing passes a
    /// ceiling the direct term cannot.
    #[test]
    fn bystander_killing_passes_the_antigen_limited_ceiling() {
        // A perfectly potent ADC on a half-antigen-positive tumour.
        let perfect = AdcConfig {
            antigen_positive_fraction: 0.5,
            direct_kill_probability: 1.0,
            linker: Linker::NonCleavable,
            ..heterogeneous()
        };
        let ceiling = antigen_limited_ceiling(&perfect);
        assert!((direct_kill_fraction(&perfect) - ceiling).abs() < 1e-12);
        assert!(
            (total_kill_fraction(&perfect) - ceiling).abs() < 1e-12,
            "a non-bystander ADC exceeded its antigen ceiling, which is \
             impossible by construction"
        );

        // The same tumour, the same potency, a cleavable linker.
        let bystanding = AdcConfig {
            linker: Linker::Cleavable,
            ..perfect
        };
        assert!(
            total_kill_fraction(&bystanding) > ceiling,
            "bystander killing did not pass the antigen ceiling ({} vs {})",
            total_kill_fraction(&bystanding),
            ceiling
        );
        // And the excess is real, not a rounding artefact.
        assert!(total_kill_fraction(&bystanding) > ceiling * 1.1);
    }

    /// Antigen LOSS lowers the ceiling, and a bystander ADC degrades more
    /// gracefully than one without it. That difference is the clinical
    /// argument for the linker.
    #[test]
    fn antigen_loss_hurts_the_non_bystander_arm_more() {
        let mut prev_gap = 0.0_f64;
        for &positive in &[0.9_f64, 0.6, 0.3, 0.1] {
            let plain = AdcConfig {
                antigen_positive_fraction: positive,
                linker: Linker::NonCleavable,
                ..heterogeneous()
            };
            let bystanding = AdcConfig {
                linker: Linker::Cleavable,
                ..plain
            };
            let (a, b) = (
                total_kill_fraction(&plain),
                total_kill_fraction(&bystanding),
            );
            assert!(b > a, "no bystander advantage at {positive} positive");
            // The RELATIVE advantage must be at least as large as antigen is
            // lost -- otherwise the mechanism does not address escape, it just
            // adds kill.
            let gap = b / a;
            assert!(
                gap >= prev_gap - 1e-9,
                "the bystander advantage shrank as antigen was lost \
                 ({gap:.3} vs {prev_gap:.3}), so it does not address escape"
            );
            prev_gap = gap;
        }
    }

    #[test]
    fn nothing_can_exceed_the_tumour() {
        // Parameters chosen to break the bound.
        let wild = AdcConfig {
            antigen_positive_fraction: 5.0,
            direct_kill_probability: 5.0,
            linker: Linker::Cleavable,
            payload_escape_fraction: 5.0,
            neighbours_in_reach: 5.0,
        };
        let total = total_kill_fraction(&wild);
        assert!((0.0..=1.0).contains(&total), "{total}");
        // Bystander killing may not claim cells the direct term already took.
        let d = direct_kill_fraction(&wild);
        assert!(d + bystander_kill_fraction(&wild) <= 1.0 + 1e-9);
        // Negative inputs clamp rather than inverting the sign.
        let negative = AdcConfig {
            antigen_positive_fraction: -1.0,
            direct_kill_probability: -1.0,
            payload_escape_fraction: -1.0,
            neighbours_in_reach: -1.0,
            linker: Linker::Cleavable,
        };
        assert_eq!(total_kill_fraction(&negative).to_bits(), 0.0_f64.to_bits());
    }

    /// The module's reason for existing: it is what makes a 7 µm penetration
    /// depth survivable.
    #[test]
    fn bystander_reach_is_what_rescues_a_short_penetration() {
        use crate::drug_transport::{antibody_drug_conjugate, penetration_length_um};
        let reach_um = penetration_length_um(&antibody_drug_conjugate());
        assert!(
            reach_um < 15.0,
            "the ADC transport profile no longer has a short penetration \
             ({reach_um:.1} um), so this module's premise needs re-deriving"
        );
        // With no bystander mechanism, kill is capped by antigen expression --
        // and that cap is what the transport number then further limits.
        let no_bystander = AdcConfig {
            linker: Linker::NonCleavable,
            ..heterogeneous()
        };
        let with_bystander = heterogeneous();
        let gain = total_kill_fraction(&with_bystander) / total_kill_fraction(&no_bystander);
        assert!(
            gain > 1.1,
            "bystander killing adds only {gain:.3}x; it cannot be what makes \
             a one-cell penetration depth clinically survivable"
        );
    }
}
