//! Adoptive cell therapy: the three barriers between an infusion and a kill.
//!
//! # Why this module exists
//!
//! `immune::adoptive_transfer_kills` models what a redirected T cell does once
//! it is touching a tumour cell, and `analysis/modality-calibration.md` fits it
//! to the published B-ALL complete-remission band of 70–94%. Both are correct
//! and neither is the interesting case, because the same source that supplies
//! that band says the result has NOT transferred to solid tumours.
//!
//! A model that reproduces the leukaemia number and has nothing to say about
//! the solid-tumour failure is describing the indication these therapies were
//! approved for rather than the setting this engine simulates. This module is
//! the difference.
//!
//! # The three barriers, named by one sentence in the corpus
//!
//! PMID 31848460 (`corpus/by-pmid/31848460.md`) lists them in order: *"Key
//! challenges relating to CAR T cells include severe toxicities, restricted
//! **trafficking to, infiltration into and activation within** tumours,
//! suboptimal **persistence** in vivo"*. Those are three sequential filters and
//! a fourth failure mode, and they are modelled as such:
//!
//! 1. **Trafficking** — the fraction of infused cells that reach the tumour at
//!    all. A leukaemia is in the blood the cells are infused into; a solid
//!    tumour is not.
//! 2. **Infiltration** — of those that arrive, the fraction that get past the
//!    rim. The same source records that the *"collective lack of efficacy …
//!    targeting several different solid tumour antigens suggests the existence
//!    of general barriers"* — general, meaning not antigen-specific.
//! 3. **Activation** — of those inside, the fraction not suppressed on
//!    arrival.
//!
//! # Why sequential filters and not a single efficiency
//!
//! Because they MULTIPLY, and that is the whole point. Three barriers at 50%
//! each leave an eighth, which is how a therapy that cures a blood cancer can
//! fail a solid one without any single step looking catastrophic. Collapsing
//! them into one number would fit the same endpoints and lose the reason —
//! and the reason is what a reader needs in order to know which step to
//! attack.
//!
//! [`delivery_efficiency`] is the product, and it is the quantity to read.
//!
//! # The asymmetry with the ADC, which is worth stating
//!
//! Antigen escape defeats both modalities (`corpus/by-pmid/31947597.md`:
//! escape and downregulation "limiting therapy effectiveness by leading to
//! antigen negative relapse"). [`crate::adc`] models a mechanism that answers
//! it — a cleavable linker kills antigen-negative neighbours. **Adoptive
//! therapy has no equivalent.** A T cell that cannot see a cell cannot kill
//! it, and there is no diffusing payload. That asymmetry is expressed here as
//! a hard ceiling rather than a parameter, because it is structural.
//!
//! # What is NOT modelled
//!
//! Toxicity, which the same corpus sentence names FIRST. Cytokine release
//! syndrome is the dose-limiting problem in the indication where this therapy
//! works, so a module that models only efficacy is describing the easy half.
//! Nothing here should be read as a dose recommendation.

use serde::{Deserialize, Serialize};

/// Barriers between an infusion and a tumour cell.
///
/// [`Default`] is the LEUKAEMIA case — every barrier fully open — because
/// that is the setting the published band is measured in, and a default that
/// silently applied solid-tumour barriers would make
/// `modality-calibration.md`'s fit describe something other than what it
/// claims.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct AdoptiveBarriers {
    /// Fraction of infused cells reaching the tumour.
    #[serde(default = "one")]
    pub trafficking: f64,
    /// Of those, the fraction penetrating past the rim.
    #[serde(default = "one")]
    pub infiltration: f64,
    /// Of those, the fraction not suppressed on arrival.
    #[serde(default = "one")]
    pub activation: f64,
    /// Per-step loss of effector function from persistent signalling.
    #[serde(default)]
    pub exhaustion_rate: f64,
    /// Fraction of the tumour expressing the target antigen. A T cell cannot
    /// kill what it cannot see, and unlike an ADC it has no payload that
    /// diffuses to neighbours.
    #[serde(default = "one")]
    pub antigen_positive_fraction: f64,
}

fn one() -> f64 {
    1.0
}

impl Default for AdoptiveBarriers {
    fn default() -> Self {
        Self {
            trafficking: 1.0,
            infiltration: 1.0,
            activation: 1.0,
            exhaustion_rate: 0.0,
            antigen_positive_fraction: 1.0,
        }
    }
}

impl AdoptiveBarriers {
    /// The barriers a solid tumour presents, as documented directions.
    ///
    /// **Every value here is a placeholder and the CALIBRATION_STATUS row says
    /// so.** What the corpus establishes is that all three barriers are real,
    /// that they are general rather than antigen-specific, and that their
    /// combined effect is large enough to erase a therapy that cures a blood
    /// cancer. It does not establish which of the three dominates, and this
    /// preset must not be read as claiming it does.
    #[must_use]
    pub fn solid_tumour() -> Self {
        Self {
            trafficking: 0.3,
            infiltration: 0.4,
            activation: 0.5,
            exhaustion_rate: 0.02,
            antigen_positive_fraction: 0.8,
        }
    }
}

/// The product of the three sequential barriers.
///
/// The quantity to read, and the reason they are modelled separately: three
/// barriers at 50% leave an eighth, which is how a therapy can fail without
/// any single step looking catastrophic.
#[must_use = "the efficiency is the function's only output"]
pub fn delivery_efficiency(b: &AdoptiveBarriers) -> f64 {
    b.trafficking.clamp(0.0, 1.0) * b.infiltration.clamp(0.0, 1.0) * b.activation.clamp(0.0, 1.0)
}

/// Effector function remaining after `steps` of persistent signalling.
///
/// The fourth failure mode — "suboptimal persistence in vivo" — and separate
/// from the three barriers because it acts on cells that already ARRIVED. A
/// model folding it into activation would predict that improving trafficking
/// fixes persistence, which it does not.
#[must_use = "the remaining function is the function's only output"]
pub fn persistence_factor(b: &AdoptiveBarriers, steps: u32) -> f64 {
    (1.0 - b.exhaustion_rate.clamp(0.0, 1.0)).powi(steps as i32)
}

/// Effectors that are both present and functional at `steps`.
#[must_use = "the effector count is the function's only output"]
pub fn effective_effectors(infused: f64, b: &AdoptiveBarriers, steps: u32) -> f64 {
    infused.max(0.0) * delivery_efficiency(b) * persistence_factor(b, steps)
}

/// The ceiling antigen expression imposes, which adoptive therapy cannot pass.
///
/// Stated as its own function because the CONTRAST is the point: [`crate::adc`]
/// has `bystander_kill_fraction`, which passes exactly this ceiling by killing
/// antigen-negative neighbours with a diffusing payload. A T cell has no such
/// mechanism, so for this modality the ceiling is structural.
#[must_use = "the ceiling is the function's only output"]
pub fn antigen_ceiling(b: &AdoptiveBarriers) -> f64 {
    b.antigen_positive_fraction.clamp(0.0, 1.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_default_is_the_leukaemia_case_with_every_barrier_open() {
        // This matters beyond tidiness: `modality-calibration.md` fits the
        // published B-ALL band with suppression at zero, and a default that
        // silently applied solid-tumour barriers would make that fit describe
        // something other than what it claims.
        let d = AdoptiveBarriers::default();
        assert_eq!(delivery_efficiency(&d).to_bits(), 1.0_f64.to_bits());
        assert_eq!(persistence_factor(&d, 500).to_bits(), 1.0_f64.to_bits());
        assert_eq!(
            effective_effectors(1000.0, &d, 200).to_bits(),
            1000.0_f64.to_bits()
        );
        assert_eq!(antigen_ceiling(&d).to_bits(), 1.0_f64.to_bits());
    }

    /// The barriers MULTIPLY, which is how a therapy fails without any single
    /// step looking catastrophic.
    #[test]
    fn three_survivable_barriers_compound_into_an_unsurvivable_one() {
        let each_half = AdoptiveBarriers {
            trafficking: 0.5,
            infiltration: 0.5,
            activation: 0.5,
            ..AdoptiveBarriers::default()
        };
        let eff = delivery_efficiency(&each_half);
        assert!(
            (eff - 0.125).abs() < 1e-12,
            "three 50% barriers should leave an eighth, got {eff}"
        );
        // No SINGLE barrier at 50% comes close to that.
        for b in [
            AdoptiveBarriers {
                trafficking: 0.5,
                ..AdoptiveBarriers::default()
            },
            AdoptiveBarriers {
                infiltration: 0.5,
                ..AdoptiveBarriers::default()
            },
            AdoptiveBarriers {
                activation: 0.5,
                ..AdoptiveBarriers::default()
            },
        ] {
            assert!(
                delivery_efficiency(&b) > eff * 3.0,
                "a single barrier is nearly as bad as three, so modelling \
                 them separately buys nothing"
            );
        }
        // And the solid-tumour preset is far below the leukaemia case.
        let solid = delivery_efficiency(&AdoptiveBarriers::solid_tumour());
        assert!(
            solid < 0.1,
            "the solid-tumour preset leaves {solid} of the infusion, which is \
             not the collapse the corpus describes"
        );
        assert!(solid > 0.0, "the preset is not a total block");
    }

    /// Persistence acts on cells that already ARRIVED, so it must be separable
    /// from the barriers -- improving trafficking must not fix exhaustion.
    #[test]
    fn persistence_is_independent_of_the_delivery_barriers() {
        let exhausting = AdoptiveBarriers {
            exhaustion_rate: 0.05,
            ..AdoptiveBarriers::solid_tumour()
        };
        // Perfect delivery, same exhaustion: persistence is unchanged.
        let perfect_delivery = AdoptiveBarriers {
            trafficking: 1.0,
            infiltration: 1.0,
            activation: 1.0,
            ..exhausting
        };
        assert_eq!(
            persistence_factor(&exhausting, 100).to_bits(),
            persistence_factor(&perfect_delivery, 100).to_bits(),
            "fixing the delivery barriers changed persistence, so the model \
             would predict that better trafficking cures exhaustion"
        );
        // But the EFFECTOR count does improve, because delivery is a real
        // separate factor.
        assert!(
            effective_effectors(1000.0, &perfect_delivery, 100)
                > effective_effectors(1000.0, &exhausting, 100)
        );
        // Monotone decay, and it must actually bite over a run.
        let mut prev = 1.0;
        for s in [0_u32, 10, 50, 200] {
            let f = persistence_factor(&exhausting, s);
            assert!(f <= prev, "not monotone at {s}");
            prev = f;
        }
        assert!(persistence_factor(&exhausting, 200) < 0.01);
        // Zero rate is the exact identity, so an unconfigured run is unmoved.
        let no_exhaustion = AdoptiveBarriers {
            exhaustion_rate: 0.0,
            ..exhausting
        };
        assert_eq!(
            persistence_factor(&no_exhaustion, 10_000).to_bits(),
            1.0_f64.to_bits()
        );
    }

    /// The asymmetry with the ADC, asserted against the other module rather
    /// than described in prose.
    #[test]
    fn adoptive_therapy_has_no_answer_to_antigen_escape_and_the_adc_does() {
        use crate::adc::{bystander_kill_fraction, total_kill_fraction, AdcConfig, Linker};

        let escaped = 0.5; // half the tumour has lost the antigen
        let t_cells = AdoptiveBarriers {
            antigen_positive_fraction: escaped,
            ..AdoptiveBarriers::default()
        };
        // A T-cell ceiling IS the antigen fraction, with no mechanism past it.
        assert!((antigen_ceiling(&t_cells) - escaped).abs() < 1e-12);

        // The ADC, on the same tumour, passes it.
        let adc = AdcConfig {
            antigen_positive_fraction: escaped,
            direct_kill_probability: 1.0,
            linker: Linker::Cleavable,
            payload_escape_fraction: 0.5,
            neighbours_in_reach: 0.6,
        };
        assert!(bystander_kill_fraction(&adc) > 0.0);
        assert!(
            total_kill_fraction(&adc) > antigen_ceiling(&t_cells),
            "the ADC no longer passes the antigen ceiling, so the asymmetry \
             this module documents is gone and both paragraphs need \
             re-deriving"
        );
    }

    #[test]
    fn everything_is_bounded_at_hostile_parameters() {
        let wild = AdoptiveBarriers {
            trafficking: 5.0,
            infiltration: -5.0,
            activation: 5.0,
            exhaustion_rate: 5.0,
            antigen_positive_fraction: 5.0,
        };
        let eff = delivery_efficiency(&wild);
        assert!((0.0..=1.0).contains(&eff), "{eff}");
        assert_eq!(
            eff.to_bits(),
            0.0_f64.to_bits(),
            "a negative barrier must clamp to zero"
        );
        assert!((0.0..=1.0).contains(&persistence_factor(&wild, 10)));
        assert!((0.0..=1.0).contains(&antigen_ceiling(&wild)));
        assert!(effective_effectors(-100.0, &wild, 10) >= 0.0);
    }
}
