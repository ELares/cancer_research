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
//! # What these four numbers can and cannot be told apart
//!
//! At a SINGLE time point they cannot. `effective_effectors` is
//! `infused · (T·I·A) · (1−r)^steps`, so at fixed `steps` any trafficking,
//! infiltration and activation with the same product are bit-identical, and
//! an exhaustion rate can be absorbed into activation exactly — verified in
//! `exhaustion_is_absorbable_into_activation_at_one_timepoint`, which
//! constructs the collapsed config and asserts the two agree to the bit.
//!
//! Only varying `steps` separates persistence from the barriers, and nothing
//! separates the three barriers from each other at all. So this module is
//! four numbers whose PRODUCT is one measurable quantity, and its value is
//! that the product decomposes into steps a reader can attack — not that the
//! decomposition is identifiable from an outcome. Fitting all four to one
//! endpoint would be fitting a single scalar four times.
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
///
/// `exhaustion_rate` is CLAMPED at both ends and both ends matter: a negative
/// rate would make `(1 - r) > 1` and this decay function would return more
/// effector function than was infused (`1.5^180` is about 5e31), which is the
/// same failure the `powi` wrap produced and is reached by a different route.
#[must_use = "the remaining function is the function's only output"]
pub fn persistence_factor(b: &AdoptiveBarriers, steps: u32) -> f64 {
    // `powi(steps as i32)` was WRONG and the bound test could not see it: the
    // cast wraps, so `u32::MAX` gave exponent -1 and this decay function
    // returned 2.0 -- more effector function than was infused -- while
    // `steps = 2^31` gave `inf`. `powf` takes the whole u32 range.
    (1.0 - b.exhaustion_rate.clamp(0.0, 1.0)).powf(f64::from(steps))
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

/// The most cells adoptive therapy can kill, whatever the dose.
///
/// This is where the ceiling is APPLIED, and it exists because for one review
/// round it was not: `antigen_ceiling` was a clamped getter that nothing read,
/// so the module documented a structural limit it did not have and a tenfold
/// infusion killed tenfold more cells against a 10%-antigen-positive tumour.
#[must_use = "the cap is the function's only output"]
pub fn max_killable(tumour_cells: f64, b: &AdoptiveBarriers) -> f64 {
    tumour_cells.max(0.0) * antigen_ceiling(b)
}

/// Kills after the antigen ceiling, given a raw effector-driven kill count.
///
/// A `min`, deliberately, and the contrast with [`crate::adc`] is the reason
/// the function exists at all. `adc::bystander_kill_fraction` PASSES the
/// antigen ceiling by killing antigen-negative neighbours with a diffusing
/// payload; a T cell cannot see them, so here the ceiling is a wall rather
/// than a coefficient. That distinction is not decorative: a multiplier can
/// be compensated by infusing more cells and a cap cannot, which is why
/// antigen escape ends the therapy instead of raising its dose.
#[must_use = "the kill count is the function's only output"]
pub fn barrier_limited_kills(raw_kills: f64, tumour_cells: f64, b: &AdoptiveBarriers) -> f64 {
    raw_kills.max(0.0).min(max_killable(tumour_cells, b))
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

    #[test]
    fn persistence_actually_reaches_the_effector_count_over_a_run() {
        // A reviewer deleted `persistence_factor` from `effective_effectors`
        // and all six tests stayed green: one ran on the default where
        // persistence is exactly 1, and another compared two configs sharing
        // an exhaustion rate and step count, so the factor cancelled on both
        // sides. Neither varied `steps`, which is the ONLY axis that can see
        // it. This one does.
        let b = AdoptiveBarriers::solid_tumour();
        let early = effective_effectors(1000.0, &b, 0);
        let late = effective_effectors(1000.0, &b, 500);
        assert!(
            late < early * 0.01,
            "500 steps of exhaustion barely moved the effector count: {early} -> {late}"
        );
        // And the shape is the decay, not merely "smaller": each step removes
        // the same FRACTION, so the ratio over equal spans is constant.
        let a = effective_effectors(1000.0, &b, 100) / effective_effectors(1000.0, &b, 50);
        let c = effective_effectors(1000.0, &b, 200) / effective_effectors(1000.0, &b, 150);
        assert!(
            (a - c).abs() < 1e-12,
            "exhaustion is not geometric: {a} vs {c}"
        );
    }

    #[test]
    fn exhaustion_is_absorbable_into_activation_at_one_timepoint() {
        // NOT a property to be proud of -- it is the identifiability limit,
        // pinned so the module cannot quietly start claiming four separable
        // axes. At fixed `steps` the four numbers are one scalar, and the
        // collapsed config agrees to the BIT.
        let b = AdoptiveBarriers::solid_tumour();
        let steps = 100;
        let collapsed = AdoptiveBarriers {
            activation: b.activation * persistence_factor(&b, steps),
            exhaustion_rate: 0.0,
            ..b
        };
        assert_eq!(
            effective_effectors(1000.0, &b, steps).to_bits(),
            effective_effectors(1000.0, &collapsed, steps).to_bits(),
            "the collapse is the point: only `steps` separates these"
        );
        // Varying steps is what breaks the tie, and it must.
        assert!(
            effective_effectors(1000.0, &b, 200)
                < effective_effectors(1000.0, &collapsed, 200) * 0.9,
            "nothing distinguishes persistence from activation at any step"
        );
    }

    #[test]
    fn the_solid_tumour_preset_is_three_real_barriers_not_one_catastrophe() {
        // A reviewer replaced the preset with 0.99/0.099/0.99 -- one collapsed
        // step and two wide open -- and every test passed, because the only
        // assertion on it was `delivery_efficiency < 0.1`, which 0.097 meets.
        // The corpus says the barriers are GENERAL rather than antigen-specific
        // (PMID 31848460), so a preset resting on a single step contradicts the
        // sentence the module is built from. No barrier may carry the collapse.
        let b = AdoptiveBarriers::solid_tumour();
        let each = [b.trafficking, b.infiltration, b.activation];
        // The criterion is a SHARE, not a band. A first version also required
        // each value in `0.2..=0.7`, which was arbitrary, lived in the same
        // file as the values it bounded, and rejected legitimate presets with
        // a false message: 0.1/0.15/0.2 is three general barriers by this
        // file's own measure and was reported as "a single catastrophe". The
        // share test alone does the work -- it correctly rejects 0.2/0.7/0.7,
        // a real catastrophe sitting INSIDE that band.
        let total = -delivery_efficiency(&b).ln();
        assert!(total > 0.0, "the preset presents no barrier at all");
        // The preset is the ONLY place the antigen ceiling and the exhaustion
        // rate can be exercised, so both must actually be set: a mutant moving
        // `antigen_positive_fraction` to 1.0 removed the ceiling from the whole
        // suite and nothing failed.
        assert!(
            (0.0..1.0).contains(&b.antigen_positive_fraction),
            "the solid-tumour preset has no antigen escape, so the ceiling is \
             unexercised everywhere"
        );
        assert!(
            b.exhaustion_rate > 0.0,
            "the solid-tumour preset has no exhaustion, so persistence is \
             unexercised everywhere"
        );
        for (name, v) in ["trafficking", "infiltration", "activation"]
            .iter()
            .zip(each)
        {
            assert!(
                (0.0..1.0).contains(&v),
                "{name} at {v} is not a fraction of the cells that got this far"
            );
            assert!(
                -v.ln() < 0.5 * total,
                "{name} at {v} carries more than half the collapse in log terms, so the \
                 preset rests on one step and contradicts the corpus sentence the module \
                 is built from (barriers GENERAL, not antigen-specific)"
            );
        }
    }

    #[test]
    fn the_antigen_ceiling_is_a_wall_a_larger_dose_cannot_climb() {
        // For one review round `antigen_ceiling` was a clamped getter nothing
        // read, so the module's structural claim was false in the most direct
        // way available: a 1000x infusion killed 1000x more cells against a
        // 10%-antigen-positive tumour. The cap is applied now, and this is the
        // test that would have caught it.
        let b = AdoptiveBarriers {
            antigen_positive_fraction: 0.1,
            ..AdoptiveBarriers::default()
        };
        let tumour = 20_000.0;
        let small = barrier_limited_kills(1_000.0, tumour, &b);
        let huge = barrier_limited_kills(1_000_000.0, tumour, &b);
        assert_eq!(
            huge.to_bits(),
            (tumour * 0.1).to_bits(),
            "the cap is not applied"
        );
        assert!(small < huge, "the cap fires below the ceiling too");
        // A thousandfold more cells buys nothing once the cap binds.
        assert_eq!(
            huge.to_bits(),
            barrier_limited_kills(1_000_000_000.0, tumour, &b).to_bits()
        );
        // The ADC contrast is the reason the wall is a `min`, and it is
        // asserted rather than described: at the SAME 10% antigen fraction a
        // cleavable-linker ADC kills strictly MORE than that fraction, because
        // its payload reaches the antigen-negative cells a T cell cannot see.
        let adc = crate::adc::AdcConfig {
            antigen_positive_fraction: 0.1,
            direct_kill_probability: 1.0,
            linker: crate::adc::Linker::Cleavable,
            payload_escape_fraction: 0.8,
            neighbours_in_reach: 0.5,
        };
        let past_the_ceiling = crate::adc::bystander_kill_fraction(&adc);
        assert!(
            past_the_ceiling > 0.0,
            "the ADC must pass the ceiling or the contrast this module draws is empty"
        );
        assert_eq!(
            huge / tumour,
            0.1,
            "adoptive therapy stops at the ceiling the ADC steps past"
        );
    }

    #[test]
    fn persistence_does_not_wrap_at_a_step_count_no_run_would_reach() {
        // `powi(steps as i32)` wrapped: `u32::MAX` gave exponent -1 and this
        // decay function returned 2.0, and `steps = 2^31` gave `inf`. A decay
        // returning more than it started with is the failure a bound test
        // exists for, and the bound test only probed `steps = 10`.
        let b = AdoptiveBarriers {
            exhaustion_rate: 0.5,
            ..AdoptiveBarriers::default()
        };
        for steps in [0u32, 1, 10, 1 << 31, u32::MAX] {
            let f = persistence_factor(&b, steps);
            assert!(
                (0.0..=1.0).contains(&f),
                "persistence_factor({steps}) = {f} is not a fraction"
            );
            assert!(
                effective_effectors(1000.0, &b, steps).is_finite(),
                "effective_effectors is not finite at {steps} steps"
            );
        }
        assert_eq!(persistence_factor(&b, 0).to_bits(), 1.0f64.to_bits());
    }

    #[test]
    fn the_serde_defaults_are_the_leukaemia_case_and_are_actually_reachable() {
        // `fn one()` could be changed to return 0.0 and nothing failed, because
        // nothing in the crate or its tests ever deserialised these barriers --
        // five `#[serde(default = "one")]` attributes and the function behind
        // them were dead. A config file omitting a field would then have
        // silently applied a total barrier instead of an open one.
        let from_empty: AdoptiveBarriers =
            serde_json::from_str("{}").expect("every field must have a default");
        assert_eq!(from_empty, AdoptiveBarriers::default());
        // And a partial config must fill only what it omits.
        let partial: AdoptiveBarriers = serde_json::from_str(r#"{"trafficking": 0.25}"#).unwrap();
        assert_eq!(partial.trafficking.to_bits(), 0.25f64.to_bits());
        assert_eq!(partial.infiltration.to_bits(), 1.0f64.to_bits());
        assert_eq!(
            partial.antigen_positive_fraction.to_bits(),
            1.0f64.to_bits()
        );
        assert_eq!(partial.exhaustion_rate.to_bits(), 0.0f64.to_bits());
    }

    #[test]
    fn the_caps_hold_at_inputs_no_caller_should_produce() {
        // Four mutants survived round 1 because the bounds test was never
        // extended to the two functions that round ADDED. Every clamp in
        // `max_killable` and `barrier_limited_kills` is exercised here, at the
        // value that makes its absence visible rather than a plausible one.
        let open = AdoptiveBarriers::default();

        // An out-of-range antigen fraction must not manufacture targets.
        let impossible = AdoptiveBarriers {
            antigen_positive_fraction: 5.0,
            ..open
        };
        assert_eq!(
            barrier_limited_kills(1000.0, 100.0, &impossible).to_bits(),
            100.0f64.to_bits(),
            "an antigen fraction above 1 killed more cells than exist"
        );
        assert_eq!(
            max_killable(100.0, &impossible).to_bits(),
            100.0f64.to_bits()
        );

        // Negative kills are not negative deaths, and a negative tumour is not
        // a negative cap.
        assert_eq!(
            barrier_limited_kills(-500.0, 100.0, &open).to_bits(),
            0.0f64.to_bits()
        );
        assert_eq!(
            barrier_limited_kills(10.0, -100.0, &open).to_bits(),
            0.0f64.to_bits()
        );
        assert_eq!(max_killable(-100.0, &open).to_bits(), 0.0f64.to_bits());

        // EVERY CLAMP, AT BOTH ENDS, ONE FACTOR AT A TIME. Nine mutants
        // survived round 2 and almost all of them for one reason: the hostile
        // config set `infiltration: -5.0`, which zeroes the whole product, so
        // a broken clamp on trafficking or activation was masked by a zero and
        // `>= 0.0` passed. A bound satisfied by zero tests nothing, so each
        // factor below is made hostile ALONE while the others stay open.
        for wild in [-5.0f64, 5.0] {
            for (name, b) in [
                (
                    "trafficking",
                    AdoptiveBarriers {
                        trafficking: wild,
                        ..open
                    },
                ),
                (
                    "infiltration",
                    AdoptiveBarriers {
                        infiltration: wild,
                        ..open
                    },
                ),
                (
                    "activation",
                    AdoptiveBarriers {
                        activation: wild,
                        ..open
                    },
                ),
            ] {
                let e = delivery_efficiency(&b);
                assert!(
                    (0.0..=1.0).contains(&e),
                    "{name} at {wild} gives delivery {e}, which is not a fraction"
                );
            }
            // A negative exhaustion rate must not manufacture effector
            // function: unclamped it gives 1.5^180, about 5e31.
            let ex = AdoptiveBarriers {
                exhaustion_rate: wild,
                ..open
            };
            for steps in [0u32, 11, 180] {
                let f = persistence_factor(&ex, steps);
                assert!(
                    (0.0..=1.0).contains(&f),
                    "exhaustion {wild} at {steps} steps gives {f}"
                );
            }
            // And a negative antigen fraction must not make a negative cap.
            let ag = AdoptiveBarriers {
                antigen_positive_fraction: wild,
                ..open
            };
            assert!((0.0..=1.0).contains(&antigen_ceiling(&ag)));
            assert!(max_killable(100.0, &ag) >= 0.0);
            assert!(barrier_limited_kills(50.0, 100.0, &ag) >= 0.0);
        }
        // `infused` is clamped too, and the old assertion could not see it:
        // it used a config whose delivery was already 0, so the product was
        // `-0.0` and `-0.0 >= 0.0` holds. Delivery is open here.
        assert_eq!(
            effective_effectors(-100.0, &open, 10).to_bits(),
            0.0f64.to_bits(),
            "a negative infusion produced a negative effector count"
        );

        // Total exhaustion at zero steps is no exhaustion, not NaN: a
        // reformulation as `(ln(1-r) * steps).exp()` returns NaN here, since
        // ln(0) is -inf and -inf * 0 is NaN. The wrap test probes steps = 0
        // only at a survivable rate, so it cannot see this corner.
        let dead = AdoptiveBarriers {
            exhaustion_rate: 1.0,
            ..open
        };
        assert_eq!(persistence_factor(&dead, 0).to_bits(), 1.0f64.to_bits());
        assert_eq!(persistence_factor(&dead, 1).to_bits(), 0.0f64.to_bits());
        assert!(effective_effectors(1000.0, &dead, 3).is_finite());
    }
}
