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
/// Bystander killing is drawn from the cells the DIRECT term left alive, which
/// is not the same set as the antigen-negative cells and was described as if it
/// were. With `direct_kill_probability < 1` some antigen-POSITIVE cells also
/// survive, so the surviving pool is `1 - dying` and only `(1 - antigen) /
/// (1 - dying)` of it is antigen-negative. Use [`bystander_kill_on_negative`]
/// for the antigen-negative part, which is the quantity an escape experiment
/// measures.
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

/// The part of the bystander kill that lands on ANTIGEN-NEGATIVE cells.
///
/// This exists because the preregistration needed a quantity an experiment can
/// produce and [`bystander_kill_fraction`] is not one. Divided by the
/// antigen-negative pool it was reported as a "share" of **216%** -- impossible
/// for a share -- because the bystander term is bounded by every surviving
/// cell, antigen-positive ones included, while the sentence beside it claimed
/// the antigen-negative pool.
///
/// Payload that escapes a dying cell does not know its neighbour's antigen
/// status, so it is apportioned across the surviving pool in proportion:
/// `bystander * (1 - antigen) / (1 - dying)`. That is bounded by the
/// antigen-negative pool by construction, which is what makes P10's
/// falsification threshold scoreable.
#[must_use = "the kill fraction is the function's only output"]
pub fn bystander_kill_on_negative(cfg: &AdcConfig) -> f64 {
    let dying = direct_kill_fraction(cfg);
    let surviving = (1.0 - dying).max(0.0);
    if surviving <= 0.0 {
        return 0.0;
    }
    let negative = (1.0 - cfg.antigen_positive_fraction.clamp(0.0, 1.0)).max(0.0);
    bystander_kill_fraction(cfg) * (negative / surviving)
}

/// [`bystander_kill_on_negative`] as a share of the antigen-negative pool.
///
/// In `[0, 1]` by construction, and the quantity P10 is registered on.
#[must_use = "the share is the function's only output"]
pub fn negative_pool_reached(cfg: &AdcConfig) -> f64 {
    let negative = (1.0 - cfg.antigen_positive_fraction.clamp(0.0, 1.0)).max(0.0);
    if negative <= 0.0 {
        return 0.0;
    }
    // NO CLAMP. A clamp here made the guard policing this bound satisfy
    // itself: reverting the numerator to the raw bystander term -- literally
    // the 216% this function exists to retract -- passed every test, because
    // the clamp capped it and the assertion then checked the clamp. The
    // apportionment is bounded by construction, so a value above one means the
    // construction is wrong and must surface rather than be flattened.
    bystander_kill_on_negative(cfg) / negative
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

// ── Drug loading, and the barrier that made it a trade-off ───────────────
//
// The module above treats an ADC as a payload with a delivery problem. The
// engineering question it cannot ask is the one that actually decides the
// design: HOW MUCH payload per antibody?
//
// More looks strictly better -- each antibody that arrives carries more drug --
// and it is not, because loading changes the antibody. A more heavily
// conjugated antibody is cleared faster, so fewer arrive. Two monotonic
// effects in opposite directions, which is the shape this chapter keeps
// finding, and here the trade-off has been MEASURED: Hamblett 2004 (PMID
// 15501986) reports that an eight-loaded conjugate clears three times faster
// than a four-loaded one and five times faster than a two-loaded one, and
// that at equal antibody dose the four-loaded conjugate performed comparably
// to the eight.
//
// The second effect is the binding-site barrier: a higher-affinity antibody
// binds the first antigen it meets and never reaches the cell behind it, so
// affinity buys retention at the cost of penetration. `scripts/validate_
// penetration.py` deliberately did NOT add this for small molecules, on the
// grounds that it is an antibody phenomenon and physically weak for a small
// drug. This is the case that refusal excluded.

/// Clearance of an ADC relative to a two-loaded conjugate, by drug-antibody
/// ratio.
///
/// Anchored on Hamblett 2004 (PMID 15501986), which measured RATIOS rather
/// than absolute rates: an eight-loaded conjugate clears three times faster
/// than a four-loaded one and five times faster than a two-loaded one. Taking
/// DAR 2 as the reference, that fixes `c(2) = 1`, `c(4) = 5/3` and `c(8) = 5`.
///
/// **A single power law cannot pass through both measurements, and finding
/// that out is worth more than the curve.** A law fitted to `c(8)/c(2) = 5`
/// predicts `c(8)/c(4) = 2.24` against a measured 3. The two ratios together
/// say the penalty ACCELERATES: clearance rises more slowly than linearly
/// between DAR 2 and 4, and faster than linearly between 4 and 8. That
/// acceleration is what produces the interior optimum below -- so a curve
/// chosen for smoothness would have removed the finding.
///
/// Interpolated log-linearly WITHIN each measured segment and extrapolated
/// with the nearest segment's exponent outside `[2, 8]`, where nothing is
/// measured and the result should not be trusted.
#[must_use]
pub fn clearance_multiplier(dar: f64) -> f64 {
    let d = dar.max(0.1);
    // c(2) = 1, c(4) = 5/3, c(8) = 5.
    let (lo, hi, c_lo, c_hi) = if d <= 4.0 {
        (2.0f64, 4.0f64, 1.0f64, 5.0f64 / 3.0)
    } else {
        (4.0f64, 8.0f64, 5.0f64 / 3.0, 5.0f64)
    };
    let exponent = (c_hi / c_lo).ln() / (hi / lo).ln();
    c_lo * (d / lo).powf(exponent)
}

/// The drug-antibody ratios Hamblett 2004 measured, and the clearance ratios
/// it reported between them.
///
/// Public so a validation can assert the crate reproduces them rather than
/// restating them in a second place.
pub const HAMBLETT_DAR_POINTS: [(f64, f64); 3] = [(2.0, 1.0), (4.0, 5.0 / 3.0), (8.0, 5.0)];

/// Payload delivered to the tumour per unit of antibody dose.
///
/// `dar / clearance(dar)`: more drug on each antibody, fewer antibodies
/// surviving to arrive. **The interior optimum is not put in.** Both factors
/// are monotonic in DAR; the ratio is not, and where it peaks follows from the
/// measured clearance ratios rather than from a choice made here.
#[must_use]
pub fn payload_delivered_per_dose(dar: f64) -> f64 {
    let d = dar.max(0.0);
    d / clearance_multiplier(d).max(f64::MIN_POSITIVE)
}

/// In-vitro potency relative to a two-loaded conjugate.
///
/// Rises monotonically with loading, because in a dish there is no clearance:
/// every conjugate reaches every cell. That is why the in-vitro ordering and
/// the in-vivo ordering DISAGREE, and reproducing the disagreement is the
/// point -- it is the same in-vitro-to-in-vivo gap this book argues about
/// everywhere else, arriving in the one place where the field measured both
/// halves.
#[must_use]
pub fn in_vitro_potency(dar: f64) -> f64 {
    (dar.max(0.0) / 2.0).max(0.0)
}

/// Penetration depth relative to a low-affinity antibody, by binding affinity.
///
/// THE BINDING-SITE BARRIER. A higher-affinity antibody is captured by the
/// first antigen it meets, so the front of the tumour consumes it and the
/// cells behind never see it. Penetration therefore FALLS as affinity rises,
/// which is the opposite of the intuition that a better binder is a better
/// drug (Thurber 2008, PMID 18541331).
///
/// UNCALIBRATED. The form is `1 / (1 + affinity/scale)`, which has the right
/// direction and the right limits and is not fitted to a measured depth.
#[must_use]
pub fn penetration_vs_affinity(relative_affinity: f64, scale: f64) -> f64 {
    let a = relative_affinity.max(0.0);
    (1.0 / (1.0 + a / scale.max(f64::MIN_POSITIVE))).clamp(0.0, 1.0)
}

/// Tumour cells reached, combining loading and the binding-site barrier.
///
/// The product of what arrives (loading against clearance) and how far it
/// spreads (affinity against the barrier). Both halves have interior optima
/// and they are not at the same place, which is why an ADC is an engineering
/// compromise rather than a maximisation.
#[must_use]
pub fn delivered_reach(dar: f64, relative_affinity: f64, barrier_scale: f64) -> f64 {
    payload_delivered_per_dose(dar) * penetration_vs_affinity(relative_affinity, barrier_scale)
}

/// The drug-antibody ratio that maximises delivered payload, by scan.
///
/// Returns `(dar, delivered)`. Scanned rather than solved because the point is
/// the shape, and because the same scan produces the figure.
#[must_use]
pub fn optimal_dar(max_dar: f64) -> (f64, f64) {
    let mut best = (0.0, f64::NEG_INFINITY);
    let n = 400;
    for i in 1..=n {
        let d = max_dar * f64::from(i) / f64::from(n);
        let v = payload_delivered_per_dose(d);
        if v > best.1 {
            best = (d, v);
        }
    }
    best
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

    /// THE BYSTANDER MECHANISM IS SELF-LIMITING, and this test used to claim
    /// the opposite.
    ///
    /// It was `antigen_loss_hurts_the_non_bystander_arm_more`, and it asserted
    /// that the relative advantage `b/a` is NON-DECREASING as antigen is lost
    /// -- with a comment saying that is what separates "addresses escape" from
    /// "just adds kill". Measured at this module's own `heterogeneous()`
    /// config the ratio is EXACTLY CONSTANT at `1 + escape*reach = 1.30` for
    /// every antigen fraction, because `bystander = dying * escape * reach`
    /// and `dying = apf * kill`, so both arms scale with `apf` and the ratio
    /// cancels. `>=` accepts a constant, and a constant ratio IS "just adds
    /// kill" -- the reading the guard existed to rule out. The name was false
    /// too: both arms are hurt in exact proportion, so neither is hurt more.
    ///
    /// What the model actually says is more interesting and is the direction
    /// worth preregistering: the payload comes from cells that TOOK UP the
    /// ADC, i.e. antigen-positive ones, and antigen escape removes exactly
    /// that population. So the bystander effect does not answer escape better
    /// as escape worsens -- it is starved by it. The share of the
    /// antigen-NEGATIVE pool a cleavable linker reaches falls from 77.1% at
    /// 0.9 positive to 2.6% at 0.1 positive. (This docstring said 90% and 1%,
    /// which were the values the RETRACTED formula produced -- the retraction
    /// reached the code and not the paragraph performing it.)
    #[test]
    fn the_bystander_effect_is_starved_by_the_escape_it_answers() {
        // THE QUANTITY IS THE ANTIGEN-NEGATIVE SHARE, and getting that wrong
        // is what this test is for. A previous version divided the whole
        // bystander kill by the antigen-negative pool and published 216% as a
        // "share" -- impossible -- because the bystander term is bounded by
        // every surviving cell, antigen-positive ones included.
        let mut prev = f64::INFINITY;
        let mut shares = Vec::new();
        for &positive in &[0.9_f64, 0.6, 0.3, 0.1] {
            let cfg = AdcConfig {
                antigen_positive_fraction: positive,
                linker: Linker::Cleavable,
                ..heterogeneous()
            };
            let share = negative_pool_reached(&cfg);
            assert!(
                (0.0..=1.0).contains(&share),
                "a SHARE of the antigen-negative pool came out at {share} for \
                 {positive} positive, which no experiment can produce"
            );
            assert!(
                share < prev,
                "the bystander arm reached a LARGER share of the \
                 antigen-negative pool at {positive} positive ({share:.4} vs \
                 {prev:.4}); the payload comes from antigen-POSITIVE cells, so \
                 losing them cannot help"
            );
            prev = share;
            shares.push(share);
        }
        // AND THE DECLINE MUST BE STEEP, not merely non-increasing. The guard
        // this replaces accepted a CONSTANT -- `>=` on a ratio that cancels --
        // and that vacuity is the whole reason the direction had to be
        // rederived. A model where bystander kill were proportional to the
        // negative pool would hold `share` flat and must fail here.
        assert!(
            shares[0] > shares[3] * 5.0,
            "the reach falls only from {:.4} to {:.4}; a near-flat decline is \
             the shape the retracted guard could not distinguish from a real \
             one",
            shares[0],
            shares[3]
        );
        // The unbounded quantity still exists and is still NOT a share: this
        // pins the distinction the retraction turns on.
        let dense = AdcConfig {
            antigen_positive_fraction: 0.9,
            linker: Linker::Cleavable,
            ..heterogeneous()
        };
        assert!(
            bystander_kill_fraction(&dense) > 1.0 - dense.antigen_positive_fraction,
            "the raw bystander term no longer exceeds the antigen-negative \
             pool, so the two quantities can no longer be confused and this \
             test's subject is gone"
        );
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
    // ── Drug loading, and the barrier ────────────────────────────────────

    #[test]
    fn the_clearance_curve_passes_through_both_measured_ratios() {
        // Hamblett 2004 (PMID 15501986) measured two ratios, not a curve. A
        // single power law fitted to one of them misses the other by a third,
        // which is why this is piecewise -- and why the acceleration it
        // encodes is a property of the DATA rather than of the fit.
        for (dar, expected) in HAMBLETT_DAR_POINTS {
            let got = clearance_multiplier(dar);
            assert!(
                (got - expected).abs() < 1e-9,
                "clearance at DAR {dar} is {got}, measured {expected}"
            );
        }
        let c8 = clearance_multiplier(8.0);
        assert!((c8 / clearance_multiplier(4.0) - 3.0).abs() < 1e-9);
        assert!((c8 / clearance_multiplier(2.0) - 5.0).abs() < 1e-9);
    }

    #[test]
    fn the_penalty_for_loading_accelerates() {
        // The shape that produces the optimum. Between DAR 2 and 4 clearance
        // rises MORE SLOWLY than linearly; between 4 and 8, faster. A single
        // exponent cannot do both, and a model that used one would have no
        // interior optimum at all.
        let low = clearance_multiplier(4.0) / clearance_multiplier(2.0);
        let high = clearance_multiplier(8.0) / clearance_multiplier(4.0);
        assert!(
            high > low,
            "the penalty did not accelerate: {low} then {high}"
        );
        assert!(
            low < 2.0,
            "doubling the load below DAR 4 more than doubled clearance"
        );
        assert!(
            high > 2.0,
            "doubling the load above DAR 4 did not more than double it"
        );
    }

    #[test]
    fn the_optimal_loading_falls_out_of_the_measured_clearance() {
        // THE FINDING, and it is not tuned: both factors in
        // `payload_delivered_per_dose` are monotonic in DAR, the ratio is not,
        // and where it peaks follows from ratios somebody else measured.
        let (dar, _) = optimal_dar(16.0);
        assert!((dar - 4.0).abs() < 0.2,
                "the optimum moved to DAR {dar}; it should sit where the                  clearance penalty starts accelerating");
        // And the consequence Hamblett reports: twice the payload per antibody
        // does NOT double what arrives.
        let ratio = payload_delivered_per_dose(8.0) / payload_delivered_per_dose(4.0);
        assert!(ratio < 1.0,
                "an eight-loaded conjugate delivered MORE than a four-loaded                  one ({ratio}), which is the opposite of the measurement");
        assert!(ratio > 0.5, "it should be comparable, not halved: {ratio}");
    }

    #[test]
    fn the_in_vitro_and_in_vivo_orderings_disagree() {
        // The reason this arm is worth a section in THIS book. In a dish
        // there is no clearance, so potency rises monotonically with loading
        // and the eight-loaded conjugate wins. In an animal the four-loaded
        // one matches it at equal antibody dose. Same molecule, opposite
        // ordering, and the model reproduces both halves.
        assert!(in_vitro_potency(8.0) > in_vitro_potency(4.0));
        assert!(in_vitro_potency(4.0) > in_vitro_potency(2.0));
        assert!(payload_delivered_per_dose(8.0) < payload_delivered_per_dose(4.0));
    }

    #[test]
    fn higher_affinity_costs_penetration() {
        // The binding-site barrier (Thurber 2008, PMID 18541331). Deliberately
        // NOT added for small molecules by `scripts/validate_penetration.py`,
        // on the grounds that it is an antibody phenomenon; this is the case
        // that refusal excluded.
        let weak = penetration_vs_affinity(0.1, 1.0);
        let strong = penetration_vs_affinity(10.0, 1.0);
        assert!(
            strong < weak,
            "affinity did not cost penetration: {weak} {strong}"
        );
        assert!(
            penetration_vs_affinity(0.0, 1.0) > 0.99,
            "a non-binding antibody should penetrate freely"
        );
        assert!(strong > 0.0, "penetration should fall, not vanish");
        // And the combined reach is worse than either factor alone suggests.
        let reach = delivered_reach(4.0, 10.0, 1.0);
        assert!(reach < payload_delivered_per_dose(4.0));
    }
}
