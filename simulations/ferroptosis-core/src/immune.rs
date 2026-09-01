//! ICD-immune cascade model.
//!
//! Models: ferroptotic death → DAMP release → DC maturation → T cell priming → tumor killing.
//! Deliberately simple (3-stage linear cascade) to show qualitative ICD differences
//! between physical modalities and pharmacologic inducers.

use serde::Serialize;

use crate::params::ImmuneParams;

/// Result of the immune cascade for one treatment condition.
#[derive(Clone, Debug, Serialize)]
pub struct ImmuneResult {
    /// Total DAMP signal from all dead cells.
    pub total_damps: f64,
    /// Average DAMP per dead cell (higher = more immunogenic death).
    pub damp_per_dead_cell: f64,
    /// Number of dead cells that contributed DAMPs.
    pub n_dead_cells: usize,
    /// Fraction of DCs activated by DAMP signal.
    pub dc_activation_fraction: f64,
    /// Number of mature DCs.
    pub mature_dcs: f64,
    /// Number of primed T cells.
    pub primed_tcells: f64,
    /// Tumor cells killed by immune response.
    pub immune_kills: f64,
    /// Whether anti-PD-1 was applied.
    pub with_anti_pd1: bool,
}

/// Calculate total DAMP release from a set of dead cells.
///
/// DAMP is proportional to lipid peroxidation at death. Key biological insight:
/// - SDT/PDT kill via runaway LP cascade → LP at death is MUCH higher than threshold
/// - RSL3 kills via slow GPX4 inhibition → LP at death is NEAR threshold
/// Therefore SDT/PDT-killed cells release MORE DAMPs per cell.
///
/// Ref: Krysko et al., Nat Rev Cancer 2012 (ICD markers)
///      Berezhnoy et al., PLoS Comput Biol 2020 (Boolean ICD model)
pub fn calculate_damp_release(dead_cell_lps: &[f64], params: &ImmuneParams) -> (f64, f64) {
    if dead_cell_lps.is_empty() {
        return (0.0, 0.0);
    }
    let total: f64 = dead_cell_lps.iter().map(|lp| lp * params.damp_per_lp).sum();
    let per_cell = total / dead_cell_lps.len() as f64;
    (total, per_cell)
}

/// Published durable-response rates for the oncolytic virus T-VEC, treated
/// and control (`oncolytic-virus`: 5,006 census articles, 201 trials).
///
/// 16% against 2% (PMID 27298410, `corpus/by-pmid/27298410.md`): "With
/// talimogene laherparepvec, the primary end point of durable response rate
/// (DRR; continuous response lasting >= 6 months) was significantly higher
/// (16% v 2%...)".
///
/// **The CONTROL arm is stored with it deliberately.** A 16% durable response
/// quoted alone reads as a modest result; quoted against 2% it is an eightfold
/// ratio. Neither number means much without the other, and a constant holding
/// only the treated arm invites the reader to supply the wrong baseline.
pub const TVEC_DURABLE_RESPONSE: (f64, f64) = (0.16, 0.02);

/// Oncolytic-virus lysis and the immunogenic death it produces.
///
/// The modality is TWO effects and the second is the one that matters
/// clinically: the virus replicates in and lyses tumour cells, and that lysis
/// is IMMUNOGENIC, converting a local infection into a systemic anti-tumour
/// response. T-VEC's own trial reports responses at UNINJECTED lesions, which
/// direct lysis cannot explain.
///
/// So this returns both, and a consumer feeds `lysed` into
/// [`immune_cascade`]'s `dead_cell_lps` exactly as it would feed ferroptotic
/// deaths — the ICD chain is shared because the biology is shared. That is
/// also why this belongs in `immune.rs`: an oncolytic virus is an
/// immunotherapy that arrives as an infection.
///
/// `infected_fraction = 0.0` returns `(0.0, 0.0)`, so an unconfigured run is
/// unmoved.
///
/// ## What is NOT modelled, and it is most of the virology
///
/// Replication kinetics, antiviral immunity clearing the virus before it
/// spreads, and the interferon response that makes some tumours permissive and
/// others resistant — all absent. `infected_fraction` is an input, not a
/// result, which means this cannot answer "will the virus spread?", only "if
/// it spreads this far, what follows?".
#[must_use = "both outputs carry the modality's two effects"]
pub fn oncolytic_lysis(
    total_tumor_cells: usize,
    infected_fraction: f64,
    lysis_efficiency: f64,
    immunogenicity: f64,
) -> (f64, f64) {
    let infected = infected_fraction.clamp(0.0, 1.0) * total_tumor_cells as f64;
    let lysed = infected * lysis_efficiency.clamp(0.0, 1.0);
    // The per-cell "quality of death" the ICD chain reads, on the same scale
    // as `lp_at_grace_end`: viral lysis is highly immunogenic, so this is a
    // fraction of the maximum rather than a small perturbation of it.
    let damp_per_cell = immunogenicity.clamp(0.0, 1.0);
    (lysed, if lysed > 0.0 { damp_per_cell } else { 0.0 })
}

/// Published complete-remission band for CD19 CAR-T in B-ALL./// Published complete-remission band for CD19 CAR-T in B-ALL.
///
/// 70–94% (PMID 32607912, `corpus/by-pmid/32607912.md`): "Targeting of the
/// CD19 antigen using CD19-specific CAR-T cells ... with complete remission
/// rates of 70–94% seen in some clinical trials."
///
/// **The same source carries the number that keeps this honest**, and it is
/// quoted in [`adoptive_transfer_kills`]'s docs rather than left out: "the use
/// of CAR-T cells in solid tumours has been less successful", and "in 30–50%
/// of cases, the response was not durable". A model reproducing only the
/// headline band would be describing the indication these therapies were
/// approved for, not the setting this engine simulates.
pub const CART_B_ALL_CR_BAND: (f64, f64) = (0.70, 0.94);

/// How the effector T cells got there — and the two differ in a way that
/// matters more than the kill term does.
///
/// Both bypass dendritic-cell priming, which is why they share
/// [`adoptive_transfer_kills`]. What separates them is PERSISTENCE, and it is
/// the clinical difference rather than a modelling convenience:
///
/// * **CAR-T** cells are transferred once and then EXPAND and persist —
///   `tisagenlecleucel` is a single infusion. The engineered population is
///   self-renewing, so effector numbers grow after delivery.
/// * A **bispecific engager** is a drug, not a cell. It redirects the
///   patient's own resident T cells only while it is present, so its effect
///   tracks the dose schedule and stops when the drug clears.
///
/// The model expresses that as a persistence multiplier per step rather than
/// as two kill formulas, because the KILLING is the same event — a T cell
/// lysing a target — and only the supply differs.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize)]
pub enum EffectorSource {
    /// Transferred, self-renewing (`car-t`: 15,358 census articles).
    CarT,
    /// Redirected while the drug is present (`bispecific-antibody`: 3,462
    /// census articles, and the joint-highest trial share in the taxonomy at
    /// 9.4%).
    BispecificEngager,
}

impl EffectorSource {
    /// Effector population after `steps`, from a starting count.
    ///
    /// CAR-T expands geometrically toward a cap; a bispecific decays with the
    /// drug. `expansion_rate` and `clearance_rate` are BOTH uncalibrated and
    /// the CALIBRATION_STATUS row says so — what is modelled here is the
    /// DIRECTION, which is documented and clinically unambiguous, not the
    /// magnitude.
    #[must_use = "the effector count is the function's only output"]
    pub fn effectors_after(self, initial: f64, steps: u32, rate: f64, cap: f64) -> f64 {
        let r = rate.max(0.0);
        match self {
            // Logistic-style approach to a cap: expansion is not unbounded,
            // and a model that let it be would predict cures the literature
            // does not report.
            EffectorSource::CarT => {
                let mut n = initial.max(0.0);
                for _ in 0..steps {
                    n += r * n * (1.0 - n / cap.max(f64::MIN_POSITIVE));
                    n = n.clamp(0.0, cap.max(0.0));
                }
                n
            }
            // Exponential clearance: the redirecting drug washes out.
            EffectorSource::BispecificEngager => {
                initial.max(0.0) * (1.0 - r.min(1.0)).powi(steps as i32)
            }
        }
    }
}

/// Kills from ADOPTIVELY TRANSFERRED or REDIRECTED effector T cells.
///
/// Covers two of the taxonomy's absent mechanisms at once, because in this
/// engine they differ in how the effector arrives and not in what it does:
/// CAR-T (15,358 census articles) transfers pre-armed T cells, and bispecific
/// antibodies (3,462, and the joint-highest trial share at 9.4%) redirect
/// resident ones. Both **bypass the DC-priming step entirely**, which is the
/// whole point of the modality — antigen presentation is exactly what they do
/// not wait for.
///
/// So this is NOT [`immune_cascade`] with a different constant. That function
/// gates on `n_dead` and on DAMP quality, and a redirected T cell needs
/// neither. `effector_cells` enters the kill term directly:
///
/// `kills = effector_cells · kill_rate · (1 − brake) · (1 − suppression)`
///
/// **The brake still applies**, and that is a claim rather than a convenience:
/// CAR-T cells express PD-1 and are suppressed in solid tumours, which is the
/// documented reason the B-ALL result has not transferred. A model exempting
/// them would predict solid-tumour efficacy the literature does not report.
///
/// `suppression` is where the TME enters — the same Treg/MDSC and exhaustion
/// machinery [`crate::immune_spatial`] already carries — so the difference
/// between the leukaemia and solid-tumour settings is a parameter here rather
/// than two models.
///
/// Returns `0.0` for zero effector cells, so an unconfigured run is unmoved.
#[must_use = "the kill count is the function's only output"]
pub fn adoptive_transfer_kills(
    effector_cells: f64,
    total_tumor_cells: usize,
    params: &ImmuneParams,
    suppression: f64,
    with_anti_pd1: bool,
) -> f64 {
    let effectors = effector_cells.max(0.0);
    if effectors == 0.0 {
        return 0.0;
    }
    let effective_brake = if with_anti_pd1 {
        params.pd1_brake * (1.0 - params.anti_pd1_efficacy)
    } else {
        params.pd1_brake
    };
    let raw = effectors
        * params.tcell_kill_rate
        * (1.0 - effective_brake).clamp(0.0, 1.0)
        * (1.0 - suppression.clamp(0.0, 1.0));
    raw.min(total_tumor_cells as f64).max(0.0)
}

/// Macrophage phagocytosis after CD47/SIRPα blockade — the "don't eat me" axis.
///
/// The taxonomy's `phagocytosis-checkpoint` mechanism (918 census articles).
/// CD47 on the tumour binds SIRPα on macrophages and dendritic cells and
/// SUPPRESSES engulfment (PMID 30320184, `corpus/by-pmid/30320184.md`: "The
/// inhibitory effect of CD47 on phagocytosis is mediated by its binding to
/// signal-regulatory protein α (SIRPα), which is expressed on macrophages and
/// DCs"). Blocking it releases the brake.
///
/// Deliberately the same SHAPE as the T-cell brake — a fraction of a rate
/// removed by a drug — because that is what the biology is, and giving it a
/// different form would assert a distinction nobody measured. What differs is
/// the effector: macrophages engulf, they do not lyse, so this is a separate
/// count rather than an addend to the T-cell kills.
///
/// `cd47_expression = 1.0` with no blockade returns exactly `0.0`: a fully
/// protected tumour is not eaten at all.
#[must_use = "the phagocytosis count is the function's only output"]
pub fn phagocytosis_kills(
    macrophages: f64,
    total_tumor_cells: usize,
    engulf_rate: f64,
    cd47_expression: f64,
    blockade_efficacy: f64,
) -> f64 {
    let residual_protection = (cd47_expression.clamp(0.0, 1.0)
        * (1.0 - blockade_efficacy.clamp(0.0, 1.0)))
    .clamp(0.0, 1.0);
    let raw = macrophages.max(0.0) * engulf_rate.max(0.0) * (1.0 - residual_protection);
    raw.min(total_tumor_cells as f64).max(0.0)
}

/// Run the DC maturation → T cell priming → killing cascade.
///
/// Model stages:
/// 1. DAMPs activate DCs: saturating Michaelis-Menten: DAMP / (DAMP + Kd)
/// 2. Activated DCs mature with probability dc_maturation_rate
/// 3. Each mature DC primes tcell_priming_rate T cells
/// 4. Each T cell kills tcell_kill_rate tumor cells
/// 5. PD-1 brake suppresses a fraction of T-cell killing
/// 6. Anti-PD-1 removes a fraction of the brake
///
/// This is a coarse-grained model showing QUALITATIVE differences, not absolute numbers.
pub fn immune_cascade(
    dead_cell_lps: &[f64],
    total_tumor_cells: usize,
    params: &ImmuneParams,
    with_anti_pd1: bool,
) -> ImmuneResult {
    let (total_damps, damp_per_dead_cell) = calculate_damp_release(dead_cell_lps, params);
    let n_dead = dead_cell_lps.len();

    // DC activation: saturating response to DAMP *per dead cell* (quality of death).
    // Using per-cell average rather than total prevents the number of dead cells from
    // dominating the activation signal, allowing SDT vs RSL3 quality differences to show.
    let dc_activation_fraction =
        damp_per_dead_cell / (damp_per_dead_cell + params.dc_activation_kd);

    // Mature DCs: activation quality × maturation rate × number of
    // antigen-presenting deaths.
    //
    // PLUS a ferroptosis-INDEPENDENT term (#728). Every activation path in
    // this engine used to be gated on `n_dead`, so with no ferroptotic death
    // the product was zero and anti-PD-1 -- which enters only through the
    // brake -- multiplied zero at every setting. Real tumours present antigen
    // without being killed first, and that is why checkpoint blockade is a
    // monotherapy. `baseline_antigenicity = 0` (the default) makes this term
    // vanish and reproduces the prior behaviour exactly.
    let baseline_presenting = params.baseline_antigenicity.max(0.0) * total_tumor_cells as f64;
    let mature_dcs = dc_activation_fraction * params.dc_maturation_rate * n_dead as f64
        + params.dc_maturation_rate * baseline_presenting;

    // T cell priming
    let primed_tcells = mature_dcs * params.tcell_priming_rate;

    // T cell killing (with PD-1 brake)
    let effective_brake = if with_anti_pd1 {
        params.pd1_brake * (1.0 - params.anti_pd1_efficacy)
    } else {
        params.pd1_brake
    };
    let kill_efficiency = 1.0 - effective_brake;
    let raw_kills = primed_tcells * params.tcell_kill_rate * kill_efficiency;

    // Cap at remaining alive tumor cells
    let remaining = total_tumor_cells.saturating_sub(n_dead) as f64;
    let immune_kills = raw_kills.min(remaining);

    ImmuneResult {
        total_damps,
        damp_per_dead_cell,
        n_dead_cells: n_dead,
        dc_activation_fraction,
        mature_dcs,
        primed_tcells,
        immune_kills,
        with_anti_pd1,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // The 2D ICD cascade drives sim-icd and sim-combo but had no coverage
    // (#295), while its 3D successor `immune_spatial` is well-tested. These pin
    // the DAMP arithmetic and the cascade's qualitative invariants.

    #[test]
    fn damp_release_empty_is_zero() {
        let p = ImmuneParams::default();
        assert_eq!(calculate_damp_release(&[], &p), (0.0, 0.0));
    }

    #[test]
    fn damp_release_total_and_per_cell() {
        // total = Σ lp·damp_per_lp; per_cell = total / n. Default damp_per_lp = 1.0.
        let p = ImmuneParams::default();
        let (total, per_cell) = calculate_damp_release(&[2.0, 4.0], &p);
        assert!((total - 6.0 * p.damp_per_lp).abs() < 1e-12, "total={total}");
        assert!(
            (per_cell - 3.0 * p.damp_per_lp).abs() < 1e-12,
            "per_cell={per_cell}"
        );
    }

    #[test]
    fn empty_cascade_produces_no_kills() {
        let p = ImmuneParams::default();
        let r = immune_cascade(&[], 1000, &p, false);
        assert_eq!(r.total_damps, 0.0);
        assert_eq!(r.n_dead_cells, 0);
        assert_eq!(r.immune_kills, 0.0);
    }

    /// No ferroptotic death means no kills AT EVERY BLOCKADE SETTING.
    ///
    /// This is the claim `analysis/modality-coverage.md` rests on when it
    /// files immunotherapy as a MODIFIER rather than a treatment: checkpoint
    /// blockade cannot be the thing under test in this engine, because it
    /// enters only through the brake, and the brake multiplies a term that is
    /// already zero when nothing died by ferroptosis.
    ///
    /// Swept rather than asserted at one setting.
    /// `empty_cascade_produces_no_kills` above tests `with_anti_pd1 = false`
    /// only, so it could not have seen a blockade path that manufactured
    /// kills of its own.
    #[test]
    fn no_ferroptotic_death_means_no_kills_at_any_blockade_setting() {
        for &efficacy in &[0.0, 0.25, 0.5, 0.75, 1.0] {
            let p = ImmuneParams {
                anti_pd1_efficacy: efficacy,
                ..ImmuneParams::default()
            };
            for &with in &[false, true] {
                let r = immune_cascade(&[], 10_000, &p, with);
                assert_eq!(
                    r.immune_kills, 0.0,
                    "anti_pd1={efficacy}, with_anti_pd1={with} produced kills \
                     with no ferroptotic death"
                );
                assert_eq!(r.total_damps, 0.0);
            }
        }
        // And with NO brake at all -- the most favourable case blockade could
        // ever reach -- still zero.
        let p = ImmuneParams {
            pd1_brake: 0.0,
            ..ImmuneParams::default()
        };
        assert_eq!(immune_cascade(&[], 10_000, &p, true).immune_kills, 0.0);
    }

    /// END-TO-END: no ferroptotic death anywhere in the chain means no immune
    /// kills, at every blockade setting.
    ///
    /// **This test exists because three rounds of static scanning could not
    /// decide the claim it makes.** `analysis/modality-coverage.md` argues
    /// that checkpoint blockade cannot be a treatment arm in this engine, and
    /// a text scan over source lines cannot establish that: every regex
    /// written for it was defeated by an ordinary Rust idiom that put the
    /// mutation and the field name on different lines. A behaviour test can,
    /// for the path it executes.
    ///
    /// It runs the REAL chain -- generate cells, run the ferroptosis engine,
    /// collect lipid peroxidation at death, feed the immune cascade -- with
    /// death made impossible, and asserts zero kills however hard the
    /// blockade is pushed. Any antigen source reachable from this path, in
    /// any idiom, makes it fail.
    ///
    /// Paired with a POSITIVE case, so it cannot pass by the chain being
    /// broken: with death possible, the same chain must produce kills, and
    /// they must respond to blockade. A test that only asserts zero is
    /// satisfied by an engine that does nothing.
    #[test]
    fn no_ferroptotic_death_end_to_end_means_no_kills_at_any_blockade() {
        use crate::biochem::sim_cell;
        use crate::cell::{gen_cell, Phenotype, Treatment};
        use crate::params::Params;
        use rand::{rngs::StdRng, SeedableRng};

        const N: usize = 400;

        fn run(
            params: &Params,
            tx: Treatment,
            immune: &ImmuneParams,
            blockade: bool,
        ) -> (usize, f64) {
            let mut lps = Vec::new();
            for i in 0..N {
                let mut rng = StdRng::seed_from_u64(0xC0FFEE + i as u64);
                let cell = gen_cell(Phenotype::Glycolytic, &mut rng);
                let (dead, lp, _, _) = sim_cell(&cell, tx, params, &mut rng);
                if dead {
                    lps.push(lp);
                }
            }
            let r = immune_cascade(&lps, N, immune, blockade);
            (lps.len(), r.immune_kills)
        }

        // NEGATIVE: death is impossible AND `baseline_antigenicity` is 0,
        // so no antigen can exist. The second half is now a real condition
        // rather than an automatic one -- see
        // `baseline_antigenicity_makes_blockade_a_treatment` below, which is
        // the whole point of the immunotherapy arm and the case where this
        // assertion MUST NOT hold.
        assert_eq!(
            ImmuneParams::default().baseline_antigenicity,
            0.0,
            "the default gained antigen presentation; every committed immune \
             number moved and this test is now about a different engine"
        );
        let no_death = Params {
            death_threshold: f64::INFINITY,
            ..Params::default()
        };
        for &efficacy in &[0.0, 0.5, 1.0] {
            let immune = ImmuneParams {
                anti_pd1_efficacy: efficacy,
                ..ImmuneParams::default()
            };
            for &blockade in &[false, true] {
                for &tx in &[Treatment::Control, Treatment::RSL3, Treatment::SDT] {
                    let (dead, kills) = run(&no_death, tx, &immune, blockade);
                    assert_eq!(dead, 0, "{tx:?} killed cells with an infinite threshold");
                    assert_eq!(
                        kills, 0.0,
                        "{tx:?} produced {kills} immune kills with no ferroptotic \
                         death (anti_pd1={efficacy}, blockade={blockade})"
                    );
                }
            }
        }

        // POSITIVE: with death possible the same chain DOES kill, and blockade
        // moves it. Without this the negative case is satisfied by a chain
        // that is simply disconnected -- and `immune_cascade` caps kills at
        // the survivors, so the rates are scaled down here to keep the
        // comparison off that ceiling. At the defaults both arms saturate at
        // the same number and the blockade term is invisible.
        let immune = ImmuneParams {
            tcell_kill_rate: 0.02,
            ..ImmuneParams::default()
        };
        let (dead, kills_off) = run(&Params::default(), Treatment::SDT, &immune, false);
        let (_, kills_on) = run(&Params::default(), Treatment::SDT, &immune, true);
        let survivors = (N - dead) as f64;
        assert!(
            dead > 0,
            "SDT killed nothing, so the negative case proves nothing"
        );
        assert!(
            kills_off > 0.0,
            "the immune chain is disconnected: {kills_off} kills"
        );
        assert!(
            kills_on < survivors,
            "both arms are at the survivor cap ({survivors}), so this cannot \
             see the blockade term at all"
        );
        assert!(
            kills_on > kills_off,
            "blockade did not raise kills ({kills_on} vs {kills_off}), so the \
             sweep in the negative case is not exercising a live path"
        );
    }

    /// The arm itself: with antigen present, blockade KILLS without ferroptosis.
    ///
    /// This is the assertion `analysis/modality-coverage.md` said the engine
    /// could not make. It filed immunotherapy as a MODIFIER -- the largest
    /// trial count in the taxonomy, and unaskable -- because every activation
    /// path was gated on ferroptotic death, so anti-PD-1 multiplied zero at
    /// every setting.
    ///
    /// Pinned in BOTH directions, because "blockade now does something" is
    /// satisfied by an engine that ignores the brake entirely: kills must be
    /// zero at zero antigenicity, positive above it, MONOTONE in it, and
    /// strictly higher with blockade than without.
    #[test]
    fn baseline_antigenicity_makes_blockade_a_treatment() {
        const N: usize = 10_000;
        let off = ImmuneParams::default();
        assert_eq!(off.baseline_antigenicity, 0.0);
        assert_eq!(
            immune_cascade(&[], N, &off, true).immune_kills,
            0.0,
            "with no antigen and no death, blockade must still kill nothing"
        );

        // Positive, monotone, and responsive to the brake.
        let mut prev = 0.0;
        for &a in &[0.001, 0.005, 0.01, 0.05] {
            let p = ImmuneParams {
                baseline_antigenicity: a,
                // Keep the kill rate off the survivor cap so the brake term
                // stays visible; at the defaults both arms saturate.
                tcell_kill_rate: 0.02,
                ..ImmuneParams::default()
            };
            let no_drug = immune_cascade(&[], N, &p, false).immune_kills;
            let drug = immune_cascade(&[], N, &p, true).immune_kills;
            assert!(
                no_drug > 0.0,
                "antigenicity {a} with NO ferroptotic death produced no kills, \
                 so blockade is still not a treatment"
            );
            assert!(
                drug > no_drug,
                "blockade did not raise kills at antigenicity {a}: {drug} vs \
                 {no_drug}"
            );
            assert!(drug < N as f64, "kills hit the population cap at {a}");
            assert!(
                drug > prev,
                "kills are not monotone in antigenicity at {a}: {drug} vs {prev}"
            );
            prev = drug;
        }

        // And it composes with ferroptotic death rather than replacing it:
        // the same antigenicity plus real deaths must kill MORE.
        let p = ImmuneParams {
            baseline_antigenicity: 0.005,
            tcell_kill_rate: 0.02,
            ..ImmuneParams::default()
        };
        let alone = immune_cascade(&[], N, &p, true).immune_kills;
        let with_death = immune_cascade(&[8.0; 50], N, &p, true).immune_kills;
        assert!(
            with_death > alone,
            "ferroptotic death no longer adds to the baseline term: \
             {with_death} vs {alone}"
        );
    }

    /// Redirected effectors do NOT wait for antigen presentation, and the
    /// brake still applies to them.
    ///
    /// Both halves are the claim. The first is the modality: CAR-T and
    /// bispecifics bypass DC priming, so a run with no ferroptotic death and
    /// no baseline antigenicity must still kill -- exactly the case
    /// `immune_cascade` returns zero for. The second is what keeps the model
    /// honest about solid tumours: CAR-T cells express PD-1 and are
    /// suppressed, which is the documented reason the 70-94% B-ALL band has
    /// not transferred, and a model exempting them would predict efficacy the
    /// literature does not report.
    #[test]
    fn adoptive_transfer_bypasses_priming_but_not_the_brake() {
        const N: usize = 10_000;
        let p = ImmuneParams::default();

        // The case the DC cascade cannot reach: no deaths, no antigen.
        assert_eq!(immune_cascade(&[], N, &p, true).immune_kills, 0.0);
        let redirected = adoptive_transfer_kills(100.0, N, &p, 0.0, false);
        assert!(
            redirected > 0.0,
            "redirected effectors produced no kills without prior death, so \
             the modality is indistinguishable from the DC cascade"
        );

        // Zero effectors is bit-zero, so an unconfigured run is unmoved.
        assert_eq!(
            adoptive_transfer_kills(0.0, N, &p, 0.0, false).to_bits(),
            0.0_f64.to_bits()
        );

        // The brake applies, and blockade lifts it.
        let blocked = adoptive_transfer_kills(100.0, N, &p, 0.0, true);
        assert!(
            blocked > redirected,
            "anti-PD-1 did not raise redirected kills ({blocked} vs \
             {redirected}); CAR-T cells are not brake-exempt"
        );

        // Suppression is the leukaemia-vs-solid-tumour difference, and it
        // must be able to erase the effect entirely rather than shade it.
        let suppressed = adoptive_transfer_kills(100.0, N, &p, 0.9, false);
        assert!(
            suppressed < redirected * 0.2,
            "{suppressed} vs {redirected}"
        );
        assert_eq!(
            adoptive_transfer_kills(100.0, N, &p, 1.0, true).to_bits(),
            0.0_f64.to_bits(),
            "total suppression must zero the kill, or the solid-tumour \
             setting is unreachable"
        );

        // Monotone in effector count, and capped at the population.
        let mut prev = 0.0;
        for &e in &[1.0_f64, 10.0, 100.0, 1_000.0] {
            let k = adoptive_transfer_kills(e, N, &p, 0.0, false);
            assert!(k > prev, "not monotone at {e}");
            prev = k;
        }
        assert!(adoptive_transfer_kills(1e9, N, &p, 0.0, true) <= N as f64);
    }

    /// A fully CD47-protected tumour is not eaten at all, and blockade is what
    /// changes that.
    #[test]
    fn phagocytosis_is_gated_by_the_do_not_eat_me_signal() {
        const N: usize = 10_000;
        // Full protection, no drug: bit-zero. A small nonzero here would make
        // the checkpoint look leaky and the drug look optional.
        assert_eq!(
            phagocytosis_kills(500.0, N, 0.01, 1.0, 0.0).to_bits(),
            0.0_f64.to_bits()
        );
        // No protection at all: the macrophages engulf freely.
        let unprotected = phagocytosis_kills(500.0, N, 0.01, 0.0, 0.0);
        assert!(unprotected > 0.0);
        // Blockade recovers the unprotected rate exactly at full efficacy --
        // the drug removes the brake, it does not add an effect of its own.
        let fully_blocked = phagocytosis_kills(500.0, N, 0.01, 1.0, 1.0);
        assert!(
            (fully_blocked - unprotected).abs() < 1e-12,
            "{fully_blocked}"
        );
        // Monotone in blockade, and bounded.
        let mut prev = 0.0;
        for &eff in &[0.0_f64, 0.25, 0.5, 0.75, 1.0] {
            let k = phagocytosis_kills(500.0, N, 0.01, 1.0, eff);
            assert!(k >= prev, "not monotone at efficacy {eff}");
            prev = k;
        }
        assert!(phagocytosis_kills(1e9, N, 1.0, 0.0, 1.0) <= N as f64);
    }

    /// The published CAR-T band is quoted WITH the caveat from the same
    /// sentence, so the constant cannot be read as a target this engine meets.
    /// CAR-T persists and expands; a bispecific washes out. The DIRECTIONS
    /// are the model, and they are opposite.
    #[test]
    fn the_two_effector_sources_differ_in_persistence_not_in_killing() {
        let (initial, cap) = (100.0_f64, 10_000.0_f64);

        let cart_early = EffectorSource::CarT.effectors_after(initial, 5, 0.3, cap);
        let cart_late = EffectorSource::CarT.effectors_after(initial, 40, 0.3, cap);
        assert!(cart_early > initial, "CAR-T did not expand: {cart_early}");
        assert!(cart_late > cart_early, "CAR-T expansion is not sustained");
        assert!(cart_late <= cap, "expansion escaped the cap: {cart_late}");
        // Bounded, not unbounded: a model without the cap predicts cures the
        // literature does not report.
        let runaway = EffectorSource::CarT.effectors_after(initial, 500, 0.9, cap);
        assert!(runaway <= cap && (runaway - cap).abs() < cap * 0.01);

        let bs_early = EffectorSource::BispecificEngager.effectors_after(initial, 5, 0.3, cap);
        let bs_late = EffectorSource::BispecificEngager.effectors_after(initial, 40, 0.3, cap);
        assert!(bs_early < initial, "the engager did not clear: {bs_early}");
        assert!(bs_late < bs_early, "clearance is not sustained");

        // The two must move in OPPOSITE directions over the same interval,
        // which is the whole distinction; a shared sign would make the enum
        // decorative.
        assert!(
            cart_late > initial && bs_late < initial,
            "CAR-T {cart_late} and engager {bs_late} did not diverge from \
             {initial}"
        );

        // Zero rate is the identity for BOTH, so an unconfigured run is
        // unmoved whichever source it names.
        for src in [EffectorSource::CarT, EffectorSource::BispecificEngager] {
            assert_eq!(
                src.effectors_after(initial, 20, 0.0, cap).to_bits(),
                initial.to_bits(),
                "{src:?} moved at rate 0"
            );
            assert_eq!(src.effectors_after(0.0, 20, 0.5, cap), 0.0);
        }

        // And the KILL term is genuinely shared: same effector count, same
        // kills, whichever source produced it.
        let p = ImmuneParams::default();
        let k = adoptive_transfer_kills(250.0, 10_000, &p, 0.1, false);
        assert!(k > 0.0);
        assert_eq!(
            adoptive_transfer_kills(250.0, 10_000, &p, 0.1, false).to_bits(),
            k.to_bits()
        );
    }

    /// Oncolytic lysis feeds the SAME ICD chain ferroptotic death does, and
    /// that sharing is the modelling claim.
    #[test]
    fn oncolytic_lysis_drives_the_shared_icd_chain() {
        const N: usize = 10_000;
        let p = ImmuneParams::default();

        // Unconfigured is bit-zero on both outputs.
        let (l0, d0) = oncolytic_lysis(N, 0.0, 0.9, 0.8);
        assert_eq!(l0.to_bits(), 0.0_f64.to_bits());
        assert_eq!(d0.to_bits(), 0.0_f64.to_bits());

        // Lysis is monotone in infection and bounded by the population.
        let mut prev = 0.0;
        for &f in &[0.01_f64, 0.1, 0.5, 1.0] {
            let (l, _) = oncolytic_lysis(N, f, 0.9, 0.8);
            assert!(l > prev, "lysis not monotone at infected fraction {f}");
            assert!(l <= N as f64);
            prev = l;
        }

        // THE POINT: lysed cells drive the immune cascade exactly as
        // ferroptotic deaths do, so an oncolytic run produces immune kills
        // through the shared chain rather than through a parallel one.
        let (lysed, quality) = oncolytic_lysis(N, 0.2, 0.9, 0.8);
        let lps: Vec<f64> = vec![quality * 10.0; lysed as usize];
        let immune = immune_cascade(&lps, N, &p, false);
        assert!(
            immune.immune_kills > 0.0,
            "viral lysis produced no immune kills, so the modality's systemic \
             effect -- responses at UNINJECTED lesions -- is unreachable"
        );
        // And blockade raises it, since the chain is the same one.
        let with_drug = immune_cascade(&lps, N, &p, true);
        assert!(with_drug.immune_kills >= immune.immune_kills);

        // Immunogenicity must MOVE the quality signal, or the second effect
        // is decorative.
        let (_, low) = oncolytic_lysis(N, 0.2, 0.9, 0.1);
        let (_, high) = oncolytic_lysis(N, 0.2, 0.9, 0.9);
        assert!(
            high > low,
            "immunogenicity does not change the death quality"
        );
    }

    /// The T-VEC constant carries its CONTROL arm, because 16% alone reads as
    /// a modest result and 16% against 2% is an eightfold ratio.
    #[test]
    fn the_oncolytic_band_keeps_its_control_arm() {
        let (treated, control) = TVEC_DURABLE_RESPONSE;
        assert!(treated > control, "{treated} vs {control}");
        assert!(
            treated / control > 5.0,
            "the treated/control ratio is {:.1}, which is not the published \
             effect",
            treated / control
        );
        assert!(
            treated < 0.5,
            "a durable response rate of {treated} is not T-VEC's"
        );
    }

    #[test]
    fn the_cart_band_is_the_leukaemia_one_and_the_docs_say_so() {
        let (lo, hi) = CART_B_ALL_CR_BAND;
        assert!(lo > 0.0 && hi > lo && hi <= 1.0, "{lo}-{hi}");
        let src = include_str!("immune.rs");
        for caveat in ["solid tumours has been less successful", "not durable"] {
            assert!(
                src.contains(caveat),
                "the CAR-T band no longer carries its own caveat: {caveat:?}"
            );
        }
    }

    #[test]
    fn activation_tracks_per_cell_quality_not_count() {
        // The model's central claim: DC activation responds to DAMP-PER-DEATH
        // (kill quality), so a few high-LP deaths out-activate many low-LP
        // deaths. This is the SDT-vs-RSL3 asymmetry the 2D cascade encodes.
        let p = ImmuneParams::default();
        let few_high = immune_cascade(&[10.0], 10_000, &p, false); // 1 cell, high LP
        let many_low = immune_cascade(&vec![1.0; 50], 10_000, &p, false); // 50 cells, low LP
        assert!(
            few_high.dc_activation_fraction > many_low.dc_activation_fraction,
            "quality should beat quantity: {} !> {}",
            few_high.dc_activation_fraction,
            many_low.dc_activation_fraction
        );
    }

    #[test]
    fn anti_pd1_increases_killing() {
        // Removing part of the PD-1 brake raises kill efficiency; with ample
        // remaining tumor (uncapped), anti-PD-1 must yield more immune kills.
        let p = ImmuneParams::default();
        let lps = vec![5.0; 20];
        let base = immune_cascade(&lps, 1_000_000, &p, false);
        let treated = immune_cascade(&lps, 1_000_000, &p, true);
        assert!(
            treated.immune_kills > base.immune_kills,
            "anti-PD-1 {} should exceed baseline {}",
            treated.immune_kills,
            base.immune_kills
        );
    }

    #[test]
    fn immune_kills_capped_at_remaining_tumor() {
        let p = ImmuneParams::default();
        // total_tumor == n_dead ⇒ zero alive remaining ⇒ no immune kills possible.
        let r = immune_cascade(&[8.0; 30], 30, &p, true);
        assert_eq!(r.immune_kills, 0.0, "no alive cells left to kill");
        // And in general immune_kills never exceeds remaining.
        let r2 = immune_cascade(&[8.0; 30], 35, &p, true);
        assert!(
            r2.immune_kills <= 5.0 + 1e-9,
            "kills {} > remaining 5",
            r2.immune_kills
        );
    }
}
