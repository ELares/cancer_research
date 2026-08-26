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

    // Mature DCs: activation quality × maturation rate × number of antigen-presenting deaths
    let mature_dcs = dc_activation_fraction * params.dc_maturation_rate * n_dead as f64;

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

        // NEGATIVE: death is impossible, so no antigen can exist.
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
