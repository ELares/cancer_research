//! Cell types, phenotypes, treatments, and cell generation.

use rand::prelude::*;
use rand_distr::Normal;
use serde::{Deserialize, Serialize};

/// Sample from a normal distribution. Does NOT clamp — callers must apply
/// `.max(threshold)` if a positive value is required.
pub fn norm(rng: &mut StdRng, mean: f64, sd: f64) -> f64 {
    Normal::new(mean, sd).unwrap().sample(rng)
}

/// A single cell's biochemical state.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Cell {
    /// Labile iron pool (µM). Ref: 0.2-1.5 normal, 2-6 overloaded
    /// (Kakhlon & Cabantchik, Free Radic Biol Med 2002)
    pub iron: f64,
    /// Glutathione (mM). Ref: 1-10 mM intracellular
    /// (Forman et al., Free Radic Biol Med 2009)
    pub gsh: f64,
    /// GPX4 activity (relative, 1.0 = normal). Ref: kcat ~40/s
    /// (Ursini et al., Free Radic Biol Med 1995)
    pub gpx4: f64,
    /// FSP1/DHODH activity (relative). GPX4-independent CoQ10 pathway.
    /// (Bersuker et al., Nature 2019; Mao et al., Nature 2021)
    pub fsp1: f64,
    /// Basal mitochondrial ROS production (relative).
    /// OXPHOS cells: ~2-3× higher due to active ETC (Murphy, Biochem J 2009)
    pub basal_ros: f64,
    /// Lipid unsaturation: PUFA content determines peroxidation susceptibility.
    /// OXPHOS cells have more mitochondrial membranes = more target.
    /// (Yang et al., Cell 2016 — PUFA requirement for ferroptosis)
    pub lipid_unsat: f64,
    /// NRF2 transcriptional activity. Master regulator of antioxidant response.
    /// Drives GSH synthesis (via GCL/GSS), GPX4 expression.
    /// (Dodson et al., Free Radic Biol Med 2019)
    pub nrf2: f64,
    /// Per-cell MUFA carrying capacity (#270). `None` ⇒ the global
    /// `Params::scd_mufa_max`. The spheroid model sets this radially (rim-high,
    /// core-low) so position-dependent MUFA is **durable**:
    /// `update_mufa_protection` relaxes toward a steady state that scales with
    /// this cap, instead of every cell converging to the single uniform M_ss.
    /// `None` (the default) keeps every non-spheroid path byte-identical.
    ///
    /// When `Some`, this is the effective cap — it **replaces** (and may exceed)
    /// the global `Params::scd_mufa_max`, since `update_mufa_protection` clamps
    /// to the per-cell value. The old global 0.25 ceiling on `SpheroidConfig`'s
    /// MUFA values therefore no longer applies; a config can now set a
    /// supra-global per-cell cap deliberately.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub mufa_cap: Option<f64>,
    /// Per-cell SCD1/MUFA accumulation RATE (#363). `None` ⇒ the global
    /// `Params::scd_mufa_rate`. SCD1/MUFA enrichment kinetics differ by
    /// phenotype (e.g. a drug-tolerant persister remodels lipids toward MUFA at a
    /// different rate than a proliferating glycolytic cell), so a consumer can set
    /// this per-phenotype via [`crate::phenotype_mufa::PhenotypeMufaConfig`] to give the
    /// acute-versus-established MUFA build-up (`mufa_acute_start`, #339) a
    /// phenotype-specific time constant instead of one shared `scd_mufa_rate`.
    ///
    /// `None` (the default) keeps every path byte-identical (the call sites fall
    /// back to `params.scd_mufa_rate`, exactly as before). Independent of
    /// [`Cell::mufa_cap`]: the rate controls how fast MUFA protection accumulates,
    /// the cap controls the steady state it saturates toward.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub mufa_rate: Option<f64>,
}

/// Treatment modalities.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum Treatment {
    Control,
    RSL3,
    SDT,
    PDT,
}

/// Cell phenotypes.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum Phenotype {
    Glycolytic,
    OXPHOS,
    Persister,
    PersisterNrf2,
    /// Non-tumor stromal cells (for spatial model). Low iron, high GSH, inert.
    Stromal,
}

/// Generate a cell with stochastic parameters for the given phenotype.
/// Identical to the original v3 simulation for the four base phenotypes.
pub fn gen_cell(pheno: Phenotype, rng: &mut StdRng) -> Cell {
    match pheno {
        Phenotype::Glycolytic => Cell {
            iron: norm(rng, 1.0, 0.25).max(0.3),
            gsh: norm(rng, 5.0, 1.0).max(1.5),
            gpx4: norm(rng, 1.0, 0.12).max(0.4),
            fsp1: norm(rng, 1.0, 0.12).max(0.4),
            basal_ros: norm(rng, 0.2, 0.05).max(0.05),
            lipid_unsat: norm(rng, 1.0, 0.12).max(0.5),
            nrf2: norm(rng, 1.0, 0.12).max(0.4),
            mufa_cap: None,
            mufa_rate: None,
        },
        Phenotype::OXPHOS => Cell {
            iron: norm(rng, 2.8, 0.6).max(0.8),
            gsh: norm(rng, 4.0, 0.8).max(1.0),
            gpx4: norm(rng, 1.0, 0.12).max(0.4),
            fsp1: norm(rng, 1.0, 0.12).max(0.4),
            basal_ros: norm(rng, 0.5, 0.12).max(0.1),
            lipid_unsat: norm(rng, 1.6, 0.2).max(0.7),
            nrf2: norm(rng, 1.2, 0.15).max(0.5),
            mufa_cap: None,
            mufa_rate: None,
        },
        Phenotype::Persister => Cell {
            iron: norm(rng, 1.5, 0.3).max(0.5),
            gsh: norm(rng, 4.8, 0.8).max(1.8),
            gpx4: norm(rng, 0.7, 0.15).max(0.15),
            fsp1: norm(rng, 0.15, 0.06).max(0.01),
            basal_ros: norm(rng, 0.25, 0.06).max(0.05),
            lipid_unsat: norm(rng, 1.4, 0.15).max(0.6),
            nrf2: norm(rng, 0.7, 0.15).max(0.2),
            mufa_cap: None,
            mufa_rate: None,
        },
        Phenotype::PersisterNrf2 => Cell {
            iron: norm(rng, 2.8, 0.6).max(0.8),
            gsh: norm(rng, 7.0, 1.2).max(3.0),
            gpx4: norm(rng, 1.3, 0.15).max(0.5),
            fsp1: norm(rng, 0.2, 0.08).max(0.02),
            basal_ros: norm(rng, 0.5, 0.12).max(0.1),
            lipid_unsat: norm(rng, 1.6, 0.2).max(0.7),
            nrf2: norm(rng, 3.0, 0.4).max(1.5),
            mufa_cap: None,
            mufa_rate: None,
        },
        Phenotype::Stromal => Cell {
            iron: norm(rng, 0.3, 0.08).max(0.1),
            gsh: norm(rng, 8.0, 1.0).max(4.0),
            gpx4: norm(rng, 1.5, 0.15).max(0.8),
            fsp1: norm(rng, 1.0, 0.12).max(0.4),
            basal_ros: norm(rng, 0.1, 0.03).max(0.02),
            lipid_unsat: norm(rng, 0.6, 0.1).max(0.3),
            nrf2: norm(rng, 1.0, 0.12).max(0.4),
            mufa_cap: None,
            mufa_rate: None,
        },
    }
}

/// Recovery rates for persister phenotype transitions post-chemotherapy.
/// Used by the vulnerability window simulation.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct RecoveryRates {
    /// FSP1 re-expression half-time (days). Slowest — epigenetic.
    pub fsp1_half_recovery_days: f64,
    /// GPX4 re-expression half-time (days). Transcriptional.
    pub gpx4_half_recovery_days: f64,
    /// NRF2 re-activation half-time (days).
    pub nrf2_half_recovery_days: f64,
    /// GSH resynthesis half-time (days). Fastest — metabolic.
    pub gsh_half_recovery_days: f64,
}

impl Default for RecoveryRates {
    fn default() -> Self {
        RecoveryRates {
            fsp1_half_recovery_days: 7.0,
            gpx4_half_recovery_days: 3.0,
            nrf2_half_recovery_days: 5.0,
            gsh_half_recovery_days: 1.0,
        }
    }
}

/// Apply time-dependent recovery to a persister cell.
/// Returns new parameter means (not a full Cell — caller generates stochastic cell from these).
pub fn recovered_persister_means(days: f64, rates: &RecoveryRates) -> (f64, f64, f64, f64) {
    // Exponential recovery: fraction recovered = 1 - exp(-ln(2) * t / t_half)
    let frac = |t_half: f64| -> f64 { 1.0 - (-(2.0_f64.ln()) * days / t_half).exp() };

    // Persister baseline → Glycolytic normal targets
    let fsp1 = 0.15 + (1.0 - 0.15) * frac(rates.fsp1_half_recovery_days);
    let gpx4 = 0.7 + (1.0 - 0.7) * frac(rates.gpx4_half_recovery_days);
    let nrf2 = 0.7 + (1.0 - 0.7) * frac(rates.nrf2_half_recovery_days);
    let gsh = 4.8 + (5.0 - 4.8) * frac(rates.gsh_half_recovery_days);

    (fsp1, gpx4, nrf2, gsh)
}

/// Generate a persister cell at a given number of days post-chemo withdrawal.
pub fn gen_recovered_persister(days: f64, rates: &RecoveryRates, rng: &mut StdRng) -> Cell {
    let (fsp1_mean, gpx4_mean, nrf2_mean, gsh_mean) = recovered_persister_means(days, rates);
    Cell {
        iron: norm(rng, 1.5, 0.3).max(0.5),
        gsh: norm(rng, gsh_mean, 0.8).max(1.8),
        gpx4: norm(rng, gpx4_mean, 0.15).max(0.15),
        fsp1: norm(rng, fsp1_mean, 0.06).max(0.01),
        basal_ros: norm(rng, 0.25, 0.06).max(0.05),
        lipid_unsat: norm(rng, 1.4, 0.15).max(0.6),
        nrf2: norm(rng, nrf2_mean, 0.15).max(0.2),
        mufa_cap: None,
        mufa_rate: None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // The post-chemo persister-recovery means/cell generator drives sim-combo
    // and sim-window's vulnerability-window result but had no coverage (#295).
    // This is a separate formula from persister.rs's competing-rate step. Pin
    // the exponential-recovery endpoints, the half-time midpoint, monotonicity,
    // and the stochastic generator's determinism + floors.

    #[test]
    fn means_at_day_zero_are_persister_baseline() {
        let (fsp1, gpx4, nrf2, gsh) = recovered_persister_means(0.0, &RecoveryRates::default());
        assert!((fsp1 - 0.15).abs() < 1e-12, "fsp1={fsp1}");
        assert!((gpx4 - 0.7).abs() < 1e-12, "gpx4={gpx4}");
        assert!((nrf2 - 0.7).abs() < 1e-12, "nrf2={nrf2}");
        assert!((gsh - 4.8).abs() < 1e-12, "gsh={gsh}");
    }

    #[test]
    fn means_converge_to_normal_targets() {
        let (fsp1, gpx4, nrf2, gsh) = recovered_persister_means(1.0e6, &RecoveryRates::default());
        assert!((fsp1 - 1.0).abs() < 1e-6, "fsp1={fsp1}");
        assert!((gpx4 - 1.0).abs() < 1e-6, "gpx4={gpx4}");
        assert!((nrf2 - 1.0).abs() < 1e-6, "nrf2={nrf2}");
        assert!((gsh - 5.0).abs() < 1e-6, "gsh={gsh}");
    }

    #[test]
    fn gpx4_is_halfway_recovered_at_its_half_time() {
        // At days == t_half the recovered fraction is exactly 0.5, so
        // gpx4 = 0.7 + (1.0 - 0.7)·0.5 = 0.85, independent of the other rates.
        let rates = RecoveryRates::default();
        let (_, gpx4, _, _) = recovered_persister_means(rates.gpx4_half_recovery_days, &rates);
        assert!(
            (gpx4 - 0.85).abs() < 1e-9,
            "gpx4 at t_half = {gpx4}, expected 0.85"
        );
    }

    #[test]
    fn means_are_monotonic_in_time() {
        let r = RecoveryRates::default();
        let (f0, g0, n0, s0) = recovered_persister_means(0.0, &r);
        let (f1, g1, n1, s1) = recovered_persister_means(1.0, &r);
        let (f10, g10, n10, s10) = recovered_persister_means(10.0, &r);
        assert!(f0 < f1 && f1 < f10, "fsp1 not increasing: {f0},{f1},{f10}");
        assert!(g0 < g1 && g1 < g10, "gpx4 not increasing: {g0},{g1},{g10}");
        assert!(n0 < n1 && n1 < n10, "nrf2 not increasing: {n0},{n1},{n10}");
        assert!(s0 < s1 && s1 < s10, "gsh not increasing: {s0},{s1},{s10}");
    }

    #[test]
    fn gen_recovered_persister_is_deterministic() {
        let r = RecoveryRates::default();
        let a = gen_recovered_persister(2.0, &r, &mut StdRng::seed_from_u64(7));
        let b = gen_recovered_persister(2.0, &r, &mut StdRng::seed_from_u64(7));
        // Same seed ⇒ identical draws on every field.
        assert_eq!(a.iron, b.iron);
        assert_eq!(a.gsh, b.gsh);
        assert_eq!(a.gpx4, b.gpx4);
        assert_eq!(a.fsp1, b.fsp1);
        assert_eq!(a.nrf2, b.nrf2);
        assert_eq!(a.basal_ros, b.basal_ros);
        assert_eq!(a.lipid_unsat, b.lipid_unsat);
    }

    #[test]
    fn gen_recovered_persister_respects_floors_and_caps() {
        // Sweep many seeds at the days where floors are most reachable (day 0
        // and day 1): the `.max(..)` floors must never be violated, and a
        // recovered persister carries no per-cell MUFA cap. NOTE: only the
        // basal_ros (0.05) and fsp1 (0.01) floors actually bind for some seed
        // here — the iron/gsh/gpx4/lipid_unsat/nrf2 means sit several sigma
        // above their floors at every recovery time, so those clamps are
        // defensive. The stronger guard that the generator USES the recovery
        // means correctly is `gen_recovered_persister_means_track_formula`.
        let r = RecoveryRates::default();
        for &days in &[0.0, 1.0] {
            for seed in 0..300u64 {
                let c = gen_recovered_persister(days, &r, &mut StdRng::seed_from_u64(seed));
                assert!(c.iron >= 0.5, "iron {} < 0.5 (seed {seed})", c.iron);
                assert!(c.gsh >= 1.8, "gsh {} < 1.8 (seed {seed})", c.gsh);
                assert!(c.gpx4 >= 0.15, "gpx4 {} < 0.15 (seed {seed})", c.gpx4);
                assert!(c.fsp1 >= 0.01, "fsp1 {} < 0.01 (seed {seed})", c.fsp1);
                assert!(
                    c.basal_ros >= 0.05,
                    "basal_ros {} < 0.05 (seed {seed})",
                    c.basal_ros
                );
                assert!(
                    c.lipid_unsat >= 0.6,
                    "lipid_unsat {} < 0.6 (seed {seed})",
                    c.lipid_unsat
                );
                assert!(c.nrf2 >= 0.2, "nrf2 {} < 0.2 (seed {seed})", c.nrf2);
                assert!(
                    c.mufa_cap.is_none(),
                    "recovered persister should carry no MUFA cap"
                );
            }
        }
    }

    #[test]
    fn gen_recovered_persister_means_track_formula() {
        // The load-bearing guard: the generator must actually DRAW from
        // recovered_persister_means(days) — a regression that swapped fields,
        // ignored the means, or used a wrong target would pass the floor sweep
        // but fail here. At day 5 every mean sits far above its floor, so the
        // `.max()` truncation introduces no bias and the sample means converge
        // to the formula. Also checks the days-independent fields (iron 1.5,
        // basal_ros 0.25, lipid_unsat 1.4).
        let r = RecoveryRates::default();
        let days = 5.0;
        let (fsp1_m, gpx4_m, nrf2_m, gsh_m) = recovered_persister_means(days, &r);
        let n = 4000;
        let (mut s_iron, mut s_gsh, mut s_gpx4, mut s_fsp1) = (0.0, 0.0, 0.0, 0.0);
        let (mut s_basal, mut s_lipid, mut s_nrf2) = (0.0, 0.0, 0.0);
        for seed in 0..n as u64 {
            let c = gen_recovered_persister(days, &r, &mut StdRng::seed_from_u64(seed));
            s_iron += c.iron;
            s_gsh += c.gsh;
            s_gpx4 += c.gpx4;
            s_fsp1 += c.fsp1;
            s_basal += c.basal_ros;
            s_lipid += c.lipid_unsat;
            s_nrf2 += c.nrf2;
        }
        let nf = n as f64;
        // Recovery-driven fields track recovered_persister_means(days).
        assert!(
            (s_gpx4 / nf - gpx4_m).abs() < 0.03,
            "gpx4 mean {} vs {gpx4_m}",
            s_gpx4 / nf
        );
        assert!(
            (s_fsp1 / nf - fsp1_m).abs() < 0.03,
            "fsp1 mean {} vs {fsp1_m}",
            s_fsp1 / nf
        );
        assert!(
            (s_nrf2 / nf - nrf2_m).abs() < 0.03,
            "nrf2 mean {} vs {nrf2_m}",
            s_nrf2 / nf
        );
        assert!(
            (s_gsh / nf - gsh_m).abs() < 0.05,
            "gsh mean {} vs {gsh_m}",
            s_gsh / nf
        );
        // Days-independent fields track their fixed sampling means.
        assert!(
            (s_iron / nf - 1.5).abs() < 0.03,
            "iron mean {}",
            s_iron / nf
        );
        assert!(
            (s_basal / nf - 0.25).abs() < 0.03,
            "basal_ros mean {}",
            s_basal / nf
        );
        assert!(
            (s_lipid / nf - 1.4).abs() < 0.03,
            "lipid_unsat mean {}",
            s_lipid / nf
        );
    }
}
