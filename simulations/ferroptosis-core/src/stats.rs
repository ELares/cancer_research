//! Statistics, results aggregation, and parallel execution.

use rand::prelude::*;
use rayon::prelude::*;
use serde::Serialize;

use crate::biochem::sim_cell;
use crate::cell::{gen_cell, Phenotype, Treatment};
use crate::params::Params;

/// Aggregated result for one condition (phenotype × treatment).
#[derive(Serialize, Clone, Debug)]
pub struct SimResult {
    pub phenotype: String,
    pub treatment: String,
    pub n_cells: usize,
    pub n_dead: usize,
    pub death_rate: f64,
    pub ci_low: f64,
    pub ci_high: f64,
    pub mean_lipid_perox: f64,
    pub mean_gsh_final: f64,
    pub mean_gpx4_final: f64,
}

/// Wilson score confidence interval for binomial proportion.
///
/// Returns (0.0, 0.0) if `n == 0`. Panics in debug mode if `k > n`.
#[must_use]
pub fn wilson_ci(n: usize, k: usize) -> (f64, f64) {
    debug_assert!(k <= n, "wilson_ci: k ({k}) > n ({n})");
    if n == 0 {
        return (0.0, 0.0);
    }
    let (nf, p, z) = (n as f64, k as f64 / n as f64, 1.96);
    let d = 1.0 + z * z / nf;
    let c = (p + z * z / (2.0 * nf)) / d;
    let s = z * ((p * (1.0 - p) / nf + z * z / (4.0 * nf * nf)).sqrt()) / d;
    ((c - s).max(0.0), (c + s).min(1.0))
}

/// Run a single condition (phenotype × treatment) with n cells in parallel.
/// Uses the same RNG seeding as the original v3 for bitwise reproducibility.
pub fn run_condition(
    pheno: Phenotype,
    tx: Treatment,
    params: &Params,
    n: usize,
    pname: &str,
    tname: &str,
) -> SimResult {
    let outcomes: Vec<(bool, f64, f64, f64)> = (0..n)
        .into_par_iter()
        .map(|i| {
            let mut cell_rng = StdRng::seed_from_u64(i as u64 * 2);
            let mut sim_rng = StdRng::seed_from_u64(i as u64 * 2 + 1);
            let cell = gen_cell(pheno, &mut cell_rng);
            sim_cell(&cell, tx, params, &mut sim_rng)
        })
        .collect();

    let dead = outcomes.iter().filter(|(d, _, _, _)| *d).count();
    let dr = dead as f64 / n as f64;
    let (cl, ch) = wilson_ci(n, dead);

    SimResult {
        phenotype: pname.to_string(),
        treatment: tname.to_string(),
        n_cells: n,
        n_dead: dead,
        death_rate: dr,
        ci_low: cl,
        ci_high: ch,
        mean_lipid_perox: outcomes.iter().map(|(_, l, _, _)| l).sum::<f64>() / n as f64,
        mean_gsh_final: outcomes.iter().map(|(_, _, g, _)| g).sum::<f64>() / n as f64,
        mean_gpx4_final: outcomes.iter().map(|(_, _, _, p)| p).sum::<f64>() / n as f64,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // The Wilson-score CI feeds ci_low/ci_high in sim-original, sim-window,
    // sim-invivo, sim-tissue-pk, sim-combo-mech, sim-tumor-pk, and the Python
    // bindings, but had no coverage (#295). A wrong z-constant or formula slip
    // would silently skew every reported interval; these pin the invariants.

    #[test]
    fn n_zero_returns_zero_interval() {
        assert_eq!(wilson_ci(0, 0), (0.0, 0.0));
    }

    #[test]
    fn interval_is_clamped_to_unit() {
        // All-dead and none-dead must not escape [0, 1]. The `.max(0.0)` clamp
        // genuinely fires at k=0 (the unclamped lower bound is tiny negative
        // fp-noise); the `.min(1.0)` clamp is defensive (the Wilson upper bound
        // is exactly 1.0 at k=n, never overshooting for valid inputs).
        let (lo0, hi0) = wilson_ci(10, 0);
        assert!(lo0 >= 0.0 && hi0 <= 1.0, "k=0: ({lo0}, {hi0})");
        assert_eq!(lo0, 0.0, "k=0 lower bound clamps to 0");
        let (lo1, hi1) = wilson_ci(10, 10);
        assert!(lo1 >= 0.0 && hi1 <= 1.0, "k=n: ({lo1}, {hi1})");
        assert_eq!(hi1, 1.0, "k=n upper bound clamps to 1");
    }

    #[test]
    fn invariants_hold_across_input_grid() {
        // Broad sweep: for every (n, k), the interval must satisfy
        // 0 <= lo <= p <= hi <= 1. A sign error, swapped bound, or bad clamp
        // that the point-cases miss would surface somewhere in the grid.
        for n in [1usize, 2, 5, 17, 100, 999, 5000] {
            for k in 0..=n {
                let p = k as f64 / n as f64;
                let (lo, hi) = wilson_ci(n, k);
                assert!(0.0 <= lo, "n={n} k={k}: lo {lo} < 0");
                assert!(lo <= p + 1e-12, "n={n} k={k}: lo {lo} > p {p}");
                assert!(p <= hi + 1e-12, "n={n} k={k}: p {p} > hi {hi}");
                assert!(hi <= 1.0, "n={n} k={k}: hi {hi} > 1");
                assert!(lo <= hi, "n={n} k={k}: lo {lo} > hi {hi}");
            }
        }
    }

    #[test]
    fn interval_contains_point_estimate() {
        for (n, k) in [(10, 3), (100, 50), (1000, 137), (50, 1)] {
            let p = k as f64 / n as f64;
            let (lo, hi) = wilson_ci(n, k);
            assert!(
                lo <= p && p <= hi,
                "p={p} not in ({lo}, {hi}) for n={n}, k={k}"
            );
        }
    }

    #[test]
    fn symmetric_around_half_at_p_half() {
        // Wilson center is exactly p when the z-correction is symmetric; at
        // p=0.5 the interval midpoint must be 0.5.
        for n in [10, 100, 1000] {
            let (lo, hi) = wilson_ci(n, n / 2);
            let mid = (lo + hi) / 2.0;
            assert!((mid - 0.5).abs() < 1e-9, "n={n}: midpoint {mid} != 0.5");
            assert!(
                lo < 0.5 && 0.5 < hi,
                "n={n}: 0.5 not strictly inside ({lo}, {hi})"
            );
        }
    }

    #[test]
    fn lower_bound_monotonic_in_k() {
        // More successes ⇒ higher lower bound (point estimate rises).
        let n = 200;
        let mut prev = -1.0;
        for k in [0, 20, 50, 100, 150, 200] {
            let (lo, _) = wilson_ci(n, k);
            assert!(lo >= prev, "ci_low not monotonic at k={k}: {lo} < {prev}");
            prev = lo;
        }
    }

    #[test]
    fn interval_narrows_as_n_grows() {
        // Same proportion (0.5), larger n ⇒ tighter interval.
        let (lo_small, hi_small) = wilson_ci(10, 5);
        let (lo_big, hi_big) = wilson_ci(1000, 500);
        assert!(
            (hi_big - lo_big) < (hi_small - lo_small),
            "width(1000) {} not < width(10) {}",
            hi_big - lo_big,
            hi_small - lo_small
        );
    }

    #[test]
    fn known_value_n100_k50() {
        // z=1.96, p=0.5, n=100: center 0.5, half-width via the closed form.
        let (lo, hi) = wilson_ci(100, 50);
        // Reference values from the Wilson formula with z=1.96.
        assert!((lo - 0.4038).abs() < 1e-3, "lo={lo}");
        assert!((hi - 0.5962).abs() < 1e-3, "hi={hi}");
    }

    /// Audit follow-up: `wilson_ci`'s `k <= n` `debug_assert` documents a
    /// debug-mode panic contract, but no test exercised it. Feeding k > n must
    /// panic in debug builds. Gated to debug builds, where the assert lives (it
    /// is compiled out under `--release`, so it would not panic there).
    #[test]
    #[cfg(debug_assertions)]
    #[should_panic(expected = "wilson_ci")]
    fn wilson_ci_panics_when_k_exceeds_n() {
        let _ = wilson_ci(5, 9);
    }
}
