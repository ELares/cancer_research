//! Therapy-induced senescence as a ferroptosis program (#341).
//!
//! `analysis/principle-resistance-tradeoff.md` lists therapy-induced senescence
//! (TIS) as a primary resistance/escape route, but the suite had no senescence
//! program. This adds an off-by-default one that biases a configurable fraction
//! of tumor cells into a senescence state with two independently-toggleable,
//! well-supported axes: iron accumulation and raised antioxidant/GPX4 defenses.
//!
//! ## Direction is axis-dependent, and it CORRECTS the one-sided "resist" framing
//!
//! Whether senescent cancer cells RESIST or are SENSITIVE to ferroptosis depends
//! on WHICH node of the pathway the therapy hits, and the verified literature is
//! explicit about the split:
//!   - **Upstream triggers** (cystine/system-xc deprivation, iron/erastin):
//!     senescent cells RESIST. They accumulate iron (~30x) but compartmentalize
//!     it in ferritin/lysosomes so it is not Fenton-available, and they raise
//!     antioxidant capacity (Masaldan et al., Redox Biol 2018, PMID 28888202;
//!     Feng et al., Aging 2024, PMID 38683121; Loo et al., Nat Commun 2025, PMID
//!     40731111; Machii et al., FEBS Lett 2026, PMID 42003248).
//!   - **Direct GPX4 inhibition** (RSL3 / ML162-like): senescent cells are
//!     SENSITIVE, a senolytic vulnerability. They are "primed" (high labile
//!     Fe2+, high ROS) and GPX4-DEPENDENT, so removing GPX4 kills them
//!     selectively (D'Ambrosio et al. / Gil, Nat Cell Biol 2026, PMID 42032311;
//!     the resistance papers above ALSO report senescent cells stay sensitive to
//!     direct GPX4 inhibitors even while resisting cystine deprivation).
//!
//! So the issue's one-sided "senescent cells resist acute ferroptosis" framing
//! is incomplete: it holds for upstream triggers but REVERSES under the
//! GPX4-inhibition node, and this model's primary inducer (RSL3) IS a GPX4
//! inhibitor. We let the evidence lead by making the program genuinely
//! BIDIRECTIONAL rather than forcing one sign: the iron and defense axes are
//! independent, so a consumer can express the senolytic-vulnerable state (set
//! `iron_mul > 1`, leave the defense multipliers near `1` ⇒ more ferroptosis
//! under RSL3, D'Ambrosio/Gil 2026) OR the resistant state (set `iron_mul ≈ 1`
//! for compartmentalized iron and raise the defense multipliers ⇒ less
//! ferroptosis, Masaldan 2018 / Loo 2025). The unit test demonstrates BOTH
//! directions under the same RSL3 inducer. The NET sign is set by which axis the
//! applied therapy defeats, which is exactly the contested biology; the model
//! reproduces the contest rather than hiding it behind a single default.
//!
//! ## Model
//!
//! For each senescent cell: `cell.iron *= iron_mul`, `cell.gpx4 *= gpx4_mul`,
//! `cell.nrf2 *= nrf2_mul`, `cell.fsp1 *= fsp1_mul`. The NET ferroptosis outcome
//! emerges from the biochem + the applied treatment (the GPX4 crutch is large
//! absent a GPX4 inhibitor, small once RSL3 knocks GPX4 down ~92%, leaving the
//! raised iron to dominate). Non-proliferation is implicit: the 180-step window
//! has no division, so a senescence growth-arrest needs no extra term. The
//! magnitudes are UNCALIBRATED placeholders; the axis structure (iron-up,
//! defense-up, net set by the therapy node) is the result, not the numbers.
//!
//! A fraction of tumor cells is marked via an INDEPENDENT `StdRng(seed)`, so the
//! cell grid's RNG stream is untouched. [`SenescenceConfig::default`] (identity,
//! `fraction == 0`) marks nothing and perturbs nothing, so the matrix stays
//! byte-identical. The SASP (senescence-associated secretory phenotype)
//! immune coupling is deferred (it is the issue's "optional" part).

use crate::cell::Cell;
use crate::grid::TumorGrid3D;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

/// Senescence-program configuration. `default()` is identity (`fraction == 0`,
/// all multipliers `1.0`) so it is byte-identical. The two axes (iron
/// accumulation vs antioxidant/GPX4 defense) are independent so a consumer can
/// express either the senolytic-vulnerable (iron-heavy) or the resistant
/// (defense-heavy, compartmentalized-iron) state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SenescenceConfig {
    /// Fraction of tumor cells driven into the senescence program. `0.0` ⇒
    /// identity ⇒ byte-identical.
    pub fraction: f64,
    /// Labile-iron multiplier. `> 1` ⇒ the primed, high-Fe2+ senolytic-
    /// vulnerable state (Masaldan 2018; D'Ambrosio 2026). `≈ 1` ⇒ iron
    /// compartmentalized (not Fenton-available), the resistant case.
    pub iron_mul: f64,
    /// GPX4 multiplier. `> 1` ⇒ a GPX4 crutch / dependency: it raises baseline
    /// defense but is exactly what direct GPX4 inhibition (RSL3) removes,
    /// producing the senolytic vulnerability.
    pub gpx4_mul: f64,
    /// NRF2 (antioxidant setpoint) multiplier. `> 1` ⇒ raised GSH/GPX4
    /// regeneration capacity.
    pub nrf2_mul: f64,
    /// FSP1 (GPX4-independent backup) multiplier. `> 1` ⇒ raised backup defense.
    pub fsp1_mul: f64,
}

impl Default for SenescenceConfig {
    fn default() -> Self {
        // Identity: no cells senescent, no perturbation ⇒ byte-identical.
        SenescenceConfig {
            fraction: 0.0,
            iron_mul: 1.0,
            gpx4_mul: 1.0,
            nrf2_mul: 1.0,
            fsp1_mul: 1.0,
        }
    }
}

impl SenescenceConfig {
    /// Literature-motivated placeholder encoding the broadly-documented
    /// senescent state: iron accumulation (Masaldan 2018 PMID 28888202) AND
    /// raised antioxidant/GPX4 capacity (Loo 2025 PMID 40731111). With these
    /// magnitudes the antioxidant axis happens to outweigh the raised iron under
    /// RSL3, so the NET is resistant, but that is NOT a fixed property of the
    /// program: it is parameter- and therapy-node-dependent (the unit test drives
    /// both senolytic and resistant nets from explicit iron- vs defense-dominant
    /// configs). UNCALIBRATED magnitudes; calibrate vs senescent-vs-non-senescent
    /// ferroptosis dose-response under GPX4 inhibition (D'Ambrosio/Gil 2026 PMID
    /// 42032311) and under cystine deprivation (Loo 2025 PMID 40731111).
    pub fn literature() -> Self {
        SenescenceConfig {
            fraction: 0.2,
            iron_mul: 2.5,
            gpx4_mul: 1.3,
            nrf2_mul: 1.3,
            fsp1_mul: 1.3,
        }
    }

    /// True when the config applies no perturbation (`fraction == 0` or all
    /// multipliers are `1.0`).
    pub fn is_identity(&self) -> bool {
        self.fraction == 0.0
            || (self.iron_mul == 1.0
                && self.gpx4_mul == 1.0
                && self.nrf2_mul == 1.0
                && self.fsp1_mul == 1.0)
    }
}

/// Apply the senescence program to ONE cell: scale its labile iron and its
/// antioxidant/GPX4 defenses. Pure, deterministic, no RNG. Multipliers are
/// floored at `0`.
pub fn apply_senescence_to_cell(cell: &mut Cell, cfg: &SenescenceConfig) {
    cell.iron *= cfg.iron_mul.max(0.0);
    cell.gpx4 *= cfg.gpx4_mul.max(0.0);
    cell.nrf2 *= cfg.nrf2_mul.max(0.0);
    cell.fsp1 *= cfg.fsp1_mul.max(0.0);
}

/// Mark a `fraction` of tumor cells senescent and apply the program to them.
/// Uses an INDEPENDENT `StdRng(seed)` so the cell grid's RNG stream is
/// untouched. Returns the per-cell senescence mask (`true` = senescent). No-op
/// (empty effect, all-`false` mask) when `cfg.is_identity()`, so the matrix
/// stays byte-identical when the layer is off.
pub fn apply_senescence_program_3d(
    grid: &mut TumorGrid3D,
    cfg: &SenescenceConfig,
    seed: u64,
) -> Vec<bool> {
    let n = grid.cells.len();
    let mut mask = vec![false; n];
    if cfg.is_identity() {
        return mask;
    }
    let mut rng = StdRng::seed_from_u64(seed);
    let frac = cfg.fraction.clamp(0.0, 1.0);
    for (idx, is_sen) in mask.iter_mut().enumerate() {
        if !grid.cells[idx].is_tumor {
            continue;
        }
        if rng.gen::<f64>() < frac {
            *is_sen = true;
            apply_senescence_to_cell(&mut grid.cells[idx].cell, cfg);
        }
    }
    mask
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::biochem::sim_cell;
    use crate::cell::{gen_cell, Phenotype, Treatment};
    use crate::params::Params;

    fn grid() -> TumorGrid3D {
        TumorGrid3D::generate(40, 40, 40, 20.0, 42)
    }

    #[test]
    fn identity_is_a_no_op() {
        let mut g = grid();
        let baseline = grid();
        let mask = apply_senescence_program_3d(&mut g, &SenescenceConfig::default(), 7);
        assert!(mask.iter().all(|&s| !s), "identity marks no cells");
        for (a, b) in baseline.cells.iter().zip(g.cells.iter()) {
            assert_eq!(a.cell.iron, b.cell.iron);
            assert_eq!(a.cell.gpx4, b.cell.gpx4);
        }
        assert!(SenescenceConfig::default().is_identity());
        assert!(!SenescenceConfig::literature().is_identity());
    }

    #[test]
    fn marking_is_tumor_only_and_deterministic() {
        let mut a = grid();
        let mut b = grid();
        let cfg = SenescenceConfig::literature();
        let ma = apply_senescence_program_3d(&mut a, &cfg, 123);
        let mb = apply_senescence_program_3d(&mut b, &cfg, 123);
        assert_eq!(ma, mb, "same seed ⇒ same mask");
        let n_sen = ma.iter().filter(|&&s| s).count();
        assert!(n_sen > 0, "literature fraction marks some cells");
        // Only tumor cells are ever senescent.
        for (idx, &is_sen) in ma.iter().enumerate() {
            if is_sen {
                assert!(
                    a.cells[idx].is_tumor,
                    "a senescent cell must be a tumor cell"
                );
            }
        }
        // A different seed gives a different layout.
        let mut c = grid();
        let mc = apply_senescence_program_3d(&mut c, &cfg, 999);
        assert_ne!(ma, mc, "different seed ⇒ different senescent set");
    }

    /// #341: under the model's RSL3 (a DIRECT GPX4 inhibitor), the senescence
    /// program's raised labile iron dominates the residual GPX4 crutch (RSL3
    /// removes ~92% of GPX4), so a senescent cell is a SENOLYTIC target: it
    /// peroxidizes under RSL3. The verified literature is genuinely BIDIRECTIONAL
    /// and the net is parameter-sensitive, so rather than force one sign we prove
    /// the model can express BOTH literature-supported directions under the same
    /// RSL3 (GPX4-inhibition) inducer:
    ///   - an iron-dominant senescent state (high labile Fe2+, no extra backup)
    ///     is a SENOLYTIC target: more lipid peroxidation (D'Ambrosio/Gil 2026,
    ///     PMID 42032311; the GPX4-dependency node);
    ///   - a defense-dominant senescent state (compartmentalized iron, raised
    ///     antioxidant capacity) RESISTS: less lipid peroxidation (Masaldan 2018
    ///     PMID 28888202; Loo 2025 PMID 40731111; the upstream-trigger node).
    /// The axis structure is the result; the net sign is set by which axis the
    /// applied therapy defeats, which is exactly the contested biology.
    #[test]
    fn senescence_program_expresses_both_senolytic_and_resistant_directions() {
        let final_lp = |cfg: &SenescenceConfig| -> f64 {
            let mut gen_rng = StdRng::seed_from_u64(42);
            let mut cell = gen_cell(Phenotype::OXPHOS, &mut gen_rng);
            apply_senescence_to_cell(&mut cell, cfg);
            let mut rng = StdRng::seed_from_u64(7);
            let (_dead, lp, _, _) = sim_cell(&cell, Treatment::RSL3, &Params::default(), &mut rng);
            lp
        };
        let base = final_lp(&SenescenceConfig::default());
        // Iron-dominant (high labile Fe2+, GPX4 crutch removed by RSL3, no backup
        // boost) ⇒ senolytic: MORE peroxidation.
        let senolytic = SenescenceConfig {
            fraction: 1.0,
            iron_mul: 3.0,
            gpx4_mul: 1.0,
            nrf2_mul: 1.0,
            fsp1_mul: 1.0,
        };
        // Defense-dominant (iron compartmentalized ⇒ not raised, antioxidant
        // capacity up) ⇒ resist: LESS peroxidation.
        let resistant = SenescenceConfig {
            fraction: 1.0,
            iron_mul: 1.0,
            gpx4_mul: 1.5,
            nrf2_mul: 1.5,
            fsp1_mul: 1.5,
        };
        assert!(
            final_lp(&senolytic) > base,
            "iron-dominant senescent state should be a senolytic target under RSL3: \
             senolytic={} vs base={base}",
            final_lp(&senolytic)
        );
        assert!(
            final_lp(&resistant) < base,
            "defense-dominant senescent state should resist RSL3: resistant={} vs base={base}",
            final_lp(&resistant)
        );
    }
}
