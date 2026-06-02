//! Clonal heterogeneity (#242).
//!
//! Real tumors are genetic mosaics: 4–8+ subclones with measurably different
//! ferroptosis vulnerabilities coexist, spatially organized into patches.
//! High-mesenchymal (ZEB1+) subclones are highly ferroptosis-vulnerable;
//! epithelial subclones are resistant, and the between-subclone variance often
//! exceeds the between-treatment variance in real drug screens (Conrad et al.,
//! Nat Rev Mol Cell Biol 2018; Marusyk et al., Cancer Cell 2014; Heindl et al.,
//! Nat Methods 2019 for spatial clonal mosaics).
//!
//! ## Design: separate field + independent RNG (byte-identity)
//!
//! [`assign_subclones_3d`] returns a consumer-owned `Vec<u8>` of per-cell
//! subclone ids (`0` = non-tumor, `1..=k` = subclone) computed by Voronoi
//! assignment from `k` seed points. It draws from its **own** `StdRng`, so it
//! never perturbs [`TumorGrid3D::generate`](crate::grid::TumorGrid3D::generate)'s
//! RNG stream — the generated cell grid (phenotypes, persister clusters) is
//! bit-for-bit unchanged whether or not subclones are assigned. The consumer
//! applies [`SubclonePerturbation`]s as RNG-neutral setup mutations (like the
//! `oxygen`/`ph` gradients). A `k = 1` identity config is therefore a no-op
//! and byte-identical to having no clonal model.

use crate::biochem::CellState;
use crate::cell::{gen_cell, Cell, Treatment};
use crate::grid::{GridCell, TumorGrid3D, TUMOR_RADIUS_FRACTION};
use crate::params::Params;
use rand::prelude::*;

/// Per-subclone perturbation of the baseline ferroptosis parameters, applied
/// by the consumer as one-time setup mutations (RNG-neutral). All-identity
/// (the [`SubclonePerturbation::identity`] value) leaves the cell unchanged.
#[derive(Clone, Copy, Debug)]
pub struct SubclonePerturbation {
    /// Multiplier on labile iron (>1 ⇒ more Fenton ROS ⇒ more vulnerable).
    pub iron_mul: f64,
    /// Multiplier on antioxidant/GPX4 capacity (<1 ⇒ less repair ⇒ more
    /// vulnerable). The consumer scales **both** the initial `state.gpx4`
    /// reserve **and** the static `cell.nrf2` setpoint that GPX4 relaxes toward
    /// (`gpx4_target = nrf2 · gpx4_nrf2_target_multiplier`), so the axis is
    /// **durable** across the run rather than an early-window transient that
    /// relaxes back to the shared NRF2 setpoint (#266; the GPX4-transience
    /// finding from the #265 review). Because `cell.nrf2` is the master
    /// antioxidant regulator it also drives GSH resynthesis, so this axis is
    /// deliberately "general antioxidant capacity" (NRF2-low subclone), not
    /// GPX4-reserve in isolation — biologically coherent for a subclone
    /// identity.
    pub gpx4_mul: f64,
    /// Multiplier on oxidizable-PUFA content `cell.lipid_unsat` — the MUFA
    /// membrane-remodeling axis (<1 ⇒ MUFA-enriched ⇒ less peroxidizable lipid
    /// ⇒ more resistant). Scales a **static** `Cell` field so the effect is
    /// durable across steps; perturbing the homeostatically-reset
    /// `state.mufa_protection` instead would be silently overwritten on step 1
    /// under the default params (where `scd_mufa_max == 0`).
    pub lipid_unsat_mul: f64,
}

impl SubclonePerturbation {
    /// The no-op perturbation: iron ×1, GPX4 ×1, lipid_unsat ×1.
    pub fn identity() -> Self {
        SubclonePerturbation {
            iron_mul: 1.0,
            gpx4_mul: 1.0,
            lipid_unsat_mul: 1.0,
        }
    }

    /// True when applying this perturbation changes nothing.
    pub fn is_identity(&self) -> bool {
        self.iron_mul == 1.0 && self.gpx4_mul == 1.0 && self.lipid_unsat_mul == 1.0
    }

    /// Apply this perturbation to a tumor cell + its **already-initialized**
    /// state, as one RNG-neutral setup mutation (the `oxygen`/`ph` consumer
    /// pattern). Scales, in order:
    /// - `cell.iron` ← `iron_mul` (static ⇒ durable)
    /// - `cell.lipid_unsat` ← `lipid_unsat_mul` (static ⇒ durable; MUFA axis)
    /// - `state.gpx4` ← `gpx4_mul` (the **initial** reserve)
    /// - `cell.nrf2` ← `gpx4_mul` (the **static setpoint** GPX4 relaxes toward,
    ///   #266 — this is what makes the antioxidant axis durable instead of an
    ///   early-window transient; it also scales GSH resynthesis, see the
    ///   [`gpx4_mul`](Self::gpx4_mul) field doc).
    ///
    /// [`identity`](Self::identity) is a no-op (all multipliers 1.0), so a K=1
    /// identity config leaves the cell byte-identical. Lives here (rather than
    /// inline in the consumer) so the full set of scaled fields — including the
    /// durable `nrf2` axis — is unit-testable.
    pub fn apply(&self, cell: &mut Cell, state: &mut CellState) {
        cell.iron *= self.iron_mul;
        cell.lipid_unsat *= self.lipid_unsat_mul;
        state.gpx4 *= self.gpx4_mul; // initial reserve
        cell.nrf2 *= self.gpx4_mul; // durable setpoint (#266) — also scales GSH resynthesis
    }
}

/// Clonal-heterogeneity configuration: one [`SubclonePerturbation`] per
/// subclone (`perturbations.len() == k`; subclone ids are `1..=k`).
#[derive(Clone, Debug)]
pub struct ClonalConfig {
    pub perturbations: Vec<SubclonePerturbation>,
    /// Per-step probability that a dead tumor site is repopulated by a living
    /// Moore-neighbor's subclone — spatial **clonal expansion** (#266 item 3).
    /// Resistant subclones, having more living neighbors, win more vacated
    /// sites and grow their territory. `0.0` (the default for every
    /// constructor) ⇒ no repopulation ⇒ the static-map behavior, byte-identical.
    /// Set via [`with_repopulation`](Self::with_repopulation).
    pub repopulation_rate: f64,
}

impl ClonalConfig {
    /// Number of subclones (`k`).
    pub fn k(&self) -> usize {
        self.perturbations.len()
    }

    /// A single identity subclone (`k = 1`, no perturbation). A run with this
    /// config is byte-identical to one with no clonal model — the K=1
    /// byte-identity guarantee.
    pub fn single_identity() -> Self {
        ClonalConfig {
            perturbations: vec![SubclonePerturbation::identity()],
            repopulation_rate: 0.0,
        }
    }

    /// Builder: enable spatial clonal expansion at the given per-step
    /// repopulation probability (#266 item 3). `0.0` leaves it off.
    pub fn with_repopulation(mut self, rate: f64) -> Self {
        self.repopulation_rate = rate;
        self
    }

    /// Literature-informed 4-subclone table (placeholders pending calibration).
    /// Spans the mesenchymal⇄epithelial ferroptosis-vulnerability axis. With the
    /// durable `gpx4_mul` (#266), the GPX4/antioxidant ordering below now holds
    /// for the whole run, not just the early window:
    /// - 1 high-mesenchymal (ZEB1+): iron-loaded, antioxidant-low ⇒ most vulnerable.
    /// - 2 intermediate-mesenchymal: mildly vulnerable.
    /// - 3 intermediate-epithelial: mildly resistant (MUFA-enriched).
    /// - 4 epithelial: antioxidant-high, MUFA-enriched (low PUFA) ⇒ most resistant.
    ///
    /// **Re-check pending (#266 calibration).** These multipliers were chosen
    /// for qualitative direction under the *old transient* GPX4 axis. The
    /// durable axis amplifies their effect at the resistant end — a `gpx4_mul`
    /// > 1 now raises both the GPX4 setpoint *and* GSH resynthesis for the whole
    /// run (compounding on phenotypes that already start NRF2-high, e.g.
    /// `PersisterNrf2`), so the effective between-subclone spread is wider than
    /// when these values were picked. Treat the current spread as illustrative,
    /// not intentional, until calibrated against drug-screen kill-rate spreads
    /// (Conrad 2018; Viswanathan 2017) — tracked as item 2 of #266.
    pub fn literature_4() -> Self {
        ClonalConfig {
            perturbations: vec![
                SubclonePerturbation {
                    iron_mul: 1.5,
                    gpx4_mul: 0.6,
                    lipid_unsat_mul: 1.0,
                },
                SubclonePerturbation {
                    iron_mul: 1.2,
                    gpx4_mul: 0.85,
                    lipid_unsat_mul: 1.0,
                },
                SubclonePerturbation {
                    iron_mul: 0.9,
                    gpx4_mul: 1.1,
                    lipid_unsat_mul: 0.9,
                },
                SubclonePerturbation {
                    iron_mul: 0.7,
                    gpx4_mul: 1.3,
                    lipid_unsat_mul: 0.8,
                },
            ],
            repopulation_rate: 0.0,
        }
    }

    /// True when every subclone perturbation is the identity (no effect).
    pub fn is_identity(&self) -> bool {
        self.perturbations.iter().all(|p| p.is_identity())
    }
}

/// Assign each grid cell to a subclone via Voronoi tessellation from `k` seed
/// points placed uniformly in the tumor sphere. Returns a `Vec<u8>` of length
/// `grid.cells.len()`: `0` for non-tumor cells, `1..=k` for tumor cells (the
/// id of the nearest seed). Deterministic given `(grid dims, k, seed)`.
///
/// Uses an **independent** `StdRng(seed)` so it never advances the RNG used by
/// [`TumorGrid3D::generate`], preserving byte-identity of the cell grid.
///
/// # Panics
/// If `k == 0` or `k > 255` (ids must fit in `u8`, `0` reserved for stroma).
pub fn assign_subclones_3d(grid: &TumorGrid3D, k: usize, seed: u64) -> Vec<u8> {
    assert!(
        k >= 1 && k <= u8::MAX as usize,
        "k must be in 1..=255, got {k}"
    );
    let mut rng = StdRng::seed_from_u64(seed);
    let center = (
        grid.rows as f64 / 2.0,
        grid.cols as f64 / 2.0,
        grid.layers as f64 / 2.0,
    );
    let tumor_radius = (grid.rows.min(grid.cols).min(grid.layers) as f64) * TUMOR_RADIUS_FRACTION;

    // Seed points uniformly distributed in the tumor sphere (cbrt radial
    // sampling to avoid center bias — same convention as generate's clusters).
    let seeds: Vec<(f64, f64, f64)> = (0..k)
        .map(|_| {
            let dist = rng.gen::<f64>().cbrt() * tumor_radius * 0.9;
            let theta = rng.gen::<f64>() * std::f64::consts::TAU;
            let cos_phi = 2.0 * rng.gen::<f64>() - 1.0;
            let sin_phi = (1.0 - cos_phi * cos_phi).sqrt();
            (
                center.0 + dist * cos_phi,
                center.1 + dist * sin_phi * theta.cos(),
                center.2 + dist * sin_phi * theta.sin(),
            )
        })
        .collect();

    (0..grid.cells.len())
        .map(|idx| {
            if !grid.cells[idx].is_tumor {
                return 0u8;
            }
            let (r, c, l) = grid.coords(idx);
            let (rf, cf, lf) = (r as f64, c as f64, l as f64);
            let mut best = 0usize;
            let mut best_d2 = f64::INFINITY;
            for (si, &(sr, sc, sl)) in seeds.iter().enumerate() {
                let d2 = (rf - sr).powi(2) + (cf - sc).powi(2) + (lf - sl).powi(2);
                if d2 < best_d2 {
                    best_d2 = d2;
                    best = si;
                }
            }
            (best + 1) as u8
        })
        .collect()
}

/// Spatial clonal **expansion** (#266 item 3): repopulate dead tumor sites from
/// living Moore-26 neighbors, so resistant subclones (which have more living
/// neighbors) grow their territory over the run. The revived site becomes a
/// fresh living cell of a randomly-chosen living neighbor's subclone (its
/// phenotype + that subclone's [`SubclonePerturbation`]), and `subclone_ids` is
/// updated so the expansion shows in the per-subclone census.
///
/// **Deterministic + order-independent**: a start-of-pass living snapshot is
/// taken first, all repopulations are decided against it, then applied — so a
/// site revived this pass cannot itself seed another this step. Each dead site
/// uses a per-site RNG mixed from `(seed, idx, step)`; pass a condition-derived
/// `seed` (like the immune-kill RNG) so it never touches the setup-seed streams.
///
/// No-op (returns 0, mutates nothing) when `cfg.repopulation_rate <= 0.0`, so a
/// consumer that doesn't opt in keeps the static-map behavior, byte-identical.
/// Returns the number of sites repopulated this step.
pub fn repopulate_dead_sites_3d(
    grid: &mut TumorGrid3D,
    subclone_ids: &mut [u8],
    cfg: &ClonalConfig,
    params: &Params,
    seed: u64,
    step: u32,
) -> usize {
    if cfg.repopulation_rate <= 0.0 {
        return 0;
    }
    let n = grid.cells.len();
    // Phase 1: snapshot living tumor sites (start of pass).
    let living: Vec<bool> = (0..n)
        .map(|i| grid.cells[i].is_tumor && !grid.cells[i].state.dead)
        .collect();
    // Phase 2: decide + build each repopulation against the snapshot.
    let mut repop: Vec<(usize, u8, GridCell)> = Vec::new();
    for idx in 0..n {
        // Only DEAD tumor sites are candidates for repopulation.
        if !grid.cells[idx].is_tumor || living[idx] {
            continue;
        }
        let mut rng = StdRng::seed_from_u64(
            seed ^ (idx as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15)
                ^ (step as u64).wrapping_mul(0x8000_0000_0000_01B3),
        );
        if rng.gen::<f64>() >= cfg.repopulation_rate {
            continue;
        }
        // Living tumor Moore-26 neighbors (from the snapshot).
        let (r, c, l) = grid.coords(idx);
        let (nbrs, cnt) = grid.neighbors(r, c, l);
        let donors: Vec<usize> = nbrs[..cnt]
            .iter()
            .map(|&(nr, nc, nl)| grid.flat_index(nr, nc, nl))
            .filter(|&ni| living[ni])
            .collect();
        if donors.is_empty() {
            continue; // isolated dead site — nothing to repopulate from
        }
        let donor = donors[rng.gen_range(0..donors.len())];
        let donor_id = subclone_ids[donor];
        let pheno = grid.cells[donor].phenotype;
        // A fresh living daughter cell of the donor subclone.
        let mut cell = gen_cell(pheno, &mut rng);
        let mut state = CellState::from_cell(&cell, Treatment::Control, params, &mut rng);
        cfg.perturbations[(donor_id - 1) as usize].apply(&mut cell, &mut state);
        repop.push((
            idx,
            donor_id,
            GridCell {
                cell,
                phenotype: pheno,
                state,
                is_tumor: true,
                extra_iron: 0.0,
                lp_at_death: 0.0,
                newly_dead: false,
            },
        ));
    }
    // Phase 3: apply.
    let count = repop.len();
    for (idx, id, gc) in repop {
        grid.cells[idx] = gc;
        subclone_ids[idx] = id;
    }
    count
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn identity_perturbation_is_noop() {
        assert!(SubclonePerturbation::identity().is_identity());
        assert!(ClonalConfig::single_identity().is_identity());
    }

    #[test]
    fn literature_4_is_not_identity_and_spans_vulnerability() {
        let c = ClonalConfig::literature_4();
        assert_eq!(c.k(), 4);
        assert!(!c.is_identity());
        // Subclone 1 is more vulnerable than subclone 4 on every axis: more
        // iron, less GPX4, and more oxidizable PUFA (higher lipid_unsat).
        let v = c.perturbations[0];
        let r = c.perturbations[3];
        assert!(
            v.iron_mul > r.iron_mul
                && v.gpx4_mul < r.gpx4_mul
                && v.lipid_unsat_mul > r.lipid_unsat_mul
        );
    }

    #[test]
    fn assignment_is_deterministic_given_seed() {
        let g = TumorGrid3D::generate(20, 20, 20, 20.0, 42);
        let a = assign_subclones_3d(&g, 4, 7);
        let b = assign_subclones_3d(&g, 4, 7);
        assert_eq!(a, b);
        // A different seed generally moves at least one boundary cell.
        let c = assign_subclones_3d(&g, 4, 8);
        assert_ne!(a, c);
    }

    #[test]
    fn ids_are_zero_for_stroma_and_in_range_for_tumor() {
        let g = TumorGrid3D::generate(20, 20, 20, 20.0, 42);
        let ids = assign_subclones_3d(&g, 4, 7);
        assert_eq!(ids.len(), g.cells.len());
        for (idx, &id) in ids.iter().enumerate() {
            if g.cells[idx].is_tumor {
                assert!((1..=4).contains(&id), "tumor cell id {id} out of 1..=4");
            } else {
                assert_eq!(id, 0, "non-tumor cell must be subclone 0");
            }
        }
    }

    #[test]
    fn k1_assigns_every_tumor_cell_to_subclone_one() {
        let g = TumorGrid3D::generate(20, 20, 20, 20.0, 42);
        let ids = assign_subclones_3d(&g, 1, 7);
        for (idx, &id) in ids.iter().enumerate() {
            assert_eq!(id, u8::from(g.cells[idx].is_tumor));
        }
    }

    #[test]
    #[should_panic(expected = "k must be in 1..=255")]
    fn k_zero_panics() {
        let g = TumorGrid3D::generate(10, 10, 10, 20.0, 42);
        let _ = assign_subclones_3d(&g, 0, 7);
    }

    /// #266 durable GPX4 axis: the consumer scales BOTH the initial
    /// `state.gpx4` reserve and the static `cell.nrf2` setpoint. This test pins
    /// the durability the fix is for — that the GPX4/antioxidant difference
    /// persists late in the run rather than relaxing back to the shared NRF2
    /// setpoint.
    ///
    /// Three variants of the same baseline cell, run to step 160 under Control
    /// (no treatment ⇒ the cell stays healthy and `gpx4` relaxes cleanly toward
    /// `nrf2 · gpx4_nrf2_target_multiplier`):
    /// - **identity**: no perturbation.
    /// - **transient** (the OLD behavior): scale only the initial `state.gpx4`.
    /// - **durable** (#266): scale the initial `state.gpx4` AND `cell.nrf2`.
    ///
    /// The transient knock relaxes back UP toward identity's setpoint, while the
    /// durable knock holds a lower setpoint — so the durable axis keeps a much
    /// larger late-run differentiation.
    #[test]
    fn gpx4_axis_persists_late_only_when_nrf2_is_scaled() {
        use crate::biochem::{sim_cell_step, CellState};
        use crate::cell::{gen_cell, Phenotype, Treatment};
        use crate::params::Params;
        use rand::prelude::*;

        let params = Params::default();
        let mul = 0.5; // a "GPX4-low / antioxidant-low" subclone
        let steps = 160;

        // One baseline cell, cloned so the only difference is the perturbation.
        let mut gen_rng = StdRng::seed_from_u64(1);
        let base = gen_cell(Phenotype::OXPHOS, &mut gen_rng);

        // Run a variant to step `steps`, returning the late-run gpx4. Each run
        // uses the SAME step-RNG seed so the only difference is the scaling.
        let run = |scale_gpx4_init: bool, scale_nrf2: bool| -> f64 {
            let mut cell = base.clone();
            if scale_nrf2 {
                cell.nrf2 *= mul; // durable setpoint (#266)
            }
            let mut init_rng = StdRng::seed_from_u64(7);
            let mut state = CellState::from_cell(&cell, Treatment::Control, &params, &mut init_rng);
            if scale_gpx4_init {
                state.gpx4 *= mul; // initial reserve (transient on its own)
            }
            let mut step_rng = StdRng::seed_from_u64(99);
            for step in 0..steps {
                sim_cell_step(&mut state, &cell, &params, step, 0.0, &mut step_rng);
            }
            state.gpx4
        };

        let identity = run(false, false);
        let transient = run(true, false); // old behavior: initial-only
        let durable = run(true, true); // #266: initial + nrf2 setpoint

        // The durable knock holds a meaningfully lower late-run gpx4 (its
        // setpoint is halved), well below the identity cell.
        assert!(
            durable < 0.75 * identity,
            "durable axis should stay well below identity late: durable={durable}, identity={identity}"
        );
        // The transient knock relaxed back UP, ending above the durable one.
        assert!(
            transient > durable,
            "transient (initial-only) knock should relax back above the durable one: transient={transient}, durable={durable}"
        );
        // Quantitatively: the durable axis maintains a far larger late-run
        // differentiation from identity than the transient one does.
        let gap_durable = identity - durable;
        let gap_transient = identity - transient;
        assert!(
            gap_durable > 2.0 * gap_transient,
            "durable differentiation should dwarf the transient one: gap_durable={gap_durable}, gap_transient={gap_transient}"
        );
    }

    /// `apply` scales every axis the consumer relies on — crucially including
    /// the durable `cell.nrf2` setpoint (#266). This pins the exact field set so
    /// dropping the `nrf2` scaling (the durability fix) fails a fast unit test,
    /// which a sim-level kill comparison can't catch (the initial-gpx4 knock
    /// alone already produces a kill differential). Identity must be a no-op.
    #[test]
    fn apply_scales_all_axes_including_durable_nrf2() {
        use crate::biochem::CellState;
        use crate::cell::{gen_cell, Phenotype, Treatment};
        use crate::params::Params;
        use rand::prelude::*;

        let params = Params::default();
        let mut rng = StdRng::seed_from_u64(3);
        let mut cell = gen_cell(Phenotype::Glycolytic, &mut rng);
        let mut state = CellState::from_cell(&cell, Treatment::Control, &params, &mut rng);
        let (iron0, lipid0, nrf20, gpx40) = (cell.iron, cell.lipid_unsat, cell.nrf2, state.gpx4);

        let p = SubclonePerturbation {
            iron_mul: 2.0,
            gpx4_mul: 0.5,
            lipid_unsat_mul: 0.25,
        };
        p.apply(&mut cell, &mut state);

        assert!((cell.iron - iron0 * 2.0).abs() < 1e-12, "iron scaled");
        assert!(
            (cell.lipid_unsat - lipid0 * 0.25).abs() < 1e-12,
            "lipid_unsat scaled"
        );
        assert!(
            (state.gpx4 - gpx40 * 0.5).abs() < 1e-12,
            "initial gpx4 reserve scaled"
        );
        // The #266 durable axis: nrf2 (setpoint) scaled by gpx4_mul.
        assert!(
            (cell.nrf2 - nrf20 * 0.5).abs() < 1e-12,
            "durable nrf2 setpoint scaled by gpx4_mul (#266)"
        );

        // Identity is a no-op.
        let mut c2 = gen_cell(Phenotype::Glycolytic, &mut rng);
        let mut s2 = CellState::from_cell(&c2, Treatment::Control, &params, &mut rng);
        let (i2, n2, g2) = (c2.iron, c2.nrf2, s2.gpx4);
        SubclonePerturbation::identity().apply(&mut c2, &mut s2);
        assert_eq!(c2.iron, i2);
        assert_eq!(c2.nrf2, n2);
        assert_eq!(s2.gpx4, g2);
    }

    // ===== Spatial clonal expansion / repopulation (#266 item 3) =====

    fn killed_plane(seed: u64) -> (TumorGrid3D, Vec<u8>, Vec<usize>) {
        let mut g = TumorGrid3D::generate(20, 20, 20, 20.0, 42);
        let ids = assign_subclones_3d(&g, 4, 7);
        // Kill the central r=10 plane of tumor cells (living neighbors remain
        // at r=9/11), leaving dead sites with donors to repopulate from.
        let mut killed = Vec::new();
        for idx in 0..g.cells.len() {
            let (r, _, _) = g.coords(idx);
            if r == 10 && g.cells[idx].is_tumor {
                g.cells[idx].state.dead = true;
                killed.push(idx);
            }
        }
        let _ = seed;
        (g, ids, killed)
    }

    #[test]
    fn repopulate_is_noop_when_rate_zero() {
        let (mut g, mut ids, killed) = killed_plane(0);
        let ids_before = ids.clone();
        let cfg = ClonalConfig::literature_4(); // repopulation_rate = 0
        let n = repopulate_dead_sites_3d(&mut g, &mut ids, &cfg, &Params::default(), 1, 0);
        assert_eq!(n, 0, "no repopulation when rate is 0");
        assert_eq!(ids, ids_before, "subclone map unchanged");
        assert!(
            killed.iter().all(|&i| g.cells[i].state.dead),
            "dead sites stay dead"
        );
    }

    #[test]
    fn repopulate_revives_dead_sites_and_is_deterministic() {
        let cfg = ClonalConfig::literature_4().with_repopulation(1.0);
        let (mut g, mut ids, killed) = killed_plane(0);
        assert!(!killed.is_empty(), "the killed plane is non-empty");
        let n = repopulate_dead_sites_3d(&mut g, &mut ids, &cfg, &Params::default(), 5, 0);
        assert!(
            n > 0,
            "rate 1.0 repopulates dead sites that have living neighbors"
        );
        // Revived sites are alive again with a valid subclone id (1..=4).
        let revived = killed.iter().filter(|&&i| !g.cells[i].state.dead).count();
        assert_eq!(revived, n, "revived count matches the return value");
        for &i in &killed {
            if !g.cells[i].state.dead {
                assert!(
                    (1..=4).contains(&ids[i]),
                    "revived site inherits a subclone"
                );
            }
        }
        // Deterministic: same seeds ⇒ identical decisions.
        let (mut g2, mut ids2, _) = killed_plane(0);
        let n2 = repopulate_dead_sites_3d(&mut g2, &mut ids2, &cfg, &Params::default(), 5, 0);
        assert_eq!(n, n2, "deterministic repopulation count");
        assert_eq!(ids, ids2, "deterministic subclone map after repopulation");
    }
}
