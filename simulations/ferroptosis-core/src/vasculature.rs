//! Explicit 3D tumor vasculature (#191).
//!
//! The 2D simulations (and the 3D `oxygen::radial_o2_field`) use "distance
//! from the tumor edge" as a proxy for vasculature: the spheroid surface is
//! the only O2/drug source, so supply decays monotonically inward. Real 3D
//! tumors instead carry **internal** vessels that penetrate the volume, so
//! oxygenation is patchy: well-supplied near a vessel, hypoxic in the gaps
//! between them. Inter-vessel distance is ~100–200 µm in well-vascularized
//! tumors and ~300–500 µm in poorly-vascularized ones (Vaupel; Krogh-cylinder
//! supply geometry).
//!
//! This module implements **Option A** from #191: random vessel seed points in
//! the tumor volume, with each cell's supply set by the Krogh-style decay
//! `exp(-distance_to_nearest_vessel / λ)`. It is a drop-in alternative to
//! `oxygen::radial_o2_field` — same `Vec<f64>` per-cell-factor output, same
//! "non-tumor cells = 1.0 (well-perfused bulk)" convention — so the consumer
//! multiplies it into `cell.basal_ros` (O2) or a drug-availability field
//! identically. The same proximity factor supplies both O2 and drug.
//!
//! ## Design: independent RNG (byte-identity)
//!
//! [`place_vessels_3d`] draws vessel positions from its **own** `StdRng`, so it
//! never perturbs [`TumorGrid3D::generate`](crate::grid::TumorGrid3D::generate)'s
//! stream — the cell grid is bit-for-bit unchanged whether or not vessels are
//! placed. A consumer that doesn't opt into the vessel model keeps using
//! `radial_o2_field` and stays byte-identical.

use crate::grid::{TumorGrid3D, TUMOR_RADIUS_FRACTION};
use rand::prelude::*;

/// How vessel positions are laid out.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum VesselTopology {
    /// Uniform-random seed points in the tumor volume (the original #191 model).
    #[default]
    Random,
    /// Fractal-branching vascular tree (#268): trunks enter from the periphery
    /// and bifurcate inward with high, tumor-like variability — a hierarchical
    /// but chaotic network with dead ends and perfusion gaps, unlike the smooth
    /// coverage of uniform-random points. Motivated by the well-documented
    /// fractal/irregular architecture of tumor vasculature (Baish & Jain,
    /// *Cancer Res* 2000, PMID 10919633): tumor vessels are disorganized, with
    /// high variability in segment length and branching angle and a higher,
    /// space-filling fractal dimension (~1.89) than normal vasculature (~1.70).
    ///
    /// Note: only the central-sphere / spheroid placement path honors this; the
    /// patient-scale slab path ([`place_vessels_in_slab_3d`], #240/#272) always
    /// uses uniform-in-box placement and ignores `topology`.
    Fractal,
}

/// Vessel-network configuration. `inter_vessel_um` is the target mean spacing
/// between vessel seed points; the vessel count is derived from it and the
/// tumor volume in [`place_vessels_3d`]. The Krogh decay length λ is supplied
/// per call (it matches the condition's O2 reference λ, like `radial_o2_field`).
#[derive(Clone, Copy, Debug)]
pub struct VasculatureConfig {
    /// Target mean inter-vessel spacing (µm). Smaller ⇒ denser ⇒ better-oxygenated.
    pub inter_vessel_um: f64,
    /// Layout: uniform-random (default) or fractal-branching tree (#268).
    pub topology: VesselTopology,
}

impl VasculatureConfig {
    /// Well-vascularized tumor (~150 µm inter-vessel spacing).
    pub fn well_vascularized() -> Self {
        VasculatureConfig {
            inter_vessel_um: 150.0,
            topology: VesselTopology::Random,
        }
    }

    /// Poorly-vascularized tumor (~400 µm inter-vessel spacing).
    pub fn poorly_vascularized() -> Self {
        VasculatureConfig {
            inter_vessel_um: 400.0,
            topology: VesselTopology::Random,
        }
    }

    /// Switch this config to the fractal-branching topology (#268).
    ///
    /// Footgun: this is a no-op for slab geometry — [`place_vessels_in_slab_3d`]
    /// (#240/#272) always scatters vessels uniform-in-box regardless of
    /// `topology`. Only the central-sphere / spheroid path branches on it.
    pub fn with_fractal(mut self) -> Self {
        self.topology = VesselTopology::Fractal;
        self
    }
}

/// Vessel count from tumor volume and target spacing (`n ≈ volume / spacing³`,
/// cubic packing, floored at 1). Shared by [`place_vessels_3d`] (sphere) and
/// [`place_vessels_in_slab_3d`] (box) so the two placements agree on density.
fn derive_vessel_count(grid: &TumorGrid3D, cfg: &VasculatureConfig) -> usize {
    debug_assert!(
        cfg.inter_vessel_um.is_finite() && cfg.inter_vessel_um > 0.0,
        "inter_vessel_um must be finite and positive, got {}",
        cfg.inter_vessel_um
    );
    let cell_um3 = grid.cell_size_um.powi(3);
    let n_tumor = grid.cells.iter().filter(|gc| gc.is_tumor).count();
    let tumor_volume_um3 = n_tumor as f64 * cell_um3;
    (tumor_volume_um3 / cfg.inter_vessel_um.powi(3))
        .round()
        .max(1.0) as usize
}

/// Place vessel seed points uniformly in the tumor sphere. The count is an
/// approximation from the tumor volume and target inter-vessel spacing
/// (`n ≈ tumor_volume / inter_vessel³`, assuming cubic packing, floored at 1).
/// Returns positions in **lattice (cell) coordinates**. Deterministic given
/// `(grid dims, cfg, seed)`.
///
/// Uses an **independent** `StdRng(seed)` so it never advances the RNG used by
/// [`TumorGrid3D::generate`], preserving byte-identity of the cell grid.
pub fn place_vessels_3d(
    grid: &TumorGrid3D,
    cfg: &VasculatureConfig,
    seed: u64,
) -> Vec<(f64, f64, f64)> {
    let n_vessels = derive_vessel_count(grid, cfg);

    let mut rng = StdRng::seed_from_u64(seed);
    let center = (
        grid.rows as f64 / 2.0,
        grid.cols as f64 / 2.0,
        grid.layers as f64 / 2.0,
    );
    let tumor_radius = (grid.rows.min(grid.cols).min(grid.layers) as f64) * TUMOR_RADIUS_FRACTION;

    // Uniform-in-sphere sampling (cbrt radial avoids center bias — same
    // convention as generate's persister clusters / clonal seeds).
    (0..n_vessels)
        .map(|_| {
            let dist = rng.gen::<f64>().cbrt() * tumor_radius * 0.95;
            let theta = rng.gen::<f64>() * std::f64::consts::TAU;
            let cos_phi = 2.0 * rng.gen::<f64>() - 1.0;
            let sin_phi = (1.0 - cos_phi * cos_phi).sqrt();
            (
                center.0 + dist * cos_phi,
                center.1 + dist * sin_phi * theta.cos(),
                center.2 + dist * sin_phi * theta.sin(),
            )
        })
        .collect()
}

/// Place vessel seed points uniformly across the **whole grid box**
/// (uniform-in-box), rather than the central tumor sphere of
/// [`place_vessels_3d`]. For a patient-scale slab (#240, an all-tumor block,
/// not a sphere) vessels should pervade the entire block so deep tissue
/// *throughout* — not just a central pocket — can have focal well-perfused
/// regions (#272 slab+vasculature coupling). Count matches `place_vessels_3d`
/// (same volume/spacing rule) so the two agree on density; only the spatial
/// distribution differs (box vs sphere).
///
/// Returns positions in **lattice (cell) coordinates**. Deterministic given
/// `(grid dims, cfg, seed)`, with an **independent** `StdRng(seed)` so it never
/// advances the grid-generation RNG (byte-identity preserved).
pub fn place_vessels_in_slab_3d(
    grid: &TumorGrid3D,
    cfg: &VasculatureConfig,
    seed: u64,
) -> Vec<(f64, f64, f64)> {
    // Precondition: an all-tumor grid (`TumorGrid3D::generate_slab`). Uniform-in-
    // box sampling over the whole grid only places vessels in tumor space when
    // every cell IS tumor; on a sphere grid (with stroma) it would scatter
    // vessels into non-tumor space while `derive_vessel_count` sizes the count
    // from the tumor volume, under-vascularizing the actual tumor. Use
    // `place_vessels_3d` for sphere grids. (cheap O(cells) check, debug-only.)
    debug_assert!(
        grid.cells.iter().all(|gc| gc.is_tumor),
        "place_vessels_in_slab_3d expects an all-tumor (slab) grid; use place_vessels_3d for sphere grids"
    );
    let n_vessels = derive_vessel_count(grid, cfg);

    let mut rng = StdRng::seed_from_u64(seed);
    (0..n_vessels)
        .map(|_| {
            (
                rng.gen::<f64>() * grid.rows as f64,
                rng.gen::<f64>() * grid.cols as f64,
                rng.gen::<f64>() * grid.layers as f64,
            )
        })
        .collect()
}

/// Place vessels as a **fractal-branching tree** (#268): trunks enter from the
/// tumor periphery and bifurcate inward (breadth-first), returning sampled points
/// along the branches (lattice coordinates) for [`vessel_supply_field`]. Unlike
/// the uniform-random [`place_vessels_3d`], this produces a hierarchical-but-
/// chaotic network — points cluster along branches with avascular gaps and dead
/// ends between them — matching the documented irregular/fractal architecture of
/// tumor vasculature (see [`VesselTopology::Fractal`]; Baish & Jain 2000). The
/// total point count is capped at the same `inter_vessel_um` target as
/// [`place_vessels_3d`] — but note this matches the raw vessel-point COUNT, NOT
/// effective spatial coverage. Because these points are 1-cell-spaced along a
/// few branches (rather than spread through the volume), the fractal network
/// covers far less unique territory at equal count, which is exactly why it
/// leaves more avascular tissue. So the higher hypoxic fraction vs random is a
/// clustering-coverage effect, not a controlled "same density, different
/// topology" result — read it qualitatively. Breadth-first growth + the cap
/// keep the multi-trunk network balanced.
///
/// Uses an **independent** `StdRng(seed)`, so it never advances
/// [`TumorGrid3D::generate`]'s stream — byte-identity preserved. Deterministic
/// given `(grid dims, cfg, seed)`. Tuning constants (branch angle, length ratio,
/// dead-end rate) are not calibrated; they encode the *qualitative* tumor-vessel
/// chaos (high angle/length variability), a #268 follow-up to ground against
/// micro-CT morphometry.
pub fn place_vessels_fractal_3d(
    grid: &TumorGrid3D,
    cfg: &VasculatureConfig,
    seed: u64,
) -> Vec<(f64, f64, f64)> {
    use std::collections::VecDeque;
    use std::f64::consts::TAU;

    let target = derive_vessel_count(grid, cfg);
    let mut rng = StdRng::seed_from_u64(seed);
    let center = (
        grid.rows as f64 / 2.0,
        grid.cols as f64 / 2.0,
        grid.layers as f64 / 2.0,
    );
    let radius = (grid.rows.min(grid.cols).min(grid.layers) as f64) * TUMOR_RADIUS_FRACTION;

    // Tumor vasculature is disorganized — high variability in branch angle and
    // length (Baish & Jain 2000). Base bifurcation ~35° with large jitter;
    // branch length shrinks ~0.8/generation; ~12% of branches dead-end early
    // (perfusion gaps). Trunk count scales with the point-count target at
    // ~1 feeding trunk per 30 vessel-points, clamped to [2, 64] — uncalibrated,
    // chosen so small grids still get ≥2 entry points and patient-scale grids a
    // realistic handful of feeding trunks rather than hundreds.
    let n_trunks = ((target as f64 / 30.0).round() as usize).clamp(2, 64);
    const STEP: f64 = 1.0; // point spacing along a branch (cells)
    const BASE_ANGLE: f64 = 0.61; // ~35°
    const LENGTH_RATIO: f64 = 0.80;
    const MIN_LEN: f64 = 3.0;
    const MAX_DEPTH: u32 = 14;
    const DEAD_END_PROB: f64 = 0.12;

    let norm = |v: (f64, f64, f64)| -> (f64, f64, f64) {
        let n = (v.0 * v.0 + v.1 * v.1 + v.2 * v.2).sqrt().max(1e-9);
        (v.0 / n, v.1 / n, v.2 / n)
    };
    let inside = |p: (f64, f64, f64)| -> bool {
        let (dx, dy, dz) = (p.0 - center.0, p.1 - center.1, p.2 - center.2);
        (dx * dx + dy * dy + dz * dz).sqrt() <= radius
    };
    // A random unit vector perpendicular to `d` (Gram-Schmidt + random roll).
    let perp = |d: (f64, f64, f64), rng: &mut StdRng| -> (f64, f64, f64) {
        let a = if d.0.abs() < 0.9 {
            (1.0, 0.0, 0.0)
        } else {
            (0.0, 1.0, 0.0)
        };
        let adot = a.0 * d.0 + a.1 * d.1 + a.2 * d.2;
        let u = norm((a.0 - adot * d.0, a.1 - adot * d.1, a.2 - adot * d.2));
        let w = (
            d.1 * u.2 - d.2 * u.1,
            d.2 * u.0 - d.0 * u.2,
            d.0 * u.1 - d.1 * u.0,
        );
        let roll = rng.gen::<f64>() * TAU;
        (
            u.0 * roll.cos() + w.0 * roll.sin(),
            u.1 * roll.cos() + w.1 * roll.sin(),
            u.2 * roll.cos() + w.2 * roll.sin(),
        )
    };

    // Trunk roots on the surface (0.95R), each pointing inward.
    let mut queue: VecDeque<((f64, f64, f64), (f64, f64, f64), f64, u32)> = VecDeque::new();
    for _ in 0..n_trunks {
        let theta = rng.gen::<f64>() * TAU;
        let cos_phi = 2.0 * rng.gen::<f64>() - 1.0;
        let sin_phi = (1.0 - cos_phi * cos_phi).sqrt();
        let outward = (cos_phi, sin_phi * theta.cos(), sin_phi * theta.sin());
        let root = (
            center.0 + outward.0 * radius * 0.95,
            center.1 + outward.1 * radius * 0.95,
            center.2 + outward.2 * radius * 0.95,
        );
        queue.push_back((root, (-outward.0, -outward.1, -outward.2), radius * 0.5, 0));
    }

    let mut points: Vec<(f64, f64, f64)> = Vec::with_capacity(target);
    while let Some((start, dir, length, depth)) = queue.pop_front() {
        let dir = norm(dir);
        let n_steps = (length / STEP).round().max(1.0) as usize;
        let mut p = start;
        for _ in 0..n_steps {
            let np = (p.0 + dir.0 * STEP, p.1 + dir.1 * STEP, p.2 + dir.2 * STEP);
            if !inside(np) {
                break;
            }
            points.push(np);
            p = np;
            if points.len() >= target {
                break; // exact density cap, mid-segment
            }
        }
        if points.len() >= target {
            break; // density cap (BFS ⇒ the network stays balanced when capped)
        }
        // Stop branching: short segment, depth cap, or a random dead end.
        if length < MIN_LEN || depth >= MAX_DEPTH || rng.gen::<f64>() < DEAD_END_PROB {
            continue;
        }
        let child_len = length * LENGTH_RATIO * (0.8 + 0.4 * rng.gen::<f64>());
        // One bifurcation PLANE per branch point (a single shared perpendicular
        // `axis`), with the two daughters placed on OPPOSITE sides of the parent
        // via `sign` — a genuine Y-split. Each daughter keeps its own angle
        // jitter so the split stays tumor-irregular rather than perfectly
        // symmetric. (Re-rolling `axis` per daughter would make `sign` inert,
        // collapsing the Y into two independent random forward branches.)
        let axis = perp(dir, &mut rng);
        for sign in [1.0_f64, -1.0] {
            let ang = BASE_ANGLE * (0.4 + 1.2 * rng.gen::<f64>()); // high variability
            let child = norm((
                dir.0 * ang.cos() + sign * axis.0 * ang.sin(),
                dir.1 * ang.cos() + sign * axis.1 * ang.sin(),
                dir.2 * ang.cos() + sign * axis.2 * ang.sin(),
            ));
            queue.push_back((p, child, child_len, depth + 1));
        }
    }

    // vessel_supply_field requires ≥1 vessel; a degenerate tiny grid could yield none.
    if points.is_empty() {
        points.push(center);
    }
    points
}

/// Uniform-grid spatial index over vessel positions for **exact** nearest-vessel
/// distance queries (#268). Vessels are binned into a coarse lattice (sized for
/// ≈1 vessel/bin); a query expands outward in Chebyshev shells and stops as soon
/// as the best distance found is provably ≤ the closest distance any *unsearched*
/// bin could hold. It therefore returns the **same** `nearest_d2` as the
/// brute-force min (an exact acceleration, not an approximation — `min` over the
/// same finite set is order-independent), keeping [`vessel_supply_field`]
/// byte-identical while cutting it from `O(cells × vessels)` to roughly
/// `O(cells)` with a bounded per-query neighborhood.
///
/// Assumes a **physical** inter-vessel spacing (≳ cell size): the bin edge is
/// floored at 1 lattice unit, so a non-physical sub-cell spacing
/// (`inter_vessel_um < cell_size_um`) degenerates to per-cell bins — `nb` grows
/// to the grid dims and `bins` holds ~`dims³` (mostly empty) vectors. Correct,
/// but memory-heavy in that regime; real vasculature spacing (~100–500 µm vs a
/// ~20 µm cell) keeps `bin` at several cells and `nb³` small.
struct VesselIndex {
    /// Bin edge length in lattice units.
    bin: f64,
    /// Bins per axis `(rows, cols, layers)`.
    nb: (usize, usize, usize),
    /// Vessels grouped by flat bin index (row-major over `nb`).
    bins: Vec<Vec<(f64, f64, f64)>>,
}

impl VesselIndex {
    /// Bin index of a lattice coordinate along one axis, clamped to `[0, n-1]`.
    /// Returned as `isize` so the query can offset it by shell deltas; `build`
    /// casts to `usize`. Shared so bin assignment and querying can't drift.
    fn bin_index(coord: f64, bin: f64, n: usize) -> isize {
        ((coord / bin).floor() as isize).clamp(0, n as isize - 1)
    }

    fn build(grid: &TumorGrid3D, vessels: &[(f64, f64, f64)]) -> Self {
        // Target ≈1 vessel/bin: bin volume ≈ total_cells / n_vessels, so the
        // edge ≈ cbrt(that). Clamp to ≥1 lattice unit. Keeps the typical query
        // shell radius small (the nearest vessel is usually 0–2 shells away).
        let total = (grid.rows * grid.cols * grid.layers).max(1) as f64;
        let bin = (total / vessels.len() as f64).cbrt().max(1.0);
        let nb = (
            ((grid.rows as f64 / bin).ceil() as usize).max(1),
            ((grid.cols as f64 / bin).ceil() as usize).max(1),
            ((grid.layers as f64 / bin).ceil() as usize).max(1),
        );
        let mut bins = vec![Vec::new(); nb.0 * nb.1 * nb.2];
        for &v in vessels {
            let br = Self::bin_index(v.0, bin, nb.0) as usize;
            let bc = Self::bin_index(v.1, bin, nb.1) as usize;
            let bl = Self::bin_index(v.2, bin, nb.2) as usize;
            bins[(br * nb.1 + bc) * nb.2 + bl].push(v);
        }
        VesselIndex { bin, nb, bins }
    }

    /// Exact squared distance to the nearest vessel from integer lattice
    /// `(r, c, l)`. Bit-for-bit equal to the brute-force min over all vessels.
    fn nearest_d2(&self, r: usize, c: usize, l: usize) -> f64 {
        let (rf, cf, lf) = (r as f64, c as f64, l as f64);
        let cb = (
            Self::bin_index(rf, self.bin, self.nb.0),
            Self::bin_index(cf, self.bin, self.nb.1),
            Self::bin_index(lf, self.bin, self.nb.2),
        );
        // Largest shell needed to reach every bin from this cell's bin.
        let r_max = (cb.0.max(self.nb.0 as isize - 1 - cb.0))
            .max(cb.1.max(self.nb.1 as isize - 1 - cb.1))
            .max(cb.2.max(self.nb.2 as isize - 1 - cb.2));
        let mut best = f64::INFINITY;
        let mut shell = 0isize;
        while shell <= r_max {
            // Bins at Chebyshev distance == shell (the not-yet-searched ring).
            // This rescans the full (2·shell+1)³ index box and `continue`s the
            // interior, so reaching shell R is O(R⁴) index iterations (most
            // skipped). Negligible at the ≈1-vessel/bin sizing — the early-stop
            // below almost always fires by R≤2 — so a boundary-only enumeration
            // isn't worth the complexity unless a pathologically sparse grid
            // ever profiles hot.
            for dr in -shell..=shell {
                for dc in -shell..=shell {
                    for dl in -shell..=shell {
                        if dr.abs().max(dc.abs()).max(dl.abs()) != shell {
                            continue;
                        }
                        let (br, bc, bl) = (cb.0 + dr, cb.1 + dc, cb.2 + dl);
                        if br < 0
                            || bc < 0
                            || bl < 0
                            || br >= self.nb.0 as isize
                            || bc >= self.nb.1 as isize
                            || bl >= self.nb.2 as isize
                        {
                            continue;
                        }
                        let fi =
                            ((br as usize) * self.nb.1 + bc as usize) * self.nb.2 + bl as usize;
                        for &(vr, vc, vl) in &self.bins[fi] {
                            let d2 = (rf - vr).powi(2) + (cf - vc).powi(2) + (lf - vl).powi(2);
                            if d2 < best {
                                best = d2;
                            }
                        }
                    }
                }
            }
            // After fully searching shell `shell`, any unsearched vessel sits in
            // a bin at Chebyshev ≥ shell+1, hence at lattice distance ≥
            // shell*bin. If the best found is already within that, no closer
            // vessel can exist — stop. (Conservative ⇒ exact.)
            if best.is_finite() && best <= (shell as f64 * self.bin).powi(2) {
                break;
            }
            shell += 1;
        }
        best
    }
}

/// Per-cell supply factor from the explicit vessel network: `exp(-d/λ)` where
/// `d` is the distance (µm) to the **nearest** vessel. Drop-in replacement for
/// [`crate::oxygen::radial_o2_field`]: returns a `Vec<f64>` of length
/// `grid.cells.len()`, non-tumor cells = `1.0` (well-perfused bulk), tumor
/// cells clamped to `[0, 1]`. Supplies both O2 (× `basal_ros`) and drug.
///
/// Distances are computed in lattice units and scaled by `grid.cell_size_um`.
///
/// **Cost**: nearest-vessel is resolved through a [`VesselIndex`] uniform-grid
/// spatial index (#268) — roughly `O(cells)` with a bounded per-query
/// neighborhood, **exact** (byte-identical to the former brute-force min). This
/// replaces the old `O(tumor_cells × vessels)` scan, which was quadratic-ish at
/// scale (a well-vascularized 200³ ≈ 34B evals, minutes) and was the blocker for
/// vasculature on patient-scale (#240) grids — worst of all on an all-tumor
/// **slab** (#272: every cell is tumor AND the box-volume vessel count is higher).
///
/// # Panics
/// If `vessels` is empty (no source ⇒ undefined supply); callers pass the
/// output of [`place_vessels_3d`], which is floored at 1 vessel.
pub fn vessel_supply_field(
    grid: &TumorGrid3D,
    vessels: &[(f64, f64, f64)],
    lambda_um: f64,
) -> Vec<f64> {
    assert!(!vessels.is_empty(), "vessel_supply_field needs ≥1 vessel");
    debug_assert!(
        lambda_um.is_finite() && lambda_um > 0.0,
        "vessel_supply_field: lambda_um must be finite and positive, got {lambda_um}"
    );
    let cell_size = grid.cell_size_um;
    let index = VesselIndex::build(grid, vessels);
    (0..grid.cells.len())
        .map(|idx| {
            if !grid.cells[idx].is_tumor {
                return 1.0;
            }
            let (r, c, l) = grid.coords(idx);
            let nearest_d2 = index.nearest_d2(r, c, l);
            let dist_um = nearest_d2.sqrt() * cell_size;
            (-dist_um / lambda_um).exp().clamp(0.0, 1.0)
        })
        .collect()
}

/// Fraction of tumor cells whose supply factor is below `threshold` (the
/// hypoxic fraction). Used to compare the explicit-vessel field against the
/// edge-distance proxy (#191 AC: irregular vasculature shifts the hypoxic
/// fraction relative to a smooth radial gradient).
pub fn hypoxic_fraction(grid: &TumorGrid3D, field: &[f64], threshold: f64) -> f64 {
    let (hyp, tot) = grid
        .cells
        .iter()
        .zip(field)
        .fold((0usize, 0usize), |(h, t), (gc, &f)| {
            if gc.is_tumor {
                (h + usize::from(f < threshold), t + 1)
            } else {
                (h, t)
            }
        });
    if tot > 0 {
        hyp as f64 / tot as f64
    } else {
        0.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::oxygen::radial_o2_field;
    use std::time::Instant;

    fn grid() -> TumorGrid3D {
        TumorGrid3D::generate(40, 40, 40, 20.0, 42)
    }

    /// Perf smoke (ignored — wall-clock, environment-dependent): the indexed
    /// `vessel_supply_field` on a large well-vascularized slab should finish in
    /// well under a second, vs the brute force which is `cells × vessels` evals
    /// (≈ minutes at this scale). Run with `cargo test -- --ignored --nocapture`.
    #[test]
    #[ignore]
    fn spatial_index_is_fast_at_scale() {
        let dim = 100;
        let g = TumorGrid3D::generate_slab(dim, dim, dim, 20.0, 42);
        let vessels = place_vessels_in_slab_3d(&g, &VasculatureConfig::well_vascularized(), 7);
        let t = Instant::now();
        let field = vessel_supply_field(&g, &vessels, 120.0);
        let ms = t.elapsed().as_secs_f64() * 1e3;
        let brute_evals = (g.cells.len() as u128) * (vessels.len() as u128);
        eprintln!(
            "indexed vessel_supply_field: {dim}³ all-tumor, {} vessels, {} cells → {ms:.1} ms \
             (brute force would be {brute_evals} distance evals)",
            vessels.len(),
            g.cells.len()
        );
        assert_eq!(field.len(), g.cells.len());
    }

    #[test]
    fn placement_is_deterministic_and_nonempty() {
        let g = grid();
        let a = place_vessels_3d(&g, &VasculatureConfig::well_vascularized(), 7);
        let b = place_vessels_3d(&g, &VasculatureConfig::well_vascularized(), 7);
        assert_eq!(a, b);
        assert!(!a.is_empty());
        // Denser config places more vessels than a sparse one.
        let sparse = place_vessels_3d(&g, &VasculatureConfig::poorly_vascularized(), 7);
        assert!(
            a.len() > sparse.len(),
            "well={}, poor={}",
            a.len(),
            sparse.len()
        );
    }

    #[test]
    fn fractal_placement_is_deterministic_and_gappier_than_random() {
        let g = grid();
        let cfg = VasculatureConfig::well_vascularized().with_fractal();
        let a = place_vessels_fractal_3d(&g, &cfg, 7);
        let b = place_vessels_fractal_3d(&g, &cfg, 7);
        assert_eq!(a, b, "fractal placement must be deterministic");
        assert!(!a.is_empty());
        let target = derive_vessel_count(&g, &cfg);
        let random = place_vessels_3d(&g, &VasculatureConfig::well_vascularized(), 7);
        // COUNT-PARITY GUARD: the fractal tree must fill to essentially the same
        // point count as the random placer (which returns exactly `target`).
        // This is the load-bearing invariant for the "gappier" claim below — if
        // the fractal network were merely SPARSER, a higher hypoxic fraction
        // would be a count artifact, not a topology (clustering) effect. A loose
        // bound (e.g. 0.5×) would let that confound slip through silently.
        assert!(
            a.len() <= target + 1 && a.len() as f64 >= 0.9 * random.len() as f64,
            "fractal count {} must be within 10% of random count {} (target {}) \
             so the hypoxic comparison stays count-controlled",
            a.len(),
            random.len(),
            target
        );
        // The KEY structural property: at (near-)equal point count, a fractal
        // network clusters its points along a few branches, leaving avascular
        // GAPS, so it leaves MORE tumor cells hypoxic than uniform-random
        // placement spread through the volume (perfusion holes from clustering).
        let lambda = 100.0;
        let hyp_fractal = hypoxic_fraction(&g, &vessel_supply_field(&g, &a, lambda), 0.2);
        let hyp_random = hypoxic_fraction(&g, &vessel_supply_field(&g, &random, lambda), 0.2);
        assert!(
            hyp_fractal > hyp_random,
            "fractal network should leave more avascular gaps than random: \
             fractal hypoxic={hyp_fractal:.3} (n={}), random hypoxic={hyp_random:.3} (n={})",
            a.len(),
            random.len()
        );
    }

    #[test]
    fn slab_placement_is_deterministic_and_fills_the_box() {
        // #272: slab-uniform placement spreads vessels across the WHOLE box,
        // unlike the sphere placement which confines them to the central
        // ~0.4-radius sphere. Same count (same volume/spacing rule), but some
        // vessels land outside that sphere — the property the coupling needs so
        // deep tissue throughout the slab gets perfused, not just the center.
        // Uses an all-tumor SLAB grid (the function's precondition); the central
        // sphere it's compared against is `place_vessels_3d`'s confinement.
        let g = TumorGrid3D::generate_slab(40, 40, 40, 20.0, 42);
        let cfg = VasculatureConfig::well_vascularized();
        let a = place_vessels_in_slab_3d(&g, &cfg, 7);
        let b = place_vessels_in_slab_3d(&g, &cfg, 7);
        assert_eq!(a, b, "slab placement must be deterministic");
        // Density agrees with the sphere placement (shared count rule).
        assert_eq!(a.len(), place_vessels_3d(&g, &cfg, 7).len());

        let center = (
            g.rows as f64 / 2.0,
            g.cols as f64 / 2.0,
            g.layers as f64 / 2.0,
        );
        let sphere_r = (g.rows.min(g.cols).min(g.layers) as f64) * TUMOR_RADIUS_FRACTION * 0.95;
        let outside = a
            .iter()
            .filter(|&&(r, c, l)| {
                let d = ((r - center.0).powi(2) + (c - center.1).powi(2) + (l - center.2).powi(2))
                    .sqrt();
                d > sphere_r
            })
            .count();
        assert!(
            outside > 0,
            "slab-uniform placement should reach beyond the central sphere (r={sphere_r:.1}); \
             {outside}/{} vessels outside",
            a.len()
        );
        // All vessels stay within the box.
        for &(r, c, l) in &a {
            assert!(r >= 0.0 && r <= g.rows as f64, "row {r} out of box");
            assert!(c >= 0.0 && c <= g.cols as f64, "col {c} out of box");
            assert!(l >= 0.0 && l <= g.layers as f64, "layer {l} out of box");
        }
    }

    #[test]
    fn spatial_index_matches_brute_force_bit_for_bit() {
        // #268: the uniform-grid index is an EXACT acceleration, so
        // `vessel_supply_field` stays byte-identical to the former brute force.
        // Recompute the field by brute force and assert bit-for-bit equality
        // (`to_bits`) across configs spanning bin sizes/occupancy and both vessel
        // placements (sphere and slab-uniform-in-box).
        let brute = |g: &TumorGrid3D, vessels: &[(f64, f64, f64)], lambda: f64| -> Vec<f64> {
            let cell_size = g.cell_size_um;
            (0..g.cells.len())
                .map(|idx| {
                    if !g.cells[idx].is_tumor {
                        return 1.0;
                    }
                    let (r, c, l) = g.coords(idx);
                    let (rf, cf, lf) = (r as f64, c as f64, l as f64);
                    let mut nd2 = f64::INFINITY;
                    for &(vr, vc, vl) in vessels {
                        let d2 = (rf - vr).powi(2) + (cf - vc).powi(2) + (lf - vl).powi(2);
                        if d2 < nd2 {
                            nd2 = d2;
                        }
                    }
                    (-(nd2.sqrt() * cell_size) / lambda).exp().clamp(0.0, 1.0)
                })
                .collect()
        };
        // (dims, inter_vessel_um, seed, slab?)
        let cases = [
            (30usize, 150.0, 1u64, false),
            (40, 300.0, 2, false),
            (24, 100.0, 3, false),
            (32, 200.0, 4, true),
        ];
        for (dims, spacing, seed, slab) in cases {
            let cfg = VasculatureConfig {
                inter_vessel_um: spacing,
                topology: VesselTopology::Random,
            };
            let (g, vessels) = if slab {
                let g = TumorGrid3D::generate_slab(dims, dims, dims, 20.0, seed);
                let v = place_vessels_in_slab_3d(&g, &cfg, seed);
                (g, v)
            } else {
                let g = TumorGrid3D::generate(dims, dims, dims, 20.0, seed);
                let v = place_vessels_3d(&g, &cfg, seed);
                (g, v)
            };
            let lambda = 120.0;
            let indexed = vessel_supply_field(&g, &vessels, lambda);
            let reference = brute(&g, &vessels, lambda);
            assert_eq!(indexed.len(), reference.len());
            for (i, (a, b)) in indexed.iter().zip(&reference).enumerate() {
                assert_eq!(
                    a.to_bits(),
                    b.to_bits(),
                    "indexed≠brute at cell {i} (dims={dims}, spacing={spacing}, slab={slab}, \
                     vessels={}): {a} vs {b}",
                    vessels.len()
                );
            }
        }
    }

    #[test]
    fn supply_is_one_for_stroma_and_in_range_for_tumor() {
        let g = grid();
        let v = place_vessels_3d(&g, &VasculatureConfig::well_vascularized(), 7);
        let field = vessel_supply_field(&g, &v, 100.0);
        assert_eq!(field.len(), g.cells.len());
        for (idx, &f) in field.iter().enumerate() {
            if g.cells[idx].is_tumor {
                assert!((0.0..=1.0).contains(&f), "tumor supply {f} out of [0,1]");
            } else {
                assert_eq!(f, 1.0, "stroma must be well-perfused (1.0)");
            }
        }
    }

    #[test]
    fn vessel_field_oxygenates_the_core_unlike_the_edge_proxy() {
        // The key #191 difference: the edge-distance proxy makes ALL deep-core
        // cells hypoxic (supply decays monotonically from the surface), whereas
        // explicit internal vessels reach the core, so some deep cells are
        // well-supplied. So among deep-core tumor cells, the mean vessel supply
        // exceeds the mean edge-proxy supply at matched λ — the irregular,
        // non-radial oxygenation the model is meant to capture.
        let g = grid();
        let lambda = 100.0;
        let edge = radial_o2_field(&g, lambda);
        let v = place_vessels_3d(&g, &VasculatureConfig::well_vascularized(), 7);
        let vessel = vessel_supply_field(&g, &v, lambda);

        let tumor_radius_um =
            (g.rows.min(g.cols).min(g.layers) as f64) * TUMOR_RADIUS_FRACTION * g.cell_size_um;
        let deep_threshold = 0.5 * tumor_radius_um;
        // Mean supply over deep-core tumor cells (depth > half the radius).
        let deep_mean = |field: &[f64]| -> f64 {
            let (sum, n) = (0..g.cells.len()).fold((0.0_f64, 0usize), |(s, n), idx| {
                let (r, c, l) = g.coords(idx);
                if g.cells[idx].is_tumor && g.radial_depth_um(r, c, l) > deep_threshold {
                    (s + field[idx], n + 1)
                } else {
                    (s, n)
                }
            });
            if n > 0 {
                sum / n as f64
            } else {
                0.0
            }
        };
        let (vessel_core, edge_core) = (deep_mean(&vessel), deep_mean(&edge));
        assert!(
            vessel_core > edge_core,
            "explicit vessels must oxygenate the deep core more than the edge proxy: \
             vessel_core_mean={vessel_core:.4}, edge_core_mean={edge_core:.4}"
        );
    }

    #[test]
    #[should_panic(expected = "needs ≥1 vessel")]
    fn empty_vessels_panics() {
        let g = grid();
        let _ = vessel_supply_field(&g, &[], 100.0);
    }
}
