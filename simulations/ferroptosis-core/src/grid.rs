//! Tumor grids for the spatial simulation.
//!
//! Provides tumor architecture generation, neighbor iteration, and iron
//! diffusion. The 2D model ([`TumorGrid`], 8-Moore neighbors, circular
//! tumor) is the established default used by `sim-spatial` and `sim-tme`.
//! The 3D model ([`TumorGrid3D`], 26-Moore neighbors, spherical tumor)
//! was added as foundational infrastructure for the spheroid-validation
//! series (#185–#197); analytics (`depth_kill_curve`, `death_heatmap`
//! analogs) land with the binary that actually uses them (#194).

use ndarray::Array2;
use rand::prelude::*;
use serde::Serialize;

use crate::cell::{gen_cell, Cell, Phenotype};
use crate::biochem::CellState;

/// A single cell in the spatial grid.
#[derive(Clone, Debug, Serialize)]
pub struct GridCell {
    pub cell: Cell,
    pub phenotype: Phenotype,
    pub state: CellState,
    pub is_tumor: bool,
    /// Extra iron from neighbor deaths, accumulated between timesteps.
    pub extra_iron: f64,
    /// LP at death (for DAMP calculation).
    pub lp_at_death: f64,
    /// Whether this cell just died this step (for diffusion).
    pub newly_dead: bool,
}

/// 2D tumor grid.
pub struct TumorGrid {
    pub cells: Vec<GridCell>,
    pub rows: usize,
    pub cols: usize,
    pub cell_size_um: f64,
}

impl TumorGrid {
    /// Generate a heterogeneous tumor grid.
    ///
    /// Layout:
    /// - Circular tumor centered in grid
    /// - Core (inner 60%): Glycolytic 80%, OXPHOS 15%, Persister 5%
    /// - Periphery (outer 40%): OXPHOS 70%, Glycolytic 20%, Persister 8%, PersisterNrf2 2%
    /// - 10-20 persister cluster seeds scattered throughout
    /// - Border: Stromal cells
    pub fn generate(rows: usize, cols: usize, cell_size_um: f64, seed: u64) -> Self {
        let mut rng = StdRng::seed_from_u64(seed);
        let center_r = rows as f64 / 2.0;
        let center_c = cols as f64 / 2.0;
        let tumor_radius = (rows.min(cols) as f64) * 0.45;
        let core_radius = tumor_radius * 0.6;

        // Generate persister cluster centers
        let n_clusters = 10 + (rng.gen_range(0..11) as usize);
        let cluster_centers: Vec<(f64, f64)> = (0..n_clusters)
            .map(|_| {
                let angle = rng.gen::<f64>() * std::f64::consts::TAU;
                let dist = rng.gen::<f64>().sqrt() * tumor_radius * 0.9;
                (center_r + angle.sin() * dist, center_c + angle.cos() * dist)
            })
            .collect();
        let cluster_radius = 4.0; // cells

        let mut cells = Vec::with_capacity(rows * cols);

        for r in 0..rows {
            for c in 0..cols {
                let dr = r as f64 - center_r;
                let dc = c as f64 - center_c;
                let dist = (dr * dr + dc * dc).sqrt();

                // Check if in a persister cluster
                let in_cluster = cluster_centers.iter().any(|(cr, cc)| {
                    let d = ((r as f64 - cr).powi(2) + (c as f64 - cc).powi(2)).sqrt();
                    d <= cluster_radius
                });

                let (phenotype, is_tumor) = if dist > tumor_radius {
                    (Phenotype::Stromal, false)
                } else if in_cluster {
                    // Persister clusters are pure persister
                    (Phenotype::Persister, true)
                } else if dist <= core_radius {
                    // Core: mostly glycolytic
                    let roll: f64 = rng.gen();
                    if roll < 0.80 {
                        (Phenotype::Glycolytic, true)
                    } else if roll < 0.95 {
                        (Phenotype::OXPHOS, true)
                    } else {
                        (Phenotype::Persister, true)
                    }
                } else {
                    // Periphery: mostly OXPHOS
                    let roll: f64 = rng.gen();
                    if roll < 0.70 {
                        (Phenotype::OXPHOS, true)
                    } else if roll < 0.90 {
                        (Phenotype::Glycolytic, true)
                    } else if roll < 0.98 {
                        (Phenotype::Persister, true)
                    } else {
                        (Phenotype::PersisterNrf2, true)
                    }
                };

                let cell = gen_cell(phenotype, &mut rng);
                // State will be initialized later when treatment is applied
                let state = CellState {
                    gsh: cell.gsh,
                    gpx4: cell.gpx4,
                    fsp1: cell.fsp1,
                    mufa_protection: 0.0,
                    lp: 0.0,
                    dead: false,
                    death_step: None,
                    exo_ros_peak: 0.0,
                };

                cells.push(GridCell {
                    cell,
                    phenotype,
                    state,
                    is_tumor,
                    extra_iron: 0.0,
                    lp_at_death: 0.0,
                    newly_dead: false,
                });
            }
        }

        TumorGrid {
            cells,
            rows,
            cols,
            cell_size_um,
        }
    }

    /// Access cell at (row, col).
    #[inline]
    pub fn get(&self, r: usize, c: usize) -> &GridCell {
        &self.cells[r * self.cols + c]
    }

    /// Mutable access to cell at (row, col).
    #[inline]
    pub fn get_mut(&mut self, r: usize, c: usize) -> &mut GridCell {
        &mut self.cells[r * self.cols + c]
    }

    /// Return indices of Moore neighborhood (8-neighbors) for cell (r, c).
    /// Respects grid boundaries (no wrapping).
    /// Returns (array, count) — use `&result[..count]` to iterate.
    /// Zero-allocation: uses a stack-allocated fixed-size array.
    pub fn neighbors(&self, r: usize, c: usize) -> ([(usize, usize); 8], usize) {
        let mut result = [(0usize, 0usize); 8];
        let mut count = 0;
        for dr in [-1_i64, 0, 1] {
            for dc in [-1_i64, 0, 1] {
                if dr == 0 && dc == 0 {
                    continue;
                }
                let nr = r as i64 + dr;
                let nc = c as i64 + dc;
                if nr >= 0 && nr < self.rows as i64 && nc >= 0 && nc < self.cols as i64 {
                    result[count] = (nr as usize, nc as usize);
                    count += 1;
                }
            }
        }
        (result, count)
    }

    /// Distribute iron from newly dead cells to their living neighbors.
    /// Each living neighbor receives `neighbor_iron_fraction` of the released iron.
    pub fn diffuse_iron(&mut self, iron_per_death: f64, neighbor_fraction: f64) {
        // Collect positions of newly dead cells
        let dead_positions: Vec<(usize, usize)> = (0..self.rows)
            .flat_map(|r| (0..self.cols).map(move |c| (r, c)))
            .filter(|&(r, c)| self.get(r, c).newly_dead)
            .collect();

        // Distribute iron to living neighbors
        for (r, c) in dead_positions {
            let (neighbors, count) = self.neighbors(r, c);
            for &(nr, nc) in &neighbors[..count] {
                let neighbor = self.get_mut(nr, nc);
                if !neighbor.state.dead {
                    neighbor.extra_iron += iron_per_death * neighbor_fraction;
                }
            }
            // Clear newly_dead flag
            self.get_mut(r, c).newly_dead = false;
        }
    }

    /// Count cells by phenotype and alive/dead status.
    pub fn census(&self) -> GridCensus {
        let mut census = GridCensus::default();
        for gc in &self.cells {
            if !gc.is_tumor {
                census.stromal += 1;
                continue;
            }
            census.total_tumor += 1;
            if gc.state.dead {
                census.total_dead += 1;
                match gc.phenotype {
                    Phenotype::Glycolytic => census.glycolytic_dead += 1,
                    Phenotype::OXPHOS => census.oxphos_dead += 1,
                    Phenotype::Persister => census.persister_dead += 1,
                    Phenotype::PersisterNrf2 => census.persister_nrf2_dead += 1,
                    Phenotype::Stromal => {}
                }
            }
            match gc.phenotype {
                Phenotype::Glycolytic => census.glycolytic += 1,
                Phenotype::OXPHOS => census.oxphos += 1,
                Phenotype::Persister => census.persister += 1,
                Phenotype::PersisterNrf2 => census.persister_nrf2 += 1,
                Phenotype::Stromal => {}
            }
        }
        census
    }
}

/// Summary statistics of grid composition.
#[derive(Default, Debug, serde::Serialize)]
pub struct GridCensus {
    pub total_tumor: usize,
    pub total_dead: usize,
    pub stromal: usize,
    pub glycolytic: usize,
    pub glycolytic_dead: usize,
    pub oxphos: usize,
    pub oxphos_dead: usize,
    pub persister: usize,
    pub persister_dead: usize,
    pub persister_nrf2: usize,
    pub persister_nrf2_dead: usize,
}

/// Compute death rate by depth (row).
/// Returns Vec of (depth_um, death_rate, total_cells) for each row.
pub fn depth_kill_curve(grid: &TumorGrid) -> Vec<(f64, f64, usize)> {
    (0..grid.rows)
        .map(|r| {
            let depth_um = r as f64 * grid.cell_size_um;
            let mut total = 0usize;
            let mut dead = 0usize;
            for c in 0..grid.cols {
                let gc = grid.get(r, c);
                if gc.is_tumor {
                    total += 1;
                    if gc.state.dead {
                        dead += 1;
                    }
                }
            }
            let rate = if total > 0 {
                dead as f64 / total as f64
            } else {
                0.0
            };
            (depth_um, rate, total)
        })
        .collect()
}

/// Export death map as a 2D array (0 = alive/stromal, 1 = dead tumor, 2 = alive tumor).
pub fn death_heatmap(grid: &TumorGrid) -> Array2<u8> {
    Array2::from_shape_fn((grid.rows, grid.cols), |(r, c)| {
        let gc = grid.get(r, c);
        if !gc.is_tumor {
            0
        } else if gc.state.dead {
            1
        } else {
            2
        }
    })
}

// =====================================================================
// 3D tumor grid (#185 foundational infrastructure for spheroid modeling)
// =====================================================================

/// 3D tumor grid — additive analog of [`TumorGrid`] with the same
/// `GridCell` element type, spherical tumor geometry, and 26-Moore
/// neighbors.
///
/// **Design choice:** This is a new struct rather than extending
/// `TumorGrid` with an optional z-dimension. The decision is documented
/// in #185 — 2D physics, 2D test snapshots, and all 2D consumers
/// (sim-spatial, sim-tme) keep their existing `TumorGrid` semantics
/// bit-for-bit. The two types share `GridCell` and `GridCensus`; any
/// common abstraction can land as a trait when patterns emerge.
///
/// **Memory:** `GridCell` is ~150–170 B; 100³ ≈ 170 MB, 200³ ≈ 1.4 GB.
/// Tests use ≤ 20³ (≈ 1.4 MB). Realistic 3D simulation grids land with
/// the binary that constructs them (#194 `sim-spatial-3d`).
///
/// **Storage:** flat `Vec<GridCell>` indexed as
/// `r * cols * layers + c * layers + l` — row-major (C-order) with
/// strides `(cols*layers, layers, 1)`. The `l` axis has stride 1, so
/// `for r { for c { for l { ... } } }` is cache-friendly. Extends
/// `TumorGrid`'s row-major `r * cols + c` to a third dimension.
pub struct TumorGrid3D {
    pub cells: Vec<GridCell>,
    pub rows: usize,
    pub cols: usize,
    pub layers: usize,
    pub cell_size_um: f64,
}

impl TumorGrid3D {
    /// Generate a heterogeneous spherical tumor grid (3D analog of
    /// [`TumorGrid::generate`]).
    ///
    /// Layout:
    /// - Spherical tumor centered in `(rows/2, cols/2, layers/2)` with
    ///   radius = `min(rows, cols, layers) * 0.45`
    /// - Core (inner 60% of tumor radius): Glycolytic 80%, OXPHOS 15%, Persister 5%
    /// - Periphery (outer 40%): OXPHOS 70%, Glycolytic 20%, Persister 8%, PersisterNrf2 2%
    /// - 10–20 persister cluster seeds scattered uniformly by volume
    /// - Outside tumor sphere: Stromal cells
    ///
    /// Deterministic from `seed`.
    ///
    /// **Geometry:** cell positions are lattice-centered (`r` etc. are
    /// treated as continuous coords with no `+ 0.5` voxel offset), so
    /// the cell at `(rows/2, cols/2, layers/2)` sits exactly at the
    /// sphere center when dimensions are even. Same convention as
    /// `TumorGrid::generate`.
    ///
    /// **Memory:** `GridCell` is ~150–170 B; a `100³` grid is ~170 MB,
    /// `200³` is ~1.4 GB. Callers should size against available RAM
    /// before invoking; `debug_assert!` guards against accidental
    /// allocations above ~1 G cells (≈ 170 GB, almost certainly a typo).
    pub fn generate(
        rows: usize,
        cols: usize,
        layers: usize,
        cell_size_um: f64,
        seed: u64,
    ) -> Self {
        // Memory guard — catches "rows = 1000" typos in tests and CI.
        // Active only in debug builds; release callers are expected to
        // size against their own hardware. 10^9 cells × 170 B ≈ 170 GB
        // is well past any realistic spheroid sim and almost certainly
        // an off-by-orders-of-magnitude error.
        debug_assert!(
            rows.saturating_mul(cols).saturating_mul(layers) <= 1_000_000_000,
            "TumorGrid3D::generate: rows*cols*layers = {} exceeds 1e9 — likely a typo (use saturating mul to avoid wraparound)",
            rows.saturating_mul(cols).saturating_mul(layers)
        );
        let mut rng = StdRng::seed_from_u64(seed);
        let center_r = rows as f64 / 2.0;
        let center_c = cols as f64 / 2.0;
        let center_l = layers as f64 / 2.0;
        let tumor_radius = (rows.min(cols).min(layers) as f64) * 0.45;
        let core_radius = tumor_radius * 0.6;

        // Generate persister cluster centers uniformly distributed by
        // *volume* inside the tumor. Spherical-coordinate sampling:
        //   θ ~ U(0, 2π)              (azimuth)
        //   cos(φ) ~ U(-1, 1)         (polar; uniform-area on sphere)
        //   r = cbrt(U(0, 1)) * R     (uniform-volume in ball)
        // Without the cbrt, samples bias toward the center (analogous to
        // 2D where sqrt(rand) gives uniform-area; cbrt(rand) is 3D's
        // uniform-volume).
        //
        // Axis convention: `r` is the polar (vertical) axis, `(c, l)` form
        // the azimuthal plane. So a cluster at φ = 0 sits on the r axis
        // (highest row offset), and θ rotates within the (c, l) plane.
        // Choice is arbitrary; documented here because the asymmetric
        // r vs c/l usage is non-obvious from the math.
        let n_clusters = 10 + (rng.gen_range(0..11) as usize);
        let cluster_centers: Vec<(f64, f64, f64)> = (0..n_clusters)
            .map(|_| {
                let theta = rng.gen::<f64>() * std::f64::consts::TAU;
                let cos_phi = 1.0 - 2.0 * rng.gen::<f64>();
                let sin_phi = (1.0 - cos_phi * cos_phi).sqrt();
                let dist = rng.gen::<f64>().cbrt() * tumor_radius * 0.9;
                (
                    center_r + dist * cos_phi,
                    center_c + dist * sin_phi * theta.cos(),
                    center_l + dist * sin_phi * theta.sin(),
                )
            })
            .collect();
        // TODO(#194): cluster_radius=4.0 cells gives ~268 cells per
        // cluster in 3D (4/3·π·4³) vs ~50 in 2D (π·4²) — so persister
        // clusters are ~5× larger relative to tumor volume than the
        // 2D analog. A volume-equivalent retune (≈ cluster_radius =
        // 4·(50/268)^(1/3) ≈ 2.4 cells, giving ~50 cells/cluster) keeps
        // the 2D↔3D census proportions comparable. Left at 4.0 for v1
        // to match the 2D code shape literally; sim-spatial-3d (#194)
        // is the natural place to calibrate against experimental data.
        let cluster_radius = 4.0;

        let mut cells = Vec::with_capacity(rows * cols * layers);

        for r in 0..rows {
            for c in 0..cols {
                for l in 0..layers {
                    let dr = r as f64 - center_r;
                    let dc = c as f64 - center_c;
                    let dl = l as f64 - center_l;
                    let dist = (dr * dr + dc * dc + dl * dl).sqrt();

                    let in_cluster = cluster_centers.iter().any(|(cr, cc, cl)| {
                        let d = ((r as f64 - cr).powi(2)
                            + (c as f64 - cc).powi(2)
                            + (l as f64 - cl).powi(2))
                        .sqrt();
                        d <= cluster_radius
                    });

                    let (phenotype, is_tumor) = if dist > tumor_radius {
                        (Phenotype::Stromal, false)
                    } else if in_cluster {
                        (Phenotype::Persister, true)
                    } else if dist <= core_radius {
                        let roll: f64 = rng.gen();
                        if roll < 0.80 {
                            (Phenotype::Glycolytic, true)
                        } else if roll < 0.95 {
                            (Phenotype::OXPHOS, true)
                        } else {
                            (Phenotype::Persister, true)
                        }
                    } else {
                        let roll: f64 = rng.gen();
                        if roll < 0.70 {
                            (Phenotype::OXPHOS, true)
                        } else if roll < 0.90 {
                            (Phenotype::Glycolytic, true)
                        } else if roll < 0.98 {
                            (Phenotype::Persister, true)
                        } else {
                            (Phenotype::PersisterNrf2, true)
                        }
                    };

                    let cell = gen_cell(phenotype, &mut rng);
                    let state = CellState {
                        gsh: cell.gsh,
                        gpx4: cell.gpx4,
                        fsp1: cell.fsp1,
                        mufa_protection: 0.0,
                        lp: 0.0,
                        dead: false,
                        death_step: None,
                        exo_ros_peak: 0.0,
                    };
                    cells.push(GridCell {
                        cell,
                        phenotype,
                        state,
                        is_tumor,
                        extra_iron: 0.0,
                        lp_at_death: 0.0,
                        newly_dead: false,
                    });
                }
            }
        }

        TumorGrid3D {
            cells,
            rows,
            cols,
            layers,
            cell_size_um,
        }
    }

    /// Access cell at (row, col, layer).
    #[inline]
    pub fn get(&self, r: usize, c: usize, l: usize) -> &GridCell {
        &self.cells[r * self.cols * self.layers + c * self.layers + l]
    }

    /// Mutable access to cell at (row, col, layer).
    #[inline]
    pub fn get_mut(&mut self, r: usize, c: usize, l: usize) -> &mut GridCell {
        &mut self.cells[r * self.cols * self.layers + c * self.layers + l]
    }

    /// Return indices of 3D Moore neighborhood (26-neighbors) for cell
    /// (r, c, l). Respects grid boundaries (no wrapping). Zero-allocation:
    /// uses a stack-allocated fixed-size array. Use `&result[..count]`
    /// to iterate.
    ///
    /// Counts by position type:
    /// - Interior cell: 26 neighbors (3³ − 1)
    /// - Face cell (one coordinate at the boundary): 17 (2·3·3 − 1)
    /// - Face-edge cell (two coordinates at boundary): 11 (2·2·3 − 1)
    /// - Corner cell (all three at boundary): 7 (2³ − 1)
    pub fn neighbors(
        &self,
        r: usize,
        c: usize,
        l: usize,
    ) -> ([(usize, usize, usize); 26], usize) {
        let mut result = [(0usize, 0usize, 0usize); 26];
        let mut count = 0;
        for dr in [-1_i64, 0, 1] {
            for dc in [-1_i64, 0, 1] {
                for dl in [-1_i64, 0, 1] {
                    if dr == 0 && dc == 0 && dl == 0 {
                        continue;
                    }
                    let nr = r as i64 + dr;
                    let nc = c as i64 + dc;
                    let nl = l as i64 + dl;
                    if nr >= 0
                        && nr < self.rows as i64
                        && nc >= 0
                        && nc < self.cols as i64
                        && nl >= 0
                        && nl < self.layers as i64
                    {
                        result[count] = (nr as usize, nc as usize, nl as usize);
                        count += 1;
                    }
                }
            }
        }
        (result, count)
    }

    /// Distribute iron from newly dead cells to their living neighbors.
    /// Each living neighbor receives `neighbor_fraction` of the released
    /// iron.
    ///
    /// **Note for callers:** 2D `TumorGrid::diffuse_iron` uses
    /// `neighbor_fraction = 0.1` against 8 Moore neighbors, so up to 80%
    /// of released iron is distributed. In 3D with 26 Moore neighbors,
    /// the same `0.1` would distribute up to 260% — non-physical. The
    /// natural 3D analog is `0.1 * 8 / 26 ≈ 0.0308`, but the actual
    /// calibration is left to the caller (sim-spatial-3d, #194) so they
    /// can choose based on their experimental targets.
    pub fn diffuse_iron(&mut self, iron_per_death: f64, neighbor_fraction: f64) {
        // Hoist dimensions so the closures don't try to capture `&mut self`.
        let (rows, cols, layers) = (self.rows, self.cols, self.layers);
        let dead_positions: Vec<(usize, usize, usize)> = (0..rows)
            .flat_map(|r| (0..cols).flat_map(move |c| (0..layers).map(move |l| (r, c, l))))
            .filter(|&(r, c, l)| self.get(r, c, l).newly_dead)
            .collect();

        for (r, c, l) in dead_positions {
            let (neighbors, count) = self.neighbors(r, c, l);
            for &(nr, nc, nl) in &neighbors[..count] {
                let neighbor = self.get_mut(nr, nc, nl);
                if !neighbor.state.dead {
                    neighbor.extra_iron += iron_per_death * neighbor_fraction;
                }
            }
            self.get_mut(r, c, l).newly_dead = false;
        }
    }

    /// Count cells by phenotype and alive/dead status. Reuses
    /// [`GridCensus`] from the 2D model.
    pub fn census(&self) -> GridCensus {
        let mut census = GridCensus::default();
        for gc in &self.cells {
            if !gc.is_tumor {
                census.stromal += 1;
                continue;
            }
            census.total_tumor += 1;
            if gc.state.dead {
                census.total_dead += 1;
                match gc.phenotype {
                    Phenotype::Glycolytic => census.glycolytic_dead += 1,
                    Phenotype::OXPHOS => census.oxphos_dead += 1,
                    Phenotype::Persister => census.persister_dead += 1,
                    Phenotype::PersisterNrf2 => census.persister_nrf2_dead += 1,
                    Phenotype::Stromal => {}
                }
            }
            match gc.phenotype {
                Phenotype::Glycolytic => census.glycolytic += 1,
                Phenotype::OXPHOS => census.oxphos += 1,
                Phenotype::Persister => census.persister += 1,
                Phenotype::PersisterNrf2 => census.persister_nrf2 += 1,
                Phenotype::Stromal => {}
            }
        }
        census
    }
}

#[cfg(test)]
mod tests_3d {
    use super::*;

    /// Index round-trip at several positions including corners and the
    /// last cell. Catches off-by-one in `r * cols * layers + c * layers + l`.
    #[test]
    fn index_round_trip() {
        let g = TumorGrid3D::generate(5, 7, 11, 20.0, 42);
        // Probe diverse positions.
        let probes: &[(usize, usize, usize)] = &[
            (0, 0, 0),
            (4, 6, 10), // last cell
            (2, 3, 5),  // interior
            (0, 6, 0),  // edge
            (4, 0, 10), // opposite corner
        ];
        for &(r, c, l) in probes {
            let flat = r * g.cols * g.layers + c * g.layers + l;
            assert_eq!(
                &g.cells[flat] as *const _,
                g.get(r, c, l) as *const _,
                "(r={r}, c={c}, l={l}) → flat={flat} doesn't round-trip"
            );
        }
        // Total cell count matches.
        assert_eq!(g.cells.len(), 5 * 7 * 11);
    }

    /// 26-Moore neighbor counts by position type. The exact counts catch
    /// boundary off-by-ones in any direction.
    #[test]
    fn neighbor_counts_at_boundary_types() {
        let g = TumorGrid3D::generate(5, 5, 5, 20.0, 42);

        // Interior (all three coords strictly inside): 26
        let (_, n_interior) = g.neighbors(2, 2, 2);
        assert_eq!(n_interior, 26, "interior cell should have 26 neighbors");

        // Face (one coord on boundary, others interior): 17 = 2·3·3 − 1
        let (_, n_face) = g.neighbors(0, 2, 2);
        assert_eq!(n_face, 17, "face cell should have 17 neighbors");
        let (_, n_face_b) = g.neighbors(2, 4, 2);
        assert_eq!(n_face_b, 17, "opposite-face cell should also have 17");

        // Face-edge (two coords on boundary, one interior): 11 = 2·2·3 − 1
        let (_, n_edge) = g.neighbors(0, 0, 2);
        assert_eq!(n_edge, 11, "face-edge cell should have 11 neighbors");

        // Corner (all three on boundary): 7 = 2³ − 1
        let (_, n_corner) = g.neighbors(0, 0, 0);
        assert_eq!(n_corner, 7, "corner cell should have 7 neighbors");
        let (_, n_corner_b) = g.neighbors(4, 4, 4);
        assert_eq!(n_corner_b, 7, "opposite corner should also have 7");
    }

    /// Neighbor positions are always within bounds and never include the
    /// cell itself. Catches an off-by-one that would include or skip the
    /// center.
    #[test]
    fn neighbors_are_valid_and_exclude_self() {
        let g = TumorGrid3D::generate(5, 5, 5, 20.0, 42);
        let (neighbors, count) = g.neighbors(2, 2, 2);
        for &(nr, nc, nl) in &neighbors[..count] {
            assert!(nr < g.rows && nc < g.cols && nl < g.layers);
            assert!(
                !(nr == 2 && nc == 2 && nl == 2),
                "neighbors must not include self"
            );
        }
    }

    /// `diffuse_iron` distributes to all valid neighbors of a newly-dead
    /// interior cell. Picks a cell guaranteed to be tumor (center) so
    /// neighbors are tumor too.
    #[test]
    fn diffuse_iron_distributes_to_all_26_interior_neighbors() {
        let mut g = TumorGrid3D::generate(10, 10, 10, 20.0, 42);

        // Pick an interior cell. Sphere center may or may not be alive;
        // either way, mark a known interior position as newly-dead and
        // verify diffusion fans out to all 26 neighbors that aren't
        // already-dead.
        let (r, c, l) = (5, 5, 5);
        g.get_mut(r, c, l).newly_dead = true;
        // Ensure no neighbors are already-dead before the call.
        let (positions, count) = g.neighbors(r, c, l);
        for &(nr, nc, nl) in &positions[..count] {
            g.get_mut(nr, nc, nl).state.dead = false;
            g.get_mut(nr, nc, nl).extra_iron = 0.0;
        }
        assert_eq!(count, 26, "interior cell must have 26 neighbors");

        let iron_per_death = 2.0;
        let neighbor_fraction = 0.0308; // 3D-natural fraction per docstring
        let expected = iron_per_death * neighbor_fraction;
        g.diffuse_iron(iron_per_death, neighbor_fraction);

        for &(nr, nc, nl) in &positions[..count] {
            // Assert the *magnitude* — not just presence — so a bug that
            // applied neighbor_fraction twice or swapped iron_per_death
            // for a constant is caught. IEEE-exact: only one multiply, no
            // libm-dependent ops, so strict equality is safe.
            assert_eq!(
                g.get(nr, nc, nl).extra_iron,
                expected,
                "neighbor ({nr},{nc},{nl}) extra_iron should equal iron_per_death × neighbor_fraction"
            );
        }
        assert!(!g.get(r, c, l).newly_dead, "newly_dead flag should be cleared after diffusion");
    }

    /// Generate is deterministic from seed: same seed → same per-cell
    /// assignment (not just same aggregate census).
    ///
    /// The original version only compared `census()` counts, which would
    /// pass even if phenotype assignments were permuted across the cells
    /// vec (e.g. a bug that swapped the order of cluster checks but kept
    /// the same totals). This stronger version walks every cell and
    /// compares `(phenotype, is_tumor)` pairwise — a true bit-level
    /// determinism check on `generate`'s output.
    #[test]
    fn generate_is_deterministic() {
        let g1 = TumorGrid3D::generate(20, 20, 20, 20.0, 42);
        let g2 = TumorGrid3D::generate(20, 20, 20, 20.0, 42);

        assert_eq!(g1.cells.len(), g2.cells.len());
        for (i, (a, b)) in g1.cells.iter().zip(g2.cells.iter()).enumerate() {
            assert_eq!(
                a.phenotype, b.phenotype,
                "phenotype divergence at flat index {i}"
            );
            assert_eq!(
                a.is_tumor, b.is_tumor,
                "is_tumor divergence at flat index {i}"
            );
        }
    }

    /// Spherical tumor: census has both tumor and stromal cells, and
    /// the center of a sufficiently large grid is tumor while a corner
    /// is stromal.
    #[test]
    fn spherical_geometry_sanity() {
        let g = TumorGrid3D::generate(20, 20, 20, 20.0, 42);
        let census = g.census();
        assert!(census.total_tumor > 0, "should have some tumor cells");
        assert!(census.stromal > 0, "should have some stromal cells");
        // Total = volume of grid
        assert_eq!(census.total_tumor + census.stromal, 20 * 20 * 20);
        // Center cell is inside the tumor sphere (tumor_radius = 20*0.45 = 9.0)
        assert!(g.get(10, 10, 10).is_tumor, "center cell should be tumor");
        // Corner is well outside (distance from center ≈ 17.3 > 9.0)
        assert!(!g.get(0, 0, 0).is_tumor, "corner cell should be stromal");
    }

    /// `diffuse_iron` correctly handles a corner cell (count=7, not 26).
    /// Catches a hypothetical "iterate `&neighbors[..26]` unconditionally"
    /// bug — the interior test alone wouldn't notice that, since at
    /// interior cells count *is* 26. We override `is_tumor` on the corner
    /// and its 7 neighbors so the diffusion path isn't gated by phenotype.
    #[test]
    fn diffuse_iron_distributes_to_only_7_corner_neighbors() {
        let mut g = TumorGrid3D::generate(10, 10, 10, 20.0, 42);

        let (r, c, l) = (0, 0, 0); // grid corner
        // Force the corner to be a fresh tumor cell that just died.
        let gc = g.get_mut(r, c, l);
        gc.is_tumor = true;
        gc.state.dead = true;
        gc.newly_dead = true;
        gc.extra_iron = 0.0;

        // Record neighbor identity + force them to be alive tumor.
        let (positions, count) = g.neighbors(r, c, l);
        assert_eq!(count, 7, "corner should have 7 neighbors");
        for &(nr, nc, nl) in &positions[..count] {
            let nb = g.get_mut(nr, nc, nl);
            nb.is_tumor = true;
            nb.state.dead = false;
            nb.extra_iron = 0.0;
        }

        let iron_per_death = 2.0;
        let neighbor_fraction = 0.0308;
        let expected = iron_per_death * neighbor_fraction;
        g.diffuse_iron(iron_per_death, neighbor_fraction);

        // All 7 valid neighbors received exactly the right magnitude
        // (catches a "iterate ..26 unconditionally" bug — out-of-bounds
        // would panic — AND a "wrong multiplier" bug).
        for &(nr, nc, nl) in &positions[..count] {
            assert_eq!(
                g.get(nr, nc, nl).extra_iron,
                expected,
                "corner neighbor ({nr},{nc},{nl}) extra_iron should equal expected"
            );
        }
        assert!(!g.get(r, c, l).newly_dead, "newly_dead should be cleared");
    }

    /// `layers = 1` collapses TumorGrid3D to a single-slab grid: index
    /// math becomes `r * cols + c` (same as 2D), and neighbor counts
    /// fall back to 8-Moore in the (r,c) plane since dl ∈ {0} only.
    /// Locks down the useful "TumorGrid3D-as-2D-fallback" emergent
    /// property — anything that breaks this is a regression in the
    /// general index/neighbor formulas.
    #[test]
    fn layers_eq_1_degenerates_to_2d() {
        let g = TumorGrid3D::generate(5, 5, 1, 20.0, 42);

        // Index math: r * 5 * 1 + c * 1 + 0 == r * 5 + c
        for r in 0..5 {
            for c in 0..5 {
                let flat = r * g.cols * g.layers + c * g.layers;
                assert_eq!(flat, r * 5 + c, "layers=1 index should collapse to 2D");
            }
        }

        // Interior cell at (2, 2, 0): 8 neighbors (Moore in (r,c) plane;
        // dl ∈ {0} only since layers=1).
        let (_, n_interior) = g.neighbors(2, 2, 0);
        assert_eq!(
            n_interior, 8,
            "interior cell with layers=1 should have 8 Moore neighbors (2D fallback)"
        );

        // Corner at (0, 0, 0): 3 neighbors ({0,1}×{0,1} - self = 3).
        let (_, n_corner) = g.neighbors(0, 0, 0);
        assert_eq!(
            n_corner, 3,
            "corner with layers=1 should have 3 neighbors (2D-corner fallback)"
        );

        // Edge at (0, 2, 0): 5 neighbors ({0,1}×{1,2,3} - self = 5).
        let (_, n_edge) = g.neighbors(0, 2, 0);
        assert_eq!(
            n_edge, 5,
            "edge with layers=1 should have 5 neighbors (2D-edge fallback)"
        );
    }

    /// Census on a fresh grid has no dead cells.
    #[test]
    fn fresh_grid_has_no_dead_cells() {
        let g = TumorGrid3D::generate(10, 10, 10, 20.0, 42);
        let census = g.census();
        assert_eq!(census.total_dead, 0);
        assert_eq!(census.glycolytic_dead, 0);
        assert_eq!(census.oxphos_dead, 0);
        assert_eq!(census.persister_dead, 0);
        assert_eq!(census.persister_nrf2_dead, 0);
    }
}
