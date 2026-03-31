//! 2D tumor grid for the spatial simulation.
//!
//! Provides tumor architecture generation, neighbor iteration, and iron diffusion.

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
