//! Dump the steady-state reaction-diffusion field for an external-solver
//! cross-check (#408). Builds an all-tumor cubic grid (uniform consumption, to
//! match a uniform-decay BioFVM scenario), runs `reaction_diffusion_solve` with
//! Dirichlet vessel sources at integer voxel coordinates, and prints the field as
//! `r,c,l,c` CSV (grid index order). Used by scripts/validate_rd_vs_biofvm.py.
//!
//! Usage: cargo run --release --example rd_field_dump -- N h lambda  r0 c0 l0 [r1 c1 l1 ...]

use ferroptosis_core::grid::TumorGrid3D;
use ferroptosis_core::reaction_diffusion::{reaction_diffusion_solve, ReactionDiffusionConfig};
use std::env;

fn main() {
    let a: Vec<String> = env::args().collect();
    let n: usize = a[1].parse().unwrap();
    let h: f64 = a[2].parse().unwrap();
    let lambda: f64 = a[3].parse().unwrap();
    let mut vessels: Vec<(f64, f64, f64)> = Vec::new();
    let mut i = 4;
    while i + 2 < a.len() {
        vessels.push((
            a[i].parse().unwrap(),
            a[i + 1].parse().unwrap(),
            a[i + 2].parse().unwrap(),
        ));
        i += 3;
    }

    let mut grid = TumorGrid3D::generate(n, n, n, h, 7);
    for cell in grid.cells.iter_mut() {
        cell.is_tumor = true; // uniform consumption everywhere (matches uniform-decay BioFVM)
    }
    let mut cfg = ReactionDiffusionConfig::new(lambda);
    cfg.max_iters = 200_000;
    cfg.tol = 1e-10;
    let sol = reaction_diffusion_solve(&grid, &vessels, &cfg);

    println!("r,c,l,val");
    for r in 0..n {
        for c in 0..n {
            for l in 0..n {
                let v = sol.field[grid.flat_index(r, c, l)];
                println!("{},{},{},{:.8}", r, c, l, v);
            }
        }
    }
    eprintln!(
        "rust_rd converged={} iters={} residual={:.3e}",
        sol.converged, sol.iters, sol.residual
    );
}
