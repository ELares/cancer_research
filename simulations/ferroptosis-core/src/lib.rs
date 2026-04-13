//! # ferroptosis-core
//!
//! Embeddable ferroptosis biochemistry engine for cancer simulation.
//!
//! This library models the ferroptosis cell death pathway at single-cell resolution:
//! ROS generation, GSH depletion, GPX4/FSP1 repair, lipid peroxidation, and
//! death threshold crossing. It supports both full single-cell simulations and
//! single-timestep updates for embedding in spatial or multi-scale frameworks.
//!
//! ## Key entry points
//!
//! - [`biochem::sim_cell`] — full 180-step ferroptosis simulation for one cell
//! - [`biochem::sim_cell_step`] — single timestep (for spatial model interleaving)
//! - [`cell::gen_cell`] — generate a cell with stochastic phenotype-specific parameters
//! - [`params::Params`] — all biochemistry rate constants (`default()` for 2D, `invivo()` for 3D)
//!
//! ## Modules
//!
//! | Module | Purpose |
//! |--------|---------|
//! | [`cell`] | Phenotypes, treatments, stochastic cell generation |
//! | [`params`] | Rate constants for biochemistry, physics, immune cascade |
//! | [`biochem`] | Core simulation engine |
//! | [`stats`] | Wilson CIs, parallel Monte Carlo execution |
//! | [`physics`] | Depth-dependent energy deposition (Beer-Lambert, acoustic) |
//! | [`grid`] | 2D tumor grid with heterogeneous architecture |
//! | [`immune`] | ICD/DAMP immune cascade model |
//! | [`io`] | JSON and CSV output helpers |
//! | [`drug_transport`] | Tissue-specific drug penetration (Krogh cylinder approximation) |

pub mod cell;
pub mod params;
pub mod biochem;
pub mod stats;
pub mod physics;
pub mod grid;
pub mod immune;
pub mod io;
pub mod drug_transport;
