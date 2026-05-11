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
//! | [`photosensitizer_pk`] | Photosensitizer plasma PK and drug-light-interval scaling for PDT |
//! | [`params`] | Rate constants for biochemistry, physics, immune cascade |
//! | [`biochem`] | Core simulation engine |
//! | [`stats`] | Wilson CIs, parallel Monte Carlo execution |
//! | [`physics`] | Depth-dependent energy deposition (Beer-Lambert, acoustic; 2D + 3D dispatchers) |
//! | [`grid`] | 2D and 3D tumor grids with heterogeneous architecture |
//! | [`oxygen`] | 3D radial oxygen gradients for spheroid tumors |
//! | [`ph`] | 3D radial pH gradient + iron-release and ion-trapping modulation helpers |
//! | [`stromal`] | 3D CAF-shielded boundary detection (26-Moore) + shielded kill rate |
//! | [`immune`] | ICD/DAMP immune cascade model |
//! | [`io`] | JSON and CSV output helpers |
//! | [`drug_transport`] | Tissue-specific drug penetration (Krogh cylinder approximation) |
//! | [`tumor_pk`] | Two-compartment vascular/interstitial pharmacokinetics |

pub mod cell;
// Listed before `params` because `SpatialParams` holds a `Photosensitizer`.
pub mod photosensitizer_pk;
pub mod params;
pub mod biochem;
pub mod stats;
pub mod physics;
pub mod grid;
pub mod oxygen;
pub mod ph;
pub mod stromal;
pub mod immune;
pub mod io;
pub mod drug_transport;
pub mod tumor_pk;
