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
//! | [`immune`] | ICD/DAMP immune cascade model (dimensionless, single-event) |
//! | [`immune_spatial`] | 3D spatial DAMP diffusion + per-cell immune activation/kill primitives |
//! | [`io`] | JSON and CSV output helpers |
//! | [`drug_transport`] | Tissue-specific drug penetration (Krogh cylinder approximation) |
//! | [`tumor_pk`] | Two-compartment vascular/interstitial pharmacokinetics |
//! | [`dose_schedule`] | Time-varying drug-administration schedules (bolus / multi-dose / infusion / PK-driven) |
//! | [`persister`] | Drug-tolerant persister cells (epigenetic ferroptosis tolerance, acquire/revert) |
//! | [`clonal`] | Voronoi subclone patches with per-subclone iron/GPX4/MUFA perturbations |
//! | [`vasculature`] | Explicit 3D vessel network + per-cell distance-decayed O2/drug supply |
//! | [`spheroid`] | 3D spheroid radial biology (rim/mid/core phenotypes + GSH/iron/MUFA gradients) |
//! | [`slab`] | Patient-scale slab geometry: all-tumor block + planar depth-graded O2/drug supply |
//! | [`contact`] | Cell-cell contact-mediated ferroptosis resistance (E-cadherin/NF2-YAP; dense cells lower PUFA/iron) |
//! | [`nutrient`] | Radial nutrient gradient (glucose/glutamine) scaling the antioxidant setpoint toward the starved core |
//! | [`reaction_diffusion`] | Steady-state reaction-diffusion O2/drug supply (vessel sources + consumption) vs the exponential proxy |

pub mod cell;
pub mod clonal;
pub mod contact;
// Listed before `params` because `SpatialParams` holds a `Photosensitizer`.
pub mod biochem;
pub mod dose_schedule;
pub mod drug_transport;
pub mod grid;
pub mod ifngamma;
pub mod immune;
pub mod immune_spatial;
pub mod io;
pub mod nutrient;
pub mod oxygen;
pub mod params;
pub mod persister;
pub mod ph;
pub mod phenotype_mufa;
pub mod photosensitizer_pk;
pub mod physics;
pub mod reaction_diffusion;
pub mod senescence;
pub mod slab;
pub mod spheroid;
pub mod stats;
pub mod stromal;
pub mod tumor_pk;
pub mod vasculature;
