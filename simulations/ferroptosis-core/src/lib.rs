//! # ferroptosis-core
//!
//! Embeddable cancer-therapy simulation engine.
//!
//! The crate is NAMED for ferroptosis and began as a ferroptosis engine; it
//! now also carries radiation, ablation, oncolytic spread, ADC bystander
//! killing and adoptive-cell barriers. The ferroptosis layers remain the
//! deepest and carry the calibrated legs -- see
//! `analysis/modality-module-depth.md` for the measured gap.
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
//! | [`ablation`] | Physical ablation: HIFU thermal dose (CEM43) and irreversible electroporation, as THRESHOLD phenomena rather than dose-responses |
//! | [`adoptive`] | Adoptive cell therapy: the three sequential barriers -- trafficking, infiltration, activation -- that multiply into the solid-tumour failure, plus persistence |
//! | [`adc`] | Antibody-drug conjugates: the BYSTANDER effect, which is what makes a ~7 um penetration depth survivable, and the linker chemistry that decides whether it happens |
//! | [`oncolytic`] | Oncolytic virus SPREAD as a race between replication and antiviral clearance, gated by interferon competence -- so infection extent is DERIVED rather than assumed |
//! | [`oxygen`] | 3D radial oxygen gradients for spheroid tumors |
//! | [`radiation`] | Ionizing radiation: linear-quadratic DNA lethality and radiation-induced ferroptosis as SEPARATE channels, both O2-scaled by the OER |
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
//! | [`alox`] | ALOX lipoxygenase-isoform-specific peroxidation rate + MCFA→ACSL4 PUFA sensitization (off-by-default boosts) |
//! | [`acsl4`] | ACSL4-status biomarker stratification: tumor-intrinsic PUFA-incorporation gate (ACSL4-negative ⇒ ferroptosis-refractory) |
//! | [`ifngamma`] | IFN-γ → System Xc⁻ + ACSL4 ferroptosis-sensitization coupling (the immune-amplification return arm) |
//! | [`phenotype_mufa`] | Per-phenotype SCD1/MUFA accumulation rate + carrying-capacity multipliers |
//! | [`senescence`] | Therapy-induced-senescence ferroptosis program + SASP→immune coupling + CDK4/6-primed combination |
//! | [`repair`] | ESCRT-III membrane-repair brake on death EXECUTION (finite per-cell rescue budget delays death; off-by-default) |

pub mod acsl4;
pub mod alox;
pub mod cell;
pub mod checkpoint;
pub mod chemo;
pub mod clonal;
pub mod contact;
pub mod copper;
// Listed before `params` because `SpatialParams` holds a `Photosensitizer`.
pub mod ablation;
pub mod adc;
pub mod adoptive;
pub mod biochem;
pub mod dose_schedule;
pub mod drug_transport;
pub mod grid;
pub mod ifngamma;
pub mod immune;
pub mod immune_spatial;
pub mod io;
pub mod nutrient;
pub mod oncolytic;
pub mod oxygen;
pub mod params;
pub mod persister;
pub mod ph;
pub mod phenotype_mufa;
pub mod photosensitizer_pk;
pub mod physics;
pub mod radiation;
pub mod reaction_diffusion;
pub mod repair;
pub mod senescence;
pub mod slab;
pub mod sonodynamic;
pub mod spheroid;
pub mod stats;
pub mod stromal;
pub mod trigger_wave;
pub mod tumor_pk;
pub mod vasculature;
