//! Ferroptosis Simulation Suite — shared library.
//!
//! Provides the core biochemistry engine, cell types, physics models,
//! spatial grid, immune cascade, and I/O utilities used by all simulation binaries.

pub mod cell;
pub mod params;
pub mod biochem;
pub mod stats;
pub mod physics;
pub mod grid;
pub mod immune;
pub mod io;
