//! Phenotype-specific SCD1/MUFA accumulation rates (#363).
//!
//! Follow-up to #339, which delivered the kinetic/acute MUFA start
//! ([`crate::params::Params::mufa_acute_start`]) but used a single GLOBAL
//! `scd_mufa_rate` for the accumulation dynamics, so the acute-versus-established
//! MUFA build-up has the same time constant for every cell.
//!
//! SCD1 (stearoyl-CoA desaturase) drives the monounsaturated-fatty-acid
//! enrichment that confers a ferroptosis-resistant cell state (Magtanong et al.,
//! Cell Chem Biol 2019, PMID 30686757), and that enrichment is plausibly
//! phenotype-dependent: a drug-tolerant persister remodels its lipidome under
//! drug pressure differently than a proliferating glycolytic cell. The DIRECTION
//! is genuinely uncertain — drug-tolerant persisters are simultaneously
//! **GPX4-dependent / ferroptosis-vulnerable** (Hangauer et al., Nature 2017,
//! PMID 29088702) AND can lean on lipid remodeling (SCD1/MUFA) to survive — so
//! this module exposes the per-phenotype rate as a configurable knob rather than
//! baking in one sign.
//!
//! [`PhenotypeMufaConfig`] is a small per-phenotype multiplier on the global
//! `scd_mufa_rate`. [`apply_phenotype_mufa_rates_3d`] writes the resulting
//! per-cell rate to [`crate::cell::Cell::mufa_rate`]. The default
//! [`PhenotypeMufaConfig::identity`] (all multipliers `1.0`) leaves every cell at
//! the global rate, so a consumer that opts out keeps the production matrix
//! byte-identical. [`PhenotypeMufaConfig::literature`] is an UNCALIBRATED,
//! direction-anchored placeholder, not a fitted result; calibrate against
//! time-resolved per-phenotype MUFA lipidomics.

use crate::cell::Phenotype;
use crate::grid::TumorGrid3D;

/// Per-phenotype multiplier on the global SCD1/MUFA accumulation rate
/// (`Params::scd_mufa_rate`). `1.0` for a phenotype ⇒ that phenotype keeps the
/// global rate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PhenotypeMufaConfig {
    /// Glycolytic (proliferating rim) cells.
    pub glycolytic: f64,
    /// OXPHOS (mid-zone) cells.
    pub oxphos: f64,
    /// Drug-tolerant persister cells.
    pub persister: f64,
    /// NRF2-high persister cells.
    pub persister_nrf2: f64,
    /// Non-tumor stromal cells (not dosed; included for completeness).
    pub stromal: f64,
}

impl PhenotypeMufaConfig {
    /// Identity: every phenotype keeps the global `scd_mufa_rate`. Applying this
    /// leaves the per-cell rate equal to the global rate, so a consumer that
    /// passes `base_rate = params.scd_mufa_rate` stays byte-identical.
    pub fn identity() -> Self {
        Self {
            glycolytic: 1.0,
            oxphos: 1.0,
            persister: 1.0,
            persister_nrf2: 1.0,
            stromal: 1.0,
        }
    }

    /// UNCALIBRATED, direction-anchored PLACEHOLDER (not a fitted result). Encodes
    /// ONE plausible hypothesis: drug-tolerant persister and OXPHOS-leaning cells
    /// build SCD1/MUFA protection somewhat FASTER than proliferating glycolytic
    /// cells (a lipid-remodeling survival response; Magtanong 2019 PMID 30686757),
    /// while acknowledging the opposing pull that persisters are also
    /// GPX4-dependent/ferroptosis-vulnerable (Hangauer 2017 PMID 29088702). The
    /// magnitudes are illustrative; calibrate against time-resolved per-phenotype
    /// MUFA lipidomics before reading any number from a run that uses this.
    pub fn literature() -> Self {
        Self {
            glycolytic: 1.0,
            oxphos: 1.2,
            persister: 1.5,
            persister_nrf2: 1.5,
            stromal: 1.0,
        }
    }

    /// The rate multiplier for a given phenotype.
    pub fn rate_multiplier(&self, phenotype: Phenotype) -> f64 {
        match phenotype {
            Phenotype::Glycolytic => self.glycolytic,
            Phenotype::OXPHOS => self.oxphos,
            Phenotype::Persister => self.persister,
            Phenotype::PersisterNrf2 => self.persister_nrf2,
            Phenotype::Stromal => self.stromal,
        }
    }

    /// `true` when every multiplier is exactly `1.0`. A consumer uses this to
    /// skip the apply step entirely on the default path, keeping the production
    /// matrix byte-identical (no per-cell `mufa_rate` is set).
    pub fn is_identity(&self) -> bool {
        self.glycolytic == 1.0
            && self.oxphos == 1.0
            && self.persister == 1.0
            && self.persister_nrf2 == 1.0
            && self.stromal == 1.0
    }
}

impl Default for PhenotypeMufaConfig {
    fn default() -> Self {
        Self::identity()
    }
}

/// Set each tumor cell's per-cell SCD1/MUFA accumulation rate
/// ([`crate::cell::Cell::mufa_rate`]) to `base_rate * config.rate_multiplier(phenotype)`.
///
/// `base_rate` is the global SCD1/MUFA rate the multipliers scale (the consumer
/// passes `params.scd_mufa_rate`). Only tumor cells are touched; stromal cells
/// keep their default. Geometric, no RNG. With [`PhenotypeMufaConfig::identity`]
/// every cell is set to `Some(base_rate)`, which is behaviorally identical to the
/// `None` fallback (`params.scd_mufa_rate`) when `base_rate == params.scd_mufa_rate`;
/// for a guaranteed byte-identical default, a consumer skips this call when
/// [`PhenotypeMufaConfig::is_identity`] holds.
pub fn apply_phenotype_mufa_rates_3d(
    grid: &mut TumorGrid3D,
    base_rate: f64,
    config: &PhenotypeMufaConfig,
) {
    for gc in grid.cells.iter_mut() {
        if gc.is_tumor {
            gc.cell.mufa_rate = Some(base_rate * config.rate_multiplier(gc.phenotype));
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn identity_is_all_unit_and_detected() {
        let c = PhenotypeMufaConfig::identity();
        assert!(c.is_identity());
        for p in [
            Phenotype::Glycolytic,
            Phenotype::OXPHOS,
            Phenotype::Persister,
            Phenotype::PersisterNrf2,
            Phenotype::Stromal,
        ] {
            assert_eq!(c.rate_multiplier(p), 1.0);
        }
        assert_eq!(
            PhenotypeMufaConfig::default(),
            PhenotypeMufaConfig::identity()
        );
    }

    #[test]
    fn literature_is_non_identity_and_differentiates_phenotypes() {
        let c = PhenotypeMufaConfig::literature();
        assert!(!c.is_identity());
        // The placeholder hypothesis: persister builds MUFA protection faster
        // than glycolytic (the direction is the point, not the magnitude).
        assert!(c.rate_multiplier(Phenotype::Persister) > c.rate_multiplier(Phenotype::Glycolytic));
        assert!(c.rate_multiplier(Phenotype::OXPHOS) > c.rate_multiplier(Phenotype::Glycolytic));
    }

    #[test]
    fn apply_sets_per_cell_rate_by_phenotype() {
        let mut grid = TumorGrid3D::generate(11, 11, 11, 20.0, 42);
        let cfg = PhenotypeMufaConfig {
            glycolytic: 1.0,
            oxphos: 2.0,
            persister: 3.0,
            persister_nrf2: 3.0,
            stromal: 1.0,
        };
        let base = 0.02;
        apply_phenotype_mufa_rates_3d(&mut grid, base, &cfg);
        let mut saw_tumor = false;
        for gc in &grid.cells {
            if gc.is_tumor {
                saw_tumor = true;
                let expected = base * cfg.rate_multiplier(gc.phenotype);
                assert_eq!(
                    gc.cell.mufa_rate,
                    Some(expected),
                    "tumor cell rate must be base*multiplier for its phenotype"
                );
            } else {
                // Stromal/non-tumor cells are left at their default (None).
                assert_eq!(gc.cell.mufa_rate, None);
            }
        }
        assert!(saw_tumor, "test grid should contain tumor cells");
    }
}
