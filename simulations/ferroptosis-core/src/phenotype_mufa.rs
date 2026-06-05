//! Phenotype-specific SCD1/MUFA accumulation rates (#363).
//!
//! Follow-up to #339, which delivered the kinetic/acute MUFA start
//! ([`crate::params::Params::mufa_acute_start`]) but used a single GLOBAL
//! `scd_mufa_rate` for the accumulation dynamics, so the acute-versus-established
//! MUFA build-up has the same time constant for every cell.
//!
//! Monounsaturated-fatty-acid (MUFA) enrichment confers a ferroptosis-resistant
//! cell state (Magtanong et al., Cell Chem Biol 2019, PMID 30686757; ACSL3-
//! dependent), and SCD1 (stearoyl-CoA desaturase 1) is the desaturase that
//! synthesizes that endogenous MUFA pool and itself protects cells from
//! ferroptosis (Tesfay et al., Cancer Res 2019, PMID 31270077). That enrichment
//! is plausibly phenotype-dependent: a drug-tolerant persister remodels its
//! lipidome under drug pressure differently than a proliferating glycolytic cell.
//! The DIRECTION is genuinely uncertain — drug-tolerant persisters are
//! simultaneously **GPX4-dependent / ferroptosis-vulnerable** (Hangauer et al.,
//! Nature 2017, PMID 29088702) AND can lean on lipid remodeling (SCD1/MUFA) to
//! survive — so this module exposes the per-phenotype rate as a configurable knob
//! rather than baking in one sign.
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

/// Per-phenotype multipliers on the SCD1/MUFA dynamics: the accumulation RATE
/// (the bare phenotype fields, scaling `Params::scd_mufa_rate`) and, separately
/// (#390), the carrying CAP / steady state (the `*_cap` fields, scaling the
/// effective `Cell::mufa_cap`, or the global `Params::scd_mufa_max` when no
/// per-cell cap is set). `1.0` for a phenotype on both axes ⇒ that phenotype
/// keeps the global rate and cap. Rate sets how FAST MUFA protection builds; cap
/// sets the steady state it saturates TOWARD — biologically distinct (a
/// phenotype can build fast but plateau low, or vice versa).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PhenotypeMufaConfig {
    /// Glycolytic (proliferating rim) cells — RATE multiplier.
    pub glycolytic: f64,
    /// OXPHOS (mid-zone) cells — RATE multiplier.
    pub oxphos: f64,
    /// Drug-tolerant persister cells — RATE multiplier.
    pub persister: f64,
    /// NRF2-high persister cells — RATE multiplier.
    pub persister_nrf2: f64,
    /// Non-tumor stromal cells (not dosed; included for completeness) — RATE.
    pub stromal: f64,
    /// Glycolytic cells — CAP (steady-state) multiplier (#390).
    pub glycolytic_cap: f64,
    /// OXPHOS cells — CAP multiplier (#390).
    pub oxphos_cap: f64,
    /// Persister cells — CAP multiplier (#390).
    pub persister_cap: f64,
    /// NRF2-high persister cells — CAP multiplier (#390).
    pub persister_nrf2_cap: f64,
    /// Stromal cells — CAP multiplier (#390).
    pub stromal_cap: f64,
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
            glycolytic_cap: 1.0,
            oxphos_cap: 1.0,
            persister_cap: 1.0,
            persister_nrf2_cap: 1.0,
            stromal_cap: 1.0,
        }
    }

    /// UNCALIBRATED, direction-anchored PLACEHOLDER (not a fitted result). Encodes
    /// ONE plausible hypothesis: drug-tolerant persister and OXPHOS-leaning cells
    /// build SCD1/MUFA protection somewhat FASTER than proliferating glycolytic
    /// cells (a lipid-remodeling survival response; SCD1 protects from ferroptosis,
    /// Tesfay 2019 PMID 31270077, and MUFA enrichment confers resistance,
    /// Magtanong 2019 PMID 30686757), while acknowledging the opposing pull that
    /// persisters are also GPX4-dependent/ferroptosis-vulnerable (Hangauer 2017
    /// PMID 29088702). The
    /// magnitudes are illustrative; calibrate against time-resolved per-phenotype
    /// MUFA lipidomics before reading any number from a run that uses this.
    pub fn literature() -> Self {
        Self {
            glycolytic: 1.0,
            oxphos: 1.2,
            persister: 1.5,
            persister_nrf2: 1.5,
            stromal: 1.0,
            // CAP placeholders: the drug-tolerant lipid-remodeled persister state
            // is hypothesized to sustain a modestly HIGHER MUFA steady state, and
            // OXPHOS slightly so. Even more uncertain than the rate direction —
            // purely illustrative, calibrate against per-phenotype MUFA lipidomics.
            glycolytic_cap: 1.0,
            oxphos_cap: 1.1,
            persister_cap: 1.3,
            persister_nrf2_cap: 1.3,
            stromal_cap: 1.0,
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

    /// The cap (steady-state) multiplier for a given phenotype (#390). Applied
    /// multiplicatively on the effective cap, so `1.0` leaves the cap untouched
    /// (preserving a spheroid radial cap, or the global `scd_mufa_max`).
    pub fn cap_multiplier(&self, phenotype: Phenotype) -> f64 {
        match phenotype {
            Phenotype::Glycolytic => self.glycolytic_cap,
            Phenotype::OXPHOS => self.oxphos_cap,
            Phenotype::Persister => self.persister_cap,
            Phenotype::PersisterNrf2 => self.persister_nrf2_cap,
            Phenotype::Stromal => self.stromal_cap,
        }
    }

    /// `true` when every multiplier (rate AND cap) is exactly `1.0`. A consumer
    /// uses this to skip the apply step entirely on the default path, keeping the
    /// production matrix byte-identical (no per-cell `mufa_rate`/`mufa_cap` set).
    pub fn is_identity(&self) -> bool {
        self.glycolytic == 1.0
            && self.oxphos == 1.0
            && self.persister == 1.0
            && self.persister_nrf2 == 1.0
            && self.stromal == 1.0
            && self.glycolytic_cap == 1.0
            && self.oxphos_cap == 1.0
            && self.persister_cap == 1.0
            && self.persister_nrf2_cap == 1.0
            && self.stromal_cap == 1.0
    }
}

impl Default for PhenotypeMufaConfig {
    fn default() -> Self {
        Self::identity()
    }
}

/// Set each tumor cell's per-cell SCD1/MUFA RATE
/// ([`crate::cell::Cell::mufa_rate`]) to `base_rate * config.rate_multiplier(phenotype)`,
/// and (#390) scale its CAP ([`crate::cell::Cell::mufa_cap`]) by
/// `config.cap_multiplier(phenotype)`.
///
/// `base_rate` is the global SCD1/MUFA rate (`params.scd_mufa_rate`) and
/// `base_max` the global cap (`params.scd_mufa_max`) the multipliers scale. Only
/// tumor cells are touched; stromal cells keep their default. Geometric, no RNG.
///
/// The CAP is applied **multiplicatively on the effective cap**: when
/// `cap_multiplier == 1.0` the cap is left UNTOUCHED (preserving a spheroid
/// radial `mufa_cap`, or the implicit global `scd_mufa_max` via `None`); when it
/// differs, the cap becomes `effective_cap * cap_multiplier`, where
/// `effective_cap = cell.mufa_cap.unwrap_or(base_max)`. So run this AFTER any
/// radial-cap layer (the spheroid) for the two to compose. With
/// [`PhenotypeMufaConfig::identity`] the rate is set to `Some(base_rate)`
/// (behaviorally identical to the global fallback) and the cap is untouched; for
/// a guaranteed byte-identical default, a consumer skips this call when
/// [`PhenotypeMufaConfig::is_identity`] holds.
pub fn apply_phenotype_mufa_3d(
    grid: &mut TumorGrid3D,
    base_rate: f64,
    base_max: f64,
    config: &PhenotypeMufaConfig,
) {
    for idx in 0..grid.cells.len() {
        apply_phenotype_mufa_at_3d(grid, idx, base_rate, base_max, config);
    }
}

/// Re-apply the per-cell phenotype MUFA rate + cap to a SINGLE cell index, from
/// that cell's CURRENT phenotype. No-op if the index is not a tumor cell.
///
/// Used after clonal repopulation (#266) revives a dead site as a fresh cell:
/// `gen_cell` resets `mufa_rate`/`mufa_cap`, so a revived cell would otherwise
/// lose its phenotype-specific rate (and cap). Re-deriving them from the revived
/// cell's phenotype keeps `clonal(repopulation)` + `phenotype_mufa` coherent —
/// the analogue of [`crate::contact::apply_contact_resistance_at_3d`] (#302).
/// (A revived cell's `mufa_cap` is `None` from `gen_cell`, so a non-1.0 cap
/// multiplier scales `base_max` for it — the spheroid radial cap is not restored
/// here, matching how the contact re-application uses the geometric setup value.)
pub fn apply_phenotype_mufa_at_3d(
    grid: &mut TumorGrid3D,
    idx: usize,
    base_rate: f64,
    base_max: f64,
    config: &PhenotypeMufaConfig,
) {
    let gc = &mut grid.cells[idx];
    if gc.is_tumor {
        gc.cell.mufa_rate = Some(base_rate * config.rate_multiplier(gc.phenotype));
        let cap_mul = config.cap_multiplier(gc.phenotype);
        if cap_mul != 1.0 {
            let effective_cap = gc.cell.mufa_cap.unwrap_or(base_max);
            gc.cell.mufa_cap = Some(effective_cap * cap_mul);
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
        // #390: and saturates at a higher MUFA cap.
        assert!(c.cap_multiplier(Phenotype::Persister) > c.cap_multiplier(Phenotype::Glycolytic));
    }

    #[test]
    fn rate_only_config_leaves_cap_untouched() {
        // A config that perturbs only the RATE (all cap multipliers 1.0) must NOT
        // touch `mufa_cap` — preserving a spheroid radial cap or the `None` global
        // fallback. This is the byte-identity-relevant cap invariant.
        let mut grid = TumorGrid3D::generate(11, 11, 11, 20.0, 42);
        let cfg = PhenotypeMufaConfig {
            persister: 3.0,
            ..PhenotypeMufaConfig::identity()
        };
        apply_phenotype_mufa_3d(&mut grid, 0.02, 0.25, &cfg);
        for gc in &grid.cells {
            assert_eq!(
                gc.cell.mufa_cap, None,
                "a rate-only config must leave every mufa_cap at its prior value (None here)"
            );
        }
    }

    #[test]
    fn apply_sets_per_cell_rate_and_cap_by_phenotype() {
        let mut grid = TumorGrid3D::generate(11, 11, 11, 20.0, 42);
        let cfg = PhenotypeMufaConfig {
            glycolytic: 1.0,
            oxphos: 2.0,
            persister: 3.0,
            persister_nrf2: 3.0,
            stromal: 1.0,
            glycolytic_cap: 1.0,
            oxphos_cap: 1.5,
            persister_cap: 2.0,
            persister_nrf2_cap: 2.0,
            stromal_cap: 1.0,
        };
        let base = 0.02;
        let base_max = 0.25;
        apply_phenotype_mufa_3d(&mut grid, base, base_max, &cfg);
        let mut saw_tumor = false;
        for gc in &grid.cells {
            if gc.is_tumor {
                saw_tumor = true;
                let expected_rate = base * cfg.rate_multiplier(gc.phenotype);
                assert_eq!(
                    gc.cell.mufa_rate,
                    Some(expected_rate),
                    "tumor cell rate must be base*rate_multiplier for its phenotype"
                );
                // Cap: started at None ⇒ effective base_max, scaled by cap_mul.
                // cap_mul == 1.0 ⇒ left untouched (None); otherwise base_max*cap_mul.
                let cap_mul = cfg.cap_multiplier(gc.phenotype);
                if cap_mul == 1.0 {
                    assert_eq!(gc.cell.mufa_cap, None, "cap_mul 1.0 must leave cap None");
                } else {
                    assert_eq!(
                        gc.cell.mufa_cap,
                        Some(base_max * cap_mul),
                        "cap must be base_max*cap_multiplier for its phenotype"
                    );
                }
            } else {
                assert_eq!(gc.cell.mufa_rate, None);
                assert_eq!(gc.cell.mufa_cap, None);
            }
        }
        assert!(saw_tumor, "test grid should contain tumor cells");
    }

    #[test]
    fn cap_multiplier_composes_with_an_existing_radial_cap() {
        // #390: the cap multiplier scales the EFFECTIVE cap, so a phenotype cap
        // applied on top of a spheroid-style radial cap multiplies it (rather than
        // replacing it). Here we pre-set a per-cell cap (as the spheroid would),
        // then apply a 2x persister cap and confirm it composed.
        let mut grid = TumorGrid3D::generate(11, 11, 11, 20.0, 42);
        let radial = 0.1;
        for gc in grid.cells.iter_mut() {
            if gc.is_tumor {
                gc.cell.mufa_cap = Some(radial);
            }
        }
        let cfg = PhenotypeMufaConfig {
            persister_cap: 2.0,
            ..PhenotypeMufaConfig::identity()
        };
        apply_phenotype_mufa_3d(&mut grid, 0.02, 0.25, &cfg);
        for gc in &grid.cells {
            if gc.is_tumor {
                let expected = if gc.phenotype == Phenotype::Persister {
                    Some(radial * 2.0) // composed with the pre-set radial cap
                } else {
                    Some(radial) // cap_mul 1.0 ⇒ untouched
                };
                assert_eq!(
                    gc.cell.mufa_cap, expected,
                    "persister cap must scale the pre-set radial cap; others untouched"
                );
            }
        }
    }

    #[test]
    fn at_helper_reapplies_rate_and_cap_to_a_reset_cell() {
        // Coherence with clonal repopulation (#266/#302): a revived dead site is a
        // fresh `gen_cell` with `mufa_rate`/`mufa_cap` reset. The per-index
        // re-apply must re-derive the per-phenotype rate (and cap) from the revived
        // cell's CURRENT phenotype, mirroring the contact re-application.
        let mut grid = TumorGrid3D::generate(11, 11, 11, 20.0, 42);
        let cfg = PhenotypeMufaConfig::literature();
        let base = 0.02;
        let base_max = 0.25;
        apply_phenotype_mufa_3d(&mut grid, base, base_max, &cfg);
        let idx = grid
            .cells
            .iter()
            .position(|gc| gc.is_tumor)
            .expect("test grid should contain a tumor cell");
        // Simulate the revival reset.
        grid.cells[idx].cell.mufa_rate = None;
        grid.cells[idx].cell.mufa_cap = None;
        apply_phenotype_mufa_at_3d(&mut grid, idx, base, base_max, &cfg);
        let gc = &grid.cells[idx];
        assert_eq!(
            gc.cell.mufa_rate,
            Some(base * cfg.rate_multiplier(gc.phenotype)),
            "the at-helper must re-derive the per-phenotype rate for a reset (revived) cell"
        );
        let cap_mul = cfg.cap_multiplier(gc.phenotype);
        let expected_cap = if cap_mul == 1.0 {
            None
        } else {
            Some(base_max * cap_mul)
        };
        assert_eq!(
            gc.cell.mufa_cap, expected_cap,
            "the at-helper must re-derive the per-phenotype cap for a reset (revived) cell"
        );
    }

    /// #390 acceptance: two phenotypes saturate at DIFFERENT MUFA steady states
    /// under the same dosing when a per-phenotype CAP multiplier is enabled. A
    /// cap-only config (persister cap 2x, others 1x) is applied to a grid; a
    /// persister cell and a glycolytic cell are then each run from the same naive
    /// MUFA start under `Params::spheroid()` + Control (alive, so MUFA accumulates
    /// to its cap). The 2x-cap persister must plateau HIGHER.
    #[test]
    fn phenotype_caps_yield_different_mufa_steady_states() {
        use crate::biochem::{sim_cell_step, CellState};
        use crate::cell::Treatment;
        use crate::params::Params;
        use rand::rngs::StdRng;
        use rand::SeedableRng;

        let params = Params::spheroid(); // MUFA active (scd_mufa_rate > 0)
        let mut grid = TumorGrid3D::generate(15, 15, 15, 20.0, 7);
        // Cap-only: persister saturates at 2x base_max, glycolytic at 1x.
        let cfg = PhenotypeMufaConfig {
            persister_cap: 2.0,
            ..PhenotypeMufaConfig::identity()
        };
        apply_phenotype_mufa_3d(&mut grid, params.scd_mufa_rate, params.scd_mufa_max, &cfg);

        let steady = |pheno: Phenotype| -> f64 {
            let gc = grid
                .cells
                .iter()
                .find(|g| g.is_tumor && g.phenotype == pheno)
                .expect("test grid should contain a tumor cell of this phenotype");
            let cell = gc.cell.clone();
            let mut init_rng = StdRng::seed_from_u64(5);
            let mut state = CellState::from_cell(&cell, Treatment::Control, &params, &mut init_rng);
            state.mufa_protection = 0.0; // same naive start for both
            let mut step_rng = StdRng::seed_from_u64(99);
            for step in 0..400 {
                sim_cell_step(&mut state, &cell, &params, step, 0.0, &mut step_rng);
            }
            state.mufa_protection
        };

        let persister_ss = steady(Phenotype::Persister);
        let glyco_ss = steady(Phenotype::Glycolytic);
        assert!(
            persister_ss > glyco_ss * 1.2,
            "the 2x-cap persister must saturate at a higher MUFA steady state than glycolytic: \
             persister={persister_ss}, glycolytic={glyco_ss}"
        );
    }
}
