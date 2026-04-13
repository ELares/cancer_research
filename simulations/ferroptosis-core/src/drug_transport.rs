//! Drug penetration modeling for tissue-specific pharmacokinetics.
//!
//! Models how drug concentration drops with distance from the nearest
//! blood vessel using an exponential decay approximation of the Krogh
//! cylinder steady-state solution.
//!
//! The key equation: `C(r) = C₀ × exp(-r / λ)` where `λ = √(D/k)` is
//! the penetration length, `D` is the effective diffusion coefficient,
//! and `k` is the total clearance rate (cellular uptake + metabolism).
//!
//! # References
//!
//! - Minchinton AI, Tannock IF. "Drug penetration in solid tumours."
//!   Nature Reviews Cancer 6:583-592, 2006.
//! - Thurber GM, et al. "Antibody tumor penetration." Advanced Drug
//!   Delivery Reviews 60:1421-1434, 2008.
//! - El-Kareh AW, Secomb TW. "A mathematical model for comparison of
//!   bolus injection, continuous infusion, and liposomal delivery of
//!   doxorubicin." Neoplasia 2:325-338, 2000.

use serde::{Deserialize, Serialize};

/// Drug physicochemical and pharmacokinetic parameters.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct DrugParams {
    /// Effective diffusion coefficient in tissue (cm²/s).
    /// Depends on drug molecular weight, charge, and tissue density.
    /// Small molecules (~500 Da): 1-10 × 10⁻⁷ cm²/s.
    /// Antibodies (~150 kDa): 0.1-1 × 10⁻⁷ cm²/s.
    pub diffusion_coeff_cm2_s: f64,

    /// Cellular uptake rate (1/s). How fast cells internalize the drug.
    pub uptake_rate: f64,

    /// Extracellular metabolism/degradation rate (1/s).
    pub metabolism_rate: f64,

    /// Drug-intrinsic bioavailability at the vessel wall (normalized, 0-1).
    /// Accounts for plasma protein binding and endothelial exclusion
    /// specific to the drug molecule. Set to 1.0 for freely permeable
    /// small molecules. This is multiplied by the tissue's vascular
    /// permeability to get the interstitial concentration, so do NOT
    /// duplicate the tissue permeability factor here.
    pub vessel_wall_conc: f64,

    /// Human-readable name for output.
    pub name: &'static str,
}

/// Tissue-specific transport parameters affecting drug penetration.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct TissueParams {
    /// Mean inter-vessel distance (μm). Determines how far a drug must
    /// diffuse to reach the most remote cells. Inversely related to
    /// vascular density. Typical: 100-300μm in solid tumors.
    pub inter_vessel_distance_um: f64,

    /// Vascular permeability factor (0-1). Fraction of vessel-wall
    /// concentration that reaches the interstitium. Reduced by tight
    /// junctions (BBB), elevated interstitial fluid pressure, etc.
    pub vascular_permeability: f64,

    /// Human-readable name for output.
    pub name: &'static str,
}

/// Characteristic penetration length (μm).
///
/// This is the distance at which drug concentration drops to 1/e (~37%)
/// of its vessel-wall value. Determined by the balance between diffusion
/// (spreading the drug) and clearance (cells consuming it).
///
/// `λ = √(D / k_total)` where `k_total = uptake_rate + metabolism_rate`.
pub fn penetration_length_um(drug: &DrugParams) -> f64 {
    let k_total = drug.uptake_rate + drug.metabolism_rate;
    if k_total <= 0.0 {
        return f64::INFINITY;
    }
    // Convert D from cm²/s to μm²/s (1 cm = 10⁴ μm, so 1 cm² = 10⁸ μm²)
    let d_um2_s = drug.diffusion_coeff_cm2_s * 1e8;
    (d_um2_s / k_total).sqrt()
}

/// Drug concentration at radial distance `r_um` (μm) from the nearest vessel.
///
/// Returns normalized concentration in [0, 1]. Uses the exponential decay
/// approximation of the Krogh cylinder steady-state solution, which is valid
/// when the vessel radius is much smaller than the tissue radius (typically
/// R_vessel ≈ 5-10μm vs R_tissue ≈ 50-150μm).
///
/// `C(r) = C_vessel × permeability × exp(-r / λ)`
pub fn concentration_at_distance(r_um: f64, drug: &DrugParams, tissue: &TissueParams) -> f64 {
    let lambda = penetration_length_um(drug);
    let c0 = drug.vessel_wall_conc * tissue.vascular_permeability;
    (c0 * (-r_um / lambda).exp()).min(1.0)
}

/// Maximum distance from a vessel (half the inter-vessel distance).
pub fn max_distance_um(tissue: &TissueParams) -> f64 {
    tissue.inter_vessel_distance_um / 2.0
}

/// Compute the concentration profile across the full radial range.
///
/// Returns `n_bins` evenly spaced `(distance_um, concentration)` pairs
/// from 0 to `max_distance_um`.
pub fn concentration_profile(
    drug: &DrugParams,
    tissue: &TissueParams,
    n_bins: usize,
) -> Vec<(f64, f64)> {
    let r_max = max_distance_um(tissue);
    (0..n_bins)
        .map(|i| {
            let r = r_max * i as f64 / (n_bins - 1).max(1) as f64;
            (r, concentration_at_distance(r, drug, tissue))
        })
        .collect()
}

// ============================================================
// Drug Presets
// ============================================================

/// RSL3-like small molecule GPX4 inhibitor.
///
/// MW ~500 Da. Moderate diffusion, moderate cellular uptake.
/// Penetration length ~100-120μm in well-vascularized tissue.
pub fn rsl3_like() -> DrugParams {
    DrugParams {
        // Small molecule in tissue: D ≈ 5 × 10⁻⁷ cm²/s
        // Ref: El-Kareh & Secomb 2000 (doxorubicin range 1-8 × 10⁻⁷)
        diffusion_coeff_cm2_s: 5.0e-7,
        // Moderate uptake: cells internalize but don't trap heavily
        uptake_rate: 0.004,
        // Low extracellular metabolism for a stable small molecule
        metabolism_rate: 0.001,
        // Freely permeable small molecule
        vessel_wall_conc: 1.0,
        name: "RSL3-like",
    }
}

/// Doxorubicin-like transport profile (penetration calibration reference).
///
/// Uses doxorubicin's well-characterized transport parameters
/// (MW ~540 Da, D ≈ 3×10⁻⁷ cm²/s, high uptake from DNA trapping)
/// to validate that the exponential model produces a penetration
/// length in the published 40-80μm range (Minchinton & Tannock 2006).
///
/// **Important:** This is a transport-only reference. The cell-level
/// pharmacology still uses the RSL3/GPX4-inhibition pathway, not
/// doxorubicin's actual mechanism (DNA intercalation, topoisomerase II).
/// Comparative kill rates between this and `rsl3_like()` reflect only
/// differences in tissue penetration depth, not drug mechanism.
pub fn doxorubicin_transport_reference() -> DrugParams {
    DrugParams {
        // Doxorubicin D ≈ 3 × 10⁻⁷ cm²/s in tissue
        // Ref: El-Kareh & Secomb 2000
        diffusion_coeff_cm2_s: 3.0e-7,
        // High uptake due to DNA binding/trapping
        uptake_rate: 0.01,
        // Moderate metabolism
        metabolism_rate: 0.002,
        // Freely permeable
        vessel_wall_conc: 1.0,
        name: "Doxorubicin-transport",
    }
}

// ============================================================
// Tissue Presets
// ============================================================

/// Well-vascularized epithelial tissue (breast, lung, colorectal).
///
/// Dense capillary network, moderate permeability.
/// Inter-vessel distance ~100-150μm.
pub fn epithelial_well_vascularized() -> TissueParams {
    TissueParams {
        inter_vessel_distance_um: 120.0,
        vascular_permeability: 0.8,
        name: "Epithelial (well-vascularized)",
    }
}

/// Poorly vascularized epithelial tissue (pancreatic, some liver).
///
/// Sparse vasculature, high interstitial fluid pressure, low permeability.
/// Inter-vessel distance ~200-300μm.
/// Ref: Olive et al., Science 2009 (pancreatic desmoplasia)
pub fn epithelial_poorly_vascularized() -> TissueParams {
    TissueParams {
        inter_vessel_distance_um: 250.0,
        vascular_permeability: 0.4,
        name: "Epithelial (poorly-vascularized)",
    }
}

/// CNS/neuroectodermal tissue (glioblastoma).
///
/// Blood-brain barrier severely restricts drug entry. Even disrupted BBB
/// in tumor core has much lower permeability than systemic vasculature.
/// Inter-vessel distance moderate but permeability very low.
/// Ref: Sarkaria et al., Neuro-Oncology 2018 (BBB and drug delivery)
pub fn neuroectodermal_cns() -> TissueParams {
    TissueParams {
        inter_vessel_distance_um: 150.0,
        vascular_permeability: 0.15,
        name: "Neuroectodermal (CNS/BBB)",
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn concentration_at_vessel_wall_is_c0() {
        let drug = rsl3_like();
        let tissue = epithelial_well_vascularized();
        let c = concentration_at_distance(0.0, &drug, &tissue);
        let expected = drug.vessel_wall_conc * tissue.vascular_permeability;
        assert!((c - expected).abs() < 1e-10, "At r=0: expected {expected}, got {c}");
    }

    #[test]
    fn concentration_decays_with_distance() {
        let drug = rsl3_like();
        let tissue = epithelial_well_vascularized();
        let c_near = concentration_at_distance(10.0, &drug, &tissue);
        let c_far = concentration_at_distance(100.0, &drug, &tissue);
        assert!(c_far < c_near, "Concentration should decrease with distance");
    }

    #[test]
    fn concentration_near_zero_at_large_distance() {
        let drug = rsl3_like();
        let tissue = epithelial_well_vascularized();
        let c = concentration_at_distance(1000.0, &drug, &tissue);
        assert!(c < 0.01, "Concentration at 1mm should be negligible, got {c}");
    }

    #[test]
    fn penetration_length_scales_with_sqrt_d_over_k() {
        let drug1 = DrugParams {
            diffusion_coeff_cm2_s: 4.0e-7,
            uptake_rate: 0.004,
            metabolism_rate: 0.0,
            vessel_wall_conc: 1.0,
            name: "test1",
        };
        // 4× diffusion should give 2× penetration length (sqrt scaling)
        let drug2 = DrugParams {
            diffusion_coeff_cm2_s: 16.0e-7,
            ..drug1.clone()
        };
        let ratio = penetration_length_um(&drug2) / penetration_length_um(&drug1);
        assert!((ratio - 2.0).abs() < 0.01, "λ should scale as √D: ratio={ratio}");
    }

    #[test]
    fn doxorubicin_transport_penetration_matches_literature() {
        // Minchinton & Tannock 2006: doxorubicin penetrates ~40-80μm
        let drug = doxorubicin_transport_reference();
        let lambda = penetration_length_um(&drug);
        assert!(
            lambda > 30.0 && lambda < 120.0,
            "Doxorubicin transport λ should be ~50-80μm, got {lambda:.1}μm"
        );
    }

    #[test]
    fn bbb_reduces_effective_concentration() {
        let drug = rsl3_like();
        let normal = epithelial_well_vascularized();
        let cns = neuroectodermal_cns();
        let c_normal = concentration_at_distance(50.0, &drug, &normal);
        let c_cns = concentration_at_distance(50.0, &drug, &cns);
        assert!(
            c_cns < c_normal * 0.5,
            "BBB should substantially reduce concentration: normal={c_normal:.3}, cns={c_cns:.3}"
        );
    }

    #[test]
    fn profile_has_correct_length() {
        let drug = rsl3_like();
        let tissue = epithelial_well_vascularized();
        let profile = concentration_profile(&drug, &tissue, 50);
        assert_eq!(profile.len(), 50);
        assert!((profile[0].0 - 0.0).abs() < 1e-10, "First point should be at r=0");
        let r_max = max_distance_um(&tissue);
        assert!((profile[49].0 - r_max).abs() < 1e-10, "Last point should be at r_max");
    }

    #[test]
    fn zero_clearance_gives_infinite_penetration() {
        let drug = DrugParams {
            diffusion_coeff_cm2_s: 5.0e-7,
            uptake_rate: 0.0,
            metabolism_rate: 0.0,
            vessel_wall_conc: 1.0,
            name: "no-clearance",
        };
        let lambda = penetration_length_um(&drug);
        assert!(lambda.is_infinite(), "Zero clearance should give infinite penetration");
    }
}
