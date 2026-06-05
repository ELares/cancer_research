//! Spatial immune coupling: DAMP diffusion + activation primitives.
//!
//! Ferroptotic cell death releases damage-associated molecular patterns
//! (DAMPs) — calreticulin, HMGB1, ATP — that diffuse through the
//! extracellular space, activate nearby dendritic cells via TLR/RAGE
//! signaling, and prime T-cell-mediated tumor killing. This module
//! holds the spatial primitives that consumers (#195 sim-tme-3d, plus
//! sim-tme via the dim-agnostic helpers) compose into a full
//! immune-step orchestrator.
//!
//! **Module naming (renamed from `immune_3d` in #224).** The 3D-specific
//! grid helpers ([`diffuse_damp_3d_step`]) live here, but so do the
//! dimensionality-agnostic scalar helpers ([`dc_activation`],
//! [`immune_kill_probability`]) and the cross-dimension
//! [`DAMP_KILL_THRESHOLD`] constant — both 2D (sim-tme) and 3D
//! (sim-tme-3d) consume the same scalars. Calling the module
//! `immune_spatial` honestly covers both: the 3D-shaped grid helpers
//! keep their `_3d` suffix on the function name, the scalar helpers
//! work in either geometry.
//!
//! **Scope vs the existing [`crate::immune`] module.** That module is the
//! *dimensionless* single-event ICD cascade (one death → one DAMP burst →
//! one DC activation → one T-cell kill). This module is the *spatial*
//! complement: how DAMPs diffuse across a spheroid grid and how local
//! DAMP concentration drives per-cell kill probability. The two compose;
//! `immune` answers "what does one death contribute?" and `immune_spatial`
//! answers "where does it spread and who does it affect?"
//!
//! **The 104:1 question (issue #188).** Sim-tme's 2D model finds SDT
//! produces ~104× more immune kills than RSL3 because SDT's dense kill
//! field creates a high local DAMP concentration. The issue asks
//! whether this ratio holds in 3D. **Answering it requires a full
//! multi-step simulation** (sim-tme-3d, #195) — not a library unit
//! test. This module provides the diffusion primitive; the kill-ratio
//! comparison lands with #196 (3D validation).
//!
//! ## ⚠️ Stability requirement (critical for 3D)
//!
//! [`diffuse_damp_3d_step`] is mathematically stable only when
//! `diffusion_fraction × max_neighbor_count < 1`. In 3D with up to 26
//! Moore neighbors, that means **`diffusion_fraction < ≈ 0.038`**.
//!
//! Sim-tme's 2D default `0.08` is **unsafe in 3D** (0.08 × 26 = 2.08 > 1):
//! the source cell would lose more DAMP per step than it has, the
//! defensive `.max(0.0)` clamp would destroy mass, and the field would
//! silently produce nonsense. **A regular `assert!` (NOT
//! `debug_assert!`) enforces the stability check in both debug and
//! release** — the silent-failure mode is too pernicious for
//! debug-only catching. A consumer porting sim-tme's parameters verbatim
//! panics immediately with a clear stability message instead of
//! silently corrupting the field.
//!
//! **Suggested 3D-safe value: 0.025**, which gives the same per-step
//! total-diffusion fraction as 2D's `0.08 × 8 = 0.64` (compare to
//! `0.025 × 26 = 0.65`).
//!
//! ## API design — scratch-buffer pattern
//!
//! [`diffuse_damp_3d_step`] takes both `damp_field: &mut [f64]` AND
//! `scratch: &mut [f64]` (both length `grid.cells.len()`). The scratch
//! buffer ensures the spread step is order-independent (otherwise
//! source-ordering would bias the result). **Allocate scratch ONCE
//! before the simulation loop and reuse per step** — at 100³ × 180
//! steps, per-step allocation would be ~1.4 GB of churn.
//!
//! ## Quick example
//!
//! ```
//! use ferroptosis_core::grid::TumorGrid3D;
//! use ferroptosis_core::immune_spatial::{diffuse_damp_3d_step, dc_activation, immune_kill_probability};
//!
//! let g = TumorGrid3D::generate(10, 10, 10, 20.0, 42);
//! let n = g.cells.len();
//!
//! // Allocate state once, reuse per step.
//! let mut damp_field = vec![0.0_f64; n];
//! let mut scratch = vec![0.0_f64; n];
//!
//! // Simulate a single death at the center: DAMP burst.
//! let center_idx = g.flat_index(5, 5, 5);
//! damp_field[center_idx] = 10.0;
//!
//! // 3D-safe diffusion_fraction (see stability requirement above).
//! diffuse_damp_3d_step(&mut damp_field, &mut scratch, &g, 0.025, 0.03);
//!
//! // Per-cell immune activation and kill probability.
//! for &local_damp in &damp_field {
//!     let activation = dc_activation(local_damp, 50.0);
//!     let kill_prob = immune_kill_probability(activation, 0.02, 0.21);
//!     // consumer rolls: if rng.gen() < kill_prob { ... }
//!     let _ = (activation, kill_prob);
//! }
//! ```

use crate::grid::{TumorGrid3D, TUMOR_RADIUS_FRACTION};
use rand::prelude::*;

/// Minimum local-DAMP concentration above which a cell is eligible for
/// immune-mediated kill. Cells with `local_damp < DAMP_KILL_THRESHOLD`
/// are skipped — both as a performance optimization (don't roll RNG
/// when activation × kill_rate would round to ~0 anyway) and as a
/// floor against numerical noise in the diffused DAMP field.
///
/// Used by both sim-tme (2D) and sim-tme-3d (3D) so the kill-eligibility
/// floor stays consistent across dimensionality. The module name
/// `immune_spatial` covers both geometries — see [`immune_kill_probability`]
/// for the matching dimensionality-agnostic helper.
pub const DAMP_KILL_THRESHOLD: f64 = 0.01;

/// Maximum Moore-neighbor count in 3D (3×3×3 cube − self).
///
/// **Coupled to [`TumorGrid3D::neighbors`]**: if a future refactor of
/// `neighbors` changes the neighborhood radius (e.g. Moore-2 with 124
/// neighbors), this constant must be updated together — otherwise the
/// stability check in [`diffuse_damp_3d_step`] would silently
/// under-estimate the per-step source loss and admit unstable
/// `diffusion_fraction` values.
const MAX_3D_NEIGHBORS: usize = 26;

/// Source cells below this DAMP value skip spread + self-decrement
/// (matches sim-tme's `if local < 0.001 { continue; }` sub-threshold cutoff in
/// its DAMP-diffusion step). The cutoff exists for performance: cells
/// at this magnitude contribute negligibly to neighbors. Mass is exactly
/// preserved for sub-threshold cells (they don't lose to anyone), so
/// the field's long-run behavior is unchanged — just faster.
const DIFFUSION_SOURCE_CUTOFF: f64 = 0.001;

/// One step of DAMP diffusion + exponential clearance on a 3D spheroid
/// grid. Mutates `damp_field` in place using `scratch` to avoid
/// source-ordering bias (each step, every source spreads `share = local
/// × diffusion_fraction` to each of its up-to-26 Moore neighbors, then
/// every cell decays by `(1 − clearance_rate)`).
///
/// **Stability requirement** (see module doc): `diffusion_fraction ×
/// MAX_3D_NEIGHBORS < 1.0`. Enforced by a regular `assert!` that fires
/// in **both debug and release** — the silent-failure mode of mass
/// destruction via the defensive `.max(0.0)` clamp is too pernicious
/// for debug-only catching. Per-call cost is one multiply + one compare
/// (negligible). **Suggested 3D-safe value: `0.025`** (matches 2D's
/// per-step total diffusion of ~64%).
///
/// **Length contract**: `damp_field.len() == scratch.len() == grid.cells.len()`.
/// Validated with `assert!` (matches stromal pattern — programming-
/// contract bug deserves a clear release-mode panic).
///
/// **Cost**: O(N × 26) for N = `grid.cells.len()`. A per-step full-grid
/// sweep; the same hoist-constants-once pattern as the `oxygen`/`ph` radial
/// fields (landed in #289) applies if this ever needs tightening.
/// A source cell with DAMP below [`DIFFUSION_SOURCE_CUTOFF`] = 0.001
/// skips both spread and self-decrement (matches sim-tme's 2D
/// optimization; sub-threshold cells contribute negligibly and skipping
/// them preserves their mass exactly).
pub fn diffuse_damp_3d_step(
    damp_field: &mut [f64],
    scratch: &mut [f64],
    grid: &TumorGrid3D,
    diffusion_fraction: f64,
    clearance_rate: f64,
) {
    let n = grid.cells.len();
    assert!(
        damp_field.len() == n,
        "diffuse_damp_3d_step: damp_field.len() {} must equal grid.cells.len() {}",
        damp_field.len(),
        n
    );
    assert!(
        scratch.len() == n,
        "diffuse_damp_3d_step: scratch.len() {} must equal grid.cells.len() {}",
        scratch.len(),
        n
    );
    debug_assert!(
        diffusion_fraction.is_finite() && diffusion_fraction >= 0.0,
        "diffuse_damp_3d_step: diffusion_fraction must be finite and ≥ 0, got {diffusion_fraction}"
    );
    // **Stability check is `assert!` (not `debug_assert!`)** because the
    // failure mode is silent: violation causes the source's `local + delta`
    // to go negative, then the defensive `.max(0.0)` clamp destroys mass
    // with no panic. A consumer porting sim-tme's 2D default (0.08) would
    // get wrong DAMP fields in release with no indication. The per-call
    // cost (one multiply + one compare) is negligible.
    assert!(
        diffusion_fraction * (MAX_3D_NEIGHBORS as f64) < 1.0,
        "diffuse_damp_3d_step: diffusion_fraction × 26 = {} must be < 1 for stability. \
         Got diffusion_fraction = {}. Sim-tme's 2D default (0.08) is UNSAFE in 3D \
         (0.08 × 26 = 2.08 > 1) — use ≤ 0.038, suggested 0.025 to match 2D's per-step total.",
        diffusion_fraction * (MAX_3D_NEIGHBORS as f64),
        diffusion_fraction
    );
    debug_assert!(
        clearance_rate.is_finite() && (0.0..=1.0).contains(&clearance_rate),
        "diffuse_damp_3d_step: clearance_rate must be in [0, 1], got {clearance_rate}"
    );

    // Compute the spread into scratch.
    scratch.fill(0.0);
    for r in 0..grid.rows {
        for c in 0..grid.cols {
            for l in 0..grid.layers {
                let idx = grid.flat_index(r, c, l);
                let local = damp_field[idx];
                if local < DIFFUSION_SOURCE_CUTOFF {
                    continue;
                }
                let share = local * diffusion_fraction;
                let (neighbors, count) = grid.neighbors(r, c, l);
                for &(nr, nc, nl) in &neighbors[..count] {
                    let nidx = grid.flat_index(nr, nc, nl);
                    scratch[nidx] += share;
                }
                // Source loses exactly `share × count` (where count is the
                // ACTUAL neighbor count, not always 26 at boundaries).
                // This is what makes mass conservation hold across edges.
                scratch[idx] -= share * count as f64;
            }
        }
    }

    // Apply spread + clearance.
    for i in 0..n {
        // `.max(0.0)` is unreachable under the stability invariant
        // (`fraction × 26 < 1` guarantees `share × count < local`, so
        // `local + scratch[i] > 0` after diffusion). Kept as defensive
        // safety net for two cases: (1) caller-mutated negative inputs,
        // (2) future relaxation of the stability `assert!`. Cheap.
        damp_field[i] = (damp_field[i] + scratch[i]).max(0.0);
        damp_field[i] *= 1.0 - clearance_rate;
    }
}

/// Dendritic-cell activation as a function of local DAMP concentration.
///
/// Michaelis-Menten saturation: `activation = damp / (damp + kd)`, in
/// `[0, 1]` for non-negative inputs. At `damp = 0`: returns `0.0`
/// exactly (IEEE: `0/kd = 0`). At `damp = kd`: returns `0.5` exactly
/// (IEEE: `x/(2x) = 0.5` for finite x). For `damp >> kd`: asymptotes
/// to `1.0`.
///
/// **Pure scalar function.** No clamp; trust caller. Sim-tme uses
/// `kd = 50.0` (the `dc_activation_kd` field of `ImmuneConfig`).
///
/// **Release behavior for invalid inputs**:
///
/// | Bad input | Output |
/// |-----------|--------|
/// | `damp = NaN` | `NaN` (propagates) |
/// | `kd = NaN` | `NaN` |
/// | `damp + kd = 0` (both 0) | `NaN` (`0/0`) |
/// | `damp < 0`, `damp + kd < 0` | negative; consumer's problem |
/// | `kd < 0` | math still works but biologically meaningless |
///
/// Consumers passing finite `damp >= 0` and `kd > 0` (the normal regime)
/// don't hit any of these.
#[inline]
#[must_use = "the activation is the function's only output; ignoring it suggests a logic bug"]
pub fn dc_activation(local_damp: f64, kd: f64) -> f64 {
    local_damp / (local_damp + kd)
}

/// Per-cell immune kill probability per step.
///
/// `probability = (activation × kill_rate × (1 − effective_brake)).min(0.99)`
///
/// Dimensionality-agnostic (pure scalar — no grid). Called by both
/// sim-tme (2D) and sim-tme-3d (3D); lives in `immune_spatial` because
/// that module owns the shared spatial-immune primitives across both
/// geometries.
///
/// The `.min(0.99)` cap matches sim-tme: even at full activation with no
/// PD-1 brake, kills are never guaranteed (preserves stochasticity over
/// long simulations). At `activation = 0` or `kill_rate = 0` or
/// `effective_brake = 1`: returns `0.0` exactly.
///
/// **Pure scalar function.** Sim-tme uses `kill_rate = 0.02` and computes
/// `effective_brake = pd1_brake × (1 - anti_pd1_efficacy)` upstream.
///
/// **Lower bound**: NOT clamped at zero. For pathological inputs (e.g.,
/// `activation < 0`) the return value can be negative. Caller's
/// `rng.gen::<f64>() < probability` correctly produces zero kills in
/// that case (random in `[0, 1)` is never less than a negative), so
/// downstream semantics are safe — but the return-type contract is
/// `(-∞, 0.99]`, not `[0, 0.99]`.
///
/// **Release behavior for invalid inputs**:
///
/// | Bad input | Output |
/// |-----------|--------|
/// | any `NaN` argument | `NaN` (NaN-preserving cap; see implementation note) |
/// | `activation < 0` or `effective_brake > 1` | negative; capped only at upper bound |
///
/// **Implementation note on the cap.** Cannot use `.min(0.99)`: Rust's
/// `f64::min` treats NaN as the maximum (`NaN.min(0.99) = 0.99`), which
/// would silently convert a NaN input into a near-certain kill
/// probability. Using an explicit `if raw > 0.99` branch instead, which
/// is NaN-preserving (`NaN > 0.99` is `false`, so the `else` branch
/// returns `raw` = NaN). Caller-visible contract: NaN in → NaN out, so
/// downstream `rng.gen() < NaN` correctly produces zero kills.
#[inline]
#[must_use = "the kill probability is the function's only output; ignoring it suggests a logic bug"]
pub fn immune_kill_probability(activation: f64, kill_rate: f64, effective_brake: f64) -> f64 {
    let raw = activation * kill_rate * (1.0 - effective_brake);
    // NaN-preserving cap (cannot use raw.min(0.99) — f64::min treats
    // NaN as max). See implementation note in the rustdoc.
    if raw > 0.99 {
        0.99
    } else {
        raw
    }
}

/// T-cell exhaustion multiplier (#243, Phase 1).
///
/// Sustained cytotoxic activity in a region drives resident T cells toward a
/// dysfunctional ("exhausted") state that lowers their per-encounter kill
/// probability (Wherry, Nat Immunol 2011; Snell et al., Cell 2018). Modeled
/// as hyperbolic decay in the number of immune kills accumulated in the
/// cell's Moore neighborhood:
///
/// ```text
/// factor = 1 / (1 + exhaustion_rate · cumulative_kills)
/// ```
///
/// A consumer multiplies this into [`immune_kill_probability`]. Returns
/// exactly `1.0` (no suppression) when `exhaustion_rate == 0.0` OR
/// `cumulative_kills == 0`, so the default config (`exhaustion_rate = 0.0`)
/// is a no-op and the consumer's output stays byte-identical. Monotonically
/// non-increasing in both arguments, bounded in `(0, 1]`.
#[inline]
#[must_use]
pub fn exhaustion_factor(cumulative_kills: u32, exhaustion_rate: f64) -> f64 {
    debug_assert!(
        exhaustion_rate >= 0.0,
        "exhaustion_rate must be >= 0; a negative rate pushes the factor above 1 or negative"
    );
    1.0 / (1.0 + exhaustion_rate * cumulative_kills as f64)
}

// =====================================================================
// Treg / MDSC immunosuppressor field (#264, immune realism Phase 2)
// =====================================================================

/// Configuration for the Treg/MDSC immunosuppressor field (#264 Phase 2): a
/// second diffusing field that locally scales immune kill probability DOWN,
/// opposing the DAMP→kill effect. Sources (Treg/MDSC niches) replenish the
/// field each step; it diffuses and clears like the DAMP field.
///
/// **Off by default**: [`disabled`](Self::disabled) (`suppression_strength = 0`)
/// makes [`suppressor_kill_multiplier`] return `1.0`, so a consumer that does
/// not opt in stays byte-identical. The diffusion/clearance reuse
/// [`diffuse_damp_3d_step`], so `diffusion_fraction` is bound by the same
/// `× 26 < 1` stability requirement.
#[derive(Clone, Copy, Debug)]
pub struct SuppressorConfig {
    /// Strength of local kill suppression: `kill ×= 1/(1 + strength · field)`.
    /// `0.0` ⇒ no suppression (identity / byte-identical).
    pub suppression_strength: f64,
    /// Per-step suppressor released at each source (Treg/MDSC) cell, before
    /// diffusion. The field is clamped to `[0, 1]`.
    pub replenish_rate: f64,
    /// Fraction of suppressor shared with each Moore-26 neighbor per step
    /// (reuses [`diffuse_damp_3d_step`]; must satisfy `× 26 < 1`).
    pub diffusion_fraction: f64,
    /// Exponential clearance per step (Treg/MDSC turnover + drainage).
    pub clearance_rate: f64,
    /// Radius (µm) around each seed point within which tumor cells are marked
    /// as suppressor sources (the niche size).
    pub niche_radius_um: f64,
    /// Number of heuristic source seed points when no vessel positions are
    /// supplied. Ignored in perivascular mode (vessels are the seed points).
    pub n_sources: usize,
}

impl SuppressorConfig {
    /// Disabled: `suppression_strength = 0` ⇒ [`suppressor_kill_multiplier`]
    /// is the identity, so the consumer's output is byte-identical.
    pub fn disabled() -> Self {
        SuppressorConfig {
            suppression_strength: 0.0,
            replenish_rate: 0.0,
            diffusion_fraction: 0.025,
            clearance_rate: 0.03,
            niche_radius_um: 60.0,
            n_sources: 8,
        }
    }

    /// Literature-informed enabled config (placeholders pending calibration):
    /// moderate suppression, perivascular niches ~60 µm, slow turnover. Refs:
    /// Tauriello et al., Nature 2018 (TGFβ–Treg exclusion axis).
    pub fn enabled() -> Self {
        SuppressorConfig {
            suppression_strength: 6.0,
            replenish_rate: 0.15,
            diffusion_fraction: 0.025,
            clearance_rate: 0.03,
            niche_radius_um: 60.0,
            n_sources: 8,
        }
    }

    /// True when this config has no effect (no suppression).
    pub fn is_disabled(&self) -> bool {
        self.suppression_strength == 0.0
    }
}

/// Local immune-kill multiplier from the suppressor field: `1/(1 + strength · s)`
/// (mirrors [`exhaustion_factor`]'s form). Returns exactly `1.0` (no
/// suppression) when `suppression_strength == 0.0` OR `local_suppressor == 0.0`,
/// so the disabled config is a no-op and the consumer stays byte-identical.
/// Monotonically non-increasing in both arguments, bounded in `(0, 1]`.
#[inline]
#[must_use]
pub fn suppressor_kill_multiplier(local_suppressor: f64, suppression_strength: f64) -> f64 {
    debug_assert!(
        suppression_strength >= 0.0 && local_suppressor >= 0.0,
        "suppressor inputs must be >= 0; got strength={suppression_strength}, local={local_suppressor}"
    );
    1.0 / (1.0 + suppression_strength * local_suppressor)
}

/// Local immune-kill multiplier from the IMMUNOSUPPRESSIVE arm of ferroptosis
/// (#337): `1/(1 + strength · local_damp)` (mirrors [`suppressor_kill_multiplier`]
/// and [`exhaustion_factor`]).
///
/// The model otherwise treats ferroptotic ICD as net pro-immune (DAMP ->
/// `dc_activation` -> CD8 kill). But ferroptotic cells also co-release factors
/// that BLUNT dendritic-cell maturation and CD8 priming: extracellular GPX4
/// binding DC ZP3 (Liu et al., Cell 2026, PMID 41494530), oxidized lipids and
/// PGE2 (Kim et al., Nature 2022, PMID 36385526), so ferroptotic cells impede
/// DC-mediated anti-tumor immunity (Wiernicki et al., Nat Commun 2022, PMID
/// 35760796; review Tang et al., Immunol Rev 2024, PMID 37424139). Keying the
/// suppression on the SAME `local_damp` ferroptotic-death signal that drives
/// `dc_activation` means the two arms compose: at low death density the
/// pro-immune `activation` term dominates, but as density rises `activation`
/// saturates while this multiplier keeps falling, so the NET immune effect can
/// flip from pro- to anti-tumor.
///
/// Returns exactly `1.0` when `strength == 0.0` (the default) OR
/// `local_damp == 0.0`, so the consumer stays byte-identical. Monotonically
/// non-increasing in both arguments, bounded in `(0, 1]`. The direction is
/// timing-dependent (a small early-ferroptotic fraction can be immunogenic,
/// Efimova 2020 PMID 33188036), so this should dominate only at sustained/high
/// death density; the magnitude is an uncalibrated placeholder, the sign is the
/// result.
#[inline]
#[must_use]
pub fn ferroptotic_immunosuppression(local_damp: f64, strength: f64) -> f64 {
    debug_assert!(
        strength >= 0.0 && local_damp >= 0.0,
        "ferroptotic immunosuppression inputs must be >= 0; got strength={strength}, local_damp={local_damp}"
    );
    1.0 / (1.0 + strength.max(0.0) * local_damp)
}

/// Local immune-kill multiplier from the diffusing **SASP field** (#376).
///
/// The senescence SASP is a SECRETED, paracrine signal, so unlike the
/// cell-autonomous [`crate::senescence::sasp_immune_multiplier`] (#341, which
/// only scales a senescent cell's OWN kill probability) this modulates EVERY
/// cell exposed to the field — including adjacent NON-senescent tumor cells the
/// secretome protects or sensitizes. The consumer (sim-tme-3d) seeds the field
/// at the senescent cells (the senescence mask), diffuses it with
/// [`diffuse_damp_3d_step`], then applies this multiplier at each cell.
///
/// `strength` is SIGNED to carry the SASP's documented bidirectionality, the
/// same axis dependence as [`crate::senescence::SenescenceConfig::sasp_immune_mult`]:
/// - `strength > 0` ⇒ the immunosuppressive arm: `1/(1 + strength · field)`
///   LOWERS the local immune-kill probability (IL-1RA protecting neighboring
///   tumor cells, recruitment of surveillance-antagonizing myeloid cells — Di
///   Mitri 2014 PMID 25156255; Eggert 2016 PMID 27728804). Mirrors
///   [`suppressor_kill_multiplier`] / [`ferroptotic_immunosuppression`].
/// - `strength < 0` ⇒ the surveillance arm: `1 + (-strength) · field` RAISES it
///   (CD4 T-cell / monocyte recruitment to the senescent niche — Kang 2011 PMID
///   22080947).
/// - `strength == 0.0` (default) OR `field == 0.0` ⇒ exactly `1.0` ⇒ the consumer
///   stays byte-identical.
///
/// Bounded in `(0, ∞)` (immunosuppressive arm in `(0, 1]`; surveillance arm in
/// `[1, ∞)`). The magnitude is an UNCALIBRATED placeholder; the bidirectional
/// structure and the neighbor/bystander reach (vs the #341 cell-autonomous
/// multiplier) are the result.
#[inline]
#[must_use]
pub fn sasp_field_kill_multiplier(local_sasp: f64, strength: f64) -> f64 {
    debug_assert!(
        local_sasp >= 0.0 && local_sasp.is_finite() && strength.is_finite(),
        "sasp_field inputs must be finite and field >= 0; got field={local_sasp}, strength={strength}"
    );
    if strength == 0.0 || local_sasp == 0.0 {
        1.0
    } else if strength > 0.0 {
        1.0 / (1.0 + strength * local_sasp)
    } else {
        1.0 + (-strength) * local_sasp
    }
}

/// Per-cell boolean mask of Treg/MDSC suppressor **source** cells (#264).
///
/// Seed points are **perivascular** when `vessels` is supplied (Tregs cluster
/// at perivascular niches — Tauriello 2018), reusing the vessel lattice
/// positions; otherwise `cfg.n_sources` **heuristic** points are sampled
/// uniformly in the tumor sphere via an **independent** `StdRng(seed)` (so the
/// cell grid's RNG stream is untouched). A tumor cell is a source iff it lies
/// within `cfg.niche_radius_um` of any seed point. Non-tumor cells are never
/// sources. Deterministic given `(grid dims, cfg, vessels, seed)`.
pub fn suppressor_source_mask_3d(
    grid: &TumorGrid3D,
    cfg: &SuppressorConfig,
    vessels: Option<&[(f64, f64, f64)]>,
    seed: u64,
) -> Vec<bool> {
    // Seed points in lattice (cell) coordinates — same convention as
    // `vasculature::place_vessels_3d`.
    let seed_points: Vec<(f64, f64, f64)> = match vessels {
        Some(v) if !v.is_empty() => v.to_vec(),
        _ => {
            let mut rng = StdRng::seed_from_u64(seed);
            let center = (
                grid.rows as f64 / 2.0,
                grid.cols as f64 / 2.0,
                grid.layers as f64 / 2.0,
            );
            let tumor_radius =
                (grid.rows.min(grid.cols).min(grid.layers) as f64) * TUMOR_RADIUS_FRACTION;
            (0..cfg.n_sources.max(1))
                .map(|_| {
                    // Uniform-in-sphere (cbrt radial), matching place_vessels_3d.
                    let dist = rng.gen::<f64>().cbrt() * tumor_radius * 0.95;
                    let theta = rng.gen::<f64>() * std::f64::consts::TAU;
                    let cos_phi = 2.0 * rng.gen::<f64>() - 1.0;
                    let sin_phi = (1.0 - cos_phi * cos_phi).sqrt();
                    (
                        center.0 + dist * cos_phi,
                        center.1 + dist * sin_phi * theta.cos(),
                        center.2 + dist * sin_phi * theta.sin(),
                    )
                })
                .collect()
        }
    };
    // Niche radius in lattice units, squared (compare in cell coordinates).
    let niche_cells2 = (cfg.niche_radius_um / grid.cell_size_um).powi(2);
    (0..grid.cells.len())
        .map(|idx| {
            if !grid.cells[idx].is_tumor {
                return false;
            }
            let (r, c, l) = grid.coords(idx);
            let (rf, cf, lf) = (r as f64, c as f64, l as f64);
            seed_points.iter().any(|&(sr, sc, sl)| {
                (rf - sr).powi(2) + (cf - sc).powi(2) + (lf - sl).powi(2) <= niche_cells2
            })
        })
        .collect()
}

// =====================================================================
// Multi-checkpoint immune brake (#264, immune realism Phase 3)
// =====================================================================

/// One immune-checkpoint axis (PD-1, CTLA-4, LAG-3, or TIM-3), modeled as an
/// independent brake on T-cell killing that its inhibitor drug can lift.
#[derive(Clone, Copy, Debug)]
pub struct Checkpoint {
    /// Intrinsic brake strength ∈ [0, 1] — the fraction of kill this checkpoint
    /// suppresses when fully engaged and not drug-blocked. `0.0` ⇒ inactive.
    pub brake: f64,
    /// Fraction of this checkpoint's brake removed by its inhibitor drug
    /// ∈ [0, 1] (e.g. `anti_pd1_efficacy` for PD-1). `0.0` ⇒ no drug.
    pub drug_efficacy: f64,
}

impl Checkpoint {
    /// An inactive axis (no brake, no drug) — contributes nothing.
    pub fn inactive() -> Self {
        Checkpoint {
            brake: 0.0,
            drug_efficacy: 0.0,
        }
    }

    /// Residual brake after drug blockade: `brake · (1 − drug_efficacy)`.
    #[inline]
    pub fn residual(&self) -> f64 {
        // Out-of-range inputs push `combined_brake` outside [0,1], which makes
        // `immune_kill_probability`'s `1 − effective_brake` negative or > 1 and
        // silently produces nonsense kill probabilities (no panic). Guard like
        // `exhaustion_factor` / `suppressor_kill_multiplier` do.
        debug_assert!(
            (0.0..=1.0).contains(&self.brake) && (0.0..=1.0).contains(&self.drug_efficacy),
            "checkpoint brake and drug_efficacy must be in [0, 1]; got brake={}, drug_efficacy={}",
            self.brake,
            self.drug_efficacy
        );
        self.brake * (1.0 - self.drug_efficacy)
    }
}

/// A panel of immune checkpoints (#264 Phase 3): PD-1, CTLA-4, LAG-3, TIM-3,
/// generalizing the single PD-1 brake. Each axis is an **independent** brake on
/// the kill, so the combined brake is `1 − Π(1 − residualᵢ)` — and a panel with
/// only PD-1 active reduces *exactly* to the single-PD-1 model
/// (`SpatialImmuneConfig::effective_brake`), so a consumer that doesn't opt into
/// the panel stays byte-identical. Models anti-PD-1 + anti-CTLA-4 combinations
/// (Sharma & Allison, Cell 2015). `Copy`, so `Overrides` stays `Copy`-friendly.
///
/// **When a consumer supplies a panel, it fully replaces** the single-brake
/// `SpatialImmuneConfig::{pd1_brake, anti_pd1_efficacy}` — set the PD-1 axis
/// here (`pd1_only` / `pd1_ctla4_tumor` + `with_anti_pd1`) rather than on the
/// immune config, whose brake fields are then ignored.
#[derive(Clone, Copy, Debug)]
pub struct CheckpointPanel {
    pub pd1: Checkpoint,
    pub ctla4: Checkpoint,
    pub lag3: Checkpoint,
    pub tim3: Checkpoint,
}

impl CheckpointPanel {
    /// PD-1 only (CTLA-4/LAG-3/TIM-3 inactive). `combined_brake` then equals the
    /// single-PD-1 `pd1_brake · (1 − anti_pd1)` — the byte-identical baseline.
    pub fn pd1_only(pd1_brake: f64, anti_pd1: f64) -> Self {
        CheckpointPanel {
            pd1: Checkpoint {
                brake: pd1_brake,
                drug_efficacy: anti_pd1,
            },
            ctla4: Checkpoint::inactive(),
            lag3: Checkpoint::inactive(),
            tim3: Checkpoint::inactive(),
        }
    }

    /// A tumor expressing both PD-1 and CTLA-4 brakes, no drug yet — the
    /// substrate for the anti-PD-1 vs anti-PD-1+anti-CTLA-4 comparison. Apply
    /// drugs with [`with_anti_pd1`](Self::with_anti_pd1) /
    /// [`with_anti_ctla4`](Self::with_anti_ctla4). Placeholders pending
    /// calibration.
    pub fn pd1_ctla4_tumor() -> Self {
        CheckpointPanel {
            pd1: Checkpoint {
                brake: 0.7,
                drug_efficacy: 0.0,
            },
            ctla4: Checkpoint {
                brake: 0.5,
                drug_efficacy: 0.0,
            },
            lag3: Checkpoint::inactive(),
            tim3: Checkpoint::inactive(),
        }
    }

    /// Set the anti-PD-1 drug efficacy (fraction of the PD-1 brake removed).
    pub fn with_anti_pd1(mut self, efficacy: f64) -> Self {
        self.pd1.drug_efficacy = efficacy;
        self
    }

    /// Set the anti-CTLA-4 drug efficacy.
    pub fn with_anti_ctla4(mut self, efficacy: f64) -> Self {
        self.ctla4.drug_efficacy = efficacy;
        self
    }

    /// Set the anti-LAG-3 drug efficacy.
    pub fn with_anti_lag3(mut self, efficacy: f64) -> Self {
        self.lag3.drug_efficacy = efficacy;
        self
    }

    /// Set the anti-TIM-3 drug efficacy.
    pub fn with_anti_tim3(mut self, efficacy: f64) -> Self {
        self.tim3.drug_efficacy = efficacy;
        self
    }

    /// Combined brake from all axes acting independently:
    /// `1 − Π(1 − brakeᵢ·(1−drug_efficacyᵢ))`, in `[0, 1]`. Feeds
    /// [`immune_kill_probability`] in place of the single-PD-1
    /// `effective_brake`. Inactive axes (brake 0) contribute a factor of 1.
    ///
    /// **Assumes independence** between checkpoints (no shared/redundant
    /// signaling), so the multiplicative composition is an *upper bound* on the
    /// true combined brake — real checkpoints partially overlap, which would
    /// make dual blockade slightly less additive than modeled. An uncalibrated
    /// first-order choice; treat the combined magnitude qualitatively.
    #[inline]
    #[must_use]
    pub fn combined_brake(&self) -> f64 {
        let pass: f64 = [self.pd1, self.ctla4, self.lag3, self.tim3]
            .iter()
            .map(|c| 1.0 - c.residual())
            .product();
        1.0 - pass
    }
}

/// Dendritic-cell subset composition (#264 Phase 4).
///
/// The anti-tumor efficacy of the DAMP -> T-cell priming step depends on WHICH
/// DC subset takes up the DAMPs. **cDC1** (Batf3-dependent, cross-presenting)
/// are the rare but critical drivers of anti-tumor CD8 immunity, whereas
/// **cDC2** skew toward Th17/CD4 help and are far less effective at priming
/// tumor killing. A cDC1-poor tumor therefore converts the same DAMP signal
/// into LESS immune kill (Broz et al., Cancer Cell 2014, PMID 25446897: rare
/// activating cDC1 are critical for T-cell immunity).
///
/// The layer collapses to one uniform anti-tumor **priming-efficiency** scalar
/// ([`priming_efficiency`](Self::priming_efficiency)) that a consumer multiplies
/// into the immune kill probability (the `dc_activation` x
/// `immune_kill_probability` product). [`balanced`](Self::balanced) gives
/// efficiency `1.0` (identity, byte-identical); a cDC1-poor
/// [`literature`](Self::literature) config gives efficiency `< 1.0` (reduced
/// immune kill). The subset fractions/efficiencies are UNCALIBRATED placeholders;
/// the direction (cDC1-poor => weaker anti-tumor priming) is the claim, not the
/// number.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DcSubsetConfig {
    /// Fraction of tumor DCs that are cDC1 (cross-presenting); the rest are cDC2.
    pub cdc1_fraction: f64,
    /// Anti-tumor priming efficiency contributed by cDC1 (cross-presentation).
    pub cdc1_efficiency: f64,
    /// Anti-tumor priming efficiency contributed by cDC2 (Th17-skewing, weaker).
    pub cdc2_efficiency: f64,
}

impl DcSubsetConfig {
    /// Identity: both subsets fully effective => priming efficiency `1.0` =>
    /// the consumer's multiplier is exactly `1.0` => byte-identical.
    pub fn balanced() -> Self {
        DcSubsetConfig {
            cdc1_fraction: 1.0,
            cdc1_efficiency: 1.0,
            cdc2_efficiency: 1.0,
        }
    }

    /// Literature-informed cDC1-poor tumor (placeholders pending calibration):
    /// cDC1 are rare (~10% of tumor DCs) and cross-present efficiently, while the
    /// dominant cDC2 prime anti-tumor CD8 killing far less well. Refs: Broz et
    /// al., Cancer Cell 2014 (PMID 25446897). The fractions/efficiencies are
    /// UNCALIBRATED; only the direction (cDC1-poor => weaker priming) is claimed.
    /// Note: `cdc1_efficiency` is held at 1.0 here, so this shipped default
    /// exercises only the `cdc1_fraction` and `cdc2_efficiency` knobs (the cDC1
    /// efficacy knob exists for callers modeling impaired cDC1 cross-presentation).
    pub fn literature() -> Self {
        DcSubsetConfig {
            cdc1_fraction: 0.1,
            cdc1_efficiency: 1.0,
            cdc2_efficiency: 0.3,
        }
    }

    /// Uniform anti-tumor priming-efficiency scalar (a subset-weighted average):
    /// `cdc1_fraction · cdc1_efficiency + (1 − cdc1_fraction) · cdc2_efficiency`.
    /// A consumer multiplies this into the immune kill probability.
    #[must_use]
    pub fn priming_efficiency(&self) -> f64 {
        self.cdc1_fraction * self.cdc1_efficiency
            + (1.0 - self.cdc1_fraction) * self.cdc2_efficiency
    }

    /// True when the mix has no effect (priming efficiency exactly `1.0`), so the
    /// consumer stays byte-identical.
    ///
    /// The exact `== 1.0` compare is **intentional**: only `balanced()` (and the
    /// degenerate all-cDC1-at-full-efficiency case) reaches exactly 1.0, and that
    /// is the only identity config the matrix/snapshot path constructs. Do NOT
    /// relax this to an epsilon compare: a near-1.0 placeholder must route through
    /// the realism kill path (where its sub-1.0 scalar is applied), not slip onto
    /// the default byte-identical path.
    pub fn is_identity(&self) -> bool {
        self.priming_efficiency() == 1.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::grid::TumorGrid;

    // =============================
    // diffuse_damp_3d_step tests
    // =============================

    /// **v2 addition**: stability `assert` (fires in BOTH debug and
    /// release — silent-failure class) rejects sim-tme's 2D default
    /// (0.08), which is unsafe in 3D (0.08 × 26 = 2.08 > 1). Critical
    /// bug-class guard; per-call cost is one multiply + one compare.
    #[test]
    #[should_panic(expected = "diffusion_fraction × 26")]
    fn diffusion_fraction_stability_assertion() {
        let g = TumorGrid3D::generate(5, 5, 5, 20.0, 42);
        let mut field = vec![0.0; g.cells.len()];
        let mut scratch = vec![0.0; g.cells.len()];
        diffuse_damp_3d_step(&mut field, &mut scratch, &g, 0.08, 0.03);
    }

    /// Field length mismatch panics in release (matches stromal pattern).
    #[test]
    #[should_panic(expected = "damp_field.len()")]
    fn field_length_validation() {
        let g = TumorGrid3D::generate(5, 5, 5, 20.0, 42);
        let mut field = vec![0.0; g.cells.len() / 2];
        let mut scratch = vec![0.0; g.cells.len()];
        diffuse_damp_3d_step(&mut field, &mut scratch, &g, 0.025, 0.03);
    }

    /// Scratch length mismatch panics in release.
    #[test]
    #[should_panic(expected = "scratch.len()")]
    fn scratch_length_validation() {
        let g = TumorGrid3D::generate(5, 5, 5, 20.0, 42);
        let mut field = vec![0.0; g.cells.len()];
        let mut scratch = vec![0.0; g.cells.len() / 2];
        diffuse_damp_3d_step(&mut field, &mut scratch, &g, 0.025, 0.03);
    }

    /// A single DAMP source at an INTERIOR cell spreads to all 26
    /// neighbors after one step. With no clearance, each neighbor
    /// receives exactly `share = local × diffusion_fraction`, and the
    /// source retains `local × (1 - 26 × diffusion_fraction)`.
    #[test]
    fn interior_source_spreads_to_26_neighbors() {
        let g = TumorGrid3D::generate(10, 10, 10, 20.0, 42);
        let n = g.cells.len();
        let mut field = vec![0.0; n];
        let mut scratch = vec![0.0; n];

        let center = g.flat_index(5, 5, 5); // interior cell
        let source_value = 100.0_f64;
        let fraction = 0.025_f64;
        field[center] = source_value;

        diffuse_damp_3d_step(&mut field, &mut scratch, &g, fraction, 0.0);

        // Source retains local × (1 - 26 × fraction) = 100 × (1 - 0.65) = 35
        let expected_source = source_value * (1.0 - 26.0 * fraction);
        assert!(
            (field[center] - expected_source).abs() < 1e-9,
            "source should retain {expected_source}, got {}",
            field[center]
        );

        // Each of the 26 neighbors should have received exactly share = 2.5
        let expected_share = source_value * fraction;
        let (neighbors, count) = g.neighbors(5, 5, 5);
        assert_eq!(count, 26);
        for &(nr, nc, nl) in &neighbors[..26] {
            let nidx = g.flat_index(nr, nc, nl);
            assert!(
                (field[nidx] - expected_share).abs() < 1e-9,
                "neighbor ({nr},{nc},{nl}) should have {expected_share}, got {}",
                field[nidx]
            );
        }
    }

    /// A source on a grid FACE (one coordinate on boundary) has 17
    /// neighbors. Verifies the in-between boundary case between corner
    /// (7) and interior (26) — closes the reviewer's flagged coverage
    /// gap on partial-boundary cells. Mass conservation still holds:
    /// source loses `share × 17`, 17 neighbors gain `share` each.
    #[test]
    fn face_source_spreads_to_17_neighbors() {
        let g = TumorGrid3D::generate(10, 10, 10, 20.0, 42);
        let n = g.cells.len();
        let mut field = vec![0.0; n];
        let mut scratch = vec![0.0; n];

        // Face cell: r=0 (boundary), c=5 (interior), l=5 (interior). Per
        // grid::tests_3d::neighbor_counts_at_boundary_types this has 17.
        let face = g.flat_index(0, 5, 5);
        let source_value = 100.0_f64;
        let fraction = 0.025_f64;
        field[face] = source_value;

        diffuse_damp_3d_step(&mut field, &mut scratch, &g, fraction, 0.0);

        // Source retains 100 × (1 - 17 × 0.025) = 100 × 0.575 = 57.5
        let expected_source = source_value * (1.0 - 17.0 * fraction);
        assert!(
            (field[face] - expected_source).abs() < 1e-9,
            "face source should retain {expected_source}, got {}",
            field[face]
        );

        let (_neighbors, count) = g.neighbors(0, 5, 5);
        assert_eq!(
            count, 17,
            "face cell (one coord on boundary) should have exactly 17 Moore neighbors"
        );
    }

    /// A source at the grid CORNER has only 7 neighbors. Source loses
    /// `share × 7` (not `share × 26`); each of the 7 neighbors gains
    /// `share`. Mass is still conserved.
    #[test]
    fn corner_source_spreads_to_7_neighbors() {
        let g = TumorGrid3D::generate(10, 10, 10, 20.0, 42);
        let n = g.cells.len();
        let mut field = vec![0.0; n];
        let mut scratch = vec![0.0; n];

        let corner = g.flat_index(0, 0, 0);
        let source_value = 100.0_f64;
        let fraction = 0.025_f64;
        field[corner] = source_value;

        diffuse_damp_3d_step(&mut field, &mut scratch, &g, fraction, 0.0);

        // Source retains 100 × (1 - 7 × 0.025) = 100 × 0.825 = 82.5
        let expected_source = source_value * (1.0 - 7.0 * fraction);
        assert!(
            (field[corner] - expected_source).abs() < 1e-9,
            "corner source should retain {expected_source}, got {}",
            field[corner]
        );

        let (_neighbors, count) = g.neighbors(0, 0, 0);
        assert_eq!(
            count, 7,
            "grid corner should have exactly 7 Moore neighbors"
        );
    }

    /// Total DAMP is conserved (modulo clearance) across a diffusion step.
    /// `sum_after = sum_before × (1 - clearance_rate)` to numerical
    /// tolerance.
    #[test]
    fn total_damp_conserved_modulo_clearance() {
        let g = TumorGrid3D::generate(10, 10, 10, 20.0, 42);
        let n = g.cells.len();
        let mut field = vec![0.0; n];
        let mut scratch = vec![0.0; n];

        // Sparse, varied initial conditions.
        field[g.flat_index(5, 5, 5)] = 100.0;
        field[g.flat_index(3, 7, 2)] = 50.0;
        field[g.flat_index(8, 1, 6)] = 25.0;

        let sum_before: f64 = field.iter().sum();
        let clearance = 0.03_f64;
        diffuse_damp_3d_step(&mut field, &mut scratch, &g, 0.025, clearance);
        let sum_after: f64 = field.iter().sum();

        let expected = sum_before * (1.0 - clearance);
        let rel_error = (sum_after - expected).abs() / expected;
        assert!(
            rel_error < 1e-12,
            "mass not conserved: before={sum_before}, after={sum_after}, expected≈{expected} (rel error {rel_error:.2e})"
        );
    }

    /// Zero diffusion → only clearance acts → every cell decays by
    /// `(1 - clearance)` exactly.
    #[test]
    fn zero_diffusion_just_applies_clearance() {
        let g = TumorGrid3D::generate(10, 10, 10, 20.0, 42);
        let n = g.cells.len();
        let mut field = vec![1.0_f64; n]; // uniform
        let mut scratch = vec![0.0; n];

        diffuse_damp_3d_step(&mut field, &mut scratch, &g, 0.0, 0.03);

        let expected = 0.97_f64;
        for (i, &v) in field.iter().enumerate() {
            assert!(
                (v - expected).abs() < 1e-12,
                "cell {i} should be {expected}, got {v}"
            );
        }
    }

    /// Determinism: same inputs → same outputs.
    #[test]
    fn deterministic_same_inputs_same_outputs() {
        let g = TumorGrid3D::generate(10, 10, 10, 20.0, 42);
        let n = g.cells.len();
        let init: Vec<f64> = (0..n).map(|i| (i as f64 % 7.0) * 0.5).collect();

        let mut field1 = init.clone();
        let mut scratch1 = vec![0.0; n];
        diffuse_damp_3d_step(&mut field1, &mut scratch1, &g, 0.025, 0.03);

        let mut field2 = init.clone();
        let mut scratch2 = vec![0.0; n];
        diffuse_damp_3d_step(&mut field2, &mut scratch2, &g, 0.025, 0.03);

        assert_eq!(field1, field2);
    }

    /// **Cross-geometry** (AC #4, scope-bounded): at matched
    /// `diffusion_fraction` and an interior source, a 3D source decays
    /// **faster** than a 2D source per step because 26 > 8 neighbors take
    /// share. This is the geometric component of the issue's kill-ratio
    /// question; the full kill-ratio answer (whether SDT's 104:1
    /// differential holds in 3D) requires multi-step simulation and lands
    /// with #195/#196.
    ///
    /// Uses 3D-safe `0.025` so the test exercises real diffusion math
    /// rather than mass-destroying via the stability violation.
    #[test]
    fn three_d_source_decays_faster_than_two_d_at_matched_diffusion() {
        let fraction = 0.025_f64;
        let source = 100.0_f64;

        // 2D: 10×10 interior cell, 8 neighbors. Inline sim-tme math
        // (source of truth: sim-tme/src/main.rs:701-716).
        let g2 = TumorGrid::generate(10, 10, 20.0, 42);
        let n2 = g2.cells.len();
        let mut field2 = vec![0.0; n2];
        let center_2d = 5 * g2.cols + 5;
        field2[center_2d] = source;
        let (_, count_2d) = g2.neighbors(5, 5);
        let source_after_2d = source * (1.0 - count_2d as f64 * fraction);

        // 3D: 10×10×10 interior cell, 26 neighbors.
        let g3 = TumorGrid3D::generate(10, 10, 10, 20.0, 42);
        let n3 = g3.cells.len();
        let mut field3 = vec![0.0; n3];
        let mut scratch3 = vec![0.0; n3];
        let center_3d = g3.flat_index(5, 5, 5);
        field3[center_3d] = source;
        diffuse_damp_3d_step(&mut field3, &mut scratch3, &g3, fraction, 0.0);

        assert!(
            field3[center_3d] < source_after_2d,
            "3D source ({}) should retain LESS than 2D source ({}) at matched diffusion — \
             geometry: 26 neighbors take share vs 8. Empirical: 3D loses 26·share = {}, \
             2D loses 8·share = {}.",
            field3[center_3d],
            source_after_2d,
            26.0 * source * fraction,
            8.0 * source * fraction
        );

        // Quantitative: 3D should lose ~26/8 = 3.25× more per step.
        let loss_3d = source - field3[center_3d];
        let loss_2d = source - source_after_2d;
        let ratio = loss_3d / loss_2d;
        let expected_ratio = 26.0 / 8.0; // = 3.25
        assert!(
            (ratio - expected_ratio).abs() < 0.01,
            "loss ratio 3D/2D should be ≈ {expected_ratio} (26/8), got {ratio}"
        );
    }

    // =============================
    // dc_activation tests
    // =============================

    /// `dc_activation(0, kd) = 0` exactly (IEEE: 0/kd = 0 for kd > 0).
    #[test]
    fn dc_activation_at_zero_damp_returns_zero() {
        for &kd in &[1.0, 50.0, 1000.0] {
            assert_eq!(dc_activation(0.0, kd), 0.0);
        }
    }

    /// `dc_activation(kd, kd) = 0.5` exactly (IEEE: x/(2x) = 0.5 for
    /// finite x; even though `kd` may not be IEEE-exact, the doubling
    /// `kd + kd = 2*kd` is exact, then x/(2x) when 2x is finite is
    /// IEEE-exactly 0.5).
    #[test]
    fn dc_activation_at_damp_eq_kd_returns_half() {
        for &kd in &[1.0_f64, 50.0, 1000.0, 7.4] {
            assert_eq!(dc_activation(kd, kd), 0.5);
        }
    }

    /// `dc_activation(damp >> kd, kd) ≈ 1.0` asymptote.
    #[test]
    fn dc_activation_large_damp_approaches_one() {
        let kd = 50.0;
        assert!(dc_activation(1e9, kd) > 0.999999);
        assert!(dc_activation(1e9, kd) <= 1.0);
    }

    /// Monotone non-decreasing: more DAMP → higher activation.
    #[test]
    fn dc_activation_monotone_non_decreasing() {
        let kd = 50.0;
        let probes = [0.0, 1.0, 10.0, 50.0, 100.0, 1000.0];
        let mut prev = dc_activation(probes[0], kd);
        for &d in &probes[1..] {
            let cur = dc_activation(d, kd);
            assert!(
                cur >= prev,
                "activation decreased from {prev} to {cur} at damp={d}"
            );
            prev = cur;
        }
    }

    // =============================
    // immune_kill_probability tests
    // =============================

    /// `min(0.99)` cap: even at activation=1, rate=1, brake=0 → result
    /// is exactly 0.99. Hardcoded cap matches sim-tme.
    #[test]
    fn immune_kill_clamps_at_0_99_for_extreme_inputs() {
        let result = immune_kill_probability(1.0, 1.0, 0.0);
        assert_eq!(result, 0.99);

        // Even more extreme: activation > 1 (shouldn't happen but test the
        // clamp behavior).
        let result2 = immune_kill_probability(2.0, 1.0, 0.0);
        assert_eq!(result2, 0.99);
    }

    /// `effective_brake = 1.0` → full PD-1 suppression → zero kill.
    /// IEEE-exact: `(activation × rate × 0.0) = 0.0` for any finite
    /// activation, rate.
    #[test]
    fn immune_kill_full_brake_returns_zero() {
        for &activation in &[0.0_f64, 0.5, 1.0] {
            for &rate in &[0.01_f64, 0.02, 0.05] {
                assert_eq!(immune_kill_probability(activation, rate, 1.0), 0.0);
            }
        }
    }

    /// **Reviewer-flagged correctness regression guard**: NaN inputs
    /// must propagate to NaN output, not silently cap at 0.99.
    ///
    /// Background: Rust's `f64::min` treats NaN as the maximum, so
    /// `NaN.min(0.99) = 0.99` (not NaN). A naïve `(raw).min(0.99)` cap
    /// would turn NaN inputs (e.g., from `dc_activation(0.0, 0.0) = 0/0
    /// = NaN`) into a near-certain kill probability — silent
    /// correctness bug. The implementation uses an explicit
    /// `if raw > 0.99` branch instead (NaN comparisons return false, so
    /// the `else` branch returns the NaN unchanged).
    ///
    /// This test locks down the NaN-propagation contract so a future
    /// refactor back to `.min(0.99)` would fail loudly.
    #[test]
    fn immune_kill_propagates_nan() {
        // Direct NaN injection on each argument.
        assert!(immune_kill_probability(f64::NAN, 0.02, 0.21).is_nan());
        assert!(immune_kill_probability(0.5, f64::NAN, 0.21).is_nan());
        assert!(immune_kill_probability(0.5, 0.02, f64::NAN).is_nan());

        // Realistic scenario: dc_activation(0, 0) = 0/0 = NaN propagating
        // through the kill probability.
        let bad_activation = dc_activation(0.0, 0.0);
        assert!(bad_activation.is_nan(), "test precondition");
        let kill_prob = immune_kill_probability(bad_activation, 0.02, 0.21);
        assert!(
            kill_prob.is_nan(),
            "NaN from upstream should propagate to NaN kill probability, \
             not silently cap at 0.99. Got {kill_prob} (a regression to `.min(0.99)` \
             would yield 0.99 — a near-certain kill from invalid input)."
        );
    }

    /// Default sim-tme numerical example: activation=0.5, rate=0.02,
    /// brake=0.21 (= pd1_brake 0.7 × (1 - anti_pd1 0.7)) → ≈ 0.0079.
    /// Tight tolerance via libm.
    #[test]
    fn immune_kill_at_sim_tme_defaults_matches_expected() {
        let activation = 0.5;
        let rate = 0.02;
        let effective_brake = 0.7 * (1.0 - 0.7); // = 0.21
        let prob = immune_kill_probability(activation, rate, effective_brake);
        // Expected: 0.5 × 0.02 × (1 - 0.21) = 0.5 × 0.02 × 0.79 = 0.0079
        let expected = 0.5 * 0.02 * 0.79;
        assert!(
            (prob - expected).abs() < 1e-9,
            "expected ≈ {expected}, got {prob}"
        );
        assert!(prob < 0.99, "should not be capped at this activation");
    }

    // =============================
    // exhaustion_factor tests (#243)
    // =============================

    #[test]
    fn exhaustion_is_identity_at_zero_rate() {
        // rate 0 => no suppression for any cumulative count (byte-identity).
        for cum in [0u32, 1, 10, 1000] {
            assert_eq!(exhaustion_factor(cum, 0.0), 1.0);
        }
    }

    #[test]
    fn exhaustion_is_identity_at_zero_cumulative() {
        assert_eq!(exhaustion_factor(0, 0.5), 1.0);
    }

    #[test]
    fn exhaustion_decreases_with_cumulative_kills() {
        let rate = 0.1;
        let f1 = exhaustion_factor(1, rate);
        let f10 = exhaustion_factor(10, rate);
        assert!(f1 < 1.0 && f10 < f1);
        // Bounded in (0, 1].
        assert!(f10 > 0.0);
        // Exact hyperbolic value: 1/(1 + 0.1*10) = 0.5.
        assert!((exhaustion_factor(10, 0.1) - 0.5).abs() < 1e-12);
    }

    #[test]
    fn exhaustion_decreases_with_rate() {
        assert!(exhaustion_factor(5, 0.2) < exhaustion_factor(5, 0.05));
    }

    // ===== Suppressor field (#264 Phase 2) =====

    #[test]
    fn suppressor_multiplier_identity_when_disabled() {
        // strength 0 ⇒ identity regardless of field.
        assert_eq!(suppressor_kill_multiplier(0.9, 0.0), 1.0);
        // field 0 ⇒ identity regardless of strength.
        assert_eq!(suppressor_kill_multiplier(0.0, 6.0), 1.0);
        // both > 0 ⇒ suppression < 1, monotonic in each argument.
        assert!(suppressor_kill_multiplier(0.5, 6.0) < 1.0);
        assert!(suppressor_kill_multiplier(0.8, 6.0) < suppressor_kill_multiplier(0.2, 6.0));
        assert!(suppressor_kill_multiplier(0.5, 10.0) < suppressor_kill_multiplier(0.5, 2.0));
        // bounded in (0, 1].
        assert!(suppressor_kill_multiplier(1.0, 6.0) > 0.0);
    }

    // ===== Immunosuppressive ferroptosis (#337) =====

    #[test]
    fn ferroptotic_immunosuppression_identity_at_zero_and_suppresses() {
        // strength 0 ⇒ identity regardless of local DAMP (byte-identical).
        for damp in [0.0, 1.0, 50.0, 1000.0] {
            assert_eq!(ferroptotic_immunosuppression(damp, 0.0), 1.0);
        }
        // local DAMP 0 ⇒ identity regardless of strength.
        assert_eq!(ferroptotic_immunosuppression(0.0, 0.5), 1.0);
        // both > 0 ⇒ suppression < 1, monotone decreasing in each argument.
        assert!(ferroptotic_immunosuppression(10.0, 0.1) < 1.0);
        assert!(
            ferroptotic_immunosuppression(20.0, 0.1) < ferroptotic_immunosuppression(10.0, 0.1)
        );
        assert!(
            ferroptotic_immunosuppression(10.0, 0.2) < ferroptotic_immunosuppression(10.0, 0.1)
        );
        // Exact hyperbolic value: 1/(1 + 0.1*10) = 0.5.
        assert!((ferroptotic_immunosuppression(10.0, 0.1) - 0.5).abs() < 1e-12);
        // bounded in (0, 1].
        assert!(ferroptotic_immunosuppression(1.0e9, 1.0) > 0.0);
    }

    #[test]
    fn sasp_field_multiplier_is_signed_and_identity_at_zero() {
        // strength 0 ⇒ identity regardless of field (byte-identical).
        for field in [0.0, 1.0, 50.0, 1000.0] {
            assert_eq!(sasp_field_kill_multiplier(field, 0.0), 1.0);
        }
        // field 0 ⇒ identity regardless of strength (a cell out of SASP reach is
        // never affected, so an enabled-but-empty field stays byte-identical).
        for strength in [-2.0, -0.5, 0.5, 2.0] {
            assert_eq!(sasp_field_kill_multiplier(0.0, strength), 1.0);
        }
        // Immunosuppressive arm (strength > 0): lowers kill, hyperbolic, in (0,1).
        assert!(sasp_field_kill_multiplier(10.0, 0.1) < 1.0);
        assert!((sasp_field_kill_multiplier(10.0, 0.1) - 0.5).abs() < 1e-12); // 1/(1+1)
        assert!(
            sasp_field_kill_multiplier(20.0, 0.1) < sasp_field_kill_multiplier(10.0, 0.1),
            "monotone decreasing in field on the suppressive arm"
        );
        assert!(sasp_field_kill_multiplier(1.0e9, 1.0) > 0.0, "bounded > 0");
        // Surveillance arm (strength < 0): raises kill, linear, >= 1.
        assert!(sasp_field_kill_multiplier(10.0, -0.1) > 1.0);
        assert!((sasp_field_kill_multiplier(10.0, -0.1) - 2.0).abs() < 1e-12); // 1 + 0.1*10
        assert!(
            sasp_field_kill_multiplier(20.0, -0.1) > sasp_field_kill_multiplier(10.0, -0.1),
            "monotone increasing in field on the surveillance arm"
        );
        // The two arms straddle 1.0 at the same |strength|·field.
        assert!(
            sasp_field_kill_multiplier(10.0, 0.1) < 1.0
                && sasp_field_kill_multiplier(10.0, -0.1) > 1.0
        );
    }

    /// #376: the SASP field reaches a NON-senescent NEIGHBOR. One senescent
    /// source cell seeds the field; after the exact diffusion operator the
    /// consumer (sim-tme-3d) uses, an adjacent non-senescent cell sees a positive
    /// field, so its `sasp_field_kill_multiplier` shifts off 1.0 — the bystander
    /// coupling the cell-autonomous #341 multiplier could not express. The
    /// immunosuppressive arm LOWERS the neighbor's kill, the surveillance arm
    /// RAISES it.
    #[test]
    fn sasp_field_shifts_a_non_senescent_neighbor_kill() {
        let g = TumorGrid3D::generate(10, 10, 10, 20.0, 42);
        let n = g.cells.len();
        let mut field = vec![0.0; n];
        let mut scratch = vec![0.0; n];

        // One senescent SOURCE cell at an interior site; a Moore neighbor of it
        // is NON-senescent (it is never a source).
        let src = g.flat_index(5, 5, 5);
        let nbr = g.flat_index(5, 5, 6); // face-adjacent, distinct cell

        // Seed + diffuse a few steps exactly as the consumer does.
        let replenish = 0.15_f64;
        for _ in 0..5 {
            field[src] = (field[src] + replenish).min(1.0);
            diffuse_damp_3d_step(&mut field, &mut scratch, &g, 0.025, 0.03);
        }

        // The signal reached the non-senescent neighbor.
        assert!(
            field[nbr] > 0.0,
            "SASP field should diffuse to the non-senescent neighbor; got {}",
            field[nbr]
        );

        // Its kill multiplier shifts off 1.0, signed by the arm.
        let suppress = sasp_field_kill_multiplier(field[nbr], 0.5);
        let surveil = sasp_field_kill_multiplier(field[nbr], -0.5);
        assert!(
            suppress < 1.0,
            "immunosuppressive SASP should LOWER a neighbor's kill; got {suppress}"
        );
        assert!(
            surveil > 1.0,
            "surveillance SASP should RAISE a neighbor's kill; got {surveil}"
        );
    }

    #[test]
    fn suppressor_sources_heuristic_are_tumor_only_and_deterministic() {
        let g = TumorGrid3D::generate(30, 30, 30, 20.0, 42);
        let cfg = SuppressorConfig::enabled();
        let a = suppressor_source_mask_3d(&g, &cfg, None, 123);
        let b = suppressor_source_mask_3d(&g, &cfg, None, 123);
        assert_eq!(a, b, "deterministic given the same seed");
        let n_sources = a.iter().filter(|&&s| s).count();
        assert!(n_sources > 0, "heuristic seeding marks some sources");
        // Only tumor cells are ever sources.
        for (idx, &is_src) in a.iter().enumerate() {
            if is_src {
                assert!(g.cells[idx].is_tumor, "a source must be a tumor cell");
            }
        }
        // A different seed gives a different layout (sources move).
        let c = suppressor_source_mask_3d(&g, &cfg, None, 999);
        assert_ne!(a, c, "different seed ⇒ different niches");
    }

    #[test]
    fn suppressor_sources_perivascular_cluster_near_vessels() {
        let g = TumorGrid3D::generate(30, 30, 30, 20.0, 42);
        let cfg = SuppressorConfig::enabled();
        // One vessel at the center ⇒ sources are the central niche only.
        let vessels = vec![(15.0, 15.0, 15.0)];
        let mask = suppressor_source_mask_3d(&g, &cfg, Some(&vessels), 0);
        let niche_cells = cfg.niche_radius_um / g.cell_size_um; // 60/20 = 3 cells
        for (idx, &is_src) in mask.iter().enumerate() {
            if is_src {
                let (r, c, l) = g.coords(idx);
                let d = (((r as f64 - 15.0).powi(2)
                    + (c as f64 - 15.0).powi(2)
                    + (l as f64 - 15.0).powi(2))
                .sqrt())
                .round();
                assert!(
                    d <= niche_cells + 1.0,
                    "perivascular source at lattice dist {d} exceeds niche {niche_cells}"
                );
            }
        }
        assert!(
            mask.iter().any(|&s| s),
            "the central vessel seeds a non-empty niche"
        );
    }

    // ===== Multi-checkpoint brake (#264 Phase 3) =====

    #[test]
    fn pd1_only_panel_matches_the_single_pd1_brake() {
        // A PD-1-only panel must reproduce `effective_brake = pd1·(1−anti_pd1)`
        // exactly — the byte-identity equivalence to the single-brake model.
        assert!((CheckpointPanel::pd1_only(0.7, 0.0).combined_brake() - 0.7).abs() < 1e-12);
        let blocked = CheckpointPanel::pd1_only(0.7, 0.8).combined_brake();
        assert!((blocked - 0.7 * 0.2).abs() < 1e-12, "got {blocked}");
        // Inactive extra axes never change the brake.
        assert_eq!(
            CheckpointPanel::pd1_only(0.5, 0.3).combined_brake(),
            CheckpointPanel {
                lag3: Checkpoint::inactive(),
                tim3: Checkpoint::inactive(),
                ..CheckpointPanel::pd1_only(0.5, 0.3)
            }
            .combined_brake()
        );
    }

    #[test]
    fn dual_blockade_lowers_brake_below_anti_pd1_alone() {
        // A PD-1 + CTLA-4 tumor: anti-PD-1 alone leaves CTLA-4 braking, so the
        // combined brake stays high; adding anti-CTLA-4 lifts both → lower brake
        // → more killing (the combination-immunotherapy result, Sharma & Allison
        // 2015). Lower brake ⇒ higher kill in `immune_kill_probability`.
        let mono = CheckpointPanel::pd1_ctla4_tumor()
            .with_anti_pd1(0.8)
            .combined_brake();
        let combo = CheckpointPanel::pd1_ctla4_tumor()
            .with_anti_pd1(0.8)
            .with_anti_ctla4(0.8)
            .combined_brake();
        let untreated = CheckpointPanel::pd1_ctla4_tumor().combined_brake();
        assert!(
            combo < mono && mono < untreated,
            "dual blockade < anti-PD-1 alone < untreated: combo={combo}, mono={mono}, untreated={untreated}"
        );
        // All brakes bounded in [0, 1].
        for b in [combo, mono, untreated] {
            assert!((0.0..=1.0).contains(&b), "brake {b} out of [0,1]");
        }
    }

    #[test]
    fn combined_brake_is_one_minus_product_of_passes() {
        // Two equal independent brakes of 0.5 combine to 1 − (0.5·0.5) = 0.75,
        // NOT 1.0 — independent brakes compose multiplicatively on the pass.
        let p = CheckpointPanel {
            pd1: Checkpoint {
                brake: 0.5,
                drug_efficacy: 0.0,
            },
            ctla4: Checkpoint {
                brake: 0.5,
                drug_efficacy: 0.0,
            },
            lag3: Checkpoint::inactive(),
            tim3: Checkpoint::inactive(),
        };
        assert!((p.combined_brake() - 0.75).abs() < 1e-12);
    }

    #[test]
    fn dc_subsets_balanced_is_identity() {
        let b = DcSubsetConfig::balanced();
        assert!(b.is_identity());
        assert!((b.priming_efficiency() - 1.0).abs() < 1e-12);
        // A cDC1-poor literature mix is NOT identity and primes < 1.0.
        let lit = DcSubsetConfig::literature();
        assert!(!lit.is_identity());
        assert!(lit.priming_efficiency() < 1.0);
        assert!(lit.priming_efficiency() > 0.0);
    }

    #[test]
    fn dc_subsets_priming_is_the_subset_weighted_average() {
        // 0.1·1.0 + 0.9·0.3 = 0.37 (the literature cDC1-poor default).
        let lit = DcSubsetConfig::literature();
        assert!((lit.priming_efficiency() - 0.37).abs() < 1e-12);
        // More cDC1 (the effective subset) ⇒ higher priming efficiency.
        let cdc1_rich = DcSubsetConfig {
            cdc1_fraction: 0.8,
            ..DcSubsetConfig::literature()
        };
        assert!(cdc1_rich.priming_efficiency() > lit.priming_efficiency());
        // All-cDC1 with full efficiency is the identity (efficiency 1.0).
        let all_cdc1 = DcSubsetConfig {
            cdc1_fraction: 1.0,
            ..DcSubsetConfig::literature()
        };
        assert!((all_cdc1.priming_efficiency() - 1.0).abs() < 1e-12);
        assert!(all_cdc1.is_identity());
    }
}
