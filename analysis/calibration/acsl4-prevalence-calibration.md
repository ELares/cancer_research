# ACSL4-status prevalence calibration (TCGA PanCancer Atlas, #462)

This anchors the ACSL4-status biomarker layer (#444, `ferroptosis-core` `acsl4`
module) to real per-cancer-type expression data. The #444 layer maps a per-tumor
ACSL4 expression *status* (1.0 = wild-type baseline, `< 1` = ferroptosis-refractory
via a collapsed PUFA substrate) to a PUFA-incorporation boost, but the per-cancer-type
prevalence of ACSL4-low/negative tumors was flagged **DATA-GATED**. This closes that
gap with the only login-free public source: the cBioPortal REST API over the 32 TCGA
PanCancer Atlas studies.

- **Fetch / artifact:** `scripts/fetch_acsl4_prevalence.py` (offline contract: run
  locally, hits the cBioPortal API; CI never fetches). Committed derived artifacts:
  `analysis/calibration/acsl4_prevalence_tcga.csv` (32 cancer types) and
  `acsl4_prevalence_tcga.json`.
- **Genes:** ACSL4 (2182), with GPX4 (2879) and SLC7A11 (23657) for context.
- **Anchor:** Doll et al., Nat Chem Biol 2017, PMID 27842070 (ACSL4 dictates
  ferroptosis sensitivity).

## What the data says

### 1. Within-cohort low-ACSL4 prevalence (the usable population prior)

Using the within-study mRNA z-scores (`*_rna_seq_v2_mrna_median_all_sample_Zscores`),
the fraction of tumors in each cancer type's low-ACSL4 tail:

| Threshold | Interpretation | Across 32 cancer types |
|-----------|----------------|------------------------|
| `z < -1` ("low", ~ACSL4_LOW) | refractory-leaning | min 10.8%, **median 14.4%**, max 18.8% |
| `z < -2` ("very low", ~ACSL4_NEGATIVE) | strongly refractory-leaning | min 0%, **median 3.0%**, max 5.4% |

This is the **real, committed population prior** for #444: when simulating a patient
cohort, roughly **1 in 7 tumors** falls in the low-ACSL4 (refractory-leaning) tail
and roughly **3%** in the very-low tail, fairly uniformly across cancer types. Because
the z-scores are computed within each study, this is a within-cohort stratification
number by construction (about a normal lower tail), and the small cross-type spread is
expected, not a defect.

### 2. The calibrated data→model bridge

The shipped status constants turn out to be exactly the integer-z points of a simple
linear bridge, now added as `acsl4::status_from_zscore(z) = max(0, 1 + z/2)`:

| ACSL4 mRNA z-score | status | constant |
|--------------------|--------|----------|
| `+1` | 1.5 | `ACSL4_HIGH` |
| `0`  | 1.0 | `ACSL4_NORMAL` (baseline) |
| `-1` | 0.5 | `ACSL4_LOW` |
| `-2` | 0.0 | `ACSL4_NEGATIVE` (PUFA-collapse floor) |

So a real patient ACSL4 z-score now maps directly onto the model's status scale via
`pufa_boost_from_status(status_from_zscore(z))`. The slope (`/2`) is the placeholder
that reproduces the existing constants; the **within-cohort z interpretation** is the
data-anchored part, not the absolute status→ferroptosis magnitude (that remains the
uncalibrated #444 linear placeholder).

## The honest negative result: bulk mRNA does NOT show HCC/AML as ACSL4-low

Doll 2017 reports HCC and some AML subtypes as constitutively ACSL4-low (hence
refractory). The cross-cancer-type **raw RSEM** ranking does **not** support this at
the bulk-mRNA-median level: of the 32 cancer types, **lihc (HCC) ranks highest** (top
percentile) and **laml (AML) ranks high** (~0.81 percentile) on raw ACSL4 RSEM. The
lowest-RSEM types are uvm, prad, brca, thym, ucs.

Two reasons this is not a contradiction, both reported rather than buried:

1. **Cross-study RSEM medians are batch-confounded.** Different TCGA studies were
   sequenced and normalized in different batches, so cross-study raw-expression
   comparison is unreliable for ranking cancer *types*. This ranking is therefore
   read qualitatively and used only to test the literature claim, which it fails.
2. **The ACSL4-low refractory phenotype is a protein/subtype/functional property**,
   not a bulk-mRNA-median one. Doll 2017's evidence is protein-level and cell-line
   level; a cancer type can have high bulk ACSL4 mRNA yet contain a refractory
   ACSL4-low subtype (HCC in fact has one of the larger within-cohort very-low tails,
   4.9%, despite the highest bulk median).

**Consequence for the model.** The ACSL4 status should be set from the **within-cohort
z-score** (which the bridge does), not from a cross-type bulk-mRNA ranking. The
per-cancer-type low-ACSL4 prevalence is a genuine prior; a clean per-cancer-type
"these whole tumor types are refractory" stratification is **not** supported by bulk
TCGA mRNA and would need protein (RPPA / IHC) or ACSL4-deletion data to establish. That
data gap is now specific and flagged, not a blanket "DATA-GATED".

## Status of the #444 layer after this leg

Moves from **uncalibrated / DATA-GATED** to **partially anchored**: the status→z-score
bridge and the per-cancer-type low-ACSL4 prevalence prior are anchored to TCGA; the
status→ferroptosis-magnitude mapping and the protein-level refractory-subtype
prevalence remain uncalibrated/data-gated. The remaining cell-line ACSL4-status-vs-
dose-response meta-analysis (#444) would calibrate the magnitude.
