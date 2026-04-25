# Weight Sensitivity Pre-Registration

Pre-registered BEFORE running any alternative weighting schemes.
Commit this document before implementing or running sensitivity analysis code.

## Current baseline formula

`weight = tier_weight × citation_modifier × recency_modifier`

- Tier weights: phase3-clinical=12.0, phase2-clinical=8.0, phase1-clinical=5.0, clinical-other=3.0, preclinical-invivo=2.0, preclinical-invitro=1.0, theoretical=0.5
- Citation modifier: 1.0 + (iCite_percentile / 200.0) → range 1.0–1.5×
- Recency modifier: 0.9 + ((year − 2015) / 11) × 0.2 → range 0.9–1.1×

## Target conclusions to test

1. Immunotherapy maintains weighted rank 1 under all schemes
2. Nanoparticle drops from article-count rank 2 to weighted rank 6 or lower
3. ADC and bispecific antibody rank in the weighted top 5
4. Bottom 5 mechanisms (HIFU, frequency-therapy, microbiome, electrolysis, and one other) are stable
5. Tier weights are the dominant driver of ranking differences (clinical maturity drives reordering)

## Stability definition

- **Stable**: conclusion holds under ALL 5 alternative schemes plus the baseline (6 total)
- **Directionally stable**: rank shift is ≤2 positions across all schemes
- **Fails**: any scheme reverses the conclusion

## Alternative weighting schemes

| Scheme | Tier weights | Citation range | Recency range | What it tests |
|--------|-------------|---------------|--------------|---------------|
| Baseline | phase3=12...theoretical=0.5 | 1.0–1.5× | 0.9–1.1× | Current default |
| Equal tiers | all=1.0 | 1.0–1.5× | 0.9–1.1× | Does tier weighting drive the nanoparticle drop? |
| Tier-only | baseline tiers | 1.0× (off) | 1.0× (off) | How much do citation and recency matter? |
| Flattened | phase3=4...theoretical=0.5 | 1.0–1.5× | 0.9–1.1× | Does compressed tier spread change rankings? |
| Citation-heavy | baseline tiers | 1.0–2.0× | 0.9–1.1× | Does citation impact dominate tier? |
| No recency | baseline tiers | 1.0–1.5× | 1.0× (off) | Does recency correction matter? |

## Reporting commitment

ALL 5 conclusions above will be reported regardless of outcome. Stable conclusions, directionally stable conclusions, and failures will all be documented explicitly. No selective omission.
