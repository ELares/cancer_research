# Diagnostic-to-Therapy Matching Layer

## Purpose

This layer surfaces the translational chain: **diagnostic modality → targetable feature → intervention class**. It complements the existing mechanism, evidence-tier, and convergence layers by answering a different question: not "what mechanisms exist?" but "what diagnostic led to what therapy choice?"

## Schema

Each chain has three keyword groups:

- **diagnostic**: the test, imaging, or profiling modality (e.g., "PSMA PET", "tumor mutational burden", "neoantigen prediction")
- **feature**: the targetable characteristic revealed (e.g., "PSMA expression", "TMB-high", "neoantigen")
- **intervention**: the therapy selected based on that feature (e.g., "177Lu-PSMA", "pembrolizumab", "mRNA vaccine")

## Matching Rule

A paper must mention the **intervention** keyword PLUS at least one of (**diagnostic** or **feature**). This is more conservative than 2-of-3 and prevents false positives from papers that discuss only a diagnostic imaging technique or only a biomarker without connecting it to a therapy choice.

## Current Coverage (6 chains, 4 modalities)

| Chain | Modality | Articles |
|-------|----------|----------|
| psma-imaging-to-radioligand | Radioligand | 4 |
| sstr-imaging-to-prrt | Radioligand | 0 |
| pdl1-ihc-to-checkpoint | Checkpoint selection | 13 |
| tmb-msi-to-immunotherapy | Checkpoint selection | 33 |
| neoantigen-profiling-to-mrna-vaccine | mRNA vaccine | 79 |
| oncolytic-susceptibility-to-virotherapy | Oncolytic virus | 1 |

## What This Layer Can Do

- Identify papers that connect a specific diagnostic to a specific therapy class
- Cross-tabulate with cancer types and evidence levels to show where diagnostic-therapy chains have clinical evidence
- Surface examples of biomarker-guided therapy selection across modalities
- Provide a starting point for understanding how different intervention choices are justified

## What This Layer Cannot Do

- Detect novel biomarker-therapy pairs not in the keyword dictionary
- Distinguish between papers that study the chain experimentally vs. papers that merely discuss it in a review context
- Replace clinical decision support systems or patient-matching tools
- Capture diagnostic-therapy chains expressed in non-standard terminology
- Quantify the strength of diagnostic-therapy evidence (use the evidence-tier layer for that)

## Known Limitations

- The SSTR-to-PRRT chain has 0 matches, likely because the corpus has very few PRRT-focused papers and the keywords are specific
- The oncolytic-susceptibility chain has only 1 match because viral receptor profiling as a patient-selection strategy is uncommon in the current literature
- PD-L1 and TMB keywords are ubiquitous in a cancer corpus; the intervention-required rule prevents most false positives but some broad review papers may still match
- The neoantigen chain (79 papers) is the strongest because the diagnostic-to-vaccine pipeline is tightly coupled in the mRNA vaccine literature

## How to Extend

1. **Add a new chain**: Define diagnostic, feature, and intervention keyword groups in `DIAGNOSTIC_THERAPY_KEYWORDS` in `scripts/config.py`. Add the chain ID to `DIAGNOSTIC_THERAPY_ORDER`.
2. **Run the pipeline**: `python scripts/tag_articles.py && python scripts/build_index.py && python scripts/analyze_corpus.py`
3. **Validate**: Check the generated `analysis/diagnostic-therapy-audit.md` and spot-check tagged papers.

## Related Files

- `scripts/config.py` — keyword definitions (`DIAGNOSTIC_THERAPY_KEYWORDS`)
- `scripts/tag_articles.py` — matching function (`match_diagnostic_therapy_links`)
- `scripts/build_index.py` — index field (`diagnostic_therapy_links`)
- `scripts/analyze_corpus.py` — analysis function (`build_diagnostic_therapy_audit`)
- `analysis/diagnostic-therapy-audit.md` — generated audit output
- `analysis/diagnostic-therapy-pilot.csv` — manually reviewed pilot sample
