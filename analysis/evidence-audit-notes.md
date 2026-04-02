# Evidence Audit Notes

This note documents the main interpretation risk in the current corpus-wide evidence analysis.

## Current state

- The manuscript reports evidence-level coverage of 1,779 / 4,830 full-text articles.
- That is about 36.8% coverage.
- Any claim of "no Phase II/III evidence" should therefore be interpreted as "not detected in the current keyword-derived evidence layer" unless it has been externally verified.

## Why this matters

The repo uses evidence-level summaries in:

- [analysis/evidence-tiers.md](/Users/ezequiellares/go/src/github.com/ELares/cancer_research/analysis/evidence-tiers.md)
- [analysis/gap-analysis.md](/Users/ezequiellares/go/src/github.com/ELares/cancer_research/analysis/gap-analysis.md)
- [article/drafts/v1.md](/Users/ezequiellares/go/src/github.com/ELares/cancer_research/article/drafts/v1.md)

That is acceptable for hypothesis generation, but weak for strong absence claims.

## Highest-risk interpretation pattern

The repo is most likely to overstate absence when all three are true:

- mechanism article count is substantial
- evidence-tag coverage for that mechanism is low
- the mechanism has known external clinical activity outside the current query/taxonomy frame

Examples flagged in discussion and repo review:

- synthetic lethality outside the current PARP-heavy framing
- radioligand / theranostics
- targeted protein degradation
- phagocytosis checkpoints
- CAP and other non-core physical ROS modalities

## Recommended guardrails

- Prefer "not detected in our current keyword-based evidence analysis" over "no clinical evidence exists."
- Re-check important absence claims with external PubMed or registry verification before using them in the manuscript.
- Treat evidence tables as a coverage layer, not as ground truth.

## Follow-up

The pipeline code now includes better hooks for broader taxonomy work, but a full evidence-tagging audit still needs either:

- improved automated classification, or
- manual sampling of untagged records by mechanism

That follow-up now exists as a first-pass gold-set evaluation:

- `analysis/evidence-gold-set-v1.csv`
- `analysis/evidence-gold-labels-v1.csv`
- `analysis/evidence-gold-eval.md`

The current measured result is that the heuristic tagger behaves conservatively:

- binary evidence-detection precision is high in the sampled set
- recall is much lower, especially for `clinical-other`, `theoretical`, and unclassified preclinical records

So evidence-level absence claims should still be phrased as `not detected by the current heuristic layer` unless they have been externally verified.
