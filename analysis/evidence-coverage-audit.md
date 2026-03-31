# Evidence Coverage Audit

Evidence-level tags are present for 1779/4830 records (36.8%). 3051 records remain unclassified.

## Mechanisms Most Exposed To Overstated Absence Claims

| Mechanism | Total | Tagged for evidence | Coverage |
|---|---|---|---|
| **antibody-drug-conjugate** | 284 | 124 | 43.7% |
| **bioelectric** | 182 | 143 | 78.6% |
| **bispecific-antibody** | 247 | 146 | 59.1% |
| **car-t** | 473 | 124 | 26.2% |
| **cold-atmospheric-plasma** | 3 | 1 | 33.3% |
| **crispr** | 331 | 141 | 42.6% |
| **electrochemical-therapy** | 185 | 61 | 33.0% |
| **electrolysis** | 11 | 5 | 45.5% |
| **epigenetic** | 183 | 80 | 43.7% |
| **frequency-therapy** | 71 | 14 | 19.7% |
| **hifu** | 81 | 22 | 27.2% |
| **immunotherapy** | 2297 | 660 | 28.7% |
| **mRNA-vaccine** | 317 | 81 | 25.6% |
| **metabolic-targeting** | 274 | 74 | 27.0% |
| **microbiome** | 109 | 9 | 8.3% |
| **nanoparticle** | 515 | 223 | 43.3% |
| **oncolytic-virus** | 378 | 173 | 45.8% |
| **phagocytosis-checkpoint** | 28 | 18 | 64.3% |
| **radioligand-therapy** | 52 | 20 | 38.5% |
| **sonodynamic** | 187 | 111 | 59.4% |
| **synthetic-lethality** | 367 | 173 | 47.1% |
| **targeted-protein-degradation** | 19 | 7 | 36.8% |
| **ttfields** | 262 | 112 | 42.7% |

## Recommended Interpretation Guardrails

- Treat `0 Phase 2+` as `not detected in current keyword-derived evidence tags` unless manually verified.
- Re-check any high-priority mechanism with external PubMed or trial-registry verification before using it as a headline gap.
- Prefer coverage-aware language in the manuscript and analysis files whenever evidence tagging is below 50% for a mechanism.