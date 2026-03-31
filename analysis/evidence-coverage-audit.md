# Evidence Coverage Audit

Evidence-level tags are present for 1826/4830 records (37.8%). Of the unclassified records, 1921 are review-like and 29 are protocol-like by design; 1054 primary-study-like records remain uncategorized. Primary-study-like evidence coverage is 1826/2880 (63.4%).

## Mechanisms Most Exposed To Overstated Absence Claims

| Mechanism | Total | Tagged | Review-like | Protocol-like | Other untagged | Primary-study-like coverage |
|---|---|---|---|---|---|---|
| **immunotherapy** | 2297 | 682 | 1161 | 15 | 439 | 682/1121 (60.8%) |
| **mRNA-vaccine** | 317 | 84 | 95 | 0 | 138 | 84/222 (37.8%) |
| **electrochemical-therapy** | 185 | 63 | 47 | 2 | 73 | 63/136 (46.3%) |
| **ttfields** | 262 | 115 | 72 | 4 | 71 | 115/186 (61.8%) |
| **synthetic-lethality** | 367 | 173 | 128 | 1 | 65 | 173/238 (72.7%) |
| **nanoparticle** | 515 | 228 | 230 | 0 | 57 | 228/285 (80.0%) |
| **bispecific-antibody** | 247 | 150 | 47 | 5 | 45 | 150/195 (76.9%) |
| **car-t** | 474 | 133 | 298 | 0 | 43 | 133/176 (75.6%) |
| **metabolic-targeting** | 274 | 74 | 166 | 0 | 34 | 74/108 (68.5%) |
| **antibody-drug-conjugate** | 284 | 127 | 121 | 2 | 34 | 127/161 (78.9%) |
| **crispr** | 331 | 144 | 155 | 0 | 32 | 144/176 (81.8%) |
| **hifu** | 81 | 22 | 29 | 1 | 29 | 22/51 (43.1%) |
| **epigenetic** | 183 | 82 | 71 | 1 | 29 | 82/111 (73.9%) |
| **oncolytic-virus** | 378 | 177 | 172 | 0 | 29 | 177/206 (85.9%) |
| **sonodynamic** | 187 | 113 | 50 | 0 | 24 | 113/137 (82.5%) |
| **frequency-therapy** | 71 | 14 | 36 | 1 | 20 | 14/34 (41.2%) |
| **bioelectric** | 182 | 144 | 19 | 0 | 19 | 144/163 (88.3%) |
| **microbiome** | 109 | 10 | 84 | 2 | 13 | 10/23 (43.5%) |
| **radioligand-therapy** | 52 | 20 | 23 | 1 | 8 | 20/28 (71.4%) |
| **electrolysis** | 11 | 5 | 0 | 0 | 6 | 5/11 (45.5%) |
| **phagocytosis-checkpoint** | 28 | 18 | 5 | 0 | 5 | 18/23 (78.3%) |
| **cold-atmospheric-plasma** | 3 | 1 | 2 | 0 | 0 | 1/1 (100.0%) |
| **targeted-protein-degradation** | 19 | 7 | 12 | 0 | 0 | 7/7 (100.0%) |

## Sample Of Unclassified Primary-Study-Like Records

Illustrative examples below come from the uncategorized primary-study-like pool rather than the review/protocol bucket. These are the records most likely to affect absence claims if the evidence classifier is expanded.


### immunotherapy

- **PMID 28436963** (2017) — *A STING-activating nanovaccine for cancer immunotherapy.*
- **PMID 35488273** (2022) — *Integrated analysis of single-cell and bulk RNA sequencing data reveals a pan-cancer stemness signature predicting immunotherapy response.*
- **PMID 33473101** (2021) — *Personalized neoantigen pulsed dendritic cell vaccine for advanced lung cancer.*

### mRNA-vaccine

- **PMID 33473101** (2021) — *Personalized neoantigen pulsed dendritic cell vaccine for advanced lung cancer.*
- **PMID 34861036** (2022) — *Efficacy of a third BNT162b2 mRNA COVID-19 vaccine dose in patients with CLL who failed standard 2-dose vaccination.*
- **PMID 35661819** (2022) — *Durable complete response to neoantigen-loaded dendritic-cell vaccine following anti-PD-1 therapy in metastatic gastric cancer.*

### electrochemical-therapy

- **PMID 25442132** (2014) — *Image-guided tumor ablation: standardization of terminology and reporting criteria--a 10-year update.*
- **PMID 25179590** (2014) — *Initial assessment of safety and clinical feasibility of irreversible electroporation in the focal treatment of prostate cancer.*
- **PMID 30986263** (2019) — *Prostate cancer treatment with Irreversible Electroporation (IRE): Safety, efficacy and clinical experience in 471 treatments.*

### ttfields

- **PMID 30534421** (2018) — *Tumor treating fields increases membrane permeability in glioblastoma cells.*
- **PMID 33080774** (2020) — *Tumor Treating Fields (TTFields) Hinder Cancer Cell Motility through Regulation of Microtubule and Acting Dynamics.*
- **PMID 23095807** (2012) — *Long-term survival of patients suffering from glioblastoma multiforme treated with tumor-treating fields.*

### synthetic-lethality

- **PMID 35444283** (2022) — *CCNE1 amplification is synthetic lethal with PKMYT1 kinase inhibition.*
- **PMID 38509368** (2024) — *Transcription-replication conflicts underlie sensitivity to PARP inhibitors.*
- **PMID 33333017** (2021) — *Defective ALC1 nucleosome remodeling confers PARPi sensitization and synthetic lethality with HRD.*

## What The Current Miss-Rate Signal Likely Means

- The raw 36.8% coverage number is pessimistic because review-like and protocol-like records are intentionally excluded from evidence tagging.
- The more relevant upper-bound miss rate is the share of `other_untagged` records within the primary-study-like subset. Mechanisms with the largest remaining uncertainty are immunotherapy, mRNA-vaccine, electrochemical-therapy, TTFields, and CAR-T.
- The sampled uncategorized records are enriched for observational clinical studies, biomarker/antigen-discovery papers, and translational engineering studies that do not announce phase or preclinical status in obvious keywords.
- This means the main risk is overstating `no detected clinical evidence` for modalities with many non-phase clinical or translational papers, not silently missing large numbers of explicit Phase III trials.

## Recommended Interpretation Guardrails

- Treat `0 Phase 2+` as `not detected in current keyword-derived evidence tags` unless manually verified.
- Distinguish review/protocol exclusions from true uncategorized primary-study-like records when discussing evidence coverage.
- Re-check any high-priority mechanism with external PubMed or trial-registry verification before using it as a headline gap.
- Prefer coverage-aware language in the manuscript and analysis files whenever evidence tagging is below 50% for a mechanism.