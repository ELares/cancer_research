# Evidence Coverage Audit

Evidence-level tags are present for 1779/4830 records (36.8%). Of the unclassified records, 1921 are review-like and 29 are protocol-like by design; 1101 primary-study-like records remain uncategorized. Primary-study-like evidence coverage is 1779/2880 (61.8%).

## Mechanisms Most Exposed To Overstated Absence Claims

| Mechanism | Total | Tagged | Review-like | Protocol-like | Other untagged | Primary-study-like coverage |
|---|---|---|---|---|---|---|
| **immunotherapy** | 2297 | 660 | 1161 | 15 | 461 | 660/1121 (58.9%) |
| **mRNA-vaccine** | 317 | 81 | 95 | 0 | 141 | 81/222 (36.5%) |
| **electrochemical-therapy** | 185 | 61 | 47 | 2 | 75 | 61/136 (44.9%) |
| **ttfields** | 262 | 112 | 72 | 4 | 74 | 112/186 (60.2%) |
| **synthetic-lethality** | 367 | 173 | 128 | 1 | 65 | 173/238 (72.7%) |
| **nanoparticle** | 515 | 223 | 230 | 0 | 62 | 223/285 (78.2%) |
| **car-t** | 473 | 124 | 297 | 0 | 52 | 124/176 (70.5%) |
| **bispecific-antibody** | 247 | 146 | 47 | 5 | 49 | 146/195 (74.9%) |
| **antibody-drug-conjugate** | 284 | 124 | 121 | 2 | 37 | 124/161 (77.0%) |
| **crispr** | 331 | 141 | 155 | 0 | 35 | 141/176 (80.1%) |
| **metabolic-targeting** | 274 | 74 | 166 | 0 | 34 | 74/108 (68.5%) |
| **oncolytic-virus** | 378 | 173 | 172 | 0 | 33 | 173/206 (84.0%) |
| **epigenetic** | 183 | 80 | 71 | 1 | 31 | 80/111 (72.1%) |
| **hifu** | 81 | 22 | 29 | 1 | 29 | 22/51 (43.1%) |
| **sonodynamic** | 187 | 111 | 50 | 0 | 26 | 111/137 (81.0%) |
| **frequency-therapy** | 71 | 14 | 36 | 1 | 20 | 14/34 (41.2%) |
| **bioelectric** | 182 | 143 | 19 | 0 | 20 | 143/163 (87.7%) |
| **microbiome** | 109 | 9 | 84 | 2 | 14 | 9/23 (39.1%) |
| **radioligand-therapy** | 52 | 20 | 23 | 1 | 8 | 20/28 (71.4%) |
| **electrolysis** | 11 | 5 | 0 | 0 | 6 | 5/11 (45.5%) |
| **phagocytosis-checkpoint** | 28 | 18 | 5 | 0 | 5 | 18/23 (78.3%) |
| **cold-atmospheric-plasma** | 3 | 1 | 2 | 0 | 0 | 1/1 (100.0%) |
| **targeted-protein-degradation** | 19 | 7 | 12 | 0 | 0 | 7/7 (100.0%) |

## Sample Of Unclassified Primary-Study-Like Records

Illustrative examples below come from the uncategorized primary-study-like pool rather than the review/protocol bucket. These are the records most likely to affect absence claims if the evidence classifier is expanded.


### immunotherapy

- **PMID 33542131** (2021) — *Fecal microbiota transplant overcomes resistance to anti-PD-1 therapy in melanoma patients.*
- **PMID 37165196** (2023) — *Personalized RNA neoantigen vaccines stimulate T cells in pancreatic cancer.*
- **PMID 28436963** (2017) — *A STING-activating nanovaccine for cancer immunotherapy.*

### mRNA-vaccine

- **PMID 37165196** (2023) — *Personalized RNA neoantigen vaccines stimulate T cells in pancreatic cancer.*
- **PMID 35969778** (2022) — *Lipid nanoparticle-mediated lymph node-targeting delivery of mRNA cancer vaccine elicits robust CD8+ T cell response.*
- **PMID 33473101** (2021) — *Personalized neoantigen pulsed dendritic cell vaccine for advanced lung cancer.*

### electrochemical-therapy

- **PMID 25442132** (2014) — *Image-guided tumor ablation: standardization of terminology and reporting criteria--a 10-year update.*
- **PMID 25179590** (2014) — *Initial assessment of safety and clinical feasibility of irreversible electroporation in the focal treatment of prostate cancer.*
- **PMID 30986263** (2019) — *Prostate cancer treatment with Irreversible Electroporation (IRE): Safety, efficacy and clinical experience in 471 treatments.*

### ttfields

- **PMID 29260225** (2017) — *Effect of Tumor-Treating Fields Plus Maintenance Temozolomide vs Maintenance Temozolomide Alone on Survival in Patients With Glioblastoma: A Randomize*
- **PMID 26670971** (2015) — *Maintenance Therapy With Tumor-Treating Fields Plus Temozolomide vs Temozolomide Alone for Glioblastoma: A Randomized Clinical Trial.*
- **PMID 30534421** (2018) — *Tumor treating fields increases membrane permeability in glioblastoma cells.*

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