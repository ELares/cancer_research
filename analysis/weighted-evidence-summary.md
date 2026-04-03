# Weighted Evidence Summary

Heuristic weighting of detected evidence by tier, citation percentile, recency, and a small clinical flag modifier.

This is a ranking aid, not a formal study-quality score. It only applies to records with detected evidence tags, and it inherits the tagger’s conservative recall.

Weight formula: `tier_weight × citation_modifier × recency_modifier × clinical_modifier`, with evidence tier as the dominant term.

Multi-tag papers contribute to every mechanism they are tagged with, so scores are useful for within-lane ranking but are not independent or additive across mechanisms.

## Weighted Ranking By Mechanism

| Mechanism | Weighted score | Tagged evidence rows | Primary-study-like coverage | Avg weight per tagged row |
|---|---|---|---|---|
| **immunotherapy** | 2701.8 | 763 | 763/1121 (68.1%) | 3.54 |
| **bispecific-antibody** | 731.7 | 168 | 168/195 (86.2%) | 4.36 |
| **ttfields** | 585.0 | 135 | 135/186 (72.6%) | 4.33 |
| **antibody-drug-conjugate** | 574.4 | 140 | 140/161 (87.0%) | 4.10 |
| **oncolytic-virus** | 503.5 | 182 | 182/206 (88.3%) | 2.77 |
| **nanoparticle** | 494.7 | 228 | 228/285 (80.0%) | 2.17 |
| **synthetic-lethality** | 452.1 | 178 | 178/238 (74.8%) | 2.54 |
| **car-t** | 397.7 | 136 | 136/176 (77.3%) | 2.92 |
| **crispr** | 331.0 | 144 | 144/176 (81.8%) | 2.30 |
| **epigenetic** | 286.5 | 85 | 85/111 (76.6%) | 3.37 |
| **electrochemical-therapy** | 274.1 | 86 | 86/136 (63.2%) | 3.19 |
| **sonodynamic** | 258.9 | 113 | 113/137 (82.5%) | 2.29 |
| **bioelectric** | 244.0 | 144 | 144/163 (88.3%) | 1.69 |
| **mRNA-vaccine** | 198.2 | 67 | 67/119 (56.3%) | 2.96 |
| **metabolic-targeting** | 158.2 | 76 | 76/108 (70.4%) | 2.08 |
| **hifu** | 91.3 | 29 | 29/51 (56.9%) | 3.15 |
| **microbiome** | 53.0 | 12 | 12/23 (52.2%) | 4.41 |
| **frequency-therapy** | 49.2 | 20 | 20/34 (58.8%) | 2.46 |
| **phagocytosis-checkpoint** | 37.3 | 18 | 18/23 (78.3%) | 2.07 |
| **radioligand-therapy** | 21.9 | 8 | 8/9 (88.9%) | 2.74 |
| **targeted-protein-degradation** | 14.8 | 7 | 7/7 (100.0%) | 2.11 |
| **electrolysis** | 11.0 | 5 | 5/11 (45.5%) | 2.20 |
| **cold-atmospheric-plasma** | 2.7 | 1 | 1/1 (100.0%) | 2.75 |

## Top Weighted Studies By Mechanism


### immunotherapy

- **PMID 36814222** (2023) — `phase3-clinical` — weight 19.13 — *INTEGRATE II: randomised phase III controlled trials of regorafenib containing regimens versus standard of care in refractory Advanced Gastr*
- **PMID 38245601** (2024) — `phase3-clinical` — weight 19.06 — *Real-world analysis of teclistamab in 123 RRMM patients from Germany.*
- **PMID 36600534** (2023) — `phase3-clinical` — weight 18.76 — *Refining adjuvant treatment in endometrial cancer based on molecular features: the RAINBO clinical trial program.*

### bispecific-antibody

- **PMID 38245601** (2024) — `phase3-clinical` — weight 19.06 — *Real-world analysis of teclistamab in 123 RRMM patients from Germany.*
- **PMID 37824808** (2024) — `phase3-clinical` — weight 17.65 — *Ibrutinib-based therapy reinvigorates CD8+ T cells compared to chemoimmunotherapy: immune monitoring from the E1912 trial.*
- **PMID 41037766** (2025) — `phase3-clinical` — weight 14.28 — *Mosunetuzumab Plus Polatuzumab Vedotin in Transplant-Ineligible Refractory/Relapsed Large B-Cell Lymphoma: Primary Results of the Phase III *

### ttfields

- **PMID 29260225** (2017) — `phase3-clinical` — weight 18.54 — *Effect of Tumor-Treating Fields Plus Maintenance Temozolomide vs Maintenance Temozolomide Alone on Survival in Patients With Glioblastoma: A*
- **PMID 26670971** (2015) — `phase3-clinical` — weight 17.81 — *Maintenance Therapy With Tumor-Treating Fields Plus Temozolomide vs Temozolomide Alone for Glioblastoma: A Randomized Clinical Trial.*
- **PMID 30506499** (2019) — `phase3-clinical` — weight 17.12 — *Increased compliance with tumor treating fields therapy is prognostic for improved survival in the treatment of glioblastoma: a subgroup ana*

### antibody-drug-conjugate

- **PMID 39232496** (2024) — `phase3-clinical` — weight 18.18 — *CMG901, a Claudin18.2-specific antibody-drug conjugate, for the treatment of solid tumors.*
- **PMID 41037766** (2025) — `phase3-clinical` — weight 14.28 — *Mosunetuzumab Plus Polatuzumab Vedotin in Transplant-Ineligible Refractory/Relapsed Large B-Cell Lymphoma: Primary Results of the Phase III *
- **PMID 40460679** (2025) — `phase3-clinical` — weight 14.28 — *The association of high body mass index with the safety and efficacy of sacituzumab govitecan in patients with metastatic triple-negative br*

### oncolytic-virus

- **PMID 27298410** (2016) — `phase3-clinical` — weight 18.10 — *Talimogene Laherparepvec in Combination With Ipilimumab in Previously Untreated, Unresectable Stage IIIB-IV Melanoma.*
- **PMID 38378555** (2024) — `phase2-clinical` — weight 13.61 — *Reshaping the tumor microenvironment of cold soft-tissue sarcomas with oncolytic viral therapy: a phase 2 trial of intratumoral JX-594 combi*
- **PMID 37142291** (2023) — `phase2-clinical` — weight 13.55 — *Talimogene laherparepvec in combination with ipilimumab versus ipilimumab alone for advanced melanoma: 5-year final analysis of a multicente*

### nanoparticle

- **PMID 40515479** (2025) — `phase2-clinical` — weight 9.52 — *First-in-human phase I/II, open-label study of mRNA-2416 alone or combined with durvalumab in patients with advanced solid tumors and ovaria*
- **PMID 37165196** (2023) — `phase1-clinical` — weight 8.62 — *Personalized RNA neoantigen vaccines stimulate T cells in pancreatic cancer.*
- **PMID 37323470** (2023) — `phase1-clinical` — weight 6.74 — *A personalized mRNA vaccine has exhibited potential in the treatment of pancreatic cancer.*

### synthetic-lethality

- **PMID 36082969** (2023) — `phase3-clinical` — weight 20.69 — *Overall Survival With Maintenance Olaparib at a 7-Year Follow-Up in Patients With Newly Diagnosed Advanced Ovarian Cancer and a BRCA Mutatio*
- **PMID 36600534** (2023) — `phase3-clinical` — weight 18.76 — *Refining adjuvant treatment in endometrial cancer based on molecular features: the RAINBO clinical trial program.*
- **PMID 37552839** (2023) — `phase2-clinical` — weight 13.74 — *MRTX1719 Is an MTA-Cooperative PRMT5 Inhibitor That Exhibits Synthetic Lethality in Preclinical Models and Patients with MTAP-Deleted Cancer*

### car-t

- **PMID 38245601** (2024) — `phase3-clinical` — weight 19.06 — *Real-world analysis of teclistamab in 123 RRMM patients from Germany.*
- **PMID 39657136** (2025) — `phase3-clinical` — weight 12.98 — *Lisocabtagene maraleucel for relapsed/refractory large B-cell lymphoma: a cell therapy consortium real-world analysis.*
- **PMID 39908461** (2025) — `phase2-clinical` — weight 9.52 — *Phase 2 trial of ibrutinib and nivolumab in patients with relapsed CNS lymphomas.*

## Guardrails

- Evidence tier dominates the score. Citation percentile and recency only adjust within-tier ordering.
- These weights do not estimate true study quality, patient benefit, or sample size.
- Mechanisms with low evidence-tag coverage can still be under-ranked even if their real literature is stronger.
- The gold-set evaluation suggests the tagger is conservative, so use this report as `quality among detected evidence`, not `quality of the whole field`.