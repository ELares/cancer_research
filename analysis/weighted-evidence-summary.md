# Weighted Evidence Summary

Heuristic weighting of detected evidence by tier, citation percentile, and recency.

This is a ranking aid, not a formal study-quality score. It only applies to records with detected evidence tags, and it inherits the tagger’s conservative recall.

Weight formula: `tier_weight × citation_modifier × recency_modifier`, with evidence tier as the dominant term.

Multi-tag papers contribute to every mechanism they are tagged with, so scores are useful for within-lane ranking but are not independent or additive across mechanisms.

## Weighted Ranking By Mechanism

| Mechanism | Weighted score | Tagged evidence rows | Primary-study-like coverage | Avg weight per tagged row |
|---|---|---|---|---|
| **immunotherapy** | 2616.8 | 769 | 769/1121 (68.6%) | 3.40 |
| **bispecific-antibody** | 687.0 | 168 | 168/195 (86.2%) | 4.09 |
| **ttfields** | 568.4 | 138 | 138/186 (74.2%) | 4.12 |
| **antibody-drug-conjugate** | 557.7 | 143 | 143/161 (88.8%) | 3.90 |
| **oncolytic-virus** | 493.7 | 182 | 182/206 (88.3%) | 2.71 |
| **nanoparticle** | 492.5 | 228 | 228/285 (80.0%) | 2.16 |
| **synthetic-lethality** | 450.7 | 179 | 179/238 (75.2%) | 2.52 |
| **car-t** | 369.5 | 136 | 136/176 (77.3%) | 2.72 |
| **crispr** | 330.3 | 144 | 144/176 (81.8%) | 2.29 |
| **electrochemical-therapy** | 285.0 | 90 | 90/136 (66.2%) | 3.17 |
| **epigenetic** | 277.7 | 85 | 85/111 (76.6%) | 3.27 |
| **sonodynamic** | 258.9 | 113 | 113/137 (82.5%) | 2.29 |
| **bioelectric** | 244.0 | 144 | 144/163 (88.3%) | 1.69 |
| **mRNA-vaccine** | 192.3 | 67 | 67/119 (56.3%) | 2.87 |
| **metabolic-targeting** | 157.7 | 76 | 76/108 (70.4%) | 2.08 |
| **hifu** | 88.8 | 29 | 29/51 (56.9%) | 3.06 |
| **frequency-therapy** | 51.6 | 21 | 21/34 (61.8%) | 2.46 |
| **microbiome** | 43.7 | 12 | 12/23 (52.2%) | 3.64 |
| **phagocytosis-checkpoint** | 37.3 | 18 | 18/23 (78.3%) | 2.07 |
| **radioligand-therapy** | 21.9 | 8 | 8/9 (88.9%) | 2.74 |
| **targeted-protein-degradation** | 14.8 | 7 | 7/7 (100.0%) | 2.11 |
| **electrolysis** | 11.0 | 5 | 5/11 (45.5%) | 2.20 |
| **cold-atmospheric-plasma** | 2.7 | 1 | 1/1 (100.0%) | 2.75 |

## Top Weighted Studies By Mechanism


### immunotherapy

- **PMID 36600534** (2023) — `phase3-clinical` — weight 18.76 — *Refining adjuvant treatment in endometrial cancer based on molecular features: the RAINBO clinical trial program.*
- **PMID 38454445** (2024) — `phase3-clinical` — weight 17.80 — *S100A9+CD14+ monocytes contribute to anti-PD-1 immunotherapy resistance in advanced hepatocellular carcinoma by attenuating T cell-mediated *
- **PMID 37824808** (2024) — `phase3-clinical` — weight 17.65 — *Ibrutinib-based therapy reinvigorates CD8+ T cells compared to chemoimmunotherapy: immune monitoring from the E1912 trial.*

### bispecific-antibody

- **PMID 37824808** (2024) — `phase3-clinical` — weight 17.65 — *Ibrutinib-based therapy reinvigorates CD8+ T cells compared to chemoimmunotherapy: immune monitoring from the E1912 trial.*
- **PMID 40248696** (2025) — `phase3-clinical` — weight 12.98 — *Cadonilimab plus chemotherapy as first-line treatment for persistent, recurrent, or metastatic cervical cancer: a cost-effectiveness analysi*
- **PMID 41037766** (2025) — `phase3-clinical` — weight 12.98 — *Mosunetuzumab Plus Polatuzumab Vedotin in Transplant-Ineligible Refractory/Relapsed Large B-Cell Lymphoma: Primary Results of the Phase III *

### ttfields

- **PMID 30506499** (2019) — `phase3-clinical` — weight 17.12 — *Increased compliance with tumor treating fields therapy is prognostic for improved survival in the treatment of glioblastoma: a subgroup ana*
- **PMID 29260225** (2017) — `phase3-clinical` — weight 16.85 — *Effect of Tumor-Treating Fields Plus Maintenance Temozolomide vs Maintenance Temozolomide Alone on Survival in Patients With Glioblastoma: A*
- **PMID 32535723** (2020) — `phase3-clinical` — weight 16.53 — *Global post-marketing safety surveillance of Tumor Treating Fields (TTFields) in patients with high-grade glioma in clinical practice.*

### antibody-drug-conjugate

- **PMID 39232496** (2024) — `phase3-clinical` — weight 18.18 — *CMG901, a Claudin18.2-specific antibody-drug conjugate, for the treatment of solid tumors.*
- **PMID 41037766** (2025) — `phase3-clinical` — weight 12.98 — *Mosunetuzumab Plus Polatuzumab Vedotin in Transplant-Ineligible Refractory/Relapsed Large B-Cell Lymphoma: Primary Results of the Phase III *
- **PMID 40460679** (2025) — `phase3-clinical` — weight 12.98 — *The association of high body mass index with the safety and efficacy of sacituzumab govitecan in patients with metastatic triple-negative br*

### oncolytic-virus

- **PMID 27298410** (2016) — `phase3-clinical` — weight 16.45 — *Talimogene Laherparepvec in Combination With Ipilimumab in Previously Untreated, Unresectable Stage IIIB-IV Melanoma.*
- **PMID 38378555** (2024) — `phase2-clinical` — weight 12.38 — *Reshaping the tumor microenvironment of cold soft-tissue sarcomas with oncolytic viral therapy: a phase 2 trial of intratumoral JX-594 combi*
- **PMID 37142291** (2023) — `phase2-clinical` — weight 12.32 — *Talimogene laherparepvec in combination with ipilimumab versus ipilimumab alone for advanced melanoma: 5-year final analysis of a multicente*

### nanoparticle

- **PMID 40515479** (2025) — `phase2-clinical` — weight 8.65 — *First-in-human phase I/II, open-label study of mRNA-2416 alone or combined with durvalumab in patients with advanced solid tumors and ovaria*
- **PMID 37165196** (2023) — `phase1-clinical` — weight 7.84 — *Personalized RNA neoantigen vaccines stimulate T cells in pancreatic cancer.*
- **PMID 37323470** (2023) — `phase1-clinical` — weight 6.74 — *A personalized mRNA vaccine has exhibited potential in the treatment of pancreatic cancer.*

### synthetic-lethality

- **PMID 36082969** (2023) — `phase3-clinical` — weight 18.81 — *Overall Survival With Maintenance Olaparib at a 7-Year Follow-Up in Patients With Newly Diagnosed Advanced Ovarian Cancer and a BRCA Mutatio*
- **PMID 36600534** (2023) — `phase3-clinical` — weight 18.76 — *Refining adjuvant treatment in endometrial cancer based on molecular features: the RAINBO clinical trial program.*
- **PMID 37552839** (2023) — `phase2-clinical` — weight 12.49 — *MRTX1719 Is an MTA-Cooperative PRMT5 Inhibitor That Exhibits Synthetic Lethality in Preclinical Models and Patients with MTAP-Deleted Cancer*

### car-t

- **PMID 39908461** (2025) — `phase2-clinical` — weight 8.65 — *Phase 2 trial of ibrutinib and nivolumab in patients with relapsed CNS lymphomas.*
- **PMID 38583184** (2024) — `phase1-clinical` — weight 7.96 — *CD70-Targeted Allogeneic CAR T-Cell Therapy for Advanced Clear Cell Renal Cell Carcinoma.*
- **PMID 39735354** (2024) — `phase1-clinical` — weight 7.79 — *Recent Advances and Future Directions in Sonodynamic Therapy for Cancer Treatment.*

## Guardrails

- Evidence tier dominates the score. Citation percentile and recency only adjust within-tier ordering.
- Scores are taxonomy-overlap dependent because the same study can legitimately contribute to umbrella and subclass mechanism lanes.
- These weights do not estimate true study quality, patient benefit, or sample size.
- Mechanisms with low evidence-tag coverage can still be under-ranked even if their real literature is stronger.
- The gold-set evaluation suggests the tagger is conservative, so use this report as `quality among detected evidence`, not `quality of the whole field`.