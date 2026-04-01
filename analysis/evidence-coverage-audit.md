# Evidence Coverage Audit

Evidence-level tags are present for 2020/4830 records (41.8%). Of the unclassified records, 1921 are review-like and 29 are protocol-like by design; 860 primary-study-like records remain uncategorized. Primary-study-like evidence coverage is 2020/2880 (70.1%).

## Mechanisms Most Exposed To Overstated Absence Claims

| Mechanism | Total | Tagged | Review-like | Protocol-like | Other untagged | Primary-study-like coverage |
|---|---|---|---|---|---|---|
| **immunotherapy** | 2297 | 763 | 1161 | 15 | 358 | 763/1121 (68.1%) |
| **synthetic-lethality** | 367 | 178 | 128 | 1 | 60 | 178/238 (74.8%) |
| **nanoparticle** | 515 | 228 | 230 | 0 | 57 | 228/285 (80.0%) |
| **mRNA-vaccine** | 179 | 67 | 60 | 0 | 52 | 67/119 (56.3%) |
| **ttfields** | 262 | 135 | 72 | 4 | 51 | 135/186 (72.6%) |
| **electrochemical-therapy** | 185 | 86 | 47 | 2 | 50 | 86/136 (63.2%) |
| **car-t** | 474 | 136 | 298 | 0 | 40 | 136/176 (77.3%) |
| **metabolic-targeting** | 274 | 76 | 166 | 0 | 32 | 76/108 (70.4%) |
| **crispr** | 331 | 144 | 155 | 0 | 32 | 144/176 (81.8%) |
| **bispecific-antibody** | 247 | 168 | 47 | 5 | 27 | 168/195 (86.2%) |
| **epigenetic** | 183 | 85 | 71 | 1 | 26 | 85/111 (76.6%) |
| **sonodynamic** | 187 | 113 | 50 | 0 | 24 | 113/137 (82.5%) |
| **oncolytic-virus** | 378 | 182 | 172 | 0 | 24 | 182/206 (88.3%) |
| **hifu** | 81 | 29 | 29 | 1 | 22 | 29/51 (56.9%) |
| **antibody-drug-conjugate** | 284 | 140 | 121 | 2 | 21 | 140/161 (87.0%) |
| **bioelectric** | 182 | 144 | 19 | 0 | 19 | 144/163 (88.3%) |
| **frequency-therapy** | 71 | 20 | 36 | 1 | 14 | 20/34 (58.8%) |
| **microbiome** | 109 | 12 | 84 | 2 | 11 | 12/23 (52.2%) |
| **electrolysis** | 11 | 5 | 0 | 0 | 6 | 5/11 (45.5%) |
| **phagocytosis-checkpoint** | 28 | 18 | 5 | 0 | 5 | 18/23 (78.3%) |
| **radioligand-therapy** | 11 | 8 | 1 | 1 | 1 | 8/9 (88.9%) |
| **cold-atmospheric-plasma** | 3 | 1 | 2 | 0 | 0 | 1/1 (100.0%) |
| **targeted-protein-degradation** | 19 | 7 | 12 | 0 | 0 | 7/7 (100.0%) |

## Sample Of Unclassified Primary-Study-Like Records

Illustrative examples below come from the uncategorized primary-study-like pool rather than the review/protocol bucket. These are the records most likely to affect absence claims if the evidence classifier is expanded.


### immunotherapy

- **PMID 28436963** (2017) — *A STING-activating nanovaccine for cancer immunotherapy.*
- **PMID 35488273** (2022) — *Integrated analysis of single-cell and bulk RNA sequencing data reveals a pan-cancer stemness signature predicting immunotherapy response.*
- **PMID 31186412** (2019) — *LIF regulates CXCL9 in tumor-associated macrophages and prevents CD8+ T cell tumor-infiltration impairing anti-PD1 therapy.*

### synthetic-lethality

- **PMID 35444283** (2022) — *CCNE1 amplification is synthetic lethal with PKMYT1 kinase inhibition.*
- **PMID 38509368** (2024) — *Transcription-replication conflicts underlie sensitivity to PARP inhibitors.*
- **PMID 33333017** (2021) — *Defective ALC1 nucleosome remodeling confers PARPi sensitization and synthetic lethality with HRD.*

### nanoparticle

- **PMID 28436963** (2017) — *A STING-activating nanovaccine for cancer immunotherapy.*
- **PMID 30881202** (2018) — *Nanoparticle-assisted ultrasound: A special focus on sonodynamic therapy against cancer.*
- **PMID 37626073** (2023) — *A photo-triggered self-accelerated nanoplatform for multifunctional image-guided combination cancer immunotherapy.*

### mRNA-vaccine

- **PMID 34872567** (2021) — *Tumor antigens and immune subtypes guided mRNA vaccine development for kidney renal clear cell carcinoma.*
- **PMID 37428918** (2023) — *Comb-structured mRNA vaccine tethered with short double-stranded RNA adjuvants maximizes cellular immunity for cancer treatment.*
- **PMID 35185876** (2022) — *Identification of Tumor Antigens and Immune Subtypes of Glioblastoma for mRNA Vaccine Development.*

### ttfields

- **PMID 30534421** (2018) — *Tumor treating fields increases membrane permeability in glioblastoma cells.*
- **PMID 33080774** (2020) — *Tumor Treating Fields (TTFields) Hinder Cancer Cell Motility through Regulation of Microtubule and Acting Dynamics.*
- **PMID 26558989** (2015) — *NovoTTF™-100A System (Tumor Treating Fields) transducer array layout planning for glioblastoma: a NovoTAL™ system user study.*

## What The Current Miss-Rate Signal Likely Means

- The raw 41.8% coverage number is pessimistic because review-like and protocol-like records are intentionally excluded from evidence tagging.
- The more relevant upper-bound miss rate is the share of `other_untagged` records within the primary-study-like subset. Mechanisms with the largest remaining uncertainty are immunotherapy, mRNA-vaccine, electrochemical-therapy, TTFields, and CAR-T.
- After adding a `clinical-other` bucket, the remaining uncategorized records are still enriched for translational engineering studies, biomarker/antigen-discovery papers, and mechanistic studies that do not announce phase or preclinical status in obvious keywords.
- The main residual risk is now twofold: under-classifying ambiguous patient studies that still do not emit clear textual signals, and overstating absence when key landmark papers are missing from the local full-text archive.
- See `analysis/landmark-corpus-gaps.md` for a small manually curated shortlist of known missing papers that are important enough to change field-level interpretation.

## Recommended Interpretation Guardrails

- Treat `0 Phase 2+` as `not detected in current keyword-derived evidence tags` unless manually verified.
- Treat `clinical-other` as non-phase patient-study signal that is informative for field maturity, but not interchangeable with registrational phase evidence.
- Distinguish review/protocol exclusions from true uncategorized primary-study-like records when discussing evidence coverage.
- Re-check any high-priority mechanism with external PubMed or trial-registry verification before using it as a headline gap.
- Prefer coverage-aware language in the manuscript and analysis files whenever evidence tagging is below 50% for a mechanism.