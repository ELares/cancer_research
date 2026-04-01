# Designed Combination Audit

First-pass separation of broad multi-mechanism co-mentions from papers that look like deliberate combination studies.

This layer is heuristic. It is designed to complement the existing co-occurrence map, not replace it.

## Schema

- `co-mention-only`: multi-tagged paper without strong designed-combination language.
- `designed-combination-preclinical`: preclinical paper with explicit combination/synergy language.
- `designed-combination-clinical`: patient-study signal plus explicit combination language.
- `review-or-perspective-multi-lane`: review/prospective paper spanning multiple lanes.

## Corpus-Level Counts

- **designed-combination-clinical**: 67 (3.4% of 2+ mechanism papers)
- **designed-combination-preclinical**: 319 (16.3% of 2+ mechanism papers)
- **co-mention-only**: 621 (31.8% of 2+ mechanism papers)
- **review-or-perspective-multi-lane**: 945 (48.4% of 2+ mechanism papers)

## Highest-Count Designed Combination Lanes

| Mechanism pair | Designed-combination articles |
|---|---|
| **immunotherapy + oncolytic-virus** | 69 |
| **nanoparticle + sonodynamic** | 52 |
| **immunotherapy + nanoparticle** | 44 |
| **bispecific-antibody + immunotherapy** | 42 |
| **car-t + immunotherapy** | 33 |
| **immunotherapy + sonodynamic** | 22 |
| **epigenetic + immunotherapy** | 22 |
| **crispr + immunotherapy** | 20 |
| **immunotherapy + mRNA-vaccine** | 20 |
| **mRNA-vaccine + nanoparticle** | 13 |
| **crispr + synthetic-lethality** | 12 |
| **antibody-drug-conjugate + immunotherapy** | 10 |
| **car-t + crispr** | 10 |
| **immunotherapy + metabolic-targeting** | 10 |
| **bioelectric + nanoparticle** | 10 |

## Audited Priority Lanes

The samples below are manually reviewed examples selected from recent or highly cited records in the priority lanes discussed in issue #42.


### immunotherapy + oncolytic-virus

- **PMID 27298410** (2016) — `designed-combination-clinical` / `phase3-clinical` — *Talimogene Laherparepvec in Combination With Ipilimumab in Previously Untreated, Unresectable Stage IIIB-IV Melanoma.*
- **PMID 33232299** (2021) — `designed-combination-preclinical` / `preclinical-invitro` — *Zika virus oncolytic activity requires CD8+ T cells and is boosted by immune checkpoint blockade.*
- **PMID 37142291** (2023) — `designed-combination-clinical` / `phase2-clinical` — *Talimogene laherparepvec in combination with ipilimumab versus ipilimumab alone for advanced melanoma: 5-year final analysis of a multicenter, randomi*

### immunotherapy + mRNA-vaccine

- **PMID 36168634** (2022) — `designed-combination-preclinical` / `preclinical-invivo` — *Circular RNA cancer vaccines drive immunity in hard-to-treat malignancies.*
- **PMID 38584166** (2024) — `designed-combination-clinical` / `phase2-clinical` — *Personalized neoantigen vaccine and pembrolizumab in advanced hepatocellular carcinoma: a phase 1/2 trial.*
- **PMID 38268001** (2024) — `designed-combination-preclinical` / `preclinical-invivo` — *mRNA-based precision targeting of neoantigens and tumor-associated antigens in malignant brain tumors.*

### immunotherapy + radioligand-therapy

- **PMID 38698840** (2024) — `designed-combination-preclinical` / `preclinical-invivo` — *Designing combination therapies for cancer treatment: application of a mathematical framework combining CAR T-cell immunotherapy and targeted radionuc*
- **PMID 34063642** (2021) — `designed-combination-preclinical` / `preclinical-invivo` — *Combined Radionuclide Therapy and Immunotherapy for Treatment of Triple Negative Breast Cancer.*
- **PMID 38302933** (2024) — `review-or-perspective-multi-lane` / `protocol_like` — *A multicentric, single arm, open-label, phase I/II study evaluating PSMA targeted radionuclide therapy in adult patients with metastatic clear cell re*

### immunotherapy + sonodynamic

- **PMID 31048681** (2019) — `designed-combination-preclinical` / `preclinical-invitro` — *Checkpoint blockade and nanosonosensitizer-augmented noninvasive sonodynamic therapy combination reduces tumour growth and metastases in mice.*
- **PMID 37914681** (2023) — `designed-combination-preclinical` / `preclinical-invivo` — *Nanosensitizer-mediated augmentation of sonodynamic therapy efficacy and antitumor immunity.*
- **PMID 35568916** (2022) — `designed-combination-preclinical` / `preclinical-invivo` — *Enhancement of antitumor immunotherapy using mitochondria-targeted cancer cell membrane-biomimetic MOF-mediated sonodynamic therapy and checkpoint blo*

## Interpretation

- The designed-combination counts are materially smaller than the raw multi-tag co-occurrence totals, which confirms that convergence maps and designed-treatment maps should not be treated as interchangeable.
- Clinical combination signal is concentrated in a handful of lanes, especially immunotherapy-centered combinations. Much of the remaining multi-tag literature is still review-heavy or conceptual.
- This is a deliberately conservative first pass. The main purpose is to create a usable schema and an audited artifact before attempting more aggressive extraction.