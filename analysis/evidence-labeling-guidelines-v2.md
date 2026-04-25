# Evidence Labeling Guidelines (v2)

Stable labeling guidelines for the evidence-tier gold set. Commit this document BEFORE any v2 labeling begins. All labelers must read this document and apply these definitions consistently.

## Purpose

The gold set measures how well the automated evidence tagger performs. Consistent human labels are the ground truth. Inconsistent labels make the evaluation meaningless. These guidelines define what each tier means, with examples and edge-case rulings derived from the v1 evaluation's confusion patterns.

## Labeling Procedure

1. Read the article's **title**, **abstract**, and **methods section** (if available in the full text)
2. Identify the **highest evidence tier** present in the article's **primary research** (not in cited references or background discussion)
3. Assign exactly one tier from the list below
4. If uncertain, consult the edge cases section. If still uncertain, assign the **more conservative (lower) tier** and note the uncertainty in `gold_notes`

## Tier Definitions

### phase3-clinical
Randomized controlled trial (RCT) or pivotal Phase III trial. Multi-center, powered for primary efficacy endpoint. Examples: EF-14 trial for TTFields in glioblastoma, VISION trial for lutetium-177 PSMA.

### phase2-clinical
Single-arm Phase II, randomized Phase II, or dose-escalation Phase II trial. Reports efficacy endpoints (response rate, PFS, OS) in a defined patient cohort. Includes expansion cohorts from Phase I/II designs.

### phase1-clinical
First-in-human, dose-finding, safety/tolerability studies. Primary endpoint is safety (MTD, DLT). May report preliminary efficacy signals but not powered for them.

### clinical-other
Patient-level data without phase designation. Includes: case reports, case series, retrospective chart reviews, registry analyses, compassionate-use reports, real-world evidence studies, single-patient expanded access. The key criterion: **actual patient data is reported** (not just discussed in a review context).

### preclinical-invivo
Animal model studies using cancer models: xenograft (human tumor in immunodeficient mice), syngeneic (murine tumor in immunocompetent mice), genetically engineered mouse models (GEM), patient-derived xenografts (PDX). Must include in-vivo tumor measurement (volume, survival, imaging).

### preclinical-invitro
Cell line studies, organoid experiments, 2D/3D culture assays. Reports quantitative results (viability, IC50, gene expression, Western blot, flow cytometry) from cancer-related cell experiments. Includes patient-derived organoids (NOT patient-level evidence — the organoid is the experimental unit, not the patient).

### theoretical
Computational modeling, mathematical simulation, bioinformatics analysis, network pharmacology, molecular docking. The article presents **novel computational results** (not just a review of existing computational work). Includes machine learning models trained on cancer data, Monte Carlo simulations, and pathway analysis with novel predictions.

### none-applicable
Reviews (narrative, systematic, scoping), meta-analyses, protocols, editorials, commentaries, opinion pieces, methodology papers without cancer-specific primary data. Also includes papers that discuss multiple evidence tiers but do not report their own primary research.

## Edge Cases (from v1 confusion analysis)

These rulings address the three largest confusion patterns from the v1 evaluation (34 of 54 errors):

### Computational vs Review (15 errors in v1)
- Computational modeling paper that **generates novel predictions** → **theoretical**
- Review paper that **summarizes existing computational methods** → **none-applicable**
- Bioinformatics paper that **reanalyzes public datasets with new methods** → **theoretical**
- Commentary on computational trends → **none-applicable**

### Clinical-other vs None-applicable (10 errors in v1)
- Retrospective analysis of patient outcomes → **clinical-other**
- Case report (even N=1) with treatment data → **clinical-other**
- Registry analysis of treatment patterns → **clinical-other**
- Review paper that mentions patient outcomes from cited studies → **none-applicable**
- Expert opinion citing clinical experience → **none-applicable**

### Preclinical-invitro vs None-applicable (9 errors in v1)
- Paper reporting IC50 values in cancer cell lines → **preclinical-invitro**
- Methodology paper developing a new assay (tested on cancer cells as proof of concept) → **preclinical-invitro** (if quantitative cancer-specific results are reported)
- Methodology paper describing assay development without cancer-specific endpoints → **none-applicable**
- Nanoparticle characterization without biological testing → **none-applicable**

### Additional Edge Cases
- Meta-analysis of Phase III trials → **none-applicable** (review-like, not primary evidence)
- Phase I/II study (combined design) → **phase1-clinical** (assign lowest phase when ambiguous)
- Patient-derived organoid study → **preclinical-invitro** (organoid is the experimental unit, not the patient)
- In-silico drug screening validated in cell lines → **preclinical-invitro** (highest evidence wins)
- CRISPR screen identifying cancer targets (no therapeutic intervention tested) → **theoretical**
