# Pathway Target Audit

First-pass tracking for ferroptosis-resistance and adjacent cell-death pathway targets that were previously present in the corpus text but not modeled as a dedicated layer.

Current pathway-target coverage in the index: 179/4830 records (3.7%).

## Target Counts

Counts below are split so broad review coverage does not get conflated with pathway-centered primary-study-like signal.

| Pathway target | Total | Primary-study-like | Review-like | Protocol-like | Top mechanisms | Top cancers |
|---|---|---|---|---|---|---|
| **cuproptosis-core** | 38 | 7 | 31 | 0 | immunotherapy (25), nanoparticle (15), metabolic-targeting (7) | breast (4), liver (2), pancreatic (1) |
| **dhcr7-7dhc-axis** | 7 | 2 | 5 | 0 | immunotherapy (4), nanoparticle (2), microbiome (1) | kidney (1), breast (1), melanoma (1) |
| **dhodh-defense** | 22 | 5 | 17 | 0 | metabolic-targeting (9), immunotherapy (7), epigenetic (3) | leukemia (3), glioblastoma (3), prostate (2) |
| **disulfidptosis-core** | 87 | 28 | 59 | 0 | immunotherapy (45), metabolic-targeting (25), nanoparticle (14) | lung (6), ovarian (5), colorectal (5) |
| **fdx1-cuproptosis-axis** | 38 | 7 | 31 | 0 | immunotherapy (25), nanoparticle (15), metabolic-targeting (7) | breast (4), liver (2), pancreatic (1) |
| **mboat1-mboat2-axis** | 4 | 0 | 4 | 0 | immunotherapy (3), epigenetic (1) | colorectal (1) |
| **scd-mufa-axis** | 62 | 12 | 50 | 0 | metabolic-targeting (32), immunotherapy (27), crispr (10) | liver (9), breast (5), colorectal (3) |
| **trim25-gpx4-degradation** | 7 | 1 | 6 | 0 | immunotherapy (5), mRNA-vaccine (1), synthetic-lethality (1) | glioblastoma (1), breast (1), ovarian (1) |

## Example Articles

Examples prefer primary-study-like records when available, then fall back to the most-cited review-like articles.


### disulfidptosis-core

- **PMID 37889752** (2023) — crispr, immunotherapy — `preclinical-invitro` — *KEAP1 mutation in lung adenocarcinoma promotes immune evasion and immunotherapy resistance.*
- **PMID 36594611** (2023) — sonodynamic — `other_untagged` — *Covalent Organic Framework Nanobowls as Activatable Nanosensitizers for Tumor-Specific and Ferroptosis-Augmented Sonodynamic Therapy.*
- **PMID 37369808** (2023) — immunotherapy, mRNA-vaccine — `preclinical-invitro` — *Crosstalk of ferroptosis regulators and tumor immunity in pancreatic adenocarcinoma: novel perspective to mRNA vaccines and personalized immunotherapy*

### scd-mufa-axis

- **PMID 35859734** (2022) — untagged — `preclinical-invitro` — *Targeting stearoyl-coa desaturase enhances radiation induced ferroptosis and immunogenic cell death in esophageal squamous cell carcinoma.*
- **PMID 32660617** (2020) — untagged — `preclinical-invivo` — *Progesterone receptor membrane component 1 regulates lipid homeostasis and drives oncogenic signaling resulting in breast cancer progression.*
- **PMID 31819198** (2020) — metabolic-targeting — `other_untagged` — *Development of cancer metabolism as a therapeutic target: new pathways, patient studies, stratification and combination therapy.*

### cuproptosis-core

- **PMID 39209877** (2024) — bioelectric, nanoparticle — `preclinical-invitro` — *A singular plasmonic-thermoelectric hollow nanostructure inducing apoptosis and cuproptosis for catalytic cancer therapy.*
- **PMID 38938647** (2024) — immunotherapy — `other_untagged` — *Copper(II)-Based Nano-Regulator Correlates Cuproptosis Burst and Sequential Immunogenic Cell Death for Synergistic Cancer Immunotherapy.*
- **PMID 25482950** (2014) — bioelectric, metabolic-targeting — `preclinical-invitro` — *Targeting of two aspects of metabolism in breast cancer treatment.*

### fdx1-cuproptosis-axis

- **PMID 39209877** (2024) — bioelectric, nanoparticle — `preclinical-invitro` — *A singular plasmonic-thermoelectric hollow nanostructure inducing apoptosis and cuproptosis for catalytic cancer therapy.*
- **PMID 38938647** (2024) — immunotherapy — `other_untagged` — *Copper(II)-Based Nano-Regulator Correlates Cuproptosis Burst and Sequential Immunogenic Cell Death for Synergistic Cancer Immunotherapy.*
- **PMID 25482950** (2014) — bioelectric, metabolic-targeting — `preclinical-invitro` — *Targeting of two aspects of metabolism in breast cancer treatment.*

### dhodh-defense

- **PMID 35103292** (2022) — epigenetic — `preclinical-invitro` — *Epigenetic therapy with chidamide alone or combined with 5‑azacitidine exerts antitumour effects on acute myeloid leukaemia cells in vitro.*
- **PMID 35286311** (2022) — metabolic-targeting, synthetic-lethality — `preclinical-invitro` — *A network-based approach to integrate nutrient microenvironment in the prediction of synthetic lethality in cancer metabolism.*
- **PMID 40490790** (2025) — immunotherapy, nanoparticle — `preclinical-invitro` — *Sonodynamic therapy-boosted biomimetic nanoplatform targets ferroptosis and CD47 as vulnerabilities for cancer immunotherapy.*

## Interpretation

- `scd-mufa-axis` and `disulfidptosis-core` already have enough corpus presence to affect how the repo frames in vivo ferroptosis escape and residual-state vulnerabilities.
- `dhodh-defense`, `dhcr7-7dhc-axis`, `fdx1-cuproptosis-axis`, and `trim25-gpx4-degradation` are smaller but non-zero. They should be treated as candidate stratification or escape markers rather than ignored side notes.
- The total counts are still inflated by broad reviews and pathway-survey papers, so prioritization should use the primary-study-like column rather than the raw total alone.
- The key repo-level shift is from modality-only comparison to vulnerability-layer comparison: these targets help explain when ferroptosis logic fails, and when adjacent programs like cuproptosis or disulfidptosis may be more relevant.