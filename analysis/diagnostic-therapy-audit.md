# Diagnostic-to-Therapy Matching Audit

First-pass extraction of diagnostic → targetable feature → intervention chains. Matching requires the intervention link plus at least one of (diagnostic, feature). Chain membership is recomputed from the frozen corpus text using the current chain set, so the index itself is never mutated (#441).

- Articles with at least one diagnostic-therapy link: **240** / 4830
- Chains evaluated: 10

## Chain Counts

| Chain | Articles | Top cancer types | Top evidence levels |
|---|---|---|---|
| **psma-imaging-to-radioligand** | 4 | prostate (3), kidney (1), breast (1) | preclinical-invivo (1) |
| **sstr-imaging-to-prrt** | 0 | . | . |
| **pdl1-ihc-to-checkpoint** | 13 | lung (6), bladder (2), melanoma (2) | phase2-clinical (4), phase1-clinical (2), clinical-other (1) |
| **tmb-msi-to-immunotherapy** | 33 | colorectal (15), lung (4), melanoma (3) | preclinical-invivo (5), clinical-other (4), phase2-clinical (2) |
| **neoantigen-profiling-to-mrna-vaccine** | 79 | melanoma (18), pancreatic (16), lung (8) | preclinical-invivo (12), phase1-clinical (6), preclinical-invitro (5) |
| **oncolytic-susceptibility-to-virotherapy** | 1 | melanoma (1) | preclinical-invivo (1) |
| **her2-testing-to-trastuzumab** | 28 | breast (20), gastric (4), lung (2) | clinical-other (5), preclinical-invivo (4), phase2-clinical (3) |
| **brca-mutation-to-parp-inhibitor** | 71 | ovarian (30), breast (29), pancreatic (5) | preclinical-invivo (16), preclinical-invitro (10), clinical-other (2) |
| **egfr-mutation-to-egfr-inhibitor** | 12 | lung (12) | preclinical-invivo (4), clinical-other (4), phase1-clinical (1) |
| **kras-g12c-mutation-to-sotorasib** | 0 | . | . |

## Example Papers


### psma-imaging-to-radioligand

*Diagnostic:* psma pet, psma imaging, psma scan | *Feature:* psma expression, psma-positive, psma positive | *Intervention:* 177lu-psma, lu-psma, psma radioligand

- **PMID 35750683** (2022, 242 cites) — unclassified — immunotherapy, synthetic-lethality — *Targeting signaling pathways in prostate cancer: mechanisms and clinical trials.*
- **PMID 31640747** (2019, 44 cites) — preclinical-invivo — untagged — *Targeting of prostate-specific membrane antigen for radio-ligand therapy of triple-negative breast cancer.*
- **PMID 38302933** (2024, 16 cites) — unclassified — immunotherapy, radioligand-therapy — *A multicentric, single arm, open-label, phase I/II study evaluating PSMA targeted radionuclide therapy in adult patients*

### pdl1-ihc-to-checkpoint

*Diagnostic:* pd-l1 immunohistochemistry, pd-l1 ihc, pd-l1 staining | *Feature:* pd-l1 positive, pd-l1 high, pd-l1 expression | *Intervention:* pembrolizumab, nivolumab, atezolizumab

- **PMID 29347993** (2018, 328 cites) — phase2-clinical — untagged — *Updated efficacy of avelumab in patients with previously treated metastatic Merkel cell carcinoma after ≥1 year of follo*
- **PMID 33203644** (2021, 118 cites) — phase2-clinical — epigenetic, immunotherapy — *Entinostat plus Pembrolizumab in Patients with Metastatic NSCLC Previously Treated with Anti-PD-(L)1 Therapy.*
- **PMID 39722028** (2024, 49 cites) — unclassified — antibody-drug-conjugate, car-t, immunotherapy — *Current and future immunotherapy for breast cancer.*

### tmb-msi-to-immunotherapy

*Diagnostic:* tumor mutational burden, tmb-high, tmb-h | *Feature:* tmb-high, tmb-h, msi-h | *Intervention:* pembrolizumab, nivolumab, checkpoint inhibitor

- **PMID 37492581** (2023, 80 cites) — unclassified — immunotherapy — *Immune escape and resistance to immunotherapy in mismatch repair deficient tumors.*
- **PMID 39123202** (2024, 72 cites) — unclassified — car-t, immunotherapy — *Recent developments in immunotherapy for gastrointestinal tract cancers.*
- **PMID 37655661** (2023, 57 cites) — unclassified — immunotherapy — *Linked CD4+/CD8+ T cell neoantigen vaccination overcomes immune checkpoint blockade resistance and enables tumor regress*

### neoantigen-profiling-to-mrna-vaccine

*Diagnostic:* neoantigen prediction, neoantigen profiling, neoantigen identification | *Feature:* neoantigen, neo-antigen, tumor-specific antigen | *Intervention:* mrna vaccine, mrna cancer vaccine, personalized vaccine

- **PMID 37165196** (2023, 1153 cites) — phase1-clinical — immunotherapy, mRNA-vaccine, nanoparticle — *Personalized RNA neoantigen vaccines stimulate T cells in pancreatic cancer.*
- **PMID 33473101** (2021, 257 cites) — clinical-other — immunotherapy, mRNA-vaccine — *Personalized neoantigen pulsed dendritic cell vaccine for advanced lung cancer.*
- **PMID 38584166** (2024, 189 cites) — phase2-clinical — immunotherapy, mRNA-vaccine — *Personalized neoantigen vaccine and pembrolizumab in advanced hepatocellular carcinoma: a phase 1/2 trial.*

### oncolytic-susceptibility-to-virotherapy

*Diagnostic:* viral receptor expression, nectin-1 expression, cd46 expression | *Feature:* viral entry receptor, nectin-1, cd46 | *Intervention:* t-vec, talimogene, oncolytic herpes

- **PMID 34205379** (2021, 13 cites) — preclinical-invivo — crispr, oncolytic-virus — *Nectin-1 Expression Correlates with the Susceptibility of Malignant Melanoma to Oncolytic Herpes Simplex Virus In Vitro *

### her2-testing-to-trastuzumab

*Diagnostic:* her2 ihc, her2 immunohistochemistry, her2 fish | *Feature:* her2-positive, her2 positive, her2 overexpression | *Intervention:* trastuzumab, herceptin, pertuzumab

- **PMID 37801674** (2023, 111 cites) — phase2-clinical — antibody-drug-conjugate — *Patritumab Deruxtecan (HER3-DXd), a Human Epidermal Growth Factor Receptor 3-Directed Antibody-Drug Conjugate, in Patien*
- **PMID 31451760** (2019, 82 cites) — preclinical-invitro — antibody-drug-conjugate, crispr — *CRISPR-Cas9 screens identify regulators of antibody-drug conjugate toxicity.*
- **PMID 38398191** (2024, 40 cites) — unclassified — antibody-drug-conjugate — *Next-Generation HER2-Targeted Antibody-Drug Conjugates in Breast Cancer.*

### brca-mutation-to-parp-inhibitor

*Diagnostic:* brca testing, brca1/2 testing, brca mutation testing | *Feature:* brca mutation, brca-mutant, brca-mutated | *Intervention:* olaparib, lynparza, niraparib

- **PMID 36082969** (2023, 446 cites) — phase3-clinical — synthetic-lethality — *Overall Survival With Maintenance Olaparib at a 7-Year Follow-Up in Patients With Newly Diagnosed Advanced Ovarian Cance*
- **PMID 29533782** (2018, 298 cites) — preclinical-invivo — synthetic-lethality — *BRD4 Inhibition Is Synthetic Lethal with PARP Inhibitors through the Induction of Homologous Recombination Deficiency.*
- **PMID 32122376** (2020, 244 cites) — unclassified — synthetic-lethality — *PARP inhibitors in pancreatic cancer: molecular mechanisms and clinical applications.*

### egfr-mutation-to-egfr-inhibitor

*Diagnostic:* egfr mutation testing, egfr testing, egfr mutation analysis | *Feature:* egfr mutation, egfr-mutant, egfr-mutated | *Intervention:* erlotinib, tarceva, gefitinib

- **PMID 21856766** (2011, 414 cites) — unclassified — untagged — *Disease flare after tyrosine kinase inhibitor discontinuation in patients with EGFR-mutant lung cancer and acquired resi*
- **PMID 34548332** (2022, 68 cites) — preclinical-invivo — antibody-drug-conjugate — *EGFR Inhibition Enhances the Cellular Uptake and Antitumor-Activity of the HER3 Antibody-Drug Conjugate HER3-DXd.*
- **PMID 37057110** (2023, 24 cites) — preclinical-invivo — untagged — *CKAP4 is a potential exosomal biomarker and therapeutic target for lung cancer.*

## Interpretation

- This pilot covers 10 diagnostic-therapy chains: the original 4-modality set (radioligands, checkpoint selection, mRNA vaccines, oncolytic viruses) plus the four most clinically-deployed predictive-biomarker-to-targeted-drug chains added in #441 (HER2-to-trastuzumab, BRCA-to-PARP-inhibitor, EGFR-to-EGFR-inhibitor, KRAS-G12C-to-sotorasib).
- Counts reflect this mechanism-keyword-built corpus, not a general-oncology corpus, so the targeted-therapy chains read far lower than their true clinical literature volume (EGFR and KRAS-G12C in particular), and a chain returning zero means the corpus lacks those papers, not that the chain is unimportant.
- The matching rule (intervention required + at least one other link) is conservative; papers that discuss only a diagnostic or only a therapy without the chain are excluded.
- Chain counts depend on keyword coverage and should not be read as exhaustive. Papers using non-standard terminology for diagnostics or interventions may be missed.