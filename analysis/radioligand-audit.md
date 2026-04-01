# Radioligand Lane Audit

Audit note for the cleaned `radioligand-therapy` mechanism after removing generic theranostic spillover and adding a minimal target-level layer.

Current `radioligand-therapy` full-text count: **11**.

## Evidence Mix

- Tagged evidence records: 8
- Review-like records: 1
- Protocol-like records: 1
- Other untagged primary-study-like records: 1

## Target-Level Distinctions

- **psma**: 3 articles
- **fap**: 1 articles
- **sstr**: 1 articles

## Audited Former False Positives

These PMIDs were previously strong contamination candidates because generic theranostic language could bridge into the radioligand lane.

- **PMID 25728459**: removed from radioligand lane — *A review of low-intensity ultrasound for cancer therapy.*
- **PMID 31410214**: removed from radioligand lane — *Focused ultrasound-augmented targeting delivery of nanosonosensitizers from homogenous exosomes for enhanced sonodynamic cancer therapy.*
- **PMID 30613291**: removed from radioligand lane — *Nanosonosensitizers for Highly Efficient Sonodynamic Cancer Theranostics.*
- **PMID 40321808**: removed from radioligand lane — *A Bi2O3-TiO2 Heterojunction for Triple-Modality Cancer Theranostics.*

## Representative Retained Positives

- **PMID 41342316** (2026) — fap, psma, sstr — `review_like` — *Advancements in Targeted Radiopharmaceuticals: Innovations in Diagnosis and Therapy for Enhanced Cancer Management.*
- **PMID 38302933** (2024) — psma — `protocol_like` — *A multicentric, single arm, open-label, phase I/II study evaluating PSMA targeted radionuclide therapy in adult patients with metastatic clear cell re*
- **PMID 38446353** (2024) — psma — `clinical-other` — *Lutetium-177 Labelled Anti-PSMA Monoclonal Antibody (Lu-TLX591) Therapy for Metastatic Prostate Cancer: Treatment Toxicity and Outcomes.*
- **PMID 30911535** (2018) — target-unspecified — `other_untagged` — *Thyroid Cancer Radiotheragnostics: the case for activity adjusted 131I therapy.*
- **PMID 29468134** (2018) — target-unspecified — `clinical-other` — *Yttrium-90 microsphere selective internal radiation therapy for liver metastases following systemic chemotherapy and surgical resection for metastatic*

## Interpretation

- Generic `theranostic` phrasing is no longer sufficient by itself to create a radioligand hit. The lane now requires radionuclide-specific therapy signals or a target-plus-radionuclide pattern.
- The cleaned lane is smaller but more defensible, and it is now usable for target-level questions such as whether PSMA, FAP, or SSTR dominate the accessible local full-text archive.
- The lane is still constrained by corpus coverage. The missing VISION trial remains a known archive artifact and still limits how strong any absence claim should be.