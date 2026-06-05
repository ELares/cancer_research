# Landmark Corpus Gaps

This note records a small set of clinically important papers that are missing from the local full-text corpus and are large enough to distort gap claims if ignored.

A reusable guardrail, `scripts/landmark_coverage.py` (#345), checks each curated landmark's membership in the corpus and writes `analysis/landmark-coverage-report.md`; run it with `--recover-missing` to fetch any absent landmark's abstract from PubMed and add it as a tagged abstract record. All five PMIDs below were verified against PubMed (June 2026); VISION (34161051) was recovered via that script. Recovered records are ABSTRACT-ONLY and do NOT change the frozen full-text quantitative results (which analyze the 4,830 full-text records); they fix the coverage guardrail so absence claims are not corpus artifacts.

## Confirmed gaps

| Mechanism | PMID | Current local status | Why it matters |
|---|---|---|---|
| radioligand-therapy | 34161051 | **recovered to `corpus/abstracts/by-pmid/` (#345)**; still missing from full-text `corpus/by-pmid/` | VISION is a field-defining trial for `177Lu-PSMA-617`; its absence made `0 phase-labeled trial evidence detected` for radioligands a corpus artifact. Now a tagged abstract record (radioligand-therapy / prostate / phase3-clinical, iCite RCR 143). |
| ttfields | 40448572 | present in `corpus/abstracts/by-pmid/`, missing from `corpus/by-pmid/` | PANOVA-3 is a pivotal phase III TTFields study in locally advanced pancreatic cancer and materially affects how mature TTFields looks outside glioblastoma. |
| mRNA-vaccine | 33016924 | present in `corpus/abstracts/by-pmid/`, missing from `corpus/by-pmid/` | Phase I/II gastrointestinal neoantigen-vaccine study with strong citation impact; omission weakens the apparent patient-level depth of the vaccine field. |
| mRNA-vaccine | 36027916 | present in `corpus/abstracts/by-pmid/`, missing from `corpus/by-pmid/` | Phase I NEO-PV-01 plus chemotherapy and anti-PD-1 in NSCLC; important for combination-vaccine framing in lung cancer. |
| mRNA-vaccine | 35970920 | present in `corpus/abstracts/by-pmid/`, missing from `corpus/by-pmid/` | Phase I/II individualized heterologous adenovirus plus self-amplifying mRNA neoantigen vaccine study; relevant to platform breadth, not just one vaccine format. |

## How to use this list

- Do not treat this as a comprehensive audit of all missing important papers.
- Use it as a guardrail when writing mechanism-level absence claims from the full-text corpus.
- If a mechanism has a known gap on this list, prefer `not detected in the local full-text corpus` over stronger wording.
- If a manuscript paragraph depends on one of these areas, verify against PubMed or the trial paper directly before using the corpus count as evidence of immaturity.
