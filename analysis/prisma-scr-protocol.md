# PRISMA-ScR protocol for the cancer-therapy scoping corpus (#505)

> From the 2026 fresh-eyes honesty review. This document reframes the corpus work
> honestly: it is **automated keyword scoping of open-access literature**, not a
> systematic review, and not (yet) a fully PRISMA-ScR-compliant scoping review.
> It documents the methodology against the PRISMA-ScR (Tricco et al. 2018,
> PMID 30178033) structure so a reader can see exactly which scoping-review
> elements are met and which are not.

## Honest scope statement (read this first)

The corpus is a **keyword-built convenience corpus** retrieved by automated
PubMed queries and filtered to records with locally resolvable open-access full
text. By the methodological standards of evidence synthesis it is **not** a
systematic review and does **not yet** meet every PRISMA-ScR requirement:

- There was **no a-priori registered protocol** at collection time. This document
  is **retrospective** methodology documentation, not a pre-registration. (The
  forward-looking, registrable predictions live in `PREREGISTRATION.md`; that is
  a different artifact and covers the simulation, not the corpus.)
- There was **no dual independent screening** and **no consensus adjudication**.
  Selection was a single automated pipeline.
- Eligibility was **"what the query returned and we could get OA full text for"**,
  not a human-applied inclusion/exclusion protocol.

Stating this is the point. The phrase "corpus-level non-detection" is only a
defensible scientific statement once these limits are explicit, and a
non-detection is reported as **"not detected in the local keyword-derived
analysis"**, never as proof of absence. Where this document uses "scoping", it
means scoping in the descriptive-mapping sense, with the compliance gaps above
named, not a claim of a completed PRISMA-ScR review.

## Why a scoping review, not a systematic review

PRISMA-ScR is the right target frame because the goal is to **map** the extent,
range, and concentration of cancer-therapy research across mechanisms and cancer
types and to surface gaps and biases, not to answer a focused effectiveness
question or pool effect sizes. A systematic review would require a registered
protocol, dual screening, risk-of-bias assessment, and (often) meta-analysis,
none of which this automated pipeline performs.

---

## 1. Research question (PCC framing)

Scoping reviews use Population / Concept / Context rather than PICO.

- **Population:** cancer (any type; 22 cancer types are charted).
- **Concept:** therapeutic mechanisms and modalities (19 frozen mechanism tags,
  spanning molecular targets, drug classes, delivery platforms, device-based
  interventions, and biological strategies), together with evidence maturity and
  pathway/resistance annotations.
- **Context:** the indexed, predominantly open-access biomedical literature
  (PubMed-retrievable, OA-resolvable full text), 2001 to 2026.

**Primary question.** Across cancer-therapy mechanisms, how is the
PubMed-indexed, open-access-resolvable literature distributed by mechanism,
cancer type, and evidence tier, and which apparent gaps are attributable to
search design, taxonomy granularity, or open-access skew rather than to genuine
research absence?

---

## 2. Eligibility criteria (as actually applied)

| | Criterion |
| --- | --- |
| **Include** | Record returned by one of the 19 mechanism-specific PubMed queries AND addressing at least one of the 19 mechanisms in a cancer-treatment or cancer-biology context. |
| **Full-text subset** | Of included records, those with locally resolvable OA full text form the 4,830-record quantitative corpus. |
| **Abstract-only archive** | Included records without resolvable full text (5,585) are retained separately and used only for the open-access sensitivity check, never for full-text-derived statistics. |
| **Exclude** | No additional manual inclusion/exclusion beyond the query constraints. Records not returned by the queries are out of scope by construction. |
| **Languages** | English-biased (PubMed query + OA-resolution), a documented limitation. |
| **Publication types** | Original research, trial reports, reviews, systematic reviews, narrative reviews, and some case series are all retained; review-like articles are left unclassified for evidence tier rather than forced into one. |

---

## 3. Information sources and search strategy

- **Source:** PubMed via the NCBI E-utilities API, indexed through March 2026.
- **Enrichment:** OpenAlex (open-access status, citations, topics); PubMed Central
  and publisher OA endpoints for full text.
- **Search strategy:** the 19 mechanism-specific query sets are committed verbatim
  in [`scripts/queries.txt`](../scripts/queries.txt) (reproducibility: the exact
  strings are the search record).
- **Documented date-window asymmetry.** Newer/novel mechanisms (CRISPR, CAR-T,
  ADCs, bispecifics, mRNA vaccines, synthetic lethality, oncolytic viruses,
  epigenetic therapy, targeted protein degradation, radioligand therapy) were
  retrieved with a recent lower bound (2020/2022/2023 through 2026); historically
  broader physical, frequency, and immunotherapy modalities were retrieved with no
  date window. Every query was fetched with a per-query cap of 500 records.
- **Consequence (stated, not hidden):** cross-mechanism volume counts are a
  **within-fetch comparison, not a field census**. This is disclosed in the
  manuscript Section 3.3.1 ("Query-design bias") and noted in `queries.txt`.

---

## 4. Selection process (as actually performed)

A **single automated pipeline** (fetch, OA-resolve, enrich, keyword-tag, index),
**no dual independent screening, no consensus step**. This is the largest
departure from a compliant PRISMA-ScR selection process and is the reason the
work is labeled automated keyword scoping rather than a scoping review.
Reproducibility substitutes for inter-rater reliability here: the frozen
`corpus/INDEX.jsonl`, the pinned environment, and CI make the selection
deterministic and re-runnable, but determinism is not the same as independent
human screening, and we do not claim it is.

---

## 5. Data-charting method

Each record is charted along three axes by automated keyword matching against
curated dictionaries (`scripts/config.py`):

1. **Therapeutic mechanism(s):** one or more of 19 tags (multi-tag allowed,
   enabling convergence analysis).
2. **Cancer type(s):** one or more of 22 types.
3. **Evidence tier:** one of seven tiers (Phase III/RCT, Phase II, Phase I,
   clinical-other, preclinical in vivo, preclinical in vitro,
   theoretical/computational); review-like and protocol articles left
   unclassified.

Additional charted layers: pathway targets, tissue-of-origin, weighted evidence,
and diagnostic-to-therapy chains (recomputed on the fly from frozen text).

---

## 6. PRISMA-ScR flow (tied to the real counts)

This reproduces the flow already shown as manuscript Figure 20 (the #135
PRISMA-corpus-flow work), with the committed counts:

```
Identification
  PubMed records returned by 19 mechanism queries .......... 10,415 unique
        |
Full-text resolution
  OA full text locally resolvable ......................... 4,830 (46.4%)
  Abstract-only (retained in separate archive) ............ 5,585 (53.6%)
        |
Charting (full-text subset)
  Evidence-tier tagged .................................... 2,038 (42.2% of 4,830)
  Indexed but not evidence-classified ..................... 2,792
        |
Included for quantitative analysis ....................... 4,830 full-text
  (803 journals; 19 mechanisms; 22 cancer types; 2001-2026)
```

No records were excluded by manual screening; the only "exclusion" is the
full-text-resolution step, which is an open-access access filter, not an
eligibility judgment, and is itself quantified as the open-access bias in
Section 3.3.1.

---

## 7. Two recall types (never conflate them)

The corpus reports **two different recall numbers** that measure different
things. They are not interchangeable and neither stands in for the other:

| Recall type | Value | What it measures | Source |
| --- | --- | --- | --- |
| **Evidence-tier recall** | **55%** | Of records that truly carry an evidence-bearing tier, the fraction the tagger labels (with 96% binary evidence-presence precision). Measured on the 100-article gold set; this is the **only** recall that sample can measure. | Manuscript Sec. 3.6, Figure 9 |
| **Mechanism-presence recall** | **~90%** (90.6% volume-weighted / 89.2% macro) | Of expert-MeSH-labelled mechanism records, the fraction the tagger recovers from title+abstract alone, measured non-circularly against independent MeSH leaves. | `mechanism-recall-report.md` (#412) |

**Implication for absence claims.** Absence claims rest on **mechanism
detection**, which is high (~90%) where it can be checked without circularity, so
for those mechanisms a non-detection carries more weight. For modalities with no
discriminative MeSH descriptor (TTFields, electrolysis, bioelectric, cold
atmospheric plasma, cuproptosis, disulfidptosis, targeted protein degradation,
radioligand therapy), recall is **not MeSH-measurable** and the provisional-absence
caution still applies. The 55% evidence-tier recall constrains how strongly
**evidence-maturity** claims can be made, a separate axis from mechanism presence.

---

## 8. PRISMA-ScR 20-item checklist: applicability

Honest item-by-item status (met / partial / not met). The point is transparency
about the gaps, not a compliance badge.

| PRISMA-ScR item | Status | Note |
| --- | --- | --- |
| 1-2 Title/abstract as scoping review | Partial | Reframed to "automated keyword scoping"; manuscript language corrected (Sec. 3.2). |
| 3-4 Rationale, objectives (PCC) | Met | Section 1 above. |
| 5 Protocol/registration | **Not met** | No a-priori protocol; this doc is retrospective. |
| 6 Eligibility criteria | Partial | Stated (Sec. 2) but query-driven, not human-applied. |
| 7 Information sources | Met | PubMed + OpenAlex + OA endpoints, dated. |
| 8 Search | Met | Exact queries committed in `queries.txt`. |
| 9 Selection of sources | Partial | Automated, **no dual screening** (Sec. 4). |
| 10 Data charting | Met | Three axes + layers (Sec. 5). |
| 11 Data items | Met | Defined dictionaries (`config.py`). |
| 12 Critical appraisal | **Not met** | No per-source risk-of-bias appraisal (consistent with scoping-review norms, which make appraisal optional). |
| 13-16 Results / flow | Met | Flow (Sec. 6); bias quantification (Sec. 3.3.1). |
| 17-18 Summary, limitations | Met | Manuscript Sec. 3.6; this document. |
| 19-20 Conclusions, funding | Met | Open repo; funding/competing interests in CITATION.cff. |

**Net:** the work is a transparent, reproducible **keyword-scoping map** that
meets most charting/reporting items but fails the registration, dual-screening,
and critical-appraisal items required to claim a completed systematic review, and
falls short of full PRISMA-ScR on the registration and dual-screening items. We
report it as exactly that.

## Acceptance-criteria check (#505)

- [x] PRISMA-ScR protocol: research question (PCC), eligibility, search strategy
  (with documented date windows), selection process, data-charting method.
- [x] PRISMA flow tied to the real counts (Sec. 6, built on the Figure 20 / #135
  flow).
- [x] Activity renamed honestly as "automated keyword scoping of OA literature";
  the lone "systematic literature search" over-claim in the manuscript (Sec. 3.2)
  is corrected to match.
- [x] Mechanism-presence recall (~90%, #412) and evidence-tier recall (55%) are
  distinguished and never substituted (Sec. 7).
