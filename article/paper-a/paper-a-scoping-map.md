# A reproducible open-access keyword-scoping map of the cancer-therapy literature, with quantified search, taxonomy, and access biases

**Authors:** Ezequiel Lares (lead); _open for a co-author with evidence-synthesis / scientometrics expertise — see [`analysis/collaborator-brief.md`](../../analysis/collaborator-brief.md)._

> **STATUS: SKELETON — not a submission.** This is the Paper A scaffold called for by #523, operationalizing the two-paper split outline ([`analysis/manuscript-split-outline.md`](../../analysis/manuscript-split-outline.md), #504). Section bodies are content notes + source pointers into the existing manuscript and analysis outputs, not finished prose. **Venue and length targets below are provisional placeholders the maintainer can change.**

**Provisional target:** a research-on-research / evidence-synthesis methods venue, a scientometrics venue, or a cancer-informatics venue that publishes reproducible evidence-mapping methodology. Preprint first (a general biomedical server) so the map is citable immediately.

**Length target:** standard research/methods article, ~4,000–6,000 words main text.

---

## Lead claim (the one sentence this paper defends)

Several apparent "gaps" in the cancer-therapy literature are artifacts of **search design, taxonomy granularity, and open-access skew — not of biology** — and a fully reproducible automated keyword-scoping pipeline can *quantify each of those biases* rather than hide them.

## What this paper is NOT

Not a systematic review (no a-priori registered protocol, no dual independent screening), not a clinical-claims paper, and not a census of absolute literature sizes. It is a **methods + landscape** paper. Methodology framing inherited from [`analysis/prisma-scr-protocol.md`](../../analysis/prisma-scr-protocol.md) (#505).

---

## Abstract (draft stub, ~200 words — replace)

> Cancer-therapy literatures mature in parallel and are rarely compared on shared axes (evidence depth, resistant-state biology, pathway dependencies, delivery constraints). We present a fully reproducible, automated keyword-scoping pipeline that maps 4,830 full-text open-access cancer-therapy articles (plus a 5,585-record abstract-only archive) across 19 mechanisms and 22 cancer types, and — crucially — quantifies its own biases rather than presenting the map as ground truth. We report the tagger's measured quality with two distinct, never-conflated recall numbers (evidence-tier: 96% binary evidence-presence precision, 55% recall on a 100-article gold set; mechanism-presence: ~90% measured non-circularly against independent MeSH leaves), and we show that three classes of bias — open-access skew (98.7% OA in the full-text corpus), date-window/fetch-cap query design, and taxonomy granularity — each move the apparent rankings in a measurable, documented direction. We demonstrate that several "zero-publication gaps" collapse under coarser taxonomy groupings, and that mechanism×cancer concentration is better read as hypergeometric over-representation (BH-FDR) than as raw counts. The contribution is the reproducible, bias-aware method and the honest accounting it enables; absence is reported throughout as "not detected in the local keyword-derived analysis," never as evidence of absence. All code, queries, the frozen index, and the pinned environment are released for re-running and extension.

---

## Section outline (skeleton)

### 1. Introduction
- Cancer-therapy literatures mature in parallel; rarely compared on shared axes.
- The contribution: a reproducible map **plus** an honest accounting of what it can and cannot say.
- _Source: v1.md §1 Introduction (de-claimed: "cross-literature analysis," not "systematic")._

### 2. Methods — the reproducible pipeline
- Search strategy: 19 mechanism-specific PubMed queries (`scripts/queries.txt`); state the date-window asymmetry + 500-record fetch cap **up front** (#510).
- Retrieval + full-text resolution (PubMed Central + publisher OA endpoints); OpenAlex enrichment.
- Automated keyword tagging on three axes (mechanism / cancer type / evidence tier); **no manual article-level screening**, stated plainly.
- Reproducibility: frozen `corpus/INDEX.jsonl`, pinned `requirements-lock.txt`, CI, release manifest.
- _Source: v1.md §3.1–3.2; `scripts/{fetch_articles,tag_articles,build_index}.py`._

### 3. Tagger validation — the honesty core
- Evidence-tier tagger: 46% exact-label, **96% binary evidence-presence precision, 55% recall** (100-article gold set); per-tier precision corrected (#509).
- Mechanism-presence recall measured **separately and non-circularly** vs independent MeSH leaves: 90.6% volume-weighted / 89.2% macro (#412).
- The two recall numbers measure different things and are **never substituted** for one another.
- _Source: `analysis/evidence-coverage-audit.md`, `analysis/mechanism-recall-report.md`._

### 4. The map — descriptive results
- Mechanism×cancer matrix (now with **hypergeometric over-representation + BH-FDR**, #525, ranking by fold over expected rather than raw count), convergence map, evidence-tier distribution, growth trajectories.
- All framed descriptively; absence = "not detected in the local keyword-derived analysis."
- _Source: `analysis/mechanism-matrix.md`, `analysis/convergence-map.md`, `analysis/evidence-tiers.md`._

### 5. Bias quantification — the methodological payload
- **Open-access skew:** 98.7% OA full-text vs 29.1% OA abstract-only; immunotherapy share 34.4%→28.7%, physical class 14.7%→22.4% (`scripts/oa_bias_analysis.py`, #348).
- **Query-design bias:** date-capped + fetch-capped ⇒ within-fetch comparison only (#510).
- **Taxonomy-granularity bias:** zero-publication gaps 94 → 29–38 under collapsed groupings; the three "must-survive" conclusions hold (#133).
- _Source: `analysis/oa-bias-report.md`, `analysis/taxonomy-rerun-notes.md`._

### 6. Discussion
- What a reproducible map buys (citable, re-runnable, bias-aware) and what it cannot replace (a registered systematic review with dual screening).
- Living-review path (`scripts/living_review_update.py`, #349) for keeping the map current.

### 7. Supplementary
- Full per-mechanism OA tables; the diagnostic-to-therapy chains (ten chains, #441); the taxonomy-sensitivity preregistration.

## Figures (provisional)
- F1: Cross-literature landscape (mechanism×cancer, FDR-enriched).
- F2: Evidence-tier distribution per mechanism.
- F3: OA-skew ranking shift (full-text vs abstract-only).
- F4: Taxonomy-granularity gap-count sensitivity.

## What moves OUT of Paper A
All simulation content (→ Paper B). The corpus paper stands entirely on the pipeline + bias analyses + validation numbers.

## Source map (book → this paper)
| v1.md content | Section here |
|---|---|
| §3 corpus construction + bias analyses | §2, §5 |
| §3.6 tagger validation + #412 | §3 |
| §4 mechanism landscape / rankings | §4 |
| Diagnostic-to-therapy chains | §7 supplementary |
