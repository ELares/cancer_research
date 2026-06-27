# Manuscript split: two focused papers (#504)

> From the 2026 fresh-eyes honesty review. This is a **plan and a pair of
> outlines, not a submission**. Venue names are candidate targets with rationale,
> not acceptances. It operationalizes the venue split already argued in
> [`contribution-plan-2026.md`](contribution-plan-2026.md) Part 1 into two
> concrete, section-by-section outlines, each leading with the single claim it
> will defend.

## Why split (the one-paragraph version)

The current draft (`article/drafts/v1.md`) is a ~115 page, ~45,800 word, 11
chapter book that repeatedly self-describes as "this book", while its own header
targets a 3,000 to 5,000 word Perspective. That length/target mismatch is a desk
reject before the content is read. Worse, the book bundles two genuinely
different contributions (a literature scoping map and a simulation engine) whose
audiences, review criteria, and failure modes do not overlap. Bundling lets a
weakness in the uncalibrated simulation undermine the corpus method, and invites
reviewers to judge the simulation as if it claimed calibrated predictions, which
it does not. The fix is two short papers, each with a defended contribution, and
the book demoted to background/supplementary material.

The hard rule for both: **lead with the claim the paper defends. Do not staple
two half-papers together.** Each paper must survive its own reviewers without
leaning on the other.

---

## Paper A: a reproducible open-access keyword-scoping map of the cancer-therapy literature

**Lead claim (the one sentence the paper defends).** Several apparent "gaps" in
the cancer-therapy literature are artifacts of search design, taxonomy
granularity, and open-access skew, not of biology, and an automated, fully
reproducible keyword-scoping pipeline can quantify each of those biases rather
than hide them.

**What this paper is NOT.** It is not a systematic review (no a-priori protocol,
no dual independent screening), not a clinical-claims paper, and not a
field-census of absolute literature sizes. It is a methods + landscape paper.
See `prisma-scr-protocol.md` for the honest methodology framing it inherits.

**Target venue (candidates, by fit).** A research-on-research / evidence-synthesis
methods venue, a scientometrics venue, or a cancer-informatics venue that
publishes reproducible evidence-mapping methodology. Preprint first
(general biomedical server) so the map is citable immediately.

**Length target.** A standard research/methods article (roughly 4,000 to 6,000
words main text), not a book.

### Section outline

1. **Introduction.** Cancer-therapy literatures mature in parallel and are rarely
   compared on shared axes. The contribution is a reproducible map plus an honest
   accounting of what the map can and cannot say.
2. **Methods (the reproducible pipeline).**
   - Search strategy: 19 mechanism-specific PubMed queries (`scripts/queries.txt`),
     with the documented date-window asymmetry stated up front.
   - Retrieval and full-text resolution (PubMed Central + publisher OA endpoints);
     OpenAlex enrichment.
   - Automated keyword tagging along three axes (mechanism, cancer type, evidence
     tier). No manual article-level screening, stated plainly.
   - Reproducibility: frozen `corpus/INDEX.jsonl`, pinned environment, CI.
3. **Tagger validation (the honesty core).**
   - Evidence-tier tagger: 46% exact-label accuracy, 96% binary
     evidence-presence precision, 55% evidence-tier recall (100-article gold set).
   - Mechanism-presence recall measured separately and non-circularly against
     independent MeSH leaves: 90.6% volume-weighted / 89.2% macro (#412).
   - The two recall numbers measure different things and are never substituted
     for one another.
4. **The map (descriptive results).** Mechanism x cancer matrix, convergence map,
   evidence-tier distribution, growth trajectories. All framed descriptively;
   absence reported as "not detected in the local keyword-derived analysis".
5. **Bias quantification (the methodological payload).**
   - Open-access skew (98.7% OA full-text vs 29.1% OA abstract-only; immunotherapy
     share 34.4% to 28.7%; physical class 14.7% to 22.4%), from
     `oa_bias_analysis.py`.
   - Query-design bias (date-capped + 500-record cap = within-fetch comparison).
   - Taxonomy-granularity bias (zero-publication gaps 94 to 29-38 under collapsed
     groupings; the three "must survive" conclusions hold).
6. **Discussion.** What a reproducible map buys (citable, re-runnable, bias-aware)
   and what it cannot replace (a registered systematic review with dual
   screening). Living-review path (`living_review_update.py`).
7. **Supplementary.** Full per-mechanism OA tables, the diagnostic-to-therapy
   chains, the taxonomy-sensitivity preregistration.

### What moves OUT of Paper A

All simulation content. The corpus paper stands entirely on the pipeline + the
bias analyses + the validation numbers.

---

## Paper B: an embeddable ferroptosis simulation engine with an honest, data-anchored validation core

**Lead claim (the one sentence the paper defends).** A small, reusable,
byte-identity-disciplined ferroptosis biochemistry engine can reproduce several
**independently measured** ferroptosis phenomena (an in-vitro GPX4-inhibitor
dose-response, spheroid zone geometry, a paired plasma/tumor PK course, a
drug-penetration length, and a propagating ferroptotic trigger-wave speed), and
the engineering discipline that keeps those anchored results honest (a model
card, a per-layer calibration ledger, a calibration-regression CI gate, and a
byte-identity containment gate) is itself the contribution.

**What this paper is NOT.** It is not a paper of validated quantitative clinical
predictions. The dozens of off-by-default tumor-microenvironment realism layers
are explicitly NOT the headline. They are a single supplementary catalogue,
labeled directional / uncalibrated, used in zero reported numbers.

**Target venue (candidates, by fit).** A computational / systems-biology methods
venue that publishes reusable tools with explicit validation-status reporting.
Preprint first (quantitative-biology category).

**Length target.** A methods article (roughly 4,000 to 6,000 words main text)
plus a software/figure supplement.

### Section outline

1. **Introduction.** A reusable ferroptosis engine is useful only if its
   validation status is legible. The contribution is the tool plus the honesty
   infrastructure, not a set of clinical numbers.
2. **Engine design.** The single-cell biochemical core (the bistable
   recover-or-collapse switch), the FFI/Python bindings, and the off-by-default
   byte-identical discipline (why every layer ships inert).
3. **The data-anchored validation core (the headline).** Lead with ONLY the legs
   tied to independent data, each with its held-out metric:
   - **#330** in-vitro kill switch vs CTRPv2 GPX4-inhibitor dose-response
     (ML162 fit RMSE 0.05, ML210 held-out RMSE 0.07).
   - **#333** spheroid zone geometry vs Browning 2021 confocal structure
     (size-aware thresholds, mean boundary error 0.025 vs 0.338 fixed).
   - **#334** tumor PK vs imidazole-ketone-erastin paired plasma/tumor course
     (partition Kp 0.90; the multi-compartment finding).
   - **#335** Krogh drug-penetration form + reference length vs Primeau/Tannock
     (lambda 50 um within the measured 25 to 75 um).
   - **#464** PDT/SDT source-independent singlet-oxygen kill threshold vs
     Zhu 2015 (~0.5 mM).
   - **#482** ferroptotic trigger-wave speed vs Co 2024 (5.52 um/min baseline,
     c proportional to sqrt(iron), reproducing the DFO/loading 2.33/9.40).
   - Cross-validation against published ODE models (#344) and the Arbatskiy
     73-ODE structural cross-check (#471).
4. **Identifiability and uncertainty.** State plainly that no headline magnitude
   is point-estimable under the documented parameter uncertainty
   (`identifiability-report.md`): 11 swept parameters, 6 non-identifiable from
   the kill rate, prior-predictive intervals, the in-vitro/in-vivo ABC
   disjunction.
5. **Honesty infrastructure (the second half of the contribution).** The model
   card, `CALIBRATION_STATUS.md` per-layer ledger, the calibrate-or-cut policy
   (#501), the calibration-regression CI gate (#499), and the byte-identity
   containment gate (#253) as a reusable pattern for "broad but honest" modeling.
6. **Discussion.** What an anchored-core + uncalibrated-periphery tool is good
   for (hypothesis generation, falsifiable predictions) and what it cannot do
   (stand in for in-vivo dose-response).
7. **Supplementary catalogue.** ONE table listing every off-by-default realism
   layer with its tier, its named calibration target, and its
   used-in-any-reported-number = N status. This is where the layer catalogue
   lives, demoted from the body.

### What moves OUT of Paper B

All corpus/literature content. The uncalibrated layer catalogue moves to a single
supplementary table. The book's narrative chapters become background.

---

## How the book maps onto the two papers

| Book content | Destination |
| --- | --- |
| Ch. 3 corpus construction + bias analyses | Paper A body |
| Ch. 3.6 tagger validation + #412 mechanism recall | Paper A body (honesty core) |
| Ch. 4 mechanism landscape / rankings | Paper A body (descriptive) |
| Diagnostic-to-therapy chains | Paper A supplementary |
| Ch. 5 single-cell switch + identifiability | Paper B body (anchored core + Sec. 4) |
| Ch. 6-7 drug penetration, combinations, TME layers | Paper B: anchored legs to body, the rest to the supplementary catalogue |
| Ch. 8.4 structural limitations | Paper B discussion |
| The 30+ off-by-default realism layers | Paper B: ONE supplementary table |
| Narrative framing ("this book", mission chapters) | Background; not in either submission |

## Acceptance-criteria check (#504)

- [x] Corpus/resource paper outline with a defended contribution + candidate venues.
- [x] Simulation methods paper outline built around the data-anchored legs, with
  the uncalibrated catalogue demoted to a single supplementary table.
- [x] Each piece states an explicit lead claim at an appropriate length; the book
  becomes background.

Venue choice and submission are human actions (the issue assigns them to the
author). This document is the draftable artifact the issue asks for.

**Skeletons delivered (#523).** The two paper scaffolds built from this outline now
live at [`article/paper-a/paper-a-scoping-map.md`](../article/paper-a/paper-a-scoping-map.md)
and [`article/paper-b/paper-b-ferroptosis-engine.md`](../article/paper-b/paper-b-ferroptosis-engine.md)
— each leads with its defended claim, carries a draft abstract stub, a
section-by-section skeleton with source pointers into `v1.md` and the analysis
outputs, and a provisional (changeable) venue + length target. Filling the section
bodies and choosing the final venues remain the author's steps.
