# CLAUDE.md

## Author

Ezequiel Lares

## What This Repo Is For

This is an open cancer-research workspace. The point is to help people—not to be right about one hypothesis.

The repo exists to compare therapeutic mechanisms, evidence depth, resistant-state biology, pathway targets, and simulation-based ideas honestly, so that anyone who reads it can form their own informed view. If the evidence says a direction is weak, say so. If a new direction looks promising, surface it. Don't protect old framing at the expense of clarity.

## Guiding Principles

1. **Let the evidence lead.** The repo started around a PDT/SDT-persister-ferroptosis thesis. That's still worth evaluating, but it's one lane among many. Don't treat it as the default winner.

2. **Stay open.** New modalities, targets, pathways, and interpretations should be welcomed when the corpus or external literature supports them. The README invites anyone with curiosity to contribute—the codebase should reflect that same openness.

3. **Be honest about what we don't know.** Corpus-level non-detection is provisional, not proof of absence. Missing landmark papers distort field-level claims. Taxonomy choices inflate or deflate conclusions. Say so directly rather than burying caveats.

4. **Make it reproducible.** Scripts should be re-runnable. Analysis outputs should be regenerated, not hand-edited. Separate generated files from handwritten interpretation notes so it's clear what came from the pipeline and what came from a person.

5. **Keep it human.** This project matters because cancer takes people from their families. Technical rigor serves that mission—but so does making the work accessible, welcoming contributions, and not hiding behind jargon when plain language works.

6. **The work is more important than the paper.** Don't optimize for journal word limits or publication constraints. If a finding needs context, give context. If a decision needs explaining, explain it. Clarity and contribution to the scientific community matter more than fitting a format. A reader should understand why we made each decision, what the results yield, and why we believe they signal value.

## Current Workstreams

- manuscript drafting and revision (Parts I–V substantially complete: 112 pages, 11 chapters + 3 appendices)
- corpus fetching, enrichment, tagging, and indexing
- evidence-tier audits and coverage caveats (gold-set evaluation: 46% exact, 96% precision, 55% recall)
- taxonomy and search refinement
- pathway-target and resistant-state analysis
- diagnostic-to-therapy chain extraction (6 chains, 129 articles mapped)
- tissue-of-origin analysis layer (5 tissue categories, 62% coverage)
- simulation work: ferroptosis biochemistry, drug penetration, calibration
- ferroptosis-core library packaging for external use
- news source authentication pipeline (fetch, extract claims, verify, score, index)
- broader strategy review of alternative therapies and biological bottlenecks
- operational maturity: CI, figure traceability, archival release (Phase 2 in progress)

## Current Repo State

- local full-text corpus: 4,830 records
- abstract-only archive: 5,584 records
- mechanism taxonomy, evidence tiers, pathway-targets, biology-process tags, and resistant-state scaffolding are all active
- evidence tagging is improved but still incomplete (gold-set measured)
- tissue-of-origin and weighted-evidence layers are active
- diagnostic-therapy matching layer covers 6 chains across 4 modalities (radioligand, checkpoint, mRNA vaccine, oncolytic)
- manuscript: 112 pages (book format), 11 chapters + 3 appendices, 19 figures, 33,344 words
- simulation suite: 10 binaries (incl. sim-tumor-pk) + ferroptosis-core library (MIT licensed, 10 modules, 31 unit tests) + Python bindings + 50 Python smoke tests
- news authentication pipeline: 5 scripts (fetch, extract claims, verify against PubMed, score credibility, build claim-centric index) implementing the 3-tier source framework from analysis/news-source-criteria.md
- simulation calibration: 5 targets documented, evaluate script operational
- drug penetration module: 3 tissue types, exponential Krogh approximation
- drug combination modeling: 4 drugs, pairwise Bliss synergy scoring with pathway traces
- tumor microenvironment: oxygen gradients, spatial immune zones (DAMP diffusion, T cell activation, anti-PD-1), LP overshoot multiplier, CAF-mediated stromal protection (GSH/MUFA supply), pH gradient (iron release + drug ion trapping)
- some landmark papers are known to be missing from the local full-text archive
- content provenance manifest (PROVENANCE.yaml) documents asset licensing and redistribution rights
- pinned Python environment (requirements-lock.txt, 32 packages) and Rust toolchain (rust-toolchain.toml, 1.92.0)
- contributor guide (CONTRIBUTING.md), citation metadata (CITATION.cff), and pytest in tracked dependencies

## What To Optimize For

- claims that are traceable and caveated
- taxonomy choices that do not inflate conclusions
- language that reflects uncertainty honestly
- outputs that help compare alternatives fairly
- maintainable scripts and reproducible reruns
- a tone that invites contribution rather than gatekeeping

## What To Avoid

- assuming the repo exists only to defend PDT/SDT
- writing stronger absence claims than the evidence model supports
- confusing patient-study signal with phase-labeled trial maturity
- treating broad process coverage as proof of therapeutic depth
- letting historical framing survive after the underlying analysis changed
- making the codebase feel closed or intimidating to newcomers

## Key Files

```text
README.md                                 repo-level framing and entry points
article/drafts/v1.{md,tex,pdf}            manuscript drafts
analysis/evidence-coverage-audit.md       evidence-tier coverage and guardrails
analysis/taxonomy-rerun-notes.md          taxonomy/query caveats after reruns
analysis/pathway-target-audit.md          pathway-target coverage
analysis/landmark-corpus-gaps.md          known missing papers that distort claims
corpus/INDEX.jsonl                        master index
scripts/                                  Python pipeline
simulations/                              Rust simulation work
simulations/ferroptosis-python/           Python bindings (PyO3)
simulations/ferroptosis-ffi/              C FFI bindings (PhysiCell integration)
tags/                                     precomputed tag indexes
article/book-outline.md                   frozen book outline and chapter contracts
article/AUTHORING.md                      writing rules and heading conventions
news/                                     news source scaffolding (issue #99)
PROVENANCE.yaml                           content provenance and redistribution rights
CONTRIBUTING.md                           contributor setup, testing, and PR guide
CITATION.cff                              citation metadata (renders GitHub "Cite" button)
requirements-lock.txt                     pinned Python dependency versions
```

## Search Conventions

Prefer fast repo-native inspection first:

```bash
rg "term" scripts analysis article
rg --files corpus/by-pmid | head
sed -n '1,120p' analysis/evidence-coverage-audit.md
```

## Writing Conventions

- every strong claim should be traceable to the corpus, analysis outputs, or external verification
- separate generated outputs from handwritten interpretation notes
- use coverage-aware language such as `not detected in the local keyword-derived analysis` where appropriate
- if a known taxonomy artifact or corpus gap applies, mention it directly rather than burying it
- keep the repo open to thesis revision rather than optimizing for rhetorical neatness
