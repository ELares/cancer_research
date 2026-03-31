# Cancer Research Synthesis and Analysis

## Why This Repo Exists

Cancer research is too important to leave trapped in disconnected literatures, siloed tooling, or narrow framing.

This repository is an open attempt to map, compare, and stress-test therapeutic ideas across a large cancer-research corpus using literature analysis, tagging pipelines, generated audits, and simulation work. It started from a specific hypothesis around PDT/SDT, ferroptosis, and drug-tolerant persister cells, but the repo has grown beyond that one thesis.

The current goal is not to prove a predetermined answer. It is to build a transparent workspace where hypotheses can be tested, weakened, expanded, or replaced when the corpus and external evidence demand it.

## Current Framing

This repo currently supports several linked but distinct workstreams:

- cross-literature mapping of therapeutic mechanisms across cancer types
- evidence-tier analysis with explicit coverage caveats
- taxonomy and query refinement so field-level claims are less artifact-prone
- pathway-target and resistant-state tracking for recurrence and therapy escape
- manuscript drafting and revision
- simulation work around ferroptosis-related dynamics and escape mechanisms

The PDT/SDT-persister-ferroptosis hypothesis is still part of the project, but it is now treated as one candidate research direction inside a broader resistant-state and cross-modality analysis effort.

## What The Repo Contains

- a local corpus of full-text and abstract-only cancer research records
- Python scripts for fetching, enriching, tagging, indexing, and analyzing the corpus
- generated analysis files that summarize mechanism coverage, evidence tiers, pathway-target signals, and known gaps
- manuscript drafts and figures
- Rust simulations for ferroptosis-related modeling

## Current Priorities

The repo is currently oriented around questions like:

- are our conclusions stable under better taxonomy and evidence tagging?
- where are the biggest corpus-coverage blind spots?
- which signals are real versus query artifacts?
- are we over-centered on PDT/SDT relative to the broader field?
- which alternative modalities, targets, or resistant states deserve equal or greater attention?

## Known Limits

Several important constraints should shape how this repo is used:

- the corpus is large, but not complete
- some field-defining papers are missing from the local full-text archive
- evidence tagging is heuristic and coverage-aware, not definitive
- taxonomy choices materially affect gap counts and field comparisons
- resistant-state tagging is intentionally conservative and still sparse
- simulation outputs are useful for hypothesis support, not experimental proof

In practice, repo outputs should be read as structured non-detection, prioritization support, and hypothesis scaffolding unless externally verified.

## Key Entry Points

- [article/drafts/v1.md](article/drafts/v1.md)
  Current manuscript draft.
- [analysis/evidence-coverage-audit.md](analysis/evidence-coverage-audit.md)
  Evidence-tag coverage and interpretation guardrails.
- [analysis/taxonomy-rerun-notes.md](analysis/taxonomy-rerun-notes.md)
  What changed after taxonomy tightening and reruns.
- [analysis/pathway-target-audit.md](analysis/pathway-target-audit.md)
  First-pass pathway-target layer for ferroptosis resistance and adjacent programs.
- [analysis/landmark-corpus-gaps.md](analysis/landmark-corpus-gaps.md)
  Known missing papers large enough to distort claims.
- [scripts/](scripts/)
  Fetching, enrichment, tagging, indexing, and analysis pipeline.
- [simulations/](simulations/)
  Rust simulation code.

## Repository Structure

```text
analysis/                    generated audits, maps, gap notes, and interpretation docs
article/drafts/              manuscript drafts
article/figures/             figure outputs
article/references/          bibliography
corpus/by-pmid/              local full-text corpus
corpus/abstracts/by-pmid/    abstract-only archive
scripts/                     Python pipeline
simulations/                 Rust simulation work
tags/                        precomputed tag indexes
plans/                       planning notes
books/                       reference texts
```

## Reproduction

```bash
pip install -r requirements.txt
cp .env.example .env

python scripts/tag_articles.py
python scripts/build_index.py
python scripts/analyze_corpus.py
python scripts/generate_figures.py
```

Simulation examples:

```bash
cd simulations
cargo build --release
cargo run --release -p sim-original
cargo run --release -p sim-spatial
cargo run --release -p sim-window
cargo run --release -p sim-icd
cargo run --release -p sim-combo
```

## Contributing

Open an issue or pull request if you want to:

- challenge a conclusion
- tighten the taxonomy
- add missing landmark papers
- broaden the mechanism or pathway coverage
- improve the evidence model
- revise the manuscript
- extend the simulations or charts

The repo is most useful when it stays falsifiable and open to revision.

## License

MIT License. See [LICENSE](LICENSE).
