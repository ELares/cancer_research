# Beyond Pharmacologic Ferroptosis Inducers: Physical ROS Modalities for Drug-Tolerant Persister Cells

## Why This Exists

Too often has cancer taken from us the people we love.

There was a point in my life where I volunteered at a children's hospital and saw firsthand what this disease does to little kids and their families. It traumatized me. I wanted to help — but my mind wasn't built for medicine. Mathematics and computers came easier to me, and I went into computer science instead.

But now we have AI. And we should be using it for more than vibe-coding the next get-rich-quick app.

This repository is an attempt to do something that matters. It's a cross-literature analysis of thousands of cancer research articles, combined with Monte Carlo biochemical simulations, proposing a specific and testable idea: that physical ROS-generating modalities (PDT, SDT) should be evaluated as spatially selective ferroptosis inducers for drug-tolerant persister cells. The full analysis, data, simulations, and article are here — open, free, and reproducible.

**I want nothing in return.** If someone takes this idea, validates it in a lab, spins it off, and it works — the world benefits. That's the point. Much like the polio vaccine, breakthroughs against diseases that destroy lives should be a human right, not a revenue stream.

**The mission is to crowdsource the minds of the global community — researchers, engineers, students, anyone — amplified by AI, to work on problems that actually matter.** Not another SaaS product. Not another chatbot wrapper. Real problems. Hard problems. The kind where the payoff isn't money — it's fewer empty chairs at the dinner table.

If you have expertise in oncology, biochemistry, ferroptosis, photodynamic therapy, sonodynamic therapy, immunology, computational biology, or just have ideas — open an issue, submit a PR, or fork and run with it. Everything here is MIT licensed. Take it. Use it. Make it better.

— Ezequiel Lares

---

## The Idea

Drug-tolerant persister cells are ferroptosis-sensitive (established, PMID:41481741). The field is searching for clinical ferroptosis inducers using *only pharmacologic agents* (RSL3, erastin). Meanwhile, PDT (355 ferroptosis articles) and SDT (121) trigger ferroptosis through ROS + produce immunogenic cell death — advantages systemic drugs lack.

**The proposal**: Physical ROS modalities should be evaluated as persister-targeting tools offering (1) spatial selectivity and (2) ICD for immunotherapy synergy. SDT extends this to deep tumors where PDT can't reach.

**What's novel**: The modality-class question — should the persister field look beyond drugs? — has not been systematically framed despite 355 PDT and 121 SDT ferroptosis papers existing independently.

**Key caveat**: PDT has 40 years of development without demonstrating robust ICD-immune synergy in randomized trials.

## Computational Simulation

A Rust-based Monte Carlo simulation of the ferroptosis cascade (16M cells across 16 conditions) validates the biochemical plausibility:

| Phenotype | Control | RSL3 | SDT/PDT |
|-----------|---------|------|---------|
| Glycolytic | 0.00% | 0.00% | 87.2% |
| OXPHOS | 0.04% | 1.1% | 99.9% |
| **Persister (FSP1↓)** | 1.2% | **42.5%** | **100.0%** |
| Persister+NRF2 | 0.00% | 0.05% | 99.5% |

Key features: autocatalytic LP propagation gated by GSH/GPX4, dynamic GPX4 regulation, FSP1 as GPX4-independent repair pathway. All baselines <2%. Sensitivity: 22/22 holds (100%). RSL3 kills persisters (42.5%) — matching published biology (PMID:41481741). NRF2 rescues from RSL3 but NOT from SDT.

## Article

**Author:** Ezequiel Lares

- `article/drafts/v1.md` — Markdown draft (~12,000 words, 114 references)
- `article/drafts/v1.tex` — LaTeX with proper tables and embedded figures
- `article/drafts/v1.pdf` — Compiled PDF manuscript
- `article/references/bibliography.bib` — BibTeX bibliography (114 entries)

10+ review rounds including adversarial peer review, falsification analysis, novelty assessment, and simulation validation.

## Figures

| Figure | Content |
|--------|---------|
| Fig 1 | Publication trends 2015-2025 |
| Fig 2 | 19×22 mechanism-cancer heatmap |
| Fig 3 | Ferroptosis engagement comparison (χ²=97.3, p<10⁻²³) |
| Fig 4 | Molecular pathway overlap across modalities |
| Fig 5 | Literature disconnect between communities |
| Fig 6 | SDT ferroptosis-ICD chain evidence depth |
| Fig 7 | Monte Carlo simulation results (1M cells/condition) |

## Repository Structure

```
article/drafts/          # Manuscript (md, tex, pdf)
article/figures/         # 7 publication-quality figures (pdf + png)
article/references/      # BibTeX bibliography
corpus/by-pmid/          # 4,830 full-text articles with YAML frontmatter
corpus/abstracts/by-pmid/ # 5,584 abstract-only records
tags/                    # Pre-computed indexes
scripts/                 # Python pipeline (fetch, enrich, tag, analyze, figures)
simulations/             # Rust Monte Carlo ferroptosis simulation suite
analysis/                # Hypothesis documents and data analysis
books/                   # Reference textbooks (LFS)
plans/                   # Research plans
```

## Reproduction

```bash
# Corpus
pip install -r requirements.txt && cp .env.example .env
cd scripts/
python fetch_articles.py --query-file queries.txt --max 500
python fetch_semantic_scholar.py --mode search --max 200
python enrich_metadata.py && python tag_articles.py
python build_index.py && python analyze_corpus.py
python generate_figures.py

# Simulation
cd ../simulations && cargo build --release
cargo run --release -p sim-original        # Original single-cell model
cargo run --release -p sim-spatial         # Spatial tumor with energy physics
cargo run --release -p sim-window          # Vulnerability window dynamics
cargo run --release -p sim-icd             # ICD-immune cascade comparison
cargo run --release -p sim-combo           # Combination therapy optimizer

# PDF
cd ../article/drafts && pdflatex v1.tex && bibtex v1 && pdflatex v1.tex && pdflatex v1.tex
```

## Contributing

Open an issue. Submit a PR. Fork it and run. No permission needed — that's the point.

If you're a researcher with lab access and want to test the core hypothesis (PDT vs SDT vs RSL3 on persister cells, measuring cell death + ICD markers), please reach out. That single experiment is the decisive test.

## License

MIT License. See [LICENSE](LICENSE) for details.
