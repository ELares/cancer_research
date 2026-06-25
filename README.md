# Cancer Research Synthesis and Analysis

## Why This Exists

Too often has cancer taken from us the people we love.

There was a point in my life where I volunteered at a children's hospital and saw firsthand what this disease does to little kids and their families. It traumatized me. I wanted to help — but my mind wasn't built for medicine. Mathematics and computers came easier to me, and I went into computer science instead.

But now we have AI. And we should be using it for more than vibe-coding the next get-rich-quick app.

This repository is an attempt to do something that matters. It's a cross-literature analysis of thousands of cancer research articles, combined with Monte Carlo biochemical simulations, all open and reproducible. The goal is to let the evidence guide us—not to lock into one hypothesis from the start. The full analysis, data, simulations, and ongoing drafts are here.

**I want nothing in return.** If someone takes an idea from this work, validates it in a lab, spins it off, and it works — the world benefits. That's the point. Much like the polio vaccine, breakthroughs against diseases that destroy lives should be a human right, not a revenue stream.

**The mission is to crowdsource the minds of the global community — researchers, engineers, students, anyone — amplified by AI, to work on problems that actually matter.** Not another SaaS product. Not another chatbot wrapper. Real problems. Hard problems. The kind where the payoff isn't money — it's fewer empty chairs at the dinner table.

If you have expertise in oncology, biochemistry, ferroptosis, immunology, computational biology, or just have ideas — open an issue, submit a PR, or fork and run with it. Everything here is MIT licensed. Take it. Use it. Make it better.

— Ezequiel Lares

## What's here

- **4,830 full-text cancer research articles** across 19 mechanisms, 22 cancer types, 803 journals (2001-2026); the corpus skews toward immunotherapy (the single most-studied mechanism, ~1,685 articles), with physical and pharmacologic ferroptosis approaches a smaller, more preclinical slice
- **Python pipeline** for corpus fetching, tagging (7 tag layers), indexing, analysis, and figure generation
- **11 Rust simulation binaries**, a mechanistic claim-testing engine for cancer therapies: single-cell and spatial Monte Carlo, drug penetration across tissue types, drug combinations, tumor microenvironment (oxygen gradients, spatial immune zones, DAMP-mediated T-cell activation, stromal shielding, vasculature, clonal heterogeneity), vulnerability windows, ICD immune cascades, and tumor PK. Worked implementations include ferroptosis/RSL3 biochemistry and PDT/SDT depth physics (2D row-based and 3D radial-depth dispatchers; sim-tme-3d is the 3D-spheroid capstone consuming the full TME library stack) plus photosensitizer PK (drug-light-interval scaling, saturating distribution phase, relative singlet-O₂ yield)
- **ferroptosis-core library** (MIT, with Python bindings) — embeddable ferroptosis biochemistry engine; module list and current unit-test count in [`simulations/ferroptosis-core/README.md`](simulations/ferroptosis-core/README.md)
- **Calibration infrastructure** linking simulation parameters to published experimental data
- **[Model card](MODEL_CARD.md)** with the simulation suite's intended use, out-of-scope cases, assumptions/scope checklist, and per-layer calibration/validation status (the honest "broad but mostly uncalibrated" accounting, consolidated from [`CALIBRATION_STATUS.md`](simulations/calibration/CALIBRATION_STATUS.md))
- **Book-format manuscript (~115 pp)** with 11 chapters, 3 appendices, and 24 figures (~39,400 words), cross-referenced against all analysis outputs

Everything is organised so you can re-run the pipeline, challenge the conclusions, or extend the work in directions we haven't thought of yet.

## What we found

This work is first a **consolidation of the cancer-therapy literature**: mapping where research is concentrated, where apparent gaps are artifacts of search design rather than biology, and which mechanistic ideas can be compared on shared axes (evidence depth, resistant-state relevance, delivery constraints, tissue access). Immunotherapy dominates the corpus, and the analysis is deliberately honest about coverage limits (the evidence tagger has 96% binary evidence-presence precision but only 55% recall, so absence claims are provisional; an off-by-default MeSH-descriptor fallback lifts that recall to ~68% at ~95% precision but is not yet applied to the production corpus).

On top of that landscape, the simulations act as a **claim-testing engine**: we take specific mechanistic claims and try to validate or disprove them with reproducible, fact-grounded models. Three results that, if validated experimentally, would have translational implications:

1. **Combination synergy (ferroptosis case study).** Dual inhibition of GPX4 and FSP1 produces 1.99× Bliss synergy, because depleting both parallel repair pathways drops antioxidant defense below the autocatalytic lipid-peroxidation threshold. A general lesson about combining parallel-pathway blocks, tested in the RSL3 system.

2. **Microenvironment barriers affect drug-based and physical approaches differently.** Under simulated hypoxia, stromal shielding, and acidic pH, pharmacologic ferroptosis (RSL3) kill collapses (hypoxia 3.7% to 0.1%; stromal 3.0% to 1.5%; pH 163 to 77) while light- and ultrasound-delivered ROS (PDT/SDT) are less affected. This is one worked comparison of how mechanistically distinct modalities meet different barrier landscapes; it is directional, not a verdict. The hypoxia leg is the least certain (SDT's own oxygen-dependence is contested), and the immune-coupling amplification (a model-predicted 104× more immune kills, medium confidence) shrinks to roughly 4:1 in 3D.

3. **In-vitro-to-in-vivo penetration gap (applies to any systemic drug).** Tissue-specific delivery drops a RSL3-like drug from 40% (2D culture) to 12.1% (well-vascularized) to 2.6% (poorly-vascularized) to 1.8% (CNS/BBB), even at the blood vessel wall.

These are computational predictions with documented assumptions and caveats, not clinical claims. All parameters are documented with literature sources and confidence ratings. See the [manuscript](article/drafts/v1.md) for full context.

## Explore the work

| Directory | What you'll find |
|-----------|-----------------|
| `analysis/` | 15+ analysis outputs: evidence tiers, tissue-of-origin, diagnostic-therapy matching, combination audits, gap analysis |
| `article/drafts/` | Manuscript (v1.md + v1.tex) with 24 figures |
| `scripts/` | Python pipeline: tagging, indexing, analysis, figure generation, LaTeX generation, news authentication pipeline |
| `simulations/` | [11 Rust binaries](simulations/README.md) (each with its own README) + [ferroptosis-core library](simulations/ferroptosis-core/) + [Python bindings](simulations/ferroptosis-python/) + [calibration](simulations/calibration/) |
| `corpus/` | Full-text articles by PubMed ID + INDEX.jsonl |
| `tags/` | Precomputed tag indexes (mechanism, cancer type, tissue, evidence level, diagnostic-therapy) |
| `news/` | News source scaffolding: fetched articles, extracted claims, verification results, credibility scores |
| `tests/` | 334 Python tests (pipeline smoke + figure traceability + manuscript-inventory drift guard + depth-kill physics-constant guard + flagship-figure data guard + quantitative-figure drift guards (Figs 21/22/23) + invariant/integration + calibrate-extractor + MeSH evidence-fallback + gold-set precision-floor regression (#346) + Bliss/sim-tme/penetration prior-predictive intervals + ABC posterior (#332) + non-circular mechanism-recall (#412) + CTRPv2 calibration target + in-vitro kill-switch fit (#330) + System Xc-/erastin fit (#502) + joint multi-inducer posterior (#500) + spheroid structure validation (#333) + embedding evidence leg (#411) + RD-vs-BioFVM cross-check (#408) + dashboard data layer (#354) + tumor-PK measured-data anchor (#334) + Krogh penetration validation (#335) + spheroid size-aware zone thresholds (#333) + spheroid kill-vs-size direction (#333) + ferroptosis-python bindings) |

Start with the files in `analysis/` if you want to see what we've concluded so far—and where we're still uncertain.

## Get it running

See [CONTRIBUTING.md](CONTRIBUTING.md) for full setup instructions, or the quick version:

```bash
pip install -r requirements.txt          # or requirements-lock.txt for exact versions
cp .env.example .env

python scripts/tag_articles.py
python scripts/build_index.py
python scripts/analyze_corpus.py
python scripts/generate_figures.py
```

For the simulations (see [simulations/README.md](simulations/README.md) for all 11 binaries):

```bash
cd simulations
cargo build --release
cargo test --workspace                  # ferroptosis-core unit tests + per-binary integration tests
cargo run --release -p sim-original     # Monte Carlo ferroptosis baseline
cargo run --release -p sim-spatial      # 2D tumor with PDT/SDT depth physics
cargo run --release -p sim-tissue-pk    # drug penetration across tissue types
cargo run --release -p sim-combo-mech   # pairwise drug combination synergy
cargo run --release -p sim-tme          # tumor microenvironment (O2 gradients)
```

For the Python bindings:

```bash
cd simulations
pip install maturin
maturin develop -m ferroptosis-python/Cargo.toml --release
python -c "import ferroptosis_core as fc; print(fc.sim_batch('Persister', 'RSL3', n=1000, seed=42))"
```

For the interactive dashboard (corpus exploration + a single-cell parameter sweep):

```bash
pip install -r requirements-dashboard.txt   # optional UI deps (streamlit, pandas); not in the pinned core
streamlit run scripts/dashboard.py
```

The Corpus tab (filters, mechanism/cancer/evidence views, the mechanism x cancer
matrix) needs only the committed `corpus/INDEX.jsonl`. The Simulation-sweep tab runs
a live `ferroptosis_core.sim_batch` sweep when the bindings above are built, and
otherwise degrades to the committed prior-predictive intervals. Self-hosting: behind
auth, `streamlit run scripts/dashboard.py --server.address 0.0.0.0 --server.port 8501`.

## Philosophy

**The work is more important than the paper.** We don't optimize for journal word limits or publication formats. If a finding needs context, we give context. If a decision needs explaining, we explain it. Every result in this repo includes the reasoning chain that produced it — what we assumed, what we measured, what we're uncertain about, and why we believe the finding signals value despite those uncertainties.

We'd rather publish a longer, clearer document that a graduate student can follow end-to-end than a compressed paper that only specialists can decode. Breakthroughs against diseases that destroy lives should be accessible to anyone willing to read carefully.

## Contribute

This project is most useful when it's questioned, expanded, and corrected. You don't need to be a cancer researcher—curiosity and a willingness to look at the evidence are enough.

- Read [CONTRIBUTING.md](CONTRIBUTING.md) for setup, testing, and PR guidelines
- Open an issue with a question, a counter-example, or a missing paper. Issue templates (bug, corpus/literature contribution, simulation extension, manuscript correction) are in [.github/ISSUE_TEMPLATE/](.github/ISSUE_TEMPLATE/)
- Submit a pull request that improves the code, the corpus, or the manuscript
- See [CONTRIBUTORS.md](CONTRIBUTORS.md) for how contributions are recognized
- Fork the repo and go in a completely new direction—MIT license means you're free to do that

The model's falsifiable predictions and the experiments that would confirm or refute them are registered in [PREREGISTRATION.md](PREREGISTRATION.md), so the predictions are locked in before the calibration work that tests them.

We're not trying to steer everyone toward one answer. The goal is to build a shared space where good ideas can emerge.

## Cite this work

If you use this work in your research, GitHub renders a "Cite this repository" button from [CITATION.cff](CITATION.cff). Content provenance and redistribution rights for bundled assets are documented in [PROVENANCE.yaml](PROVENANCE.yaml).

## License

MIT License (code only). See [LICENSE](LICENSE). Bundled data assets have their own licenses — see [PROVENANCE.yaml](PROVENANCE.yaml) for details.
