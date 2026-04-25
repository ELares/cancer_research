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

- **4,830 full-text cancer research articles** across 19 mechanisms, 22 cancer types, 803 journals (2001-2026)
- **Python pipeline** for corpus fetching, tagging (7 tag layers), indexing, analysis, and figure generation
- **10 Rust simulation binaries** modeling ferroptosis biochemistry: single-cell Monte Carlo, spatial tumors, drug penetration, drug combinations, tumor microenvironment (oxygen gradients, spatial immune zones, DAMP-mediated T cell activation), vulnerability windows, ICD immune cascades, tumor PK
- **ferroptosis-core library** (MIT, with Python bindings) — embeddable ferroptosis biochemistry engine with 10 modules and 31 unit tests
- **Calibration infrastructure** linking simulation parameters to published experimental data
- **112-page book-format manuscript** with 11 chapters, 3 appendices, and 20 figures (~34,600 words), cross-referenced against all analysis outputs

Everything is organised so you can re-run the pipeline, challenge the conclusions, or extend the work in directions we haven't thought of yet.

## What we found

Three simulation findings that, if validated experimentally, would have translational implications:

1. **RSL3 + FSP1 inhibitor produces 1.99× Bliss synergy** through dual-pathway depletion — depleting both GPX4 and FSP1 repair pathways simultaneously drops antioxidant defense below the autocatalytic lipid peroxidation threshold.

2. **The tumor microenvironment selectively favors physical over pharmacologic ferroptosis inducers through four mechanisms.** Hypoxia collapses RSL3 kill from 3.7% to 0.1%; SDT maintains 87.8%. Stromal shielding halves RSL3's peripheral kill (3.0% to 1.5%); SDT barely affected. Acidic pH halves RSL3 ferroptosis kills (163→77, ion trapping dominates over iron release); SDT gains slightly (+0.8%). SDT's dense kill field generates a model-predicted 104× more immune kills than RSL3 (medium confidence). Three resistance mechanisms (hypoxia, stromal, pH) and one amplification effect (immune coupling) all favor physical modalities.

3. **Tissue-specific drug penetration creates a substantial in-vitro-to-in-vivo gap.** RSL3-like drug kill drops from 40% (2D culture) to 12.1% (well-vascularized) to 2.6% (poorly-vascularized) to 1.8% (CNS/BBB) — even at the blood vessel wall.

These are computational predictions with documented assumptions and caveats, not clinical claims. All parameters are documented with literature sources and confidence ratings. See the [manuscript](article/drafts/v1.md) for full context.

## Explore the work

| Directory | What you'll find |
|-----------|-----------------|
| `analysis/` | 15+ analysis outputs: evidence tiers, tissue-of-origin, diagnostic-therapy matching, combination audits, gap analysis |
| `article/drafts/` | Manuscript (v1.md + v1.tex) with 20 figures |
| `scripts/` | Python pipeline: tagging, indexing, analysis, figure generation, LaTeX generation, news authentication pipeline |
| `simulations/` | [10 Rust binaries](simulations/README.md) + ferroptosis-core library + [Python bindings](simulations/ferroptosis-python/) + calibration infrastructure |
| `corpus/` | Full-text articles by PubMed ID + INDEX.jsonl |
| `tags/` | Precomputed tag indexes (mechanism, cancer type, tissue, evidence level, diagnostic-therapy) |
| `news/` | News source scaffolding: fetched articles, extracted claims, verification results, credibility scores |
| `tests/` | 79 Python tests (50 pipeline smoke + 10 figure traceability + 19 invariant/integration) |

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

For the simulations (see [simulations/README.md](simulations/README.md) for all 10 binaries):

```bash
cd simulations
cargo build --release
cargo test --workspace                  # 31 unit tests
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

## Philosophy

**The work is more important than the paper.** We don't optimize for journal word limits or publication formats. If a finding needs context, we give context. If a decision needs explaining, we explain it. Every result in this repo includes the reasoning chain that produced it — what we assumed, what we measured, what we're uncertain about, and why we believe the finding signals value despite those uncertainties.

We'd rather publish a longer, clearer document that a graduate student can follow end-to-end than a compressed paper that only specialists can decode. Breakthroughs against diseases that destroy lives should be accessible to anyone willing to read carefully.

## Contribute

This project is most useful when it's questioned, expanded, and corrected. You don't need to be a cancer researcher—curiosity and a willingness to look at the evidence are enough.

- Read [CONTRIBUTING.md](CONTRIBUTING.md) for setup, testing, and PR guidelines
- Open an issue with a question, a counter-example, or a missing paper
- Submit a pull request that improves the code, the corpus, or the manuscript
- Fork the repo and go in a completely new direction—MIT license means you're free to do that

We're not trying to steer everyone toward one answer. The goal is to build a shared space where good ideas can emerge.

## Cite this work

If you use this work in your research, GitHub renders a "Cite this repository" button from [CITATION.cff](CITATION.cff). Content provenance and redistribution rights for bundled assets are documented in [PROVENANCE.yaml](PROVENANCE.yaml).

## License

MIT License (code only). See [LICENSE](LICENSE). Bundled data assets have their own licenses — see [PROVENANCE.yaml](PROVENANCE.yaml) for details.
