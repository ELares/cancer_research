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

- A local collection of cancer research papers (full-text and abstracts)
- Python scripts to fetch, tag, index, and analyze the corpus
- Generated analyses, gap notes, and a draft manuscript
- Rust simulations exploring biochemical dynamics

Everything is organised so you can re-run the pipeline, challenge the conclusions, or extend the work in directions we haven't thought of yet.

## Explore the work

The repository is structured to make it easy to dig in:

| Directory | What you'll find |
|-----------|-----------------|
| `analysis/` | Audits, gap notes, and interpretative documents |
| `article/drafts/` | Current manuscript drafts |
| `scripts/` | Python pipeline for fetching, tagging, and analysis |
| `simulations/` | Rust simulation code |
| `corpus/` | Raw text data (by PubMed ID) |
| `tags/` | Precomputed tag indexes |

Start with the files in `analysis/` if you want to see what we've concluded so far—and where we're still uncertain.

## Get it running

If you want to reproduce or modify the analysis:

```bash
pip install -r requirements.txt
cp .env.example .env

python scripts/tag_articles.py
python scripts/build_index.py
python scripts/analyze_corpus.py
python scripts/generate_figures.py
```

For the simulations:

```bash
cd simulations
cargo build --release
cargo run --release -p sim-original   # try other sim-* packages
```

## Contribute

This project is most useful when it's questioned, expanded, and corrected. You don't need to be a cancer researcher—curiosity and a willingness to look at the evidence are enough.

- Open an issue with a question, a counter-example, or a missing paper
- Submit a pull request that improves the code, the corpus, or the manuscript
- Fork the repo and go in a completely new direction—MIT license means you're free to do that

We're not trying to steer everyone toward one answer. The goal is to build a shared space where good ideas can emerge.

## License

MIT License. See [LICENSE](LICENSE).
