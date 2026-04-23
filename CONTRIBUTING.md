# Contributing

Thank you for your interest in contributing to this cancer research project. All contributions — code, analysis, manuscript corrections, corpus additions, and simulation extensions — are welcome.

## Prerequisites

- **Python 3.10+** with pip
- **Rust 1.56+** (install via [rustup.rs](https://rustup.rs/); the repo pins the tested version in `simulations/rust-toolchain.toml`)
- **LaTeX** (pdflatex + bibtex) for manuscript compilation
- **Git LFS** for the `books/` directory

## Setup

```bash
git clone https://github.com/ELares/cancer_research.git
cd cancer_research
pip install -r requirements.txt          # or requirements-lock.txt for exact versions
cd simulations && cargo build --release   # build all simulation binaries
cd ..
```

## Running Tests

```bash
# Python pipeline + news tests (50 tests)
python3 -m pytest tests/ -q

# Rust simulation tests (31+ tests)
cd simulations && cargo test --workspace
```

All tests must pass before submitting a PR.

## Code Style

- **Python**: Follow existing patterns. Use type hints, docstrings, and `tqdm` for progress bars. Scripts in `scripts/` import from `config.py` and `article_io.py`.
- **Rust**: Follow existing patterns. Use `cargo fmt` and `cargo clippy` before committing.
- **Manuscript** (`article/drafts/v1.md`): Follow the conventions in `article/AUTHORING.md` — heading levels, citation formats, narrative inflation guardrails.

## Submitting Changes

1. **Fork and branch** from `main`
2. **Make focused changes** — one issue per PR where possible
3. **Run tests** (both Python and Rust)
4. **Describe your changes** in the PR body — what changed, why, and what was tested
5. **Link to the relevant issue** if one exists

## Dependency Management

When adding a Python dependency:
1. Add it to `requirements.txt` with a minimum version (`>=X.Y.Z`)
2. Regenerate the lockfile:
   ```bash
   python3 -m venv /tmp/lock-env
   /tmp/lock-env/bin/pip install -r requirements.txt
   /tmp/lock-env/bin/pip freeze > requirements-lock.txt
   rm -rf /tmp/lock-env
   ```
3. Commit both `requirements.txt` and `requirements-lock.txt`

## Where to Start

- **Corpus contributions**: Submit a PR adding a missing landmark paper (see `analysis/landmark-corpus-gaps.md` for known gaps)
- **Simulation extensions**: See Chapter 8.4 in the manuscript for structural limitations that could become new simulation features
- **Manuscript corrections**: If you find a factual error, unclear wording, or missing caveat, open an issue or submit a PR
- **Pipeline improvements**: The scripts in `scripts/` follow a fetch → enrich → tag → index pattern. Extensions follow the same pattern.

## Project Philosophy

See `CLAUDE.md` for guiding principles: let the evidence lead, stay open, be honest about what we don't know, make it reproducible, keep it human.

## Code of Conduct

Be respectful, constructive, and focused on the mission: contributing to cancer research openly and honestly. See [GitHub's Community Guidelines](https://docs.github.com/en/site-policy/github-terms/github-community-guidelines).
