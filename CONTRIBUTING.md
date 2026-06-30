# Contributing

Thank you for your interest in contributing to this cancer research project. All contributions — code, analysis, manuscript corrections, corpus additions, and simulation extensions — are welcome.

## Prerequisites

- **Python 3.12+** with [uv](https://docs.astral.sh/uv/) (recommended) or pip
- **Rust 1.56+** (install via [rustup.rs](https://rustup.rs/); the repo pins the tested version in `simulations/rust-toolchain.toml`)
- **LaTeX** (pdflatex + bibtex) for manuscript compilation
- **Git LFS** for the `books/` directory
- **Optional:** Nix, if you want the reproducible outer dev shell in `flake.nix`

## Setup

```bash
git clone https://github.com/ELares/cancer_research.git
cd cancer_research
uv lock
uv sync --frozen
cd simulations && cargo build --release   # build all simulation binaries
cd ..
```

With Nix:

```bash
git clone https://github.com/ELares/cancer_research.git
cd cancer_research
nix develop
uv lock
uv sync --frozen
cd simulations && cargo build --release
cd ..
```

## Running Tests

```bash
# Python pipeline, news, figure traceability, invariant, and integration tests (351 tests)
uv run pytest tests/ -q

# Rust simulation tests (full workspace unit + integration suite)
cd simulations && cargo test --workspace
```

All tests must pass before submitting a PR.

## Code Style

- **Python**: Follow existing patterns. Use type hints, docstrings, and `tqdm` for progress bars. Scripts in `scripts/` import from `config.py` and `article_io.py`.
- **Rust**: Follow existing patterns. Use `cargo fmt` and `cargo clippy` before committing.
- **Manuscript** (`article/drafts/v1.md`): Follow the conventions in `article/AUTHORING.md` — heading levels, citation formats, narrative inflation guardrails.

## Opening an Issue

Issue templates live in `.github/ISSUE_TEMPLATE/`: bug report, corpus or literature contribution, simulation extension proposal, and manuscript correction. They encode the project's conventions (for example, that new simulation layers stay off by default and byte-identical), so following the matching template makes a proposal easier to act on. Blank issues are still allowed for anything that does not fit a template.

## Submitting Changes

1. **Fork and branch** from `main`
2. **Make focused changes** — one issue per PR where possible
3. **Run tests** (both Python and Rust)
4. **Describe your changes** in the PR body. The pull-request template (`.github/PULL_REQUEST_TEMPLATE.md`) has a checklist covering tests, the off-by-default byte-identical rule for simulation changes, and verifiable citations for claims
5. **Link to the relevant issue** if one exists

## Simulation layer policy (calibrate-or-cut, #501)

The simulation suite grew broad fast: most layers are off-by-default,
uncalibrated, and direction-only. To stop breadth from outrunning depth, a new
off-by-default biochem or tumor-microenvironment axis is **not** mergeable on the
"off-by-default byte-identical" rule alone. In the same PR it must also carry:

1. **A named calibration target.** A specific assay, dataset, or measured
   quantity that would fit or falsify the layer. "The direction is the result" is
   acceptable for the magnitude, but the PR must say what data would calibrate it.
2. **A row in `simulations/calibration/CALIBRATION_STATUS.md`** stating the
   layer's tier, its used-in-any-reported-number status (expected: N for an
   off-by-default layer), and that calibration target.

Off-by-default byte-identical remains necessary (it keeps the production matrix
and the manuscript numbers safe), but it is no longer sufficient. Prefer
calibrating or improving an existing layer over adding a new one; a model that
does a few things with data behind them is more credible than one that does many
with none. A layer that has no calibration path and is used in no claim is a
retire-candidate, not a default keep. See the "Calibrate-or-cut policy and
accounting" section of `CALIBRATION_STATUS.md` for the standing accounting.

## Dependency Management

When adding a Python dependency:
1. Add it to `pyproject.toml`
2. Refresh the lockfile:
   ```bash
   uv lock
   ```
3. Regenerate the compatibility exports:
   ```bash
   make export-python
   ```
4. Commit `pyproject.toml`, `uv.lock`, and any regenerated compatibility files

The top-level Python source of truth is `pyproject.toml` + `uv.lock`. The
committed `requirements-lock.txt` and `requirements-dashboard.txt` are
transitional compatibility exports for non-uv consumers and automation that has
not migrated yet. `requirements.txt` remains as a short-term bridge file for
older setups, but new dependency edits should start in `pyproject.toml`.

## Maintainer Notes

- Sync the default environment with `uv sync --frozen`
- Generate or refresh the lockfile with `uv lock`
- Add dashboard dependencies with `uv sync --frozen --group dashboard`
- Regenerate compatibility exports with `make export-python`
- Verify the main workflows with `make test` and `make reproduce`
- Verify the Docker path still matches the lockfile with `docker build -t cancer-research .`

The nested [`simulations/ferroptosis-python/pyproject.toml`](simulations/ferroptosis-python/pyproject.toml)
is intentionally separate. It remains the packaging source of truth for the
distributable `ferroptosis-core` wheel and should not be merged into the
top-level app environment.

## Where to Start

**New here?** Read the one-page [plain-language explainer](docs/EXPLAINER.md) first,
then pick up a [**good first issue**](https://github.com/ELares/cancer_research/issues?q=is%3Aissue+is%3Aopen+label%3A%22good+first+issue%22)
— these are small, self-contained tasks (deploy the dashboard demo, add a unit test,
extend the landmark registry, review a flagged duplicate) that don't need deep domain
expertise. The categories below are the broader on-ramps:

- **Corpus contributions**: Submit a PR adding a missing landmark paper (see `analysis/landmark-corpus-gaps.md` for known gaps)
- **Simulation extensions**: See Chapter 8.4 in the manuscript for structural limitations that could become new simulation features
- **Manuscript corrections**: If you find a factual error, unclear wording, or missing caveat, open an issue or submit a PR
- **Pipeline improvements**: The scripts in `scripts/` follow a fetch → enrich → tag → index pattern. Extensions follow the same pattern.

## Project Philosophy

See `CLAUDE.md` for guiding principles: let the evidence lead, stay open, be honest about what we don't know, make it reproducible, keep it human.

## Code of Conduct

Be respectful, constructive, and focused on the mission: contributing to cancer research openly and honestly. See [GitHub's Community Guidelines](https://docs.github.com/en/site-policy/github-terms/github-community-guidelines).
