# Cancer Research Synthesis and Analysis

## Why This Exists

Too often has cancer taken from us the people we love.

There was a point in my life where I volunteered at a children's hospital and saw firsthand what this disease does to little kids and their families. It traumatized me. I wanted to help — but my mind wasn't built for medicine. Mathematics and computers came easier to me, and I went into computer science instead.

But now we have AI. And we should be using it for more than vibe-coding the next get-rich-quick app.

This repository is an attempt to do something that matters. It's a cross-literature analysis of thousands of cancer research articles, combined with Monte Carlo biochemical simulations, all open and reproducible. The goal is to let the evidence guide us—not to lock into one hypothesis from the start. The full analysis, data, simulations, and ongoing drafts are here.

**I want nothing in return.** If someone takes an idea from this work, validates it in a lab, spins it off, and it works — the world benefits. That's the point. Much like the polio vaccine, breakthroughs against diseases that destroy lives should be a human right, not a revenue stream.

**The mission is to crowdsource the minds of the global community — researchers, engineers, students, anyone — amplified by AI, to work on problems that actually matter.** Not another SaaS product. Not another chatbot wrapper. Real problems. Hard problems. The kind where the payoff isn't money — it's fewer empty chairs at the dinner table.

If you have expertise in oncology, biochemistry, radiation physics, immunology and cell therapy, drug delivery, virology, computational biology — or just have ideas — open an issue, submit a PR, or fork and run with it. Everything here is MIT licensed. Take it. Use it. Make it better.

— Ezequiel Lares

## What's here

This project reads the cancer literature at census scale. It holds **5,187,265
cancer articles** — the 4,403,994 that MeSH indexes under the neoplasms tree,
plus 783,271 more recovered by text-matching because MeSH has not indexed them
yet — with **1,116,481 open-access full texts**, **10,700,928 typed entity
relations** over 2,129,080 articles, and 289 million sentences mined for
co-mention. A daily-update stream adds the newest literature on top: the latest
window carried 188,850 distinct articles, 65,966 of them new to the census.

That is the whole indexed cancer literature, not a sample of it. It exists
because most literature analyses — including this project's own, at first — are
computed over a few thousand hand-retrieved papers, and there is no way to tell
from inside such a corpus whether an apparent gap in the science is real or an
artifact of which papers were retrieved. At census scale you can check.

**Three surfaces, and it matters which one a number comes from:**

- **The census** (above): every cancer article PubMed has indexed. Built from
  the PubMed annual baseline because E-utilities cannot page past 10,000 hits,
  defined as MeSH tree C04 plus nine adjacent descriptors so that foundational
  mechanism papers outside the cancer tree are not lost — both founding FSP1
  papers are outside C04, as is 62% of all ferroptosis work.
- **An open-access full-text layer** of **1,116,481 cancer articles** (25.4% of
  the indexed census), drawn from PMC's 28 bulk packages — 737,929
  redistributable `oa_comm`, 378,552 `oa_noncomm`. Text is what an evidence
  label can be read from; the other three quarters of the census are measured
  from MeSH descriptors and NLM publication types instead.
- **A living review** that re-runs monthly and reports a dated delta.

Every prevalence figure here is a share of the census, measured on
expert-assigned MeSH descriptors rather than on this project's own retrieval.
That distinction is not cosmetic: it is the difference between describing the
literature and describing a search.

- **Python pipeline** for corpus fetching, tagging (7 tag layers), indexing, analysis, and figure generation
- **13 Rust simulation binaries**, a mechanistic claim-testing engine covering
  **ten treatment arms plus an untreated control**: ferroptosis induction
  (RSL3), photodynamic and sonodynamic therapy, ionizing radiation, cytotoxic
  chemotherapy, checkpoint blockade, adoptive cell therapy, oncolytic virus,
  thermal and electrical ablation, and antibody-drug conjugates. Around them sit single-cell and
  spatial Monte Carlo, drug penetration across tissue types, drug combinations,
  the tumour microenvironment (oxygen gradients, spatial immune zones,
  DAMP-mediated T-cell activation, stromal shielding, vasculature, clonal
  heterogeneity), vulnerability windows, ICD immune cascades and tumour PK
  — with `sim-tme-3d` as the 3D-spheroid capstone and `sim-modality-panel`
  running every applicable arm against the identical tumour from the identical
  seed. **Read the depth honestly:** the arms a modality owns outright are
  9 modules and 1562 lines against the ferroptosis engine's 26 modules and
  3,805, so the newer arms are roughly 2.4 to 3.0 times smaller — measured in
  [`analysis/modality-module-depth.md`](analysis/modality-module-depth.md),
  which reports it against itself rather than leaving you to count.
- **ferroptosis-core library** (MIT, with Python bindings) — the ferroptosis
  biochemistry engine plus the modality layers built beside it (radiation,
  ablation, oncolytic spread, ADC bystander killing, adoptive-cell barriers);
  the name is historical and the crate is now broader than it. Module list and
  current unit-test count in
  [`simulations/ferroptosis-core/README.md`](simulations/ferroptosis-core/README.md)
- **Calibration infrastructure** linking simulation parameters to published experimental data
- **[Model card](MODEL_CARD.md)** with the simulation suite's intended use, out-of-scope cases, assumptions/scope checklist, and per-layer calibration/validation status (the honest "broad but mostly uncalibrated" accounting, consolidated from [`CALIBRATION_STATUS.md`](simulations/calibration/CALIBRATION_STATUS.md))
- **Book-format manuscript (~270 pp at 6x9 trim)** with 12 chapters, 3 appendices, and 40 figures (~73,900 words), cross-referenced against all analysis outputs and indexed in [`FIGURES.yaml`](FIGURES.yaml)

## What the work is actually about

The corpus is broad. **The work is not, and that is worth saying plainly rather
than leaving you to count files.**

| | ferroptosis / physical-ROS | other therapy | method & tooling |
|---|--:|--:|--:|
| committed analyses | 21 | **1** | 123 |
| preregistered predictions | **10 of 22** | 0 | — |
| engine modules mentioning it anywhere | **36 of 40** | — | — |
| engine modules mentioning it in code | **21 of 40** | — | — |
| engine modules mentioning it in PRODUCTION code | **15 of 40** | — | — |

The module rows read `33 of 33` until 2026-08-17 — the same number on both
sides of "of", produced by counting `.rs` files without opening one. Measured,
four modules (`adoptive`, `reaction_diffusion`, `spheroid`,
`vasculature`) mention neither ferroptosis nor a physical-ROS modality anywhere
in their text — and the count has MOVED, in the direction that matters: it was
five, and `oncolytic` left the list by acquiring prose about the immune bind
rather than by anyone editing this sentence. And the `1` is
a filename marker rather than a subject measurement: the therapy bucket matches
on filenames only while the ferroptosis bucket also reads body text, so the
count moves on a rename. Applying the ferroptosis rule's body route to a
therapy vocabulary admits 50 -- but 10 of those sit in this table's own
ferroptosis column, so 50 is not a bound either. Neither rule measures subject. See
[`analysis/scope-audit.md`](analysis/scope-audit.md).

8 of 20 preregistered predictions, and 33 of 38 modules of the simulation engine, concern ferroptosis or the physical-ROS modalities (PDT and
SDT). That first count read "every falsifiable commitment this project makes" until P9 to P13 registered predictions for the modality arms, and it sat sixteen lines below a table this same change had already updated to 8 of 13 -- the file contradicted itself. One committed analysis is FILED as taking another therapy as its subject, and that filing is a filename match rather than a measurement -- see the scope audit for what the rules do and do not establish.

A narrow thesis on a broad corpus is how most good science works, and the census
above is genuinely broad. But a reader arriving at a front door that opens with
five million articles will reasonably infer the analysis is commensurate, and it
is not — so: **the census is a measuring instrument, and the thing being measured
is still mostly ferroptosis.** The largest bucket above is method work, which is
about the instrument rather than about any biology.

**That is less true than it was, and the table is the place to see by how
much.** A deliberate campaign widened the engine from 3 selectable treatment
arms to 9, took the mechanisms with no engine representation from 13 of 16 to
0 of 16, and registered the first 5 preregistered predictions that are not
about ferroptosis or physical ROS — which is why the predictions row reads
8 of 20 rather than the 8 of 8 it read for most of this project's life. None of
that makes the work broad yet: the modality arms are roughly 2 to 3 times
smaller than the ferroptosis engine by line count, every one of them is
recorded in `CALIBRATION_STATUS.md` as feeding no number in the manuscript's
QUANTITATIVE chapters — which is narrower than "no number this manuscript
reports", the phrasing used here until Chapter 6 began reporting these arms'
own figures and made it false. Of the 8 arms that have a
published calibration target at all, 4 reproduce it — 1 of those 4 is genuinely
PINNED by its observable and the other 3 merely round-trip — while 1 is
UNCONSTRAINED, 1 INADMISSIBLE, 1 constrains only a DIRECTION, and 1 is PARTLY
REFUTED: an independent study contradicts one of its three directional claims,
which is reported as prominently as the two that survive. The remaining 2 arms
have no target to fit at all. The gap is narrower and it is
now measured rather than argued.

Two consequences worth carrying into anything you read here:

- Mechanism shares are shares of the **census records carrying a
  discriminative MeSH descriptor** (165,700 across the measured mechanisms).
  Immunotherapy is 31,890 = **19.2%**. Read that against what a keyword
  retrieval reported for the same mechanism — 47.6% — and the gap is the point:
  a corpus assembled from queries about immunotherapy finds immunotherapy
  everywhere. Two cautions on the ranking: `epigenetic` leads on volume only
  because 75% of its records come from `DNA Methylation`, a descriptor carried
  by any paper that MEASURES the process, so its rank is a scope artifact and
  not a discovery; and the taxonomy reaches about 6.8% of the census at all —
  see [`analysis/atlas-taxonomy-reach.md`](analysis/atlas-taxonomy-reach.md).
- Areas with no lane at all are not absent because they were weighed and
  dismissed. **Radiotherapy was this section's worked example of that, and it
  no longer is:** it has a `Treatment::Radiation` arm, a linear-quadratic
  DNA-damage channel separate from the ferroptosis one, and a row in the
  coverage table. Of the 16 mechanisms this project's taxonomy can measure at
  census scale, **0 now have no engine representation at all**, against 13 when
  that count was first taken. What replaces the absence question is a harder
  one the coverage page now leads with — presence is not applicability. 8 of
  the 16 can be APPLIED as a treatment a run selects; the other 8 are MODIFIERS
  that only change how another treatment lands. That split read 2 and 14
  until the tier rule was checked: it decided by matching a mechanism's
  keyword against the SPELLING of a Rust enum variant, and `oncolytic` cannot
  match inside `OncolyticVirus` — so four arms this README lists above were
  reported as unselectable while `sim-modality-panel` was selecting them.
  And "representation" is a low bar in the other direction: for several
  mechanisms it is a function with no production caller at all. See
  [`analysis/modality-coverage.md`](analysis/modality-coverage.md).

Counts derived by [`scripts/scope_audit.py`](scripts/scope_audit.py); the
bucketing is listed in [`analysis/scope-audit.md`](analysis/scope-audit.md) so a
placement can be disputed.

Everything is organised so you can re-run the pipeline, challenge the conclusions, or extend the work in directions we haven't thought of yet.

## What we found

This work is first a **consolidation of the cancer-therapy literature**: mapping where research is concentrated, where apparent gaps are artifacts of search design rather than biology, and which mechanistic ideas can be compared on shared axes (evidence depth, resistant-state relevance, delivery constraints, tissue access). Immunotherapy dominates the corpus, and the analysis is deliberately honest about coverage limits (the evidence tagger has 96% binary evidence-presence precision but only 55% recall, so absence claims are provisional; an off-by-default MeSH-descriptor fallback lifts that recall to ~68% at ~95% precision but is not yet applied to the production corpus; a further off-by-default rebuild of the evidence tagger — which reads the Methods/Results sections of the stored full text instead of the abstract alone — cuts the exact-label error 2.6x on held-out records and lifts binary F1 to 0.97, and is likewise not yet applied, see [`analysis/evidence-tagger-v2.md`](analysis/evidence-tagger-v2.md)).

On top of that landscape, the simulations act as a **claim-testing engine**: we take specific mechanistic claims and try to validate or disprove them with reproducible, fact-grounded models. Six results that, if validated experimentally, would have translational implications — the first three from the ferroptosis and physical-ROS work this project started with, the last three from the multi-modality arms added since.

> **Read the numbers below as directions, not magnitudes.** [`analysis/identifiability-report.md`](analysis/identifiability-report.md) prices them: with eleven free rate constants, six non-identifiable from the kill rate, and **zero** of these outputs conditioned on data in the regime that produces them, **none is point-estimable**. Each is given with its interval below; where an interval spans the plausible range, the number is order-of-magnitude and the direction is the result.


1. **Combination synergy (ferroptosis case study).** Dual inhibition of GPX4 and FSP1 produces ~1.99× Bliss synergy — 95% prior-predictive ~1.0× to 5.2×, median ~1.35×, so the **direction** is robust and the magnitude is not point-estimable — because depleting both parallel repair pathways drops antioxidant defense below the autocatalytic lipid-peroxidation threshold. A general lesson about combining parallel-pathway blocks, tested in the RSL3 system.

2. **Microenvironment barriers affect drug-based and physical approaches differently.** Under simulated hypoxia, stromal shielding, and acidic pH, pharmacologic ferroptosis (RSL3) kill collapses (hypoxia 3.7% to 0.1%; stromal 3.0% to 1.5%; pH 163 to 77) while light- and ultrasound-delivered ROS (PDT/SDT) are less affected. This is one worked comparison of how mechanistically distinct modalities meet different barrier landscapes; it is directional, not a verdict. The hypoxia leg is the least certain: the SDT hypoxic-zone advantage brackets 0% to 86.6%, collapsing to roughly 0% if SDT's ROS is fully O2-dependent — the regime the lead clinical agent SONALA-001 occupies. The immune-coupling amplification (a model-predicted 104× more immune kills, medium confidence) shrinks to roughly 4:1 in 3D.

3. **In-vitro-to-in-vivo penetration gap (applies to any systemic drug).** Tissue-specific delivery drops a RSL3-like drug from 40% (2D culture) to 12.1% (well-vascularized) to 2.6% (poorly-vascularized) to 1.8% (CNS/BBB), even at the blood vessel wall. The **ordering** is what is parameter-robust — it held in 300 of 300 draws — not these magnitudes, whose intervals span nearly the whole range at every tissue (well-vascularized median 23%, ~[0%, 93%]; CNS/BBB median 4%, ~[0%, 77%]).

4. **What transport costs, priced on one tumour.** Run every applicable arm against the identical tumour from the identical seed and an antibody-drug conjugate kills **1.8%** where sonodynamic therapy kills **87.2%**, both driven by the same exogenous-ROS constant through the same engine. **Read that as arithmetic more than as a discovery:** the two arms differ in the transport factor and in nothing else, so the ratio prices this model's transport layer rather than measuring an antibody in tissue. Two things keep it honest. The ADC arm does NOT run RSL3's pharmacology — the GPX4-inhibition branch fires only for `Treatment::RSL3` — and the panel said otherwise until it was checked against the code. And the free drug carrying the actual RSL3 payload kills **0.00%** in the same table, LESS than the antibody-delivered arm, because the glycolytic state is where a GPX4 inhibitor has almost nothing to act on. The panel is **not a ranking**: every arm prints its own calibration tier, and what it shows is structural — which arms need ferroptotic death and which do not, which are dose-responses and which are thresholds. See [`analysis/modality-panel.md`](analysis/modality-panel.md).

5. **The same CAR-T construct collapses 633-fold between two diseases, and no single step looks catastrophic.** Three barriers the literature names in one sentence — trafficking to, infiltration into, and activation within the tumour — multiply to 6% delivery; exhaustion removes most of what is left; the antigen ceiling contributes nothing here because the kill is already far below it. That decomposition is the point: a model fitting one efficacy scalar would reproduce the same endpoint and lose the reason, which is what tells you which step to attack. Every barrier value is an uncalibrated placeholder and the direction is the result.

6. **A registered prediction that contradicts this project's own prior belief.** The ADC module was built expecting a cleavable linker's advantage to *grow* as antigen is lost — the mechanism that would answer antigen escape — and shipped a guard asserting it. The guard was vacuous: the ratio is exactly constant, because bystander kill is proportional to the dying antigen-positive population, so both arms scale together. What the model actually says is that the bystander effect is **starved by the escape it answers** — the share of the antigen-negative pool it reaches falls from 77.1% to 2.6% as antigen is lost. That is registered as **P10** with the corrected sign, alongside four other falsifiable predictions for the new arms ([`PREREGISTRATION.md`](PREREGISTRATION.md)).

These are computational predictions with documented assumptions and caveats, not clinical claims. All parameters are documented with literature sources and confidence ratings. See the [manuscript](article/drafts/v1.md) for full context.

## Explore the work

| Directory | What you'll find |
|-----------|-----------------|
| `analysis/` | 144 analysis outputs. Frozen-corpus work (evidence tiers, tissue-of-origin, diagnostic-therapy matching, combination audits, gap analysis) plus 34 census-scale ones, including the drug-by-variant map, the co-treatment layer, retraction exposure, entity-ambiguity impact, and the manuscript-vs-census re-test |
| `article/drafts/` | Manuscript (v1.md + v1.tex) with 40 figures; [`FIGURES.yaml`](FIGURES.yaml) indexes 45 entries, the extra 4 being supplementary and 1 orphan |
| `scripts/` | Python pipeline: tagging, indexing, analysis, figure generation, LaTeX generation, news authentication. 38 of them are census-scale (33 `atlas_*.py` + 5 `comention_*.py`); `scripts/atlas_pipeline.sh` documents the dependency order, which is load-bearing and fails silently if run wrong |
| `simulations/` | [13 Rust binaries](simulations/README.md) (most with their own README; `sim-modality-panel`, `sim-scale` do not yet) + [ferroptosis-core library](simulations/ferroptosis-core/) + [Python bindings](simulations/ferroptosis-python/) + [calibration](simulations/calibration/) |
| `corpus/` | Frozen full text by PubMed ID + INDEX.jsonl; `corpus/atlas/` holds the census (bulk gitignored, committed artifacts in `analysis/`); `corpus/living/` documents the frozen-versus-living split (the monthly deltas themselves are uploaded as workflow artifacts, never committed) |
| `tags/` | Precomputed tag indexes (mechanism, cancer type, tissue, evidence level, diagnostic-therapy) |
| `news/` | News source scaffolding: fetched articles, extracted claims, verification results, credibility scores |
| `tests/` | 1925 Python tests (pipeline smoke + figure traceability + calibration-status ref guard + manuscript-inventory drift guard + depth-kill physics-constant guard + flagship-figure data guard + quantitative-figure drift guards (Figs 21/22/23) + invariant/integration + calibrate-extractor + MeSH evidence-fallback + gold-set precision-floor regression (#346) + Bliss/sim-tme/penetration prior-predictive intervals + ABC posterior (#332) + non-circular mechanism-recall (#412) + CTRPv2 calibration target + in-vitro kill-switch fit (#330) + System Xc-/erastin fit (#502) + joint multi-inducer posterior (#500) + spheroid structure validation (#333) + embedding evidence leg (#411) + RD-vs-BioFVM cross-check (#408) + dashboard data layer (#354) + tumor-PK measured-data anchor (#334) + Krogh penetration validation (#335) + spheroid size-aware zone thresholds (#333) + spheroid kill-vs-size direction (#333) + gene-symbol ambiguity/FSP1 sense disambiguation (#ATLAS-AMBIG) + rare-event Poisson intervals + tail-resolution classification + ferroptosis-python bindings) |

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

For the simulations (see [simulations/README.md](simulations/README.md) for all 13 binaries):

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

The Census tab reads the committed census aggregates under `analysis/` (~62 KB
of JSON): mechanism volumes, clinical-trial shares from NLM publication types,
anatomical concentration, convergence partners and growth. Record-level browsing
of the census is deliberately not offered — it is 5,187,265 records and
gitignored, so a panel that appeared to browse it would be browsing something
else. The Corpus tab (filters, mechanism/cancer/evidence views, the mechanism x
cancer matrix) needs only the committed `corpus/INDEX.jsonl` and browses the
earlier keyword-retrieved archive, retained as a method-comparison arm: holding
descriptor scope constant across it and the census is what separates a
labelling effect from a selection effect. No figure rests on it. The Simulation-sweep tab runs
a live `ferroptosis_core.sim_batch` sweep when the bindings above are built, and
otherwise degrades to the committed prior-predictive intervals. Self-hosting: behind
auth, `streamlit run scripts/dashboard.py --server.address 0.0.0.0 --server.port 8501`.

**Live demo: https://elares.github.io/cancer_research/** — the Corpus tab runs
entirely in your browser via [stlite](https://github.com/whitphx/stlite) (Streamlit
compiled to WebAssembly/Pyodide): it executes `scripts/dashboard.py` on the committed
census aggregates and `corpus/INDEX.jsonl` client-side, with **no server and no install** (first load ~30-60 s
while Pyodide + pandas + matplotlib download, then cached). The Simulation-sweep tab
shows a read-only notice pointing to the committed prior-predictive intervals (the
compiled `ferroptosis_core` extension is not available under Pyodide). The page is
`docs/index.html`, served by GitHub Pages;
alternatively the app deploys 1-click on Streamlit Community Cloud by pointing it at
`scripts/dashboard.py`.

**New here / not a specialist?** Start with the one-page plain-language explainer:
[`docs/EXPLAINER.md`](docs/EXPLAINER.md) — what ferroptosis is, what the three
headline results mean, and the "directional, not clinical" caveat, in everyday terms.

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
