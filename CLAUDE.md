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
- simulation work: ferroptosis biochemistry, drug penetration, calibration, photosensitizer PK (drug-light-interval scaling, saturating distribution phase, relative singlet-O₂ yield, FromStr-based clap CLI integration in sim-spatial), 3D spheroid scaffolding (TumorGrid3D #185, signed radial depth + 3D energy-physics dispatcher #186, 3D radial O₂ gradient + zone-census #187, 3D radial pH gradient + iron/ion-trap helpers #190, 3D CAF-shielded boundary detection #189, 3D spatial immune coupling + DAMP diffusion #188) for the #185–#197 spheroid-validation series; sim-tme-3d TME capstone (#195); 2D-math lift + `immune_3d`→`immune_spatial` rename + JSON schema_version (#220/#224); 3D trajectory snapshots + animated axial-slice GIF renderer (#193/#238); time-varying multi-dose pharmacokinetics (`dose_schedule` module — Constant/Bolus/MultiDose/Infusion/FromPk — wired into sim-tme-3d via `--dose-sweep` + the `--snapshot=multidose` preset, with the orphaned `tumor_pk` ODE finally bridged in via `FromPk`, #239); 3D performance and scalability work (`--bench` harness + within-condition rayon parallelism, byte-identical via position-independent per-cell RNG, 3.8x to 4.9x speedup on single large grids, dense 200³ measured at ~1.29 GB so sparse storage deferred to #254, #192); full-production byte-identity regression CI guarding sim-tme-3d's default-matrix `summary.json` SHA (#253); drug-tolerant persister cells (`persister` module + `PersisterConfig`, epigenetic ferroptosis tolerance acquired under drug exposure and reverting after clearance; off-by-default identity config keeps the matrix byte-identical, with a `--snapshot=persister` render overlay, #241); T-cell exhaustion (#243 Phase 1 — `immune_spatial::exhaustion_factor` + per-cell neighborhood `cumulative_kills`; sustained local killing suppresses further immune kills via `1/(1+rate·cumulative)`; `exhaustion_rate=0` default keeps the matrix byte-identical); clonal heterogeneity (#242 — `clonal` module: Voronoi subclone patches via an independent RNG + per-subclone iron/GPX4/MUFA perturbations, `summary.json` per-subclone kill reporting, K=1-identity default keeps the matrix byte-identical, `--snapshot=clonal` subclone overlay)
- ferroptosis-core library packaging for external use
- news source authentication pipeline (fetch, extract claims, verify, score, index)
- broader strategy review of alternative therapies and biological bottlenecks
- operational maturity: Phase 2 complete — CI (#126), figure traceability (#127), archival release tooling (#131); workspace `cargo fmt --check` gate added to Rust CI (#209/#236); off-PR sim-tme-3d production byte-identity regression workflow (#253)
- manuscript integrity: Phase 3 complete — immune coupling confidence (#130), structural uncertainty qualifiers (#137), PRISMA corpus flow diagram (#135), retrieval bias subsection (#140)
- sensitivity analyses: weight-sensitivity (#128), taxonomy-sensitivity (#133), PRCC global sensitivity (#134), and O2 cycling (#138) complete — pre-registered, run, results in manuscript
- test expansion (#139) complete — 19 invariant/integration tests added (schema, weight monotonicity, tagging correctness)

## Current Repo State

- local full-text corpus: 4,830 records
- abstract-only archive: 5,584 records
- mechanism taxonomy, evidence tiers, pathway-targets, biology-process tags, and resistant-state scaffolding are all active
- evidence tagging is improved but still incomplete (gold-set measured)
- tissue-of-origin and weighted-evidence layers are active
- diagnostic-therapy matching layer covers 6 chains across 4 modalities (radioligand, checkpoint, mRNA vaccine, oncolytic)
- manuscript: 112 pages (book format), 11 chapters + 3 appendices, 20 figures, ~36,700 words
- simulation suite: 11 binaries (incl. sim-tumor-pk + sim-tme-3d) + ferroptosis-core library (MIT licensed, 18 modules including `photosensitizer_pk`, `oxygen`, `ph`, `stromal`, `immune_spatial`, `dose_schedule`, `persister`, and `clonal`; v0.7.0 adds 3D radial-depth + 3D ROS-multiplier APIs alongside the 2D path #185–#186; v0.8.0 adds 3D radial O₂ field + zone-census #187; v0.9.0 adds 3D radial pH field + iron-release/ion-trap helpers #190; v0.10.0 adds 3D CAF-shielded boundary detection + adjacent-kill-rate #189; v0.11.0 adds 3D spatial immune coupling #188; sim-tme-3d capstone binary consumes all five #195; #220/#224 lift the TME config structs + 2D depth/pH/stromal helpers into the library and rename `immune_3d`→`immune_spatial`; v0.12.0 (#239) adds the `dose_schedule` time-varying-PK module + `biochem::exo_decay_factor`; v0.13.0 (#241) adds the `persister` drug-tolerant-persister module + `PersisterConfig`; v0.14.0 (#243) adds `immune_spatial::exhaustion_factor` (T-cell exhaustion) + `TumorGrid3D::coords` + `SpatialImmuneConfig::exhaustion_rate`; v0.15.0 (#242) adds the `clonal` module (Voronoi subclones + `ClonalConfig`/`SubclonePerturbation`); current crate version 0.15.0; unit-test count tracked in CI / `cargo test --workspace`) + Python bindings + 105 Python tests (pipeline smoke + figure traceability + invariant/integration + calibrate-extractor + ferroptosis-python bindings)
- news authentication pipeline: 5 scripts (fetch, extract claims, verify against PubMed, score credibility, build claim-centric index) implementing the 3-tier source framework from analysis/news-source-criteria.md
- simulation calibration: 8 targets documented (5 original + 3 new 3D self-consistency targets from #196 covering hypoxia-RSL3 collapse, immune SDT/RSL3 ratio, stromal boundary shielding); evaluate script operational; 3D validation infrastructure in place (`3d_validation_report.md` + upgrade-to-calibration path)
- drug penetration module: 3 tissue types, exponential Krogh approximation
- drug combination modeling: 4 drugs, pairwise Bliss synergy scoring with pathway traces
- tumor microenvironment: oxygen gradients, spatial immune zones (DAMP diffusion, T cell activation, anti-PD-1, T-cell exhaustion #243), LP overshoot multiplier, CAF-mediated stromal protection (GSH/MUFA supply), pH gradient (iron release + drug ion trapping)
- some landmark papers are known to be missing from the local full-text archive
- content provenance manifest (PROVENANCE.yaml) documents asset licensing and redistribution rights
- pinned Python environment (requirements-lock.txt, 32 packages) and Rust toolchain (rust-toolchain.toml, 1.96.0)
- contributor guide (CONTRIBUTING.md), citation metadata (CITATION.cff), and pytest in tracked dependencies
- Python CI workflow (.github/workflows/python-test.yml): Linux on PR/push, macOS weekly
- Rust CI workflow (.github/workflows/cargo-test.yml): `cargo test --workspace` + `cargo fmt --all --check` gate on PR/push (fmt pinned to the 1.96.0 toolchain, #209/#236)
- sim-tme-3d production regression workflow (.github/workflows/sim-tme-3d-regression.yml): weekly + on-demand full 60³×180 run asserting `summary.json`'s SHA-256 against a checked-in hash on the pinned 1.96.0 toolchain (#253)
- figure traceability index (FIGURES.yaml) mapping all 23 figures to generators, inputs, and types
- archival release tooling (.zenodo.json metadata template, scripts/generate_release_manifest.py for SHA256 manifest + filtered archive)

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
README.md                                    repo-level framing and entry points
article/drafts/v1.{md,tex,pdf}               manuscript drafts
analysis/evidence-coverage-audit.md          evidence-tier coverage and guardrails
analysis/taxonomy-rerun-notes.md             taxonomy/query caveats after reruns
analysis/pathway-target-audit.md             pathway-target coverage
analysis/landmark-corpus-gaps.md             known missing papers that distort claims
corpus/INDEX.jsonl                           master index
scripts/                                     Python pipeline
simulations/                                 Rust simulation work
simulations/ferroptosis-python/              Python bindings (PyO3)
simulations/ferroptosis-ffi/                 C FFI bindings (PhysiCell integration)
tags/                                        precomputed tag indexes
article/book-outline.md                      frozen book outline and chapter contracts
article/AUTHORING.md                         writing rules and heading conventions
news/                                        news source scaffolding (issue #99)
PROVENANCE.yaml                              content provenance and redistribution rights
CONTRIBUTING.md                              contributor setup, testing, and PR guide
CITATION.cff                                 citation metadata (renders GitHub "Cite" button)
requirements-lock.txt                        pinned Python dependency versions
FIGURES.yaml                                 figure-to-script traceability index (23 figures)
.github/workflows/python-test.yml            Python CI (Linux PR/push, macOS weekly)
.github/workflows/cargo-test.yml             Rust CI (cargo test + cargo fmt --check gate)
.github/workflows/sim-tme-3d-regression.yml  sim-tme-3d production byte-identity regression (weekly + manual)
.zenodo.json                                 Zenodo deposit metadata template
scripts/generate_release_manifest.py         SHA256 manifest + filtered archive builder
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
