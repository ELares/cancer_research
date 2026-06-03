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

- manuscript drafting and revision (Parts I–V substantially complete: ~115 pages, 11 chapters + 3 appendices)
- corpus fetching, enrichment, tagging, and indexing
- evidence-tier audits and coverage caveats (gold-set evaluation: 46% exact, 96% precision, 55% recall)
- taxonomy and search refinement
- pathway-target and resistant-state analysis
- diagnostic-to-therapy chain extraction (6 chains, 129 articles mapped)
- tissue-of-origin analysis layer (5 tissue categories, 62% coverage)
- simulation work: ferroptosis biochemistry, drug penetration, calibration, photosensitizer PK (drug-light-interval scaling, saturating distribution phase, relative singlet-O₂ yield, FromStr-based clap CLI integration in sim-spatial), 3D spheroid scaffolding (TumorGrid3D #185, signed radial depth + 3D energy-physics dispatcher #186, 3D radial O₂ gradient + zone-census #187, 3D radial pH gradient + iron/ion-trap helpers #190, 3D CAF-shielded boundary detection #189, 3D spatial immune coupling + DAMP diffusion #188) for the #185–#197 spheroid-validation series; sim-tme-3d TME capstone (#195); 2D-math lift + `immune_3d`→`immune_spatial` rename + JSON schema_version (#220/#224); 3D trajectory snapshots + animated axial-slice GIF renderer (#193/#238); time-varying multi-dose pharmacokinetics (`dose_schedule` module — Constant/Bolus/MultiDose/Infusion/FromPk — wired into sim-tme-3d via `--dose-sweep` + the `--snapshot=multidose` preset, with the orphaned `tumor_pk` ODE finally bridged in via `FromPk`, #239); 3D performance and scalability work (`--bench` harness + within-condition rayon parallelism, byte-identical via position-independent per-cell RNG, 3.8x to 4.9x speedup on single large grids, dense 200³ measured at ~1.29 GB so sparse storage deferred to #254, #192); full-production byte-identity regression CI guarding sim-tme-3d's default-matrix `summary.json` SHA (#253); drug-tolerant persister cells (`persister` module + `PersisterConfig`, epigenetic ferroptosis tolerance acquired under drug exposure and reverting after clearance; off-by-default identity config keeps the matrix byte-identical, with a `--snapshot=persister` render overlay, #241; #262 replaces the acquire-xor-revert either-or with a competing-rate `step` so sustained sub-saturating drug reaches a sub-cap equilibrium, + persister calibration provenance); T-cell exhaustion (#243 Phase 1 — `immune_spatial::exhaustion_factor` + per-cell neighborhood `cumulative_kills`; sustained local killing suppresses further immune kills via `1/(1+rate·cumulative)`; `exhaustion_rate=0` default keeps the matrix byte-identical; #264 Phase 2 adds a Treg/MDSC suppressor field — a second diffusing field, replenished at perivascular-or-heuristic niche sources, that scales immune kill DOWN via `1/(1+strength·field)`; off-by-default byte-identical, anti-PD-1 Treg-depletion validation test + `--snapshot=suppressor` niche overlay); #264 Phase 3 generalizes the single PD-1 brake to a multi-checkpoint panel (PD-1/CTLA-4/LAG-3/TIM-3, each independently drug-modulated; combined brake 1−Π(1−residualᵢ)) to model anti-PD-1 + anti-CTLA-4 combinations, off-by-default byte-identical with a dual-blockade validation test + `--snapshot=checkpoint`); cross-layer composition validation (#278 — pairwise-distinct realism-layer RNG-seed guard + multi-layer integration tests (clonal×suppressor×checkpoints, spheroid×vasculature) asserting deterministic composition + coherent per-layer metrics, plus a `--snapshot=combined-realism` kitchen-sink preset; all-off matrix stays byte-identical); clonal heterogeneity (#242 — `clonal` module: Voronoi subclone patches via an independent RNG + per-subclone iron/GPX4/MUFA perturbations, `summary.json` per-subclone kill reporting, K=1-identity default keeps the matrix byte-identical, `--snapshot=clonal` subclone overlay; #266 makes the GPX4 axis durable — `gpx4_mul` now scales the static `cell.nrf2` setpoint too, not just the initial `state.gpx4`, so a GPX4/antioxidant-low subclone stays differentiated for the whole run instead of relaxing to the shared setpoint; #266 item 3 adds spatial clonal expansion — dead tumor sites repopulated from living Moore-neighbors so resistant subclones grow their territory, with `subclone_kills.initial_tumor` reporting the shift; off-by-default (repopulation_rate=0) byte-identical); explicit 3D vasculature (#191 — `vasculature` module: random internal vessel seeds via an independent RNG + per-cell `exp(-dist_to_nearest_vessel/λ)` unified supply replacing the edge-distance O2 proxy and scaling drug delivery on all paths (SDT/PDT exo, RSL3 constant knockdown, dosed RSL3 availability), off-by-default byte-identical, `vascular_hypoxic_fraction` reporting + `--snapshot=vasculature` O2-supply overlay; #268 replaces the brute-force nearest-vessel scan with a `VesselIndex` uniform-grid spatial index — exact/byte-identical (Chebyshev-shell search with a provable early-stop, verified bit-for-bit vs brute force), ~O(cells) instead of O(cells×vessels), so vasculature scales to patient-size grids the #272 slab coupling needs; #268 also adds `VesselTopology::{Random,Fractal}` + `place_vessels_fractal_3d` — a fractal-branching vessel tree (trunks enter from the periphery and bifurcate inward via BFS with high, tumor-like angle/length jitter and occasional dead ends, per Baish & Jain 2000 PMID 10919633, capped at the same point-count target as the random network — matching raw point count, NOT effective coverage) producing hierarchical-but-chaotic perfusion with avascular gaps; at near-equal point count the fractal network leaves a higher hypoxic fraction than uniform-random seeding (a clustering-coverage effect: 1-cell-spaced branch points cover far less unique volume than scattered points — qualitative, not density-controlled), slab geometry ignores topology, off-by-default topology=Random byte-identical, no preset overlay yet); 3D spheroid radial biochemistry (#197 — `spheroid` module: radial phenotype re-assignment (glycolytic rim / OXPHOS mid / persister core) via an independent RNG + position-dependent MUFA/GSH/iron gradients, run under the new `Params::spheroid()` partial-MUFA context, off-by-default byte-identical, `--snapshot=spheroid` phenotype overlay; #270 makes the zone thresholds VOLUME-based (compare frac³) and grounds the defaults in the Browning 2021 eLife limiting structure (necrotic core 0.73 of radius → 0.39 of volume; rim begins 0.90 of radius), fixing the radial-threshold inversion that gave a ~4%-volume core / ~71%-volume rim — geometry now literature-grounded, biochem gradient strengths still placeholders; #270 makes the radial MUFA durable via a per-cell `Cell::mufa_cap` so the rim-vs-core spread persists instead of relaxing to the uniform M_ss); patient-scale slab geometry (#240 — `slab` module + `TumorGrid3D::generate_slab`: an all-tumor block with a planar depth-graded `exp(-depth/λ)` supply (the 1-D analog of the radial O2 field, +z face vessel-proximal, −z reflective/no-flux) replacing O2 and scaling drug on all paths, `SlabConfig::patient_deep()`/`surface()` (10 mm virtual tumor), a `scale_interpretation` output string, off-by-default byte-identical, `--snapshot=slab` depth-gradient viz; addresses the in-vitro-vs-patient scale gap — a 4 mm-deep slab kills <20% of the spheroid under the same SDT, though the depth-collapse magnitude is an uncalibrated first-order Krogh approximation; #272 couples the slab with explicit vasculature (#191) — `place_vessels_in_slab_3d` scatters vessels uniform-in-box across the all-tumor block and the per-cell supply becomes `max(planar_depth, nearest_vessel)`, so internal vessels deliver drug to focal deep pockets and RAISE deep killing, i.e. a slab is less therapy-resistant at depth than the planar-only model implies wherever real vasculature reaches; off-by-default byte-identical, gated on both overrides, `--snapshot=slab-vessels` overlay)
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
- manuscript: ~115 pages (book format), 11 chapters + 3 appendices, 24 figures, ~36,700 words
- simulation suite: 11 binaries (incl. sim-tumor-pk + sim-tme-3d) + ferroptosis-core library (MIT licensed, 21 modules including `photosensitizer_pk`, `oxygen`, `ph`, `stromal`, `immune_spatial`, `dose_schedule`, `persister`, `clonal`, `vasculature`, `spheroid`, and `slab`; v0.7.0 adds 3D radial-depth + 3D ROS-multiplier APIs alongside the 2D path #185–#186; v0.8.0 adds 3D radial O₂ field + zone-census #187; v0.9.0 adds 3D radial pH field + iron-release/ion-trap helpers #190; v0.10.0 adds 3D CAF-shielded boundary detection + adjacent-kill-rate #189; v0.11.0 adds 3D spatial immune coupling #188; sim-tme-3d capstone binary consumes all five #195; #220/#224 lift the TME config structs + 2D depth/pH/stromal helpers into the library and rename `immune_3d`→`immune_spatial`; v0.12.0 (#239) adds the `dose_schedule` time-varying-PK module + `biochem::exo_decay_factor`; v0.13.0 (#241) adds the `persister` drug-tolerant-persister module + `PersisterConfig`; v0.14.0 (#243) adds `immune_spatial::exhaustion_factor` (T-cell exhaustion) + `TumorGrid3D::coords` + `SpatialImmuneConfig::exhaustion_rate`; v0.15.0 (#242) adds the `clonal` module (Voronoi subclones + `ClonalConfig`/`SubclonePerturbation`); v0.16.0 (#191) adds the `vasculature` module (explicit vessel network + `VasculatureConfig`); v0.17.0 (#197) adds the `spheroid` module (radial phenotype/MUFA/GSH/iron + `SpheroidConfig`) + `Params::spheroid()`; v0.18.0 (#240) adds the `slab` module (patient-scale all-tumor block + planar depth-graded supply, `SlabConfig`, `scale_interpretation`, `KROGH_LAMBDA_UM`) + `TumorGrid3D::generate_slab`; v0.19.0 (#270) adds per-cell durable MUFA (`Cell::mufa_cap`, threaded into `update_mufa_protection`, set radially by the spheroid so position-dependent MUFA no longer relaxes to the uniform M_ss); v0.20.0 (#264 Phase 2) adds the Treg/MDSC immunosuppressor field to `immune_spatial` (`SuppressorConfig` + `suppressor_kill_multiplier` + `suppressor_source_mask_3d`); v0.21.0 (#264 Phase 3) adds the multi-checkpoint immune brake (`CheckpointPanel` + `Checkpoint`, `combined_brake` = 1 − Π(1 − brakeᵢ·(1−drug_effᵢ)) over PD-1/CTLA-4/LAG-3/TIM-3); v0.22.0 (#262) adds the persister competing-rate `step` (acquisition + reversion act simultaneously ⇒ sub-cap equilibrium under sustained sub-saturating drug); v0.23.0 (#266 item 3) adds spatial clonal expansion (`repopulate_dead_sites_3d` + `ClonalConfig::repopulation_rate`/`with_repopulation`); v0.24.0 (#272) adds slab+vasculature coupling (`vasculature::place_vessels_in_slab_3d` uniform-in-box vessel placement, so a patient-scale slab can carry internal vessels whose proximity supply combines element-wise-MAX with the planar depth gradient); v0.25.0 (#268) adds the `VesselIndex` uniform-grid spatial index for `vessel_supply_field`'s nearest-vessel lookup (exact/byte-identical to the former brute force, ~O(cells) instead of O(cells×vessels), so vasculature scales to patient-size grids — 100³/1M cells ≈ 105 ms); v0.26.0 (#268 + spheroid zone fix) adds `VesselTopology::{Random,Fractal}` + `place_vessels_fractal_3d` (a fractal-branching vessel tree — trunks enter from the periphery and bifurcate inward with high, tumor-like variability per Baish & Jain 2000 PMID 10919633, producing hierarchical-but-chaotic perfusion with avascular gaps, capped at the same point-count target as random (raw point count, NOT effective coverage — clustered 1-cell-spaced branch points cover far less unique volume than the same number of scattered points, so the higher hypoxic fraction is a clustering-coverage effect read qualitatively, not a density-controlled result); slab geometry ignores topology; off-by-default topology=Random keeps the matrix byte-identical) and makes `spheroid::radial_phenotype` volume-based (compares `frac.powi(3)` to volume-fraction thresholds with `literature()` defaults `glycolytic_frac:0.73`/`oxphos_frac:0.39` from Browning 2021 eLife so rim/mid/core zones reflect equal-volume shells, not equal-radius bands; spheroid layer is opt-in so byte-identical); current crate version 0.26.0; unit-test count tracked in CI / `cargo test --workspace`) + Python bindings + 119 Python tests (pipeline smoke + figure traceability + depth-kill physics-constant guard + flagship-figure data guard + invariant/integration + calibrate-extractor + ferroptosis-python bindings)
- news authentication pipeline: 5 scripts (fetch, extract claims, verify against PubMed, score credibility, build claim-centric index) implementing the 3-tier source framework from analysis/news-source-criteria.md
- simulation calibration: 8 targets documented (5 original + 3 new 3D self-consistency targets from #196 covering hypoxia-RSL3 collapse, immune SDT/RSL3 ratio, stromal boundary shielding); evaluate script operational; 3D validation infrastructure in place (`3d_validation_report.md` + upgrade-to-calibration path); `simulations/calibration/CALIBRATION_STATUS.md` is the single honest per-layer accounting (calibrated / partially anchored / self-consistency / uncalibrated) + calibration roadmap, consolidating the previously-scattered caveats — the suite is broad but mostly uncalibrated (mechanistic scaffolding, excluded from the manuscript's quantitative claims)
- manuscript honesty pass (integrity audit): the 3D suite was absent from the manuscript + several §8.4 "the model cannot do X" caveats had gone stale (clonal/persister DO model adaptive evolution; suppressor/checkpoints DO model immune dampening; vasculature DOES replace the edge-O2 proxy — all uncalibrated); fixed by a consolidated §8.4 "3D simulation suite" honesty paragraph (names the suite, its uncalibrated status, the deliberate exclusion from quantitative claims, and the 104:1→4:1 immune-ratio shrinkage in 3D) + drift fixes (version 0.2.0→0.25.0, "ten modules"→21, missing sim-tme-3d/sim-tumor-pk binaries); `v1.md` edited then `v1.tex` regenerated (compiles, ~115 pp)
- manuscript scientific review (deep external+internal cross-check, `analysis/manuscript-scientific-review.md`): five-lens review (supporting lit, contradicting lit, math/physics, fresh-eyes coherence, figure gaps) with every new PMID personally verified. Math/physics confirmed sound (Beer-Lambert, Krogh, acoustic attenuation, Bliss 1.99×, ferroptosis ODEs); fixed the `iron_diffusion_coeff` mis-citation (was Jacques-2013-optics → tortuosity-reduced estimate). Key science: the central "physical ROS BYPASSES hypoxia" claim was over-extended vs the SDT field (SDT ROS is widely O2-dependent; lead clinical agent SONALA-001 is Type II/O2-dependent), so §7.1/§7.2/§7.5/§8.4/§10.1 were rebalanced (direction survives, "bypass"+magnitude do not) with 6 verified citations added on both sides (Hubbi/Dang 2026 PMID 41932308 supports the mechanism in PDAC but shows it combines GSH-defense + mito-ROS-suppression; Zou 2019 PMID 30962421 shows hypoxia SENSITIZES ccRCC to RSL3; Hayashi 2020 PMID 33288764 DAMP-balance≠quantity; Dai 2020 PMID 33311482 ferroptosis accelerates PDAC in vivo). Figure roadmap documented; figures not yet generated
- drug penetration module: 3 tissue types, exponential Krogh approximation
- drug combination modeling: 4 drugs, pairwise Bliss synergy scoring with pathway traces
- tumor microenvironment: oxygen gradients (edge-distance proxy, explicit internal vessel network #191, or patient-scale planar depth gradient #240), spatial immune zones (DAMP diffusion, T cell activation, anti-PD-1, T-cell exhaustion #243, Treg/MDSC suppressor field + multi-checkpoint PD-1/CTLA-4/LAG-3/TIM-3 brake #264), LP overshoot multiplier, CAF-mediated stromal protection (GSH/MUFA supply), pH gradient (iron release + drug ion trapping)
- some landmark papers are known to be missing from the local full-text archive
- content provenance manifest (PROVENANCE.yaml) documents asset licensing and redistribution rights
- pinned Python environment (requirements-lock.txt, 32 packages) and Rust toolchain (rust-toolchain.toml, 1.96.0)
- contributor guide (CONTRIBUTING.md), citation metadata (CITATION.cff), and pytest in tracked dependencies
- Python CI workflow (.github/workflows/python-test.yml): Linux on PR/push, macOS weekly
- Rust CI workflow (.github/workflows/cargo-test.yml): `cargo test --workspace` + `cargo fmt --all --check` gate on PR/push (fmt pinned to the 1.96.0 toolchain, #209/#236); the `test` job also builds `ferroptosis-ffi` and compiles + runs the C harness (`tests/test_ffi.c`) against the cdylib, verifying the PhysiCell C-ABI per-PR (#295 item 4)
- sim-tme-3d production regression workflow (.github/workflows/sim-tme-3d-regression.yml): weekly + on-demand full 60³×180 run asserting `summary.json`'s SHA-256 against a checked-in hash on the pinned 1.96.0 toolchain (#253)
- figure traceability index (FIGURES.yaml) mapping all 28 figures to generators, inputs, and types (#285 quantitative simulation figures: hypoxia kill-collapse, Bliss synergy, treatment-timing window, depth-kill curves, and the flagship 2×2 resistance-mechanism asymmetry — manuscript figure 24)
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
analysis/manuscript-scientific-review.md     external+internal scientific review: verified support/contradiction, math check, reframe rationale, figure roadmap
simulations/calibration/CALIBRATION_STATUS.md per-layer calibration tiers + roadmap (the honest sim accounting)
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
FIGURES.yaml                                 figure-to-script traceability index (28 figures)
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
