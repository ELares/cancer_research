# Research review — 18 September 2026

Review base: `758ba464cf7c6d1b8a0914ee2bc9cf3d4b560367` on `main`.

## What this repository can support

The repository connects three substantial systems: a literature census, a
mechanistic simulation suite, and a manuscript that interprets both. Its most
useful next work is to strengthen the connections between an input, a computed
result, and the claim made from that result. Adding another uncalibrated
mechanism would contribute less than making these connections reliable.

The corpus has several distinct denominators. The manuscript reports 4,403,994
MeSH-indexed census records plus 783,271 text-recovered records; the frozen
4,830-record retrieval is a historical comparison, not the census. The checkout
includes 10,423 full-text/abstract Markdown files and committed analyses, but
the bulk census streams, million-article open-access layer, and identity
database are acquired separately. A successful checkout does not imply that
every analysis can be recomputed from local inputs.

The Rust workspace contains 16 crates, including 13 simulation binaries.
Single-cell and 2D results underpin the legacy quantitative chapters; the 3D
engine provides a wider, largely optional set of mechanisms. Its default-output
SHA is a useful regression check, not independent biological validation. The
eight targets in `simulations/calibration/targets.yaml` likewise need their
provenance read: matching a model-derived target is self-consistency, not a fit
to an independent experiment. The more detailed evidence accounting belongs in
[`CALIBRATION_STATUS.md`](../simulations/calibration/CALIBRATION_STATUS.md).

The manuscript source is [`article/drafts/v1.md`](../article/drafts/v1.md).
LaTeX and figures are generated artifacts. The review traced selected headline
claims to the uncertainty reports, simulation equations, figure generators,
and experimental protocol. It did not re-read millions of papers, independently
validate every citation, or rerun the full census. Large local data stores were
left unchanged; ingestion regressions were reproduced with temporary fixtures.

## Reproduced issues and repairs

| Area | Evidence | Repair |
| --- | --- | --- |
| Tumor pharmacokinetics | The solver integrated a minute before recording a row labelled with the start of that minute. Zero initial tumor concentration was reported as nonzero at time zero; a future plasma input could appear in an earlier tissue sample. | Record the state at its labelled time before advancing. Check initial conditions, an analytic decoupled ODE, and causality of a delayed plasma input. |
| PK input validation | CSV parsing accepted `NaN` and infinity; normalization could silently turn invalid exposure into zero. Zero integration substeps also produced zero exposure. | Reject nonfinite measured inputs and zero substeps explicitly. |
| Living-review deduplication | The producer writes dated `index.jsonl`, while the identity builder looked only for Markdown. A successful build could miss every living-review identifier. | Index plain JSONL as metadata, alongside existing Markdown and gzip sources; test PMID, PMC ID, DOI, unchanged-file skipping, and changed-file rescanning. |
| Census report preservation | With missing census shards, `census_mechanism_profile.py` returned success and replaced published results with a zero-record report. | Fail before writing for missing or empty inputs; retain valid zero-match results, support `FERRO_ATLAS_ROOT`, and keep offline `--render-only` usable. |
| Fetch admission rate | The scheduled macOS CI run [34857341665](https://github.com/ELares/cancer_research/actions/runs/34857341665) exposed bunched requests. Reserving future slots before sleeping let late workers consume several expired slots together. A deterministic clock reproduced two consecutive grants at the same time. | Serialize admission waiting and set the next deadline from the actual wakeup. Tests exercise late and early wakeups. Network operations remain concurrent after admission; this is an admission-spacing guarantee, not a hard real-time socket-start guarantee. |
| Engine-time audit | Its declaration lookup stopped after 12 lines, so adding API documentation could silently remove a time-bound solver from the report. | Follow attached documentation to its declaration, stopping at unrelated code; test both long documentation and an unrelated declaration. Regenerate the audit and source-line provenance. |
| FFI pointer contracts | Workspace Clippy exposed public safe Rust functions that dereference caller-supplied raw pointers. Null checks do not establish validity or exclusive access. | Mark the three pointer-consuming exports `unsafe` for Rust callers and document their validity, lifetime, and aliasing contracts. C declarations and ABI remain unchanged. |
| Combination interpretation | The committed Bliss-ratio report has a 95% prior-predictive interval of `[1.000, 5.242]` and full sampled range `[0.953, 7.800]`. Those values do not support guaranteed supra-additivity. The text also called a ratio an excess and treated failure of a synergy threshold as proof against pathway independence. | Distinguish the Bliss ratio, Bliss excess, and Chou–Talalay CI; qualify the uncertainty and the inference from a failed experiment. Preserve the historical P1 preregistration and its decision thresholds, with explicit clarification in the operational protocol. |
| Immune-coupling interpretation | The manuscript and diagram compared LP near 20 versus 7.8 as per-dead-cell means, inferring a 2.6-fold DAMP ratio. The relevant spatial condition requires LP above 10 for death and permits further nonnegative LP accumulation afterward. The cited outputs do not establish the claimed conditional means. | Remove the unsupported per-dead-cell comparison. Describe the actual assumed LP-to-DAMP mapping and separate total death density from per-cell immunogenicity. Update the diagram and generated caption with the prose. |

### Scientific interpretation of the combination correction

For the illustrative single-agent kills of 0.40 and 0.037, Bliss independence
predicts `0.40 + 0.037 - 0.40 * 0.037 = 0.4222`. This already exceeds either
single-agent kill without synergy. A modeled combination kill of 0.841 gives an
observed/expected ratio of about 1.99 and an excess of about 0.419, or 41.9
percentage points. A Chou–Talalay combination index requires dose-response
information and cannot be inferred from that ratio.

Brequinar also cannot be presented as an interchangeable selective FSP1 reagent.
The revised protocol distinguishes its DHODH target from concentration-dependent
FSP1 effects and calls for target-engagement and genetic controls. The primary
sources are [Mao et al. (2021)](https://doi.org/10.1038/s41586-021-03539-7),
[Mishima et al. (2023)](https://doi.org/10.1038/s41586-023-06269-0), and the
[reply](https://doi.org/10.1038/s41586-023-06270-7). The CI interpretation is
supported by [Chou (2010)](https://pubmed.ncbi.nlm.nih.gov/20068163/).

These corrections improve the inference drawn from existing evidence. They do
not add a new experimental validation or fit the model to biological data.

## Verification

The regressions use isolated fixtures, analytic expectations, and committed
evidence artifacts. The new analytic PK test fails against the original solver;
the delayed/early wakeup tests fail against the original rate limiter.

The default 24-condition `sim-tme-3d` production result retains SHA256
`973d9443bd05851721d53d9a15e53dc828517bba4bf307eb1c4f0b1f60d27289`.
Separate seeded old/new PK comparisons changed at most two deaths per 10,000
cells in the checked spatial conditions. This is a bounded regression
observation, not a claim that time alignment is immaterial for every input.

Required local checks are the full Python suite, full Rust workspace tests,
Rust formatting and Clippy, manuscript regeneration, and release-manifest
freshness. Use the pinned Rust toolchain and Python dependencies. Build output
and temporary experiment logs belong outside the checkout.

The full local Python run passed 2,143 tests with 22 skips; the Rust workspace
passed 656 tests with one ignored. Formatting, Clippy, and the C integration
smoke test also passed. Existing Clippy and Python dependency warnings remain.

The corrected conceptual diagram was rendered with Graphviz on the local Linux
host using Ubuntu's packages extracted into a temporary user directory (no
system installation). Only its PDF/PNG pair was regenerated and visually
checked; the existing helper removed PDF dates. The manuscript LaTeX was
regenerated, but the full manuscript PDF was not compiled locally because a
TeX toolchain was unavailable.

## Next work, in order

1. **Make all census generators distinguish unavailable data from a negative
   result.** The mechanism-profile repair covers one important path; other
   generators still assume local paths. Apply the same missing-input and
   valid-zero-match fixtures, then expose one common input preflight. A
   fresh checkout must never erase a published finding because data was not
   downloaded.
2. **Measure the immune statistic the manuscript needs.** Add matched spatial
   outputs for LP at death, LP at grace-period end, DAMP per dead cell, and
   the eligible immune-kill population. Report these by treatment, alongside
   death counts and uncertainty. This would separate the density explanation
   from the assumed per-cell contribution; external DAMP measurements would
   still be needed to calibrate the latter.
3. **Prioritize independent validation over new simulation layers.** Start
   with one existing claim, a named external assay/dataset, a prespecified
   fit/evaluation split, and observable-specific error criteria. Do not count
   recovery of `targets.yaml` values or unchanged production hashes as
   held-out predictive performance. Extend calibration checks that currently
   inspect committed JSON to recompute the relevant metric against the current
   compiled engine; artifact integrity alone cannot detect a changed solver.
   The combination dose matrix is a concrete
   proposed experiment; completing its protocol is not performing it.
4. **Extend claim-to-artifact checks beyond the repaired sections.** Audit
   chapter cross-references and summaries against their current destinations;
   distinguish absolute kill, superiority to a single agent, interaction
   relative to a reference model, and clinical efficacy. The present review
   corrected selected high-impact claims, not every paragraph of the book.
5. **Make data availability part of reproducibility reporting.** Provide an
   explicit inventory of committed summaries, acquired raw streams, input
   hashes, and commands that can run offline. Keep census, frozen retrieval,
   and living-update denominators visible wherever results are compared.

The open trial-publication-gap work and long-`Retry-After` fix were already in
separate branches during this review and were not duplicated here.
