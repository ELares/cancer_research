# Research next steps

The [mission](../MISSION.md) puts the normalized, measurable cancer evidence
base first and uses simulations to test mechanisms suggested by it. The next
increments should make existing claims reproducible and independently testable.
This order follows the [research review](RESEARCH_REVIEW_2026-09-18.md); it does
not require another simulation mechanism or a larger unmeasured corpus.

1. **Calibration support is corrected; joint inference remains unresolved.**
   The regenerated fits now use one supported cohort per compound, with input
   hashes, grids, and exclusions recorded. The former 100 µM erastin target
   exceeded every recorded maximum. The corrected joint run accepted only
   4 of 40,000 draws, below the unchanged minimum of 20; its earlier interval
   claims are superseded. See the [methods and results](CALIBRATION_DOSE_SUPPORT.md).
   **Next deliverable:** compare a more efficient sampler with the fixed prior,
   target construction, reference vector, and tolerance; report independent-run
   stability, accepted counts, and held-out errors. Do not widen the criterion
   to fill a quota. A favorable result is not a completion requirement.
   **Complete when:** the joint inference meets its documented sampling minimum
   and repeated-run diagnostics, or the remaining sampling/model mismatch is
   documented without posterior claims. Independent validation is a separate
   requirement below.

2. **Census input protection is implemented in six analysis entry points.**
   The shared [streaming helper](../scripts/census_input.py) honors
   `FERRO_ATLAS_ROOT`, validates shard stride, and rejects missing or empty
   inputs before replacing reports. Valid records with zero scientific matches
   remain valid results; `--render-only` works without raw data. Regression
   coverage includes preservation of both outputs, malformed later records,
   and render failures. See the [acquisition instructions](../analysis/atlas-README.md).
   **Next deliverable:** audit remaining census scanners and adopt the helper
   where they share these input semantics, preserving their sampling behavior.
   **Complete when:** each migrated entry point distinguishes unavailable input
   from an observed zero and its offline rendering remains reproducible.

3. **Measure the immune comparison directly.** Deliver a passive report from
   matched SDT and RSL3 scenarios using their existing per-condition seeds,
   geometry, oxygen conditions, and immune configuration. Record ferroptotic deaths, death-time
   LP, release-time LP, completed releases, and deaths whose release period is
   censored by the simulation horizon. Report terminal DAMP additions
   separately: additions after the last immune update cannot explain kills
   already counted. Record actual immune-eligible living tumor cells after
   the activation delay and DAMP threshold, with an explicit denominator for
   unique cells versus cell-time opportunities, alongside total immune kills.
   **Complete when:** event counts reconcile, every conditional mean names its
   population, empty populations are explicit, and instrumentation preserves
   existing simulation outcomes and random draws. The historical 104:1 total
   immune-kill ratio supplies no per-cell DAMP ratio or replacement experimental
   threshold. Start in [sim-tme-3d](../simulations/sim-tme-3d/src/main.rs), its
   [immune eligibility code](../simulations/ferroptosis-core/src/immune_spatial.rs),
   and [headline_sensitivity.py](../scripts/headline_sensitivity.py).

4. **Validate one observable against independent raw assay replicates.** Deliver
   a validation report with a named source/accession, raw replicate measurements,
   assay endpoint and timing, tested dose range, and a prespecified evaluation
   split and error criterion. ML210 is held out by compound in the current fit,
   but shares cell lines with ML162 within the same screen; that is limited
   transfer evidence, not an independent assay validation. **Complete when:**
   parameters are frozen before evaluation, shared cell lines and experiments
   are disclosed, uncertainty uses the independent replicate unit, and the
   observed error is reported whether it passes or fails the criterion. If
   usable raw data cannot be obtained, publish the specific availability gap
   and leave validation pending. Extend the provenance in
   [calibration-targets-ctrpv2.md](../analysis/calibration/calibration-targets-ctrpv2.md)
   and the [calibration ledger](../simulations/calibration/CALIBRATION_STATUS.md).

   **Acquisition candidate checked in September 2026:**
   [Lee et al., Nature Metabolism (2024)](https://www.nature.com/articles/s42255-024-00974-4)
   supplies three normalized observations for SK-Hep1 treated with 0.1 µM
   ML210 for 24 hours, with vehicle and ferrostatin controls. This single-dose
   endpoint is a bounded candidate; the erastin endpoint uses 6 µM for 72 hours
   and a different assay. See the [acquisition findings](INDEPENDENT_ASSAY_CANDIDATE.md)
   for source links, workbook hash and cell ranges, and descriptive summaries.
   The ML210 curves remain unresolved: both workbook and figure say “nmol”,
   two columns conflict with the caption's three replicates, and zero-dose
   control encoding is implicit. Raw readings and experiment identities are
   unavailable in the reviewed materials. Model-to-assay and exposure-time
   mappings must be specified before comparison; no validation set has been
   ingested and independent raw-replicate validation remains pending.

5. **Measure evidence quality before promoting discovery candidates.** Deliver
   a bounded, independently adjudicated census sample with sentence-attributed
   claims, stable entity identifiers, and explicit abstract/full-text and census
   stream coverage. Report precision, recall where a recall denominator exists,
   intervals, and annotator agreement. Then evaluate any proposed discovery
   change against the same temporal split, seeds, and candidate pool as the
   degree/popularity and random baselines. **Complete when:** the quality report
   exposes its sampling limits and the discovery report gives paired differences
   with uncertainty; lack of improvement remains a reportable result. Predicting
   later publications and finding biologically useful overlooked connections
   need separate evaluation claims. Start with
   [evaluate_evidence_v2.py](../scripts/evaluate_evidence_v2.py),
   [atlas_contradiction_quality.py](../scripts/atlas_contradiction_quality.py),
   [atlas_discovery_eval.py](../scripts/atlas_discovery_eval.py), and the existing
   [degree-bias](../analysis/atlas-discovery-degree-bias.md) and
   [headroom](../analysis/atlas-discovery-headroom.md) reports.

The open [trial-publication-gap PR #866](https://github.com/ELares/cancer_research/pull/866)
and [long-`Retry-After` PR #867](https://github.com/ELares/cancer_research/pull/867)
are separate work already in progress. Review and integrate those changes on
their own branches rather than duplicating them in these increments.
