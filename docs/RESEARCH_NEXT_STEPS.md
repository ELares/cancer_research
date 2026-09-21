# Research next steps

The [mission](../MISSION.md) puts the normalized, measurable cancer evidence
base first and uses simulations to test mechanisms suggested by it. The next
increments should make existing claims reproducible and independently testable.
This order follows the [research review](RESEARCH_REVIEW_2026-09-18.md); it does
not require another simulation mechanism or a larger unmeasured corpus.

1. **Joint sampling passes its biological-run screens; broader geometry challenges expose a limitation.**
   The regenerated fits now use one supported cohort per compound, with input
   hashes, grids, and exclusions recorded. The former 100 µM erastin target
   exceeded every recorded maximum. The corrected joint run accepted only
   4 of 40,000 draws, below the unchanged minimum of 20; its earlier interval
   claims are superseded. See the [methods and results](CALIBRATION_DOSE_SUPPORT.md).
   **Sampling experiment completed:** a
   [plan committed before production](JOINT_SAMPLING_PLAN.md) evaluated three
   independently fitted, frozen importance proposals with the same prior,
   targets, reference vector, and final tolerance. The
   [results](../analysis/calibration/joint-importance-sampling.md) accepted
   1, 4, and 4 of 8,192 production attempts, with ESS 1.0, 3.6, and 4.0.
   Every run failed all four per-run adequacy screens; parameter and held-out
   error stability also failed. All attempts, weights, code/data provenance,
   and report reconstruction are archived. No pooled posterior is published.
   This documents an insufficient sampling design, not an empty target or
   established biological model failure.
   **Synthetic coverage prerequisite added:** the first bounded resample/move
   design passed three boundary-slab runs but failed all three separated-region
   runs. One failed run had ESS 1,101 while missing a region containing 25% of
   the target mass. Its [failure analysis](JOINT_RESAMPLE_FAILURE.md) and
   [replayable study](../analysis/calibration/proposal-synthetic-validation.md)
   are retained; the biological driver rejects this failed prerequisite.
   **Revised sampling experiment completed:** the
   [separately frozen local-move plan](JOINT_RESAMPLE_LOCAL_PLAN.md) retained all
   18 development trials, then passed all six prospective synthetic runs before
   biological production. The [three biological runs](../analysis/calibration/joint-resample-sampling.md)
   accepted 1,269, 1,323 and 1,236 of 8,192 independent attempts, with ESS 728.4,
   844.1 and 824.1. All per-run and across-run screens passed; all six islands
   reached the unchanged tolerance. The largest median span was 2.14% of prior
   width, the largest tail-quantile span 2.72%, and held-out median-RMSE span
   0.00262. Every evaluated production curve, weight, decision and cost is
   archived, and no pooled posterior is published.
   **Additional geometry challenge completed:** a
   [protocol committed before evaluation](COVERAGE_CHALLENGE_PLAN.md) tested
   the unchanged sampler on a correlated rotated box, an annular cylinder,
   and unequal disconnected balls. [Seven of nine learned runs passed](../analysis/calibration/proposal-coverage-challenges.md).
   Two rotated-box runs exceeded the maximum-weight limit; one also exceeded
   the regional-mass error limit, despite ESS above 200 in both. Every declared
   region was observed. All nine full-target oracle controls passed; all nine
   deliberately restricted controls passed usual screens but failed truth checks.
   The [interpretation](COVERAGE_CHALLENGE_RESULTS.md) separates these finite-budget
   failures from evidence of omitted support or a biological model defect.
   **Next deliverable:** develop and separately freeze a proposal revision
   addressing correlated targets, with a normalized, evaluable full density
   and a prospective evaluation that discloses these now-seen fixtures.
   Keep this completed failure study unchanged. In parallel, establish the
   model-to-assay mapping and missing independent replicate metadata below.
   **Complete when:** the revision and independent observable have fixed inputs,
   error criteria and published outcomes, including failures. Finite moment
   checks cannot establish complete within-region distributions or precise tails;
   agreeing biological runs cannot rule out a shared unseen region.

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

3. **Passive 2D and 3D immune measurements completed; the saturation explanation is corrected.**
   The [protocol committed before capture](IMMUNE_MEASUREMENT_PROTOCOL.md)
   observes canonical Control/RSL3/SDT scenarios with their existing geometry,
   oxygen configuration and treatment-specific seeds. The
   [replayable report](../analysis/immune-measurement-report.md) records all
   ferroptotic deaths, death-time LP, completed-release LP, horizon-censored
   deaths and terminal DAMP additions separately. It counts actual eligible
   living cells both as unique identities and repeated cell-step opportunities.
   Every event census reconciles, empty populations are explicit, and the
   full 24-condition production SHA is unchanged. SDT/RSL3 have 111/28 immune
   kills, 6,068/68,916 unique eligible cells, and mean completed-release LP
   19.69/17.98. These conditional populations differ; their ratios are not
   causal estimates of per-cell potency. Older sensitivity reports now correctly
   describe their final-count normalization without claiming it isolates earlier
   immune eligibility. Their numerical results are unchanged.
   The [separately frozen 2D protocol](IMMUNE_2D_MEASUREMENT_PROTOCOL.md) now
   reconstructs the historical 521/5 kills while preserving the full 33-condition
   summary. Its [report](../analysis/immune-2d-measurement-report.md) measures
   2,314,460/770,701 eligible SDT/RSL3 cell-step opportunities, mean activation
   0.038506/0.001309, and zero opportunities at or above half-maximal activation.
   The former deep-saturation explanation is unsupported in this realization.
   Mean completed-release LP is 19.123/18.266; three RSL3 releases are censored
   and recorded separately. See the [interpretation](IMMUNE_2D_MEASUREMENT_RESULTS.md).
   These observations do not isolate a causal geometry effect between 2D and 3D.
   **Next deliverable:** separately freeze a replicate study with a defined
   independent unit and observable, while resolving the assay mapping below.
   Neither single realization supplies experimental validation or a replacement
   P5 threshold. **Complete when:** the independent unit, inputs and criteria
   are fixed before evaluation, every declared outcome is published, and both
   existing archives remain reproducible and unchanged.

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
