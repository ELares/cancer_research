# Passive observation of the canonical 2D immune headline

This study measures the existing `sim-tme` run behind the manuscript's immune
headline. It adds event accounting and observations of the DAMP field during
actual opportunities for immune killing. It does not change the biochemical,
transport or immune rules, fit parameters, or validate an experimental outcome.
The [completed 3D study](../analysis/immune-measurement-report.md) remains a
separate frozen artifact.

## Prior knowledge and capture chronology

The historical baseline counts are already known: SDT produces 521 immune
kills and RSL3 produces 5, a ratio of 104.2. The existing
[20-seed study](../analysis/seed-replication-report.md) is also known: seeds
42–61 gave a median defined per-seed ratio of about 107, a range of about
49–260, and one seed with zero RSL3 immune kills. Those are count summaries;
they do not contain the death/release cohorts or eligible-cell exposure
records defined here. The 3D study's results were seen before this design.
This is a declared measurement extension of known scenarios, not an unseen
test of the 104:1 hypothesis or a new replicate study.

The order of work is fixed:

1. Build the uninstrumented `sim-tme` from a recorded source commit with the
   pinned Rust toolchain. Save its complete default 33-condition
   `tme_summary.json`, source identity and binary identity. Establish
   `simulations/sim-tme/expected_summary.sha256` from those bytes before
   evaluating the instrumented production run.
2. Implement and test the passive observer and offline accounting. Commit
   this protocol, implementation, configuration contract, tests and unchanged
   baseline hash before capturing the new production observations. Small
   synthetic and reduced-grid tests check software behavior; they are not
   evidence for the biological interpretation.
3. Build the committed instrumented binary. Its complete default 33-condition
   summary must match the uninstrumented summary byte for byte. The three
   observed legacy results must match their exact canonical matrix rows.
   A mismatch stops capture; never change the expected hash to accept it.
4. Capture the three declared arms without tuning or selecting among their
   outcomes. Publish the observations, provenance, reconciliation and all
   declared summaries, including empty populations. The archive records the
   implementation freeze and the earlier baseline provenance separately.

The pre-instrumentation baseline was obtained from source commit
`e9ec92c022ede1b00258789fecdcb9a46d9481d8`. Its complete 33-condition summary
has SHA-256
`e04f9699ddca9b84145cdecbcd8e21abb0bc760dc3d909908358d4554a5eee98`.
The [frozen golden](../simulations/sim-tme/expected_summary.sha256) also records
the original binary identity. This baseline check precedes the new production
event and activation-window measurements.

## Inputs and populations fixed before capture

Select Control, RSL3 and SDT from the baseline immune block: oxygen condition
`gradient_120um`, `immune_on`, stromal protection off, and no pH modulation.
The later stromal and pH blocks contain additional immune-on rows and are not
interchangeable with this block. Anti-PD-1 conditions are outside this study.
The sidecar identifies the arms as `immune_Control`, `immune_RSL3` and
`immune_SDT`, in that order.

Preserve the 500 × 500 grid, 20 µm cell spacing, circular tumor of radius
225 cells (4,500 µm), geometry seed 42, oxygen penetration length 120 µm and
180 steps numbered 0–179. Use `Params::default()`, the existing 2D spatial
parameters, and `SpatialImmuneConfig::for_2d()`. The existing oxygen-independent
SDT setting remains unchanged; no oxygen-dependence override is allowed.
No regrowth or additional cause of death is introduced.

Preserve the original treatment-specific runtime seeds: Control 42, RSL3
10,000,042 and SDT 20,000,042. Preserve initialization order, per-cell seed
offsets and step strides. Shared geometry does not mean shared treatment RNG
streams. Reject `FERRO_SEED`, biochemical parameter overrides and extra mode
arguments; capture requires all `FERRO_*` environment overrides to be absent.

The 2D immune defaults are DAMP-per-LP 1, diffusion fraction 0.08 across up to
eight neighbors, clearance 0.03 per step, activation Kd 50, immune kill rate
0.02 and PD-1 brake 0.7. Anti-PD-1 efficacy, exhaustion and ferroptotic
immunosuppression are zero. These are model assumptions. Record the complete
biochemical, spatial and immune configuration, geometry, seeds, RNG offsets
and strides in the sidecar and validate them against a frozen 2D configuration
contract. Keep the 3D configuration contract unchanged.

The standalone Rust mode is `--immune-measurements` and writes
`output/tme/immune_measurements.json`, separate from the legacy summary. The
sidecar declares `schema_version: 2`, `simulator: "sim-tme"` and `dimension: 2`.
It carries the unchanged legacy result for each arm and a separate
`final_damp` summary with total and peak field values. The archive driver
performs the production parity checks; invoking the Rust mode alone does not
establish old/new matrix parity.

## Event phases and denominators

1. **Ferroptotic death:** after biochemistry, record lattice identity, death
   step and LP at threshold crossing. A lattice identity denotes one cell
   because these arms have no repopulation.
2. **Completed release:** at `death_step + 5`, record the stored grace-end LP
   and DAMP actually injected before diffusion, clearance and immune killing.
   Preserve death-time LP for this same completed-release cohort.
3. **Immune opportunity:** after diffusion and clearance, immediately before
   the immune update, include each living tumor cell at step ≥60 whose local
   DAMP is ≥0.01. Include a cell killed by that update. Record per-cell first
   and last eligible steps, number of opportunities and summed local DAMP;
   also record per-step counts and sums and every actual immune-kill event
   with its cell identity, step and local DAMP. Unique eligible cells and
   repeated cell-step opportunities have separate denominators.
4. **Horizon censoring:** a death at step 174 releases at the final simulated
   step 179. A death at 175 has release scheduled at 180 and is censored, as
   are later deaths. Record their horizon LP and terminal DAMP additions
   separately. The existing terminal flush runs after the last immune update;
   those additions cannot explain kills already observed.

Report event counts, all-death mean LP, completed-cohort mean death and release
LP, and censored-cohort mean horizon LP. Report completed and terminal DAMP
injections separately from field mass before and after the terminal flush.
Report immune kills, kills per unique eligible cell and kills per eligible
cell-step. Name the population beside every conditional mean or fraction.
Empty means, maxima and fractions are undefined; empty counts and sums are
zero. Do not substitute a censored horizon value for a completed release.

The historical fraction `immune_kills / max(total_tumor - ferroptosis_kills, 1)`
uses final non-ferroptotic counts, a different population. Neither that
normalization nor the new opportunity denominators remove competing death
or treatment-selection effects.

## Activation-window measurements fixed before capture

The manuscript attributes the large 2D count ratio partly to DAMP saturation.
Observe the field where an immune draw actually occurs, rather than infer
activation from a final heatmap or an average DAMP concentration. For each
eligible cell-step with local DAMP `D`, use the existing activation expression
`a(D) = D / (D + Kd)`, with `Kd = 50`.

Accumulate the following both per eligible cell and per simulation step:

- `activation_sum`: the sum of `a(D)` over the eligible opportunities.
- `max_local_damp`: maximum local DAMP over those opportunities.
- `damp_ge_kd_opportunities`: number with `D >= Kd` (activation at least 0.5).
- `damp_ge_9kd_opportunities`: number with `D >= 9 * Kd` (activation at least 0.9).

The thresholds come from the existing response function and are fixed before
production. The report gives the total activation sum divided by total
opportunities, the maximum eligible local DAMP, and each threshold count and
its share of opportunities. Sum over opportunities directly; do not average
per-cell means equally or apply `a()` to mean DAMP. Preserve per-cell and
per-step records so their totals, maxima and threshold counts can be checked.
Report all three arms even when no opportunity reaches either threshold.

These are descriptive exposures at the model's existing immune windows. They
can show whether eligible cells reached the specified activation levels in
this run. They do not identify how much of the treatment difference is caused
by saturation, spatial density, release per death, transport or target-cell
depletion. The 2D and 3D runs also differ in size, neighborhood, diffusion and
RNG construction; their ratio difference is not an isolated dilution effect.

## Reconciliation and scientific passivity

- Ferroptotic deaths equal completed releases plus horizon-censored deaths;
  final dead tumor cells equal ferroptotic events plus immune-kill events.
  The causes do not overlap, and immune deaths produce no ferroptotic release.
- Every recorded identity belongs to the canonical tumor circle. Immune kills
  cannot exceed unique eligible cells, which cannot exceed either the tumor
  population or the number of opportunities.
- Per-step deaths, completed releases and kills agree with their event records.
  Per-cell and per-step opportunity totals agree. Their first/last-step bounds
  and death times exclude opportunities before activation or after death.
- Per-cell and per-step DAMP and activation sums agree with relative tolerance
  `1e-10` and absolute tolerance `1e-8`. Maxima agree, and threshold counts
  reconcile
  exactly: high-threshold count ≤ low-threshold count ≤ opportunity count.
- Kill-event DAMP fits within the corresponding opportunity sums. Eligible
  DAMP is a subset of the whole field mass. From a zero initial field,
  completed releases and per-step clearance reproduce pre-terminal mass;
  the stable 2D diffusion rule conserves mass, including at grid boundaries.
- Terminal additions equal the sum from censored events and the increase in
  final field mass. Cumulative injections need not equal retained field mass.
- The observer consumes no RNG and receives no mutable simulation state.
  Tests compare observed and unobserved legacy results and cell/field
  trajectories, including the activation boundary and release/horizon boundary.
  The complete old/new production summary parity is a separate required check.

## Separate capture and offline reconstruction

Use the pinned Rust 1.96.0 toolchain and a clean committed checkout. Execution
uses the existing serial 2D loops; it does not use the 3D Rayon execution mode.
The dimension-specific driver is `scripts/immune_2d_measurement_report.py`; the
archive is `analysis/immune-2d-measurements/` and the report is
`analysis/immune-2d-measurement-report.md`. Capture must select the executable
reported by the pinned Cargo build and run a private copy in temporary output
directories. It must not replace the existing 3D artifacts or an existing 2D
archive. Reports must remain outside archives and cannot alias archive files.

```bash
# After the implementation and protocol freeze, capture the declared study:
python3 scripts/immune_2d_measurement_report.py --capture

# Reconcile the committed archive and regenerate its report without Rust:
python3 scripts/immune_2d_measurement_report.py
```

After an archive exists, an independent capture requires a new `--archive`
directory and an external `--report` path, for example:

```bash
python3 scripts/immune_2d_measurement_report.py --capture \
  --archive analysis/immune-2d-measurements-replicate \
  --report analysis/immune-2d-measurements-replicate.md
```

Archive `observations.json.gz`, `baseline-summary.json`, `sources.tar.gz` and
`manifest.json`. Record source and artifact hashes, the
baseline and instrumented binary identities, toolchain, platform and execution
settings. Verification uses the frozen baseline and sources, not the current
checkout's evolving simulator. Preserve the original 2D baseline and the
existing 3D archive when producing any later independent capture.

## Interpretation fixed before the new production observations

Publish every declared arm and measurement. There is no efficacy, superiority
or saturation acceptance gate, no parameter search and no outcome-dependent
selection of populations or thresholds. Technical reconciliation and passivity
are requirements for publishing a valid capture, not biological success tests.

This single canonical realization supplies no new replicate confidence
intervals. Cells and repeated opportunities are not independent simulation
replicates. LP and DAMP remain model units; cell-steps have no newly calibrated
conversion to hours. Activation is a model function, not measured dendritic
cell maturation. The earlier 20-seed count study remains separate and does
not provide uncertainty for these newly measured event populations.

The historical P5 text and its within-1.5x experimental falsification threshold
in [PREREGISTRATION.md](../PREREGISTRATION.md) remain unchanged. A measured
model release ratio neither establishes an assay-scale immunogenicity ratio
nor replaces that threshold. Independent assay mapping and raw-replicate
validation remain pending after this measurement study.
