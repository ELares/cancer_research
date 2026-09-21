# Passive observation of the existing 3D immune comparison

This protocol adds measurements to an existing deterministic simulation. It
does not introduce a mechanism, tune parameters, replace a historical result,
or test an experimental efficacy claim. The implementation and this protocol
are committed before the production capture; the archive names that commit.

## Scope fixed before capture

Run the existing `immune_Control`, `immune_RSL3` and `immune_SDT` conditions
from `generate_conditions()`. Preserve their names, since their FNV-1a hashes
determine treatment-specific runtime seeds. All use geometry seed 42, a 60³
grid, 20 µm cells, tumor radius 540 µm, an oxygen gradient with λ=120 µm,
180 steps, constant treatment, default biochemical parameters and
`SpatialImmuneConfig::for_3d()`. These are matched settings, not shared
treatment RNG streams. No stromal shielding, pH modulation, regrowth or
additional death modality is enabled. Reject configuration overrides.

The sidecar records the complete biochemical, spatial and immune parameter
configuration, seeds and RNG salts. The default immune configuration has
DAMP-per-LP 1, diffusion fraction 0.025, clearance 0.03, activation Kd 50,
kill rate 0.02, PD-1 brake 0.7, and zero anti-PD-1 efficacy, exhaustion and
ferroptotic immunosuppression. These are model assumptions, not independently
calibrated experimental values.

## Observation points and denominators

1. **Ferroptotic death:** after the biochemical update, record lattice identity,
   death step and threshold-crossing LP. Lattice identity denotes a unique
   cell only because these conditions have no repopulation.
2. **Completed release:** at `death_step + 5`, record LP and DAMP actually
   injected before that step's diffusion/clearance and immune update. Preserve
   death-time LP for that same completed-release cohort.
3. **Immune opportunities:** after diffusion/clearance, immediately before
   killing, count living tumor cells admitted by the existing step ≥60 and
   local DAMP ≥0.01 checks. Include cells killed by that update. Archive the
   unique cell identities, first/last eligible steps, number of opportunities,
   sum of local DAMP over opportunities, per-step counts and actual kill events.
4. **Horizon censoring:** a death at step 174 releases at the final simulated
   step, 179; a death at 175 is censored because its release is scheduled at
   180. Record horizon LP separately from completed release LP. The simulator's
   terminal flush injects DAMP from censored deaths after the last immune
   update. This addition cannot explain any observed immune kill.

Report all-death LP, completed-cohort death and release LP, and censored-cohort
horizon LP with their respective counts. Empty means are undefined, never
zero. Report total immune kills, kills divided by unique eligible cells, and
kills divided by eligible cell-step opportunities. Neither denominator removes
competing-risk or treatment-selection effects. Cell-step rates are descriptive
model summaries, not fitted hazards or independent-cell efficacy estimates.

The final non-ferroptotic pool used in older sensitivity/uncertainty analyses
is a different population. Preserve that historical statistic and numerical
results while correcting claims that it isolates per-cell amplification.

## Reconciliation and passivity checks

- Ferroptotic deaths = completed releases + horizon-censored deaths.
- Final dead tumor cells = ferroptotic death events + immune kill events.
- Immune kills ≤ unique eligible cells ≤ eligible cell-step opportunities.
- Every per-step event count and DAMP sum agrees with the event ledger.
- Per-cell and per-step opportunity/DAMP totals agree; there are no eligible
  opportunities after a cell dies or before immune activation.
- Terminal DAMP addition equals the sum from censored deaths and the increase
  in final field mass. Cumulative injections need not equal field mass because
  diffusion and clearance occur between observations.
- Observing a run consumes no RNG and receives no mutable simulation state.
  Rust tests compare serialized legacy results and per-step snapshot arrays
  with observation on/off and compare observation across Rayon thread counts.
- The entire default 24-condition `summary.json` must retain its existing SHA:
  `973d9443bd05851721d53d9a15e53dc828517bba4bf307eb1c4f0b1f60d27289`.
  Each observed legacy result must also equal the corresponding matrix row.
  Never update the expected SHA to accommodate passive instrumentation.

## Capture and offline reconstruction

Commit the implementation, tests and protocol first. With Rust 1.96.0 (pinned
by the workspace), from a clean checkout and with no `FERRO_*` environment
overrides:

```bash
python3 scripts/immune_measurement_report.py --capture --threads 8
```

The command builds with `--locked`, runs the unchanged default matrix and the
three observed conditions in separate temporary directories, validates their
results, and creates `analysis/immune-measurements/`. It refuses to replace an
existing archive. For a new independent capture, supply a new `--archive`
directory and `--report` path. `CARGO_TARGET_DIR` may select a build cache.

The archive contains the production summary, compressed raw observations,
a compressed source bundle and a manifest with source commit, individual
source/artifact hashes, binary hash, toolchain, platform and thread count.
The source bundle preserves historical reconstruction even after the current
simulator changes; no historical Git objects or extraction of archive paths
are needed to verify it.

```bash
# Reconcile every event and regenerate the committed Markdown, without Rust:
python3 scripts/immune_measurement_report.py

# Focused validation, followed by the normal repository checks:
python3 -m pytest tests/test_immune_measurement_report.py -q
cd simulations
cargo test -p sim-tme-3d immune_measurements
```

## Interpretation fixed before seeing the observations

Publish all three conditions and every population, including empty populations
and terminal additions. There is no superiority gate and no parameter search.
One fixed geometry and one existing runtime seed per condition do not support
simulation-replicate confidence intervals. Model LP/DAMP units and cell-steps
have no new experimental calibration. The 3D observations cannot retrospectively
explain the historical 2D 104:1 total-kill ratio, establish DAMP potency per
dead cell, or replace a prespecified independent experimental threshold.

The next scientific step is to measure corresponding populations in the 2D
headline scenario and/or a frozen replicate study, then compare a clearly
mapped observable with independently measured assay data. This measurement
report supplies missing denominators; it does not close the validation gap.
