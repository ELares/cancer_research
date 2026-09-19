# Calibration targets within recorded dose support

The September 2026 review found that the erastin calibration evaluated every
CTRPv2 fitted curve at 100 µM, although the largest recorded maximum dose was
66 µM. The GPX4-inhibitor grids also evaluated a few individual curves outside
their recorded ranges. A fitted curve can return a number there; that does not
make the number an observation supported by the screen.

## Target construction

[`ctrp_dose_support.py`](../scripts/ctrp_dose_support.py) selects one fixed
cohort per compound and dose grid. A retained curve's recorded minimum and
maximum must contain the entire grid. Selection therefore cannot change the
cell-line denominator at each dose. Invalid units, ranges, coefficients, empty
cohorts, and unsupported grids fail before fitting or writing results.

| Use | Compound | Dose grid range (µM) | Retained / available curves |
|---|---|---:|---:|
| GPX4 training | ML162 | 0.01–10 | 793 / 795 |
| GPX4 held-out compound | ML210 | 0.01–10 | 759 / 762 |
| System Xc− training | Erastin | 0.1–30 | 785 / 795 |
| GPX4 cross-mechanism contrast | Erastin | 0.01–10 | 790 / 795 |

The targets remain pointwise medians of published fitted curves within these
ranges. They are not raw replicate measurements at every requested dose, nor
the curve of one particular median cell line. Excluded lines may differ from
retained lines; the estimates apply to the supported cohorts. Each result JSON
records the CSV hash, requested grid, per-dose support, and excluded row IDs.
The loader hashes the same bytes it parses before starting the fit.

## What changed in the point fits

| Metric | Previous target construction | Supported-cohort refit |
|---|---:|---:|
| ML162 fit RMSE | 0.0504 | 0.0513 |
| ML210 held-out-compound RMSE | 0.0689 | 0.0684 |
| Erastin fit RMSE | 0.1003 | 0.0697 |
| Erastin shared-switch RMSE | 0.1716 | 0.1101 |

The best point-fit parameter vectors are unchanged. The erastin RMSE now uses
six supported doses instead of seven doses that included an extrapolation.
These errors use different targets: the lower value is not evidence that the
model became a better biological predictor. The previous values are archived
in [the pre-correction commit](https://github.com/ELares/cancer_research/tree/7bbc2bf564cd6e7f6de86fbbd168f4dcf8cdffad/analysis/calibration).

ML210 is held out by compound within the same screen. Cell lines overlap with
the training compound. The ABC bands describe variation among accepted model
parameter draws with a fixed simulation seed; they do not include experimental
replicate variation or validate predictions on new cell lines or another assay.
The joint acceptance rule keeps its reference vector, priors, tolerance factor,
and minimum accepted count; target correction does not justify widening the
tolerance to obtain a desired result.

## Joint inference is currently underpowered

The corrected 40,000-draw joint run accepted **4 draws**, below the unchanged
minimum of **20**. Its reference distance is 0.1586 and its acceptance threshold
is 0.1745 (factor 1.10). The previous joint credible intervals and held-out
coverage claim are superseded. The four-draw quantiles and coverage remain in
JSON for audit but are not presented as usable posterior inference. Reports
and the plot withhold those intervals; downstream information and
identifiability reports mark the joint inference unassessable.

This is a failed sampling requirement, not proof that the biological model is
wrong or correct. A follow-up should assess more efficient sampling against the
same fixed targets and criterion, with independent repeated runs to check
stability. Independent assay validation is still required. The single-inducer
ABC remains a separate, provisional result; it does not fill the joint gap.

## Reproduction

Build the Python extension from the checked-out simulation source using the
repository's pinned toolchain, then run:

```bash
python scripts/calibrate_kill_switch.py
python scripts/calibrate_erastin.py
python scripts/abc_posterior.py --n-draws 2000
python scripts/abc_joint_posterior.py --n-draws 40000
python scripts/abc_posterior_information.py
python scripts/abc_acceptance_diagnostic.py --draws 20000
python scripts/identifiability_report.py
```

The point fits use 4,000 simulated cells per dose; ABC uses 2,000. The cell
simulation seed remains 42 and the ABC prior-draw seed remains 12345. This
review rebuilt the extension from simulation sources at `7bbc2bf5`; the change
does not modify Rust simulation code or production defaults. See the generated
[point-fit report](../analysis/calibration/kill-switch-calibration.md),
[erastin report](../analysis/calibration/erastin-calibration.md), and
[joint report](../analysis/calibration/joint-posterior.md) for the actual results.

The existing [headline transfer experiment](../analysis/headline-at-fitted-cascade.md)
records older parameter vectors explicitly. It remains a historical sensitivity
experiment; it has not tested the newly refitted posterior's admissibility.
The next validation work is described in the [research roadmap](RESEARCH_NEXT_STEPS.md).
