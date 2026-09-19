# Independent assay acquisition candidate

Acquisition review, 19 September 2026. Lee et al., [*Nature Metabolism* (2024),
“Selenium reduction of ubiquinone via SQOR suppresses ferroptosis”](https://www.nature.com/articles/s42255-024-00974-4)
provides normalized source observations for a possible external endpoint check.
This review does not ingest a validation dataset, fit parameters, or establish
independent raw-replicate validation.

## Source and candidate endpoint

The [Extended Data Fig. 2 source workbook](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs42255-024-00974-4/MediaObjects/42255_2024_974_MOESM11_ESM.xlsx),
`42255_2024_974_MOESM11_ESM.xlsx`, was 20,728 bytes when inspected, with SHA256:

```text
27faf22389ea08063cda6e86d18c9c0c6ef126f83d427b41f64dd87802dee491
```

The best bounded candidate is **SK-Hep1 viability after 0.1 µM ML210 for
24 hours, Extended Data Fig. 2d**. The [primary manuscript's methods and
caption](https://pmc.ncbi.nlm.nih.gov/articles/PMC11694790/) identify Hoechst/PI
staining: viable cell count is total Hoechst-positive count minus PI-positive
count, normalized to vehicle. Panel 2e instead measures CellTiter-Glo after
6 µM erastin for 72 hours. These are different assay endpoints and durations.

| Endpoint | Workbook range | Three normalized observations | Mean | Sample SD |
|---|---|---|---:|---:|
| ML210, 0.1 µM, 24 h | `2d!F4:F6` | 0.713655, 0.642109, 0.676226 | 0.677330 | 0.035786 |
| ML210 + 10 µM ferrostatin, 24 h | `2d!G4:G6` | 1.062660, 1.081940, 1.093761 | 1.079454 | 0.015699 |
| RSL3, 0.1 µM, 24 h | `2d!D4:D6` | 0.055740, 0.039368, 0.043138 | 0.046082 | 0.008574 |
| Erastin, 6 µM, 72 h | `2e!D4:D6` | 0.085446, 0.074083, 0.073743 | 0.077757 | 0.006661 |

Means and sample SDs use the three supplied values, with SD denominator
`n − 1`, rounded to six decimals. They describe these observations; they are
not uncertainty estimates for model predictions. Rows `B4:B6` are numbered
1–3 in both sheets. Vehicle controls are `2d!C4:C6` and `2e!C4:C6`; panel 2d
also supplies RSL3 + ferrostatin (`E4:E6`) and ferrostatin alone (`H4:H6`).

## Unresolved acquisition requirements

- **Curve units:** ML210 dose columns in sheets 2g/2h and the
  [published figure axes](https://www.nature.com/articles/s42255-024-00974-4/figures/8)
  say “nmol”, not a concentration. The reviewed materials do not justify
  converting these values to nM. In 2g, ML210 doses occupy `B5:B13`.
- **Replicate identity:** the caption states biological `n=3` for 2g/2h, but
  ML210 has two columns per condition: for example, vehicle `2g!C5:D13` and
  selenite `2g!E5:F13`. FIN56 in 2g has three. Neither the workbook nor the
  [reporting summary](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs42255-024-00974-4/MediaObjects/42255_2024_974_MOESM2_ESM.pdf)
  resolves this conflict or supplies experiment/date/plate identities. The
  three numbered 2d/2e observations agree with the stated count, but do not
  establish which experimental units or paired controls produced each value.
- **Curve controls:** the lowest workbook dose is `1`; the figure labels its
  leftmost point `0`. Normalization and position suggest a plotting placeholder,
  but no explicit mapping identifies the underlying vehicle wells. Do not
  interpret that row as an established 1 nM exposure.
- **Underlying readings:** the workbook supplies normalized observations, not
  raw cell counts, luminescence, plate maps, or normalization denominators.
  The reporting summary puts raw data on author request. The reviewed
  [supplementary information](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs42255-024-00974-4/MediaObjects/42255_2024_974_MOESM1_ESM.pdf)
  does not reconcile these gaps. No authors were contacted.

## Requirement before a model comparison

A source-level reproduction of the single-dose summaries is possible now.
A predictive comparison still needs a prespecified mapping from the
[generic Glycolytic phenotype and `1 − death_rate`](../scripts/calibrate_kill_switch.py)
to SK-Hep1's vehicle-normalized viable-cell count, including proliferation and
assay effects. The [single-cell engine](../simulations/ferroptosis-core/src/biochem.rs)
uses a fixed 180-step horizon; its correspondence to 24 or 72 hours has not
been established here. Freeze parameters, endpoint mapping, and error criterion
before evaluation; disclose any shared cell lines with CTRPv2 and distinguish
experiment identity from compound holdout. The
[roadmap's raw-replicate validation requirement](RESEARCH_NEXT_STEPS.md) remains
pending even if a normalized endpoint comparison becomes feasible.
