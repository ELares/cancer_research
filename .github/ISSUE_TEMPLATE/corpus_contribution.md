---
name: Corpus or literature contribution
about: Suggest a paper, dataset, or correction to the literature corpus and tagging
title: "[corpus] "
labels: subject:taxonomy-and-search
assignees: ''
---

You do not need to write any code to contribute here. Pointing us at a paper we
are missing, a mis-tagged record, or a dataset that would calibrate a model layer
is one of the most valuable things you can do. The corpus is a frozen snapshot for
reproducibility, so most additions land in the living-review layer or a future
re-freeze rather than changing the published numbers, and that is fine.

## What kind of contribution is this

- [ ] A paper the corpus is missing (especially a landmark that would change a
      field-level claim)
- [ ] A mis-tagged or mis-classified record (wrong mechanism, evidence tier,
      cancer type, or tissue)
- [ ] A dataset that would calibrate a model layer (see
      `simulations/calibration/CALIBRATION_STATUS.md` for what is data-gated)
- [ ] A taxonomy or search-query suggestion
- [ ] Something else

## Details

For a paper: the PMID, DOI, title, and one line on why it matters (what claim it
supports, contradicts, or fills a gap in).

For a mis-tagging: the PMID and what the tag should be instead, with a reason.

For a dataset: what it measures, where it lives (a public URL if possible), its
license, and which prediction or model layer it would test.

## Honesty note

If this contradicts something the repo currently claims, say so plainly. The
project optimizes for being right for patients, not for defending its current
framing, so a contribution that refutes a claim is as welcome as one that supports
it.
