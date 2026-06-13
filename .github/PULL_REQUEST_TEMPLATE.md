<!--
Thanks for contributing. This template is a guide, not a gate. Delete the parts
that do not apply. A small, clear PR with a couple of these boxes ticked is very
welcome.
-->

## What this changes

A short description of the change and why.

Closes #<!-- issue number, if this fully resolves one. Use "refs #N" for partial work. -->

## Checklist

- [ ] Tests pass locally (`python3 -m pytest tests/ -q` and, for Rust changes,
      `cd simulations && cargo test --workspace`).
- [ ] For Rust changes, `cargo fmt --all --check` is clean.

### If this touches the simulation

- [ ] The change is **off by default** and the production matrix stays
      **byte-identical** (the `sim-tme-3d` `summary.json` SHA regression test still
      passes). New mechanisms are opt-in.
- [ ] Uncalibrated parameters are labeled as such (in the module doc and
      `simulations/calibration/CALIBRATION_STATUS.md`); the claim is a direction,
      not a fitted magnitude.

### If this touches claims, the manuscript, or the corpus

- [ ] Strong claims are traceable to the corpus, an analysis output, or a verified
      external source, and any taxonomy artifact or corpus gap is flagged rather
      than buried.
- [ ] Citations are real and verifiable (PMID or DOI), not asserted from memory.
- [ ] The frozen corpus (`corpus/INDEX.jsonl`) and the manuscript's quantitative
      numbers are unchanged unless that is the explicit point of the PR.

## Notes for reviewers

Anything that needs context: a design tradeoff, a deferred follow-up, or a part
you are unsure about. Honest uncertainty is welcome here.
