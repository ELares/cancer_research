---
name: Simulation extension proposal
about: Propose a new biology layer, parameter, or analysis for the ferroptosis simulation
title: "[sim] "
labels: subject:simulation, enhancement
assignees: ''
---

Proposals to extend the model are welcome. Before writing one, it helps to know
the project's two standing rules for simulation layers, because following them
makes a proposal much easier to accept:

1. **Off by default, byte-identical.** A new layer must default to an inert state
   that leaves the production matrix output bit-for-bit unchanged (a regression
   test guards the SHA of `sim-tme-3d`'s `summary.json`). New mechanisms are
   opt-in.
2. **Direction over magnitude, honestly labeled.** Most layers are uncalibrated
   mechanistic scaffolding. The claim is the direction of an effect, not a fitted
   number, and the uncalibrated status is stated in
   `simulations/calibration/CALIBRATION_STATUS.md` and the module doc.

## The mechanism

What biology do you want to add or refine, and what is the expected direction of
its effect on ferroptosis (more sensitive, more resistant, or context-dependent)?

## Evidence

The key reference(s), ideally with PMIDs. If the direction is contested in the
literature, say so, and say which way the evidence leans.

## How it would fit

- Which module or part of the engine would it touch
  (`simulations/ferroptosis-core/src/...`)?
- What would the off-by-default identity case be?
- Is there a falsifiable prediction it would generate?

## What would calibrate it

The measurement or dataset that would move it from uncalibrated to anchored. If
that data does not exist publicly, that is worth noting too (it may become a
`help wanted` collaborator call).
