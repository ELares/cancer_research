# Collaborator brief: a registered, falsifiable, open ferroptosis program (#498)

**One page for a prospective domain co-author or external reviewer.** Seeded from
`analysis/contribution-plan-2026.md` and `PREREGISTRATION.md`.

## The ask

This is an open, MIT-licensed, reproducible cancer-research program (a
cross-literature corpus + an embeddable ferroptosis simulation engine) built so
far by a single non-domain author with AI assistance. It has **no domain
co-author, no wet-lab, and an empty contributor ledger**, so self-consistency is
its validation ceiling. The project's own scientific-review pass already caught a
real biology over-extension and a mis-cited diffusion coefficient that no
automated test caught, exactly the class of error an independent expert removes.
We are looking for one of:

- a **ferroptosis or redox biologist** (wet-lab) to run one of the cheap,
  pre-stated experiments below, or
- a **computational / systems oncologist** to stress-test the model and its
  calibration, or
- an **external reviewer** willing to read the manuscript and the calibration
  ledger critically.

What you get: **co-authorship on a registered, testable, openly licensed program**,
with a falsification-first design that makes a clean result publishable whether it
confirms or refutes.

## What the model already integrates

- A reproducible map of **4,830 full-text cancer-therapy articles** (19 mechanisms,
  22 cancer types, 803 journals) with honest coverage caveats (`analysis/`).
- An embeddable ferroptosis biochemistry engine (`ferroptosis-core`, MIT) with a
  per-layer calibration ledger that states plainly what is and is not anchored to
  data (`simulations/calibration/CALIBRATION_STATUS.md`, `MODEL_CARD.md`).
- Data-anchored legs: the in-vitro kill switch vs CTRPv2 GPX4 inhibitors (#330),
  System Xc-/erastin (#502), a joint multi-inducer in-vitro posterior (#500),
  spheroid zone geometry vs Browning 2021 (#333), tumor PK vs imidazole-ketone
  erastin (#334), Krogh penetration vs Primeau/Tannock (#335), and the ferroptotic
  trigger-wave speed vs Co 2024 (#482).

## The falsifiable predictions (each has a numeric model output + threshold)

P1 GPX4i + FSP1i synergy in FSP1-low persisters (Bliss ~1.99x; refuted if CI > 0.8).
P2 physical-ROS modalities less depth-limited than RSL3 in large spheroids.
P3 a days-timescale post-withdrawal vulnerability window with sequential defense recovery.
P4 SDT retains more efficacy than RSL3 under hypoxia (the contested keystone, direction only).
P5 dense ferroptotic kill is more immunogenic per cell than sparse kill (~4:1 in 3D).
P6 CAFs protect RSL3 more than SDT.
P7 RSL3 efficacy drops at acidic pH via ion trapping.
P8 a persister-targeting inducer has the opposite size-dependence to generic cytotoxics.

Full statements with quantitative model outputs and pre-stated numeric
falsification thresholds are in `PREREGISTRATION.md`.

## The experiment menu (cheapest first)

- **E4 / P1 (cheapest, highest leverage):** a GPX4i + FSP1i dose matrix in a
  persister-enriched line. Detailed protocol: `analysis/p1-wetlab-protocol.md`.
- **E1 / P4 (the keystone):** spheroid RSL3 vs SDT kill at measured hypoxia
  (pimonidazole or hypoxia chamber, confocal depth sectioning).
- **E2 / P6:** CAF-coculture IC50 shift, RSL3 vs SDT.
- **E3 / P5:** spheroid-supernatant DAMP + DC-maturation assay.
- **E5 / P3:** sequential defense recovery after drug withdrawal (time course).
- **E6 / P7:** RSL3 efficacy and intracellular concentration vs pH.

Each is a single, isolatable experiment with a pre-stated falsifying outcome
(`PREREGISTRATION.md` Part 2).

## Who would be a good fit (selection criteria, not a named list)

We are deliberately not naming individuals here. A strong collaborator profile:

- publishes on **ferroptosis, lipid peroxidation, GPX4/FSP1/System Xc-, or
  drug-tolerant persisters** (wet-lab), or on **tumor-microenvironment modeling /
  computational oncology** (modeling), and
- values **open, reproducible, falsification-first** work (the program is MIT
  licensed and registered), and
- is willing to either run one cheap pre-stated experiment or critically review
  the model and its calibration.

Identifying specific candidates (recent ferroptosis-methods authors, lab leads
running the relevant assays, computational-oncology groups) is a human step left to
the author, who can apply these criteria to the corpus and the current literature.

## How to engage

Open a GitHub issue or pull request on the repository, or contact the author. The
predictions are registered (`PREREGISTRATION.md`); a collaborator joins a program
where the success criteria were fixed in advance, so a clean experiment is
publishable either way.
