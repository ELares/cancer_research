# What this project is, in plain language

A one-page explainer for anyone — no biology or coding background needed. For the
technical version, see the [README](../README.md) and the
[manuscript](../article/drafts/v1.md).

---

## The one-sentence version

This is an open, free-forever attempt to use AI to (1) **map** thousands of cancer-therapy
research papers so we can see where the science is concentrated and where the "gaps" are
just artifacts of how we searched, and (2) **stress-test** specific biological ideas with
computer simulations — honestly labeling what's backed by real data and what's still just
a hypothesis.

## What is ferroptosis (the biology at the center)?

Cells normally protect themselves from a kind of rust-like chemical damage to their
membranes (called *lipid peroxidation*). **Ferroptosis** is a way of killing a cancer cell
by overwhelming those defenses until the damage runs away with itself and the cell dies.
It's interesting in cancer because some treatment-resistant "persister" cells — the ones
that survive chemo and later cause relapse — appear to be *unusually vulnerable* to it.

The project asks: if that's true, **how would you actually exploit it**, and which obstacles
(low oxygen deep in a tumor, drugs that can't physically reach the cells, the tumor's
supportive neighbor cells) would get in the way?

## What's actually in here?

1. **A literature map.** ~4,830 full-text cancer papers, auto-sorted by treatment type and
   cancer type. The useful twist: the project *measures its own blind spots* — it shows that
   some apparent "nobody has studied this" gaps are really just side effects of which journals
   are free to read and how the search was worded, not real holes in the science.
2. **A simulation engine.** A small, reusable program that models the ferroptosis chemistry
   and the tumor environment, used to ask "if this idea were true, what would we expect to
   see?"

## The three headline results (and the big caveat)

The simulations produced three directional ideas worth testing in a lab:

1. **Hitting two defenses at once works better than one.** Blocking two parallel "repair
   crews" (GPX4 and FSP1) at the same time is more than additive.
2. **Physical treatments (light/ultrasound) and chemical drugs hit different walls.** A drug
   has to physically diffuse deep into a tumor; light/ultrasound-delivered damage is limited
   by different things. So the *type* of obstacle a treatment faces depends on the treatment.
3. **Getting a drug deep into a tumor is brutally hard.** A drug that kills 40% of cells in a
   dish may reach only ~2% effectiveness behind the blood-brain barrier — before biology even
   fights back.

> **The big caveat, stated plainly:** these are **computer predictions, not medical advice
> and not validated cures.** The project is unusually honest about this — most of the
> simulation layers are explicitly labeled "we modeled the *direction* of an effect, but not
> a trustworthy number." Nothing here has been tested in a human, and several predictions are
> flagged as the project's *least* certain. The point is to generate good, falsifiable
> hypotheses for real scientists to test — not to claim a breakthrough.

## Why it's built this way

The author's view (see the [README](../README.md)) is that breakthroughs against diseases
that destroy families should be a shared human resource, not a product. Everything is
MIT-licensed and free to take, copy, and improve. The work is deliberately written so a
motivated student — not just a specialist — can follow it end to end and check it.

## How you can engage

- **Just curious?** Read the [README](../README.md), then browse the
  [analysis outputs](../analysis/). The corpus dashboard (`scripts/dashboard.py`) lets you
  explore the literature map interactively.
- **Have expertise?** Oncology, biochemistry, ferroptosis, immunology, computational biology
  — open an issue or a PR. The simulations especially benefit from people who can say "that
  assumption is wrong, here's the data."
- **Can run a wet-lab experiment?** The cheapest experiments that would confirm or *kill*
  these predictions are written up in [`PREREGISTRATION.md`](../PREREGISTRATION.md). Testing
  even one would be the single most valuable contribution.

> You don't need to be a cancer researcher. Curiosity and a willingness to look at the
> evidence — including the parts that say "we're not sure" — are enough.
