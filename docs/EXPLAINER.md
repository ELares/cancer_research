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

## What is ferroptosis (the biology this started from)?

Cells normally protect themselves from a kind of rust-like chemical damage to their
membranes (called *lipid peroxidation*). **Ferroptosis** is a way of killing a cancer cell
by overwhelming those defenses until the damage runs away with itself and the cell dies.
It's interesting in cancer because some treatment-resistant "persister" cells — the ones
that survive chemo and later cause relapse — appear to be *unusually vulnerable* to it.

The project asks: if that's true, **how would you actually exploit it**, and which obstacles
(low oxygen deep in a tumor, drugs that can't physically reach the cells, the tumor's
supportive neighbor cells) would get in the way?

## What's actually in here?

1. **A map of the whole field.** Every cancer paper PubMed has indexed — 5,187,265 of them —
   sorted by treatment type, by where in the body the cancer is, and by what kind of study it
   was. Nothing here is a sample. That matters more than it sounds: a project that picks a few
   thousand papers to study can never tell "nobody has researched this" apart from "our search
   didn't find it", and this one does not have to guess.
2. **Labels somebody else wrote.** The sorting uses the National Library of Medicine's own
   subject headings and study-type labels, applied by professional indexers who have never
   heard of this project. So when we say a treatment has few clinical trials behind it, that
   is their label being counted, not ours being trusted.
3. **An honest account of what even that cannot see.** The census has one real blind spot and
   it is worth understanding: if the National Library of Medicine has no *name* for a
   treatment, no amount of searching will count it. Tumor Treating Fields is the clearest
   case — it is FDA-approved and has completed large trials, and it is simply invisible to
   this instrument. Throughout, such treatments are reported as *not measurable* rather than
   as zero, because a zero would read as "nobody is working on it," which would be false.
4. **A simulation engine.** A small, reusable program that models the tumor environment and
   nine different ways of attacking a cancer cell — ferroptosis chemistry (the deepest and
   the one it started with), light and ultrasound, radiation, checkpoint blockade, engineered
   T cells, a cancer-killing virus, physically destroying tissue, and antibodies carrying a
   drug — used to ask "if this idea were true, what would we expect to see?" The newer arms
   are much smaller than the ferroptosis one, and the project measures that gap rather than
   glossing it.

## The headline results (and the big caveat)

The simulations produced six directional ideas worth testing in a lab. The first three come from the ferroptosis and light/ultrasound work the project began with:

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
