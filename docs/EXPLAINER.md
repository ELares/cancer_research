# What this project is, in plain language

A one-page explainer for anyone — no biology or coding background needed. For the
technical version, see the [README](../README.md) and the
[manuscript](../article/drafts/v1.md).

---

## The one-sentence version

This open project uses AI, literature analysis, and computer simulations to **map cancer
research and test biological hypotheses**, while making its evidence, assumptions, and
failures available for others to inspect and improve.

## What is ferroptosis (the biology this started from)?

Cells protect their membranes from chemical damage called *lipid peroxidation*.
**Ferroptosis** is a form of cell death in which those defenses are overwhelmed.
The project's earliest models explored this process and the obstacles a treatment
might encounter, such as low oxygen, limited drug delivery, and neighboring cells.
That remains its most developed modeling area; the literature collection is broader.

## What's actually in here?

1. **A large literature collection.** The recorded collection contains **5,187,265
   articles**: **4,403,994** selected through National Library of Medicine subject headings
   for neoplasms and adjacent topics, plus **783,271** recovered through text matching
   before subject indexing. These are defined collection streams, not a guarantee that
   every relevant paper is included. See the [scope and denominators](../analysis/europepmc-access-ceiling.md#which-denominator-these-percentages-use).
2. **Census charts and a smaller searchable archive.** The committed census charts
   describe the **4,403,994 indexed records**. Browser article search covers a separate,
   historical **4,830-record archive** assembled through keyword searches. Searching that
   archive does not search the larger collection.
3. **Measured limits to the map.** Census analyses count subject headings and publication
   types supplied by the National Library of Medicine. The project's choice of headings,
   missing labels, and unequal full-text availability still affect what can be measured.
   An unmeasurable treatment is not a zero, and publication counts do not measure treatment
   effectiveness. The [mechanism profiles](../analysis/census-mechanism-profile.md) explain
   why raw volume cannot fairly rank different mechanisms.
4. **A simulation engine.** The suite represents **ten treatment arms plus an untreated
   control**, with unequal depth and calibration. It asks how a proposed mechanism behaves
   under stated assumptions. The [model card](../MODEL_CARD.md) describes where evidence
   anchors the models and where their behavior remains exploratory.

## What the models are testing

Examples include whether blocking two cellular defenses together can outperform either
alone, how oxygen changes simulated treatment responses, and how tissue barriers affect
drug delivery. These are questions to test under specified conditions. The
[uncertainty report](../analysis/headline-uncertainty-report.md) shows why individual
simulation percentages should not be read as dependable treatment effects.

## Where the work stands

**The latest sampling study passed its computational checks.** Three independent runs
found enough suitably weighted parameter combinations that fit a fixed set of biological
measurements, and their summaries agreed within prespecified limits. These are repeated
computer runs, not repeated laboratory experiments. See the [sampling results](../analysis/calibration/joint-resample-sampling.md).

**A harder synthetic challenge exposed a limitation.** Seven of nine runs passed on
artificial targets whose answers are known. Two runs on a correlated target failed because
too much estimated weight depended on a few samples; one also misestimated a region's
share. The [challenge therefore failed its overall criterion](COVERAGE_CHALLENGE_RESULTS.md).
Even agreeing runs cannot prove that every important region has been explored.

**Independent assay validation is still pending.** An existing compound holdout shares
cell lines and a screening experiment with the fitting data. An
[external assay candidate](INDEPENDENT_ASSAY_CANDIDATE.md) has been reviewed, but raw
readings, experiment identities, and the mapping between the assay and model remain
unresolved. Computational progress does not establish biological or clinical validity.

## Why it's built this way

The author's view (see the [README](../README.md)) is that progress against diseases that
destroy families should be shared. The project's code and original work are MIT-licensed;
third-party materials retain their [source licenses](../PROVENANCE.yaml). The aim is to
make the work accessible enough for a motivated student to follow and challenge it.

## How you can engage

- **Just curious?** Explore the [project site](https://elares.github.io/cancer_research/),
  then follow a result to its [analysis report](../analysis/).
- **Have expertise?** Oncology, biochemistry, ferroptosis, immunology, computational biology
  — open an issue or a PR. The simulations especially benefit from people who can say "that
  assumption is wrong, here's the data."
- **Can help with experiments or data?** The [research roadmap](RESEARCH_NEXT_STEPS.md)
  describes the missing validation inputs, and [preregistered predictions](../PREREGISTRATION.md)
  set out hypotheses that others can challenge.

> You don't need to be a cancer researcher. Curiosity and a willingness to look at the
> evidence — including the parts that say "we're not sure" — are enough.
