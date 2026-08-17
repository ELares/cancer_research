#!/usr/bin/env python3
"""Where does the thesis sit among ALL modalities, not the three it names? (#725)

WHY
---
`atlas_thesis_position.py` measures where this project's thesis sits by
intersecting MeSH `Ferroptosis` with a hand-written list of six modalities. Its
honest headline -- the sonodynamic leg rests on roughly thirty papers in the
entire indexed cancer literature -- is one of the best results in the repo.

But a hand-written list cannot surface a leg nobody thought to look for. The
table enumerates the modalities the thesis already believes in, so a large
ferroptosis intersection with something outside that list is structurally
invisible.

This ranks the intersection against every modality in the committed
`analysis/modality-partitions.json` -- a descriptor universe built for #724 by
five independent panels and reviewed for symmetry -- rather than a list written
by the same person making the claim.

WHAT THE ANSWER IS FOR
----------------------
Two things, and the second is why it is worth doing.

  RANK. Where do the thesis's own legs sit among all modalities that intersect
  ferroptosis? A leg being small is not news -- the repo already says so. A leg
  being small while something unexamined is large is news.

  PRECEDENT. If a modality the thesis does not discuss has a much larger
  ferroptosis literature, that is either strong support for the mechanism or a
  novelty claim that needs qualifying, and the manuscript should say which.

A WITHDRAWN CLAIM THIS REPLACES. An earlier version of #725 asserted that
ferroptosis x radiotherapy is "roughly an order of magnitude larger than the
sonodynamic leg". That reproduced under no descriptor set: the bare
`Radiotherapy` descriptor gives a count BELOW the sonodynamic leg, and the
apparent 10x came from comparing a wide text-stem count against a narrow
descriptor count -- asymmetric rules, the same error corrected in #722. This
script compares descriptor to descriptor throughout, and reports the whole
ranking so no single pair can be cherry-picked.

Usage:
    python scripts/atlas_thesis_rank.py
"""

import argparse
import gzip
import json
from collections import Counter
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
ATLAS = PROJECT_ROOT / "corpus" / "atlas"
PARTITIONS = PROJECT_ROOT / "analysis" / "modality-partitions.json"
OUT_MD = PROJECT_ROOT / "analysis" / "atlas-thesis-rank.md"
OUT_JSON = PROJECT_ROOT / "analysis" / "atlas-thesis-rank.json"

FERROPTOSIS = "ferroptosis"

# The legs the thesis names, so their rank can be located rather than assumed.
THESIS_LEGS = {
    "sonodynamic therapy": "ultrasonic therapy",
    "photodynamic therapy": "photochemotherapy",
    "focused ultrasound": "high-intensity focused ultrasound ablation",
    "drug resistance": "drug resistance, neoplasm",
}

# Modalities with no usable MeSH descriptor. Named rather than silently absent,
# which is the treatment `mechanism_recall.py` already gives unmeasurable
# mechanisms -- reporting them as a zero would be a different claim.
# Modalities claimed to have no usable descriptor, WITH the candidate each
# claim rests on so the claim can be checked rather than trusted. A
# hand-written "no descriptor" list beside a census that can answer the
# question is the same shape as a hand-written list beside an exhaustive
# match, and it was wrong: "cold atmospheric plasma (no descriptor)" was
# published while `Plasma Gases` carried 474 census records and a ferroptosis
# intersection LARGER than a leg this report ranks.
CANDIDATE_DESCRIPTORS = {
    "tumour-treating fields": "electric stimulation therapy",
    "cold atmospheric plasma": "plasma gases",
}

NOT_MEASURABLE = [
    "tumour-treating fields (no descriptor; `Electric Stimulation Therapy` is broader)",
    "sonodynamic therapy specifically. `Ultrasonic Therapy` IS broader, but "
    "that is the smaller effect: measured, its precision is 90.6% (3 records "
    "of 32) while its recall is 46.0% (34 ferroptosis-SDT papers carry no "
    "such descriptor), so the count is an UNDER-estimate by roughly twofold "
    "and every ratio against it is an UPPER bound. An earlier version of this "
    "list stated the opposite; see analysis/atlas-descriptor-recall.md",
]


def modality_universe() -> set:
    """Every modality descriptor the committed partitions name."""
    d = json.loads(PARTITIONS.read_text())
    out = set()
    for spec in d.values():
        for side in ("pharmacological", "physical"):
            for x in spec.get(side, []):
                out.add(x.strip().lower())
    if not out:
        raise SystemExit(f"no descriptors in {PARTITIONS}")
    return out


def scan(universe: set) -> dict:
    inter = Counter()
    # THE THESIS LEGS ARE COUNTED WHETHER OR NOT THEY ARE IN THE UNIVERSE.
    # `drug resistance, neoplasm` is not a MODALITY, so it is legitimately
    # outside a modality universe -- but the table printed `counts.get(desc, 0)`
    # under a column headed "count", so absence rendered as a measured ZERO.
    # Its true intersection is published in the sibling atlas-thesis-position
    # artifact from the same build, which called it "the strong one".
    legs = {d.lower() for d in THESIS_LEGS.values()}
    leg_inter = Counter()
    # And the CENSUS-WIDE total for every descriptor we report, so a share can
    # be normalised by how common the descriptor is overall. Without that, a
    # ratio between an umbrella pharmacologic-action term and a specific
    # technique measures how broad the descriptors are.
    census_desc = Counter()
    cands = {d.lower() for d in CANDIDATE_DESCRIPTORS.values()}
    cand_inter = Counter()
    reported = universe | legs | cands   # every descriptor this report prints
    ferro_total = 0
    census = 0
    for f in sorted((ATLAS / "records").glob("*.jsonl.gz")):
        with gzip.open(f, "rt", encoding="utf-8") as fh:
            for line in fh:
                r = json.loads(line)
                census += 1
                mesh = {m.lower() for m in (r.get("mesh") or [])}
                for m in mesh & reported:
                    census_desc[m] += 1
                if FERROPTOSIS not in mesh:
                    continue
                ferro_total += 1
                for m in mesh & universe:
                    inter[m] += 1
                for m in mesh & legs:
                    leg_inter[m] += 1
                for m in mesh & cands:
                    cand_inter[m] += 1
    return {"census": census, "ferroptosis_total": ferro_total,
            "universe_size": len(universe),
            "intersections": [[k, v] for k, v in inter.most_common()],
            "leg_intersections": dict(leg_inter),
            "candidate_intersections": dict(cand_inter),
            "census_descriptor_totals": dict(census_desc)}


def render(d: dict) -> str:
    rows = d["intersections"]
    rank = {k: i + 1 for i, (k, _v) in enumerate(rows)}
    counts = dict(rows)
    L = ["# Where the thesis sits among all modalities", ""]
    L += ["*Generated by `scripts/atlas_thesis_rank.py`. The modality universe "
          "is the committed `analysis/modality-partitions.json`, built for #724 "
          "by five independent panels -- not a list written by the person "
          "making the claim.*", ""]

    L += [f"**{d['ferroptosis_total']:,} census articles carry MeSH "
          f"`Ferroptosis`.** Their intersection with each of "
          f"{d['universe_size']:,} modality descriptors, ranked:", ""]
    L += ["| rank | modality descriptor | ferroptosis articles | share of ferroptosis |",
          "|--:|---|--:|--:|"]
    for i, (k, v) in enumerate(rows[:20], 1):
        L.append(f"| {i} | {k} | {v:,} | "
                 f"{100*v/max(d['ferroptosis_total'],1):.2f}% |")
    L += [""]

    L += ["## Where the thesis's own legs land", ""]
    L += ["| leg | descriptor | count | rank of "
          f"{len(rows)} |", "|---|---|--:|--:|"]
    legs = d.get("leg_intersections") or {}
    for leg, desc in THESIS_LEGS.items():
        # The count comes from the LEG scan, not from the universe. An earlier
        # version used `counts.get(desc, 0)`, so a descriptor merely absent
        # from the modality universe printed as 0 under a column headed
        # "count" -- and `drug resistance, neoplasm` really is outside a
        # MODALITY universe while having a large intersection, which the
        # sibling atlas-thesis-position artifact publishes from the same build.
        c = legs.get(desc, counts.get(desc))
        r = rank.get(desc)
        cell = f"{c:,}" if c is not None else "not measured"
        rank_cell = (str(r) if r else
                     "*not a modality; outside this universe*")
        L.append(f"| {leg} | {desc} | {cell} | {rank_cell} |")
    L += [""]
    outside = [leg for leg, desc in THESIS_LEGS.items() if not rank.get(desc)]
    if outside:
        L += [f"{', '.join(outside)} " +
              ("is" if len(outside) == 1 else "are") +
              " outside the modality universe by construction -- the universe "
              "is built from therapeutic MODALITIES and these are not one. "
              "The count column is their real census intersection, measured "
              "here; an earlier version of this table printed a zero there, "
              "which reads as a measured absence rather than a scope "
              "boundary.", ""]

    top = rows[0] if rows else ("none", 0)
    sdt = counts.get(THESIS_LEGS["sonodynamic therapy"], 0)
    if sdt:
        L += [f"The largest ferroptosis-modality intersection is **{top[0]}** at "
              f"{top[1]:,} articles, which is **{top[1]/sdt:.1f}x** the "
              f"sonodynamic leg the simulation half rests on.", ""]

        # PREVALENCE NORMALISATION. The raw ratio divides an umbrella
        # pharmacologic-action descriptor by a specific technique, so it
        # partly measures how BROAD each descriptor is. Normalising by each
        # descriptor's census-wide prevalence asks a different and fairer
        # question -- is ferroptosis literature enriched or depleted in this
        # descriptor relative to the literature as a whole -- and it does not
        # give the same answer.
        tot = d.get("census_descriptor_totals") or {}
        ferro, cen = d["ferroptosis_total"], d["census"]
        base = ferro / cen if cen else None
        enrich = {}
        for name, c in ((top[0], top[1]),
                        (THESIS_LEGS["sonodynamic therapy"], sdt)):
            n_cen = tot.get(name)
            if n_cen and base:
                enrich[name] = (c / n_cen) / base
        if len(enrich) == 2:
            L += ["### The same comparison, normalised by how common each "
                  "descriptor is", ""]
            L += ["| descriptor | ferroptosis | census | share of its own "
                  "literature | vs base rate |", "|---|--:|--:|--:|--:|"]
            for name, c in ((top[0], top[1]),
                            (THESIS_LEGS["sonodynamic therapy"], sdt)):
                n_cen = tot.get(name)
                L.append(f"| {name} | {c:,} | {n_cen:,} | "
                         f"{100*c/n_cen:.3f}% | **{enrich[name]:.2f}x** |")
            L += ["",
                  f"The census base rate is {100*base:.3f}% "
                  f"({ferro:,} of {cen:,}). A value above 1 means ferroptosis "
                  f"literature is ENRICHED in that descriptor relative to the "
                  f"literature as a whole.", ""]
            a, b = enrich[top[0]], enrich[THESIS_LEGS["sonodynamic therapy"]]
            if b > a:
                L += [f"**The direction reverses.** `{top[0]}` is the larger "
                      f"COUNT and is *{a:.2f}x* its own base rate, while the "
                      f"sonodynamic descriptor is *{b:.2f}x* its own -- so the "
                      f"ferroptosis literature is relatively "
                      f"{'enriched' if b > 1 else 'less depleted'} in "
                      f"sonodynamic work and "
                      f"{'depleted' if a < 1 else 'less enriched'} in the "
                      f"umbrella term. The raw ratio above measures how broad "
                      f"the descriptors are at least as much as where "
                      f"attention goes, and it should not be read alone.", ""]
            else:
                L += [f"The direction holds under normalisation: "
                      f"`{top[0]}` is {a:.2f}x its base rate against the "
                      f"sonodynamic descriptor's {b:.2f}x.", ""]

    L += ["## Not measurable this way, and named rather than shown as zero", ""]
    for n in NOT_MEASURABLE:
        L.append(f"* {n}")
    cand = d.get("candidate_intersections") or {}
    tot2 = d.get("census_descriptor_totals") or {}
    checked = []
    for modality, desc in CANDIDATE_DESCRIPTORS.items():
        n_f, n_c = cand.get(desc, 0), tot2.get(desc, 0)
        checked.append((modality, desc, n_f, n_c))
    if checked:
        L += ["", "Each of those claims names the descriptor it rests on, and "
              "the claim is CHECKED rather than asserted:", ""]
        L += ["| modality | candidate descriptor | census | x ferroptosis |",
              "|---|---|--:|--:|"]
        for modality, desc, n_f, n_c in checked:
            L.append(f"| {modality} | `{desc}` | {n_c:,} | {n_f:,} |")
        L += [""]
        worst = min((v for _k, v in rows), default=0)
        # Only the candidates that STRICTLY exceed the smallest ranked entry
        # refute the claim. One that merely ties it is as thin as the thinnest
        # thing already shown, and saying otherwise would overstate the
        # correction in the same direction as the error.
        refuting = [c for c in checked if c[2] > worst]
        ties = [c for c in checked if 0 < c[2] <= worst]
        if refuting:
            L += ["**" + ("One of these is" if len(refuting) == 1
                          else f"{len(refuting)} of these are") +
                  " measurable after all.** " +
                  "; ".join(
                      f"`{d2}` carries {f:,} ferroptosis "
                      f"intersection{'' if f == 1 else 's'}"
                      for _m, d2, f, _c in refuting) +
                  f", against a smallest RANKED entry of {worst:,}. Listing a "
                  f"modality as having no descriptor, while a descriptor for "
                  f"it sits in the census with a LARGER intersection than "
                  f"something this report does rank, is a hand-written claim "
                  f"the data refutes.", ""]
        if ties:
            L += ["The remaining " +
                  "; ".join(f"`{d2}` ({f:,})" for _m, d2, f, _c in ties) +
                  f" does not exceed the smallest ranked entry ({worst:,}), so "
                  f"the not-measurable framing stands for "
                  f"{'it' if len(ties) == 1 else 'those'} -- the descriptor "
                  f"exists but carries no more ferroptosis literature than the "
                  f"thinnest thing already shown.", ""]
    L += [""]

    L += ["## What this does and does not establish", ""]
    L += ["* A large intersection is ATTENTION, not endorsement. It says a "
          "literature exists in which both concepts are indexed together, not "
          "that the combination works.",
          "* It does not replace the thesis-position analysis. That one asks "
          "whether the thesis's legs are thin; this one asks whether anything "
          "outside them is thick, which a hand-written list cannot answer.",
          "* Descriptor-to-descriptor throughout. An earlier version of this "
          "question compared a wide text-stem count against a narrow descriptor "
          "count and produced a ratio that reproduced under no definition.",
          "* The whole ranking is published so no single pair can be picked out "
          "to support a conclusion chosen first.",
          ""]
    return "\n".join(L) + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--render-only", action="store_true")
    args = ap.parse_args()
    if args.render_only:
        d = json.loads(OUT_JSON.read_text())
    else:
        d = scan(modality_universe())
        if d["ferroptosis_total"] == 0:
            raise SystemExit(
                "no ferroptosis articles found, which is not a finding -- it is "
                "what a descriptor-case mismatch looks like.")
        OUT_JSON.write_text(json.dumps(d, indent=1, sort_keys=True) + "\n",
                            encoding="utf-8")
    OUT_MD.write_text(render(d), encoding="utf-8")
    print(f"wrote {OUT_MD}")
    print(f"wrote {OUT_JSON}")
    for k, v in d["intersections"][:8]:
        print(f"  {v:>6,}  {k}")


if __name__ == "__main__":
    main()
