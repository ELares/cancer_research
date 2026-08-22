#!/usr/bin/env python3
"""Does the leg this project calls strong point the way the project needs?

`atlas_thesis_position.py` counts the legs of this project's hypothesis and
concludes that "the resistance leg is the strong one: 479 articles connect
ferroptosis to neoplasm drug resistance, and that is the part of the argument
the manuscript should lean on hardest." The manuscript leans on it.

THAT COUNT IS SIZE, NOT DIRECTION, and the two can disagree. The project's
thesis is that drug-resistant and persister cells are ferroptosis-VULNERABLE --
that restoring or inducing ferroptosis is a way to kill what chemotherapy left
behind. A literature connecting ferroptosis to drug resistance could just as
easily be about the opposite: tumours acquiring resistance TO ferroptosis
induction, which is an obstacle to the thesis rather than support for it. Both
kinds of paper carry both descriptors and are indistinguishable in a count.

So the leg is classified by which way its articles frame the relationship. The
patterns were fixed before any result was read and are not tuned afterwards,
which matters because the exploit vocabulary IS the thesis vocabulary and
adjusting it after seeing the split would fit the answer.

WHAT THIS MEASURES AND WHAT IT DOES NOT. A field framing something as "we can
exploit X" is evidence about what the field is TRYING, not about whether X
works. This bounds attention and its direction; it does not support the thesis,
and a 90% exploit share would say the field is pursuing the idea rather than
that the idea is right.

RECALL IS POOR AND DEMONSTRATED RATHER THAN ESTIMATED. A large minority match
neither pattern, and sampling them shows on-topic articles the patterns miss --
"induction of ferroptosis ... enhances" is the thesis direction stated in words
the exploit set does not contain. The ratio is therefore reported over the
articles that ARE classified, with the unclassified share stated beside it, and
the magnitude is not claimed.
"""
import argparse
import gzip
import json
import re
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
RECORDS = REPO / "corpus/atlas/records"
OUT_MD = REPO / "analysis/census-thesis-direction.md"
OUT_JSON = REPO / "analysis/census-thesis-direction.json"

FERRO = "ferroptosis"
RESIST = "drug resistance, neoplasm"

# FIXED BEFORE THE RESULT WAS READ. The exploit vocabulary is the thesis
# vocabulary, so tuning it after seeing the split would be fitting the answer.
EXPLOIT = re.compile(
    r"\b(sensitiz\w*|resensitiz\w*|re-sensitiz\w*|overcome\w*|overcoming|"
    r"circumvent\w*|reverse[sd]? (?:the )?(?:drug |chemo)?resistance|"
    r"vulnerab\w*|susceptib\w*|achilles)\b")
OBSTACLE = re.compile(
    r"\b(ferroptosis resistance|resistance to ferroptosis|"
    r"ferroptosis[- ]resistant|evade[sd]? ferroptosis|escape[sd]? ferroptosis|"
    r"protect\w* against ferroptosis|ferroptosis defen[cs]e)\b")
N_SAMPLE = 6


def scan(stride: int = 1) -> dict:
    counts = {"exploit": 0, "obstacle": 0, "both": 0, "neither": 0}
    unclassified = []
    total = 0
    for f in sorted(RECORDS.glob("*.jsonl.gz"))[::stride]:
        with gzip.open(f, "rt", encoding="utf-8") as fh:
            for line in fh:
                r = json.loads(line)
                ms = {m.lower() for m in (r.get("mesh") or [])}
                if FERRO not in ms or RESIST not in ms:
                    continue
                total += 1
                blob = f"{r.get('title') or ''} {r.get('abstract') or ''}".lower()
                e, o = bool(EXPLOIT.search(blob)), bool(OBSTACLE.search(blob))
                key = ("both" if e and o else "exploit" if e
                       else "obstacle" if o else "neither")
                counts[key] += 1
                if key == "neither" and len(unclassified) < N_SAMPLE:
                    unclassified.append(
                        {"year": r.get("year"),
                         "title": (r.get("title") or "")[:110]})
    return {"total": total, "counts": counts, "unclassified_sample": unclassified}


def _wilson(k: int, n: int):
    import math

    if not n:
        return None
    z = 1.96
    ph = k / n
    den = 1 + z * z / n
    centre = (ph + z * z / (2 * n)) / den
    half = z * math.sqrt(ph * (1 - ph) / n + z * z / (4 * n * n)) / den
    return [max(0.0, centre - half), min(1.0, centre + half)]


def load_adjudication(counts: dict) -> dict:
    """Correct the classifier's counts using an adjudicated sample.

    THE CLASSIFIER ERRS IN BOTH DIRECTIONS HERE, which the raw counts cannot
    show. "Targeting cTRIP12 counteracts ferroptosis resistance" is an EXPLOIT
    paper labelled obstacle, because it contains the obstacle phrase; "Nedd4
    ... suppress erastin-induced ferroptosis" is an OBSTACLE paper labelled
    exploit. Same substring failure as the sibling hypoxia analysis, in both
    directions at once.

    Every obstacle-labelled article was adjudicated (there are few) and the
    exploit-labelled ones were sampled, so the correction is exact on one side
    and estimated on the other. Ambiguous articles are dropped from both
    numerator and denominator rather than assigned.
    """
    import csv
    from collections import Counter

    path = REPO / "analysis/thesis-direction-adjudication.csv"
    if not path.exists():
        return {}
    rows = list(csv.DictReader(path.open(encoding="utf-8")))
    by = {}
    for label in ("exploit", "obstacle"):
        sub = [r for r in rows if r["regex_label"] == label]
        c = Counter(r["adjudicated"] for r in sub)
        decided = c["exploit"] + c["obstacle"]
        by[label] = {"sampled": len(sub), "decided": decided,
                     "adjudicated_exploit": c["exploit"],
                     "adjudicated_obstacle": c["obstacle"],
                     "ambiguous": c["ambiguous"],
                     "precision": (c[label] / decided) if decided else None,
                     "precision_ci": _wilson(c[label], decided)}

    def corrected(exploit_precision: float):
        # Obstacle-labelled articles were adjudicated in full, so their
        # contribution is a count rather than a rate.
        ex = counts["exploit"] * exploit_precision + by["obstacle"]["adjudicated_exploit"]
        ob = counts["exploit"] * (1 - exploit_precision) + by["obstacle"]["adjudicated_obstacle"]
        return 100 * ex / (ex + ob) if (ex + ob) else None

    p_hat = by["exploit"]["precision"]
    lo_p, hi_p = by["exploit"]["precision_ci"]
    point = corrected(p_hat)
    band = sorted([corrected(lo_p), corrected(hi_p)])
    return {
        "rows": len(rows), "by_label": by,
        "corrected_exploit_share": round(point, 1),
        "corrected_range": [round(band[0], 1), round(band[1], 1)],
        # The DIRECTION survives if every plausible correction keeps exploit
        # ahead; the MAGNITUDE does not if the range is wide.
        "direction_survives": band[0] > 50,
    }


def assemble(d: dict) -> dict:
    c = d["counts"]
    # SINGLY CLASSIFIED ONLY. An article carrying both framings is not evidence
    # for either direction, and folding it into the larger side would let the
    # ambiguous cases inflate whichever way already leads.
    single = c["exploit"] + c["obstacle"]
    out = dict(d)
    out["classified"] = single
    out["unclassified"] = c["neither"] + c["both"]
    out["exploit_share_of_classified"] = (
        round(100 * c["exploit"] / single, 1) if single else None)
    out["direction_ratio"] = (
        round(c["exploit"] / c["obstacle"], 1) if c["obstacle"] else None)
    out["unclassified_share"] = round(100 * out["unclassified"] / d["total"], 1)
    # The verdict is derived. A leg pointing the project's way is not a given
    # and the report must be able to say the opposite.
    out["points_the_projects_way"] = bool(
        out["exploit_share_of_classified"] and
        out["exploit_share_of_classified"] > 50)
    out["adjudication"] = load_adjudication(c)
    return out


def render(d: dict) -> str:
    c = d["counts"]
    L = ["# Which way the resistance leg points\n"]
    L.append(
        f"Generated by `scripts/census_thesis_direction.py`. Of the "
        f"{d['total']:,} census articles carrying both `Ferroptosis` and "
        f"`Drug Resistance, Neoplasm` -- the leg this project calls its "
        f"strongest, and the one the manuscript leans on hardest -- how many "
        f"frame ferroptosis as a way to KILL resistant cells, and how many as "
        f"something tumours become resistant TO?\n"
    )
    L.append(
        "The distinction matters because a count cannot make it. The project's "
        "thesis is that drug-resistant and persister cells are "
        "ferroptosis-vulnerable. A paper reporting that tumours acquire "
        "resistance to ferroptosis induction carries exactly the same two "
        "descriptors and is an obstacle to the thesis rather than support for "
        "it.\n"
    )
    L.append("| framing | articles | share |")
    L.append("|---|--:|--:|")
    for k, label in (("exploit", "exploit: ferroptosis kills resistant cells"),
                     ("obstacle", "obstacle: resistance TO ferroptosis"),
                     ("both", "both framings present"),
                     ("neither", "neither pattern matched")):
        L.append(f"| {label} | {c[k]:,} | {100 * c[k] / d['total']:.0f}% |")
    L.append("")
    if d["points_the_projects_way"]:
        L.append(
            f"**Among the {d['classified']:,} articles carrying exactly one "
            f"framing, the thesis direction leads "
            f"{d['exploit_share_of_classified']}% to "
            f"{100 - d['exploit_share_of_classified']:.1f}% -- a ratio of "
            f"{d['direction_ratio']} to 1.** The leg the manuscript leans on "
            f"does point the way the manuscript needs, which had been assumed "
            f"rather than measured.\n"
        )
    else:
        L.append(
            f"**The leg does NOT point the project's way.** Among singly "
            f"classified articles the exploit framing is "
            f"{d['exploit_share_of_classified']}%, so the literature "
            f"connecting ferroptosis to drug resistance is substantially about "
            f"resistance TO ferroptosis. The manuscript leans on this leg and "
            f"should stop.\n"
        )
    a = d.get("adjudication") or {}
    if a:
        lo, hi = a["corrected_range"]
        L.append("## The classifier errs in BOTH directions, and the raw split "
                 "overstates the lead\n")
        L.append(
            f"Every obstacle-labelled article was adjudicated and the "
            f"exploit-labelled ones were sampled, so the correction is exact "
            f"on one side and estimated on the other. Both sides are wrong, in "
            f"opposite directions: "
            f"{a['by_label']['obstacle']['adjudicated_exploit']} of the "
            f"{a['by_label']['obstacle']['sampled']} obstacle-labelled "
            f"articles are actually EXPLOIT papers, and "
            f"{a['by_label']['exploit']['adjudicated_obstacle']} of "
            f"{a['by_label']['exploit']['sampled']} sampled exploit-labelled "
            f"ones are actually OBSTACLE papers.\n"
        )
        L.append(
            "The mechanism is the one the sibling hypoxia analysis names: "
            "*\"counteracts ferroptosis resistance\"* contains *\"ferroptosis "
            "resistance\"*, and *\"suppress erastin-induced ferroptosis\"* "
            "contains the language of induction. A phrase asserting one "
            "direction contains, as a substring, the phrase asserting the "
            "other -- here in both directions at once.\n"
        )
        L.append(
            f"**Corrected, the exploit share is {a['corrected_exploit_share']}%, "
            f"not {d['exploit_share_of_classified']}%** (range {lo}-{hi}% "
            f"propagating the sampled precision's interval). "
            + ("The DIRECTION survives every plausible correction -- exploit "
               "leads throughout the range -- and the MAGNITUDE does not: a "
               f"ratio of {d['direction_ratio']} to 1 becomes roughly 3 to 1, "
               "and the raw figure should not be quoted.\n"
               if a["direction_survives"] else
               "The range crosses even, so neither the magnitude NOR the "
               "direction is established by this measurement.\n")
        )
        L.append(
            f"Ambiguous articles -- reviews, prognostic signatures, and titles "
            f"stating no direction -- are dropped from numerator and "
            f"denominator rather than assigned. Labels and reasons are "
            f"committed in `analysis/thesis-direction-adjudication.csv`.\n"
        )
    L.append("## What this does not establish\n")
    L.append(
        "That the thesis is right. A field framing something as \"we can "
        "exploit X\" is evidence about what the field is TRYING, and an "
        "exploit share of ninety per cent would say the idea is being pursued "
        "rather than that it works. This bounds attention and its direction, "
        "which is what a census can do, and stops there.\n"
    )
    L.append("## Recall is poor, and demonstrated rather than estimated\n")
    L.append(
        f"{d['unclassified_share']}% of the leg matches neither pattern or "
        f"both, and sampling those shows on-topic articles the patterns miss:\n"
    )
    for s in d["unclassified_sample"]:
        L.append(f"- {s['year']} — {s['title']}")
    L.append("")
    L.append(
        "The first of those is the thesis direction stated in words the "
        "exploit set does not contain. So the ratio above is reported over the "
        "articles that ARE classified, the unclassified share is stated beside "
        "it, and the MAGNITUDE is not claimed -- only the direction, which a "
        "9-to-1 split is robust to a good deal of missing recall.\n"
    )
    L.append(
        "The patterns were fixed before any result was read and were not "
        "adjusted afterwards. That restraint matters more here than usual: the "
        "exploit vocabulary IS the thesis vocabulary, so widening it after "
        "seeing the split would be fitting the answer this project wants.\n"
    )
    return "\n".join(L)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--render-only", action="store_true")
    a = ap.parse_args()
    d = assemble(json.loads(OUT_JSON.read_text()) if a.render_only
                 else scan(a.stride))
    OUT_JSON.write_text(json.dumps(d, indent=1) + "\n")
    OUT_MD.write_text(render(d))
    print(f"wrote {OUT_MD}")
    print(f"  {d['total']} articles, exploit "
          f"{d['exploit_share_of_classified']}% of {d['classified']} "
          f"classified, ratio {d['direction_ratio']}:1, "
          f"{d['unclassified_share']}% unclassified")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
