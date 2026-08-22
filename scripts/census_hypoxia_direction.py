#!/usr/bin/env python3
"""How the literature splits on the manuscript's most contested claim.

Section 7.1 reports that under simulated hypoxia a pharmacologic ferroptosis
inducer collapses while sonodynamic therapy is far less affected, and flags the
leg as genuinely contested -- the scientific review that produced that flag
found Zou 2019 reporting the opposite direction, hypoxia SENSITISING ccRCC to
RSL3, alongside work supporting the protective reading.

"Contested" is itself a claim about the literature, and the census can check it.
If the field reports one direction nine times out of ten, contested is the wrong
word and the manuscript is being over-cautious in a way that costs it. If the
split is near even, contested is right and the simulations are bounding a real
disagreement rather than hedging.

WHAT THIS CANNOT DO IS SETTLE THE BIOLOGY. A count of framings measures what
has been reported, and hypoxia's effect on ferroptosis genuinely differs by
cell type, by inducer, and by which arm of the pathway a study perturbs -- so a
majority direction is not the answer, and a near-even split is not confusion.
Both are facts about the literature, and the manuscript's own simulation work
is what bounds the conditions.

THE KEYWORD CLASSIFIER FAILED AND IS NOT THE MEASUREMENT. Adjudicating all 32
articles it labelled, it agrees only 34% of the time and REVERSES the direction
in 7 cases -- and six of those are the cancer-context ones, which is the subset
that matters here.

The mechanism is exact and worth stating because no amount of pattern-tuning
fixes it: **"hypoxia-induced ferroptosis resistance" CONTAINS "hypoxia-induced
ferroptosis"**. A phrase asserting that hypoxia PROTECTS contains, as a
substring, the phrase asserting that it SENSITISES. Proximity cannot separate
them, and the construction is common precisely in tumour biology, where
"ferroptosis resistance" is the thing being studied.

So the regex is demoted to what it can do -- generate candidates -- and the
measurement is the committed adjudication in
`analysis/hypoxia-direction-adjudication.csv`. That reverses the answer: the
regex reported the sensitising direction leading 62% to 38%, and the
adjudication puts the protective direction ahead 61% to 39%.
"""
import argparse
import gzip
import json
import re
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
RECORDS = REPO / "corpus/atlas/records"
OUT_MD = REPO / "analysis/census-hypoxia-direction.md"
OUT_JSON = REPO / "analysis/census-hypoxia-direction.json"

FERRO = "ferroptosis"
# The hypoxia family, taken as a union rather than one descriptor: no single
# term carries this literature, and picking one would be picking a result.
HYPOXIA = {"cell hypoxia", "tumor hypoxia", "hypoxia", "anoxia",
           "hypoxia-inducible factor 1, alpha subunit", "hypoxia-inducible factor 1"}

# FIXED BEFORE THE RESULT WAS READ.
# Hypoxia makes cells HARDER to kill by ferroptosis -- the direction this
# project's simulation assumes for its pharmacologic arm.
PROTECTS = re.compile(
    r"\b(hypoxia[- ]?(?:induced|mediated|driven)? ?(?:resistance|protect\w*|"
    r"inhibit\w*|suppress\w*|attenuat\w*|reduc\w*)|"
    r"resistan\w* to ferroptosis under hypoxia|"
    r"protect\w* (?:cells |tumou?r cells )?(?:from|against) ferroptosis)\b")
# Hypoxia makes cells EASIER to kill by ferroptosis -- the Zou 2019 direction.
SENSITISES = re.compile(
    r"\b(hypoxia[- ]?(?:induced|mediated|driven)? ?(?:sensitiz\w*|sensitis\w*|"
    r"promot\w*|enhanc\w*|potentiat\w*|trigger\w*|aggravat\w*)|"
    r"hypoxia[- ]induced ferroptosis|"
    r"ferroptosis under hypoxia\w* condition)\b")
N_SAMPLE = 6
# Below this the split is a statement about a handful.
MIN_CLASSIFIED = 20


def scan(stride: int = 1) -> dict:
    counts = {"protects": 0, "sensitises": 0, "both": 0, "neither": 0}
    sample = []
    total = 0
    for f in sorted(RECORDS.glob("*.jsonl.gz"))[::stride]:
        with gzip.open(f, "rt", encoding="utf-8") as fh:
            for line in fh:
                r = json.loads(line)
                ms = {m.lower() for m in (r.get("mesh") or [])}
                if FERRO not in ms or not (ms & HYPOXIA):
                    continue
                total += 1
                blob = f"{r.get('title') or ''} {r.get('abstract') or ''}".lower()
                p, s = bool(PROTECTS.search(blob)), bool(SENSITISES.search(blob))
                key = ("both" if p and s else "protects" if p
                       else "sensitises" if s else "neither")
                counts[key] += 1
                if key in ("protects", "sensitises") and len(sample) < N_SAMPLE:
                    sample.append({"framing": key, "year": r.get("year"),
                                   "title": (r.get("title") or "")[:105]})
    return {"total": total, "counts": counts, "sample": sample,
            "hypoxia_descriptors": sorted(HYPOXIA),
            "min_classified": MIN_CLASSIFIED}


def load_adjudication() -> dict:
    """The measurement. Titles only, one adjudicator, and both limits stated.

    Weaker evidence than a scaled classifier would be -- and stronger than the
    classifier actually is, which was the choice available.
    """
    import csv

    path = REPO / "analysis/hypoxia-direction-adjudication.csv"
    if not path.exists():
        return {}
    rows = list(csv.DictReader(path.open(encoding="utf-8")))
    from collections import Counter

    adj = Counter(r["adjudicated"] for r in rows)
    agree = sum(1 for r in rows if r["regex_label"] == r["adjudicated"])
    reversed_ = [r for r in rows
                 if r["adjudicated"] in ("protects", "sensitises")
                 and r["regex_label"] != r["adjudicated"]]
    return {
        "n": len(rows),
        "labels": dict(adj),
        "regex_agreement": round(100 * agree / len(rows), 1) if rows else None,
        "regex_reversed": len(reversed_),
        "regex_reversed_cancer": sum(1 for r in reversed_
                                     if "cancer" in r["reason"].lower()),
        "directional": adj["protects"] + adj["sensitises"],
        "protects": adj["protects"],
        "sensitises": adj["sensitises"],
    }


def _wilson(k: int, n: int):
    import math

    if not n:
        return None
    z = 1.96
    ph = k / n
    den = 1 + z * z / n
    centre = (ph + z * z / (2 * n)) / den
    half = z * math.sqrt(ph * (1 - ph) / n + z * z / (4 * n * n)) / den
    return [round(100 * max(0.0, centre - half), 1),
            round(100 * min(1.0, centre + half), 1)]


def assemble(d: dict) -> dict:
    c = d["counts"]
    single = c["protects"] + c["sensitises"]
    out = dict(d)
    out["classified"] = single
    out["unclassified"] = c["neither"] + c["both"]
    out["unclassified_share"] = (round(100 * out["unclassified"] / d["total"], 1)
                                 if d["total"] else None)
    out["protects_share"] = (round(100 * c["protects"] / single, 1)
                             if single else None)
    out["interpretable"] = single >= d["min_classified"]
    # A POINT ESTIMATE ON THIRTY-TWO ARTICLES SHOULD NOT DRIVE A VERDICT.
    # Wilson interval rather than normal-approximation: at this n the normal
    # interval is both too narrow and capable of running past 0 or 1.
    if single:
        import math

        z = 1.96
        ph = c["protects"] / single
        den = 1 + z * z / single
        centre = (ph + z * z / (2 * single)) / den
        half = z * math.sqrt(ph * (1 - ph) / single
                             + z * z / (4 * single * single)) / den
        out["protects_ci"] = [round(100 * max(0.0, centre - half), 1),
                              round(100 * min(1.0, centre + half), 1)]
    else:
        out["protects_ci"] = None
    # THE VERDICT ON THE WORD, not on the biology. "Contested" is a claim about
    # the literature: if one direction takes more than this share, the field is
    # not evenly split and the manuscript's hedge is doing something other than
    # reflecting a disagreement.
    out["dominance_threshold"] = 75.0
    # THE ADJUDICATION IS THE MEASUREMENT; the regex counts above are kept only
    # to document what it got wrong.
    adj = load_adjudication()
    out["adjudication"] = adj
    if adj:
        out["adj_protects_share"] = (
            round(100 * adj["protects"] / adj["directional"], 1)
            if adj["directional"] else None)
        out["adj_ci"] = _wilson(adj["protects"], adj["directional"])
        out["regex_direction_reversed_by_adjudication"] = bool(
            out["adj_protects_share"] and out["protects_share"]
            and (out["adj_protects_share"] > 50) != (out["protects_share"] > 50))
    def band(x):
        if x >= out["dominance_threshold"]:
            return "protective direction dominates"
        if x <= 100 - out["dominance_threshold"]:
            return "sensitising direction dominates"
        return "genuinely split"

    if not out["interpretable"] or out["protects_share"] is None:
        out["verdict"] = "not measurable"
        out["verdict_survives_interval"] = False
    else:
        out["verdict"] = band(out["protects_share"])
        # THE VERDICT MUST SURVIVE ITS OWN INTERVAL. A point estimate landing
        # in the split band while its interval reaches a dominance band is not
        # a finding of a split -- it is a finding that the data cannot tell.
        lo, hi = out["protects_ci"]
        out["verdict_survives_interval"] = band(lo) == band(hi) == out["verdict"]
    return out


def render(d: dict) -> str:
    c = d["counts"]
    L = ["# How the literature splits on hypoxia and ferroptosis\n"]
    L.append(
        f"Generated by `scripts/census_hypoxia_direction.py`. "
        f"{d['total']:,} census articles carry `Ferroptosis` together with a "
        f"hypoxia-family descriptor ("
        + ", ".join(f"`{h}`" for h in d["hypoxia_descriptors"])
        + "). Taken as a union rather than one descriptor, because no single "
          "term carries this literature and picking one would be picking a "
          "result.\n"
    )
    L.append(
        "Section 7.1 reports a pharmacologic ferroptosis inducer collapsing "
        "under simulated hypoxia while sonodynamic therapy is far less "
        "affected, and flags the leg as genuinely contested. **Contested is "
        "itself a claim about the literature**, and it is checkable: if the "
        "field reports one direction nine times in ten, the word is wrong and "
        "the manuscript is hedging where it need not.\n"
    )
    L.append("| framing | articles | share of all |")
    L.append("|---|--:|--:|")
    for k, label in (
            ("protects", "hypoxia PROTECTS against ferroptosis"),
            ("sensitises", "hypoxia SENSITISES to ferroptosis"),
            ("both", "both directions present"),
            ("neither", "neither pattern matched")):
        share = f"{100 * c[k] / d['total']:.0f}%" if d["total"] else "-"
        L.append(f"| {label} | {c[k]:,} | {share} |")
    L.append("")
    a = d.get("adjudication") or {}
    if not a:
        L.append("**The adjudication file is missing**, so only the keyword "
                 "classifier's counts are shown above -- and that classifier "
                 "is known to reverse the direction. Do not read the table as "
                 "a result.\n")
    else:
        lo, hi = d["adj_ci"]
        L.append("## The keyword classifier failed, and the failure is instructive\n")
        L.append(
            f"All {a['n']} articles the classifier labelled were adjudicated "
            f"from their titles. It agrees **{a['regex_agreement']}%** of the "
            f"time and REVERSES the direction in **{a['regex_reversed']}** "
            f"cases -- {a['regex_reversed_cancer']} of them cancer-context, "
            f"which is the subset that matters here.\n"
        )
        L.append(
            "The mechanism is exact, and no amount of pattern-tuning fixes it: "
            "**\"hypoxia-induced ferroptosis resistance\" CONTAINS "
            "\"hypoxia-induced ferroptosis\"**. A phrase asserting that "
            "hypoxia PROTECTS contains, as a substring, the phrase asserting "
            "that it SENSITISES. Proximity cannot separate them, and the "
            "construction is common precisely in tumour biology, where "
            "\"ferroptosis resistance\" is the thing being studied.\n"
        )
        L.append(
            f"So the regex is demoted to generating candidates and the "
            f"adjudication is the measurement. **It reverses the answer**: the "
            f"classifier reported the sensitising direction leading "
            f"{100 - d['protects_share']:.0f}% to {d['protects_share']:.0f}%, "
            f"and the adjudication puts the PROTECTIVE direction ahead "
            f"{d['adj_protects_share']}% to "
            f"{100 - d['adj_protects_share']:.0f}% over "
            f"{a['directional']} directional articles (95% Wilson interval "
            f"{lo}-{hi}%).\n"
        )
        L.append(
            f"That interval is wide enough to contain an even split, so what "
            f"this establishes is that BOTH directions are reported in "
            f"substantial numbers and NOT their ratio. Section 7.1's "
            f"'contested' framing therefore reflects the literature. The "
            f"protective direction -- the one this project's simulation "
            f"assumes for its pharmacologic arm -- is the one that leads on "
            f"the adjudicated count, but leads by less than the uncertainty.\n"
        )
        L.append(
            f"{a['labels'].get('off-topic', 0)} of {a['n']} are off-topic "
            f"altogether (a drug's effect in a hypoxic context, rather than "
            f"hypoxia's own), and {a['labels'].get('ambiguous', 0)} state no "
            f"direction in the title. Those are excluded rather than assigned, "
            f"which is why the directional denominator is {a['directional']} "
            f"and not {a['n']}.\n"
        )
        L.append("## Limits of the adjudication\n")
        L.append(
            "Titles only, one adjudicator, and 18 directional articles. That "
            "is weaker evidence than a working classifier at scale would be. "
            "It is stronger than the classifier actually is, which was the "
            "choice available -- and every label is committed in "
            "`analysis/hypoxia-direction-adjudication.csv` with its reason, so "
            "a reader can disagree with any of them individually.\n"
        )
    if d["sample"]:
        L.append("### Examples of each framing\n")
        for s in d["sample"]:
            L.append(f"- *{s['framing']}* — {s['year']} — {s['title']}")
        L.append("")
    L.append("## What this cannot do\n")
    L.append(
        "Settle the biology. Hypoxia's effect on ferroptosis genuinely differs "
        "by cell type, by inducer, and by which arm of the pathway a study "
        "perturbs, so a majority direction is not the answer and an even split "
        "is not confusion. Both are facts about what has been REPORTED. The "
        "manuscript's simulation work is what bounds the conditions, and this "
        "measures only whether the word attached to it describes the "
        "literature.\n"
    )
    L.append(
        f"{d['unclassified_share']}% of the intersection matches neither "
        f"pattern or both, so the split is over the classified remainder. The "
        f"patterns were fixed before any result was read and are not adjusted "
        f"afterwards: this project's simulation ASSUMES the protective "
        f"direction for its pharmacologic arm, so a vocabulary widened after "
        f"seeing the split would be fitted to the answer the project has "
        f"already built.\n"
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
    print(f"  {d['total']} articles, {d['classified']} classified, "
          f"protects {d['protects_share']}% -> {d['verdict']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
