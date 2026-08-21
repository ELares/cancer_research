#!/usr/bin/env python3
"""How long between a mechanism's literature appearing and a trial appearing.

A census makes one translational question answerable that a sample cannot: for
each therapeutic mechanism, how many years pass between the literature starting
and a CLINICAL TRIAL being indexed for it. Both ends come from NLM -- the
articles from MeSH descriptors or from the mechanism's own vocabulary, the
trials from publication types -- so neither end is a judgement made here.

THE CONFOUND IS LARGER THAN THE EFFECT AND IS THE REASON THIS RUNS TWO ARMS.
A MeSH descriptor has an introduction date. `Ferroptosis` became a descriptor in
2020, so no article carries it before then however long the science existed, and
a lag measured from a descriptor's first appearance is partly a measurement of
when NLM minted the term. So every mechanism is measured twice:

  MeSH arm  -- descriptors, expert-assigned, precise, and blind to anything
               before the descriptor existed.
  TEXT arm  -- this project's own keyword vocabulary over title and abstract
               only, which can see a concept from the moment authors named it.

The GAP BETWEEN THE TWO FIRST-YEARS is the confound, measured per mechanism
rather than assumed uniform. Where the arms agree, a lag means what it appears
to mean; where the text arm starts much earlier, the MeSH lag is compressed by
the descriptor's own history.

The text arm also reaches three mechanisms MeSH cannot express at all --
TTFields, bioelectric modulation, cold atmospheric plasma -- which are reported
as *not measurable* everywhere else in this project. They are measurable here,
on the text arm alone, and that is stated on every row rather than left for a
reader to infer from a missing cell.

WHAT A LAG IS NOT. It is not time-to-approval, and it is not evidence that a
mechanism translated well: a trial being indexed says a trial happened, not that
it succeeded. Mechanisms with no trial yet are RIGHT-CENSORED and are reported
as censored rather than as a large lag, because a mechanism that never reaches
a trial would otherwise score the same as one that reached it slowly.
"""
import argparse
import gzip
import json
import re
import statistics
from collections import Counter, defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
RECORDS = REPO / "corpus/atlas/records"
MECH_MAP = REPO / "analysis/mesh-mechanism-map.yaml"
OUT_MD = REPO / "analysis/census-translation-lag.md"
OUT_JSON = REPO / "analysis/census-translation-lag.json"

TRIAL_TYPES = {
    "Clinical Trial", "Randomized Controlled Trial", "Controlled Clinical Trial",
    "Clinical Trial, Phase I", "Clinical Trial, Phase II",
    "Clinical Trial, Phase III", "Clinical Trial, Phase IV",
    "Pragmatic Clinical Trial", "Adaptive Clinical Trial",
}
# THE THRESHOLDS MUST MATCH, and a first version's did not. Requiring 5
# articles to call a literature started while accepting 1 trial as the event
# produced trials BEFORE their own literature -- TTFields showed a 2007 trial
# against a 2013 start -- because a mechanism's first trial routinely arrives
# before its fifth paper. That is an asymmetric admission rule, and the lag it
# yields is not a duration.
#
# A trial record IS an article record, so at a matched threshold of 1 the first
# trial cannot precede the first article: the ordering holds by construction
# rather than by luck, and a negative lag becomes impossible instead of being
# silently dropped.
FIRST_MIN = 1
# The threshold the first version used, kept as a ROBUSTNESS column rather than
# as the definition. A single stray match can start the clock early, and the
# distance between the two starts is how much that matters per mechanism.
STABLE_MIN = 5
# Above this many years between the two thresholds, ONE early match is carrying
# the start year and the lag beside it is a lower bound rather than a duration.
FRAGILE_AT = 5
LAST_FULL_YEAR = 2025


def _load_mesh():
    import yaml

    mp = yaml.safe_load(MECH_MAP.read_text(encoding="utf-8"))["mechanisms"]
    # Keys LOWERCASED. The two vocabularies disagree on case -- config spells it
    # `mRNA-vaccine` and the MeSH map `mrna-vaccine` -- which produced two rows
    # for one mechanism, each with one empty arm, looking like two mechanisms
    # that each happened to be measurable on only one side. The same case
    # mismatch once zeroed this exact mechanism in atlas_landscape.
    return {k.lower(): {x.lower() for x in v["descriptors"]}
            for k, v in mp.items() if v["descriptors"]}


def _load_text():
    import sys

    sys.path.insert(0, str(REPO / "scripts"))
    import config

    # WORD-BOUNDED. Unbounded substrings dated `car-t` to 1955 by matching
    # inside "s(car t)issue" and `electrolysis` to 1953 via "li(echt)enstein".
    # A first-appearance statistic is decided entirely by the EARLIEST match, so
    # unlike a prevalence estimate it has no error averaging: one false positive
    # in four million records sets the answer.
    out = {}
    for k, terms in config.MECHANISM_KEYWORDS.items():
        if terms:
            out[k.lower()] = re.compile(
                "|".join(rf"\b{re.escape(t.lower())}\b" for t in terms))
    return out


def scan(stride: int = 1) -> dict:
    mesh = _load_mesh()
    text = _load_text()
    # One alternation over EVERY keyword as a prefilter. Most records match no
    # mechanism, and running 25 patterns over 4.4M records would dominate the
    # runtime for nothing.
    prefilter = re.compile("|".join(p.pattern for p in text.values()))

    arms = {"mesh": defaultdict(Counter), "text": defaultdict(Counter)}
    trials = {"mesh": defaultdict(Counter), "text": defaultdict(Counter)}
    n = 0
    for f in sorted(RECORDS.glob("*.jsonl.gz"))[::stride]:
        with gzip.open(f, "rt", encoding="utf-8") as fh:
            for line in fh:
                r = json.loads(line)
                n += 1
                y = r.get("year")
                if not isinstance(y, int) or y > LAST_FULL_YEAR:
                    continue
                is_trial = bool(set(r.get("pub_types") or []) & TRIAL_TYPES)

                ms = {m.lower() for m in (r.get("mesh") or [])}
                if ms:
                    for k, d in mesh.items():
                        if ms & d:
                            arms["mesh"][k][y] += 1
                            if is_trial:
                                trials["mesh"][k][y] += 1

                # Title and abstract ONLY. Folding MeSH into the text arm would
                # make the two arms partly the same instrument, and the gap
                # between them is the whole measurement.
                blob = f"{r.get('title') or ''} {r.get('abstract') or ''}".lower()
                if blob and prefilter.search(blob):
                    for k, pat in text.items():
                        if pat.search(blob):
                            arms["text"][k][y] += 1
                            if is_trial:
                                trials["text"][k][y] += 1
    return {
        "census": n,
        "first_min": FIRST_MIN,
        "stable_min": STABLE_MIN,
        "last_full_year": LAST_FULL_YEAR,
        "mesh_measurable": sorted(mesh),
        "text_measurable": sorted(text),
        "arms": {a: {k: {str(y): c for y, c in sorted(v.items())}
                     for k, v in arms[a].items()} for a in arms},
        "trials": {a: {k: {str(y): c for y, c in sorted(v.items())}
                       for k, v in trials[a].items()} for a in trials},
    }


def _first_year(series: dict, threshold: int):
    for y in sorted(series, key=int):
        if series[y] >= threshold:
            return int(y)
    return None


def assemble(d: dict) -> dict:
    rows = []
    for k in sorted(set(d["arms"]["mesh"]) | set(d["arms"]["text"])):
        row = {"mechanism": k,
               "mesh_measurable": k in d["mesh_measurable"],
               "text_measurable": k in d["text_measurable"]}
        for arm in ("mesh", "text"):
            series = d["arms"][arm].get(k, {})
            tseries = d["trials"][arm].get(k, {})
            start = _first_year(series, d["first_min"])
            trial = _first_year(tseries, d["first_min"])
            row[f"{arm}_start"] = start
            row[f"{arm}_stable_start"] = _first_year(series, d["stable_min"])
            row[f"{arm}_first_trial"] = trial
            row[f"{arm}_articles"] = sum(series.values())
            row[f"{arm}_trials"] = sum(tseries.values())
            # RIGHT-CENSORED, not a large lag. A mechanism that has not reached
            # a trial must not score like one that reached it slowly.
            row[f"{arm}_censored"] = bool(start and not trial)
            if start and trial:
                # Guaranteed non-negative: a trial record is an article record,
                # so at a matched threshold the ordering cannot invert. An
                # assertion rather than a silent None, because if this ever
                # fails the thresholds have drifted apart again.
                assert trial >= start, (k, arm, start, trial)
                row[f"{arm}_lag"] = trial - start
            else:
                row[f"{arm}_lag"] = None
            # How much a single stray match moves the clock, per mechanism.
            row[f"{arm}_start_fragility"] = (
                row[f"{arm}_stable_start"] - start
                if (start and row[f"{arm}_stable_start"]) else None)
        # NOT a measure of when a descriptor was minted, which is what a first
        # version called it. A negative value means the DESCRIPTOR is older or
        # broader than the term -- `Ultrasonic Therapy` runs from 1955 while
        # the word "sonodynamic" is recent -- so this mixes introduction date
        # with descriptor breadth and is named for what it measures.
        row["arm_start_gap"] = (
            row["mesh_start"] - row["text_start"]
            if (row["mesh_start"] and row["text_start"]) else None)
        rows.append(row)

    lags = [r["text_lag"] for r in rows if r["text_lag"] is not None]
    delays = [r["arm_start_gap"] for r in rows
              if r["arm_start_gap"] is not None]
    out = dict(d)
    out["rows"] = rows
    out["median_text_lag"] = statistics.median(lags) if lags else None
    out["median_arm_start_gap"] = statistics.median(delays) if delays else None
    # The rows a reader can lean on: both starts stable under the 5-article
    # threshold. Everything else is a lower bound set by one early match.
    out["robust"] = sorted(
        r["mechanism"] for r in rows
        if r["text_lag"] is not None
        and (r["text_start_fragility"] or 0) <= FRAGILE_AT)
    out["fragile"] = sorted(
        r["mechanism"] for r in rows
        if r["text_lag"] is not None
        and (r["text_start_fragility"] or 0) > FRAGILE_AT)
    rl = [r["text_lag"] for r in rows
          if r["text_lag"] is not None
          and (r["text_start_fragility"] or 0) <= FRAGILE_AT]
    out["median_robust_lag"] = statistics.median(rl) if rl else None
    out["mesh_older_than_text"] = sorted(
        r["mechanism"] for r in rows
        if r["arm_start_gap"] is not None and r["arm_start_gap"] < 0)
    out["censored_text"] = sorted(r["mechanism"] for r in rows if r["text_censored"])
    out["text_only"] = sorted(r["mechanism"] for r in rows
                              if r["text_measurable"] and not r["mesh_measurable"])
    return out


def render(d: dict) -> str:
    L = ["# From first literature to first trial\n"]
    L.append(
        "**Headline: the measurement mostly does not work, and why it does not "
        "is the useful part.** A census can in principle date a mechanism's "
        "literature and its first trial; doing so needs a vocabulary that is "
        "precise at its EARLIEST match, and neither instrument here is. The "
        "table is published with every row flagged for which failure applies "
        "to it.\n"
    )
    L.append(
        f"Generated by `scripts/census_translation_lag.py` over {d['census']:,} "
        f"census records. A mechanism STARTS in its first year with an article "
        f"and its first trial is the first year an article carries an NLM "
        f"clinical-trial publication type. Both use the SAME threshold, which "
        f"is what makes the lag a duration: a trial record is also an article "
        f"record, so at a matched threshold the first trial cannot precede the "
        f"first article. Years after {d['last_full_year']} are excluded as "
        f"incomplete.\n"
    )
    L.append(
        "Every mechanism is measured twice. The MeSH arm uses expert-assigned "
        "descriptors and is blind to anything before the descriptor existed; "
        "the text arm uses this project's keyword vocabulary over title and "
        "abstract only. **The gap between the two start years is a measurement "
        "of MeSH's own history, not of the field**, and it is reported per "
        "mechanism rather than assumed uniform.\n"
    )
    L.append("| mechanism | text start | text 1st trial | lag | fragility | "
             "MeSH start | MeSH 1st trial | lag | arm gap |")
    L.append("|---|--:|--:|--:|--:|--:|--:|--:|--:|")
    for r in sorted(d["rows"], key=lambda r: (r["text_lag"] is None,
                                              r["text_lag"] or 0)):
        def cell(v, censored=False):
            if censored:
                return "*censored*"
            return "-" if v is None else str(v)
        L.append(
            f"| {r['mechanism']}"
            f"{'' if r['mesh_measurable'] else ' †'} | "
            f"{cell(r['text_start'])} | "
            f"{cell(r['text_first_trial'], r['text_censored'])} | "
            f"{cell(r['text_lag'])} | "
            f"{cell(r['text_start_fragility'])} | "
            f"{cell(r['mesh_start'])} | "
            f"{cell(r['mesh_first_trial'], r['mesh_censored'])} | "
            f"{cell(r['mesh_lag'])} | "
            f"{cell(r['arm_start_gap'])} |")
    L.append("")
    if d["text_only"]:
        L.append(
            "† measurable on the TEXT arm only -- MeSH has no descriptor for "
            + ", ".join(f"`{m}`" for m in d["text_only"])
            + ", so these are reported as *not measurable* everywhere else in "
              "this project. Here they are measurable, on one arm, and that is "
              "marked on the row rather than left to be inferred from an empty "
              "cell.\n"
        )
    L.append("## What the two arms disagree about\n")
    L.append(
        f"The **arm gap** is the MeSH start minus the text start, and it is NOT "
        f"a measure of when a descriptor was minted -- an earlier draft called "
        f"it that and the sign refuted it. A POSITIVE gap means the text sees "
        f"the concept first, which is the descriptor-introduction effect. A "
        f"NEGATIVE gap means the DESCRIPTOR is older or broader than the term: "
        f"`Ultrasonic Therapy` runs from the 1950s while the word "
        f"\"sonodynamic\" is recent, so the descriptor arm starts decades "
        f"earlier and is counting something wider. The column mixes both "
        f"effects, which is why it is named for what it measures rather than "
        f"for what it was meant to measure. Median {d['median_arm_start_gap']} "
        f"years; read the per-mechanism column, not the median.\n"
    )
    if d["mesh_older_than_text"]:
        L.append(
            f"{len(d['mesh_older_than_text'])} mechanism(s) start EARLIER on "
            f"the descriptor arm: "
            + ", ".join(f"`{m}`" for m in d["mesh_older_than_text"])
            + ". For these the descriptor is the broader instrument and its lag "
              "is measured over a wider literature than the term names.\n"
        )
    L.append(
        "**Fragility** is how many years later the start moves if a mechanism "
        f"must reach {d['stable_min']} articles rather than 1. A large value "
        "means one early match is carrying the start year, so the lag beside "
        "it is an upper bound on the true duration.\n"
    )
    n_rob = len(d["robust"])
    n_all = n_rob + len(d["fragile"])
    L.append(
        f"**Only {n_rob} of {n_all} mechanisms support the measurement at all.** "
        f"Their median lag from first literature to first indexed trial is "
        f"{d['median_robust_lag']} years: "
        + ", ".join(f"`{m}`" for m in d["robust"]) + f". The median over all "
        f"{n_all} is {d['median_text_lag']} years and should not be quoted -- "
        f"it averages durations together with numbers set by a single early "
        f"false positive.\n"
    )
    L.append(
        f"That {n_rob}-of-{n_all} IS the result. The question is answerable in "
        f"principle from a census and is not answerable with the vocabulary "
        f"this project has, and the reason is specific enough to act on.\n"
    )
    if d["fragile"]:
        L.append(
            f"The other {len(d['fragile'])} have a start year carried by one "
            f"early match. Their lags are UPPER BOUNDS -- a false positive "
            f"sets the start EARLIER than the truth, which makes the computed "
            f"lag LONGER, so the real duration is shorter than the number "
            f"shown. (An earlier draft said lower bounds, reasoning from the "
            f"direction of the error without carrying it through the "
            f"subtraction.) The mechanisms: "
            + ", ".join(f"`{m}`" for m in d["fragile"]) + ".\n"
        )
    L.append("## Two instrument defects found and fixed, and one that remains\n")
    L.append(
        "A first version dated `car-t` to **1947** and `electrolysis` to "
        "**1950**. Both were unbounded substring matches: `car t` fires inside "
        "\"s(car t)issue\" and `echt` inside \"li(echt)enstein\". Word "
        "boundaries moved `car-t` to 2000, a 53-year correction. A first "
        "version also used a 5-article threshold for the start and a 1-trial "
        "threshold for the event, which put trials BEFORE their own literature "
        "-- TTFields showed a 2007 trial against a 2013 start. Matching the "
        "thresholds makes the ordering hold by construction, since a trial "
        "record is also an article record.\n"
    )
    L.append(
        "WHAT REMAINS IS POLYSEMY, which no boundary fixes. `electrolysis` "
        "still starts in 1952 on a paper about cosmetic hair removal and iris "
        "cysts -- a real use of the word, a different subject. `cuproptosis` "
        "starts in 1982 on `copper ionophore`, a genuine term applied to "
        "disulfiram four decades before cuproptosis was named. Those rows are "
        "wrong in a way the fragility column flags but does not repair.\n"
    )
    L.append(
        "The general point is worth more than the table. **A first-appearance "
        "statistic has no error averaging.** A prevalence estimate from an "
        "82.5%-precise vocabulary is off by a predictable fraction; a MINIMUM "
        "computed from the same vocabulary is decided entirely by its single "
        "worst false positive across four million records. The same instrument "
        "supports one statistic and not the other.\n"
    )
    if d["censored_text"]:
        L.append(
            f"{len(d['censored_text'])} mechanism(s) have a literature and no "
            f"indexed trial: " + ", ".join(f"`{m}`" for m in d["censored_text"])
            + ". These are RIGHT-CENSORED and carry no lag. Scoring them as a "
              "large lag would make a mechanism that never reached a trial "
              "indistinguishable from one that reached it slowly, and the "
              "median would silently include them.\n"
        )
    L.append("## What a lag is not\n")
    L.append(
        "It is not time to approval, and it is not evidence that a mechanism "
        "translated WELL -- an indexed trial says a trial happened, not that it "
        "worked, and a fast lag can mean a low barrier to a first-in-human "
        "study rather than a strong result. The text arm also carries this "
        "project's keyword vocabulary, whose mechanism precision is 82.5%, so "
        "an early stray match can start the clock early -- which is what the "
        "fragility column measures per mechanism rather than assuming a single "
        "threshold fixes it.\n"
    )
    return "\n".join(L)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--render-only", action="store_true")
    a = ap.parse_args()
    if a.render_only:
        d = assemble(json.loads(OUT_JSON.read_text()))
    else:
        d = assemble(scan(a.stride))
    OUT_JSON.write_text(json.dumps(d, indent=1) + "\n")
    OUT_MD.write_text(render(d))
    print(f"wrote {OUT_MD}")
    print(f"  median text lag {d['median_text_lag']}  median arm gap "
          f"{d['median_arm_start_gap']}  censored {len(d['censored_text'])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
