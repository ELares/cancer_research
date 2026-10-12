#!/usr/bin/env python3
"""Measure a candidate mechanism keyword BEFORE it enters the taxonomy.

WHY THIS EXISTS
---------------
`config.py` carries a comment recording what happens when a keyword is added
on judgement alone. The bare phrase "membrane potential" was added to the
`bioelectric` mechanism, tagged 182 articles, and was later audited at 13.3%
precision -- in this corpus the phrase almost always appears as MITOCHONDRIAL
membrane potential, a JC-1 apoptosis assay readout carried by any paper that
ran the assay and unrelated to bioelectric modulation as a therapy. Roughly 24
of the 182 were genuine. The audit happened after the fact and the correction
was a deletion.

Nothing stopped that at the time because the keyword list had no instrument.
This is the instrument, and it runs against the corpus the keyword would be
applied to rather than against intuition.

WHAT IT MEASURES, AND WHAT IT CANNOT
------------------------------------
TOPICAL ABOUTNESS, not precision, and the distinction is the whole point of
the name. For each candidate it reports the share of matching records whose
TITLE also carries the concept, on the reasoning that a paper ABOUT a
mechanism names it where it states its subject, while a paper that merely
mentions it does so once in an abstract. That is a proxy. A paper can be about
a mechanism its title does not name, and a title can name a concept the paper
only surveys, so the number bounds nothing on its own.

What decides precision is a human reading a sample, so a deterministic sample
is emitted for exactly that, following `comention_judge_prep.py`: judging is
reproducible or it is an anecdote.

THE CONTROLS ARE LOAD-BEARING
-----------------------------
A measure with no calibration point cannot be read. Two candidates are
deliberately included to fail: a bare organism word that the infection
literature owns, and a generic term whose cancer-adjacent uses are mostly
something else. If those score as highly as the specific phrases, the measure
is not discriminating and no row in the table means anything -- the same role
`random` plays in the discovery leaderboard.
"""
from __future__ import annotations

import argparse
import json
import random
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts"))
from corpus_retrieval_index import DB, connect  # noqa: E402
from corpus_retrieval_search import count, phrase_query  # noqa: E402

OUT_MD = REPO / "analysis" / "mechanism-keyword-precision.md"
OUT_JSON = REPO / "analysis" / "mechanism-keyword-precision.json"
SAMPLE_TSV = REPO / "analysis" / "mechanism-keyword-sample.tsv"

SEED = 20261011
SAMPLE_PER_CANDIDATE = 10
DEPTH = 200          # records inspected per candidate for the title share
ABOUTNESS_FLOOR = 0.30   # below this a candidate is reported as REFUSED

# The proposed mechanism and every phrase considered for it. `control: True`
# marks a candidate included BECAUSE it is expected to score badly.
CANDIDATES = [
    {"phrase": "engineered bacteria", "control": False},
    {"phrase": "bacterial therapy", "control": False},
    {"phrase": "oncolytic bacteria", "control": False},
    {"phrase": "tumor targeting bacteria", "control": False},
    {"phrase": "bacteria mediated tumor", "control": False},
    {"phrase": "attenuated salmonella", "control": False},
    {"phrase": "clostridium novyi", "control": False},
    {"phrase": "bacterial outer membrane vesicle", "control": False},
    {"phrase": "magnetotactic bacteria", "control": False},
    {"phrase": "programmable bacteria", "control": False},
    # CONTROLS, expected to score LOW. The first pair were chosen wrong and
    # the measurement said so, which is what controls are for.
    #
    # `salmonella typhimurium` and `probiotic` were picked as organism words
    # the infection and gut-health literatures own. In an UNRESTRICTED corpus
    # that would be right. This corpus is MeSH tree C04 plus an adjacent
    # extension -- cancer papers only -- so that literature was never in it,
    # and the subject restriction had already done the discriminating, leaving
    # the control nothing to pull in. The scores that first showed this are NOT
    # quoted here: they were measured over a smaller index and over the
    # rank-ordered head that `inspect_pool` now refuses, so any figure written
    # into this comment would contradict the artifact. The live numbers are in
    # the table, which is generated.
    #
    # The control a cancer-restricted corpus needs is a microbe that is
    # PROMINENT in cancer and is not a therapy. `helicobacter pylori` is the
    # best available: it CAUSES gastric cancer and carries a large literature
    # here, so a measure that cannot separate "microbe as therapy" from
    # "microbe as cause" will score it high. `16s rrna` is a sequencing METHOD
    # that appears throughout microbiome-and-cancer work for the same reason.
    {"phrase": "helicobacter pylori", "control": True},
    {"phrase": "16s rrna", "control": True},
    # Retained as the ORIGINAL, wrongly-chosen controls so the record of why
    # they do not work survives in the artifact rather than only in this
    # comment.
    {"phrase": "salmonella typhimurium", "control": True},
    {"phrase": "probiotic", "control": True},
]

# The concept a title must carry for a match to count as on-subject. Kept
# SEPARATE from the candidate phrase: testing whether the title contains the
# candidate would make every single-word candidate score 100% by construction,
# which is a guard computing its own expectation.
# A HAND-WRITTEN LIST, AND IT WAS ALREADY WRONG ONCE. The first version
# omitted `helicobacter`, so the `helicobacter pylori` control scored 0%
# microbe-in-title over 13,739 records -- not because those titles lack a
# microbe, which they obviously do not, but because the pattern did not know
# the word. A 0% that means "the detector is ignorant" is indistinguishable in
# the table from a 0% that means "the control worked", and the error ran in the
# direction that flattered the measure. `test_every_candidate_is_matchable`
# now fails if any candidate's own organism word is unreachable by this
# pattern, so the gap cannot reopen silently.
CONCEPT_RE = re.compile(
    r"\b(bacteri\w*|bacillus|salmonella|clostridi\w*|listeria|probiotic\w*|"
    r"microb\w*|microbiot\w*|helicobacter|lactobacill\w*|bifidobacteri\w*|"
    r"e\.?\s?coli|escherichia|magnetosome\w*|16s)\b", re.I)
THERAPY_RE = re.compile(
    r"\b(tumou?r\w*|cancer\w*|oncolog\w*|neoplas\w*|carcinom\w*|"
    r"immunotherap\w*|antitumou?r)\b", re.I)


def inspect_pool(c, match: str, depth: int, rng: random.Random) -> list[dict]:
    """A UNIFORM sample of the pool, deliberately not the rank-ordered head.

    THIS IS THE CIRCULARITY THE FIRST VERSION SHIPPED. It inspected
    `search(c, m, limit=DEPTH)`, which defaults to W_TITLE_HEAVY, so bm25
    selected the documents whose match fell in the TITLE -- and then
    `aboutness` measured whether the TITLE carried the concept. Because a
    candidate phrase CONTAINS its own concept word (the "bacteria" in
    "engineered bacteria" is itself a CONCEPT_RE alternative), the ranking was
    partly computing the number the measure reports.

    Measured over this index, rank-ordered top 200 against a uniform sample of
    the same pool: `helicobacter pylori` 43% against 24%, `probiotic` 27%
    against 12%, `engineered bacteria` 63% against 52%. The inflation scales
    with pool size, and the controls have the LARGEST pools, so the bias
    pushed `best_control` up and the verdict toward refusing -- the direction
    that made the refusal look earned. CONCEPT_RE was deliberately kept
    separate from the candidate phrase to avoid exactly this, and the ranking
    reintroduced it through the back door.

    Rows are sorted before sampling because the unsorted order is a function of
    how the index was built, which would make the draw depend on shard scan
    order rather than on the seed.
    """
    rows = [{"pmid": p or "", "title": t or "", "year": y or ""}
            for p, t, y in c.execute(
                "SELECT pmid, title, year FROM doc WHERE doc MATCH ?", (match,))]
    rows.sort(key=lambda r: (r["pmid"], r["title"]))
    if len(rows) <= depth:
        return rows
    return rng.sample(rows, depth)


def aboutness(rows: list[dict]) -> dict:
    """Share of matches whose TITLE names both a microbial agent and cancer.

    BOTH are required. A title naming only a microbe is as likely to be an
    infection paper as a therapy one, and this corpus is cancer-restricted at
    the record level, not at the title level, so the cancer half is not
    redundant.
    """
    n = len(rows)
    if not n:
        return {"n": 0, "microbe_in_title": 0.0, "both_in_title": 0.0}
    micro = sum(1 for r in rows if CONCEPT_RE.search(r["title"] or ""))
    both = sum(1 for r in rows
               if CONCEPT_RE.search(r["title"] or "")
               and THERAPY_RE.search(r["title"] or ""))
    return {"n": n,
            "microbe_in_title": micro / n,
            "both_in_title": both / n}


def _separation(real: list[dict], ctrl: list[dict]) -> dict:
    """Standard error of the worst-real minus best-control gap.

    Normal approximation to the difference of two binomial proportions, which
    is the right order of magnitude here and is named rather than hidden: the
    shares come from samples of at most DEPTH records, so a seven-point gap
    between a 27% estimate over 168 and a 20% estimate over 200 is under two
    standard errors and is NOT a clean separation, however the `>` reads.
    """
    if not real or not ctrl:
        return {"separation_margin": 0.0, "separation_se": 0.0,
                "separation_exceeds_noise": False}
    lo = min(real, key=lambda r: r["both_in_title"])
    hi = max(ctrl, key=lambda r: r["both_in_title"])
    margin = lo["both_in_title"] - hi["both_in_title"]
    var = 0.0
    for r in (lo, hi):
        n = max(r["inspected"], 1)
        pr = r["both_in_title"]
        var += pr * (1.0 - pr) / n
    se = var ** 0.5
    return {"separation_margin": margin,
            "separation_se": se,
            # Two standard errors, stated before the result as the bar.
            "separation_exceeds_noise": bool(se > 0 and margin > 2.0 * se),
            "separation_worst_real": lo["phrase"],
            "separation_worst_real_n": lo["inspected"],
            "separation_best_control": hi["phrase"],
            "separation_best_control_n": hi["inspected"]}


def assemble(db: Path = DB) -> dict:
    c = connect(db)
    total = c.execute("SELECT COUNT(*) FROM doc").fetchone()[0]
    if total == 0:
        raise SystemExit("index is empty -- run corpus_retrieval_index.py first")

    rows_out, sample = [], []
    for cand in CANDIDATES:
        m = phrase_query([cand["phrase"]])
        pool = count(c, m)
        # Seeded PER CANDIDATE from the phrase, so one candidate's draw does
        # not depend on how many candidates precede it in the list.
        rng = random.Random(f"{SEED}:{cand['phrase']}")
        hits = inspect_pool(c, m, DEPTH, rng)
        a = aboutness(hits)
        verdict = ("control" if cand["control"]
                   else ("admit" if a["both_in_title"] >= ABOUTNESS_FLOOR else "refuse"))
        rows_out.append({
            "phrase": cand["phrase"], "control": cand["control"],
            "pool": pool, "inspected": a["n"],
            "microbe_in_title": a["microbe_in_title"],
            "both_in_title": a["both_in_title"],
            "verdict": verdict,
            "example_titles": [(r["title"] or "")[:90] for r in hits[:3]],
        })
        # The judging sample comes out of the SAME uniform draw. Taking it from
        # a rank-ordered head would hand a human reader the most title-like
        # matches and call the resulting precision the keyword's.
        pick = list(hits)
        rng.shuffle(pick)
        for r in pick[:SAMPLE_PER_CANDIDATE]:
            sample.append({"phrase": cand["phrase"], "pmid": r["pmid"],
                           "year": r["year"], "title": r["title"], "verdict": ""})
    c.close()

    admitted = [r for r in rows_out if r["verdict"] == "admit"]
    ctrl = [r for r in rows_out if r["control"]]
    real = [r for r in rows_out if not r["control"]]
    return {
        "docs_indexed": total,
        "seed": SEED,
        "depth": DEPTH,
        # Named in the data so a reader does not have to infer it from the
        # code, and so a future change away from uniform sampling shows up in
        # the diff of the artifact rather than only in the diff of the script.
        "sampling": "uniform random draw from the whole pool, seeded per phrase",
        "aboutness_floor": ABOUTNESS_FLOOR,
        "candidates": rows_out,
        "sample": sample,
        "admitted": [r["phrase"] for r in admitted],
        "refused": [r["phrase"] for r in rows_out if r["verdict"] == "refuse"],
        # The calibration check, computed rather than asserted.
        "control_separates": (
            min((r["both_in_title"] for r in real), default=0.0)
            > max((r["both_in_title"] for r in ctrl), default=1.0)),
        "worst_real": min((r["both_in_title"] for r in real), default=0.0),
        "best_control": max((r["both_in_title"] for r in ctrl), default=0.0),
        # IS THE SEPARATION BIGGER THAN THE NOISE? `control_separates` is a
        # bare `>`, which calls a 1-record gap a win. Both shares are binomial
        # estimates from samples of at most DEPTH, so the gap has a standard
        # error and a strict inequality on its own overclaims. Reported so the
        # page cannot announce discrimination it cannot support -- the same
        # discipline `config.py` uses when it records 13.3% precision with its
        # 95% CI and n rather than alone.
        **_separation(real, ctrl),
    }


def render(d: dict) -> str:
    L: list[str] = []
    w = L.append
    w("# Measuring a mechanism keyword before adding it")
    w("")
    w(f"Over **{d['docs_indexed']:,} indexed titles and abstracts**, "
      f"up to {d['depth']} records per candidate drawn uniformly from its "
      f"pool. NOT the top {d['depth']} by rank: ranking by title weight and "
      f"then scoring whether the title carries the concept is circular, and "
      f"the first version of this page did exactly that.")
    w("")
    w("## Does the measure discriminate?")
    w("")
    if d["control_separates"] and d.get("separation_exceeds_noise"):
        w(f"**Yes.** Every real candidate scores above every control "
          f"({d['worst_real']:.0%} against {d['best_control']:.0%}), and the "
          f"{d['separation_margin']:.0%} gap is more than twice its standard "
          f"error ({d['separation_se']:.1%}), so the column below is reading "
          f"something.")
    elif d["control_separates"]:
        # The bar was set before the result: a strict `>` on two sampled shares
        # calls a one-record gap a win. This branch exists because the
        # measurement landed in it -- removing the ranking circularity moved
        # the verdict from "no" to a "yes" the sample sizes do not support.
        w(f"**Not decidably.** Every real candidate does score above every "
          f"control ({d['worst_real']:.0%} against {d['best_control']:.0%}), "
          f"so the ordering is right, but the gap is "
          f"{d['separation_margin']:.0%} against a standard error of "
          f"{d['separation_se']:.1%} -- under two SE, so it is not separable "
          f"from sampling noise.")
        w("")
        w(f"Both ends are binomial estimates over small samples: "
          f"`{d['separation_worst_real']}` at {d['worst_real']:.0%} over "
          f"{d['separation_worst_real_n']} records against "
          f"`{d['separation_best_control']}` at {d['best_control']:.0%} over "
          f"{d['separation_best_control_n']}. Reading the `>` alone would "
          f"announce a calibrated measure on a narrow gap between two noisy "
          f"numbers, which is how a measure gets trusted for the wrong reason. "
          f"Treat the verdicts below as ordered but not validated, and note "
          f"that the candidate establishing the lower bound is itself REFUSED "
          f"by the {d['aboutness_floor']:.0%} floor.")
    else:
        w(f"**No.** The worst real candidate ({d['worst_real']:.0%}) does not "
          f"beat the best control ({d['best_control']:.0%}), so this measure "
          f"is not separating a specific phrase from a generic one and **no "
          f"verdict in the table below should be acted on.**")
        w("")
        # The explanation is DERIVED from which control actually binds, so it
        # cannot go stale when the control set changes. An earlier version
        # hard-wrote the reason the FIRST pair of controls failed (an organism
        # word cannot pull in infection literature from a cancer-only corpus)
        # and kept printing it after a different control became the binding
        # one, for a different reason.
        binding = max((r for r in d["candidates"] if r["control"]),
                      key=lambda r: r["both_in_title"], default=None)
        beaten = [r["phrase"] for r in d["candidates"]
                  if not r["control"] and binding
                  and r["both_in_title"] < binding["both_in_title"]]
        if binding:
            n_real = sum(1 for r in d["candidates"] if not r["control"])
            w(f"The binding control is `{binding['phrase']}` at "
              f"{binding['both_in_title']:.0%} over {binding['pool']:,} "
              f"records. It sits above {len(beaten)} of the {n_real} real "
              f"candidates ({', '.join(f'`{b}`' for b in beaten)}), so "
              f"{n_real - len(beaten)} DO clear it.")
            w("")
            w(f"The decision rule is the strict one and was fixed in the "
              f"docstring before the result: the WORST real candidate must "
              f"beat the BEST control, because a measure that ranks a cause "
              f"above some cures is not one to admit a therapy keyword on. "
              f"Reporting only the {n_real - len(beaten)} that clear it would "
              f"be choosing the flattering half of the same table.")
            w("")
        # The two reasons name the controls they are about by READING which
        # control binds and which are merely retained, rather than hard-writing
        # phrases. An earlier version wrote "`salmonella typhimurium` and
        # `probiotic` are not controls here at all" three lines below a comment
        # boasting that the explanation was derived -- while the table beside it
        # labelled both *(control)* and both counted toward `best_control`. The
        # page contradicted itself about its own method.
        slack = [r["phrase"] for r in d["candidates"]
                 if r["control"] and binding
                 and r["both_in_title"] < binding["both_in_title"]]
        w(f"Two distinct reasons, and the second is the one that cannot be "
          f"tuned away. First, this corpus is MeSH tree C04 plus an adjacent "
          f"extension, so it holds cancer papers only: an organism word cannot "
          f"drag in the infection literature, because that literature was "
          f"never admitted. That is why "
          + (", ".join(f"`{x}`" for x in slack) or "the organism-word controls")
          + f" do not bind here -- they are retained in the table as a record "
          f"of a control choice that did not work, not because they "
          f"discriminate. Second, and decisively, a microbe that CAUSES cancer "
          f"produces titles naming both a microbe and a cancer just as reliably "
          f"as one that treats it. *Helicobacter pylori* and gastric carcinoma "
          f"is the clearest case in medicine. Title co-occurrence cannot "
          f"separate cause from cure, and this mechanism is a therapy class, so "
          f"the measure is answering a different question from the one asked.")
        w("")
        w("The consequence is acted on rather than noted: **no keyword was "
          "added to the taxonomy on this evidence.** The sample below exists "
          "because the remaining route is reading.")
    w("")
    w("## Candidates")
    w("")
    w("| phrase | pool | inspected | microbe in title | + cancer in title | verdict |")
    w("|---|--:|--:|--:|--:|---|")
    for r in d["candidates"]:
        tag = " *(control)*" if r["control"] else ""
        w(f"| `{r['phrase']}`{tag} | {r['pool']:,} | {r['inspected']:,} | "
          f"{r['microbe_in_title']:.0%} | {r['both_in_title']:.0%} | "
          f"{r['verdict']} |")
    w("")
    w(f"Inspected records are a {d['sampling']}, capped at {d['depth']}. A "
      f"share over a small pool is a coarse number: where `inspected` is well "
      f"below {d['depth']} the pool itself was that small, and one record "
      f"moves the percentage by whole points.")
    w("")
    w(f"Admitted at the {d['aboutness_floor']:.0%} floor: "
      + (", ".join(f"`{p}`" for p in d["admitted"]) or "none"))
    w("")
    w(f"Refused: " + (", ".join(f"`{p}`" for p in d["refused"]) or "none"))
    w("")
    w("## What this does not establish")
    w("")
    w("**The taxonomy is unchanged by this page either way.** A phrase marked "
      "`admit` has cleared a proxy at a stated floor, which is a reason to "
      "read the sample, not a decision to tag articles with it. Admission is a "
      "separate change with its own audit, because the keyword this instrument "
      "exists to prevent -- `membrane potential`, 13.3% precision over 182 "
      "articles -- would have cleared a floor like this one too.")
    w("")
    w("Title aboutness is a PROXY for topical precision, not a measurement of "
      "it. A paper can be about a mechanism its title does not name, and a "
      "title can name a concept the paper only surveys. Precision is decided "
      "by reading, so a deterministic sample of "
      f"{SAMPLE_PER_CANDIDATE} records per candidate is written to "
      f"`{SAMPLE_TSV.name}` with an empty verdict column; nothing here fills "
      "it in.")
    w("")
    w("Generated by `scripts/mechanism_keyword_precision.py`.")
    return "\n".join(L) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--db", default=str(DB))
    ap.add_argument("--stdout", action="store_true")
    a = ap.parse_args()
    d = assemble(Path(a.db))
    md = render(d)
    if a.stdout:
        print(md)
        return 0
    OUT_MD.write_text(md)
    OUT_JSON.write_text(json.dumps(d, indent=2) + "\n")
    with SAMPLE_TSV.open("w") as f:
        f.write("phrase\tpmid\tyear\ttitle\tverdict\n")
        for s in d["sample"]:
            t = (s["title"] or "").replace("\t", " ")
            f.write(f"{s['phrase']}\t{s['pmid']}\t{s['year']}\t{t}\t\n")
    print(f"wrote {OUT_MD.name}, {OUT_JSON.name}, {SAMPLE_TSV.name}")
    print(f"  discriminates: {d['control_separates']}  "
          f"admitted: {len(d['admitted'])}  refused: {len(d['refused'])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
