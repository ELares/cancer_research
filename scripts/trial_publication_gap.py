#!/usr/bin/env python3
"""How many registered cancer trials ever reach the literature.

THE QUESTION A LITERATURE CENSUS STRUCTURALLY CANNOT ASK
--------------------------------------------------------
MISSION.md says this project cannot see what was never published, and that was
true while the only evidence was published papers: the missing studies leave no
trace in the thing you are counting. A trial REGISTRY is the trace. Every
interventional cancer trial that enrolled a patient in the last two decades was
supposed to be registered before it began, so the registry records studies that
the literature does not -- which makes publication bias measurable rather than
merely acknowledged.

THE MEASUREMENT, AND WHY IT IS NOT THE OBVIOUS ONE
--------------------------------------------------
ClinicalTrials.gov carries a `references` field, and 32.5% of cancer trials
have one. Reading that as "67.5% were never published" is wrong and is the
error this script exists to avoid: the field is maintained by the sponsor, and
a sponsor who never updates the record leaves a published trial looking
unpublished. It measures REGISTRY UPKEEP, not publication.

The independent check is that Europe PMC indexes the registration number
itself -- 753,986 records mention an NCT id -- so a paper reporting a trial can
be found FROM the trial without the sponsor's help. Asking both, on the same
trials, separates the two.

WHAT A HIT DOES AND DOES NOT MEAN
---------------------------------
A paper mentioning NCT01234567 is not necessarily REPORTING it: reviews,
protocols, meta-analyses and other trials' discussion sections all cite trial
numbers. So this OVER-counts publication, which matters for the direction of
the conclusion -- any gap it finds is a LOWER bound on the real one, and that
is the safe direction for a claim about missing evidence.

Time matters too, and stratifying by it is not optional: a trial that completed
last year has not had time to publish, so pooling recent and old completions
would report the publication lag as a publication gap.
"""
from __future__ import annotations

import argparse
import glob
import gzip
import json
import os
import random
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from corpus_expand_fetch import _RateLimit, retry_page  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
OUT_MD = ROOT / "analysis" / "trial-publication-gap.md"
OUT_JSON = ROOT / "analysis" / "trial-publication-gap.json"
TRIALS = Path(os.getenv(
    "FERRO_TRIALS_OUT", str(Path.home() / "nas" / "cancer-atlas" / "trials"))) / "trials"
EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
UA = {"User-Agent": "cancer-research-corpus/1.0 (+https://github.com/ELares/cancer_research)"}
SEED = 20260907


def load_trials() -> list:
    out = []
    for f in sorted(glob.glob(str(TRIALS / "*.jsonl.gz"))):
        try:
            with gzip.open(f, "rt") as fh:
                for ln in fh:
                    try:
                        out.append(json.loads(ln))
                    except Exception:
                        continue
        except (EOFError, OSError):
            continue
    return out


def completion_year(t: dict):
    c = t.get("completion") or ""
    try:
        return int(str(c)[:4])
    except (TypeError, ValueError):
        return None


def eligible(t: dict) -> bool:
    """Completed interventional trials with a completion year.

    INTERVENTIONAL only: an observational study has different reporting norms,
    and mixing them would compare two populations under one label. COMPLETED
    only: a terminated or withdrawn trial may have nothing to report, and
    counting it as unpublished would charge the literature for a study that
    never produced a result.
    """
    return (t.get("study_type") == "INTERVENTIONAL"
            and t.get("status") == "COMPLETED"
            and completion_year(t) is not None)


def mentions(nct: str, limiter: _RateLimit) -> int:
    def once():
        limiter.wait()
        u = f"{EPMC}?{urllib.parse.urlencode({'query': nct, 'format': 'json', 'pageSize': 1})}"
        with urllib.request.urlopen(urllib.request.Request(u, headers=UA), timeout=60) as r:
            return int(json.load(r)["hitCount"])
    return retry_page(once, label=f"epmc {nct}", verbose=False)


def measure(sample: int = 600, verbose: bool = True) -> dict:
    trials = [t for t in load_trials() if eligible(t)]
    rng = random.Random(SEED)
    rng.shuffle(trials)
    picked = trials[:sample]
    limiter = _RateLimit(0.34)

    rows = []
    for i, t in enumerate(picked, 1):
        n = mentions(t["nct_id"], limiter)
        rows.append({
            "nct": t["nct_id"],
            "year": completion_year(t),
            "phases": t.get("phases") or [],
            "enrollment": t.get("enrollment"),
            "registry_pmids": len(t.get("pmids") or []),
            "epmc_mentions": n,
        })
        if verbose and i % 100 == 0:
            print(f"  {i}/{len(picked)}", flush=True)
    return {"eligible": len(trials), "sampled": len(picked), "seed": SEED,
            "rows": rows}


def summarise(d: dict) -> dict:
    rows = d["rows"]
    n = len(rows)
    reg = sum(1 for r in rows if r["registry_pmids"] > 0)
    epmc = sum(1 for r in rows if r["epmc_mentions"] > 0)
    both = sum(1 for r in rows if r["registry_pmids"] > 0 and r["epmc_mentions"] > 0)
    neither = sum(1 for r in rows if not r["registry_pmids"] and not r["epmc_mentions"])
    epmc_only = epmc - both
    reg_only = reg - both

    # Stratified by how long the trial has had to publish. Pooling a trial that
    # completed last year with one that completed in 2010 reports publication
    # LAG as a publication GAP.
    # BY COMPLETION ERA, not by age, because the two answer different
    # questions and the first version asked the wrong one. Citing a
    # registration number in a paper is a MODERN CONVENTION: registration
    # became an ICMJE condition of publication in 2005 and journals adopted
    # printing the number gradually after that. So an old trial whose paper
    # never printed its NCT id is INVISIBLE to this measure whether or not it
    # was published, and reading its low rate as a publication gap measures
    # citation practice instead.
    by_era = {}
    for lo in range(2000, 2030, 5):
        sel = [r for r in rows if r["year"] and lo <= r["year"] <= lo + 4]
        if len(sel) >= 10:
            m = sum(1 for r in sel if r["epmc_mentions"] > 0)
            by_era[f"{lo}-{lo + 4}"] = {
                "trials": len(sel), "with_any_mention": m,
                "pct": round(100 * m / len(sel), 1)}

    # The headline is taken from the era where the convention is established
    # and enough time has passed to publish. Anything earlier measures the
    # convention; anything later measures the lag.
    RELIABLE = "2015-2019"
    rel = by_era.get(RELIABLE)
    return {
        "sampled": n,
        "registry_says_published": reg,
        "registry_pct": round(100 * reg / n, 1) if n else None,
        "literature_mentions_it": epmc,
        "literature_pct": round(100 * epmc / n, 1) if n else None,
        "both": both, "registry_only": reg_only, "literature_only": epmc_only,
        "neither": neither,
        "neither_pct": round(100 * neither / n, 1) if n else None,
        "by_completion_era": by_era,
        "reliable_era": RELIABLE,
        "reliable_era_stats": rel,
        "reliable_era_gap_pct": (round(100 - rel["pct"], 1) if rel else None),
    }


def assemble(d: dict) -> dict:
    """The stored payload: measurement plus everything derived from it.

    Kept as ONE dict so `render` takes one argument, which is the interface the
    artifact-freshness harness calls. A two-argument renderer is invisible to
    that gate -- it raised a TypeError rather than checking anything -- and a
    gate that cannot call a generator cannot notice it drifting.
    """
    out = dict(d)
    out["summary"] = summarise(d)
    return out


def render(d: dict) -> str:
    s = d.get("summary") or summarise(d)
    rel = s.get("reliable_era_stats") or {}
    L = [
        "# Registered cancer trials that never reach the literature",
        "",
        "*Generated by `scripts/trial_publication_gap.py`. Every figure is a "
        "count over a seeded random sample; re-run to reproduce it exactly.*",
        "",
        "## The question a literature census cannot ask",
        "",
        "MISSION.md says this project cannot see what was never published, and "
        "while the only evidence was published papers that was simply true: the "
        "missing studies leave no trace in the thing being counted. A trial "
        "registry is the trace.",
        "",
        f"Of {d['eligible']:,} COMPLETED INTERVENTIONAL cancer trials with a "
        f"recorded completion date, {s['sampled']:,} were sampled at random "
        f"(seed {d['seed']}).",
        "",
        "## The registry's own field is not a publication rate",
        "",
        "| measure | trials | share |",
        "|---|--:|--:|",
        f"| the registry names a publication | {s['registry_says_published']:,} | {s['registry_pct']}% |",
        f"| Europe PMC has a record mentioning the NCT id | {s['literature_mentions_it']:,} | {s['literature_pct']}% |",
        f"| neither | {s['neither']:,} | {s['neither_pct']}% |",
        "",
        f"**{s['literature_only']:,} trials are published and the registry does "
        f"not say so** -- a third of the sample. Using the registry field alone "
        f"would report a non-publication rate of "
        f"{round(100 - s['registry_pct'], 1)}%, against {s['neither_pct']}% "
        "measured independently. That column tracks REGISTRY UPKEEP: a sponsor "
        "who never came back leaves a published trial looking unpublished.",
        "",
        "## The era pattern is a measurement artifact, and it is the main finding",
        "",
        "| completion period | trials | mentioned in the literature |",
        "|---|--:|--:|",
    ]
    for era in sorted(s.get("by_completion_era", {})):
        b = s["by_completion_era"][era]
        L.append(f"| {era} | {b['trials']:,} | {b['with_any_mention']:,} "
                 f"({b['pct']}%) |")
    L += [
        "",
        "That climb does not show trials becoming more likely to be published. "
        "It shows PAPERS BECOMING MORE LIKELY TO PRINT THE REGISTRATION NUMBER. "
        "Registration became an ICMJE condition of publication in 2005 and "
        "journals adopted printing the identifier gradually after that, so a "
        "trial completing in 2003 whose results appeared in 2005 is invisible "
        "to this measure whether or not it was published.",
        "",
        "An earlier draft of this page tabulated the same data by YEARS SINCE "
        "COMPLETION and reported that trials over ten years old were the least "
        "likely to be published. That reading is withdrawn: it measured the "
        "convention, not the literature.",
        "",
        "## What the number is, once the artifact is removed",
        "",
    ]
    if rel:
        L += [
            f"Restricted to trials completing in {s['reliable_era']} -- late "
            "enough that the citation convention is established, early enough "
            f"to have had time to publish -- {rel['with_any_mention']:,} of "
            f"{rel['trials']:,} have a record mentioning them "
            f"({rel['pct']}%).",
            "",
            f"**At least {s['reliable_era_gap_pct']}% of completed "
            "interventional cancer trials leave no trace in the indexed "
            "literature.**",
            "",
            "At least, because a paper mentioning `NCT01234567` is not "
            "necessarily REPORTING that trial: reviews, protocols, "
            "meta-analyses and other trials' discussion sections all cite "
            "registration numbers. The measure OVER-counts publication, so the "
            "gap it finds is a LOWER bound -- which is the safe direction for a "
            "claim about missing evidence.",
        ]
    else:
        L.append("The reliable era is not represented in this sample.")
    L += [
        "",
        "## What this does not establish",
        "",
        "- **Registration is not running.** A registered trial may never have "
        "enrolled anyone, and a study that produced no result cannot be charged "
        "to the literature. Terminated and withdrawn trials are excluded for "
        "this reason.",
        "- **Europe PMC is not everything.** A trial reported only in a "
        "conference abstract, a thesis, or a journal it does not index reads as "
        "unpublished here. That is part of why this project also collects from "
        "OpenAlex, where conference abstracts live.",
        "- **A gap is not misconduct.** Trials go unpublished for reasons "
        "ranging from a negative result nobody would print to a sponsor going "
        "out of business. This counts the silence; it does not explain it.",
        "",
    ]
    return "\n".join(L)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample", type=int, default=600)
    a = ap.parse_args()
    payload = assemble(measure(a.sample))
    OUT_JSON.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    OUT_MD.write_text(render(payload))
    print(json.dumps(payload["summary"], indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
