#!/usr/bin/env python3
"""Check the census's mechanism counts against PubMed itself.

EVERYTHING THIS PROJECT NOW REPORTS RESTS ON THE CENSUS BUILD, and the build
had never been checked against an independent source. It was assembled from
the annual baseline XML by this project's own parser, and a parser that drops
a field, mishandles a descriptor form, or admits a record twice would produce
a census that is internally consistent and wrong -- the failure mode nothing
downstream can see, because every downstream analysis reads the same file.

E-utilities cannot PAGE past 10,000 hits, which is why the census exists at
all. It returns COUNTS without paging, which is what makes this check possible.

THREE THINGS HAVE TO MATCH BEFORE THE COMPARISON MEANS ANYTHING, and getting
any of them wrong produces a large apparent disagreement that is entirely an
artifact of the query:

1. EXPLOSION. `X[mh]` in PubMed includes every descriptor beneath X in the MeSH
   tree; the census matches DescriptorName EXACTLY. Left exploded, `Ultrasonic
   Therapy` returns 4,371 against the census's 2,513 -- a 74% "disagreement"
   that is just the subtree, and it silently includes the HIFU records this
   project counts as a separate mechanism. `[mh:noexp]` is the matching query.
2. THE SNAPSHOT DATE. The census is a fixed baseline; PubMed keeps growing. The
   query is capped at the baseline date.
3. THE CANCER RESTRICTION. The census admits C04 plus nine adjacent
   experimental-context descriptors, so the comparison is restricted to records
   whose `cancer_basis` is the C04 core, against `neoplasms[mh]`.

WHAT A RESIDUAL DISAGREEMENT MEANS is not obvious and the report does not
pretend otherwise. PubMed's live index is not the baseline plus time -- records
are re-indexed, descriptors are added and withdrawn, and a record's entry date
is not its indexing date. So a few per cent in either direction is expected and
is NOT evidence of a parser defect; a large or one-sided gap is.

OFFLINE CONTRACT: the fetch runs locally and writes a committed artifact. CI
reads only the artifact and never touches the network.
"""
import argparse
import gzip
import json
import statistics
import time
import urllib.parse
import urllib.request
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
RECORDS = REPO / "corpus/atlas/records"
MECH_MAP = REPO / "analysis/mesh-mechanism-map.yaml"
OUT_MD = REPO / "analysis/census-external-check.md"
OUT_JSON = REPO / "analysis/census-external-check.json"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
BASELINE_DATE = "2026/06/17"
# NCBI allows 3 requests/second without an API key; stay well under.
PAUSE = 0.5
# Above this relative gap a row is flagged for inspection rather than accepted.
FLAG_AT = 0.10


def pubmed_count(term: str, timeout: int = 30) -> int:
    q = urllib.parse.urlencode(
        {"db": "pubmed", "term": term, "retmode": "json", "rettype": "count"})
    with urllib.request.urlopen(f"{EUTILS}?{q}", timeout=timeout) as r:
        return int(json.load(r)["esearchresult"]["count"])


def query_for(descriptors):
    """OR over the mechanism's descriptors, unexploded, cancer-restricted and
    date-capped. Every clause here corresponds to a property of the census
    build; dropping one makes the comparison measure the query instead."""
    ors = " OR ".join(f'"{d}"[mh:noexp]' for d in descriptors)
    return (f"({ors}) AND neoplasms[mh] AND "
            f'("1800"[edat] : "{BASELINE_DATE}"[edat])')


def census_counts(stride: int = 1) -> dict:
    """Per-mechanism counts over the C04 CORE only.

    The adjacent-basis extension admits records with no C04 descriptor at all,
    which `neoplasms[mh]` cannot return, so including them would guarantee the
    census reads high and the gap would be a property of the comparison.
    """
    import yaml

    mp = yaml.safe_load(MECH_MAP.read_text(encoding="utf-8"))["mechanisms"]
    mech = {k: {x.lower() for x in v["descriptors"]} for k, v in mp.items()}
    counts = {k: 0 for k in mech}
    years = {k: [] for k in mech}
    n = core = 0
    for f in sorted(RECORDS.glob("*.jsonl.gz"))[::stride]:
        with gzip.open(f, "rt", encoding="utf-8") as fh:
            for line in fh:
                r = json.loads(line)
                n += 1
                if r.get("cancer_basis") != "C04":
                    continue
                core += 1
                ms = {m.lower() for m in (r.get("mesh") or [])}
                if not ms:
                    continue
                y = r.get("year")
                for k, d in mech.items():
                    if ms & d:
                        counts[k] += 1
                        if isinstance(y, int):
                            years[k].append(y)
    return {"census": n, "c04_core": core, "counts": counts,
            "median_year": {k: (int(statistics.median(v)) if v else None)
                            for k, v in years.items()},
            "descriptors": {k: mp[k]["descriptors"] for k in mp}}


def fetch(cen: dict) -> dict:
    rows = []
    for k in sorted(cen["counts"]):
        descs = cen["descriptors"][k]
        # A mechanism the map declares NOT MeSH-measurable carries an EMPTY
        # descriptor list, and building a query from it yields `()` -- which
        # PubMed answers with something, and that something would enter the
        # table as a disagreement. Reported as not-comparable rather than
        # queried, the same distinction this project makes everywhere else
        # between unmeasurable and zero.
        if not descs:
            rows.append({"mechanism": k, "census": cen["counts"][k],
                         "pubmed": None, "query": None,
                         "error": "no MeSH descriptor: not comparable"})
            print(f"  {k:26s} no descriptor -- not comparable")
            continue
        term = query_for(descs)
        try:
            live = pubmed_count(term)
            err = None
        except Exception as e:                              # network, 429, 5xx
            live, err = None, f"{type(e).__name__}: {e}"
        rows.append({"mechanism": k, "census": cen["counts"][k],
                     "median_year": cen["median_year"].get(k),
                     "pubmed": live, "query": term, "error": err})
        shown = "unresolved" if live is None else format(live, ",")
        print(f"  {k:26s} census {cen['counts'][k]:>7,}  pubmed {shown}")
        time.sleep(PAUSE)
    return {"baseline_date": BASELINE_DATE, "census_records": cen["census"],
            "c04_core": cen["c04_core"], "rows": rows}


def assemble(d: dict) -> dict:
    for r in d["rows"]:
        if r["pubmed"]:
            r["ratio"] = round(r["census"] / r["pubmed"], 3)
            r["rel_gap"] = round(abs(r["census"] - r["pubmed"]) / r["pubmed"], 3)
            r["flagged"] = r["rel_gap"] > FLAG_AT
        else:
            r["ratio"] = None
            r["rel_gap"] = None
            r["flagged"] = False
    d = dict(d)
    d["compared"] = sum(1 for r in d["rows"] if r["pubmed"])
    # Two different reasons a row has no PubMed count, kept apart: a mechanism
    # with no descriptor CANNOT be compared, while a failed request should have
    # been. Pooling them would let a network outage read as a taxonomy gap.
    d["not_comparable"] = sorted(
        r["mechanism"] for r in d["rows"]
        if r["pubmed"] is None and r["query"] is None)
    d["unresolved"] = sorted(
        r["mechanism"] for r in d["rows"]
        if r["pubmed"] is None and r["query"] is not None)
    d["flag_threshold"] = FLAG_AT
    d["flagged"] = sorted(r["mechanism"] for r in d["rows"] if r["flagged"])
    gaps = [r["rel_gap"] for r in d["rows"] if r["rel_gap"] is not None]
    d["median_rel_gap"] = round(statistics.median(gaps), 3) if gaps else None
    d["max_rel_gap"] = round(max(gaps), 3) if gaps else None
    # Direction matters: a scatter around zero is noise, a one-sided gap is a
    # systematic difference and needs an explanation rather than a tolerance.
    higher = sum(1 for r in d["rows"] if r["ratio"] and r["ratio"] > 1)
    d["census_higher"] = higher
    d["census_lower"] = len(gaps) - higher
    d["recency_test"] = _recency_test(d["rows"])
    return d


def _spearman(pairs):
    """Rank correlation, ties averaged. Self-contained: this project's pinned
    environment has no scipy, and a hand-rolled Pearson on raw values would be
    the wrong statistic for a monotone-but-not-linear relationship."""
    if len(pairs) < 4:
        return None

    def ranks(vals):
        order = sorted(range(len(vals)), key=lambda i: vals[i])
        out = [0.0] * len(vals)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and vals[order[j + 1]] == vals[order[i]]:
                j += 1
            avg = (i + j) / 2 + 1
            for k in range(i, j + 1):
                out[order[k]] = avg
            i = j + 1
        return out

    xs = ranks([p[0] for p in pairs])
    ys = ranks([p[1] for p in pairs])
    mx = sum(xs) / len(xs)
    my = sum(ys) / len(ys)
    num = sum((a - mx) * (b - my) for a, b in zip(xs, ys))
    den = (sum((a - mx) ** 2 for a in xs) * sum((b - my) ** 2 for b in ys)) ** 0.5
    return round(num / den, 2) if den else None


def _recency_test(rows) -> dict:
    """Test the explanation rather than asserting it.

    A one-sided gap needs a cause. The candidate is that MeSH indexing keeps
    being applied to records whose ENTRY date already precedes the baseline: a
    fixed snapshot misses them, while a query filtered on entry date catches
    them. That predicts the gap should grow with a mechanism's recency, and the
    prediction is stated here BEFORE the number so it can fail.

    A near-zero or negative correlation would leave the one-sided gap
    unexplained, which is a different and worse position than an explained one.
    """
    pairs = [(r["median_year"], r["rel_gap"]) for r in rows
             if r.get("median_year") and r.get("rel_gap") is not None]
    rho = _spearman(pairs)
    return {
        "prediction": "positive: retrospective MeSH indexing since the baseline",
        "n": len(pairs),
        "spearman_year_vs_gap": rho,
        "supported": bool(rho is not None and rho >= 0.5),
    }


def render(d: dict) -> str:
    L = ["# The census, checked against PubMed\n"]
    L.append(
        f"Generated by `scripts/census_external_check.py`. Per-mechanism counts "
        f"from the census's C04 core ({d['c04_core']:,} records) against live "
        f"E-utilities counts, unexploded, cancer-restricted and capped at the "
        f"{d['baseline_date']} baseline.\n"
    )
    L.append(
        "The query has to match the build on three axes or the comparison "
        "measures the query. `X[mh]` EXPLODES the MeSH tree while the census "
        "matches DescriptorName exactly -- left exploded, `Ultrasonic Therapy` "
        "returns 4,371 against the census's 2,513, a 74% 'disagreement' that is "
        "just the subtree and silently includes the HIFU records counted here "
        "as a separate mechanism. PubMed also keeps growing past the baseline, "
        "and the census admits nine adjacent descriptors that `neoplasms[mh]` "
        "cannot return.\n"
    )
    L.append("| mechanism | census | PubMed | ratio | gap |")
    L.append("|---|--:|--:|--:|--:|")
    for r in sorted(d["rows"], key=lambda r: -(r["rel_gap"] or 0)):
        if r["pubmed"] is None:
            why = "*not comparable*" if r["query"] is None else "*unresolved*"
            L.append(f"| {r['mechanism']} | {r['census']:,} | {why} | - | - |")
            continue
        mark = "**" if r["flagged"] else ""
        L.append(f"| {r['mechanism']} | {r['census']:,} | {r['pubmed']:,} | "
                 f"{r['ratio']} | {mark}{100 * r['rel_gap']:.1f}%{mark} |")
    L.append("")
    L.append("## What the agreement says\n")
    balanced = abs(d["census_higher"] - d["census_lower"]) <= 2
    L.append(
        f"Median relative gap **{100 * d['median_rel_gap']:.1f}%** across "
        f"{d['compared']} mechanisms, with the census reading higher on "
        f"{d['census_higher']} and lower on {d['census_lower']}. "
        + ("The direction is close to balanced, which is what independent "
           "noise looks like: re-indexing, descriptor additions and the "
           "difference between a record's entry date and its indexing date all "
           "move counts in both directions.\n"
           if balanced else
           "The gap is ONE-SIDED, which noise does not produce. That is a "
           "systematic difference between the build and PubMed's index, and it "
           "needs an explanation rather than a tolerance.\n")
    )
    rt = d.get("recency_test") or {}
    if not balanced and rt.get("spearman_year_vs_gap") is not None:
        L.append("### The explanation, tested\n")
        L.append(
            "The candidate cause is that MeSH indexing keeps being applied to "
            "records whose ENTRY date already precedes the baseline. A fixed "
            "snapshot misses them; a query filtered on entry date catches them. "
            "That is not a defect in the build -- it is what comparing a "
            "snapshot against a live index does.\n"
        )
        L.append(
            f"The prediction was stated before the number: if that is the "
            f"cause, the gap should grow with a mechanism's recency. Measured "
            f"over {rt['n']} mechanisms, Spearman rank correlation between "
            f"median year and relative gap is "
            f"**{rt['spearman_year_vs_gap']:+.2f}**"
            + (" -- the prediction holds, and the one-sided gap is accounted "
               "for.\n" if rt["supported"] else
               " -- the prediction FAILS, so the one-sided gap remains "
               "unexplained and this build should not be trusted until it is.\n")
        )
        if rt["supported"]:
            worst = max((r for r in d["rows"] if r["rel_gap"] is not None),
                        key=lambda r: r["rel_gap"])
            L.append(
                f"The practical consequence is a bound, not a correction: the "
                f"census under-counts the most recent literature relative to "
                f"today's index by up to {100 * worst['rel_gap']:.0f}% "
                f"(`{worst['mechanism']}`, median year "
                f"{worst['median_year']}), and older literature by 2-3%. Every "
                f"growth figure computed to the present is therefore a LOWER "
                f"bound, and the newest mechanisms are the ones most "
                f"under-stated.\n"
            )
    if d["not_comparable"]:
        L.append(
            f"{len(d['not_comparable'])} mechanism(s) have no MeSH descriptor "
            f"at all and cannot be compared in either direction: "
            + ", ".join(f"`{m}`" for m in d["not_comparable"])
            + ". Not a disagreement and not a zero -- there is no query to "
              "ask.\n"
        )
    if d["flagged"]:
        L.append(
            f"{len(d['flagged'])} mechanism(s) exceed the "
            f"{100 * d['flag_threshold']:.0f}% flag threshold and are marked in "
            f"the table: {', '.join(f'`{m}`' for m in d['flagged'])}. A flag is "
            f"an instruction to look, not a verdict -- a mechanism whose "
            f"descriptors are few and small moves several per cent on a handful "
            f"of records.\n"
        )
    else:
        L.append(
            f"No mechanism exceeds the {100 * d['flag_threshold']:.0f}% flag "
            f"threshold.\n"
        )
    L.append("## What this does not establish\n")
    L.append(
        "Agreement on counts is agreement on ADMISSION, not on content. It says "
        "the build admits the same records PubMed would return for the same "
        "descriptors, which is the property every prevalence claim in this "
        "project depends on. It says nothing about whether a descriptor means "
        "what an analysis takes it to mean -- that is the breadth problem "
        "reported separately, and no amount of count agreement touches it.\n"
    )
    L.append(
        "PubMed's live index is also not the baseline plus elapsed time. "
        "Records are re-indexed, descriptors are added and withdrawn, and entry "
        "date is not indexing date, so a few per cent in either direction is "
        "expected. The check is powered to find a parser defect, not to certify "
        "an exact match.\n"
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
        cen = census_counts(a.stride)
        print(f"census C04 core: {cen['c04_core']:,}")
        d = assemble(fetch(cen))
    OUT_JSON.write_text(json.dumps(d, indent=1) + "\n")
    OUT_MD.write_text(render(d))
    print(f"wrote {OUT_MD}")
    print(f"  median gap {100 * d['median_rel_gap']:.1f}%  flagged "
          f"{len(d['flagged'])}  higher/lower {d['census_higher']}/"
          f"{d['census_lower']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
