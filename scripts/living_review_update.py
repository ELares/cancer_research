#!/usr/bin/env python3
"""Living-review (PRISMA-LSR) incremental corpus update (#349).

The manuscript corpus is a FROZEN snapshot. This re-runs the committed mechanism
queries (`scripts/queries.txt`) against PubMed for a recent window, diffs the
results against the frozen corpus, and emits (1) a dated incremental index of the
NEW records and (2) a delta changelog (new articles per query, new landmark
detections via the mechanism tagger), WITHOUT mutating the frozen corpus. A
scheduled GitHub Action (`.github/workflows/living-review.yml`) runs this so the
consolidation does not go stale, per PRISMA for Living Systematic Reviews
(PRISMA-LSR, Brignardello-Petersen et al.).

Separation of concerns (load-bearing): the frozen manuscript corpus lives under
`corpus/by-pmid/`, `corpus/abstracts/by-pmid/`, and `corpus/INDEX.jsonl`; the
living review writes ONLY under `corpus/living/<date>/` and
`analysis/living-review/<date>.md`. It never touches the frozen files, so the
manuscript's numbers are reproducible and the living index is clearly an
addendum.

The PUBMED-touching steps are isolated from the PURE delta/changelog logic so the
latter is unit-testable without the network (see `tests/test_living_review.py`).

Usage:
  python scripts/living_review_update.py --since 2025-01-01          # full run
  python scripts/living_review_update.py --since 2025-01-01 --dry-run  # counts only
"""

import argparse
import json
import sys
from datetime import date, datetime, timezone
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts"))

QUERIES_FILE = REPO / "scripts" / "queries.txt"
LIVING_DIR = REPO / "corpus" / "living"
REPORT_DIR = REPO / "analysis" / "living-review"
FROZEN_INDEX = REPO / "corpus" / "INDEX.jsonl"
FROZEN_PMID_DIRS = [
    REPO / "corpus" / "by-pmid",
    REPO / "corpus" / "abstracts" / "by-pmid",
]


# --------------------------------------------------------------------------- #
# Pure logic (no network, no filesystem writes) — unit-tested.
# --------------------------------------------------------------------------- #
def load_queries(text):
    """Parse a queries.txt-format string into a list of query strings (one per
    non-blank, non-comment line)."""
    return [
        line.strip()
        for line in text.splitlines()
        if line.strip() and not line.strip().startswith("#")
    ]


def windowed_query(query, since_iso):
    """Append a PubMed publication-date lower bound to a query so esearch returns
    only records published on/after `since_iso` (YYYY-MM-DD). Open upper bound."""
    y, m, d = since_iso.split("-")
    lo = f"{y}/{m}/{d}"
    return f'({query}) AND ("{lo}"[Date - Publication] : "3000"[Date - Publication])'


def compute_delta(found_by_query, frozen_pmids):
    """Given {query: [pmids found]} and the set of frozen-corpus PMIDs, return
    (new_by_query, all_new) where new_by_query is {query: [pmids not in frozen]}
    (order-preserving, de-duplicated within a query) and all_new is the ordered
    de-duplicated union across all queries (a PMID counts as new only once)."""
    new_by_query = {}
    all_new = []
    seen = set()
    for query, pmids in found_by_query.items():
        fresh = []
        local_seen = set()
        for p in pmids:
            if p in frozen_pmids or p in local_seen:
                continue
            local_seen.add(p)
            fresh.append(p)
            if p not in seen:
                seen.add(p)
                all_new.append(p)
        new_by_query[query] = fresh
    return new_by_query, all_new


_LANDMARK_PUBTYPES = (
    "clinical trial, phase iii",
    "clinical trial, phase iv",
    "randomized controlled trial",
    "practice guideline",
    "meta-analysis",
)


def is_landmark(record):
    """Heuristic 'landmark' flag for the changelog: a high-evidence PubMed
    publication type (phase III/IV trial, RCT, guideline, meta-analysis) or an
    unusually highly-cited new record (when a consumer has enriched
    `cited_by_count`). Records are metadata dicts (pmid, title, year, pub_types,
    mechanisms, ...)."""
    pub_types = [str(t).lower() for t in (record.get("pub_types") or [])]
    if any(pt in _LANDMARK_PUBTYPES for pt in pub_types):
        return True
    if (record.get("cited_by_count") or 0) >= 100:
        return True
    return False


def format_changelog(run_date, since_iso, new_by_query, records, n_frozen):
    """Render the PRISMA-LSR delta changelog (markdown). `records` is the list of
    tagged metadata dicts for the new PMIDs. Pure string assembly."""
    total_new = len(records)
    by_pmid = {r["pmid"]: r for r in records}
    lines = []
    lines.append(f"# Living-review update — {run_date}\n")
    lines.append(
        f"PRISMA-LSR incremental update (#349). Re-ran the {len(new_by_query)} "
        f"committed mechanism queries against PubMed for records published since "
        f"`{since_iso}`, diffed against the frozen corpus ({n_frozen:,} PMIDs). "
        f"**The frozen manuscript corpus is unchanged**; the new records below live "
        f"only under `corpus/living/{run_date}/`.\n"
    )
    lines.append(f"**New records this window: {total_new}.**\n")

    landmarks = [r for r in records if is_landmark(r)]
    if landmarks:
        lines.append(f"## New landmark detections ({len(landmarks)})\n")
        for r in sorted(landmarks, key=lambda r: -(r.get("cited_by_count") or 0)):
            mech = ", ".join(r.get("mechanisms") or []) or "untagged"
            lines.append(
                f"- PMID {r['pmid']} ({r.get('year', '?')}) — {r.get('title', '')[:120]} "
                f"[{mech}; {r.get('evidence_level', 'n/a')}; cites {r.get('cited_by_count', 0)}]"
            )
        lines.append("")

    lines.append("## New articles per query\n")
    lines.append("| query | new |")
    lines.append("|---|--:|")
    for query, pmids in new_by_query.items():
        short = query if len(query) <= 70 else query[:67] + "..."
        lines.append(f"| {short} | {len(pmids)} |")
    lines.append("")

    if records:
        lines.append("## New records\n")
        lines.append("| pmid | year | mechanisms | title |")
        lines.append("|---|--:|---|---|")
        for r in sorted(records, key=lambda r: -int(r.get("year") or 0)):
            mech = ", ".join(r.get("mechanisms") or []) or "—"
            lines.append(
                f"| {r['pmid']} | {r.get('year', '?')} | {mech} | {r.get('title', '')[:90]} |"
            )
        lines.append("")
    _ = by_pmid  # reserved for future evidence-shift diffing
    return "\n".join(lines) + "\n"


# --------------------------------------------------------------------------- #
# IO layer (network + filesystem) — thin wrappers over the existing pipeline.
# --------------------------------------------------------------------------- #
def load_frozen_pmids():
    """Every PMID already in the frozen corpus (full-text + abstract dirs + the
    committed INDEX), so the delta never re-ingests an existing record."""
    pmids = set()
    for d in FROZEN_PMID_DIRS:
        if d.exists():
            pmids.update(p.stem for p in d.glob("*.md"))
    if FROZEN_INDEX.exists():
        for line in FROZEN_INDEX.read_text().splitlines():
            line = line.strip()
            if line:
                try:
                    pmids.add(str(json.loads(line)["pmid"]))
                except (json.JSONDecodeError, KeyError):
                    continue
    return pmids


def tag_record(meta):
    """Attach mechanism tags to a metadata dict via the committed tagger keywords,
    so the changelog reports new articles per mechanism with the same definitions
    the frozen corpus used."""
    from config import MECHANISM_KEYWORDS
    from tag_articles import get_searchable_text, has_cancer_context, match_keywords

    # Reconstruct the minimal markdown body the tagger reads the abstract from
    # (`get_searchable_text` extracts it from a `## Abstract` section).
    body = f"## Abstract\n\n{meta.get('abstract', '')}\n"
    text = get_searchable_text(meta, body, include_full_text=False)
    mechanisms = match_keywords(text, MECHANISM_KEYWORDS) if has_cancer_context(text) else []
    meta = dict(meta)
    meta["mechanisms"] = mechanisms
    return meta


def main():
    parser = argparse.ArgumentParser(description="Living-review incremental update (#349)")
    parser.add_argument(
        "--since", required=True, help="Publication-date lower bound (YYYY-MM-DD)"
    )
    parser.add_argument("--max", type=int, default=200, help="Max PMIDs per query")
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Query + diff + changelog only; do not fetch metadata or write the index",
    )
    args = parser.parse_args()
    datetime.strptime(args.since, "%Y-%m-%d")  # validate
    run_date = datetime.now(timezone.utc).date().isoformat()

    from fetch_articles import fetch_pubmed_metadata, pubmed_search

    queries = load_queries(QUERIES_FILE.read_text())
    frozen = load_frozen_pmids()
    print(f"Frozen corpus: {len(frozen):,} PMIDs; {len(queries)} queries since {args.since}")

    found_by_query = {}
    for q in queries:
        pmids = pubmed_search(windowed_query(q, args.since), max_results=args.max)
        found_by_query[q] = pmids
        print(f"  [{len(pmids):>4}] {q[:70]}")

    new_by_query, all_new = compute_delta(found_by_query, frozen)
    print(f"New (not in frozen corpus): {len(all_new)}")

    records = []
    if all_new and not args.dry_run:
        records = [tag_record(m) for m in fetch_pubmed_metadata(all_new)]
        out_dir = LIVING_DIR / run_date
        out_dir.mkdir(parents=True, exist_ok=True)
        with (out_dir / "index.jsonl").open("w") as fh:
            for r in records:
                fh.write(json.dumps(r) + "\n")
        print(f"Wrote {len(records)} records to {out_dir / 'index.jsonl'}")
    elif all_new:
        # Dry-run: still tag from the metadata so the changelog is informative,
        # but write nothing to the corpus.
        records = [tag_record(m) for m in fetch_pubmed_metadata(all_new)]

    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    report = REPORT_DIR / f"{run_date}.md"
    report.write_text(
        format_changelog(run_date, args.since, new_by_query, records, len(frozen))
    )
    print(f"Wrote changelog {report}")


if __name__ == "__main__":
    main()
