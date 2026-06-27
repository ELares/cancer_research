"""Rank candidate wet-lab / domain collaborators from the corpus itself (#541).

Track-A accelerator for #498 (collaborator recruitment). The collaborator brief
deliberately names no individuals; this script produces a REPRODUCIBLE, first-pass
SEED LIST of senior (last) authors who publish on the ferroptosis biology the
project's predictions depend on, ranked by recency, frequency, journal breadth,
and citation impact — drawn ONLY from authorship already present in the public
corpus (no fabricated names or affiliations).

OUTPUT IS LOCAL-ONLY BY DESIGN. It writes to `local/` (gitignored), so a ranked
list of named researchers is never committed to the public repo. Run it yourself:

    python scripts/rank_collaborator_candidates.py
    # -> local/collaborator-candidates.csv + .md

HONEST SCOPE (read before using the output):
- This is a SEED LIST for a human to vet, NOT a ranking of people. Author-name
  strings are not disambiguated (no ORCID): two different "Wang Wei" collapse into
  one row, and one person publishing under name variants splits into several.
- Affiliations are not parsed (the corpus frontmatter does not carry structured
  affiliations), so the list tells you WHO publishes on the relevant biology, not
  WHERE they are. Contacting anyone is a human step that requires checking the
  current affiliation, conflicts, and willingness.
- The corpus is open-access-skewed and keyword-scoped (see the manuscript's
  coverage caveats), so absence from this list means "not surfaced here," not
  "does not work on ferroptosis."
- Last-author = senior/PI is a heuristic; some fields order authors differently.
"""

import csv
import re
from collections import defaultdict
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
PMID_DIR = REPO_ROOT / "corpus" / "by-pmid"
OUT_DIR = REPO_ROOT / "local"
RECENT_YEAR = 2022  # "recent" threshold for the recency signal

# Ferroptosis-relevant matcher: the biology the P1-P8 predictions depend on
# (GPX4/FSP1 parallel repair, System Xc-, the iron/lipid-peroxidation axis,
# persister/drug-tolerance). Word-boundary where the token is short/ambiguous.
RELEVANCE = re.compile(
    r"\bferroptos|\bGPX4\b|\bFSP1\b|\bAIFM2\b|\bSLC7A11\b|\bxCT\b|\bACSL4\b|"
    r"\berastin\b|\bRSL3\b|\bML162\b|\bML210\b|\blipid peroxidation\b|"
    r"\blabile iron\b|drug-tolerant persister|\bDHODH\b|\bGCH1\b",
    re.IGNORECASE,
)


def parse_frontmatter(text: str) -> dict:
    """Minimal YAML-frontmatter parse (authors block list + scalar fields)."""
    m = re.match(r"^---\n(.*?)\n---", text, re.S)
    if not m:
        return {}
    fm = {"authors": []}
    in_authors = False
    for line in m.group(1).splitlines():
        if line.startswith("authors:"):
            in_authors = True
            continue
        if in_authors:
            ma = re.match(r"\s*-\s+(.*\S)\s*$", line)
            if ma:
                fm["authors"].append(ma.group(1).strip())
                continue
            in_authors = False
        ms = re.match(r"^(\w+):\s*(.*)$", line)
        if ms:
            fm[ms.group(1)] = ms.group(2).strip().strip("'\"")
    return fm


def collect():
    by_author = defaultdict(lambda: {
        "papers": 0, "recent": 0, "latest_year": 0, "citations": 0,
        "journals": set(), "pmids": [],
    })
    n_scanned = n_relevant = 0
    for fp in sorted(PMID_DIR.glob("*.md")):
        try:
            text = fp.read_text(encoding="utf-8", errors="ignore")
        except OSError:
            continue
        n_scanned += 1
        if not RELEVANCE.search(text):
            continue
        n_relevant += 1
        fm = parse_frontmatter(text)
        authors = fm.get("authors") or []
        if not authors:
            continue
        senior = authors[-1]  # last author = senior/PI heuristic
        try:
            year = int(fm.get("year") or 0)
        except ValueError:
            year = 0
        try:
            cites = int(float(fm.get("cited_by_count") or 0))
        except ValueError:
            cites = 0
        rec = by_author[senior]
        rec["papers"] += 1
        rec["citations"] += cites
        rec["latest_year"] = max(rec["latest_year"], year)
        if year >= RECENT_YEAR:
            rec["recent"] += 1
        if fm.get("journal"):
            rec["journals"].add(fm["journal"])
        rec["pmids"].append(fp.stem)
    return by_author, n_scanned, n_relevant


def score(rec: dict) -> float:
    """Transparent composite: reward recent + frequent + broad + cited. Recency
    is weighted highest (an active lab matters more than a historical one)."""
    return (
        3.0 * rec["recent"]
        + 1.0 * rec["papers"]
        + 0.5 * len(rec["journals"])
        + 0.01 * rec["citations"]
    )


def main():
    by_author, n_scanned, n_relevant = collect()
    rows = []
    for author, rec in by_author.items():
        rows.append({
            "senior_author": author,
            "score": round(score(rec), 2),
            "ferroptosis_papers": rec["papers"],
            "recent_papers_2022plus": rec["recent"],
            "latest_year": rec["latest_year"],
            "distinct_journals": len(rec["journals"]),
            "total_citations": rec["citations"],
            "pmids": " ".join(sorted(rec["pmids"])),
        })
    # Keep only authors with >=2 relevant papers (a single paper is too noisy for
    # a seed list); rank by score then recency.
    rows = [r for r in rows if r["ferroptosis_papers"] >= 2]
    rows.sort(key=lambda r: (-r["score"], -r["latest_year"]))

    OUT_DIR.mkdir(exist_ok=True)
    csv_path = OUT_DIR / "collaborator-candidates.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()) if rows else
                           ["senior_author", "score", "ferroptosis_papers",
                            "recent_papers_2022plus", "latest_year",
                            "distinct_journals", "total_citations", "pmids"])
        w.writeheader()
        w.writerows(rows)

    md_path = OUT_DIR / "collaborator-candidates.md"
    lines = [
        "# Collaborator candidate seed list (LOCAL-ONLY, not committed)",
        "",
        "> First-pass, reproducible seed list mined from the corpus (#541). **A list "
        "of authors to vet, NOT a ranking of people.** Names are not disambiguated, "
        "affiliations are not parsed, and the corpus is OA-skewed. Verify current "
        "affiliation, conflicts, and willingness before any outreach (a human step). "
        "See the docstring of `scripts/rank_collaborator_candidates.py` for the full "
        "scope caveats.",
        "",
        f"Scanned {n_scanned} corpus articles; {n_relevant} matched the "
        "ferroptosis-relevance filter; "
        f"{len(rows)} senior authors with >=2 relevant papers.",
        "",
        "| Rank | Senior author | Score | Ferro papers | Recent (>=2022) | Latest | Journals | Citations |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for i, r in enumerate(rows[:40], 1):
        lines.append(
            f"| {i} | {r['senior_author']} | {r['score']} | "
            f"{r['ferroptosis_papers']} | {r['recent_papers_2022plus']} | "
            f"{r['latest_year']} | {r['distinct_journals']} | {r['total_citations']} |"
        )
    md_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"Scanned {n_scanned} articles, {n_relevant} ferroptosis-relevant.")
    print(f"Wrote {len(rows)} candidates to:")
    print(f"  {csv_path.relative_to(REPO_ROOT)}")
    print(f"  {md_path.relative_to(REPO_ROOT)}")
    print("(local/ is gitignored — these named lists are not committed.)")


if __name__ == "__main__":
    main()
