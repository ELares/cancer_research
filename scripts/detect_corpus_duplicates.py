"""Detect preprint/published (and other) duplicate article pairs in the corpus (#535).

DOI-based dedup misses preprint/published twins: the bioRxiv/medRxiv version and the
journal version carry different DOIs and different journals, so both end up in the
corpus. Each twin double-counts the work AND inflates two different open-access tiers
(the preprint as `green`, the published version as `gold`/`hybrid`/`bronze`), which
biases the open-access-skew analysis the manuscript relies on (#348).

This is a GUARD, not a remover: per the maintainer's call the records are kept, and
this script simply flags the twins (a reproducible detector + a committed audit) so
the duplication is visible and the OA-bias numbers can be read with it in mind.

    python scripts/detect_corpus_duplicates.py   # writes analysis/corpus-duplicate-audit.md
"""

import json
import re
from collections import defaultdict
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
INDEX = REPO_ROOT / "corpus" / "INDEX.jsonl"
REPORT = REPO_ROOT / "analysis" / "corpus-duplicate-audit.md"

PREPRINT_RE = re.compile(r"biorxiv|medrxiv|arxiv|research square|preprints?\b|ssrn", re.I)


def norm_title(t: str) -> str:
    return re.sub(r"[^a-z0-9]", "", (t or "").lower())


def main():
    recs = [json.loads(line) for line in INDEX.read_text(encoding="utf-8").splitlines()]
    by_title = defaultdict(list)
    for r in recs:
        key = norm_title(r.get("title"))
        if len(key) >= 20:  # ignore very short/empty titles to avoid spurious hits
            by_title[key].append(r)

    groups = [g for g in by_title.values() if len(g) > 1]
    # Split into preprint/published twins vs other same-title collisions.
    preprint_twins, other = [], []
    for g in groups:
        if any(PREPRINT_RE.search(r.get("journal", "") or "") for r in g):
            preprint_twins.append(g)
        else:
            other.append(g)

    def fmt_row(r):
        return (f"| {r['pmid']} | {(r.get('journal') or '')[:40]} | {r.get('year')} | "
                f"{r.get('oa_status', '')} | {'preprint' if PREPRINT_RE.search(r.get('journal','') or '') else 'published'} |")

    lines = [
        "# Corpus duplicate / preprint-twin audit (#535)",
        "",
        "Records sharing a near-identical title (a guard, not a remover — both records "
        "are kept). Preprint/published twins double-count the work and inflate two "
        "different open-access tiers, biasing the #348 OA-skew analysis; same-title "
        "collisions without a preprint are listed separately and may be genuine "
        "duplicates, errata, or distinct articles that happen to share a title.",
        "",
        f"- Preprint/published twin groups: **{len(preprint_twins)}**",
        f"- Other same-title groups: **{len(other)}**",
        "",
        "## Preprint/published twins",
        "",
    ]
    for g in sorted(preprint_twins, key=lambda g: norm_title(g[0].get("title"))):
        lines.append(f"**{(g[0].get('title') or '')[:90]}**")
        lines.append("")
        lines.append("| PMID | Journal | Year | OA status | Kind |")
        lines.append("|---|---|---|---|---|")
        lines += [fmt_row(r) for r in g]
        lines.append("")
    if other:
        lines += ["## Other same-title collisions (review individually)", "",
                  "| PMID | Journal | Year | OA status | Kind |", "|---|---|---|---|---|"]
        for g in sorted(other, key=lambda g: norm_title(g[0].get("title"))):
            lines += [fmt_row(r) for r in g]
        lines.append("")

    REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"{len(preprint_twins)} preprint/published twin groups, {len(other)} other "
          f"same-title groups. Wrote {REPORT.relative_to(REPO_ROOT)}.")


if __name__ == "__main__":
    main()
