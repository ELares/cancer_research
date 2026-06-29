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

# Hand-verified verdicts for "other" same-title collisions (#567), keyed by the
# frozenset of the group's PMIDs (as strings). Each was checked against NCBI
# E-utilities (publication type + authors + DOI), not assumed — see the audit.
# Groups not in this map render as "not yet reviewed" so a future collision is
# flagged rather than silently presumed benign.
VERDICTS = {
    frozenset({"38487722", "35433483"}): (
        "**Not a content duplicate — two distinct corrigenda.** Both records are PubMed "
        "type *Published Erratum* (corrigenda) to the same original article (tumor "
        "treating fields + mild hyperthermia for pancreatic cancer), same authors "
        "(Bai, Pfeifer, Gross, De La Torre) but different DOIs "
        "(10.3389/fonc.2022.889215 and 10.3389/fonc.2024.1343421), published two years "
        "apart. The shared title is the corrigendum-inherits-original-title convention. "
        "Both are non-research item types, so they do not inflate the #348 OA "
        "research-article bias analysis. Keep both (the #535 guard-not-remover policy)."
    ),
    frozenset({"19997112", "19997110"}): (
        "**Not duplicates — two independent letters to the editor.** Both records are "
        "PubMed type *Comment / Letter* responding to the same source article "
        "(\"High-intensity-focused ultrasound ... the first UK series\"), but by "
        "different authors (Eggener, Gonzalgo & Yossepowitch vs. Clark) with different "
        "DOIs (10.1038/sj.bjc.6605455 and 10.1038/sj.bjc.6605453) — distinct "
        "correspondence sharing a 'Regarding:' title. Non-research item types, so no "
        "OA-bias inflation. Keep both."
    ),
}


def norm_title(t: str) -> str:
    return re.sub(r"[^a-z0-9]", "", (t or "").lower())


def verdict_for(group) -> str:
    """Return the hand-verified verdict for a same-title group, or a not-reviewed
    placeholder so a newly-appearing collision is visibly flagged (#567)."""
    return VERDICTS.get(
        frozenset(str(r.get("pmid")) for r in group),
        "_Verdict: not yet reviewed — check the publication types/DOIs on NCBI and add "
        "an entry to `VERDICTS` in `scripts/detect_corpus_duplicates.py`._",
    )


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
        lines += [
            "## Other same-title collisions (reviewed individually, #567)",
            "",
            "Each group below carries a hand-verified verdict (publication type / authors "
            "/ DOI checked on NCBI). These are kept per the guard-not-remover policy.",
            "",
        ]
        for g in sorted(other, key=lambda g: norm_title(g[0].get("title"))):
            lines.append(f"**{(g[0].get('title') or '')[:90]}**")
            lines.append("")
            lines.append("| PMID | Journal | Year | OA status | Kind |")
            lines.append("|---|---|---|---|---|")
            lines += [fmt_row(r) for r in g]
            lines.append("")
            lines.append(verdict_for(g))
            lines.append("")

    REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"{len(preprint_twins)} preprint/published twin groups, {len(other)} other "
          f"same-title groups. Wrote {REPORT.relative_to(REPO_ROOT)}.")


if __name__ == "__main__":
    main()
