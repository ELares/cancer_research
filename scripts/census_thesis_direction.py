#!/usr/bin/env python3
"""Classify resistance-leg candidates and bind adjudication to their cohort.

The fixed vocabulary screens titles and abstracts for two framings: exploiting
ferroptosis against resistant cells, and resistance to ferroptosis. A substring
match cannot establish either biological direction. Fresh scans therefore
report lexical candidates until every singly classified candidate has a
completed, identity-verified adjudication from the selected cohort.

Historical counts and sampled corrections remain available through
``--render-only``. Their original titles-only labels cannot verify article
identity, cohort coverage, or sampling representativeness and are never applied
to a new scan. The vocabulary is retained without tuning to either result.
"""
import argparse
import json
import re
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from census_adjudication import (  # noqa: E402
    AdjudicationCohort, candidate_csv, read_adjudication,
)
from census_input import atlas_root, iter_census_records  # noqa: E402

REPO = Path(__file__).resolve().parent.parent
RECORDS = atlas_root() / "records"
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
    cohort = AdjudicationCohort("thesis-direction", ("exploit", "obstacle"))
    unclassified = []
    total = 0
    for r in iter_census_records(RECORDS, stride):
        ms = {m.lower() for m in (r.get("mesh") or [])}
        if FERRO not in ms or RESIST not in ms:
            continue
        total += 1
        blob = f"{r.get('title') or ''} {r.get('abstract') or ''}".lower()
        e, o = bool(EXPLOIT.search(blob)), bool(OBSTACLE.search(blob))
        key = ("both" if e and o else "exploit" if e
               else "obstacle" if o else "neither")
        counts[key] += 1
        cohort.add(r, key)
        if key == "neither" and len(unclassified) < N_SAMPLE:
            unclassified.append(
                {"year": r.get("year"), "title": (r.get("title") or "")[:110]})
    return {"total": total, "counts": counts,
            "unclassified_sample": unclassified, "cohort": cohort.as_dict()}


def assemble(d: dict, *, adjudication_rows=None) -> dict:
    """Derive lexical summaries and, if supplied, complete verified decisions.

    ``adjudication_rows`` comes from ``read_adjudication`` for this scan's
    cohort. Embedded historical decisions are deliberately not carried into
    this calculation; offline rendering reads their snapshot directly.
    """
    c = d["counts"]
    # SINGLY CLASSIFIED ONLY. Articles carrying both vocabularies are outside
    # the candidate adjudication denominator, alongside unmatched articles.
    single = c["exploit"] + c["obstacle"]
    out = dict(d)
    out["adjudication"] = {}
    out.pop("points_the_projects_way", None)
    out["classified"] = single
    out["unclassified"] = c["neither"] + c["both"]
    out["exploit_share_of_classified"] = (
        round(100 * c["exploit"] / single, 1) if single else None)
    out["direction_ratio"] = (
        round(c["exploit"] / c["obstacle"], 1) if c["obstacle"] else None)
    out["unclassified_share"] = (
        round(100 * out["unclassified"] / d["total"], 1) if d["total"] else None)
    if adjudication_rows is not None:
        counts = Counter(row["adjudicated"] for row in adjudication_rows)
        decided = counts["exploit"] + counts["obstacle"]
        direction = (None if not decided else "tie"
                     if counts["exploit"] == counts["obstacle"] else "exploit"
                     if counts["exploit"] > counts["obstacle"] else "obstacle")
        out["adjudication"] = {
            "mode": "complete-current-cohort",
            "cohort_sha256": d["cohort"]["cohort_sha256"],
            "rows": len(adjudication_rows),
            "decisions": adjudication_rows,
            "counts": {label: counts[label]
                       for label in ("exploit", "obstacle", "ambiguous")},
            "decided": decided,
            "exploit_share_of_decided": (
                round(100 * counts["exploit"] / decided, 1) if decided else None),
            "direction": direction,
        }
    return out


def render(d: dict) -> str:
    c = d["counts"]
    L = ["# Which way the resistance leg points\n"]
    L.append(
        "Generated by `scripts/census_thesis_direction.py`. The selected "
        f"indexed census contains {d['total']:,} articles carrying both "
        "`Ferroptosis` and `Drug Resistance, Neoplasm`. The fixed title and "
        "abstract patterns screen for exploiting ferroptosis against resistant "
        "cells and for resistance to ferroptosis. These are lexical candidate "
        "labels; matching a phrase does not establish the article's direction.\n"
    )
    L.append("| lexical candidate label | articles | share of intersection |")
    L.append("|---|--:|--:|")
    for k, label in (("exploit", "exploit pattern only"),
                     ("obstacle", "obstacle pattern only"),
                     ("both", "both patterns"),
                     ("neither", "neither pattern")):
        share = f"{100 * c[k] / d['total']:.0f}%" if d["total"] else "not applicable"
        L.append(f"| {label} | {c[k]:,} | {share} |")
    L.append("")
    if not d["total"]:
        L.append("No articles match the descriptor intersection in this "
                 "selection; no direction comparison is available.\n")
    elif not d["classified"]:
        L.append("No singly classified candidates are available for a "
                 "direction comparison.\n")
    else:
        L.append(
            f"Of {d['classified']:,} singly classified candidates, "
            f"{d['exploit_share_of_classified']}% match only the exploit "
            "pattern. This percentage describes pattern matches and does not "
            "establish which biological framing leads.\n"
        )

    a = d.get("adjudication") or {}
    if a.get("mode") == "complete-current-cohort":
        L.append("## Complete adjudication of the selected candidates\n")
        L.append(
            f"All {a['rows']:,} singly classified candidates have validated "
            "record identities, regex labels and completed decisions for "
            "this selected cohort. Articles matching both or neither pattern "
            "remain outside this adjudication denominator.\n"
        )
        L.append("| adjudicated framing | articles |")
        L.append("|---|--:|")
        for label in ("exploit", "obstacle", "ambiguous"):
            L.append(f"| {label} | {a['counts'][label]:,} |")
        L.append("")
        if not a["decided"]:
            L.append("There are no decided adjudications; no direction "
                     "comparison is available.\n")
        else:
            L.append(
                f"Among {a['decided']:,} decided candidates, the adjudicated "
                f"exploit share is {a['exploit_share_of_decided']}%. "
                "Ambiguous decisions are excluded from this denominator. "
                "This is a direct count of the adjudicated candidates, without "
                "extrapolation to unmatched articles.\n"
            )
            if a["direction"] == "tie":
                L.append("The adjudicated exploit and obstacle counts are "
                         "equal in this candidate set.\n")
            else:
                L.append(f"The adjudicated {a['direction']} framing leads "
                         "within this candidate set.\n")
    elif a:
        lo, hi = a["corrected_range"]
        L.append("## Historical conditional correction\n")
        L.append(
            "This report renders the embedded historical snapshot. The "
            "original titles-only adjudication CSV does not establish a "
            "verified article-identity link or complete coverage for that "
            "cohort, and the exploit sample's representativeness is "
            "unverified. These labels are not applied to fresh scans.\n"
        )
        L.append(
            "The historical calculation assumed that every obstacle-labelled "
            "article was adjudicated and that the exploit sample's precision "
            "could be applied to the remaining exploit-labelled articles. "
            f"Its labels record {a['by_label']['obstacle']['adjudicated_exploit']} "
            f"exploit decisions among {a['by_label']['obstacle']['sampled']} "
            "obstacle-labelled rows and "
            f"{a['by_label']['exploit']['adjudicated_obstacle']} obstacle "
            f"decisions among {a['by_label']['exploit']['sampled']} "
            "exploit-labelled rows, illustrating errors in both directions.\n"
        )
        L.append(
            f"The preserved conditional calculation gives an exploit share "
            f"of **{a['corrected_exploit_share']}%** (range **{lo}-{hi}%**), "
            "using the sampled precision's interval. This range captures "
            "only that calculation's sampling term. It does not establish "
            "uncertainty from cohort identity, sampling representativeness, "
            "labeling errors or missed articles, and cannot establish the "
            "direction of the whole resistance literature.\n"
        )
        L.append(
            "The original counts and calculation remain in "
            "`analysis/census-thesis-direction.json`; the original labels "
            "remain in `analysis/thesis-direction-adjudication.csv`. Neither "
            "artifact is relabeled by this rendering.\n"
        )
    else:
        L.append("## Adjudication required for direction\n")
        L.append(
            "No validated adjudication is attached to this selected cohort. "
            "The regex candidates alone do not establish a biological "
            "direction or a corrected share. Export candidates with "
            "`--export-candidates PATH`, complete each decision and reason, "
            "then rescan the same input and stride with `--adjudication PATH`.\n"
        )

    L.append("## What this does not establish\n")
    L.append(
        "A field framing something as exploitable is evidence about what "
        "the field is TRYING. Even a fully adjudicated framing count does not "
        "establish that inducing ferroptosis kills resistant cells or that "
        "the project's thesis is right.\n"
    )
    L.append("## Recall remains unmeasured\n")
    if d["unclassified_share"] is not None:
        L.append(
            f"{d['unclassified_share']}% of the intersection matches neither "
            "pattern or both. These articles are excluded from the singly "
            "classified candidate denominator. Their biological direction "
            "has not been established by this screen.\n"
        )
    samples = d.get("unclassified_sample", [])
    if samples:
        L.append("Examples matching neither pattern:\n")
        for sample in samples:
            L.append(f"- {sample['year']} — {sample['title']}")
        L.append("")
    L.append(
        "The patterns were fixed before any result was read and were not "
        "adjusted afterwards. The exploit vocabulary is the thesis vocabulary; "
        "widening it after seeing the split would risk fitting the answer. "
        "Candidate coverage does not measure recall over the whole intersection.\n"
    )
    return "\n".join(L)


def _same_path(left: Path, right: Path) -> bool:
    return left.resolve() == right.resolve() or (
        left.exists() and right.exists() and left.samefile(right))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--render-only", action="store_true")
    ap.add_argument("--output-dir", type=Path)
    decisions = ap.add_mutually_exclusive_group()
    decisions.add_argument("--adjudication", type=Path)
    decisions.add_argument("--export-candidates", type=Path)
    a = ap.parse_args()
    if a.stride < 1:
        ap.error("--stride must be a positive integer")
    if a.render_only and (a.adjudication or a.export_candidates):
        ap.error("--render-only cannot be combined with --adjudication or "
                 "--export-candidates")
    out_json = a.output_dir / OUT_JSON.name if a.output_dir else OUT_JSON
    out_md = a.output_dir / OUT_MD.name if a.output_dir else OUT_MD
    if _same_path(out_json, out_md):
        ap.error("JSON and Markdown reports must not alias each other")
    if a.adjudication and any(
            _same_path(a.adjudication, report) for report in (out_json, out_md)):
        ap.error("--adjudication must not alias a JSON or Markdown report")
    if a.export_candidates and any(
            _same_path(a.export_candidates, report) for report in (out_json, out_md)):
        ap.error("--export-candidates must not alias a JSON or Markdown report")

    export_text = None
    json_text = None
    if a.render_only:
        d = json.loads(out_json.read_text(encoding="utf-8"))
    else:
        raw = scan(a.stride)
        rows = (read_adjudication(a.adjudication, raw["cohort"],
                                 ("exploit", "obstacle", "ambiguous"))
                if a.adjudication else None)
        d = assemble(raw, adjudication_rows=rows)
        if a.export_candidates:
            export_text = candidate_csv(raw["cohort"])
        json_text = json.dumps(d, indent=1) + "\n"
    markdown = render(d)
    # Parsing, coverage validation, serialization and rendering all finish
    # before an existing report or candidate file can be replaced.
    for path in (out_md, out_json) if json_text is not None else (out_md,):
        path.parent.mkdir(parents=True, exist_ok=True)
    if export_text is not None:
        a.export_candidates.parent.mkdir(parents=True, exist_ok=True)
        a.export_candidates.write_text(export_text, encoding="utf-8")
    if json_text is not None:
        out_json.write_text(json_text, encoding="utf-8")
    out_md.write_text(markdown, encoding="utf-8")
    print(f"wrote {out_md}")
    print(f"  {d['total']} intersection articles, "
          f"{d['classified']} singly classified candidates")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
