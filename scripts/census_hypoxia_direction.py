#!/usr/bin/env python3
"""Generate hypoxia-direction candidates and report explicitly bound labels.

The historical title-only adjudication documents a failed keyword classifier:
"hypoxia-induced ferroptosis resistance" contains "hypoxia-induced ferroptosis".
Its numerical evidence remains available as an offline snapshot, but neither
its CSV nor its stored report identifies the complete historical candidate set.
Fresh scans therefore require a completed candidate export bound to the selected
cohort before any adjudicated measurement can be reported.
"""
import argparse
from collections import Counter
import json
import re
from pathlib import Path

from atlas_baseline import atlas_root
from census_adjudication import AdjudicationCohort, candidate_csv, read_adjudication
from census_input import iter_census_records

REPO = Path(__file__).resolve().parent.parent
RECORDS = atlas_root() / "records"
OUT_MD = REPO / "analysis/census-hypoxia-direction.md"
OUT_JSON = REPO / "analysis/census-hypoxia-direction.json"

FERRO = "ferroptosis"
# No single hypoxia descriptor carries this literature.
HYPOXIA = {"cell hypoxia", "tumor hypoxia", "hypoxia", "anoxia",
           "hypoxia-inducible factor 1, alpha subunit", "hypoxia-inducible factor 1"}

# Historical patterns are retained as candidate generators, not measurements.
PROTECTS = re.compile(
    r"\b(hypoxia[- ]?(?:induced|mediated|driven)? ?(?:resistance|protect\w*|"
    r"inhibit\w*|suppress\w*|attenuat\w*|reduc\w*)|"
    r"resistan\w* to ferroptosis under hypoxia|"
    r"protect\w* (?:cells |tumou?r cells )?(?:from|against) ferroptosis)\b")
SENSITISES = re.compile(
    r"\b(hypoxia[- ]?(?:induced|mediated|driven)? ?(?:sensitiz\w*|sensitis\w*|"
    r"promot\w*|enhanc\w*|potentiat\w*|trigger\w*|aggravat\w*)|"
    r"hypoxia[- ]induced ferroptosis|"
    r"ferroptosis under hypoxia\w* condition)\b")
N_SAMPLE = 6
DECISIONS = {"protects", "sensitises", "off-topic", "ambiguous"}


def scan(stride: int = 1) -> dict:
    counts = {"protects": 0, "sensitises": 0, "both": 0, "neither": 0}
    sample = []
    total = 0
    cohort = AdjudicationCohort("hypoxia-direction", ("protects", "sensitises"))
    for record in iter_census_records(RECORDS, stride):
        ms = {m.lower() for m in (record.get("mesh") or [])}
        if FERRO not in ms or not (ms & HYPOXIA):
            continue
        total += 1
        blob = f"{record.get('title') or ''} {record.get('abstract') or ''}".lower()
        p, s = bool(PROTECTS.search(blob)), bool(SENSITISES.search(blob))
        key = ("both" if p and s else "protects" if p
               else "sensitises" if s else "neither")
        counts[key] += 1
        cohort.add(record, key)
        if key in ("protects", "sensitises") and len(sample) < N_SAMPLE:
            sample.append({"framing": key, "year": record.get("year"),
                           "title": (record.get("title") or "")[:105]})
    return {"total": total, "counts": counts, "sample": sample,
            "hypoxia_descriptors": sorted(HYPOXIA), "cohort": cohort.as_dict()}


def _summarize_adjudication(rows: list[dict], cohort: dict) -> dict:
    labels = Counter(row["adjudicated"] for row in rows)
    agree = sum(row["regex_label"] == row["adjudicated"] for row in rows)
    reversed_ = [row for row in rows
                 if row["adjudicated"] in ("protects", "sensitises")
                 and row["regex_label"] != row["adjudicated"]]
    return {
        "mode": "complete-current-cohort",
        "cohort_sha256": cohort["cohort_sha256"],
        "decisions": rows,
        "n": len(rows),
        "labels": dict(labels),
        "regex_agreement": round(100 * agree / len(rows), 1) if rows else None,
        "regex_reversed": len(reversed_),
        "directional": labels["protects"] + labels["sensitises"],
        "protects": labels["protects"],
        "sensitises": labels["sensitises"],
    }


def assemble(d: dict, *, adjudication_rows: list[dict] | None = None) -> dict:
    """Assemble a fresh measurement using only explicitly validated rows.

    Offline snapshots are passed directly to render(), never through here.
    """
    c = d["counts"]
    single = c["protects"] + c["sensitises"]
    out = dict(d)
    for key in ("adjudication", "adj_protects_share", "adj_ci", "protects_ci",
                "dominance_threshold", "min_classified", "interpretable",
                "verdict", "verdict_survives_interval", "regex_reversed_cancer",
                "regex_direction_reversed_by_adjudication"):
        out.pop(key, None)
    out["classified"] = single
    out["unclassified"] = c["neither"] + c["both"]
    out["unclassified_share"] = (round(100 * out["unclassified"] / d["total"], 1)
                                 if d["total"] else None)
    out["protects_share"] = round(100 * c["protects"] / single, 1) if single else None
    out["adjudication"] = {}
    out["adj_protects_share"] = None
    out["regex_direction_reversed_by_adjudication"] = False
    if adjudication_rows is not None:
        adj = _summarize_adjudication(adjudication_rows, d["cohort"])
        out["adjudication"] = adj
        out["adj_protects_share"] = (
            round(100 * adj["protects"] / adj["directional"], 1)
            if adj["directional"] else None)
        if out["adj_protects_share"] is not None and out["protects_share"] is not None:
            out["regex_direction_reversed_by_adjudication"] = (
                (out["adj_protects_share"] - 50) * (out["protects_share"] - 50) < 0)
    return out


def render(d: dict) -> str:
    c = d["counts"]
    a = d.get("adjudication") or {}
    historical = bool(a and a.get("mode") != "complete-current-cohort")
    lines = ["# How the literature splits on hypoxia and ferroptosis\n"]
    lines.append(
        f"Generated by `scripts/census_hypoxia_direction.py`. "
        f"{d['total']:,} census articles carry `Ferroptosis` together with a "
        "hypoxia-family descriptor ("
        + ", ".join(f"`{h}`" for h in d["hypoxia_descriptors"])
        + "). The descriptor family is taken as a union.\n")
    if historical:
        lines.append(
            "**Historical snapshot; record linkage is not verifiable.** The "
            "stored title-only adjudication has no record identifiers or "
            "complete candidate inventory. Its numerical evidence is preserved "
            "below, but the available artifacts cannot establish which census "
            "records received those labels. Equal aggregate counts and a few "
            "truncated title examples do not establish that linkage.\n")
    lines.append(
        "Section 7.1 describes the pharmacologic ferroptosis response to "
        "simulated hypoxia as contested. This report records candidate "
        "framings and, when available, adjudicated directions. Candidate "
        "counts alone cannot establish whether that literature description "
        "is warranted.\n")
    lines.extend(["| regex candidate framing | articles | share of all |",
                  "|---|--:|--:|"])
    for key, label in (("protects", "protective pattern"),
                       ("sensitises", "sensitising pattern"),
                       ("both", "both patterns"),
                       ("neither", "neither pattern")):
        share = f"{100 * c[key] / d['total']:.0f}%" if d["total"] else "-"
        lines.append(f"| {label} | {c[key]:,} | {share} |")
    lines.append("")
    if not a:
        lines.append(
            "**No adjudications are bound to this selected cohort.** The "
            "table contains regex candidates only; no adjudicated direction "
            "or classifier accuracy is measured. Export candidates with "
            "`--export-candidates PATH`, complete the labels and reasons, "
            "then provide that file with `--adjudication PATH`.\n")
    else:
        lines.append("## Historical adjudication\n" if historical
                     else "## Adjudication of the selected candidates\n")
        if historical:
            lines.append(
                f"The embedded historical summary contains {a['n']} "
                "title-only labels. Its reported classifier agreement is "
                f"**{a['regex_agreement']}%**, with **{a['regex_reversed']}** "
                "direction reversals; the stored reasons mark "
                f"{a['regex_reversed_cancer']} of those as cancer-context. "
                "These are historical summary values, not a verified accuracy "
                "estimate for a newly selected census.\n")
        else:
            agreement = (f"{a['regex_agreement']}%" if a["regex_agreement"] is not None
                         else "undefined because there are no candidates")
            lines.append(
                f"All {a['n']} singly classified candidates have validated "
                "labels and reasons bound to the selected records and cohort "
                f"`{a['cohort_sha256']}`. Classifier agreement is "
                f"{agreement}, with {a['regex_reversed']} direction reversals. "
                "The completed rows are embedded in the JSON report.\n")
        share = d.get("adj_protects_share")
        if a["directional"] and share is not None:
            lines.append(
                f"The {'historical ' if historical else ''}adjudication "
                f"labels {a['protects']} articles protective and "
                f"{a['sensitises']} sensitising: a protective share of "
                f"**{share}%** over {a['directional']} directional articles.\n")
            if d.get("regex_direction_reversed_by_adjudication"):
                lines.append(
                    "The adjudication reverses the recorded majority "
                    "direction: the protective share changes from "
                    f"{d['protects_share']}% in the regex counts to {share}% "
                    "in the adjudication.\n")
            if historical and d.get("adj_ci") is not None:
                lo, hi = d["adj_ci"]
                lines.append(
                    f"The historical calculation reports a 95% Wilson "
                    f"interval of {lo}-{hi}%. It does not account for "
                    "candidate selection, missed articles, adjudication "
                    "errors, or unverified record linkage.\n")
                if lo <= 50 <= hi:
                    lines.append(
                        "That interval is wide enough to contain an even split. "
                        "These counts do not establish their ratio in the broader "
                        "literature or validate the manuscript's framing.\n")
            else:
                lines.append(
                    "This is a direct count of the completed candidate "
                    "adjudication, without a sampling interval or "
                    "extrapolation to the broader literature.\n")
            if a["protects"] and a["sensitises"]:
                lines.append("Both directions occur in the recorded labels.\n")
        else:
            lines.append(
                "There are no directional adjudications, so the protective "
                "share is undefined.\n")
        lines.append(
            f"{a['labels'].get('off-topic', 0)} of {a['n']} labels are "
            f"off-topic and {a['labels'].get('ambiguous', 0)} are ambiguous. "
            "They are excluded from the directional denominator, leaving "
            f"{a['directional']} directional articles.\n")
        lines.append("## Limits of the adjudication\n")
        if historical:
            lines.append(
                f"Titles only, one adjudicator, and {a['directional']} "
                "directional articles. Historical labels and reasons remain "
                "in `analysis/hypoxia-direction-adjudication.csv`; some titles "
                "are abbreviated. The legacy file cannot be automatically "
                "applied to a fresh cohort. Re-rendering reads this embedded "
                "snapshot without reloading that CSV or reconstructing "
                "missing record identities.\n")
        else:
            lines.append(
                "Complete coverage applies only to the selected singly "
                "classified candidates. It does not include records matching "
                "both patterns or neither pattern, establish independent "
                "raters, or demonstrate representativeness of the literature.\n")
    lines.append("## Why the regex only generates candidates\n")
    lines.append(
        "The regex is demoted to generating candidates because "
        "**\"hypoxia-induced ferroptosis resistance\" CONTAINS "
        "\"hypoxia-induced ferroptosis\"**. A phrase asserting protection "
        "contains a substring that the classifier reads as sensitisation. "
        "The retained historical labels document this failure.\n")
    if d["sample"]:
        lines.append("### Examples of regex candidates, not adjudicated directions\n")
        for sample in d["sample"]:
            lines.append(f"- *{sample['framing']}* — {sample['year']} — {sample['title']}")
        lines.append("")
    lines.append("## What this cannot do\n")
    lines.append(
        "Settle the biology. A count of reported framings does not establish "
        "how hypoxia affects ferroptosis in a particular cell type, inducer, "
        "or pathway perturbation. It does not validate the simulation.\n")
    if d["total"]:
        lines.append(
            f"{d['unclassified_share']}% of the intersection matches neither "
            "pattern or both, and is outside the singly classified candidate "
            "set. The original patterns are retained unchanged.\n")
    else:
        lines.append(
            "The selected census has no hypoxia/ferroptosis intersection "
            "articles, so no direction can be measured.\n")
    return "\n".join(lines)


def _aliases(left: Path, right: Path) -> bool:
    return (left.resolve() == right.resolve()
            or left.exists() and right.exists() and left.samefile(right))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--render-only", action="store_true")
    ap.add_argument("--output-dir", type=Path)
    labels = ap.add_mutually_exclusive_group()
    labels.add_argument("--adjudication", type=Path)
    labels.add_argument("--export-candidates", type=Path)
    args = ap.parse_args()
    if args.render_only and (args.adjudication or args.export_candidates):
        ap.error("--render-only cannot be combined with adjudication or candidate export")
    out_json = (args.output_dir / OUT_JSON.name if args.output_dir else OUT_JSON)
    out_md = (args.output_dir / OUT_MD.name if args.output_dir else OUT_MD)
    if _aliases(out_json, out_md):
        ap.error("JSON and Markdown report outputs must not alias each other")
    if args.adjudication and any(
            _aliases(args.adjudication, path) for path in (out_json, out_md)):
        ap.error("adjudication input must not alias either report output")
    if args.export_candidates and any(
            _aliases(args.export_candidates, path) for path in (out_json, out_md)):
        ap.error("candidate export cannot overwrite either report output")
    if args.render_only:
        data = json.loads(out_json.read_text(encoding="utf-8"))
        markdown = render(data)
        out_md.parent.mkdir(parents=True, exist_ok=True)
        out_md.write_text(markdown, encoding="utf-8")
    else:
        raw = scan(args.stride)
        rows = (read_adjudication(args.adjudication, raw["cohort"], DECISIONS)
                if args.adjudication else None)
        data = assemble(raw, adjudication_rows=rows)
        serialized = json.dumps(data, indent=1) + "\n"
        markdown = render(data)
        export = candidate_csv(raw["cohort"]) if args.export_candidates else None
        # Complete every read, validation and render before changing artifacts.
        out_json.parent.mkdir(parents=True, exist_ok=True)
        out_md.parent.mkdir(parents=True, exist_ok=True)
        if args.export_candidates:
            args.export_candidates.parent.mkdir(parents=True, exist_ok=True)
            args.export_candidates.write_text(export, encoding="utf-8")
        out_json.write_text(serialized, encoding="utf-8")
        out_md.write_text(markdown, encoding="utf-8")
    print(f"wrote {out_md}")
    print(f"  {data['total']} intersection articles, {data['classified']} regex candidates")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
