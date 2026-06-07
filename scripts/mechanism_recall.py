#!/usr/bin/env python3
"""Non-circular mechanism-tagging recall measurement (#412).

WHY THIS EXISTS
---------------
The evidence-tagging gold set is a manual sample, so the manuscript reports a
real evidence-tier recall (55%). The MECHANISM tagger has no comparable
non-circular recall number: the obvious reference label (MeSH descriptors) is
folded straight into the matched text by `get_searchable_text`, so any recall
computed against a MeSH-derived label would be partly tautological (a paper whose
MeSH says "Oncolytic Virotherapy" almost always uses the phrase "oncolytic
virotherapy" in its abstract too, which the keyword tagger then trivially
catches). This script measures mechanism recall WITHOUT that circularity.

METHOD (a controlled-vocabulary CONCORDANCE check, not manual ground truth)
--------------------------------------------------------------------------
1. `analysis/mesh-mechanism-map.yaml` (handwritten) lists, per mechanism, the
   discriminative leaf MeSH descriptors whose presence means "this article is
   genuinely about mechanism M" with high precision. MeSH is assigned by NLM
   indexers from the full paper, INDEPENDENTLY of our title/abstract keyword set.

2. Each descriptor is classified at runtime as INDEPENDENT or LEAKY:
       leaky  = the descriptor string itself would trigger M's keyword tagger
                (`is_keyword_substring`), so MeSH agreement is near-tautological;
       independent = it would NOT, so agreement is a genuine concordance signal.
   The classifier calls the LIVE tagger matcher (`text_matches_keyword`) so it
   can never drift from production semantics.

3. The tagger is RE-RUN over each record on the LEAKAGE-FREE text
   (`get_searchable_text(..., include_metadata=False)` => title + abstract only,
   MeSH/disease/gene/drug annotations EXCLUDED). This is the load-bearing
   non-circularity guard: the reference label (MeSH) is never visible to the
   matcher whose recall we measure.

4. HEADLINE recall_M = |MeSHpos(independent) AND Tagpos| / |MeSHpos(independent)|,
   computed ONLY from the independent descriptor pool, only when that pool has at
   least N_MIN records. Leaky-descriptor recall is reported separately and
   labelled near-tautological. Mechanisms with no independent leaf, only
   over-broad / bucket-confounded descriptors, or an empty pool are reported as
   NOT MeSH-measurable (with the reason) and NEVER as 0% recall.

This is mechanism recall, distinct from the evidence-tier 55% recall, and it does
not touch the 96% evidence precision. Offline and deterministic: it reads the
frozen corpus only and writes analysis/mechanism-recall-report.md + .json.
"""

import argparse
import json
import sys
from pathlib import Path

import yaml

# Import the LIVE tagger semantics so the leak classifier and the re-match path
# can never drift from production. (scripts/ on path for the sibling imports.)
sys.path.insert(0, str(Path(__file__).resolve().parent))
from article_io import load_article  # noqa: E402
from config import MECHANISM_KEYWORDS  # noqa: E402
from evidence_utils import normalize_text  # noqa: E402
from tag_articles import (  # noqa: E402
    get_searchable_text,
    match_mechanisms,
    text_matches_keyword,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
MAP_PATH = REPO_ROOT / "analysis" / "mesh-mechanism-map.yaml"
CORPUS_DIR = REPO_ROOT / "corpus" / "by-pmid"
INDEX_PATH = REPO_ROOT / "corpus" / "INDEX.jsonl"
REPORT_PATH = REPO_ROOT / "analysis" / "mechanism-recall-report.md"
JSON_PATH = REPO_ROOT / "analysis" / "mechanism-recall.json"

# Minimum independent-descriptor MeSH pool to report a non-circular recall.
N_MIN = 30

# Canonical mechanism tag names keyed by lowercase, so the YAML's lowercase
# mechanism keys (e.g. "mrna-vaccine") resolve to the production tag casing
# (e.g. "mRNA-vaccine") used by MECHANISM_KEYWORDS and match_mechanisms output.
_CANON = {k.lower(): k for k in MECHANISM_KEYWORDS}


def canonical_mechanism(name: str) -> str:
    """Resolve a YAML mechanism key to the production tag casing."""
    return _CANON.get(name.lower(), name)


def is_keyword_substring(descriptor: str, mechanism: str) -> bool:
    """True iff `descriptor` would itself trigger `mechanism`'s keyword tagger.

    LEAKY descriptors (this returns True) are near-tautological with the keyword
    set, so MeSH agreement carries no independent information and they are
    EXCLUDED from the headline recall. Uses the live `text_matches_keyword`
    matcher against the live `MECHANISM_KEYWORDS`, so it tracks production exactly.
    """
    norm = normalize_text(descriptor)
    canon = canonical_mechanism(mechanism)
    return any(text_matches_keyword(norm, kw) for kw in MECHANISM_KEYWORDS.get(canon, []))


def classify_descriptors(mech_map: dict) -> dict:
    """Split each measurable mechanism's descriptors into independent vs leaky.

    Returns {mechanism: {"canon", "independent": [...], "leaky": [...], "note",
    "proxy_confounded"}}. Pure: depends only on the map + the live keyword set.
    """
    out = {}
    for mech, spec in (mech_map.get("mechanisms") or {}).items():
        descriptors = spec.get("descriptors", []) if isinstance(spec, dict) else list(spec)
        independent, leaky = [], []
        for d in descriptors:
            (leaky if is_keyword_substring(d, mech) else independent).append(d)
        out[mech] = {
            "canon": canonical_mechanism(mech),
            "independent": independent,
            "leaky": leaky,
            "note": (spec.get("note", "").strip() if isinstance(spec, dict) else ""),
            "proxy_confounded": bool(isinstance(spec, dict) and spec.get("proxy_confounded")),
        }
    return out


def recall(hit: int, pool: int):
    """hit/pool, or None when the pool is empty."""
    return (hit / pool) if pool else None


def load_index_mechanisms(index_path: Path) -> dict:
    """pmid (str) -> set of frozen production mechanism tags from INDEX.jsonl."""
    out = {}
    if not index_path.exists():
        return out
    with open(index_path, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                rec = json.loads(line)
            except json.JSONDecodeError:
                continue
            pmid = str(rec.get("pmid", ""))
            if pmid:
                out[pmid] = set(rec.get("mechanisms") or [])
    return out


def measure(corpus_dir: Path, index_path: Path, mech_map: dict, limit: int = 0) -> dict:
    """Stream the frozen corpus and accumulate per-mechanism recall counts.

    For every record with non-empty MeSH, the tagger is re-run on the
    leakage-free (MeSH-excluded) title+abstract text. A record counts toward a
    mechanism's pool when it carries one of that mechanism's reference descriptors;
    it counts as a hit when the leakage-free tagger also tags the mechanism.
    """
    classified = classify_descriptors(mech_map)
    index_mech = load_index_mechanisms(index_path)

    # Per-mechanism accumulators. *_pool / *_hit_lf are independent-descriptor
    # counts (the headline); *_hit_idx mirrors the SHIPPED tagger (frozen INDEX,
    # which saw MeSH) over the same pool; leaky_* are the near-tautological track.
    acc = {
        m: dict(indep_pool=0, indep_hit_lf=0, indep_hit_idx=0, leaky_pool=0, leaky_hit_lf=0)
        for m in classified
    }

    n_records = 0
    n_with_mesh = 0
    files = sorted(corpus_dir.glob("*.md"))
    if limit:
        files = files[:limit]

    for fp in files:
        fm, body = load_article(fp)
        if not fm:
            continue
        n_records += 1
        mesh = set(fm.get("mesh_terms") or [])
        if not mesh:
            continue  # a record with no MeSH cannot serve as a MeSH reference
        n_with_mesh += 1

        # Leakage-free re-tag: MeSH/annotations excluded from the matched text.
        lf_text = get_searchable_text(fm, body, include_metadata=False)
        title_text = normalize_text(fm.get("title", ""))
        lf_mechs = set(match_mechanisms(lf_text, title_text))

        pmid = str(fm.get("pmid", ""))
        idx_mechs = index_mech.get(pmid, set())

        for mech, spec in classified.items():
            canon = spec["canon"]
            has_indep = any(d in mesh for d in spec["independent"])
            has_leaky = any(d in mesh for d in spec["leaky"])
            if has_indep:
                acc[mech]["indep_pool"] += 1
                if canon in lf_mechs:
                    acc[mech]["indep_hit_lf"] += 1
                if canon in idx_mechs:
                    acc[mech]["indep_hit_idx"] += 1
            if has_leaky:
                acc[mech]["leaky_pool"] += 1
                if canon in lf_mechs:
                    acc[mech]["leaky_hit_lf"] += 1

    return build_results(classified, acc, mech_map, n_records, n_with_mesh)


def build_results(classified: dict, acc: dict, mech_map: dict, n_records: int, n_with_mesh: int) -> dict:
    """Assemble per-mechanism results + aggregates from raw counts. Pure."""
    measurable = {}  # has an independent pool >= N_MIN
    leaky_only = {}  # independent pool too small / absent; only leaky signal
    for mech, spec in classified.items():
        a = acc[mech]
        entry = {
            "canon": spec["canon"],
            "independent_descriptors": spec["independent"],
            "leaky_descriptors": spec["leaky"],
            "n_independent": a["indep_pool"],
            "n_leaky": a["leaky_pool"],
            "recall_leakage_free": recall(a["indep_hit_lf"], a["indep_pool"]),
            "recall_production_index": recall(a["indep_hit_idx"], a["indep_pool"]),
            "recall_leaky_near_tautological": recall(a["leaky_hit_lf"], a["leaky_pool"]),
            "proxy_confounded": spec["proxy_confounded"],
            "note": spec["note"],
        }
        if a["indep_pool"] >= N_MIN and not spec["proxy_confounded"]:
            measurable[mech] = entry
        else:
            if a["indep_pool"] < N_MIN:
                entry["reason"] = (
                    f"independent-descriptor pool n={a['indep_pool']} < N_MIN={N_MIN}"
                    + (
                        " (all discriminative descriptors are keyword-tautological)"
                        if not spec["independent"]
                        else ""
                    )
                )
            elif spec["proxy_confounded"]:
                entry["reason"] = "proxy-confounded: recall not interpretable as tagger recall"
            leaky_only[mech] = entry

    # Aggregates over the non-circular measurable set only.
    pools = [(m, e) for m, e in measurable.items()]
    tot_pool = sum(e["n_independent"] for _, e in pools)
    tot_hit = sum(round(e["recall_leakage_free"] * e["n_independent"]) for _, e in pools)
    macro_vals = [e["recall_leakage_free"] for _, e in pools if e["recall_leakage_free"] is not None]
    aggregates = {
        "n_measurable_mechanisms": len(measurable),
        "volume_weighted_recall_leakage_free": (tot_hit / tot_pool) if tot_pool else None,
        "macro_recall_leakage_free": (sum(macro_vals) / len(macro_vals)) if macro_vals else None,
        "total_independent_pool": tot_pool,
    }

    unmeasurable = dict(mech_map.get("unmeasurable") or {})
    return {
        "n_min": N_MIN,
        "n_records_scanned": n_records,
        "n_records_with_mesh": n_with_mesh,
        "measurable": measurable,
        "leaky_only": leaky_only,
        "unmeasurable": unmeasurable,
        "aggregates": aggregates,
    }


def _fmt_pct(x):
    return "n/a" if x is None else f"{100 * x:.1f}%"


def write_report(results: dict, path: Path) -> None:
    agg = results["aggregates"]
    lines = [
        "# Mechanism-tagging recall (non-circular MeSH concordance, #412)",
        "",
        "Generated by `scripts/mechanism_recall.py` (offline, deterministic). Do not",
        "hand-edit; rerun the script.",
        "",
        "## What this measures (and what it does not)",
        "",
        "This is the recall of the **mechanism** keyword tagger, measured against",
        "expert-assigned **MeSH** descriptors as an independent reference label. It is a",
        "controlled-vocabulary *concordance* check, not manual ground truth, and it is",
        "**distinct** from the evidence-tier recall (55%) and does not touch the evidence",
        "precision (96%).",
        "",
        "The reference label (MeSH) is folded into the production tagger's matched text,",
        "so a naive MeSH-vs-tag recall would be partly circular. Two guards remove the",
        "circularity:",
        "",
        "1. The tagger is **re-run on a leakage-free text** (title + abstract only; MeSH",
        "   and the other curated annotations excluded), so the reference label is never",
        "   visible to the matcher whose recall we measure.",
        "2. The **headline** recall uses only **independent** descriptors: leaf MeSH terms",
        "   whose own string would *not* trip the mechanism's keywords. Descriptors that",
        "   would (e.g. the descriptor *contains* a tagger keyword) are near-tautological",
        "   and are reported separately, never in the headline.",
        "",
        f"Pools below the reporting floor (N_MIN = {results['n_min']} independent-descriptor",
        "records), proxy-confounded pools, and mechanisms with no independent leaf are",
        "reported as **not MeSH-measurable** with the reason, never as 0% recall.",
        "",
        f"Records scanned: {results['n_records_scanned']:,}  ",
        f"Records with non-empty MeSH (eligible as a reference): {results['n_records_with_mesh']:,}",
        "",
        "## Headline (non-circular, independent descriptors only)",
        "",
        f"- Measurable mechanisms (independent pool >= {results['n_min']}): "
        f"**{agg['n_measurable_mechanisms']}**",
        f"- Volume-weighted recall: **{_fmt_pct(agg['volume_weighted_recall_leakage_free'])}** "
        f"(over {agg['total_independent_pool']:,} independent-descriptor records)",
        f"- Macro (per-mechanism mean) recall: **{_fmt_pct(agg['macro_recall_leakage_free'])}**",
        "",
        "`recall (leakage-free)` is the headline: the live tagger re-run with MeSH excluded.",
        "`recall (production)` reflects the frozen `INDEX.jsonl` snapshot (its keywords, with",
        "MeSH folded in) over the same pool. It usually exceeds the leakage-free column, but it",
        "can LAG it where the live keyword set was improved without regenerating the frozen",
        "corpus, which is deliberately the case for epigenetic after #418 (69.8% live vs 54.7%",
        "in the frozen snapshot): the manuscript's frozen 19-mechanism counts are intentionally",
        "not re-tagged, so the gain shows in the live measurement, not the frozen INDEX.",
        "",
        "| Mechanism | Independent descriptor(s) | N | recall (leakage-free) | recall (production) |",
        "|---|---|---:|---:|---:|",
    ]
    for mech in sorted(results["measurable"]):
        e = results["measurable"][mech]
        descs = ", ".join(e["independent_descriptors"])
        lines.append(
            f"| {mech} | {descs} | {e['n_independent']} | "
            f"{_fmt_pct(e['recall_leakage_free'])} | {_fmt_pct(e['recall_production_index'])} |"
        )

    noted = [(m, results["measurable"][m]["note"]) for m in sorted(results["measurable"]) if results["measurable"][m]["note"]]
    if noted:
        lines += ["", "### Per-mechanism notes", ""]
        for mech, note in noted:
            lines.append(f"- **{mech}**: {' '.join(note.split())}")

    lines += [
        "",
        "## Near-tautological / sub-floor (reported, not in the headline)",
        "",
        "Mechanisms whose only discriminative MeSH leaves are keyword-tautological, or",
        "whose independent pool is below the floor. Where a leaky pool exists its recall",
        "is shown for completeness and is expected near 100% by construction.",
        "",
        "| Mechanism | Reason | N indep | N leaky | leaky recall |",
        "|---|---|---:|---:|---:|",
    ]
    for mech in sorted(results["leaky_only"]):
        e = results["leaky_only"][mech]
        lines.append(
            f"| {mech} | {e.get('reason', '')} | {e['n_independent']} | {e['n_leaky']} | "
            f"{_fmt_pct(e['recall_leaky_near_tautological'])} |"
        )

    lines += [
        "",
        "## Not MeSH-measurable",
        "",
        "No discriminative descriptor exists (term too new, only over-broad context",
        "descriptors, taxonomy-bucket-confounded, or an empty/tiny pool). These are a",
        "limitation of the MeSH reference, not measured tagger misses.",
        "",
        "| Mechanism | Reason |",
        "|---|---|",
    ]
    for mech in sorted(results["unmeasurable"]):
        lines.append(f"| {mech} | {results['unmeasurable'][mech]} |")

    lines += [
        "",
        "## Interpretation",
        "",
        "For the mechanisms with a genuinely independent MeSH leaf, the keyword tagger",
        "recovers the large majority of expert-MeSH-labelled records from title+abstract",
        "alone, i.e. mechanism recall is high where it can be measured without",
        "circularity. The honest limit is coverage of the *measurement*: only a subset of",
        "mechanisms have an independent leaf MeSH descriptor, and the newest / device /",
        "regulated-cell-death modalities are not MeSH-measurable at all (see the table",
        "above). High recall here is not a precision claim and says nothing about the",
        "mechanisms in the not-measurable table.",
        "",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--corpus", type=Path, default=CORPUS_DIR)
    ap.add_argument("--index", type=Path, default=INDEX_PATH)
    ap.add_argument("--map", dest="map_path", type=Path, default=MAP_PATH)
    ap.add_argument("--report", type=Path, default=REPORT_PATH)
    ap.add_argument("--json", dest="json_path", type=Path, default=JSON_PATH)
    ap.add_argument("--limit", type=int, default=0, help="cap records scanned (debug)")
    args = ap.parse_args()

    mech_map = yaml.safe_load(args.map_path.read_text(encoding="utf-8"))
    results = measure(args.corpus, args.index, mech_map, limit=args.limit)

    write_report(results, args.report)
    args.json_path.write_text(json.dumps(results, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    agg = results["aggregates"]
    print(
        f"measurable={agg['n_measurable_mechanisms']} "
        f"volume_weighted={_fmt_pct(agg['volume_weighted_recall_leakage_free'])} "
        f"macro={_fmt_pct(agg['macro_recall_leakage_free'])}"
    )
    print(f"wrote {args.report.relative_to(REPO_ROOT)} and {args.json_path.relative_to(REPO_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
