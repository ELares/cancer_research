#!/usr/bin/env python3
"""Which analyses still NEED the retrieved corpus, measured rather than argued.

The plan for retiring the 4,830-article corpus was: build the census analogues
first, then measure what still needs the corpus, then cut with evidence. The
analogues exist. This is the measurement.

THE RULE, stated before the result. A consumer needs the corpus only if it
reads something the census cannot supply. So the question is decided per FIELD,
not per script, against three things the census does carry:

  1. Census record fields -- pmid, title, abstract, mesh, mesh_major, mesh_ui,
     pub_types, journal, year, doi, pmcid, cancer_basis. Anything here is
     available for 4,403,994 articles instead of 4,830.
  2. The open-access full-text layer -- 1,116,481 records with `text`, keyed by
     pmid and pmcid and carrying a licence class. THIS IS THE FINDING THAT
     MOVES THE MOST CONSUMERS: "needs full text" stopped meaning "needs the
     corpus" when that layer landed, because it holds roughly 230 times more
     full text than the corpus does.
  3. NLM's own labels -- publication types and check tags, which replace this
     project's evidence tagger for study design.

What the census genuinely cannot supply is this project's OWN annotations:
`mechanisms`, `cancer_types`, `evidence_level`, `pathway_targets`,
`biology_processes`, `diagnostic_therapy_links`. A consumer reading those is
corpus-bound -- but that is not automatically a reason to keep it, because for
several the right question is whether the ANALYSIS is still wanted, not whether
the input still exists.

THE CLASSES ARE DELIBERATELY NOT "KEEP" AND "CUT". A dependency measurement
that also issued verdicts would be making an editorial decision under cover of
an arithmetic one. It reports what each consumer reads and what that implies is
possible; what to retire is a decision for a person.
"""
import argparse
import ast
import json
import re
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
SCRIPTS = REPO / "scripts"
OUT_MD = REPO / "analysis/corpus-dependency-audit.md"
OUT_JSON = REPO / "analysis/corpus-dependency-audit.json"

# Fields the census record carries, read from a real record rather than typed.
CENSUS_FIELDS_FALLBACK = {
    "pmid", "title", "abstract", "mesh", "mesh_major", "mesh_ui", "pub_types",
    "journal", "year", "doi", "pmcid", "cancer_basis", "adjacent_ui", "cancer_ui",
}
FULLTEXT_FIELDS = {"pmid", "pmcid", "licence", "text"}

# Annotations only this project produces. A consumer reading one of these
# cannot be pointed at the census without the annotation being rebuilt.
PROJECT_ANNOTATIONS = {
    "mechanisms", "cancer_types", "evidence_level", "pathway_targets",
    "biology_processes", "diagnostic_therapy_links", "evidence_signals",
}
# Full-text access in the corpus goes through these SYMBOLS, matched as names
# in the parsed tree rather than as substrings. A substring scan put three
# false consumers in the table -- including this file, which names every marker
# in its own constants, and `atlas_baseline.py`, whose only mention is a
# docstring stating that it does NOT touch the corpus. A detector that reads
# prose about a thing as use of that thing is measuring the wrong document.
FULLTEXT_SYMBOLS = {"get_searchable_text", "load_article", "extract_abstract",
                    "PMID_DIR", "ABSTRACT_PMID_DIR"}
CORPUS_MARKERS = ("INDEX.jsonl", "by-pmid", "corpus/abstracts")
# The OTHER route into the corpus, and the one the first fix broke. Most of the
# pipeline imports its paths from `config` and never writes a literal, so a
# code-only string scan excluded `tag_articles.py`, `fetch_articles.py` and
# five more -- the scripts that BUILD the corpus. Narrowing the scan to code
# fixed a false-positive class and created a false-negative one, in the
# direction that would justify deleting a real consumer. A file is a consumer
# if EITHER route fires.
CORPUS_PATH_SYMBOLS = {"PMID_DIR", "ABSTRACT_PMID_DIR", "INDEX_FILE",
                       "CORPUS_DIR", "ABSTRACTS_DIR"}
# Not consumers: this file (it names every marker as a constant) and the module
# that DEFINES the paths, which is a declaration rather than a read.
NOT_CONSUMERS = {"corpus_dependency_audit.py", "config.py"}


def census_fields() -> set:
    """Read the schema off a real census record where one is available."""
    import gzip

    rec_dir = REPO / "corpus/atlas/records"
    shards = sorted(rec_dir.glob("*.jsonl.gz")) if rec_dir.is_dir() else []
    if not shards:
        return set(CENSUS_FIELDS_FALLBACK)
    with gzip.open(shards[len(shards) // 2], "rt", encoding="utf-8") as fh:
        return set(json.loads(fh.readline()).keys())


def code_only(src: str) -> str:
    """The source with comments and docstrings removed.

    A marker inside a docstring is prose ABOUT the corpus, and at least one
    script's only mention is a line promising it leaves the corpus alone.
    """
    try:
        tree = ast.parse(src)
    except SyntaxError:
        return src
    spans = []
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.FunctionDef, ast.AsyncFunctionDef,
                             ast.ClassDef)):
            doc = ast.get_docstring(node, clean=False)
            if doc and node.body and isinstance(node.body[0], ast.Expr):
                e = node.body[0]
                if hasattr(e, "lineno") and hasattr(e, "end_lineno"):
                    spans.append((e.lineno, e.end_lineno))
    lines = src.split("\n")
    drop = set()
    for a, b in spans:
        drop.update(range(a - 1, b))
    kept = [("" if i in drop else re.sub(r"#.*$", "", ln))
            for i, ln in enumerate(lines)]
    return "\n".join(kept)


def symbols_used(src: str) -> set:
    """Every Name id and Attribute attr in the parsed tree."""
    out = set()
    try:
        tree = ast.parse(src)
    except SyntaxError:
        return out
    for node in ast.walk(tree):
        if isinstance(node, ast.Name):
            out.add(node.id)
        elif isinstance(node, ast.Attribute):
            out.add(node.attr)
    return out


def fields_read(src: str) -> set:
    """Record fields a source file reads, from subscripts and .get() calls.

    Deliberately over-inclusive: it collects every string used as a key, and
    the caller intersects against known field sets. Missing a field would
    understate a dependency, which is the direction that matters here.
    """
    out = set()
    try:
        tree = ast.parse(src)
    except SyntaxError:
        return out
    for node in ast.walk(tree):
        if isinstance(node, ast.Subscript) and isinstance(node.slice, ast.Constant) \
                and isinstance(node.slice.value, str):
            out.add(node.slice.value)
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute) \
                and node.func.attr == "get" and node.args \
                and isinstance(node.args[0], ast.Constant) \
                and isinstance(node.args[0].value, str):
            out.add(node.args[0].value)
    return out


def classify(name: str, src: str, cfields: set, route: str = "") -> dict:
    reads = fields_read(src)
    annotations = sorted(reads & PROJECT_ANNOTATIONS)
    uses_fulltext = bool(symbols_used(src) & FULLTEXT_SYMBOLS)
    census_only = sorted(reads & cfields)
    if annotations:
        cls = "needs project annotations"
        why = ("reads " + ", ".join(f"`{a}`" for a in annotations)
               + ", which the census does not carry: NLM assigns descriptors and "
                 "publication types, not this project's mechanism, cancer-type or "
                 "evidence tags")
    elif uses_fulltext:
        cls = "full text, available at census scale"
        why = ("reads corpus full text, which the open-access layer supplies for "
               "1,116,481 census records -- roughly 230 times more than the corpus "
               "holds, keyed by the same pmid")
    elif census_only:
        cls = "census-supplied fields only"
        why = ("reads only " + ", ".join(f"`{f}`" for f in census_only[:6])
               + ", every one of which the census record carries")
    else:
        cls = "no record fields resolved"
        why = ("touches the corpus but no record field was resolved statically; "
               "read it before drawing a conclusion")
    return {"script": name, "class": cls, "why": why,
            "route": route,
            "annotations": annotations,
            "uses_fulltext": uses_fulltext,
            "census_fields_read": census_only}


def scan() -> dict:
    cfields = census_fields()
    rows = []
    excluded = []
    for f in sorted(SCRIPTS.glob("*.py")):
        src = f.read_text(encoding="utf-8")
        if f.name in NOT_CONSUMERS:
            if any(m in src for m in CORPUS_MARKERS):
                excluded.append(f.name)
            continue
        # Two routes: a path literal in CODE (prose about the corpus is not use
        # of it), or a path symbol imported from config and actually used.
        by_literal = any(m in code_only(src) for m in CORPUS_MARKERS)
        by_symbol = bool(symbols_used(src) & CORPUS_PATH_SYMBOLS)
        if not (by_literal or by_symbol):
            if any(m in src for m in CORPUS_MARKERS):
                excluded.append(f.name)
            continue
        route = ("path literal" if by_literal else "") + \
                (" + " if by_literal and by_symbol else "") + \
                ("config path symbol" if by_symbol else "")
        rows.append(classify(f.name, src, cfields, route))
    return {
        "excluded_mentions_only": sorted(excluded),
        "census_fields": sorted(cfields),
        "fulltext_fields": sorted(FULLTEXT_FIELDS),
        "project_annotations": sorted(PROJECT_ANNOTATIONS),
        "consumers": rows,
    }


def assemble(d: dict) -> dict:
    counts = {}
    for r in d["consumers"]:
        counts[r["class"]] = counts.get(r["class"], 0) + 1
    out = dict(d)
    out["by_class"] = counts
    out["n_consumers"] = len(d["consumers"])
    out["corpus_bound"] = sorted(r["script"] for r in d["consumers"]
                                 if r["annotations"])
    return out


def render(d: dict) -> str:
    L = ["# What still needs the retrieved corpus\n"]
    L.append(
        f"Generated by `scripts/corpus_dependency_audit.py`. "
        f"{d['n_consumers']} scripts read `corpus/INDEX.jsonl`, "
        f"`corpus/by-pmid/` or `corpus/abstracts/`.\n"
    )
    L.append(
        "The rule is applied per FIELD, not per script: a consumer needs the "
        "corpus only where it reads something the census cannot supply. The "
        "census record carries "
        + ", ".join(f"`{f}`" for f in d["census_fields"])
        + "; the open-access layer supplies full text for 1,116,481 records; "
        "and NLM supplies publication types and check tags.\n"
    )
    L.append("| class | scripts |")
    L.append("|---|--:|")
    for k, v in sorted(d["by_class"].items(), key=lambda kv: -kv[1]):
        L.append(f"| {k} | {v} |")
    L.append("")
    L.append(
        f"**{len(d['corpus_bound'])} of {d['n_consumers']} are genuinely "
        f"corpus-bound**, and all for the same reason: they read annotations "
        f"only this project produces. Everything else reads fields the census "
        f"carries, or full text the census's own open-access layer carries in "
        f"far greater quantity.\n"
    )
    L.append(
        "That is the finding worth acting on. \"Needs full text\" stopped "
        "meaning \"needs the corpus\" when the open-access layer landed, and "
        "several consumers in the middle class were kept on the corpus for a "
        "reason that expired.\n"
    )
    for cls in sorted(d["by_class"], key=lambda k: -d["by_class"][k]):
        L.append(f"## {cls}\n")
        for r in d["consumers"]:
            if r["class"] == cls:
                L.append(f"- `{r['script']}` — {r['why']}")
        L.append("")
    ex = d.get("excluded_mentions_only") or []
    if ex:
        L.append("## Named the corpus without reading it\n")
        L.append(
            f"{len(ex)} file(s) mention a corpus path only in a docstring, a "
            f"comment, or as a path DEFINITION, and are not consumers: "
            + ", ".join(f"`{x}`" for x in ex)
            + ". A first version of this audit counted them, which put "
            "`atlas_baseline.py` in the table on the strength of a line "
            "promising it leaves the corpus alone -- and put this file in it "
            "too, since it names every marker as a constant.\n"
        )
    L.append("## What this does not decide\n")
    L.append(
        "Nothing here says an analysis should be kept or cut. A dependency "
        "measurement that also issued verdicts would be making an editorial "
        "decision under cover of an arithmetic one. For several corpus-bound "
        "consumers the real question is whether the ANALYSIS is still wanted "
        "-- a gold-set evaluation of a tagger no manuscript figure now depends "
        "on is a different case from one that still calibrates something -- and "
        "that question is not answerable by reading imports.\n"
    )
    L.append(
        "The static read is also over-inclusive by design: it collects every "
        "string used as a subscript or `.get()` key and intersects against "
        "known field sets, so a field accessed through a variable is missed. "
        "Missing one UNDERSTATES a dependency, which is the safer direction "
        "for a document that could be used to justify deleting something.\n"
    )
    return "\n".join(L)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--render-only", action="store_true")
    a = ap.parse_args()
    if a.render_only:
        d = assemble(json.loads(OUT_JSON.read_text()))
    else:
        d = assemble(scan())
    OUT_JSON.write_text(json.dumps(d, indent=1) + "\n")
    OUT_MD.write_text(render(d))
    print(f"wrote {OUT_MD}")
    for k, v in sorted(d["by_class"].items(), key=lambda kv: -kv[1]):
        print(f"  {v:3d}  {k}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
