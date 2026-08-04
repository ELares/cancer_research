#!/usr/bin/env python3
"""Atlas: measure entity ambiguity and stop the graph guessing (#ATLAS-AMBIG).

WHY
---
`atlas_entity_audit.py` found one bad symbol the hard way: PubTator3 maps
`FSP1` to NCBI Gene 51062 = ATL1 (atlastin GTPase 1), not 84883 = AIFM2, the
ferroptosis suppressor the manuscript's headline synergy claim rests on. That
audit checked 53 hand-listed symbols and reported the mismatch. It did not ask
the obvious next question: **how many other FSP1s are in this graph?**

This script asks it. The answer changes how the graph may be queried.

The mechanism is in `atlas_graph.build_index`: a surface form's identifier is
`alias[s].most_common(1)[0][0]` -- a majority vote. When a symbol is genuinely
ambiguous the vote still returns one identifier, silently, with no signal that
anything was discarded. `resolve("FSP1")` returns atlastin and looks exactly
like a correct answer.

WHAT IT MEASURES
----------------
All three entity types, because the answer differs sharply between them: genes
sit on a contested form for 28.3% of mentions, chemicals 2.0%, diseases 2.7%.
Genes are roughly ten times worse, and structurally so -- MeSH is a curated
vocabulary with one preferred term per concept, while NCBI itself lists `FSP1`
as an official alias of three different genes.

Contested forms fall into three classes that must not be conflated, because
only one is damaging:

* **Species ambiguity** -- `GAPDH` maps to human 2597 and mouse 14433. Same
  gene, different organism. For a literature map this is benign and often
  desirable.
* **Hierarchical granularity** -- `GBM` maps to Glioblastoma and to Glioma, one
  nested under the other in the MeSH tree. Merging loses specificity but does
  not attribute the biology to an unrelated concept.
* **Sense collision** -- `ER` maps to EREG (epiregulin) and ESR1 (estrogen
  receptor); `PC` maps to Prostatic and Pancreatic Neoplasms. Unrelated
  concepts. Collapsing them corrupts the edge.

Reporting the undecomposed contested share as an error rate would overstate the
damage several-fold, so this script never reports it alone.

WHAT IT CHANGES
---------------
Sense-colliding symbols go into a committed blocklist. `atlas_graph.resolve`
consults it and returns `None` instead of a plausible-looking wrong gene, so an
analysis fails loudly rather than quietly attributing ferroptosis biology to a
spastic-paraplegia gene. Where the cancer-domain sense is genuinely
unambiguous, `DOMAIN_SENSE` records it with its justification; a caller opts in
explicitly. `FSP1` is deliberately NOT given a domain default -- it is
genuinely context-dependent (see `atlas_disambiguate.py`).

Requires network (NCBI E-utilities for genes, NLM's MeSH identifier
service for chemicals and diseases) to name the identifiers. Derived artifacts
are committed so every downstream consumer stays offline.

Usage:
    python scripts/atlas_ambiguity.py
    python scripts/atlas_ambiguity.py --top 400 --min-share 0.20
"""

import argparse
import collections
import gzip
import json
import sys
import time
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from atlas_baseline import atlas_root  # noqa: E402
from config import PROJECT_ROOT  # noqa: E402

OUT = PROJECT_ROOT / "analysis" / "atlas-ambiguity.md"
RAW = PROJECT_ROOT / "analysis" / "atlas-ambiguity.json"

# A surface form must appear this often before its ambiguity is worth reporting.
MIN_MENTIONS = 3
# The runner-up identifier must hold at least this share for the form to count
# as genuinely contested rather than carrying a long tail of stray annotations.
MIN_RUNNER_UP_SHARE = 0.20

# Cancer-domain senses that are not in genuine doubt, each with the reason it is
# safe. A caller opts in explicitly via resolve(..., allow_domain_sense=True);
# nothing applies these silently. Symbols whose correct sense depends on the
# individual paper -- FSP1 above all -- are deliberately absent.
DOMAIN_SENSE = {
    # Every share below is MEASURED, not asserted: scripts/atlas_domain_sense.py
    # samples cancer papers carrying the symbol and counts how many declare each
    # sense in their own words. See analysis/atlas-domain-sense-validation.md.
    "er": ("2099", "ESR1", "In the oncology literature 'ER' is the estrogen "
                           "receptor (ER-positive/ER-negative disease): 98.7% of "
                           "451 declaring papers, with epiregulin declared zero "
                           "times. The majority vote returns EREG (epiregulin)."),
    "psa": ("354", "KLK3", "'PSA' in cancer is prostate-specific antigen, KLK3: "
                           "528 of 528 declaring papers. Here the majority vote "
                           "ALREADY returns KLK3, so the form stays blocklisted "
                           "only because NPEPPS is a genuinely different gene. An "
                           "earlier version of this note claimed the vote returned "
                           "NPEPPS; that was asserted rather than measured, and "
                           "was wrong."),
    "cox-2": ("5743", "PTGS2", "'COX-2' in cancer pharmacology is "
                               "prostaglandin-endoperoxide synthase 2, the "
                               "celecoxib target: 541 of 541 declaring papers, "
                               "with the mitochondrial cytochrome c oxidase "
                               "subunit II declared zero times. The majority vote "
                               "returns the mitochondrial one."),
    "p62": ("8878", "SQSTM1", "'p62' in autophagy and redox signalling is "
                              "sequestosome-1: 97.4% of 272 declaring papers. The "
                              "majority vote returns NUP62, a nuclear pore "
                              "protein."),
    "p21": ("1026", "CDKN1A", "'p21' as a cell-cycle inhibitor is CDKN1A: 89.6% "
                              "of 336 declaring papers, the lowest share of the "
                              "five because RAS p21 and a histone pseudogene are "
                              "both real minority senses. The majority vote "
                              "returns H3P16, the pseudogene."),
}


# Symbols this project queries. Checked for ambiguity whatever their volume
# rank, because relevance to this repo and mention count in the wider
# literature are unrelated -- see classify().
try:
    from atlas_entity_audit import SYMBOLS as ALWAYS_CHECK
except Exception:  # keep the scan runnable if the audit module moves
    ALWAYS_CHECK = ["GPX4", "FSP1", "AIFM2", "ACSL4", "SLC7A11", "NCOA4"]


def ncbi_gene(ids: list) -> dict:
    """id -> (symbol, organism), batched against E-utilities."""
    out = {}
    for i in range(0, len(ids), 180):
        batch = ids[i:i + 180]
        url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
               "?db=gene&retmode=json&id=" + ",".join(batch))
        try:
            res = json.load(urllib.request.urlopen(url, timeout=120))["result"]
        except Exception as exc:  # network is best-effort; unresolved is reported
            print(f"  ! esummary batch failed: {exc}", file=sys.stderr)
            continue
        for g in batch:
            e = res.get(g, {})
            if e:
                out[g] = (e.get("name", "?"),
                          e.get("organism", {}).get("scientificname", "?"))
        time.sleep(0.4)
    return out


def mesh_label_and_tree(ui: str) -> tuple:
    """MeSH descriptor UI -> (label, tree numbers), via NLM's identifier service."""
    try:
        r = json.load(urllib.request.urlopen(
            f"https://id.nlm.nih.gov/mesh/{ui}.json", timeout=30))
    except Exception:
        return ("?", [])
    label = r.get("label", {})
    label = label.get("@value") if isinstance(label, dict) else label
    trees = r.get("treeNumber", [])
    if isinstance(trees, (str, dict)):
        trees = [trees]
    out = []
    for t in trees:
        t = t.get("@id", "") if isinstance(t, dict) else str(t)
        out.append(t.rsplit("/", 1)[-1])
    return (label or "?", out)


def hierarchically_related(t1: list, t2: list) -> bool:
    """True when one descriptor sits under the other in the MeSH tree.

    Glioblastoma (C04.557.470.670.200) is a descendant of Glioma
    (C04.557.470.670): merging them loses specificity but does not attribute
    the biology to an unrelated concept, which is a different and much milder
    failure than a genuine sense collision.
    """
    return any(a != b and (a.startswith(b + ".") or b.startswith(a + "."))
               for a in t1 for b in t2)


def scan(path: Path, min_mentions: int, min_share: float) -> tuple:
    """Surface form -> identifier distribution, and the contested subset."""
    alias = collections.defaultdict(collections.Counter)
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            gid, mention = parts[2], parts[3]
            # unmapped placeholders are not a competing sense; build_index
            # skips them too, and counting them inflates the contested rate
            if gid in ("-", "None", ""):
                continue
            # PubTator packs co-referent surface forms into one pipe-joined
            # field; each is a distinct way the paper named the entity.
            for surface in mention.split("|"):
                surface = surface.strip().lower()
                if surface:
                    alias[surface][gid] += 1

    kept = {a: c for a, c in alias.items() if sum(c.values()) >= min_mentions}
    contested = {}
    for a, c in kept.items():
        if len(c) < 2:
            continue
        total = sum(c.values())
        runner_up = c.most_common(2)[1][1]
        if runner_up / total >= min_share:
            contested[a] = {"total": total, "dist": dict(c.most_common(6))}
    return kept, contested


def classify(contested: dict, top: int) -> tuple:
    """Split the contested forms into species vs sense.

    Volume rank alone is not enough. `FSP1` -- the collision that motivated this
    whole layer -- ranks 2181st by mention count and would fall outside any
    sane top-N, so the symbols this project actually queries are always checked
    regardless of how often the wider literature mentions them.
    """
    ranked = sorted(contested.items(), key=lambda kv: -kv[1]["total"])[:top]
    seen = {a for a, _ in ranked}
    for sym in ALWAYS_CHECK:
        key = sym.lower()
        if key in contested and key not in seen:
            ranked.append((key, contested[key]))
            seen.add(key)
    ids = set()
    for _a, d in ranked:
        for g in sorted(d["dist"], key=lambda k: -d["dist"][k])[:2]:
            ids.add(g)
    print(f"  resolving {len(ids):,} identifiers against NCBI ...", flush=True)
    info = ncbi_gene(sorted(ids))

    species, sense, unresolved = [], [], []
    for a, d in ranked:
        pair = sorted(d["dist"], key=lambda k: -d["dist"][k])[:2]
        n1, o1 = info.get(pair[0], ("?", "?"))
        n2, o2 = info.get(pair[1], ("?", "?"))
        row = {"surface": a, "total": d["total"],
               "top": {"id": pair[0], "symbol": n1, "organism": o1},
               "runner_up": {"id": pair[1], "symbol": n2, "organism": o2}}
        if "?" in (n1, n2):
            unresolved.append(row)
        elif n1.upper() == n2.upper():
            species.append(row)
        else:
            sense.append(row)
    return species, sense, unresolved


def classify_mesh(contested: dict, top: int) -> tuple:
    """Split contested MeSH forms into hierarchical vs genuine sense collision."""
    ranked = sorted(contested.items(), key=lambda kv: -kv[1]["total"])[:top]
    uis = set()
    for _a, d in ranked:
        for m in sorted(d["dist"], key=lambda k: -d["dist"][k])[:2]:
            uis.add(m.split(":")[-1])
    print(f"  resolving {len(uis):,} MeSH descriptors against NLM ...", flush=True)
    info = {}
    for ui in sorted(uis):
        info[ui] = mesh_label_and_tree(ui)
        time.sleep(0.15)

    hier, sense, unresolved = [], [], []
    for a, d in ranked:
        pair = sorted(d["dist"], key=lambda k: -d["dist"][k])[:2]
        u1, u2 = pair[0].split(":")[-1], pair[1].split(":")[-1]
        (l1, t1), (l2, t2) = info.get(u1, ("?", [])), info.get(u2, ("?", []))
        row = {"surface": a, "total": d["total"],
               "top": {"id": pair[0], "label": l1},
               "runner_up": {"id": pair[1], "label": l2}}
        if "?" in (l1, l2) or not (t1 and t2):
            unresolved.append(row)
        elif hierarchically_related(t1, t2):
            hier.append(row)
        else:
            sense.append(row)
    return hier, sense, unresolved




def render(results: dict, top: int) -> str:
    g = results["gene"]
    lines = [
        "# Entity ambiguity in the census graph (#ATLAS-AMBIG)", "",
        "Generated by `scripts/atlas_ambiguity.py`. The entity audit found one bad",
        "symbol (`FSP1` -> atlastin). This measures how many others there are, across",
        "all three entity types.", "",
        "## Why a majority vote hides this", "",
        "`atlas_graph.build_index` resolves a surface form with",
        "`alias[s].most_common(1)[0][0]`. When a form is genuinely ambiguous the vote",
        "still returns exactly one identifier, with no signal that anything was",
        "discarded, so a wrong answer is indistinguishable from a right one.", "",
        "## Genes are the outlier", "",
        "| entity type | surface forms | contested | share of mentions on a contested form |",
        "|---|---|---|---|",
    ]
    for kind in ("gene", "chemical", "disease"):
        r = results[kind]
        lines.append(
            f"| {kind} | {r['forms']:,} | {r['contested']:,} "
            f"({100*r['contested']/max(1,r['forms']):.1f}%) | {r['mention_share']:.1f}% |")
    lines += [
        "", "Genes are roughly **ten times** worse than chemicals or diseases by",
        "mention share, and the reason is structural rather than accidental. MeSH is",
        "a curated controlled vocabulary with one preferred term per concept. Gene",
        "symbols are not: NCBI itself lists `FSP1` as an official alias of three",
        "different genes, so the ambiguity is in the reference data, not only in the",
        "extraction.", "",
        "## The decomposition that matters", "",
        "The contested share is **not** an error rate, and reporting it as one would",
        "overstate the damage several-fold. Contested forms fall into three classes",
        "with very different consequences:", "",
        "| class | what it means | consequence |", "|---|---|---|",
        "| species ambiguity | human vs mouse *GAPDH* -- the same gene | benign; a literature map usually wants them merged |",
        "| hierarchical granularity | *Glioblastoma* vs *Glioma* -- one nests under the other | lossy, not wrong; specificity is lost, biology is not misattributed |",
        "| sense collision | *EREG* vs *ESR1* for `ER` -- unrelated concepts | damaging; the edge is attributed to the wrong biology |",
        "",
        f"Measured on the {top} highest-volume contested forms per type, resolved",
        "against NCBI and NLM:", "",
        "| entity type | benign class | sense collision | unresolved |",
        "|---|---|---|---|",
    ]
    for kind, benign_name in (("gene", "species"), ("chemical", "hierarchical"),
                              ("disease", "hierarchical")):
        r = results[kind]
        d = r["benign"] + r["sense"]
        lines.append(
            f"| {kind} | {r['benign']} {benign_name} ({100*r['benign']/max(1,d):.1f}%) "
            f"| {r['sense']} ({100*r['sense']/max(1,d):.1f}%) | {r['unresolved']} |")

    lines += ["", "## The damaging class, by mention volume", ""]
    for kind in ("gene", "chemical", "disease"):
        rows = results[kind]["sense_rows"]
        if not rows:
            continue
        lines += [f"### {kind}", "",
                  "| surface form | mentions | majority vote returns | runner-up |",
                  "|---|---|---|---|"]
        for r in sorted(rows, key=lambda r: -r["total"])[:15]:
            t, u = r["top"], r["runner_up"]
            tl = t.get("symbol") or t.get("label", "?")
            ul = u.get("symbol") or u.get("label", "?")
            to = f" ({t['organism']})" if t.get("organism") else ""
            uo = f" ({u['organism']})" if u.get("organism") else ""
            lines.append(f"| `{r['surface']}` | {r['total']:,} | {tl}{to} | {ul}{uo} |")
        lines.append("")

    lines += [
        "Four are worth stating plainly, because they touch the most-queried",
        "concepts in the cancer literature:", "",
        "* `ER` resolves to **EREG (epiregulin)**, not ESR1. In oncology `ER` is the",
        "  estrogen receptor -- ER-positive and ER-negative disease is one of the",
        "  field's primary axes.",
        "* `COX-2` resolves to **mitochondrial cytochrome c oxidase subunit II**, not",
        "  PTGS2, the celecoxib target the cancer literature means.",
        "* `p21` resolves to **H3P16, a pseudogene**, not CDKN1A.",
        "* `PC` splits across **Prostatic Neoplasms** and **Pyruvate Carboxylase",
        "  Deficiency**, with Pancreatic Neoplasms close behind -- and prostate versus",
        "  pancreatic cancer is about as consequential a confusion as this graph can",
        "  make.",
        "* `NO` resolves to **Nobelium**, the synthetic actinide, rather than nitric",
        "  oxide. That one is worth keeping in mind whenever a majority vote looks",
        "  authoritative.", "",
        "## What changed as a result", "",
        "Sense-colliding forms are written to the blocklist in this script's JSON",
        "output. `atlas_graph.resolve` consults it and returns `None` rather than a",
        "plausible-looking wrong entity, so an analysis built on an ambiguous form",
        "fails loudly instead of quietly reporting the wrong biology.", "",
        "Where the cancer-domain sense is genuinely not in doubt, `DOMAIN_SENSE`",
        "records it together with the reason it is safe, and a caller opts in",
        "explicitly. Nothing applies those silently.", "",
        "`FSP1` is deliberately given **no** domain default. It is not a symbol with",
        "one right answer that PubTator got wrong -- it is genuinely",
        "context-dependent, and the plurality of its cancer mentions are correctly",
        "S100A4. Resolving it needs the per-paper evidence in",
        "`scripts/atlas_disambiguate.py`.", "",
        "## Limits", "",
        f"* Only the top {top} forms per type are resolved, so the class split is",
        "  measured on those and extrapolated to none of the others. The tail is not",
        "  characterised.",
        "* The split uses the top two identifiers only. A form contested three ways",
        "  (`FSP1` and `PC` both are) is classified on its top two.",
        "* **The MeSH class split is measured on far fewer forms than it scanned.**",
        "  Supplementary Concept Records (the `C`-prefixed identifiers) carry no tree",
        "  number, so any pair involving one cannot be tested for nesting and is",
        "  reported unresolved rather than guessed at: 285 of the chemical forms and",
        "  136 of the disease forms land there. The chemical benign-vs-sense split",
        "  therefore rests on roughly a quarter of the forms examined, and should be",
        "  read as indicative. The **contested rate** in the table above does not have",
        "  this problem -- it is computed over every surface form, and it is the",
        "  genes-are-the-outlier comparison that carries the argument.",
        "* Hierarchical relatedness is detected by MeSH tree nesting, which misses",
        "  near-synonyms placed in different branches: *Small Cell Lung Carcinoma*",
        "  (C04.588...) and *Carcinoma, Small Cell* (C04.557...) are not nested, so",
        "  `SCLC` is counted as a sense collision when it is closer to a duplicate.",
        "  That makes the sense-collision count for diseases an over-estimate.",
        "* Species ambiguity is called benign for a literature map. That is a",
        "  judgement about this use, not a general one: a study of mouse-versus-human",
        "  biology would need exactly the distinction being collapsed here.",
    ]
    return "\n".join(lines) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--top", type=int, default=400,
                    help="how many highest-volume contested forms to resolve per type")
    ap.add_argument("--min-share", type=float, default=MIN_RUNNER_UP_SHARE)
    ap.add_argument("--render-only", action="store_true",
                    help="rebuild the report from the committed JSON without "
                         "re-scanning (the scan takes ~40 min, almost all of it "
                         "sequential MeSH lookups)")
    args = ap.parse_args()

    if args.render_only:
        raw = json.loads(RAW.read_text())
        results = raw["by_type"]
        OUT.write_text(render(results, raw.get("resolved_top_per_type", args.top)))
        print(f"re-rendered {OUT} from {RAW}")
        return 0

    ents = atlas_root() / "entities"
    results, blocklist = {}, []
    for kind in ("gene", "chemical", "disease"):
        f = ents / f"{kind}.tsv.gz"
        if not f.exists():
            print(f"missing {f}; run scripts/atlas_relations.py first", file=sys.stderr)
            return 1
        print(f"scanning {kind} ...", flush=True)
        kept, contested = scan(f, MIN_MENTIONS, args.min_share)
        cm = sum(d["total"] for d in contested.values())
        am = sum(sum(c.values()) for c in kept.values())
        if kind == "gene":
            benign, sense, unres = classify(contested, args.top)
        else:
            benign, sense, unres = classify_mesh(contested, args.top)
        print(f"  forms {len(kept):,}, contested {len(contested):,}, "
              f"benign {len(benign)}, sense {len(sense)}, unresolved {len(unres)}",
              flush=True)
        results[kind] = {"forms": len(kept), "contested": len(contested),
                         "mention_share": 100 * cm / max(1, am),
                         "benign": len(benign), "sense": len(sense),
                         "unresolved": len(unres), "sense_rows": sense,
                         "benign_rows": benign, "unresolved_rows": unres}
        blocklist += [r["surface"] for r in sense]

    OUT.write_text(render(results, args.top))
    RAW.write_text(json.dumps({
        "min_mentions": MIN_MENTIONS,
        "min_runner_up_share": args.min_share,
        "resolved_top_per_type": args.top,
        "by_type": {k: {kk: vv for kk, vv in v.items()} for k, v in results.items()},
        # kept for backwards compatibility with the gene-only first version
        "species_ambiguity": results["gene"]["benign_rows"],
        "sense_collision": results["gene"]["sense_rows"],
        "unresolved": results["gene"]["unresolved_rows"],
        "blocklist": sorted(set(blocklist)),
        "domain_sense": {k: {"id": v[0], "symbol": v[1], "why": v[2]}
                         for k, v in DOMAIN_SENSE.items()},
    }, indent=2) + "\n")
    print(f"wrote {OUT}\nwrote {RAW}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
