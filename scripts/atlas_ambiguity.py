#!/usr/bin/env python3
"""Atlas: measure gene-symbol ambiguity and stop the graph guessing (#ATLAS-AMBIG).

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
Two populations that must not be conflated, because only one is damaging:

* **Species ambiguity** -- `GAPDH` maps to human 2597 and mouse 14433. Same
  gene, different organism. For a literature map this is benign and often
  desirable.
* **Sense collision** -- `ER` maps to EREG (epiregulin) and ESR1 (estrogen
  receptor). Different genes. Collapsing them corrupts the edge.

Measured over the census: of the 400 highest-volume contested surface forms,
85.4% are species ambiguity and 14.6% are sense collisions. Reporting the
undecomposed "28.3% of mentions sit on a contested form" as an error rate would
overstate the damage by roughly sevenfold, so this script never reports it
alone.

WHAT IT CHANGES
---------------
Sense-colliding symbols go into a committed blocklist. `atlas_graph.resolve`
consults it and returns `None` instead of a plausible-looking wrong gene, so an
analysis fails loudly rather than quietly attributing ferroptosis biology to a
spastic-paraplegia gene. Where the cancer-domain sense is genuinely
unambiguous, `DOMAIN_SENSE` records it with its justification; a caller opts in
explicitly. `FSP1` is deliberately NOT given a domain default -- it is
genuinely context-dependent (see `atlas_disambiguate.py`).

Requires network (NCBI E-utilities) to name the identifiers. Derived artifacts
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
    "er": ("2099", "ESR1", "In the oncology literature 'ER' is the estrogen "
                           "receptor (ER-positive/ER-negative disease). The "
                           "majority vote returns EREG (epiregulin)."),
    "psa": ("354", "KLK3", "'PSA' in cancer is prostate-specific antigen, "
                           "KLK3. The majority vote returns NPEPPS."),
    "cox-2": ("5743", "PTGS2", "'COX-2' in cancer pharmacology is "
                               "prostaglandin-endoperoxide synthase 2, the "
                               "celecoxib target, not mitochondrial cytochrome "
                               "c oxidase subunit II."),
    "p62": ("8878", "SQSTM1", "'p62' in autophagy and redox signalling is "
                              "sequestosome-1. The majority vote returns "
                              "NUP62, a nuclear pore protein."),
    "p21": ("1026", "CDKN1A", "'p21' as a cell-cycle inhibitor is CDKN1A. The "
                              "majority vote returns H3P16, a pseudogene."),
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


def scan(path: Path, min_mentions: int, min_share: float) -> tuple:
    """Surface form -> identifier distribution, and the contested subset."""
    alias = collections.defaultdict(collections.Counter)
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            gid, mention = parts[2], parts[3]
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


def render(kept, contested, species, sense, unresolved, top) -> str:
    decided = len(species) + len(sense)
    multi = sum(1 for c in kept.values() if len(c) > 1)
    contested_mentions = sum(d["total"] for d in contested.values())
    all_mentions = sum(sum(c.values()) for c in kept.values())

    lines = [
        "# Gene-symbol ambiguity in the census graph (#ATLAS-AMBIG)", "",
        "Generated by `scripts/atlas_ambiguity.py`. The entity audit found one",
        "bad symbol (`FSP1` -> atlastin). This measures how many others there are.", "",
        "## Why a majority vote hides this", "",
        "`atlas_graph.build_index` resolves a surface form with",
        "`alias[s].most_common(1)[0][0]`. When a symbol is genuinely ambiguous the",
        "vote still returns exactly one identifier, with no signal that anything",
        "was discarded, so a wrong answer is indistinguishable from a right one.", "",
        "## Scale", "",
        f"| | count |", "|---|---|",
        f"| surface forms with >= {MIN_MENTIONS} mentions | {len(kept):,} |",
        f"| ... mapping to more than one identifier | {multi:,} |",
        f"| ... genuinely contested (runner-up >= {int(MIN_RUNNER_UP_SHARE*100)}%) | {len(contested):,} |",
        f"| mentions sitting on a contested form | {contested_mentions:,} of {all_mentions:,} ({100*contested_mentions/max(1,all_mentions):.1f}%) |",
        "", "## The decomposition that matters", "",
        "That last row is **not** an error rate, and reporting it as one would",
        "overstate the damage by roughly sevenfold. Contested forms split into two",
        "populations with opposite consequences:", "",
        f"Of the {top} highest-volume contested forms, with the top two identifiers",
        "resolved against NCBI's own record:", "",
        "| class | n | share | consequence |", "|---|---|---|---|",
        f"| species ambiguity (same symbol, different organism) | {len(species)} | "
        f"{100*len(species)/max(1,decided):.1f}% | benign -- same gene, and a literature map usually wants them merged |",
        f"| sense collision (genuinely different genes) | {len(sense)} | "
        f"{100*len(sense)/max(1,decided):.1f}% | damaging -- the edge is attributed to the wrong biology |",
        f"| unresolved | {len(unresolved)} | | not counted in either share |",
        "", "## The damaging class, by mention volume", "",
        "These are the ones an analysis is most likely to query, so a silent",
        "substitution here is most likely to reach a conclusion.", "",
        "| surface form | mentions | majority vote returns | runner-up |",
        "|---|---|---|---|",
    ]
    for r in sorted(sense, key=lambda r: -r["total"])[:25]:
        t, u = r["top"], r["runner_up"]
        lines.append(f"| `{r['surface']}` | {r['total']:,} | {t['symbol']} "
                     f"({t['organism']}) | {u['symbol']} ({u['organism']}) |")

    lines += [
        "", "Three are worth stating plainly, because they touch the most-queried",
        "concepts in the cancer literature:", "",
        "* `ER` resolves to **EREG (epiregulin)**, not ESR1. In oncology `ER` is the",
        "  estrogen receptor -- ER-positive and ER-negative disease is one of the",
        "  field's primary axes.",
        "* `COX-2` resolves to **mitochondrial cytochrome c oxidase subunit II**, not",
        "  PTGS2, the celecoxib target the cancer literature means.",
        "* `p21` resolves to **H3P16, a pseudogene**, not CDKN1A.", "",
        "## What changed as a result", "",
        "Sense-colliding forms are written to the blocklist in this script's JSON",
        "output. `atlas_graph.resolve` consults it and returns `None` rather than a",
        "plausible-looking wrong gene, so an analysis built on an ambiguous symbol",
        "fails loudly instead of quietly reporting the wrong biology.", "",
        "Where the cancer-domain sense is genuinely not in doubt, `DOMAIN_SENSE`",
        "records it together with the reason it is safe, and a caller opts in",
        "explicitly. Nothing applies those silently.", "",
        "`FSP1` is deliberately given **no** domain default. It is not a symbol with",
        "one right answer that PubTator got wrong -- it is genuinely",
        "context-dependent, and the majority of its cancer mentions are correctly",
        "S100A4. Resolving it needs the per-paper evidence in",
        "`scripts/atlas_disambiguate.py`.", "",
        "## Limits", "",
        "* Only the top-volume forms are resolved against NCBI, so the species/sense",
        f"  split is measured on {top} forms and extrapolated to none of the others.",
        "  The tail is not characterised.",
        "* The split uses the top two identifiers only. A form contested three ways",
        "  (`FSP1` is one) is classified on its top two.",
        "* Species ambiguity is called benign for a literature map. That is a",
        "  judgement about this use, not a general one: a study of mouse-versus-human",
        "  biology would need exactly the distinction being collapsed here.",
    ]
    return "\n".join(lines) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--top", type=int, default=400,
                    help="how many highest-volume contested forms to resolve")
    ap.add_argument("--min-share", type=float, default=MIN_RUNNER_UP_SHARE)
    args = ap.parse_args()

    genes = atlas_root() / "entities" / "gene.tsv.gz"
    if not genes.exists():
        print(f"missing {genes}; run scripts/atlas_relations.py first", file=sys.stderr)
        return 1

    print(f"scanning {genes} ...", flush=True)
    kept, contested = scan(genes, MIN_MENTIONS, args.min_share)
    print(f"  {len(kept):,} surface forms, {len(contested):,} contested", flush=True)

    species, sense, unresolved = classify(contested, args.top)
    print(f"  species {len(species)}, sense {len(sense)}, unresolved {len(unresolved)}",
          flush=True)

    OUT.write_text(render(kept, contested, species, sense, unresolved, args.top))
    RAW.write_text(json.dumps({
        "min_mentions": MIN_MENTIONS,
        "min_runner_up_share": args.min_share,
        "surface_forms": len(kept),
        "contested": len(contested),
        "resolved_top": args.top,
        "species_ambiguity": species,
        "sense_collision": sense,
        "unresolved": unresolved,
        # the blocklist atlas_graph.resolve consults
        "blocklist": sorted(r["surface"] for r in sense),
        "domain_sense": {k: {"id": v[0], "symbol": v[1], "why": v[2]}
                         for k, v in DOMAIN_SENSE.items()},
    }, indent=2) + "\n")
    print(f"wrote {OUT}\nwrote {RAW}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
