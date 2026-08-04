#!/usr/bin/env python3
"""Is each alias form a NAME of the entity it resolves to? (#628)

WHY
---
The co-mention layer matches full text against an alias map built from
PubTator's annotations. `analysis/comention-regression.md` measured its
precision by hand at roughly 42% and found the failure mode to be generic
English words -- `treatment`, `effects`, `as`, `left` -- resolving to specific
descriptors. Hand judging bounds that on samples of tens.

`analysis/comention/authority-labels.tsv.gz` makes it answerable for the whole
map. For every identifier the map resolves to, it carries the names NLM and NCBI
give it: a MeSH label, or a gene's official symbol plus its description and
listed aliases. A form that matches none of them is not a name of that entity by
anyone's account.

WHAT IT FINDS
-------------
Only about 54% of the layer's census mentions sit on a form that is a name of
what it resolves to. The form-level figure is higher, near 60%, and the gap runs
the WRONG way for the layer: non-name forms are used MORE than name forms,
because generic English words are common words.

The discriminator this measurement suggests works on MeSH and not on genes, and
the reason is that genes did not need it. A gene symbol is already a specific
string, so the gene subset of the judged samples is high precision before any
filtering; MeSH descriptors are where generic words land.

Usage:
    python scripts/comention_name_check.py
"""

import csv
import json
import math
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from atlas_baseline import atlas_root  # noqa: E402
from atlas_comention import build_alias_map  # noqa: E402
from atlas_graph import load_index  # noqa: E402
from build_label_source import load_table  # noqa: E402
from config import PROJECT_ROOT  # noqa: E402

OUT = PROJECT_ROOT / "analysis" / "comention-name-check.md"
RAW = PROJECT_ROOT / "analysis" / "comention-name-check.json"
JUDGED = PROJECT_ROOT / "analysis" / "comention"
SAMPLES = ["abstract-visible-judgements.csv",
           "abstract-visible-heldout-judgements.csv",
           "body-only-judgements.csv", "corroborated-judgements.csv"]


def bag(s):
    return frozenset(w for w in re.split(r"[^a-z0-9]+", (s or "").lower()) if w)


def wilson(k, n, z=1.96):
    if not n:
        return (0.0, 0.0)
    p = k / n
    den = 1 + z * z / n
    c = p + z * z / (2 * n)
    m = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return (max(0.0, (c - m) / den), min(1.0, (c + m) / den))


# Normalisation for the word-bag comparison. Strict equality is what the first
# measurement used; it fails on plurals and on the cancer/tumour/neoplasm
# synonym set, which is most of MeSH's oncology vocabulary.
_PLURAL = [("ies", "y"), ("ses", "se")]
_SYNONYM = {"cancer": "neoplasm", "tumor": "neoplasm", "tumour": "neoplasm",
            "carcinoma": "neoplasm", "neoplasia": "neoplasm",
            "malignancy": "neoplasm"}


# Dropped before comparing. MeSH stores inverted renderings -- `Cancer of
# Breast` rather than `Breast Cancer` -- so a word bag matches the literature's
# form only once these are gone. Worth 3.3 points of cancer vocabulary, and
# worth nothing at all before entry terms were in the table, which is why it was
# measured twice.
_STOPWORD = {"of", "the", "and", "in", "with", "for", "to", "or", "by", "an"}


def norm_bag(s):
    """Word bag with plurals stemmed, 1-char tokens dropped, and the
    cancer/tumour family collapsed onto `neoplasm`."""
    out = set()
    for w in re.split(r"[^a-z0-9]+", (s or "").lower()):
        if len(w) < 2:
            continue
        w = _SYNONYM.get(w, w)
        for a, b in _PLURAL:
            if w.endswith(a):
                w = w[:-len(a)] + b
                break
        else:
            if w.endswith("s") and not w.endswith("ss"):
                w = w[:-1]
        w = _SYNONYM.get(w, w)
        if w not in _STOPWORD:
            out.add(w)
    return frozenset(out)


def c04_cost(root, alias, labels, support, fn):
    """What the rule would do to the cancer vocabulary the layer exists for.

    The judged samples measure precision. They cannot say whether the rule
    deletes the concepts the layer is FOR, and MeSH tree C04 is that list.
    """
    f = root / "mesh" / "c04-descriptors.tsv"
    if not f.exists():
        return None
    c04 = {"MESH:" + line.split("\t")[0]
           for line in f.read_text().splitlines()
           if line and not line.startswith("#")}
    by_ident = {}
    for form, ident in alias.items():
        if ident in c04:
            by_ident.setdefault(ident, []).append(form)
    kept = tot = dead = 0
    for ident, forms in by_ident.items():
        names = labels.get(ident, [])
        surv = [x for x in forms if any(fn(n) == fn(x) for n in names)]
        tot += sum(support.get(x, 0) for x in forms)
        kept += sum(support.get(x, 0) for x in surv)
        if not surv:
            dead += 1
    return {"descriptors_reachable": len(by_ident),
            "mass_retained": kept / max(1, tot), "descriptors_killed": dead}


def unreachable_analysis(root, alias, labels) -> dict:
    """Why can the layer not reach most of the identifiers PubTator annotates?

    `Apoptosis` is unreachable because `apoptosis` resolves to a cortical
    malformation descriptor instead. That looked like the tip of something, so
    this counts the class -- and finds it nearly empty, which is the result.
    """
    import collections
    import gzip as _gzip

    reachable = set(alias.values())
    mentions = collections.Counter()
    for kind in ("gene", "chemical", "disease"):
        f = root / "entities" / f"{kind}.tsv.gz"
        if not f.exists():
            continue
        with _gzip.open(f, "rt", errors="replace") as fh:
            for line in fh:
                p_ = line.rstrip("\n").split("\t")
                if len(p_) >= 3 and p_[2] and p_[2] != "-":
                    mentions[p_[2]] += 1
    miss = sorted(((i, c) for i, c in mentions.items() if i not in reachable),
                  key=lambda kv: -kv[1])[:3000]
    owner = {bag(f): i for f, i in alias.items()}

    unused = same = different = unlabelled = 0
    examples = []
    for ident, c in miss:
        names = labels.get(ident)
        if not names:
            unlabelled += 1
            continue
        w = owner.get(bag(names[0]))
        if w is None:
            unused += 1
        elif bag(labels.get(w, [""])[0]) == bag(names[0]):
            same += 1
        else:
            different += 1
            if len(examples) < 8:
                examples.append((names[0], labels.get(w, ["?"])[0], c))
    n = len(miss) - unlabelled
    return {"annotated": len(mentions), "reachable": len(reachable),
            "unreachable": len(mentions) - len(reachable),
            "unreachable_mentions": sum(c for i, c in mentions.items()
                                        if i not in reachable),
            "total_mentions": sum(mentions.values()),
            "examined": n, "name_unused": unused, "name_to_same": same,
            "name_to_different": different,
            "examples": sorted(examples, key=lambda e: -e[2])}


def judged_rows():
    rows = []
    for f in SAMPLES:
        p = JUDGED / f
        if not p.exists():
            continue
        for r in csv.DictReader(p.open()):
            r["_canon"] = r.get("entity") or r.get("surface_form") or ""
            r["_v"] = r.get("verdict_v2") or r.get("verdict") or ""
            if r["_v"] in ("TP", "FP"):
                rows.append(r)
    return rows


def main() -> int:
    labels = load_table()
    if not labels:
        print("no authority table; run scripts/build_label_source.py", file=sys.stderr)
        return 1
    idx = load_index(atlas_root())
    alias, _ = build_alias_map(idx)
    support = idx.get("alias_support") or {}

    f_name = f_other = m_name = m_other = 0
    for form, ident in alias.items():
        names = labels.get(ident)
        if not names:
            continue
        w = support.get(form, 0)
        if any(bag(n) == bag(form) for n in names):
            f_name += 1
            m_name += w
        else:
            f_other += 1
            m_other += w
    ft, mt = f_name + f_other, m_name + m_other

    rows = judged_rows()
    u = unreachable_analysis(atlas_root(), alias, labels)

    def evaluate_form(cmp_fn):
        """What a filter sees: the form alone, before any sentence exists."""
        sub = [r for r in rows if r["identifier"].startswith("MESH:")
               and labels.get(r["identifier"])]
        keep = []
        for r in sub:
            span = (r.get("matched_span") or "").split("|")[0]
            if span and any(cmp_fn(n) == cmp_fn(span) for n in labels[r["identifier"]]):
                keep.append(r)
        tp = sum(1 for r in sub if r["_v"] == "TP")
        ktp = sum(1 for r in keep if r["_v"] == "TP")
        fp = len(sub) - tp
        return {"n": len(sub), "kept": len(keep), "kept_tp": ktp,
                "kept_precision": ktp / len(keep) if keep else 0.0,
                "fp_removed": (fp - (len(keep) - ktp)) / fp if fp else 0.0,
                "tp_removed": (tp - ktp) / tp if tp else 0.0}

    impl, norm = evaluate_form(bag), evaluate_form(norm_bag)
    support = idx.get("alias_support") or {}
    c04_strict = c04_cost(atlas_root(), alias, labels, support, bag)
    c04_norm = c04_cost(atlas_root(), alias, labels, support, norm_bag)

    def evaluate(sel):
        sub = [r for r in rows if sel(r) and labels.get(r["identifier"])]
        keep = []
        for r in sub:
            names = labels[r["identifier"]]
            cands = {bag(r["_canon"])}
            span = (r.get("matched_span") or "").split("|")[0]
            if span:
                cands.add(bag(span))
            if any(bag(n) in cands for n in names):
                keep.append(r)
        tp = sum(1 for r in sub if r["_v"] == "TP")
        ktp = sum(1 for r in keep if r["_v"] == "TP")
        fp = len(sub) - tp
        return {"n": len(sub), "base": tp / len(sub) if sub else 0.0,
                "kept": len(keep), "kept_tp": ktp,
                "kept_precision": ktp / len(keep) if keep else 0.0,
                "kept_ci": wilson(ktp, len(keep)),
                "fp_removed": (fp - (len(keep) - ktp)) / fp if fp else 0.0,
                "tp_removed": (tp - ktp) / tp if tp else 0.0}

    allr = evaluate(lambda r: True)
    gene = evaluate(lambda r: not r["identifier"].startswith(("MESH:", "OMIM:")))
    mesh = evaluate(lambda r: r["identifier"].startswith("MESH:"))

    L = [
        "# Is each alias form a name of what it resolves to? (#628)", "",
        "Generated by `scripts/comention_name_check.py`.", "",
        "Hand judging bounded this layer's precision on samples of tens. The",
        "committed authority table (`analysis/comention/authority-labels.tsv.gz`,",
        f"{len(labels):,} identifiers) makes the same question answerable for the",
        "whole alias map: for each form, is it a name of the entity it resolves to",
        "according to NLM or NCBI?", "",
        "## The whole map", "",
        "| | forms | census mentions |", "|---|---|---|",
        f"| a name of the entity it resolves to | {f_name:,} ({100*f_name/ft:.1f}%) | "
        f"{m_name:,} ({100*m_name/mt:.1f}%) |",
        f"| not a name of it | {f_other:,} ({100*f_other/ft:.1f}%) | "
        f"{m_other:,} ({100*m_other/mt:.1f}%) |", "",
        f"**About {100*m_name/mt:.0f}% of the layer's mentions sit on a form that is",
        "a name of what it resolves to.** The mention figure is LOWER than the form",
        f"figure ({100*m_name/mt:.1f}% against {100*f_name/ft:.1f}%), which is the",
        "unfavourable direction: non-name forms are used MORE than name forms,",
        "because generic English words are common words.", "",
        "This is a property of the alias map, not a precision estimate. A form can",
        "be a real name and still be the wrong entity in context, and a form that is",
        "not a listed name can still be correct. It sits close to the hand-judged",
        "precision, and the two are independent measurements of related things,",
        "which is worth noting and not worth over-reading.", "",
        "## The one identifier that shows what this means", "",
        "The alias map has **no route to MeSH `Apoptosis` (D017209) at all**. The",
        "form `apoptosis` resolves only to D065703, *Malformations of Cortical",
        "Development, Group I*. So every apoptosis co-mention in the layer is filed",
        "under a cortical malformation descriptor, and no query for the correct",
        "identifier can reach any of them. That is not ambiguity between two senses;",
        "it is the right answer being unreachable.", "",
        "## Where the discriminator works, and where it does not", "",
        f"Applied to the {len(rows)} hand-judged mentions, split by namespace:", "",
        "| | judged | base precision | keeps | kept precision | FPs removed | TPs removed |",
        "|---|---|---|---|---|---|---|",
    ]
    for name, d in (("all", allr), ("MeSH identifiers", mesh), ("gene identifiers", gene)):
        L.append(
            f"| {name} | {d['n']} | {100*d['base']:.1f}% | {d['kept']} | "
            f"{100*d['kept_precision']:.1f}% [{100*d['kept_ci'][0]:.0f}, "
            f"{100*d['kept_ci'][1]:.0f}] | {100*d['fp_removed']:.0f}% | "
            f"{100*d['tp_removed']:.0f}% |")

    L += [
        "", "## Why the layer cannot reach most identifiers, and why that is mostly fine",
        "",
        f"PubTator annotates **{u['annotated']:,}** distinct identifiers in the",
        f"census. The alias map resolves to **{u['reachable']:,}** of them, so",
        f"{100*u['unreachable']/u['annotated']:.0f}% are unreachable, carrying",
        f"{100*u['unreachable_mentions']/u['total_mentions']:.0f}% of all annotated",
        "mentions.", "",
        "The `Apoptosis` case above suggested that might be silent corruption: an",
        "entity locked out because a competitor won its name. Counting it over the",
        f"{u['examined']:,} most-mentioned unreachable identifiers says otherwise.", "",
        "| why it is unreachable | count | share |", "|---|---|---|",
        f"| its own name is never used as a surface form | {u['name_unused']:,} | "
        f"{100*u['name_unused']/u['examined']:.1f}% |",
        f"| its name resolves to an entity carrying the SAME name | "
        f"{u['name_to_same']:,} | {100*u['name_to_same']/u['examined']:.1f}% |",
        f"| its name resolves to a DIFFERENT entity | {u['name_to_different']:,} | "
        f"{100*u['name_to_different']/u['examined']:.1f}% |", "",
        "The middle row is the ortholog and synonym class already measured in",
        "`analysis/atlas-ambiguity.md`: the human gene wins the symbol and the mouse",
        "one becomes unreachable, which for a cancer-literature question is a",
        "reasonable outcome rather than a fault.", "",
        "**And the last row is not what it looks like either.** Its largest cases are",
        "the same biological entity holding both a MeSH descriptor and an NCBI gene",
        "id, with the gene winning the form:", "",
        "| loses its name | to | mentions |", "|---|---|---|",
    ] + [f"| {a} | {b} | {c:,} |" for a, b, c in u["examples"]] + [
        "", "`Prostate-Specific Antigen` and `KLK3` are not competitors. Neither are",
        "`Insulin` and `INS`. This is cross-namespace redundancy, and the gene is",
        "arguably the better target of the two.", "",
        "**So the class that `Apoptosis` belongs to is close to empty.** That is a",
        "negative result about a hypothesis raised in the previous section, and it",
        "matters because it bounds how far that finding should be read: apoptosis",
        "co-mentions really are filed under a cortical malformation descriptor, and",
        "that really is unrecoverable, but it is an outlier rather than a symptom.",
        "The reason most identifiers are unreachable is that nobody writes their",
        "name, or that a near-identical entity holds it.", "",
        "## Two costs this measurement did not carry, and both are decisive", "",
        "The rows above score a mention as a name match if EITHER the identifier's",
        "canonical form OR the span that fired matches an authority name. That is a",
        "fair diagnostic, and it is not what a filter would see. A rule inside",
        "`build_alias_map` decides per FORM, before any sentence exists, so only the",
        "form is available to it.", "",
        "| what is compared | keeps | kept precision | FPs removed | TPs removed |",
        "|---|---|---|---|---|",
        f"| canonical form OR span (reported above) | {mesh['kept']} | "
        f"{100*mesh['kept_precision']:.1f}% | {100*mesh['fp_removed']:.0f}% | "
        f"{100*mesh['tp_removed']:.0f}% |",
        f"| **form only, as a filter would** | {impl['kept']} | "
        f"{100*impl['kept_precision']:.1f}% | {100*impl['fp_removed']:.0f}% | "
        f"**{100*impl['tp_removed']:.0f}%** |",
        f"| form only, normalised (recommended) | {norm['kept']} | "
        f"{100*norm['kept_precision']:.1f}% | {100*norm['fp_removed']:.0f}% | "
        f"{100*norm['tp_removed']:.0f}% |", "",
        "At the insertion point the rule is CLEANER than reported -- it removes every",
        "span-bearing false positive -- and costs half again as many true positives,",
        f"{100*impl['tp_removed']:.0f}% rather than {100*mesh['tp_removed']:.0f}%. The",
        "favourable number is the one this document originally quoted.", "",
    ] + ([] if not (c04_strict and c04_norm) else [
        "### And it would delete most of the cancer vocabulary", "",
        "Precision on judged mentions cannot see this. MeSH tree C04 is the cancer",
        "definition the whole census is built on, so what the rule does to C04 is",
        "what it does to the layer's purpose:", "",
        "| rule | C04 census mass retained | C04 descriptors left with no route |",
        "|---|---|---|",
        f"| form only, strict | **{100*c04_strict['mass_retained']:.1f}%** | "
        f"{c04_strict['descriptors_killed']} of {c04_strict['descriptors_reachable']} |",
        f"| form only, normalised | {100*c04_norm['mass_retained']:.1f}% | "
        f"{c04_norm['descriptors_killed']} of {c04_norm['descriptors_reachable']} |", "",
        f"**The strict rule discards {100*(1-c04_strict['mass_retained']):.0f}% of the",
        "cancer vocabulary's mentions and leaves",
        f"{c04_strict['descriptors_killed']} cancer descriptors unreachable by any form",
        "at all.** `Neoplasms`, `Breast Neoplasms` and `Lung Neoplasms` are among them,",
        "because MeSH writes `Breast Neoplasms` where the literature writes",
        "`breast cancer`.", "",
        "Normalising plurals and collapsing the cancer/tumour/carcinoma family onto",
        "`neoplasm` recovers most of that at NO precision cost -- the same 100%",
        "false-positive removal, the same kept precision, and the true-positive cost",
        f"falls from {100*impl['tp_removed']:.0f}% to {100*norm['tp_removed']:.0f}%.",
        "**Any implementation should normalise; the strict form measured above should",
        "not be built.**", "",
        "Two changes got it there, and their order matters because each was worth",
        "little without the other:", "",
        "* **MeSH entry terms in the authority table.** A preferred label is not how",
        "  the literature writes a concept -- MeSH says `Breast Neoplasms`, papers say",
        "  `breast cancer`. Adding every term of every concept the descriptor carries",
        "  took C04 retention from 24.7% to 35.0% strict, 64.2% to 69.9% normalised.",
        "* **Dropping stopwords.** MeSH stores INVERTED renderings, `Cancer of Breast`",
        "  rather than `Breast Cancer`, so the bag carries an `of` the literature's",
        "  form does not. Measured BEFORE entry terms this was worth 1.0 point and",
        "  looked not worth having; measured after, it is worth 3.3 and takes",
        "  retention to 73.2%. A normalisation is only as good as the reference it",
        "  compares against.", "",
    ]) + [
        "**It works on MeSH and does nothing useful on genes**, and the reason is",
        "that genes did not need it. A gene symbol is already a specific string, so",
        f"the gene subset is {100*gene['base']:.0f}% precise BEFORE any filtering,",
        f"against {100*mesh['base']:.0f}% for MeSH. On genes the rule removes",
        f"{100*gene['fp_removed']:.0f}% of false positives while still costing",
        f"{100*gene['tp_removed']:.0f}% of true ones -- it has nothing to gain and",
        "something to lose. On MeSH it removes",
        f"{100*mesh['fp_removed']:.0f}% of false positives and lifts precision from",
        f"{100*mesh['base']:.0f}% to {100*mesh['kept_precision']:.0f}%.", "",
        "So the honest recommendation is narrower than the rule first looked: apply",
        "it to MeSH identifiers only. Getting gene labels did not make the rule work",
        "for genes -- it made it TESTABLE for genes, and the test says not to.", "",
        "The gene row rests on "
        f"{gene['n']} judged mentions, so 'does nothing useful' is weakly supported",
        "and 'is not needed' is the safer reading of it.", "",
        "## What the label table changed about the rule itself", "",
        "Earlier measurements compared a form against a single preferred label and",
        "reported a 60% true-positive cost. Comparing against every name an",
        "authority lists -- a gene's symbol, description and aliases -- drops that to",
        f"{100*allr['tp_removed']:.0f}%, because `xCT` really is a name of SLC7A11 and",
        "`PHGPx` one of GPX4. The rule did not change; the reference did.", "",
        "## Limits", "",
        "* Names, not senses. `FSP1` is a listed alias of ATL1, so this check calls",
        "  that resolution a name match. It is still the wrong gene in a ferroptosis",
        "  paper, which is what `atlas_disambiguate.py` exists for. This measures",
        "  whether a form is a name of an entity, never whether it is the right",
        "  entity here.",
        "* 244 alias forms resolve to identifiers absent from the table (76",
        "  withdrawn gene ids and their forms); they are excluded rather than",
        "  guessed.",
        "* Ortholog pairs count as name matches, since `Mmp2` and `MMP2` share a word",
        "  bag. That is the benign species class from `analysis/atlas-ambiguity.md`",
        "  and inflates the name-match rate slightly.",
        "* The judged mentions were drawn for a different purpose and the gene subset",
        "  is small; the namespace split is a signal, not a settled result.",
    ]

    OUT.write_text("\n".join(L) + "\n")
    RAW.write_text(json.dumps({
        "identifiers_with_labels": len(labels),
        "forms": {"name": f_name, "other": f_other, "share": f_name / ft},
        "mentions": {"name": m_name, "other": m_other, "share": m_name / mt},
        "discriminator": {"all": allr, "mesh": mesh, "gene": gene},
        "unreachable": {k: v for k, v in u.items() if k != "examples"},
        "at_insertion_point": {"strict": impl, "normalised": norm},
        "c04_cost": {"strict": c04_strict, "normalised": c04_norm},
    }, indent=2) + "\n")
    print(f"forms {100*f_name/ft:.1f}% names, mentions {100*m_name/mt:.1f}% names")
    print(f"discriminator: mesh removes {100*mesh['fp_removed']:.0f}% of FPs, "
          f"gene removes {100*gene['fp_removed']:.0f}%")
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
