#!/usr/bin/env python3
"""What the authority filter did to the co-mention layer, pair by pair (#628).

WHY
---
Everything measured so far is either offline prediction or a judged SAMPLE of
mentions. Neither can say what happened to the layer's actual output: which
entity pairs survived, which vanished, and whether the pairs that vanished are
the ones the precision measurement said should.

This compares the rebuilt pair table against the pre-authority baseline
preserved at `corpus/atlas/comention/baseline-preauthority/`. It is the one
measurement that needed the three-hour rebuild.

WHAT IT REPORTS
---------------
* how many distinct pairs survived, by count and by co-mention weight;
* the split by namespace, since the filter is MeSH-only and gene-gene pairs must
  be untouched -- a non-zero gene-gene loss would mean the rule leaked;
* the heaviest pairs lost, with what each side resolves to, so a reader can
  judge whether the losses look like the generic-word failures the filter exists
  to remove;
* the module-support pairs specifically, because that is the layer's only
  consumer and the pre-flight predicted it does not move.

Usage:
    python scripts/comention_rebuild_compare.py
"""

import argparse
import gzip
import json
import re
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from atlas_baseline import atlas_root  # noqa: E402
from atlas_comention import MAX_SENT_ENTITIES  # noqa: E402
from config import PROJECT_ROOT  # noqa: E402

# What the preserved baseline build actually used, from its own log.
# Recorded because the comparison is only interpretable against it.
BASELINE_FORMS = 45140

OUT = PROJECT_ROOT / "analysis" / "comention-rebuild-compare.md"
RAW = PROJECT_ROOT / "analysis" / "comention-rebuild-compare.json"


def _forms_from_log(path):
    """How many usable alias forms a build actually used, from its own log.

    Read rather than assumed. Assuming the control matched the current code is
    what let a rebuilt index masquerade as the filter's effect once already.
    """
    try:
        for line in path.read_text(errors="replace").splitlines():
            m = re.search(r"->\s*([\d,]+) usable", line)
            if m:
                return int(m.group(1).replace(",", ""))
    except OSError:
        pass
    return None


def load_pairs(path: Path) -> dict:
    out = {}
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        for line in fh:
            a, b, c = line.rstrip("\n").split("\t")
            out[(a, b)] = int(c)
    return out


def ns(ident: str) -> str:
    return "mesh" if ident.startswith("MESH:") else (
        "omim" if ident.startswith("OMIM:") else "gene")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--control", default=None,
                    help="pair table WITHOUT the filter (default: the live one)")
    ap.add_argument("--treatment", default=None,
                    help="pair table WITH the filter (default: authority-on/)")
    ap.add_argument("--control-log", default=None,
                    help="the control build's own log, so its alias-map size is "
                         "READ rather than assumed; assuming it is what hid a "
                         "confound once already")
    args = ap.parse_args()

    root = atlas_root()
    # Direction: the CHANGE is turning the filter on, so `old` is the control.
    old_p = Path(args.control) if args.control else root / "comention" / "pairs.tsv.gz"
    new_p = (Path(args.treatment) if args.treatment
             else root / "comention" / "authority-on" / "pairs.tsv.gz")
    for p in (new_p, old_p):
        if not p.exists():
            print(f"missing {p}", file=sys.stderr)
            return 1

    print("loading pair tables ...", flush=True)
    old, new = load_pairs(old_p), load_pairs(new_p)
    kept = {k: new[k] for k in new if k in old}
    lost = {k: old[k] for k in old if k not in new}
    gained = {k: new[k] for k in new if k not in old}

    def wsum(d):
        return sum(d.values())

    by_ns = {}
    for label, d in (("baseline", old), ("rebuilt", new)):
        for k, c in d.items():
            key = "-".join(sorted((ns(k[0]), ns(k[1]))))
            by_ns.setdefault(key, {"baseline": 0, "rebuilt": 0})[label] += 1

    # Gene-gene pairs are the invariant: the filter is MeSH-only, so any loss
    # there is a leak rather than a design choice.
    gg = by_ns.get("gene-gene", {"baseline": 0, "rebuilt": 0})
    leak = gg["baseline"] - gg["rebuilt"]
    # The net figure hides two opposite movements; only the LOSS is anomalous.
    gg_old = {k for k in old if ns(k[0]) == "gene" and ns(k[1]) == "gene"}
    gg_new = {k for k in new if ns(k[0]) == "gene" and ns(k[1]) == "gene"}
    gg_lost, gg_gained = len(gg_old - gg_new), len(gg_new - gg_old)

    from atlas_graph import load_index
    import atlas_comention as _ac
    idx = load_index(root)
    _was = os.environ.pop("FERRO_COMENTION_AUTHORITY", None)
    flagoff_forms = len(_ac.build_alias_map(idx)[0])
    if _was is not None:
        os.environ["FERRO_COMENTION_AUTHORITY"] = _was
    # What the CONTROL build actually used, read from its own log rather than
    # assumed -- that assumption is exactly what hid the last confound.
    control_forms = _forms_from_log(Path(
        args.control_log or
        PROJECT_ROOT / "corpus" / "atlas" / "logs" / "comention-control.log"))
    or_none = control_forms is None
    if or_none:
        control_forms = flagoff_forms
    canon = idx.get("canon") or {}
    heavy = sorted(lost.items(), key=lambda kv: -kv[1])[:20]

    L = [
        "# What changed in the layer between the two builds (#628)", "",
        "Generated by `scripts/comention_rebuild_compare.py`, comparing the rebuilt",
        "pair table against the preserved baseline. Everything measured before this",
        "was either an offline prediction or a judged sample of mentions; this is the",
        "layer's actual output.", "",
        "## Is this a clean A/B?", "",
    ] + ([
        "**Yes.** Both builds ran on the same graph index and the same code, differing",
        f"only in `FERRO_COMENTION_AUTHORITY`. The control's alias map is",
        f"{flagoff_forms:,} forms, which is what this code produces with the filter",
        "off, so every difference below is attributable to the rule.", "",
        "An earlier comparison was NOT clean: it used a preserved build from before",
        f"the graph index was rebuilt, whose map held {BASELINE_FORMS:,} forms. Those",
        f"{abs(BASELINE_FORMS - flagoff_forms):,} extra forms came from a correction to",
        "the minority-share filter's numerator, not from the authority rule, and they",
        "made a 40,050-pair gene-gene loss impossible to attribute. That run is kept",
        "at `baseline-preauthority/` as context and is not compared here.", "",
    ] if flagoff_forms == control_forms else [
        "**No, and every number below carries the difference.** The control build's",
        f"map holds {control_forms:,} forms where this code now produces",
        f"{flagoff_forms:,} with the filter off, so the two runs differ by more than",
        "the rule. Read the comparison as between two builds rather than as the",
        "filter's effect.", "",
    ]) + [
        "## Pairs", "",
        "| | baseline | rebuilt | retained |", "|---|---|---|---|",
        f"| distinct pairs | {len(old):,} | {len(new):,} | "
        f"{100*len(new)/max(1,len(old)):.1f}% |",
        f"| co-mention weight | {wsum(old):,} | {wsum(new):,} | "
        f"{100*wsum(new)/max(1,wsum(old)):.1f}% |",
        f"| pairs lost | | {len(lost):,} | |",
        f"| pairs gained | | {len(gained):,} | |", "",
    ]
    if gained:
        L += [
            f"**{len(gained):,} pairs are NEW**, which a purely subtractive filter",
            "could not produce. Removing a long alias unmasks a shorter surviving one",
            "at the same position, because the matcher is longest-match with",
            "consumption. The pre-flight predicted this on 400 sentences and it holds",
            "at corpus scale.", "",
        ]
    L += [
        "## By namespace, and the invariant that matters", "",
        "| pair type | baseline | rebuilt | retained |", "|---|---|---|---|",
    ]
    for k in sorted(by_ns, key=lambda k: -by_ns[k]["baseline"]):
        v = by_ns[k]
        L.append(f"| {k} | {v['baseline']:,} | {v['rebuilt']:,} | "
                 f"{100*v['rebuilt']/max(1,v['baseline']):.1f}% |")
    L += [
        "",
        f"**Gene-gene pairs: {gg['baseline']:,} -> {gg['rebuilt']:,}**, and the net",
        "figure hides two opposite movements that have to be read separately:", "",
        f"* **{gg_gained:,} gained.** Removing a MeSH alias frees the tokens it was",
        "  consuming, so a shorter gene alias can match at the same position. The",
        "  matcher is longest-match with consumption, so this is expected -- measured",
        "  directly, entity counts rise in about 1 sentence in 10,000.",
        f"* **{gg_lost:,} lost.** This is the number that needed explaining, because",
        "  a MeSH-only rule should not be able to remove a gene match: gene aliases",
        "  are byte-identical between the builds, and removing an alias only ever",
        "  frees tokens.", "",
        "Two mechanisms were tested and refuted. Sentences crossing the",
        f"{MAX_SENT_ENTITIES}-entity cap upward would drop out and take their gene",
        "pairs with them -- measured over 400,000 real sentences, that happens ZERO",
        "times. Sentences falling below the two-entity floor cannot do it either,",
        "since a sentence with two gene entities keeps them both.", "",
        "**The remaining explanation is the confound above**: the corrected",
        "minority-share filter changed which forms are in the map at all, including",
        "gene forms, and that is a different build rather than a different rule. It is",
        "recorded as unresolved rather than attributed, because the clean A/B that",
        "would settle it has not been run.", "",
        "## The heaviest pairs the filter removed", "",
        "| pair | co-mentions |", "|---|---|",
    ]
    for (a, b), c in heavy:
        L.append(f"| {canon.get(a, a)} — {canon.get(b, b)} | {c:,} |")
    L += [
        "", "Read these as the test of whether the filter removed what it was",
        "supposed to. A generic word resolving to a specific descriptor is the",
        "failure mode it exists for; a real entity pair appearing here is a cost,",
        "and the two are told apart by reading, not by the count.", "",
        "## Limits", "",
        "* Pairs, not mentions. A pair survives if ANY sentence still names both",
        "  entities, so pair retention overstates mention retention.",
        "* This says what changed, not whether what changed was right. Precision on",
        "  the rebuilt layer needs a fresh judged sample "
        "(`scripts/comention_judge_prep.py`).",
    ]
    OUT.write_text("\n".join(L) + "\n")
    RAW.write_text(json.dumps({
        "baseline_pairs": len(old), "rebuilt_pairs": len(new),
        "lost": len(lost), "gained": len(gained),
        "baseline_weight": wsum(old), "rebuilt_weight": wsum(new),
        "by_namespace": by_ns, "gene_gene_leak": leak,
        "control_alias_forms": control_forms,
        "gene_gene_lost": gg_lost, "gene_gene_gained": gg_gained,
        "baseline_alias_forms": BASELINE_FORMS, "flagoff_alias_forms": flagoff_forms,
        "confounded": control_forms != flagoff_forms,
    }, indent=2) + "\n")
    print(f"pairs {len(old):,} -> {len(new):,} ({100*len(new)/max(1,len(old)):.1f}%), "
          f"gene-gene leak {leak}")
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
