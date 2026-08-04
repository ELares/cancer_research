#!/usr/bin/env python3
"""Prepare a reproducible co-mention judging sample (#628).

WHY
---
The 140 verdicts behind `analysis/comention-regression.md` were produced by a
sampling seed and a span-recovery step that existed only in one session's shell
history. The numbers are committed and the verdicts are committed, but the
SCAFFOLD they were judged on was not, so nobody could re-derive them or draw a
comparable new sample. That is a provenance hole in the measurement this project
leans on hardest for the co-mention layer.

It also blocks the one criterion #628 still needs. Every judged mention so far
was used to SELECT the authority-name discriminator, so none of them can
validate it; that requires a held-out sample drawn the same way and disjoint
from the judged ones. `--exclude` is how this script draws it.

WHAT IT DOES
------------
Reads the strata dumps written by `atlas_comention_audit.py --dump-strata` and
`--dump-abstract-visible`, samples deterministically from a NAMED seed, and
attaches the two things a verdict cannot be checked without:

  matched_span    the span that actually fired, recovered with the build's own
                  longest-match-with-consumption matcher rather than an n-gram
                  scan, so it cannot list forms that never matched
  authority_name  what the identifier denotes according to NLM, because the
                  index's `canon` is the most frequent surface form and is
                  therefore circular

The output carries an empty `verdict` column and is judged by hand. Nothing here
assigns a verdict.

LIMITS
------
NCBI Gene and OMIM identifiers have no MeSH descriptor, so `authority_name` is
blank for them and the judge has only the span and the sentence. That is a real
gap, not a formatting one -- roughly a quarter of the corroborated stratum.

Usage:
    python scripts/atlas_comention_audit.py --dump-strata /tmp/strata \\
        --dump-abstract-visible /tmp/strata/abstract-visible.jsonl
    python scripts/comention_judge_prep.py --stratum body-only --n 40 \\
        --seed 20260804 --out analysis/comention/body-only-judgements.csv
    # a held-out sample, disjoint from what has already been judged:
    python scripts/comention_judge_prep.py --stratum body-only --n 40 \\
        --seed 20260805 --exclude analysis/comention/body-only-judgements.csv \\
        --out /tmp/body-only-heldout.csv
"""

import argparse
import csv
import json
import random
import sys
import time
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from atlas_baseline import atlas_root  # noqa: E402
from atlas_comention import build_alias_map, matched_forms  # noqa: E402
from atlas_graph import load_index  # noqa: E402
from config import PROJECT_ROOT  # noqa: E402

COLUMNS = ["n", "pmid", "identifier", "entity", "matched_span",
           "authority_name", "verdict", "sentence"]
CACHE = PROJECT_ROOT / "analysis" / "comention" / ".mesh-labels.json"


def mesh_label(ident: str, cache: dict) -> str:
    """NLM's label for a MeSH descriptor, cached so a re-run is offline."""
    if ident in cache:
        return cache[ident]
    if not ident.startswith("MESH:"):
        cache[ident] = ""
        return ""
    try:
        url = f"https://id.nlm.nih.gov/mesh/{ident.split(':', 1)[1]}.json"
        with urllib.request.urlopen(url, timeout=30) as f:
            j = json.load(f)
        lbl = j.get("label")
        lbl = lbl.get("@value") if isinstance(lbl, dict) else lbl
    except Exception:
        lbl = ""
    time.sleep(0.34)
    cache[ident] = lbl or ""
    return cache[ident]


def already_judged(paths) -> set:
    """(pmid, identifier) pairs that have been judged, so a held-out sample
    can be drawn disjoint from them."""
    seen = set()
    for p in paths or []:
        for r in csv.DictReader(Path(p).open()):
            seen.add((r["pmid"], r["identifier"]))
    return seen


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stratum", required=True,
                    choices=["body-only", "corroborated", "abstract-visible"])
    ap.add_argument("--dumps", default="/tmp/strata",
                    help="directory written by --dump-strata")
    ap.add_argument("--n", type=int, default=40)
    ap.add_argument("--seed", type=int, required=True,
                    help="recorded in the output so the draw is reproducible")
    ap.add_argument("--exclude", nargs="*",
                    help="judged CSVs whose mentions must not be re-drawn")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    src = Path(args.dumps) / f"{args.stratum}.jsonl"
    if not src.exists():
        print(f"missing {src}; run scripts/atlas_comention_audit.py --dump-strata "
              f"{args.dumps} first", file=sys.stderr)
        return 1
    rows = [json.loads(l) for l in src.read_text().splitlines() if l.strip()]

    skip = already_judged(args.exclude)
    pool = [r for r in rows if (r["pmid"], r["identifier"]) not in skip]
    if skip:
        print(f"excluded {len(rows) - len(pool)} already-judged mentions")
    if len(pool) < args.n:
        print(f"pool has only {len(pool)} mentions, fewer than the {args.n} "
              f"requested", file=sys.stderr)
        return 1

    # Sort before shuffling. The dump is written in encounter order, which is an
    # artifact of shard traversal; sorting makes the draw depend on the seed and
    # the data, not on how the file happened to be produced.
    pool.sort(key=lambda r: (r["pmid"], r["identifier"], r["sentence"][:80]))
    random.Random(args.seed).shuffle(pool)
    sample = pool[:args.n]

    idx = load_index(atlas_root())
    alias, _ = build_alias_map(idx)
    cache = json.loads(CACHE.read_text()) if CACHE.exists() else {}

    out = []
    for i, r in enumerate(sample, 1):
        spans = matched_forms(r["sentence"], r["identifier"], alias)
        out.append({
            "n": i, "pmid": r["pmid"], "identifier": r["identifier"],
            "entity": r.get("entity", ""),
            "matched_span": "|".join(sorted(spans, key=len, reverse=True)),
            "authority_name": mesh_label(r["identifier"], cache),
            "verdict": "",
            "sentence": r["sentence"].replace("\n", " "),
        })
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    CACHE.write_text(json.dumps(cache, indent=0, sort_keys=True) + "\n")

    dest = Path(args.out)
    dest.parent.mkdir(parents=True, exist_ok=True)
    with dest.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLUMNS)
        w.writeheader()
        w.writerows(out)

    no_span = sum(1 for r in out if not r["matched_span"])
    no_label = sum(1 for r in out if not r["authority_name"])
    print(f"wrote {dest}: {len(out)} mentions from the {args.stratum} stratum, "
          f"seed {args.seed}")
    print(f"  without a recovered span : {no_span} "
          f"(the alias may have been filtered since the sample was built)")
    print(f"  without an NLM label     : {no_label} "
          f"(gene/OMIM identifiers have no MeSH descriptor)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
