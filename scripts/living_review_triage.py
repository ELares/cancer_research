#!/usr/bin/env python3
"""Turn a living-review delta into an answer, instead of an unread artifact.

WHY THIS EXISTS
---------------
`.github/workflows/living-review.yml` re-runs the committed mechanism queries
against PubMed every month and uploads a dated delta as a build artifact. It
deliberately never commits: the frozen corpus is frozen, and a human decides
whether anything in the window matters.

Nobody decided. Two windows had accumulated unread -- 1,855 new records in the
August one alone, with 24 landmark detections -- because the output is a 4.5 MB
zip and the question it answers ("does any of this change a committed claim?")
takes an hour of reading to reach. A monitoring mechanism whose output nobody
opens is not monitoring anything.

This reads a delta and answers that question on one page.

THE CONFOUND IT EXISTS TO PREVENT, found on the first delta triaged
--------------------------------------------------------------------
The obvious triage is to compare the new window's mechanism mix against the
frozen corpus and report what moved. Run naively that reports `cuproptosis` at
0.0% of the frozen corpus against 4.0% of the new window -- an apparently
explosive emergence.

It is an artifact. The frozen `corpus/INDEX.jsonl` was tagged with an OLDER
mechanism list than the living review uses. 38 frozen records name cuproptosis
somewhere -- 33 of them under `pathway_targets`, a tag axis maintained
separately -- while carrying `mechanisms` values like `['immunotherapy']`;
disulfidptosis is worse, at 87. The keyword list gained both after the freeze
(#347, added as a scaffold with the frozen results deliberately untouched), so
the frozen side never had the chance to be tagged for them. The corpus already
RECOGNISED these concepts; only the mechanism vocabulary is stale.

So every mechanism is classified before any share is compared:

  COMPARABLE   present in the frozen tag vocabulary -- a share change is a
               real signal about the literature;
  NOT COMPARABLE  absent from it, so the frozen share is structurally zero and
               the comparison is meaningless in that direction.

The not-comparable ones are reported separately, never as movement.

WHAT IT WILL NOT DO
-------------------
Decide. It flags what a human should look at and states what it cannot settle:
a five-week window against a decade-long corpus is not a like-for-like base
rate even for comparable mechanisms, query date-windows differ between mechanisms
(several are restricted to 2020:2026 or later), and PubMed indexing lag
under-represents the newest literature in both directions.

Usage:
    python scripts/living_review_triage.py --delta path/to/index.jsonl
    python scripts/living_review_triage.py --delta ... --out analysis/...md
"""

import argparse
import collections
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from config import PROJECT_ROOT  # noqa: E402

FROZEN = PROJECT_ROOT / "corpus" / "INDEX.jsonl"
# Publication types the repo already treats as high-evidence elsewhere.
HIGH_EVIDENCE = {"Randomized Controlled Trial", "Clinical Trial, Phase III",
                 "Clinical Trial, Phase II", "Meta-Analysis"}


def load_jsonl(p: Path) -> list:
    return [json.loads(l) for l in p.open() if l.strip()]


def mech_counts(rows: list):
    c = collections.Counter()
    for r in rows:
        for m in (r.get("mechanisms") or []):
            c[m] += 1
    return c


def frozen_vocabulary(rows: list) -> set:
    """Mechanisms the frozen index was actually TAGGED with.

    Not the current config's mechanism list -- the vocabulary as it stood when
    the corpus was tagged. A mechanism outside this set has a structurally zero
    frozen share, and comparing against it measures the tagging vintage rather
    than the literature.
    """
    return set(mech_counts(rows))


# The frozen index carries NO abstract and NO mesh_terms -- only a title plus
# several independently-maintained tag axes. Searching `abstract` there returns
# nothing and looks like a small count rather than a missing field, which is how
# the first version of this undercounted cuproptosis 38 -> 5.
FROZEN_TEXT_FIELDS = ("title", "pathway_targets", "biology_processes",
                      "cancer_types", "resistant_states")


def untagged_mentions(rows: list, mech: str) -> tuple:
    """Frozen records naming a mechanism they are not tagged with, and where.

    The evidence that a zero frozen share is a VOCABULARY gap rather than an
    absence: if the corpus holds records about the thing and simply never
    labelled them under `mechanisms`, the zero says nothing about the
    literature. Returns the count and the per-field breakdown, because "38 of
    them under pathway_targets" is a much stronger statement than a bare 38 --
    it shows the corpus already recognised the concept on a DIFFERENT axis.
    """
    needle = mech.replace("-", " ").lower()
    where = collections.Counter()
    n = 0
    for r in rows:
        if mech in (r.get("mechanisms") or []):
            continue
        hit = None
        for f in FROZEN_TEXT_FIELDS:
            v = r.get(f)
            if v is None:
                continue
            blob = (" ".join(v) if isinstance(v, list) else str(v)).lower()
            if needle in blob or mech.lower() in blob:
                hit = f
                break
        if hit:
            n += 1
            where[hit] += 1
    return n, where


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--delta", required=True,
                    help="index.jsonl from a living-review artifact")
    ap.add_argument("--out", default=None, help="write a report here")
    ap.add_argument("--top", type=int, default=12)
    args = ap.parse_args()

    delta_p = Path(args.delta)
    if not delta_p.exists():
        print(f"no delta at {delta_p}", file=sys.stderr)
        return 1
    live = load_jsonl(delta_p)
    froz = load_jsonl(FROZEN)

    lc, fc = mech_counts(live), mech_counts(froz)
    vocab = frozen_vocabulary(froz)
    l_tot, f_tot = sum(lc.values()) or 1, sum(fc.values()) or 1

    comparable, incomparable = [], []
    for m in sorted(set(lc) | set(fc), key=lambda k: -lc.get(k, 0)):
        row = {"mechanism": m, "frozen": fc.get(m, 0), "living": lc.get(m, 0),
               "frozen_share": fc.get(m, 0) / f_tot,
               "living_share": lc.get(m, 0) / l_tot}
        if m in vocab:
            comparable.append(row)
        else:
            n, where = untagged_mentions(froz, m)
            row["frozen_untagged_mentions"] = n
            row["frozen_untagged_where"] = dict(where)
            incomparable.append(row)

    # High-evidence records in the window, which is what would actually move a
    # maturity claim rather than a volume one.
    high = [r for r in live
            if set(r.get("pub_types") or []) & HIGH_EVIDENCE]
    high_by_mech = collections.Counter()
    for r in high:
        for m in (r.get("mechanisms") or []):
            high_by_mech[m] += 1

    moved = sorted((r for r in comparable if r["frozen"] >= 20),
                   key=lambda r: -(r["living_share"] / max(r["frozen_share"], 1e-9)))

    L = [
        f"# Living-review triage — {delta_p.parent.name}", "",
        "Generated by `scripts/living_review_triage.py`. The living-review",
        "workflow uploads a delta every month and never commits; this reads one",
        "and says what a human should look at.", "",
        f"**{len(live):,} new records.** {len(high):,} carry a high-evidence",
        "publication type (randomised trial, phase II/III, meta-analysis).", "",
        "## The comparison that is NOT valid, and why it is reported separately", "",
        f"{len(incomparable)} mechanism(s) tagged in this window are absent from the",
        "frozen corpus's tag vocabulary entirely. Their frozen share is",
        "structurally zero, so any 'growth' against it measures when the keyword",
        "list changed rather than what the literature did.", "",
    ]
    if incomparable:
        L += ["| mechanism | new records | frozen `mechanisms` tag | frozen records naming it elsewhere |",
              "|---|--:|--:|---|"]
        for r in incomparable:
            w = r.get("frozen_untagged_where") or {}
            detail = ", ".join(f"{n} in `{f}`" for f, n in
                               sorted(w.items(), key=lambda kv: -kv[1])) or "none"
            L.append(f"| {r['mechanism']} | {r['living']} | {r['frozen']} | "
                     f"**{r['frozen_untagged_mentions']}** ({detail}) |")
        L += ["",
              "The last column is the evidence that these zeros are a vocabulary",
              "gap and not an absence. The frozen corpus already recognises these",
              "concepts on OTHER tag axes -- it simply never labelled them under",
              "`mechanisms`, because that keyword list gained them after the",
              "corpus was frozen. Note the frozen index carries no abstract and no",
              "MeSH terms, so this searches its title and tag axes only; a naive",
              "search of `abstract` there returns nothing and reads as a small",
              "count rather than a missing field.", ""]
    else:
        L += ["None this window -- every tagged mechanism is in the frozen",
              "vocabulary, so the shares below are all comparable.", ""]

    L += [
        "## Comparable mechanisms, by share change", "",
        "Restricted to mechanisms with at least 20 frozen records, because a",
        "ratio on a handful of records is noise.", "",
        "| mechanism | frozen share | window share | new records |",
        "|---|--:|--:|--:|",
    ]
    for r in moved[:args.top]:
        L.append(f"| {r['mechanism']} | {100*r['frozen_share']:.1f}% | "
                 f"{100*r['living_share']:.1f}% | {r['living']} |")

    L += [
        "", "## Where the high-evidence records landed", "",
        "A volume shift is weaker evidence than a maturity shift. These are the",
        "mechanisms picking up trial-grade records in this window:", "",
        "| mechanism | high-evidence records |", "|---|--:|",
    ] + [f"| {m} | {n} |" for m, n in high_by_mech.most_common(10)] + [
        "", "## What this cannot settle", "",
        "* A five-week window against a decade of corpus is not a like-for-like",
        "  base rate, even for comparable mechanisms.",
        "* Query date-windows differ by mechanism -- several are restricted to",
        "  2020:2026 or later -- so the frozen shares are not on one clock either.",
        "* PubMed indexing lag under-represents the newest literature, biasing",
        "  every window's counts down.", "",
        "So read this as a triage: it says where to look, not what is true.", "",
    ]
    text = "\n".join(L) + "\n"
    if args.out:
        Path(args.out).write_text(text, encoding="utf-8")
        print(f"wrote {args.out}")
    else:
        print(text)
    print(f"  {len(live):,} records, {len(high):,} high-evidence, "
          f"{len(incomparable)} not-comparable mechanism(s)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
