#!/usr/bin/env python3
"""What happens to a modality the taxonomy has no name for? (#723)

THE QUESTION, AND WHY THE OBVIOUS ANSWER IS WRONG
--------------------------------------------------
Radiotherapy has no mechanism tag, no query and no engine term, and the easy
conclusion is that the corpus does not contain it. That conclusion is wrong, and
the way it is wrong is the finding.

Radiotherapy IS in the frozen corpus. What happens to it is that it arrives
attached to a modality the taxonomy CAN name -- radiotherapy plus checkpoint
blockade, radiotherapy plus tumour-treating fields -- and the article is then
filed under the partner. The modality without a name does not become an
untagged article. It becomes someone else's article.

That is a different failure from a coverage gap and it has a different fix. A
coverage gap is closed by retrieving more. This is closed by naming, and until
it is named every co-occurrence measurement in the project attributes the whole
paper to the partner.

WHAT THIS MEASURES
------------------
For each candidate untagged modality:
  * how often it is the SUBJECT (title-anchored, high precision)
  * how often it is substantially discussed (full text, high recall)
  * for the subject-anchored set, what tags those articles actually carry

The third is the point. A high untagged count would mean a plain coverage gap; a
high OTHER-TAG count means attribution to a partner, which is worse because it
is invisible in every downstream count rather than merely missing.

WHAT IT DOES NOT DO
-------------------
It does not add a tag. The candidate keyword sets here are measured for
precision, not installed, because a lane that fires broadly would push its
articles' existing partner tags around and change numbers the manuscript quotes.
Deciding to install one is a separate act with its own review.

Usage:
    python scripts/atlas_untagged_partner.py
"""

import json
import re
from collections import Counter
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
INDEX = PROJECT_ROOT / "corpus" / "INDEX.jsonl"
FULLTEXT = PROJECT_ROOT / "corpus" / "by-pmid"
OUT_MD = PROJECT_ROOT / "analysis" / "atlas-untagged-partner.md"
OUT_JSON = PROJECT_ROOT / "analysis" / "atlas-untagged-partner.json"

# Modalities with no mechanism tag. Patterns are deliberately specific: a loose
# pattern would inflate the subject count, which is the number the argument
# rests on.
CANDIDATES = {
    "radiotherapy": re.compile(
        r"radiotherap|radiation therapy|chemoradi|\birradiat|"
        r"stereotactic body|\bSBRT\b|\bIMRT\b|brachytherap|"
        r"proton therapy|carbon.ion|radiosensiti", re.I),
    "chemotherapy": re.compile(
        r"chemotherap|cytotoxic agent|platinum-based|\bcisplatin\b|"
        r"\bcarboplatin\b|\bdoxorubicin\b|\bpaclitaxel\b|\bgemcitabine\b|"
        r"\b5-fluorouracil\b|\bfolfox\b|\bfolfirinox\b", re.I),
    "surgery": re.compile(
        r"\bsurgical resection\b|\bsurgery\b|\bresection\b|\bmastectom|"
        r"\bhepatectom|\bpneumonectom|\blymphadenectom|\bcytoreducti", re.I),
}

# Substantial discussion, not a passing citation. Chosen before looking at the
# result and stated so the number can be re-derived at a different threshold.
STRONG_MENTIONS = 5


def load_index() -> list:
    return [json.loads(l) for l in INDEX.read_text().splitlines() if l.strip()]


def scan() -> dict:
    idx = load_index()
    by_pmid = {str(r.get("pmid")): r for r in idx}
    out = {"corpus_size": len(idx), "strong_threshold": STRONG_MENTIONS,
           "modalities": {}}

    fulltext_counts = {name: Counter() for name in CANDIDATES}
    for p in sorted(FULLTEXT.glob("*.md")):
        text = p.read_text(errors="ignore")
        pmid = p.stem
        for name, pat in CANDIDATES.items():
            n = len(pat.findall(text))
            if n:
                fulltext_counts[name][pmid] = n

    for name, pat in CANDIDATES.items():
        subject = [r for r in idx if pat.search(r.get("title") or "")]
        strong = [p for p, n in fulltext_counts[name].items()
                  if n >= STRONG_MENTIONS]
        tags = Counter()
        untagged = 0
        for r in subject:
            ms = r.get("mechanisms") or []
            if not ms:
                untagged += 1
            for m in ms:
                tags[m] += 1
        out["modalities"][name] = {
            "subject_titled": len(subject),
            "strong_fulltext": len(strong),
            "any_fulltext": len(fulltext_counts[name]),
            "subject_untagged": untagged,
            "subject_tagged_as_other": len(subject) - untagged,
            # A LIST of pairs, not a dict: the artifact is written with
            # sort_keys=True, which reorders a dict alphabetically and would
            # make `antibody-drug-conjugate` read as the dominant partner
            # instead of `immunotherapy`. Rank is the finding here, so it must
            # survive serialisation.
            "partner_tags": [[k, v] for k, v in tags.most_common(8)],
        }
    return out


def render(d: dict) -> str:
    n = d["corpus_size"]
    L = ["# What happens to a modality the taxonomy cannot name", ""]
    L += ["*Generated by `scripts/atlas_untagged_partner.py`. Every figure is "
          "recomputed.*", ""]

    L += [f"Radiotherapy, chemotherapy and surgery have no mechanism tag in this "
          f"project. The easy conclusion is that the {n:,}-article frozen corpus "
          f"does not contain them. It does, and what happens to them instead is "
          f"the finding.", ""]

    L += ["| modality | titled about it | 5+ mentions in full text | any mention |",
          "|---|--:|--:|--:|"]
    for m, s in d["modalities"].items():
        L.append(f"| {m} | {s['subject_titled']:,} ({100*s['subject_titled']/n:.1f}%) "
                 f"| {s['strong_fulltext']:,} ({100*s['strong_fulltext']/n:.1f}%) "
                 f"| {s['any_fulltext']:,} ({100*s['any_fulltext']/n:.1f}%) |")
    L += [""]

    L += ["## They are not untagged. They are filed under their partner.", ""]
    L += ["For the articles that are titled about each modality -- the "
          "high-precision set -- here is what the taxonomy actually recorded:", ""]
    L += ["| modality | titled | carry NO tag | carry a DIFFERENT modality's tag |",
          "|---|--:|--:|--:|"]
    for m, s in d["modalities"].items():
        L.append(f"| {m} | {s['subject_titled']:,} | {s['subject_untagged']:,} | "
                 f"**{s['subject_tagged_as_other']:,}** |")
    L += [""]
    for m, s in d["modalities"].items():
        if not s["partner_tags"]:
            continue
        top = ", ".join(f"`{k}` {v}" for k, v in s["partner_tags"][:4])
        L.append(f"* **{m}** is recorded as: {top}")
    L += [""]

    L += ["A modality with no name does not become an untagged article. It "
          "becomes **someone else's article**. Every co-occurrence, prevalence "
          "and capture figure the project computes attributes the whole paper "
          "to the partner that happened to have a tag.", ""]

    L += ["## Why this is a different problem from a coverage gap", ""]
    L += ["A coverage gap is closed by retrieving more literature. This is "
          "closed by naming, and it is worse in one specific way: a missing "
          "article is absent from a count, which a careful reader can suspect, "
          "whereas a misattributed article is PRESENT in the wrong count, which "
          "nothing signals.",
          "",
          "It also means the corpus is a poor witness to its own breadth. The "
          "immunotherapy share includes combination papers whose other half the "
          "taxonomy cannot see.", ""]

    L += ["## What this does not do", ""]
    L += ["* It does not install a tag. The patterns here are measured, not "
          "adopted: a new lane would move articles' existing partner tags and "
          "change figures the manuscript quotes, so adopting one is a separate "
          "act with its own review.",
          "* The title-anchored count is high precision and low recall by "
          "design; the full-text count is the reverse. Both are given because "
          "neither alone answers 'how much of this is here'.",
          f"* The {d['strong_threshold']}-mention threshold was fixed before "
          "the result was seen. It is a judgement, and the any-mention column "
          "is there so a reader can apply their own.",
          "* Presence in the corpus is not evidence about the field. The corpus "
          "was retrieved with 33 keyword queries, none of them about these "
          "modalities, so what is here arrived attached to something else -- "
          "which is exactly the mechanism this measures.",
          ""]
    return "\n".join(L) + "\n"


def main():
    d = scan()
    if not d["modalities"] or all(v["subject_titled"] == 0
                                  for v in d["modalities"].values()):
        raise SystemExit(
            "no modality had a single titled article, which is not a finding -- "
            "it is what a broken index path or pattern looks like.")
    OUT_JSON.write_text(json.dumps(d, indent=1, sort_keys=True) + "\n",
                        encoding="utf-8")
    OUT_MD.write_text(render(d), encoding="utf-8")
    print(f"wrote {OUT_MD}")
    print(f"wrote {OUT_JSON}")
    for m, s in d["modalities"].items():
        print(f"  {m:14s} titled {s['subject_titled']:>4}  "
              f"of which filed under another modality "
              f"{s['subject_tagged_as_other']:>4}")


if __name__ == "__main__":
    main()
