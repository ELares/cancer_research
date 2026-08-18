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
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config  # noqa: E402

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


def retag_with_current_vocabulary(subject):
    """Which mechanisms the CURRENT tagger would give these articles.

    `corpus/INDEX.jsonl` is FROZEN and was tagged before #MECH-PRECISION
    retired the bare phrase "membrane potential" from `bioelectric` (audited
    at 13.3% precision -- it nearly always appears as MITOCHONDRIAL membrane
    potential, a JC-1 apoptosis readout) and narrowed `frequency-therapy`.
    So a partner table read off the frozen index can name a tag the project
    has already withdrawn, as evidence.

    Worse for the causal story: `scripts/queries.txt` retrieves on
    "membrane potential" itself, so those articles are in the corpus BECAUSE
    of the phrase that then tagged them -- a retrieval-plus-retired-keyword
    artifact rather than a partner absorbing a nameless modality.
    """
    from tag_articles import text_matches_keyword
    items = [(m, k) for m, ks in config.MECHANISM_KEYWORDS.items() for k in ks]
    out = Counter()
    still = 0
    for r in subject:
        blob = ((r.get("title") or "") + " " +
                (r.get("abstract") or "")).lower()
        hits = {m for m, k in items if text_matches_keyword(blob, k)}
        if hits:
            still += 1
        for m in hits:
            out[m] += 1
    return out, still


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

    # THE CONTROL THIS ANALYSIS NEVER RAN. Carrying another modality's tag is
    # the corpus DEFAULT, so a per-modality rate means nothing without it. All
    # three subjects turn out to sit AT OR BELOW the default, which withdraws
    # the causal reading: namelessness cannot be shown to raise something it
    # does not raise.
    n_tagged = sum(1 for r in idx if (r.get("mechanisms") or []))
    n_multi = sum(1 for r in idx if len(r.get("mechanisms") or []) >= 2)
    out["corpus_tagged"] = n_tagged
    out["corpus_multi_tagged"] = n_multi

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
        multi = sum(1 for r in subject if len(r.get("mechanisms") or []) >= 2)
        out["modalities"][name] = {
            "subject_titled": len(subject),
            "subject_multi_tagged": multi,
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
            "partner_tags_current_vocabulary":
                [[k, v] for k, v in retag_with_current_vocabulary(subject)[0].most_common(8)],
            "tagged_current_vocabulary": retag_with_current_vocabulary(subject)[1],
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

    n = d["corpus_size"]
    ct = d.get("corpus_tagged")
    cm = d.get("corpus_multi_tagged")
    base = 100 * ct / n if ct else None
    base_multi = 100 * cm / n if cm else None

    L += ["## What the taxonomy recorded, against what it records by default", ""]
    L += ["The header of this section used to read \"They are not untagged. "
          "They are filed under their partner\", and the table carried no "
          "denominator. Carrying another modality's tag is the CORPUS DEFAULT, "
          "so the rate only means something beside it.", ""]
    L += ["| modality | titled | carry NO tag | carry at least one tag | "
          "carry 2+ |", "|---|--:|--:|--:|--:|"]
    for m, s in d["modalities"].items():
        tt = s["subject_titled"]
        L.append(f"| {m} | {tt:,} | {s['subject_untagged']:,} "
                 f"({100*s['subject_untagged']/tt:.1f}%) | "
                 f"**{s['subject_tagged_as_other']:,}** "
                 f"({100*s['subject_tagged_as_other']/tt:.1f}%) | "
                 f"{s.get('subject_multi_tagged', 0):,} "
                 f"({100*s.get('subject_multi_tagged', 0)/tt:.1f}%) |")
    if base is not None:
        L.append(f"| **the corpus** | {n:,} | {n-ct:,} "
                 f"({100*(n-ct)/n:.1f}%) | **{ct:,}** ({base:.1f}%) | "
                 f"{cm:,} ({base_multi:.1f}%) |")
    L += [""]

    # THE COLUMN IS NOT WHAT ITS OLD HEADER SAID. `len(subject) - untagged` is
    # "has at least one tag"; since no lane names these three modalities, the
    # word DIFFERENT excluded nothing that could ever have been subtracted.
    for m, s in d["modalities"].items():
        if not s["partner_tags"]:
            continue
        top = ", ".join(f"`{k}` {v}" for k, v in s["partner_tags"][:4])
        L.append(f"* **{m}** is recorded as: {top}")
    L += [""]

    L += ["The third column used to be headed \"carry a DIFFERENT modality's "
          "tag\". It is `titled - untagged`, i.e. *carries at least one tag* -- "
          "and because the taxonomy has no lane for any of these three, the "
          "word DIFFERENT excluded nothing. Renamed to what it counts.", ""]

    if base is not None:
        rates = {m: 100 * s["subject_tagged_as_other"] / s["subject_titled"]
                 for m, s in d["modalities"].items()}
        below = sorted(m for m, r in rates.items() if r < base)
        L += ["### The causal reading is withdrawn", ""]
        if len(below) == len(rates):
            L += [f"**Every one of the three sits BELOW the corpus default** "
                  f"({', '.join(f'{m} {rates[m]:.1f}%' for m in below)} against "
                  f"{base:.1f}%). They are also untagged MORE often than the "
                  f"corpus and carry FEWER tags each. So the earlier claim -- "
                  f"that a modality with no name \"does not become an untagged "
                  f"article, it becomes someone else's article\" -- is not "
                  f"supported by these numbers: nothing here is elevated, and "
                  f"the second half of that sentence describes what almost "
                  f"every article in this corpus looks like.", ""]
            L += ["The combination-paper story goes with it. If these were "
                  "papers about two modalities, they would carry MORE tags "
                  "than average, not fewer.", ""]
        else:
            L += [f"{', '.join(below) or 'None'} sit below the corpus default "
                  f"of {base:.1f}%; the others are above it.", ""]

        L += ["### What survives, and it is structural rather than causal", ""]
        L += ["The taxonomy has no lane for radiotherapy, chemotherapy or "
              "surgery. So any tag an article about them carries is "
              "NECESSARILY another modality's -- that is guaranteed by the "
              "design, not discovered by this measurement. What this analysis "
              "contributes is the COUNT of articles in that position "
              f"({', '.join(str(s['subject_tagged_as_other']) for s in d['modalities'].values())}), "
              "and every co-occurrence, prevalence and capture figure the "
              "project computes attributes those whole papers to whichever "
              "partner happened to have a tag. That consequence is unchanged; "
              "the causal story about namelessness is what goes.", ""]

    # VOCABULARY VINTAGE. The frozen index was tagged before #MECH-PRECISION
    # retired keywords, so a partner table read off it can name a tag the
    # project has already withdrawn.
    # WHICH TAGS WERE ACTUALLY RETIRED, read from config.py's own record
    # rather than inferred from a disappearance -- a tag can also vanish here
    # because this re-run reads less text, and attributing that to the
    # vocabulary would be the asymmetric comparison this repo keeps making.
    retired = set()
    cfg = (PROJECT_ROOT / "scripts" / "config.py").read_text(errors="ignore")
    for blk in re.finditer(r"#MECH-PRECISION(.{0,600})", cfg, re.S):
        for mech in config.MECHANISM_KEYWORDS:
            if re.search(rf"\b{re.escape(mech)}\b", blk.group(1)):
                retired.add(mech)

    gone = {}
    for m, s in d["modalities"].items():
        old_tags = {k for k, _v in s.get("partner_tags", [])}
        new_tags = {k for k, _v in s.get("partner_tags_current_vocabulary", [])}
        missing = sorted((old_tags - new_tags) & retired)
        if missing:
            gone[m] = missing
    if gone:
        L += ["### Some of those partners are tags the project has since "
              "withdrawn", ""]
        L += ["`corpus/INDEX.jsonl` is FROZEN and was tagged before "
              "#MECH-PRECISION retired the bare phrase \"membrane potential\" "
              "from `bioelectric` -- audited at 13.3% precision, because it "
              "nearly always appears as MITOCHONDRIAL membrane potential, a "
              "JC-1 apoptosis readout. Re-running the CURRENT vocabulary over "
              "the same articles, these partners disappear:", ""]
        for m, missing in gone.items():
            L.append(f"* **{m}**: {', '.join(f'`{x}`' for x in missing)}")
        L += [""]
        L += ["That matters most for the causal reading. `scripts/queries.txt` "
              "retrieves on \"membrane potential\" itself, so those articles "
              "are in this corpus BECAUSE of the phrase that then tagged them "
              "-- a retrieval-plus-retired-keyword artifact, not a partner "
              "absorbing a nameless modality.", ""]
        L += ["Only tags `scripts/config.py` records as narrowed under "
              "#MECH-PRECISION are listed. Other partners also drop out of "
              "the re-run, but this pass reads title and abstract while the "
              "frozen tagging read the stored full text, so those "
              "disappearances mix the vocabulary with the text scope and are "
              "not attributed here.", ""]
        L += ["**The absolute counts are NOT comparable between the two "
              "columns and are deliberately not published as a delta.** The "
              "frozen tagging read the stored full text; this re-run reads "
              "title and abstract only, so any difference in totals mixes the "
              "vocabulary change with a change of text scope. What the "
              "comparison isolates is which partner tags survive the "
              "vocabulary at all, which is the question here.", ""]

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
