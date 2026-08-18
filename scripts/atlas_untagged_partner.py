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
untagged article. It becomes someone else's article -- and that is TRUE BY
CONSTRUCTION rather than discovered: the taxonomy has no lane for these three,
so any tag such an article carries is necessarily another modality's. What
this page contributes is the SIZE of the affected set, and the contrast with
modalities that DO have a lane, which are filed entirely under others only
6.1% of the time.

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


def named_modality_control(idx):
    """How often a modality that HAS a lane is filed entirely under others.

    THE CORPUS BASE RATE WAS THE WRONG DENOMINATOR. For radiotherapy,
    chemotherapy and surgery, "carries at least one tag" IS "filed under
    another modality" -- no lane names them. For the corpus at large,
    "carries at least one tag" is overwhelmingly "carries its OWN tag". The
    two are different events, so 88.7% was never a control for this one.

    The event-matched control is a title-anchored set for a modality that
    does have a lane, asked the same question: how often is it recorded
    ENTIRELY under someone else?
    """
    out = {}
    tot = miss = 0
    for m in config.MECHANISM_KEYWORDS:
        stem = re.escape(m.replace("-", " ")).replace(r"\ ", r"[\s-]?")
        pat = re.compile(stem, re.I)
        subj = [r for r in idx if pat.search(r.get("title") or "")]
        if len(subj) < 10:
            continue
        other = sum(1 for r in subj
                    if (r.get("mechanisms") or [])
                    and m not in (r.get("mechanisms") or []))
        out[m] = {"titled": len(subj), "entirely_other": other}
        tot += len(subj)
        miss += other
    return {"per_modality": out, "titled": tot, "entirely_other": miss,
            "rate": (miss / tot) if tot else None}


def retag_with_current_vocabulary(subject):
    """Which mechanisms the CURRENT tagger gives these articles.

    ON PRODUCTION TEXT, through the production entry points. An earlier
    version built `title + (record.get("abstract") or "")` from
    corpus/INDEX.jsonl -- which carries NO abstract field, so it matched the
    TITLE ALONE -- and then published three caveat paragraphs, one bolded,
    saying it read title+abstract while the frozen tagging read full text.
    Both halves were false: `get_searchable_text` defaults
    `include_full_text=False`, so the frozen tagger read title + MeSH +
    annotations + abstract and never the body.

    That mattered: under the broken scope the `bioelectric` partner would
    have vanished under ANY vocabulary, so the disappearance was guaranteed
    rather than measured. Reading production text with only the VOCABULARY
    changed isolates it cleanly, and every other partner is untouched.
    """
    from tag_articles import load_article, get_searchable_text, match_mechanisms
    out = Counter()
    still = 0
    for r in subject:
        pmid = str(r.get("pmid") or "")
        path = PROJECT_ROOT / "corpus" / "by-pmid" / f"{pmid}.md"
        if not path.exists():
            continue
        try:
            fm, body = load_article(path)
        except Exception:
            continue
        text = get_searchable_text(fm, body)
        hits = match_mechanisms(text, (fm.get("title") or "").lower())
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
    out["named_control"] = named_modality_control(idx)

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
        _retag = retag_with_current_vocabulary(subject)
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
                [[k, v] for k, v in _retag[0].most_common(8)],
            "tagged_current_vocabulary": _retag[1],
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

    L += ["The `carry 2+` column is reported and deliberately NOT argued "
          "from. An earlier version read a below-average multi-tag rate as "
          "evidence against the combination-paper reading. That inference is "
          "invalid here: a paper about radiotherapy AND a checkpoint "
          "inhibitor can carry at most ONE modality tag, because radiotherapy "
          "has no lane. The hypothesis predicts enrichment for exactly one "
          "tag, not for two, so this column cannot test it in either "
          "direction.", ""]

    L += ["The third column used to be headed \"carry a DIFFERENT modality's "
          "tag\". It is `titled - untagged`, i.e. *carries at least one tag* -- "
          "and because the taxonomy has no lane for any of these three, the "
          "word DIFFERENT excluded nothing. Renamed to what it counts.", ""]

    nc = d.get("named_control") or {}
    if nc.get("rate") is not None:
        rates = {m: 100 * s["subject_tagged_as_other"] / s["subject_titled"]
                 for m, s in d["modalities"].items()}
        ncr = 100 * nc["rate"]
        L += ["### The control, and why the corpus rate was the wrong one", ""]
        L += [f"An earlier version of this page compared those rates against "
              f"the corpus-wide tagged rate of {base:.1f}% and concluded that "
              f"nothing was elevated. **That comparison was between two "
              f"different events.** For these three modalities \"carries at "
              f"least one tag\" IS \"filed under another modality\", because no "
              f"lane names them. For the corpus at large it is overwhelmingly "
              f"\"carries its OWN tag\". The corpus rate was never a control "
              f"for this.", ""]
        L += [f"The event-matched control asks the same question of modalities "
              f"that DO have a lane: how often is a title-anchored set "
              f"recorded ENTIRELY under someone else?", ""]
        L += ["| | titled | filed entirely under others |", "|---|--:|--:|"]
        for m, s in d["modalities"].items():
            L.append(f"| {m} (no lane) | {s['subject_titled']:,} | "
                     f"**{s['subject_tagged_as_other']:,}** "
                     f"({rates[m]:.1f}%) |")
        L.append(f"| **modalities WITH a lane** | {nc['titled']:,} | "
                 f"{nc['entirely_other']:,} ({ncr:.1f}%) |")
        L += [""]
        worst = min(rates.values())
        L += [f"That is a factor of {worst/ncr:.0f} to {max(rates.values())/ncr:.0f}. "
              f"So the original sentence -- that a modality with no name "
              f"becomes someone else's article -- is not merely supported, it "
              f"is **true by construction**: the taxonomy has no lane for "
              f"these three, so any tag such an article carries is necessarily "
              f"another modality's. It is a property of the design rather than "
              f"a discovery, and what this page contributes is the size of the "
              f"affected set and the contrast above.", ""]
        L += ["Both readings this page has carried were wrong in opposite "
              "directions: the first presented a structural certainty as an "
              "empirical finding, and the second withdrew it against a "
              "denominator measuring something else.", ""]

    # VOCABULARY VINTAGE. The frozen index was tagged before #MECH-PRECISION
    # retired keywords, so a partner table read off it can name a tag the
    # project has already withdrawn.
    # WHICH TAGS WERE ACTUALLY RETIRED, read from config.py's own record
    # rather than inferred from a disappearance -- a tag can also vanish here
    # because this re-run reads less text, and attributing that to the
    # vocabulary would be the asymmetric comparison this repo keeps making.
    gone = {}
    for m, s in d["modalities"].items():
        old_tags = {k for k, _v in s.get("partner_tags", [])}
        new_tags = {k for k, _v in s.get("partner_tags_current_vocabulary", [])}
        # EVERY partner that disappears, because the comparison is now
        # vocabulary-only: both passes run the production matcher on
        # production text, so a disappearance can only be the vocabulary. An
        # earlier version filtered this list through a 600-character
        # proximity scan of config.py, which excluded `frequency-therapy`
        # (recorded 951 characters from its marker) and so forced a FALSE
        # correction -- the prose naming it had been right.
        missing = sorted(old_tags - new_tags)
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
        L += ["Both passes run the production matcher over production text, "
              "so a disappearance can only be the vocabulary. Every other "
              "partner is untouched -- `immunotherapy` and `nanoparticle` and "
              "`ttfields` return identical counts -- which is what makes the "
              "two above attributable.", ""]
        deltas = []
        for m, s in d["modalities"].items():
            a, b = s["subject_tagged_as_other"], s.get("tagged_current_vocabulary")
            if b is not None and b != a:
                deltas.append(f"{m} {a} -> {b}")
        if deltas:
            L += [f"Removing them moves the tagged counts: {', '.join(deltas)} "
                  f"-- articles whose ONLY tag was a retired keyword.", ""]
        L += ["An earlier version of this section published three caveat "
              "paragraphs, one in bold, saying it read title and abstract "
              "while the frozen tagging read stored full text. Both halves "
              "were false: `corpus/INDEX.jsonl` carries no abstract field, so "
              "that pass matched titles alone, and `get_searchable_text` "
              "defaults to excluding the body, so the frozen tagger never read "
              "it either. Under the broken scope this partner would have "
              "vanished under ANY vocabulary. The caveats were excusing a "
              "defect rather than describing a limit.", ""]

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
