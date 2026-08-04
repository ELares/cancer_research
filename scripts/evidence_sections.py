#!/usr/bin/env python3
"""Section-scoped full-text extraction for evidence-tier tagging (#TAGGER-V2).

WHY THIS EXISTS
---------------
The production evidence tagger reads title + MeSH + entity annotations +
abstract. It never reads the 216 MB of full text the corpus already stores, so
it misses evidence-design statements that only ever appear in a Methods
section, and its binary evidence-detection recall sits at 55%.

Naively switching full text on lifts recall to ~97% but *drops* precision from
96% to 89%, and the damage is one-directional: preclinical articles get
promoted into clinical tiers. The cause is that Introduction and Discussion
sections describe OTHER people's studies. PMID 40700574 is a mouse study whose
Discussion says "The most recent mRNA vaccine to undergo a phase 3 clinical
trial ...", and that single sentence promoted it to phase3-clinical.

The project's own labeling guideline (analysis/evidence-labeling-guidelines-v2.md)
already states the correct rule: take the tier from the article's primary
research, "not in cited references or background discussion". This module
implements that rule mechanically by splitting the '## Full Text' blob into
sections and classifying each one:

  SELF  - the article describes its OWN design/results (Methods, Results, ...)
  CITED - the article discusses OTHER work (Introduction, Discussion, ...)
  DROP  - boilerplate (References, Funding, Ethics statements, ...)

Evidence-tier detection reads SELF only.

HOW SECTIONS ARE FOUND
----------------------
Corpus full text has no markdown heading syntax; headings survive as bare,
blank-line-surrounded short lines ("Materials and methods", "Discussion").
Measured on a 600-article random sample, this shape finds headings in 98.7% of
articles. Headings vary wildly across 803 journals, so unrecognised headings
INHERIT the enclosing class -- a "Cell culture" or "Tumor xenograft model"
subsection stays SELF without needing to be enumerated. Text before the first
recognised heading is SELF (title/abstract front matter).

Offline, deterministic, stdlib only.
"""

import re

SELF = "self"
CITED = "cited"
DROP = "drop"

# Headings whose content is the article's OWN work.
_SELF_PAT = re.compile(
    r"^(materials?( and | & )methods?|methods?|methodology|experimental( section| procedures?| methods?| design)?"
    r"|patients? and methods?|subjects? and methods?|methods? and materials?"
    r"|study design|trial design|study population|participants?|patients?|subjects?|eligibility|enrol{1,2}ment"
    r"|results?( and discussion)?|findings"
    r"|statistical analys[ei]s|statistics|data analysis"
    r"|animals?|mice|cell culture|cell lines?|in vitro|in vivo|xenografts?"
    r"|immunohistochemistry|flow cytometry|western blot(ting| analysis)?|rna sequencing|qpcr|elisa"
    r"|synthesis|characterization|simulations?|model(ing|ling)?|computational (methods?|analysis)"
    r")$"
)

# Headings whose content describes OTHER people's work. Excluded from tier
# detection -- this is the precision guard that makes full text safe to read.
_CITED_PAT = re.compile(
    r"^(introduction|background|discussions?|conclusions?|concluding remarks"
    r"|perspectives?|future (directions?|perspectives?|work)|outlook"
    r"|related work|literature review|state of the art|prior work"
    r"|significance|implications)$"
)

# Boilerplate: carries no evidence signal and adds noise.
_DROP_PAT = re.compile(
    r"^(references?|bibliography|acknowledge?ments?|acknowledgements?"
    r"|funding( information| statement| sources?)?|grant support"
    r"|(competing|conflicts?( of)?) interests?.*|conflict of interest.*|disclosures?.*"
    r"|authors?.? contributions?|author information|contributions"
    r"|ethics.*|consent.*|institutional review board.*|informed consent.*"
    r"|availability of data.*|data availability.*|supplementary.*|supporting information.*"
    r"|abbreviations?|publisher.?s note|peer review.*|reporting summary"
    r"|generative ai statement|google scholar|declarations?|orcid|keywords?"
    r"|copyright|license|footnotes?|editor.?s note|correction|erratum"
    r")$"
)

# A Methods-like heading specifically (used by callers that need to know whether
# the article documents an experimental design at all).
_METHODS_PAT = re.compile(
    r"^(materials?( and | & )methods?|methods?|methodology"
    r"|experimental( section| procedures?| methods?| design)?"
    r"|patients? and methods?|subjects? and methods?|methods? and materials?"
    r"|study design|trial design)$"
)

_FULL_TEXT_MARKER = "## Full Text"
_STRIP = re.compile(r"[^a-z0-9&\s'’-]")
_NUMPREFIX = re.compile(r"^\s*\d+(\.\d+)*[.)]?\s+")
_MAX_HEADING_CHARS = 70
_MAX_HEADING_WORDS = 7


def normalize_heading(s: str) -> str:
    """Lowercase, drop a leading '2.1 ' style number, strip punctuation."""
    s = _NUMPREFIX.sub("", s.strip()).strip().lower()
    s = _STRIP.sub("", s).strip()
    return re.sub(r"\s+", " ", s)


def classify_heading(heading: str):
    """Return SELF/CITED/DROP for a recognised heading, else None."""
    h = normalize_heading(heading)
    if not h:
        return None
    if _DROP_PAT.match(h):
        return DROP
    if _CITED_PAT.match(h):
        return CITED
    if _SELF_PAT.match(h):
        return SELF
    return None


def _is_heading_line(lines, j) -> bool:
    """Blank-line-surrounded short line that does not read like a sentence."""
    s = lines[j].strip()
    if not s or len(s) > _MAX_HEADING_CHARS:
        return False
    if s[0] in "#|-*>![](":
        return False
    if s.endswith((".", ",", ";")):
        return False
    if not 1 <= len(s.split()) <= _MAX_HEADING_WORDS:
        return False
    if not (s[0].isupper() or s[0].isdigit()):
        return False
    prev = lines[j - 1].strip() if j > 0 else ""
    nxt = lines[j + 1].strip() if j + 1 < len(lines) else ""
    return prev == "" and nxt == ""


def split_sections(body: str):
    """Return [(heading, klass, text)] for the '## Full Text' portion of an article.

    Returns [] when the article has no full-text section (abstract-only records).
    """
    i = body.find(_FULL_TEXT_MARKER)
    if i < 0:
        return []
    lines = body[i + len(_FULL_TEXT_MARKER):].split("\n")
    out = []
    head, klass, buf = "", SELF, []
    for j, line in enumerate(lines):
        if _is_heading_line(lines, j):
            k = classify_heading(line)
            if k is not None:
                if buf:
                    out.append((head, klass, "\n".join(buf)))
                head, klass, buf = line.strip(), k, []
                continue
        buf.append(line)
    if buf:
        out.append((head, klass, "\n".join(buf)))
    return out


def self_text(body: str) -> str:
    """Full text restricted to SELF sections; '' when there is no full text."""
    return "\n".join(t for _, k, t in split_sections(body) if k == SELF)


def has_methods_section(body: str, min_chars: int = 400) -> bool:
    """True when the article carries a substantive Methods-like section.

    NOTE: measured on the v1 gold set this is only a weak signal for
    'is this primary research' (26% precision as a lone none-applicable
    predictor), so it must be combined with other evidence, never used alone.
    """
    for head, _klass, text in split_sections(body):
        if _METHODS_PAT.match(normalize_heading(head)) and len(text) >= min_chars:
            return True
    return False


def section_coverage(body: str) -> dict:
    """Diagnostic: characters per class."""
    out = {SELF: 0, CITED: 0, DROP: 0}
    for _head, klass, text in split_sections(body):
        out[klass] += len(text)
    return out
