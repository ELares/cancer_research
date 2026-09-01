"""Two generators may not hand-write a QUANTITY into prose they also derive.

SCOPE, FIRST, BECAUSE THE TITLE WOULD OTHERWISE OVERSTATE IT
------------------------------------------------------------
This rule is ENFORCED on two files. It is not repo-wide.
`test_the_rule_is_honest_about_how_little_it_covers` measures and prints how
many other scripts of the same shape would fail it, so the gap is a reported
number rather than an implication. Those are candidates, not confirmed defects:
a hand-written figure that is correct today is still a figure that cannot stay
correct, which is the property being enforced here.

WHY THIS FILE EXISTS
--------------------
This is the defect this repository produces faster than it fixes. Three review
rounds on two atlas generators each found new instances IN THE FIXES FOR THE
PREVIOUS ROUND:

  * "a 40-row sample judged from title and abstract" -- never run here
  * "an independent 36-pair sweep put the boundary at 2.98/3.82" -- never run
    here, cited in three places, justifying a threshold on 25,443 rows
  * "a heuristic on a nine-pair panel" -- eight lines below a 37-row table
  * "every sequential pair below 1.3" -- the measured value is 2.98
  * "ten rsids are refused under one gene" -- the true count is zero
  * "tens or hundreds vs nought to a handful" -- three of six rows are 6, 4, 1

Patching them one at a time has not worked, because each fix is written in the
same voice that produced the defect. So the rule is mechanical: A GENERATOR MAY
NOT HAND-WRITE A QUANTITY. If a sentence needs a figure it must interpolate one,
which makes it impossible for the sentence to outlive the measurement.

WHAT COUNTS AS A QUANTITY, AND WHY THE TEST IS SHAPED THIS WAY
---------------------------------------------------------------
Not "any digit". A digit is usually an IDENTIFIER here -- `rs77375493`,
`p.V600E`, `SOLAR-1`, `MAP2K1`, `CodeBreaK 300` -- and banning those would ban
the worked examples the reports are built on. What goes stale is a number
MODIFYING A COUNTABLE NOUN, a comparative pointing at a number, or a bare range.
Those three shapes cover every instance listed above and none of the
identifiers.

The detector validates itself against planted samples, because a scan returning
zero because it is broken looks exactly like a clean file --
tests/test_no_vacuous_assertions.py records the same requirement.
"""

import ast
import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent

# ENFORCED for these. The rule is not repo-wide and this file must not imply
# that it is: `all_generators()` below measures how many other scripts of the
# same shape would fail it, and a test reports that number so the scope of the
# exemption is visible rather than assumed.
GENERATORS = [
    REPO_ROOT / "scripts" / "atlas_combination_gaps.py",
    REPO_ROOT / "scripts" / "atlas_variant_drug_map.py",
    REPO_ROOT / "scripts" / "manuscript_vs_census.py",
    # enrolled with the analysis it supports: a page whose subject is a
    # sentence outliving the measurement beside it has no business
    # hand-writing its own figures
    REPO_ROOT / "scripts" / "atlas_descriptor_recall.py",
    # enrolled after it shipped "four arms reach through `immune.rs`" beside
    # its own table row naming six consumers -- a hand-written count of a
    # countable noun, in the one page whose entire subject is a hand-written
    # count outliving the thing it counted
    REPO_ROOT / "scripts" / "modality_module_depth.py",
]


def all_generators():
    """Every script that renders prose into an analysis/*.md it also computes."""
    out = []
    for p in sorted((REPO_ROOT / "scripts").glob("*.py")):
        t = p.read_text(errors="ignore")
        if 'analysis"' not in t or ".md" not in t:
            continue
        if "write_text" not in t and "OUT_MD" not in t:
            continue
        try:
            ast.parse(t)
        except SyntaxError:
            continue        # test_every_script_at_least_parses reports these
        out.append(p)
    return out

# NOT a list of measured nouns. A closed allowlist is always one synonym from
# being stepped around -- the first version omitted `key`, the very noun the fix
# it shipped alongside introduced, so "ten refused keys" passed. The rule is
# inverted: a number adjacent to ANY ordinary lowercase word is a quantity, and
# a short EXEMPT list carries the collocations that are structural rather than
# measured. That fails safe: an unlisted noun is a false positive an author must
# consciously exempt, not a silent pass.
COUNTABLE = r"[a-z][a-z-]{2,}"
# Words that follow a number without counting anything.
# Words that are never the noun a number is counting, so the retry must not
# land on one. Without this the skip-past-exempt walk reached "this", "them"
# and "where" and reported them as counted nouns.
NEVER_THE_NOUN = {
    "this", "that", "these", "those", "them", "they", "it", "its", "their",
    "there", "here", "where", "when", "which", "who", "what", "such",
    "different", "other", "same", "own", "each", "every", "any", "some",
    "both", "either", "neither", "more", "most", "less", "least", "very",
    "only", "just", "still", "also", "then", "so", "because", "while",
    "however", "rather", "instead", "again", "already", "never", "always",
}
EXEMPT_AFTER = {
    "and", "or", "to", "of", "in", "at", "on", "is", "was", "were", "per",
    "than", "as", "the", "a", "an", "for", "with", "by", "from", "that",
    "which", "relates", "sits", "inside", "above", "below", "gives", "holds",
    "reads", "means", "denotes", "carries", "would", "will", "can", "may",
    # Structural rather than measured: these count a rule, a design, or a named
    # variant's alleles, none of which the corpus can change under the prose.
    # Each is an EXEMPTION, so an unlisted noun still fails and the author must
    # add it deliberately -- unlike an allowlist of measured nouns, which any
    # synonym walks past.
    "tests", "test", "things", "substitutions", "site", "refuses", "states",
    "classes", "shapes", "reasons", "arms", "halves", "outcomes",
    "branches", "steps", "passes", "rules", "checks",
    # NOT exempt: "tests" and "preceding". Those count the generator's OWN
    # design, and the justification for this list -- things the corpus cannot
    # change under the prose -- is true of the corpus and false of the code. A
    # fourth agreement test would stale "the three tests apply" silently.
}
# THREE and above. `one`, `two`, `both` and `several` are grammatical rather
# than measured -- "pairs resting on ONE paper" defines a category and cannot go
# stale, while "a nine-pair panel" and "ten rsids" are counts of data that did.
# Every historical instance of this defect used three or more, or a vague
# magnitude.
# Built from morphemes rather than typed out, because a list that stops at
# twelve can be outrun by counting higher: "a forty-pair panel" evaded the first
# version, and this file's own prose writes "thirty-seven".
_UNITS = ("three four five six seven eight nine ten eleven twelve thirteen "
          "fourteen fifteen sixteen seventeen eighteen nineteen").split()
_TENS = ("twenty thirty forty fifty sixty seventy eighty ninety").split()
_BIG = ("hundred thousand million dozen tens hundreds thousands dozens").split()
SPELLED = "|".join(
    sorted(_UNITS + _TENS + _BIG + [f"{t}-{u}" for t in _TENS for u in _UNITS],
           key=len, reverse=True))

QUANTITY_SHAPES = [
    # 40-row sample, 36-pair sweep, 37 pairs, 25,443 rows
    # Digits three and above; 1 and 2 are grammatical for the same reason.
    (rf"\b(?!(?:1|2)\b)\d[\d,]*(?:\.\d+)?\s*[-–]?\s*(?:{COUNTABLE})\b",
     "a number modifying a noun"),
    # A percentage is this repository's commonest way of stating a measurement
    # and no rule mentioned `%`, so "the relation holds in 96.0% of the cases"
    # was invisible.
    # No trailing \b: `%` is not a word character, so a boundary after it can
    # never match against a space, and the rule silently caught only the
    # spelled-out "percent" forms. The planted control found this.
    (r"\b\d+(?:\.\d+)?\s*(?:%|percent\b|percentage points?\b|pp\b)",
     "a percentage"),
    # A bare decimal threshold. `\d[\d,]*` never spanned the point, so the rule
    # tested the FRACTIONAL digits against the next word and caught "2.98
    # exactly" while missing "2.98 and the flag" -- it was not reading the
    # number it flagged.
    (r"\b\d+\.\d+\b", "a bare decimal threshold"),
    # nine-pair panel, ten rsids, tens of papers
    (rf"\b(?:{SPELLED})\s*[-–]?\s*(?:{COUNTABLE})\b",
     "a spelled quantity modifying a noun"),
    # below 1.3, above 3.7, at least 30, more than 2
    (r"\b(?:below|above|at least|at most|more than|fewer than|under|over|"
     r"exceeds|beyond|up to|as many as)\s+\d[\d,]*(?:\.\d+)?",
     "a comparative pointing at a number"),
    # 2.98/3.82, 7/7 vs 1/7, 8 of 10
    (r"\b\d+(?:\.\d+)?\s*(?:/|vs\.?|versus|against|to|and|of)\s*\d+(?:\.\d+)?\b",
     "a bare numeric range or ratio"),
    # tens or hundreds ... nought to a handful
    (r"\b(?:tens|hundreds|thousands|dozens)\b[^.]{0,70}"
     r"\b(?:handful|nought|none|zero|few)\b", "a vague magnitude comparison"),
]
COMPILED = [(re.compile(p, re.I), why) for p, why in QUANTITY_SHAPES]

# Identifier shapes that would otherwise trip the noun rule. Each must match a
# WHOLE token, never a substring, so `p.V600E` is exempt and "V600E in 181
# papers" is not.
IDENT_TOKEN = re.compile(
    r"^(?:rs\d+|MESH:[A-Z]?\d+|p\.[A-Z]\d+[A-Z*]|[cgmn]\.[+-]?\d+[ACGT]>[ACGT]|"
    r"Chr\d+\S*|[A-Z][A-Za-z]*\d[A-Za-z0-9]*|[A-Z]{2,}[- ]?\d+|codon-?\d+|"
    r"CodeBreaK|KRYSTAL-\d|NEJ\d+|AG\d+-\d+|CAPItello-\d+|SOLAR-\d|FLAURA\d)$")


def _prose_runs(path: Path):
    """RUNS of adjacent string literals in a list, joined, with line numbers.

    Scanning one literal at a time made coverage a coin flip. Both generators
    wrap prose across adjacent elements of an `L` list every six to ten words at
    an arbitrary column, so a number and its noun share a literal only by luck:
    the retracted "in tens or hundreds of papers ... in nought to a handful"
    passes a per-literal scan purely by wrapping after "papers,".

    A run ends at an f-string, because an interpolated figure is the sanctioned
    form and a sentence containing one is derived.
    """
    tree = ast.parse(path.read_text())
    docstrings = set()
    for node in ast.walk(tree):
        body = getattr(node, "body", None)
        if isinstance(body, list) and body and isinstance(body[0], ast.Expr) \
                and isinstance(body[0].value, ast.Constant) \
                and isinstance(body[0].value.value, str):
            docstrings.add(id(body[0].value))

    runs, seen = [], set()
    for node in ast.walk(tree):
        # Lists only. Prose is built with `L = [...]`; the panels are LISTS OF
        # TUPLES whose inner tuples hold identifiers, and joining those produced
        # "673 COMBI-d" -- a gene id beside a trial name, which is data, not a
        # sentence. Tuple elements are marked seen so the loose-literal pass
        # below does not pick them up either.
        if isinstance(node, ast.Tuple):
            # Suppress JOINING only. Marking tuple strings `seen` skipped them
            # entirely, and six live prose glosses ride inside the collision
            # tuples -- so the exclusion written for `673 COMBI-d` was hiding
            # shipped sentences. They fall through to the loose-literal pass.
            continue
        if not isinstance(node, ast.List):
            continue
        cur = []
        for el in list(node.elts) + [None]:
            if isinstance(el, ast.Constant) and isinstance(el.value, str) \
                    and id(el) not in docstrings:
                cur.append(el)
                seen.add(id(el))
                continue
            if cur:
                runs.append((cur[0].lineno, " ".join(c.value for c in cur)))
            cur = []
    # literals outside any list still count on their own
    for n in ast.walk(tree):
        if isinstance(n, ast.Constant) and isinstance(n.value, str) \
                and id(n) not in docstrings and id(n) not in seen:
            runs.append((n.lineno, n.value))
    return runs


# A section reference names a place in a document; it is not a measurement and
# cannot go stale against the corpus. Stripped as a PHRASE before tokenising,
# because `8.2` on its own is indistinguishable from a threshold.
SECTION_REF = re.compile(r"\b(?:section|§|chapter)\s*\d+(?:\.\d+)*", re.I)


def _strip_identifiers(text: str) -> str:
    text = SECTION_REF.sub(" ", text)
    return " ".join("" if IDENT_TOKEN.match(t.strip("`*|,.()[]"))
                    else t for t in text.split())


def _matches(text: str):
    """Every quantity shape in a piece of prose, exemptions applied.

    ONE definition, used by the prose scan AND the docstring scan. The two had
    separate loops and the docstring one skipped the exemption filter, so
    structural phrases the prose scan correctly ignored fired there instead --
    which is the same two-call-sites-disagree defect this guard exists to stop
    the generators committing.
    """
    cleaned = _strip_identifiers(text)
    out = []
    for rx, why in COMPILED:
        for m in rx.finditer(cleaned):
            frag = m.group(0).strip()
            if "modifying a noun" in why:
                tail = re.split(r"[\s-]+", frag)[-1].lower()
                if tail in EXEMPT_AFTER:
                    # Do NOT drop the sentence. `finditer` resumes past the
                    # match, so an exempt first word ("40 of the regimens")
                    # masked the real noun entirely. Re-test from just after
                    # the number, skipping any run of exempt words.
                    words = [t for t in cleaned[m.start():].split() if t]
                    j, skipped = 1, 0
                    # At most two exempt words. Walking further reached
                    # pronouns several clauses away and called them counted.
                    while j < len(words) and skipped < 2 and \
                            words[j].strip(",.;:").lower() in EXEMPT_AFTER:
                        j += 1
                        skipped += 1
                    if 0 < skipped and j < len(words):
                        nxt = words[j].strip(",.;:()`*")
                        low = nxt.lower()
                        if re.fullmatch(COUNTABLE, nxt) and \
                                low not in EXEMPT_AFTER and \
                                low not in NEVER_THE_NOUN:
                            out.append((why, f"{words[0]} ... {nxt}"))
                    continue
            out.append((why, frag))
    return out


def offending(path: Path):
    bad = []
    for lineno, text in _prose_runs(path):
        flat = " ".join(text.split())
        if len(flat.split()) <= 1:
            continue
        # EVERY rule that matches, not the first: joining a run means one run
        # can hold several defects, and stopping at the first hid the rest.
        for why, frag in _matches(flat):
            bad.append((lineno, why, frag, flat[:110]))
    return bad


def test_no_generator_hand_writes_a_quantity_into_its_prose():
    """The rule, in one assertion.

    A failure here is not a style complaint. Every literal it reports is a
    sentence that can outlive the number it describes, which is how this
    repository shipped a boundary of 1.3 beside a measured 2.98 and a
    "nine-pair panel" beside a table of thirty-seven rows.
    """
    offenders = {g.name: offending(g) for g in GENERATORS}
    offenders = {k: v for k, v in offenders.items() if v}
    assert not offenders, (
        "these generator literals hand-write a quantity, so the sentence can "
        "outlive the measurement:\n"
        + "\n".join(f"  {name}:{ln}  [{why}: {frag!r}]\n      {txt}"
                    for name, hits in offenders.items()
                    for ln, why, frag, txt in hits)
        + "\n\nInterpolate the figure from the result dict. If the flagged "
          "token is genuinely an identifier, extend IDENT_TOKEN with a pattern "
          "that matches the WHOLE token.")


def test_the_detector_finds_every_shape_it_claims_to():
    """A scan returning zero because it is broken looks like a clean file."""
    import tempfile
    planted = (
        'def render(r):\n'
        '    return ["a heuristic on a nine-pair panel, not a classifier",\n'
        '            "every sequential pair scores below 1.3 with no overlap",\n'
        '            "an independent 36-pair sweep settled it",\n'
        '            "ten rsids are refused under one gene",\n'
        '            "the boundary is 2.98/3.82 on this panel",\n'
        '            "tens or hundreds against nought to a handful",\n'
        '            f"derived: {r} rows and {r} pairs",\n'
        '            "| rs77375493 | p.V600E | MAP2K1 | SOLAR-1 |",\n'
        '            "Pemetrexed is tied to EGFR L858R and KRAS G12C",\n'
        '            "sotorasib + panitumumab (CodeBreaK 300)"]\n')
    with tempfile.NamedTemporaryFile("w", suffix=".py", delete=False) as fh:
        fh.write(planted)
        tmp = Path(fh.name)
    try:
        hits = offending(tmp)
        # The MATCHED fragments, not the truncated run text. Joining runs made
        # the whole planted list one run, so the reported text is clipped and a
        # check against it silently stopped testing most of the shapes.
        frags = " | ".join(f for _, _, f, _ in hits)
        runtext = " | ".join(t for _, _, _, t in hits)
        for must in ("nine-pair", "below 1.3", "36-pair", "ten rsids",
                     "2.98/3.82", "tens or hundreds"):
            assert must in frags, f"the detector missed {must!r}"
        for must_not in ("derived:", "rs77375493", "EGFR L858R", "CodeBreaK"):
            assert must_not not in frags and must_not not in runtext, (
                f"the detector flagged {must_not!r}, which is an interpolation "
                "or an identifier and is the sanctioned form")
    finally:
        tmp.unlink()


def test_every_script_at_least_parses():
    """Found by accident, and nothing in the suite would have found it on purpose.

    `scripts/evaluate_evidence_v2.py` carried a string literal broken across two
    lines with no continuation and had not parsed since the commit that
    introduced it. The file cannot be imported or run, and the analysis it
    generates cannot be regenerated, but no test imports it and no CI job
    compiles it, so it sat green.

    This is the cheapest possible guard and it belongs beside the others here:
    a generator that cannot parse cannot keep any of its figures fresh either.
    """
    bad = []
    for f in sorted((REPO_ROOT / "scripts").rglob("*.py")):
        try:
            ast.parse(f.read_text(errors="ignore"))
        except SyntaxError as e:
            bad.append(f"{f.relative_to(REPO_ROOT)}:{e.lineno}  {e.msg}")
    assert not bad, (
        "these scripts do not parse, so they cannot be run or regenerated:\n  "
        + "\n  ".join(bad))


def test_the_rule_is_honest_about_how_little_it_covers():
    """ADVISORY IN EFFECT, BLOCKING ON THE HONESTY OF THE CLAIM.

    This file enforces the rule on two generators. The same detector, run over
    every script that renders prose into an analysis/*.md it also computes,
    finds the shape widely. That is the real size of the class, and printing it
    here keeps the file from reading as a repo-wide guarantee it does not
    provide.

    It deliberately does NOT assert the uncovered count is zero. Every literal
    it counts is a candidate, not a confirmed defect: many may be correct today
    and simply unable to stay correct. Auditing them is a separate piece of
    work, and failing CI on all of them would force that work to happen inside
    an unrelated change.
    """
    gens = all_generators()
    covered = {g.resolve() for g in GENERATORS}
    uncovered = [(g.name, len(offending(g))) for g in gens
                 if g.resolve() not in covered]
    dirty = [(n, h) for n, h in uncovered if h]
    # "literals" was the wrong unit and not comparable across commits: a run is
    # a GROUP of adjacent literals, and every matching rule is now reported, so
    # the total counts (run, rule) matches. Naming it wrongly here, in the file
    # that enforces derived units, is the defect this file is about.
    print(f"\n  rule enforced on {len(covered)} of {len(gens)} generators; "
          f"{len(dirty)} of the remaining {len(uncovered)} carry a "
          f"hand-written quantity ({sum(h for _, h in dirty)} "
          f"(run, rule) matches, not distinct literals)")
    for n, h in sorted(dirty, key=lambda r: -r[1])[:10]:
        print(f"    {h:>3}  {n}")
    assert len(covered) == len(GENERATORS), "GENERATORS holds a duplicate"
    for g in GENERATORS:
        assert g.exists(), f"{g.name} is enforced but does not exist"
        assert g.resolve() in {x.resolve() for x in gens}, (
            f"{g.name} is enforced but is not a generator by the same "
            "definition this test uses, so the coverage figure it prints is "
            "measuring a different population")


def test_the_docstring_scan_would_notice_if_it_stopped_scanning():
    """It had no positive control, so a broken scan looked like clean docstrings.

    That is exactly the failure mode the detector self-test exists to prevent,
    and the commit that rewrote this scan from module-only to every-docstring
    did not extend the planted sample to cover it. Verified: mutating the
    collection to gather nothing left the suite green.
    """
    import tempfile
    planted = (
        '"""A module docstring holding a forty-pair panel."""\n'
        'def f():\n'
        '    """A function docstring holding ten rsids and 96.0% of cases."""\n'
        '    return 1\n')
    with tempfile.NamedTemporaryFile("w", suffix=".py", delete=False) as fh:
        fh.write(planted)
        tmp = Path(fh.name)
    try:
        tree = ast.parse(tmp.read_text())
        docs = [(getattr(n, "name", "<module>"), ast.get_docstring(n))
                for n in ast.walk(tree)
                if isinstance(n, (ast.Module, ast.FunctionDef, ast.ClassDef))
                and ast.get_docstring(n)]
        assert len(docs) == 2, (
            f"the docstring collection found {len(docs)} of 2 planted "
            "docstrings, so it is not walking function definitions")
        hits = [f for _, d in docs for _, f in _matches(" ".join(d.split()))]
        for must in ("forty-pair", "ten rsids", "96.0%"):
            assert any(must in h for h in hits), (
                f"the docstring scan missed {must!r} in a planted docstring")
    finally:
        tmp.unlink()


def test_the_module_docstrings_carry_no_quantity_either():
    """A docstring cannot interpolate, so it may carry no measurement at all.

    This is where the retractions did not reach: the 40-row sample was removed
    from the report and left standing in the docstring that generates the
    retraction, and the nine-pair boundary survived the same way.
    """
    offenders = {}
    for g in GENERATORS:
        tree = ast.parse(g.read_text())
        docs = []
        for node in ast.walk(tree):
            d = ast.get_docstring(node) if isinstance(
                node, (ast.Module, ast.FunctionDef, ast.AsyncFunctionDef,
                       ast.ClassDef)) else None
            if d:
                docs.append((getattr(node, "name", "<module>"), d))
        bad = []
        # Every docstring, not only the module's. A function docstring cannot
        # interpolate either, and three were shipping unscanned -- one of them
        # added by the commit that introduced this guard, in the docstring of
        # the helper written to replace the very phrase it repeated.
        for owner, doc in docs:
            # Join the docstring's lines too: a wrapped sentence in a docstring
            # hides a quantity exactly as one in a list does.
            flat_doc = " ".join(doc.split())
            for why, frag in _matches(flat_doc):
                i = flat_doc.find(frag.split()[0])
                bad.append(f"{owner}(): [{why}: {frag!r}] "
                           f"...{flat_doc[max(0, i - 40):i + 70]}...")
        if bad:
            offenders[g.name] = bad
    assert not offenders, (
        "these docstring lines carry a quantity a docstring cannot keep fresh:\n"
        + "\n".join(f"  {n}: {t}" for n, ts in offenders.items() for t in ts)
        + "\n\nState the property without the number and let the report carry "
          "the measurement.")
