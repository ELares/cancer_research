"""No generator may hand-write a QUANTITY into prose it also derives.

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

GENERATORS = [
    REPO_ROOT / "scripts" / "atlas_combination_gaps.py",
    REPO_ROOT / "scripts" / "atlas_variant_drug_map.py",
]

# Nouns whose count is a measurement in these documents.
COUNTABLE = (r"row|pair|paper|rsid|regimen|spelling|sweep|panel|sample|entry|"
             r"twin|collapse|refusal|triple|variant|gene|drug|mutation|record|"
             r"article|combination|finding|guard|instance")
# THREE and above. `one`, `two`, `both` and `several` are grammatical rather
# than measured -- "pairs resting on ONE paper" defines a category and cannot go
# stale, while "a nine-pair panel" and "ten rsids" are counts of data that did.
# Every historical instance of this defect used three or more, or a vague
# magnitude.
SPELLED = (r"three|four|five|six|seven|eight|nine|ten|eleven|twelve|"
           r"dozen|tens|hundreds|thousands")

QUANTITY_SHAPES = [
    # 40-row sample, 36-pair sweep, 37 pairs, 25,443 rows
    # Digits three and above; 1 and 2 are grammatical for the same reason.
    (rf"\b(?!(?:1|2)\b)\d[\d,]*\s*[-–]?\s*(?:{COUNTABLE})s?\b",
     "a number modifying a countable noun"),
    # nine-pair panel, ten rsids, tens of papers
    (rf"\b(?:{SPELLED})\s*[-–]?\s*(?:{COUNTABLE})s?\b", "a spelled quantity modifying a countable noun"),
    # below 1.3, above 3.7, at least 30, more than 2
    (r"\b(?:below|above|at least|at most|more than|fewer than|under|over|"
     r"exceeds|beyond|up to|as many as)\s+\d", "a comparative pointing at a number"),
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
    r"Chr\d+\S*|[A-Z][A-Za-z]*\d[A-Za-z0-9]*|[A-Z]{2,}[- ]?\d+|"
    r"CodeBreaK|KRYSTAL-\d|NEJ\d+|AG\d+-\d+|CAPItello-\d+|SOLAR-\d|FLAURA\d)$")


def _prose_literals(path: Path):
    """Non-f-string literals outside docstrings, with line numbers.

    An f-string is a JoinedStr, so interpolating exempts a sentence
    automatically. Docstrings are excluded here and checked separately, under a
    rule they can actually satisfy: they cannot interpolate at all.
    """
    tree = ast.parse(path.read_text())
    docstrings = set()
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.FunctionDef, ast.AsyncFunctionDef,
                             ast.ClassDef)):
            body = getattr(node, "body", None)
            if body and isinstance(body[0], ast.Expr) and \
                    isinstance(body[0].value, ast.Constant) and \
                    isinstance(body[0].value.value, str):
                docstrings.add(id(body[0].value))
    return [(n.lineno, n.value) for n in ast.walk(tree)
            if isinstance(n, ast.Constant) and isinstance(n.value, str)
            and id(n) not in docstrings]


def _strip_identifiers(text: str) -> str:
    return " ".join("" if IDENT_TOKEN.match(t.strip("`*|,.()[]"))
                    else t for t in text.split())


def offending(path: Path):
    bad = []
    for lineno, text in _prose_literals(path):
        flat = " ".join(text.split())
        if len(flat.split()) <= 1:
            continue
        cleaned = _strip_identifiers(flat)
        for rx, why in COMPILED:
            m = rx.search(cleaned)
            if m:
                bad.append((lineno, why, m.group(0).strip(), flat[:110]))
                break
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
        frags = " | ".join(t for _, _, _, t in hits)
        for must in ("nine-pair", "below 1.3", "36-pair", "ten rsids",
                     "2.98/3.82", "tens or hundreds"):
            assert must in frags, f"the detector missed {must!r}"
        for must_not in ("derived:", "rs77375493", "EGFR L858R", "CodeBreaK"):
            assert must_not not in frags, (
                f"the detector flagged {must_not!r}, which is an interpolation "
                "or an identifier and is the sanctioned form")
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
        doc = ast.get_docstring(ast.parse(g.read_text())) or ""
        bad = []
        for line in doc.splitlines():
            flat = " ".join(line.split())
            if len(flat.split()) <= 1:
                continue
            cleaned = _strip_identifiers(flat)
            for rx, why in COMPILED:
                m = rx.search(cleaned)
                if m:
                    bad.append(f"[{why}: {m.group(0).strip()!r}] {flat[:100]}")
                    break
        if bad:
            offenders[g.name] = bad
    assert not offenders, (
        "these docstring lines carry a quantity a docstring cannot keep fresh:\n"
        + "\n".join(f"  {n}: {t}" for n, ts in offenders.items() for t in ts)
        + "\n\nState the property without the number and let the report carry "
          "the measurement.")
