"""The manuscript's module-support paragraph must not drift from its artifact.

Section 8.4 now states how much literature stands behind each simulation layer,
and every figure in it comes from `analysis/atlas-module-support.json`. The
manuscript spells numbers as words, so nothing about the paragraph would break
if the underlying measurement moved -- which is exactly the failure mode this
repo has spent a session removing from its generated documents. The manuscript
is hand-written and cannot be regenerated, so a guard is the only option.

The paragraph is also the one place a reader meets these numbers without the
surrounding caveats, so the caveats are pinned too: an eleven-of-twenty count
quoted without the exposure finding would tell a reader that nine layers model
something the literature rejects, which is the reading the measurement refutes.
"""

import json
import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
MS = REPO_ROOT / "article" / "drafts" / "v1.md"
RAW = REPO_ROOT / "analysis" / "atlas-module-support.json"

WORDS = {6: "six", 7: "seven", 8: "eight", 9: "nine", 11: "eleven",
         17: "seventeen", 20: "twenty"}


KEY = "How much literature stands behind these mechanisms"


def _para() -> str:
    """Both paragraphs of the passage: the counts, then the correction."""
    txt = MS.read_text()
    assert KEY in txt, "the module-support passage is missing from the manuscript"
    start = txt.index("**" + KEY)
    end = txt.index("`analysis/atlas-module-support.md`.", start)
    return txt[start:end]


def _body() -> str:
    """The passage WITHOUT its bold heading.

    The heading contains the phrase "is wrong", so a guard searching the whole
    passage for it passes while the body's actual correction is deleted --
    proved by mutation. Claims must be checked where they are argued, not where
    they are announced.

    The first version of this helper was a NO-OP: `_para()` started at the key,
    which is AFTER the opening `**`, so `startswith("**")` was always false and
    the heading was never stripped. It passed only because a stricter literal
    elsewhere happened to catch the case, and its own docstring said it was
    safe. `_para()` now anchors on the `**` so the split has something to cut.
    """
    para = _para()
    return para.split("**", 2)[-1]


def test_the_counts_match_the_artifact():
    d, para = json.loads(RAW.read_text()), _para()
    assert f"{WORDS[d['corroborated']]} are corroborated" in para.lower(), (
        f"the manuscript's corroborated count is not {d['corroborated']}")
    assert f"{WORDS[d['n_claims']]} such pairs" in para.lower(), (
        f"the manuscript's denominator is not {d['n_claims']}")
    # The layers-only figure is the one a reader takes away about the LAYERS,
    # and it is the one the first draft omitted -- quoting 11/20 as though it
    # described them, when three of the twenty are core-engine mechanisms and
    # three of the four best-corroborated rows.
    assert (f"{WORDS[d['corroborated_layers_only']]} of "
            f"{WORDS[d['n_claims_layers_only']]}") in para.lower(), (
        f"the manuscript does not give the layers-only count "
        f"({d['corroborated_layers_only']} of {d['n_claims_layers_only']})")
    assert (f"{WORDS[d['zero_explained_by_exposure']]} of the "
            f"{WORDS[d['zero_relation']]} zeros") in para.lower(), (
        "the manuscript's below-the-line count is not the artifact's")


def test_the_correlation_matches_the_artifact():
    d, para = json.loads(RAW.read_text()), _para()
    m = re.search(r"rank correlation ([\d.]+)", para)
    assert m, "the manuscript states no correlation"
    assert abs(float(m.group(1)) - d["spearman_weaker_degree_vs_relations"]) < 0.005, (
        f"manuscript says {m.group(1)}, artifact says "
        f"{d['spearman_weaker_degree_vs_relations']:.4f}")


def test_the_spread_and_base_rate_are_the_artifact_s():
    """Three figures the first guard left unpinned, proved by mutation.

    'six to nearly three thousand', the GPX4 robustness claim, and the base
    rate all survived every original assertion -- the base rate was not even a
    key in the JSON, so a graph rebuild could move it and leave the manuscript
    silently wrong.
    """
    d, para = json.loads(RAW.read_text()), _para()
    assert f"{WORDS[d['weaker_min']]} to nearly" in para.lower(), (
        f"the stated lower bound is not the artifact's {d['weaker_min']}")
    assert round(d["weaker_max"], -3) == 3000, (
        f"'nearly three thousand' no longer describes {d['weaker_max']}")
    m = re.search(r"rising to ([\d.]+) when every pair containing GPX4", para)
    assert m, "the GPX4 robustness claim is gone or reworded"
    assert abs(float(m.group(1)) - d["spearman_excluding_gpx4"]) < 0.005, (
        f"manuscript says {m.group(1)}, artifact says "
        f"{d['spearman_excluding_gpx4']:.4f}")
    assert float(m.group(1)) > d["spearman_weaker_degree_vs_relations"], (
        "the manuscript claims the association RISES without GPX4; it does not")
    assert 0.44 <= d["below_floor_base_rate"] < 0.46, (
        f"the base rate moved to {d['below_floor_base_rate']:.3f}; "
        "'forty-five per cent' is no longer accurate")


def test_the_refutation_is_scoped_to_what_was_measured():
    """The correction holds for seven rows, not for all nine.

    The first draft said the obvious reading of 'the remaining nine' is wrong.
    For two of them nothing is established either way -- the paragraph's own
    caveat concedes it -- so claiming refutation for nine overstated by two.
    """
    body = _body().lower()
    assert "wrong for most of them" in body, (
        "the passage claims the reading is refuted for every uncorroborated "
        "row; it was measured for seven of nine")
    assert "nothing is established in either direction" in body, (
        "the two unresolved rows are quietly counted with the seven")


def test_the_exposure_floor_is_described_consistently():
    """'at least a hundred partners' must remain true of the measured floor."""
    d, para = json.loads(RAW.read_text()), _para()
    assert "at least a hundred partners" in para, (
        "the manuscript no longer describes the exposure floor")
    assert 100 <= d["exposure_floor"] < 200, (
        f"the floor moved to {d['exposure_floor']}; 'at least a hundred' is no "
        "longer an accurate rounding and the sentence must be rewritten")


def test_the_reading_that_the_measurement_refutes_is_pre_empted():
    """The count must never appear without the correction.

    Eleven of twenty, quoted alone, tells a reader that nine layers model
    something the literature rejects. That is the reading the exposure result
    refutes, and it is the reason the paragraph exists.
    """
    body = _body().lower()
    assert "is wrong for most of them, and measurably so" in body, (
        "the passage announces a correction in its heading but the body no "
        "longer makes one, so the count stands uncorrected where it is argued")
    assert "what has been studied, not about what is true" in body


def test_both_limits_survive():
    """The two things that stop this going further must stay stated."""
    para = _para()
    assert "minimum of an eleven-point sample" in para, (
        "the floor is presented as a detection limit")
    # DERIVED, not spelled. This assertion held the literal "forty-five per
    # cent" while the figure it guards lives in the artifact, so a graph
    # rebuild moved the base rate to forty-six and the guard failed for the
    # right reason but with the wrong remedy on offer -- edit the guard, or
    # edit the manuscript? Deriving it removes the question.
    d = json.loads(RAW.read_text())
    pct = round(100 * d["below_floor_base_rate"])
    tens = {40: "forty", 50: "fifty", 30: "thirty", 60: "sixty"}
    word = (tens[pct // 10 * 10] if pct % 10 == 0
            else f"{tens[pct // 10 * 10]}-{WORDS[pct % 10]}")
    assert f"{word} per cent" in para, (
        f"the manuscript's base rate is not the artifact's {pct}% "
        f"(expected the words '{word} per cent')")
    assert "inverts the ranking" in para, (
        "the manuscript names or implies specific genuine zeros without the "
        "cross-measure disclosure")
