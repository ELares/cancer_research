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

WORDS = {6: "six", 7: "seven", 9: "nine", 11: "eleven", 20: "twenty"}


def _para() -> str:
    txt = MS.read_text()
    key = "How much literature stands behind each layer"
    assert key in txt, "the module-support paragraph is missing from the manuscript"
    start = txt.index(key)
    return txt[start:txt.index("\n\n", start)]


def _body() -> str:
    """The paragraph WITHOUT its bold heading.

    The heading already contains the phrase "is wrong", so a guard searching the
    whole paragraph for it passes while the body's actual correction is deleted
    -- proved by mutation. Claims must be checked where they are argued, not
    where they are announced.
    """
    para = _para()
    return para.split("**", 2)[-1] if para.startswith("**") else para


def test_the_counts_match_the_artifact():
    d, para = json.loads(RAW.read_text()), _para()
    assert f"{WORDS[d['corroborated']]} of {WORDS[d['n_claims']]}" in para.lower(), (
        f"the manuscript's corroborated count is not {d['corroborated']} of "
        f"{d['n_claims']}")
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
    assert "is wrong, and measurably so" in body, (
        "the paragraph announces a correction in its heading but the body no "
        "longer makes one, so the count stands uncorrected where it is argued")
    assert "what has been studied, not about what is true" in body


def test_both_limits_survive():
    """The two things that stop this going further must stay stated."""
    para = _para()
    assert "minimum of an eleven-point sample" in para, (
        "the floor is presented as a detection limit")
    assert "forty-five per cent" in para, "the base rate that bounds it is gone"
    assert "inverts the ranking" in para, (
        "the manuscript names or implies specific genuine zeros without the "
        "cross-measure disclosure")
