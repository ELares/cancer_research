"""Every census figure the manuscript quotes, tied to the artifact that produced it.

Fourteen PRs of corrections have moved numbers in the manuscript, several of
them more than once: the sonodynamic trial share was 4.54% and is 0.29%, the
thesis-direction lead was 89.1% and is 74.6%. Each correction was propagated by
hand. This checks the result rather than trusting that.

DECLARED, NOT INFERRED. A sweep asserting that every number in the manuscript
appears in some artifact would flag years, section numbers, PMIDs, simulation
outputs and quoted historical figures -- and this repo has already learned from
`audit_heading_only_assertions.py` that a check at roughly half precision
trains exemption-adding rather than fixing. So each entry names a claim, the
artifact path that produces it, and how it is rendered in prose.

WHAT THIS CANNOT CATCH: a number that is correct but attached to the wrong
sentence. `census-mechanism-profile.json` can say 2,513 while the manuscript
attributes it to the wrong mechanism, and both would pass. This is a freshness
guard, not a reading.
"""
import json
import re
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
MANUSCRIPT = REPO / "article/drafts/v1.md"


def _get(artifact: str, path: str):
    """Fetch a value by dotted path, with [i] indexing for lists."""
    d = json.loads((REPO / "analysis" / artifact).read_text())
    for part in path.split("."):
        m = re.fullmatch(r"([^\[]*)\[(\d+)\]", part)
        if m:
            if m.group(1):
                d = d[m.group(1)]
            d = d[int(m.group(2))]
        else:
            d = d[part]
    return d


def _comma(v):
    return f"{int(v):,}"


def _pct(v):
    return f"{v}%"


def _pct_rounded(v):
    """Accept either the artifact's precision or a sensibly rounded rendering.

    Prose rounds, and it should: 11 of 18 articles is 61.1%, and quoting a
    decimal on eighteen observations is false precision. The guard should not
    force the manuscript to be worse-written than it ought to be -- but it must
    not accept an arbitrary nearby number either, so only the exact value and
    its integer rounding pass.
    """
    return [f"{v}%", f"{round(v)}%"]


# (label, artifact, dotted path, formatter)
CLAIMS = [
    ("census records", "census-evidence-design.json", "census", _comma),
    ("classifiable records", "census-evidence-design.json", "classifiable", _comma),
    ("trial publication types", "census-evidence-design.json", "classes.trial", _comma),
    ("full-text ceiling", "census-fulltext-ceiling.json", "ceiling_records", _comma),
    ("ceiling share", "census-fulltext-ceiling.json",
     "ceiling_share_of_undetermined", _pct),
    ("mechanism growth", "census-mechanism-growth.json", "union_growth", str),
    ("field growth", "census-mechanism-growth.json", "field_growth", str),
    ("chains matched", "census-diagnostic-chains.json", "census_matched", _comma),
    ("matrix universe", "census-mechanism-cancer-matrix.json", "universe", _comma),
    ("matrix zero cells", "census-mechanism-cancer-matrix.json", "n_zero", str),
    ("matrix cells", "census-mechanism-cancer-matrix.json", "n_cells", str),
    ("thesis exploit share, corrected", "census-thesis-direction.json",
     "adjudication.corrected_exploit_share", _pct),
    ("hypoxia protective share", "census-hypoxia-direction.json",
     "adj_protects_share", _pct_rounded),
]


@pytest.mark.parametrize("label,artifact,path,fmt", CLAIMS,
                         ids=[c[0] for c in CLAIMS])
def test_a_quoted_census_figure_matches_its_artifact(label, artifact, path, fmt):
    want = fmt(_get(artifact, path))
    wants = want if isinstance(want, list) else [want]
    txt = " ".join(MANUSCRIPT.read_text().split())
    assert any(w in txt for w in wants), (
        f"the manuscript does not carry {label} as any of {wants} from "
        f"{artifact}:{path}. Either the artifact was regenerated and the "
        f"manuscript not updated, or the claim was reworded and this entry "
        f"should follow it.")


def test_the_sonodynamic_correction_is_present_on_both_arms():
    """The figure this campaign corrected most consequentially.

    The descriptor arm's 4.54% may still appear -- it is real, and the
    manuscript explains what it is measuring -- but only alongside the
    corrected 0.29%, never alone.
    """
    txt = " ".join(MANUSCRIPT.read_text().split())
    d = json.loads((REPO / "analysis/census-modality-comparison.json").read_text())
    sono = next(r for r in d["rows"] if r["modality"] == "sonodynamic")
    assert f"{sono['text_trial_share']}%" in txt, (
        "the corrected sonodynamic trial share is missing")
    if f"{sono['mesh_trial_share']}%" in txt:
        assert "98% borrowed" in txt or "only 2 of those 114" in txt, (
            "the descriptor-arm figure appears without the correction that "
            "explains it is measuring ultrasound hyperthermia")


def test_no_superseded_figure_survives_anywhere():
    """Values this campaign replaced, banned outright.

    Each was quoted in the manuscript at some point and each was corrected;
    a reappearance means a propagation missed a site.
    """
    txt = " ".join(MANUSCRIPT.read_text().split())
    superseded = {
        "a ratio of 8.2 to 1": "thesis-direction lead before adjudication",
        "leads 89% to 11%": "thesis-direction lead before adjudication",
        "a factor of 1.6 below HIFU": "sonodynamic gap before the descriptor fix",
        "2,038 of 4,830": "evidence-tag coverage over the retrieved corpus",
        # NOT banned outright: the taxonomy-sensitivity finding is a real
        # historical result and Section 3.6 cites it as the reason the count
        # stopped being trusted. What must not survive is that figure standing
        # as the answer, so the ban is on it appearing WITHOUT the census
        # result that superseded it -- checked separately below.
    }
    found = [f"{k!r} ({why})" for k, why in superseded.items() if k in txt]
    # The taxonomy-sensitivity gap count may be cited as history, but never
    # without the census result that replaced it.
    if "94 (22.5% of cells)" in txt:
        assert "6 empty cells of 288" in txt, (
            "the corpus-scale gap count appears without the census result that "
            "superseded it, so it reads as the current answer")
    assert not found, (
        "the manuscript still carries superseded figures:\n  "
        + "\n  ".join(found))


def test_the_claim_list_is_declared_and_reads_the_artifacts():
    """A guard comparing prose to a hardcoded number goes stale with it."""
    assert len(CLAIMS) >= 10
    src = Path(__file__).read_text()
    assert "DECLARED, NOT INFERRED" in src
    for label, artifact, path, _ in CLAIMS:
        assert (REPO / "analysis" / artifact).exists(), (
            f"{label} names a missing artifact: {artifact}")
        _get(artifact, path)   # raises if the path has moved
