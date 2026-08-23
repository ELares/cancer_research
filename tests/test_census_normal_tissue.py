"""Guards for the measurement that bears on this project's selectivity premise.

The premise -- ferroptosis inducers are attractive because normal cells resist
them -- is the most load-bearing assumption the simulation half rests on, and
it is the one an analysis is most tempted to protect. So the guards run in the
direction that protects the reader from the author: they pin the corrections
this analysis made against itself, and they refuse the two overclaims available
in each direction.

OVERCLAIM ONE, toward the finding: reading 575 organ-toxicity articles as a
refutation of the premise. It is not one. Cisplatin injures kidneys by several
mechanisms and ferroptosis being among them says nothing directly about a
targeted GPX4 inhibitor at a therapeutic dose.

OVERCLAIM TWO, toward the project: quietly dropping where those articles sit.
558 of 575 enter through the adjacent extension rather than the cancer tree,
which is mostly a labelling effect -- and saying so is what stops the rate
comparison being read as a charge of neglect against cancer researchers.

OFFLINE: reads only committed artifacts.
"""
import csv
import json
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
JSON = REPO / "analysis/census-normal-tissue.json"
MD = REPO / "analysis/census-normal-tissue.md"
ADJ = REPO / "analysis/normal-tissue-adjudication.csv"
CTRP = REPO / "analysis/calibration/ctrpv2_ferroptosis_curves.csv"


@pytest.fixture(scope="module")
def d():
    return json.loads(JSON.read_text())


def test_the_toxicity_descriptors_exclude_the_substring_false_friends(d):
    """`Cytotoxicity, Immunologic` and `T-Lymphocytes, Cytotoxic` contain
    "toxic" and mean the opposite -- the ability to kill a target cell, which
    is the INTENDED effect. A regex on "toxic" sweeps them in."""
    named = set(d["by_descriptor"])
    for bad in ("cytotoxicity, immunologic", "t-lymphocytes, cytotoxic",
                "natural cytotoxicity triggering receptor 1"):
        assert bad not in named, (
            f"{bad!r} is counted as organ toxicity; it is not")
    # The set itself must be explicit literals, not a pattern: a pattern is
    # what produced the false friends in the first place.
    import ast
    tree = ast.parse((REPO / "scripts/census_normal_tissue.py").read_text())
    node = next(n for n in tree.body
                if isinstance(n, ast.Assign)
                and getattr(n.targets[0], "id", None) == "ORGAN_TOX")
    assert isinstance(node.value, ast.Set), "ORGAN_TOX is not a literal set"
    literals = [e.value for e in node.value.elts if isinstance(e, ast.Constant)]
    assert len(literals) == len(node.value.elts), (
        "ORGAN_TOX contains a computed element, so what it matches cannot be "
        "read off the source")
    assert named <= set(literals), "a counted descriptor is not in ORGAN_TOX"


def test_direction_comes_from_the_adjudication_not_the_keyword_arm(d):
    """The inhibitor count is presence; the verdicts are the measurement."""
    rows = list(csv.DictReader(ADJ.open()))
    assert len(rows) == d["adjudicated_n"] == d["sample_n"]
    from collections import Counter
    assert d["verdicts"] == dict(Counter(r["adjudicated"] for r in rows).most_common())
    assert d["harm_share"] == pytest.approx(
        100 * d["verdicts"].get("harm", 0) / len(rows), abs=0.05)
    md = " ".join(MD.read_text().split())
    assert "cannot be read as papers proposing to block ferroptosis" in md, (
        "the report presents the inhibitor count without the note that the "
        "same string serves the opposite role as a specificity control")


def test_the_predicted_confound_is_reported_as_absent_and_misdirected(d):
    """The prediction was written down first and its direction was wrong.

    A file that silently drops a prediction it got backwards is worth less
    than one that keeps it.
    """
    rows = list(csv.DictReader(ADJ.open()))
    named = sum(1 for r in rows if r["regex_inhibitor"] == "True")
    probe = sum(1 for r in rows
                if r["regex_inhibitor"] == "True" and r["adjudicated"] == "probe")
    assert d["inhibitor_named_in_sample"] == named
    assert d["inhibitor_is_probe"] == probe
    md = " ".join(MD.read_text().split())
    if probe == 0:
        assert "its direction was predicted wrong" in md, (
            "the report does not record that the anticipated confound failed "
            "to appear, so a reader cannot tell the caveat was tested")
        assert "understates it" in md


def test_the_stream_split_is_derived_and_its_labelling_cause_stated(d):
    """The 20x rate gap must not read as a charge of neglect."""
    for k, v in d["basis_all"].items():
        assert d["basis_rate"][k] == pytest.approx(
            100 * d["basis_toxicity"].get(k, 0) / v, abs=0.02)
    c04, adj = d["basis_rate"]["C04"], d["basis_rate"]["adjacent"]
    assert d["basis_rate_ratio"] == pytest.approx(adj / c04, abs=0.1)
    assert sum(d["basis_toxicity"].values()) == d["organ_toxicity_articles"]
    assert sum(d["basis_all"].values()) == d["ferroptosis_articles"]
    md = " ".join(MD.read_text().split())
    assert "mostly a labelling effect and should not be read as neglect" in md
    assert "indexed under `Cardiotoxicity`" in md, (
        "the labelling mechanism is asserted without the example that shows it")


def test_the_first_framing_is_retracted_in_the_report(d):
    """It called these papers cancer-patient organ damage and most are not."""
    md = " ".join(MD.read_text().split())
    assert "corrects this analysis's own first framing" in md
    assert "Nile tilapia" in md, (
        "the retraction states the error without the exhibit that establishes "
        "it, which is what makes a retraction checkable")


def test_it_refuses_to_read_as_a_refutation(d):
    md = " ".join(MD.read_text().split())
    assert "It does not refute the selectivity premise" in md
    assert "says nothing directly about RSL3" in md
    for overclaim in ("the premise is false", "disproves", "refutes the claim",
                      "normal cells do not resist"):
        assert overclaim not in md.lower()


def test_the_layer_freeze_target_is_reported_as_unnamed(d):
    """CONTRIBUTING.md requires a named calibration target before a new
    phenotype lands, and an article count is not one. The report must say so
    rather than let a literature finding read as licence to build."""
    md = " ".join(MD.read_text().split())
    assert "the target remains unnamed" in md
    assert "an article count is not a dose-response" in md
    # And the claim about the committed curves must be true of the file.
    header = next(csv.reader(CTRP.open()))
    for col in header:
        assert not any(w in col.lower() for w in ("normal", "tissue", "tumor",
                                                  "tumour", "malignan")), (
            f"CTRPv2 curves carry a {col!r} column, so the report's claim that "
            "they hold no normal-tissue annotation is stale")


def test_the_manuscript_states_the_selectivity_assumption_as_one(d):
    """Section 8.4 now names the premise the model cannot check."""
    txt = " ".join((REPO / "article/drafts/v1.md").read_text().split())
    assert "nowhere argued because it is nowhere stated" in txt
    assert f"{d['organ_toxicity_articles']:,} of the {d['ferroptosis_articles']:,}" in txt
    # It must keep the CAF correction, which is the part most easily lost: the
    # Stromal phenotype reads like a normal-tissue compartment and is not one.
    assert "tumour-resident and parameterised for shielding" in txt
    assert "not a therapeutic index" in txt
    # And it must not over-read the literature it now cites.
    assert "does not refute the premise" in txt
    assert "quantitative question about therapeutic window" in txt


def test_the_no_normal_tissue_column_claim_is_still_true():
    """A handwritten claim ABOUT the document, so it needs checking.

    Section 8.4 asserts that no Chapter 5-7 results table carries a
    normal-tissue column and no preregistered prediction is scored on one. Both are the kind
    of sentence that silently goes false the moment somebody adds the column --
    which is the point at which the paragraph should be rewritten rather than
    left standing.
    """
    import re

    md = (REPO / "article/drafts/v1.md").read_text()
    start = md.index("## Chapter 5")
    end = md.index("## Chapter 8")
    body = md[start:end]
    # Table header rows only: "| A | B | C |"
    headers = [ln for ln in body.splitlines()
               if ln.strip().startswith("|") and "---" not in ln]
    bad = re.compile(
        r"\bnormal (tissue|cell)|\bhealthy\b|\btherapeutic index\b|\bstromal\b",
        re.I)
    hits = [h for h in headers if bad.search(h)]
    assert not hits, (
        "a Chapter 5-7 table now carries a normal-tissue column, so Section "
        f"8.4's claim that none does is false: {hits[:2]}")
    prereg = (REPO / "PREREGISTRATION.md").read_text()
    part1 = prereg[prereg.index("## Part 1"):]
    part1 = part1[:part1.index("## Part 2")] if "## Part 2" in part1 else part1
    assert not bad.search(part1), (
        "a preregistered prediction now scores on normal tissue, so Section "
        "8.4's claim is stale")
