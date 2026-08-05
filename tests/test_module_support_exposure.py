"""A zero in the module-support table must not be read as a verdict on the claim.

`atlas-module-support.md` counts how many of 20 simulation-module claims the
census corroborates, and for a long time it reported the shortfall as a flat
number -- 11 of 20 -- as though the nine zeros meant the same thing.

They do not. A pair can only be asserted if BOTH its entities are written about,
and across these claims the weaker entity's partner count spans three orders of
magnitude. Seven of the nine zeros sit below the exposure of the WEAKEST claim
that does have support, so for those the census could not have corroborated the
claim whatever the biology is. That is a fact about what has been studied.

The two that sit above the line are the only rows where a zero says anything
about the claim itself, and they are the ones worth a reader's attention.

These guards pin the structure, not the numbers: the section exists, its
classification agrees with the artifact, and the distinction it draws is real.
"""

import json
import re
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "scripts"))

DOC = REPO_ROOT / "analysis" / "atlas-module-support.md"


def _section() -> str:
    txt = DOC.read_text()
    if "## A zero is usually about EXPOSURE" not in txt:
        return ""
    start = txt.index("## A zero is usually about EXPOSURE")
    end = txt.index("## Detail", start)
    return txt[start:end]


def test_the_exposure_section_is_present():
    assert _section(), (
        "the module-support document reports a corroboration count without "
        "saying that most of its zeros are exposure artifacts")


def test_every_zero_row_is_classified():
    """No zero may be left uncategorised -- that is the old flat reading."""
    sec = _section()
    rows = re.findall(r"^\| (\w+) \| `([^`]+)` - `([^`]+)` \| (\d+) \(", sec, re.M)
    assert rows, "the classification table is missing"
    for mod, a, b, weaker in rows:
        line = next(ln for ln in sec.split("\n") if f"| {mod} | `{a}`" in ln)
        assert line.rstrip().endswith("yes |") or line.rstrip().endswith("no |"), (
            f"{mod} is listed without a below-the-range verdict")


def test_entity_degree_counts_both_endpoints():
    """The degree definition itself, checked independently of its own output.

    The classification test recomputes the floor USING `_entity_degree`, so it
    agrees with any wrong definition of degree. This asserts the invariant that
    definition must satisfy: every edge contributes to exactly two entities, so
    the degrees sum to twice the edge count. A version that counted one side
    survived every other guard here.
    """
    from atlas_baseline import atlas_root
    from atlas_graph import load_index
    from atlas_module_support import _entity_degree

    idx = load_index(atlas_root())
    deg = _entity_degree(idx)
    n_edges = len(idx["edges"])
    self_pairs = sum(1 for k in idx["edges"] if k[0] == k[1])
    assert sum(deg.values()) == 2 * n_edges - self_pairs, (
        "degree does not count each edge's endpoints exactly once each. A "
        "version counting one side halves every exposure number; a version "
        "using `for i in key` double-counts self-pairs, which AIFM2 is one of "
        "and which reaches two of the twenty claims")


def test_the_classification_matches_the_measured_threshold():
    """The split must follow from the data, not from a written-in cut.

    Recomputed here from the graph so a hand-edited table cannot disagree with
    the measurement it claims to summarise.
    """
    from atlas_baseline import atlas_root
    from atlas_graph import load_index, resolve, support
    from atlas_module_support import CLAIMS, _entity_degree

    idx = load_index(atlas_root())
    deg = _entity_degree(idx)
    weaker, totals = {}, {}
    for mod, a, b, _pmid, _claim in CLAIMS:
        ia, ib = resolve(idx, a), resolve(idx, b)
        if not ia or not ib:
            continue
        s = support(idx, a, b)
        weaker[(mod, a, b)] = min(deg.get(ia, 0), deg.get(ib, 0))
        totals[(mod, a, b)] = s["n_articles"] if s else 0

    supported = [w for k, w in weaker.items() if totals[k] > 0]
    assert supported, "no claim has support; the threshold is undefined"
    floor = min(supported)

    sec = _section()
    for (mod, a, b), w in weaker.items():
        if totals[(mod, a, b)] > 0:
            continue
        line = next((ln for ln in sec.split("\n") if f"| {mod} | `{a}`" in ln), None)
        assert line, f"{mod} has zero support but is not in the table"
        expect = "| no |" if w >= floor else "| yes |"
        assert line.rstrip().endswith(expect.strip()), (
            f"{mod} (weaker degree {w}, floor {floor}) is marked as the "
            f"opposite of what the graph says")


def test_the_reported_correlation_equals_the_recomputed_one():
    """rho must be RECOMPUTED, not merely plausible.

    The first version asserted only `> 0.5`, so a hand-edited 0.99 passed, and
    so did a broken estimator whose real-data output was 0.8595. Pin it to the
    artifact, which the manifest gate in turn pins to the tree.
    """
    raw = json.loads((REPO_ROOT / "analysis" / "atlas-module-support.json").read_text())
    sec = _section()
    m = re.search(r"rho = ([\d.]+)\*\*", sec)
    assert m, "the section states no correlation"
    assert abs(float(m.group(1)) - raw["spearman_weaker_degree_vs_relations"]) < 0.005, (
        f"the document says rho={m.group(1)}, the artifact says "
        f"{raw['spearman_weaker_degree_vs_relations']:.4f}")
    assert raw["spearman_weaker_degree_vs_relations"] > 0.5, (
        "exposure no longer predicts support, so the section's central claim "
        "does not hold and its reading of the zeros must be revisited")


def test_the_rank_estimators_handle_ties():
    """Ties are the branch that matters here and the branch that was untested.

    Nine of twenty outcomes are tied at zero, so tie-averaging is load-bearing;
    the original unit checks used tie-free inputs only, and a mutation that
    replaced the tie-average with a raw index survived them.
    """
    from atlas_module_support import _kendall, _spearman

    assert _spearman([1, 2, 3, 4], [1, 2, 3, 4]) == 1.0
    assert _spearman([1, 2, 3, 4], [4, 3, 2, 1]) == -1.0
    # Perfectly tied ranks on one side must give exactly zero.
    assert abs(_spearman([1, 1, 2, 2], [1, 2, 1, 2])) < 1e-12
    assert abs(_spearman([5, 5, 5, 5], [1, 2, 3, 4])) < 1e-12
    # THREE tied groups with an exact known answer. The cases above do not
    # discriminate: with only two distinct groups a rank variable is effectively
    # binary, and correlation is invariant to which two values the groups get --
    # so a tie handler that assigns the group's FIRST index instead of its mean
    # passes them. With three groups the scale is distorted and the value moves.
    # xs=[1,1,2,3,3] vs ys=[1..5] gives exactly 9/sqrt(90).
    assert abs(_spearman([1, 1, 2, 3, 3], [1, 2, 3, 4, 5]) - 9 / 90 ** 0.5) < 1e-9
    assert abs(_kendall([1, 2, 3, 4], [1, 2, 3, 4]) - 1.0) < 1e-12
    assert abs(_kendall([1, 1, 2, 2], [1, 2, 1, 2])) < 1e-12


def test_the_section_does_not_name_a_measure_dependent_exception_set():
    """The first draft named two rows as 'the interesting ones'. They are not.

    Running the same procedure on the pair-level co-mention column already in
    the table inverts the correlation and returns a disjoint pair, so naming
    either set as the genuine unasserted links reports an artifact of the
    measure chosen. The section must disclose that rather than pick a side.
    """
    sec = _section()
    assert "does not replicate across measures" in sec, (
        "the section names exceptions without disclosing that the pair-level "
        "measure in the same table gives a disjoint answer")
    assert "no overlap" in sec
    # The phrase may appear ONLY inside the sentence that retracts it. Banning
    # it outright would forbid the document from saying what it got wrong.
    for para in sec.split("\n\n"):
        if "the interesting rows" in para:
            assert "first draft" in para, (
                "the retracted 'interesting rows' framing is asserted rather "
                "than quoted in its retraction")


def test_the_undetectable_claim_is_not_made():
    """'could not have corroborated whatever the biology is' was false.

    45% of all asserted pairs in the graph sit below the line, so it is not a
    detectability limit. The base rate must be stated, not the modal claim.
    """
    sec = _section()
    assert "whatever the biology is" not in sec or "That is false" in sec, (
        "the section asserts the census could not have found these pairs")
    assert re.search(r"\*\*\d+% of all asserted pairs\*\*", sec), (
        "the base rate that refutes the modal reading is not reported")


def test_the_floor_is_disclosed_as_a_sample_minimum():
    """It is the min of 11 points, set by a 1-article row."""
    sec = _section()
    assert "sample minimum, not a threshold" in sec
    assert "Leave-one-out" in sec


def test_the_coupling_caveat_is_stated():
    """The correlation is partly mechanical and the document must say so."""
    sec = _section()
    assert "share a term" in sec or "not fully independent" in sec, (
        "the section reports a correlation between two quantities that share a "
        "term without disclosing it")
    assert "not evidence FOR a claim" in sec, (
        "a reader could take a high partner count as support for the claim")
