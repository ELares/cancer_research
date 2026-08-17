"""Guards for the site-assignment gate on burden weighting (#729).

WHAT THIS IS
------------
#729 wants literature-per-death per cancer site. A reviewer's gate was that
site-assignment completeness must be measured and reported BEFORE any ratio,
because a per-site rate is only as good as the assignment underneath it.

Measured: 57.8% of the census is assignable, the per-site spread is 6x, and 8.2%
of assigned articles touch more than one site. The gate passes with conditions,
and the conditions are the deliverable.

THE FAILURE MODES THIS GUARDS
------------------------------
1. ANSWERING THE QUESTION IT WAS SUPPOSED TO GATE. If this script ever computes
   a burden ratio, it has answered the question before establishing whether it
   can be answered. A guard fails if mortality data appears here.

2. AN UNAUDITABLE DENOMINATOR. The site map is shallow and checkable on purpose.
   A deep subtree walk would raise assignability at the cost of a mapping nobody
   can review, and an unauditable denominator is worse than a conservative one.

3. HIDING THE DOUBLE COUNT. Multi-site articles are counted once per site, so
   the per-site column sums past the assigned total. That is correct for "how
   much literature touches this site" and wrong for a partition, and a ratio
   built on the wrong one would be silently inconsistent.
"""

import json
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPT = REPO_ROOT / "scripts" / "atlas_site_coverage.py"
MD = REPO_ROOT / "analysis" / "atlas-site-coverage.md"
JSON_OUT = REPO_ROOT / "analysis" / "atlas-site-coverage.json"


def _doc():
    return json.loads(JSON_OUT.read_text())


def test_it_does_not_compute_the_ratio_it_gates():
    """Answering first would defeat the purpose of gating."""
    src, md = SCRIPT.read_text(), MD.read_text()
    # IDENTIFIERS, not text. The word GLOBOCAN legitimately appears in the
    # report's caveat about site-definition mismatch, which lives in a string
    # literal inside the renderer -- so both a whole-file check and a
    # docstring-stripped check fail on prose that is doing the right thing.
    # What must not exist is mortality data being USED: a name, an import or an
    # attribute. Third time today a guard has tripped on the sentence that
    # disclaims the thing it forbids.
    import ast
    tree = ast.parse(src)
    used = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Name):
            used.add(node.id.lower())
        elif isinstance(node, ast.Attribute):
            used.add(node.attr.lower())
        elif isinstance(node, (ast.Import, ast.ImportFrom)):
            used.add((getattr(node, "module", "") or "").lower())
            for a in node.names:
                used.add(a.name.lower())
    for banned in ("globocan", "mortality", "deaths", "per_death"):
        hits = [u for u in used if banned in u]
        assert not hits, (
            f"{hits} used as identifiers in the gate script; this measures "
            "whether a burden ratio is computable and must not compute one")
    assert "does NOT compute a burden ratio" in src or \
        "NOT compute a burden ratio" in md, (
        "the report no longer says it is a gate rather than the analysis")


def test_the_assignment_arithmetic_holds():
    d = _doc()
    n, a = d["census"], d["assigned"]
    assert 0 < a < n, f"assigned {a:,} of {n:,} is not a proper subset"
    assert d["no_mesh"] + a <= n
    assert d["multi_site"] <= a, "more multi-site articles than assigned ones"
    # the per-site column double-counts by design; it must exceed the assigned
    # total, or the multi-site figure is not what it says
    assert sum(d["sites"].values()) >= a, (
        "the per-site counts sum below the assigned total, which is impossible "
        "if multi-site articles are counted once per site")


def test_the_spread_is_reported_because_a_ratio_divides_into_it():
    """Total coverage can look fine while cross-site comparison is unusable."""
    d, md = _doc(), MD.read_text()
    vals = list(d["sites"].values())
    assert len(vals) >= 10, f"only {len(vals)} sites; the spread is not meaningful"
    lo, hi = min(vals), max(vals)
    assert f"{lo:,} to {hi:,}" in md, (
        "the report does not state its own measured spread, which is what a "
        "burden ratio divides mortality into")
    assert "factor of" in md


def test_the_double_count_is_disclosed():
    d, md = _doc(), MD.read_text()
    share = 100 * d["multi_site"] / max(d["assigned"], 1)
    assert f"{share:.1f}% of" in md, (
        "the multi-site share is not the one the artifact supports")
    assert "wrong for a partition" in md, (
        "the report no longer distinguishes 'touches this site' from a "
        "partition, so a ratio could be built on the wrong denominator")


def test_the_site_map_stays_shallow_and_auditable():
    """A deep walk would raise coverage at the cost of reviewability."""
    src = SCRIPT.read_text()
    assert "SITES = {" in src, "the site map is no longer a literal"
    import importlib.util
    spec = importlib.util.spec_from_file_location("sc", SCRIPT)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    assert len(m.SITES) >= 15, f"only {len(m.SITES)} sites mapped"
    for site, descs in m.SITES.items():
        assert 1 <= len(descs) <= 6, (
            f"{site} maps {len(descs)} descriptors; the list is meant to stay "
            "short enough for a reader to check")
        assert all(d == d.lower() for d in descs), (
            f"{site} has a non-lowercase descriptor; matching is case-folded "
            "and a stray capital would silently never match")


def test_the_geography_caveat_survives():
    """A literature-per-death ratio partly measures where science is funded."""
    md = MD.read_text()
    assert "where science is funded" in md, (
        "the report no longer states that mortality and publication counts "
        "have different geography, which is the caveat that stops the eventual "
        "ratio being read as a neglect verdict")


def test_an_empty_assignment_refuses_to_render():
    src = SCRIPT.read_text()
    assert 'if d["assigned"] == 0:' in src
    assert "is not a finding" in src and "raise SystemExit" in src
