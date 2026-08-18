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


def _mod():
    import importlib.util
    spec = importlib.util.spec_from_file_location("sc", SCRIPT)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


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
    assert d["adjacent_basis"] + a <= n + d["multi_site"]
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


def test_the_shallow_map_stays_readable_and_the_deep_map_is_committed():
    """The shallow list must stay small enough to read; the deep one must be
    ON DISK. A cap on the shallow list is fine -- but the cap is what made
    DEPTH non-uniform (one descriptor for `stomach`, four for `brain/CNS`),
    and the page's old answer to that was a sentence nobody had measured.
    """
    m, src = _mod(), SCRIPT.read_text()
    assert "SITES = {" in src, "the site map is no longer a literal"
    assert len(m.SITES) >= 15, f"only {len(m.SITES)} sites mapped"
    for site, descs in m.SITES.items():
        assert 1 <= len(descs) <= 6, (
            f"{site} has {len(descs)} descriptors; the shallow list exists to "
            "be readable in one screen")
        assert all(d == d.lower() for d in descs)
    # the deep map is derived from the COMMITTED cancer definition, and is
    # itself committed, which is the whole answer to "nobody can audit it"
    tsv = REPO_ROOT / "analysis" / "site-descriptor-map.tsv"
    assert tsv.exists(), "the resolved deep map is not committed"
    rows = [l.split("\t") for l in tsv.read_text().splitlines()
            if l and not l.startswith("#")]
    assert len(rows) >= 200, f"the deep map has only {len(rows)} rows"
    by_site = {}
    for site, variant, desc in rows:
        by_site.setdefault(site, {})[desc] = variant
    assert set(by_site) == set(m.SITES), (
        "the committed map covers different sites than SITES")
    # regenerate-and-diff: the file must be what the rules produce
    dm = m.deep_map()
    for site in m.SITES:
        assert set(by_site[site]) == dm[site]["deep"], (
            f"{site}: the committed map is not what DEEP_RULES resolves to -- "
            "re-run the generator")
        strict = {d for d, v in by_site[site].items() if v == "deep+strict"}
        assert strict == dm[site]["strict"], f"{site}: strict variant is stale"
        # and every deep list must CONTAIN the shallow one, or the two columns
        # beside each other are not the same site
        assert m.SITES[site] <= dm[site]["deep"], (
            f"{site}: the deep list drops "
            f"{sorted(m.SITES[site] - dm[site]['deep'])}, so the deep column "
            "is not a superset of the shallow one it is compared against")
    # every deep descriptor must be a real C04 descriptor, not an invention
    c04 = {x.lower() for x in m.c04_labels()}
    for site, descs in by_site.items():
        unknown = set(descs) - c04
        assert not unknown, (
            f"{site} maps to {sorted(unknown)}, which are not in the committed "
            "C04 descriptor file")


def test_the_depth_is_measured_rather_than_traded_away_in_prose():
    """"a deep subtree walk would raise assignability, at the cost of a mapping
    nobody can audit" was asserted in four places and computed in none.
    """
    d, md, src = _doc(), MD.read_text(), SCRIPT.read_text()
    var = d.get("variants") or {}
    for k in ("shallow", "deep", "strict"):
        assert var.get(k, {}).get("assigned"), f"the {k} variant was not run"
    assert var["deep"]["assigned"] > var["shallow"]["assigned"], (
        "the deep list assigns no more than the shallow one, so the page's "
        "measured gain needs rewriting rather than restating")
    assert var["shallow"]["assigned"] == d["assigned"]
    for k in ("deep", "strict"):
        assert f"{100*var[k]['assigned']/d['census']:.1f}%" in md, (
            f"the {k} assignability is measured and not rendered")
    assert "nobody can audit" in md and "NEITHER HALF HAD BEEN MEASURED" in md, (
        "the unmeasured tradeoff is neither retracted nor replaced")
    # the per-site ratios must be rendered, since they are the finding
    for s, c in d["sites"].items():
        cd = var["deep"]["sites"].get(s, 0)
        assert f"{cd/max(c,1):.2f}x" in md, f"{s}'s depth ratio is not rendered"


def test_the_remainder_is_decomposed_rather_than_narrated():
    """"much cancer literature is about biology, methods or cancer in general"
    was a story about 1.86M articles, and it is wrong for a large share.
    """
    d, md = _doc(), MD.read_text()
    u = d.get("unassigned") or {}
    assert u.get("total") == d["census"] - d["assigned"]
    parts = ("same_sites_deeper", "generic_neoplasms", "no_c04_descriptor")
    for k in parts:
        assert u.get(k) is not None, f"the unassigned pile is not split by {k}"
        assert f"{u[k]:,}" in md, f"{k} is measured and not rendered"
    assert u["same_sites_deeper"] + u["generic_neoplasms"] <= u["total"]
    # the load-bearing claim: a large share is NOT site-less
    share = u["same_sites_deeper"] / u["total"]
    assert share > 0.10, (
        f"only {100*share:.1f}% of the unassigned are the same sites deeper; "
        "the page's correction of the old narrative needs re-checking")
    assert "narrated rather than measured" in md
    assert u.get("top_descriptors"), "the residue is not characterised at all"


def test_the_dead_no_mesh_row_is_gone_and_the_real_exclusions_are_stated():
    """That row could only ever print 0: the stream admits on a MeSH match."""
    d, md = _doc(), MD.read_text()
    assert "no_mesh" not in d, (
        "the artifact still carries a count that the admission rule fixes at "
        "zero, which reads as a property of the literature")
    assert "carry no MeSH at all" in md and "COULD NOT HAVE BEEN ANYTHING ELSE" in md
    ex = d.get("excluded_streams") or {}
    assert ex.get("text_matched_no_mesh"), (
        "the MeSH-less census stream this denominator excludes is not counted")
    both = d["census"] + ex["text_matched_no_mesh"]
    assert f"{100*d['assigned']/both:.1f}%" in md, (
        "assignability over BOTH census streams is not stated, so the "
        "headline share is quoted against a denominator chosen by exclusion")
    assert d.get("adjacent_basis"), "adjacent-basis records are not counted"
    assert f"{d['adjacent_basis']:,}" in md


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
