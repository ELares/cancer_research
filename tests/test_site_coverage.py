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
    assert sum(c for _s, c in _mod()._pairs(d["sites"])) >= a, (
        "the per-site counts sum below the assigned total, which is impossible "
        "if multi-site articles are counted once per site")


def test_the_spread_is_reported_because_a_ratio_divides_into_it():
    """Total coverage can look fine while cross-site comparison is unusable."""
    d, md = _doc(), MD.read_text()
    vals = [c for _s, c in _mod()._pairs(d["sites"])]
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
    # regenerate-and-diff: the file must be what the tree walk produces
    dm = m.deep_map()
    for site in m.SITES:
        assert set(by_site[site]) == dm[site]["deep"], (
            f"{site}: the committed map is not what the tree walk resolves to")
        assert m.SITES[site] <= dm[site]["deep"], (
            f"{site}: the deep list does not contain its own shallow roots")
    # every deep descriptor must be a real C04 descriptor, not an invention
    c04 = set(m.c04_labels().values())
    for site, descs in by_site.items():
        unknown = set(descs) - c04
        assert not unknown, (
            f"{site} maps to {sorted(unknown)}, which are not in the committed "
            "C04 descriptor file")
    # and the count the report quotes must match the file it points at
    assert len(rows) == _doc()["variants"]["deep"]["n_descriptors"], (
        f"the committed map has {len(rows)} rows and the report says "
        f"{_doc()['variants']['deep']['n_descriptors']} descriptors")


def test_membership_is_the_tree_not_a_name_match():
    """The first version matched descriptor NAMES and shipped five traps.

    `ganGLIOn cysts`, `paraganGLIOma`, `ganglioneuroma` reached brain/CNS;
    the benign salivary `adenoLYMPHOMA` reached lymphoma; a lung disease
    reached cervix/uterus via `lymphangioLEIOMYOMAtosis`. A tree walk cannot
    do that, because a descriptor either sits under a node or it does not.
    """
    m = _mod()
    src = SCRIPT.read_text()
    assert "DEEP_RULES" not in src and "STRICT_EXCLUSIONS" not in src, (
        "the hand-written per-site name rules are back")
    dm = m.deep_map()
    tree = m.c04_tree()
    lab = m.c04_labels()
    by_label = {v: k for k, v in lab.items()}
    traps = {
        "brain/CNS": ["ganglion cysts", "paraganglioma", "ganglioneuroma",
                      "neuroectodermal tumors, primitive, peripheral"],
        "lymphoma": ["adenolymphoma", "multiple myeloma"],
        "cervix/uterus": ["lymphangioleiomyomatosis"],
        "head and neck": ["craniopharyngioma"],
    }
    for site, bad in traps.items():
        for x in bad:
            if x not in by_label:
                continue
            assert x not in dm[site]["deep"], (
                f"{x!r} is under {site} again -- it is not beneath that "
                "site's tree nodes, so only a name match could put it there")
    # every placement must be justified by an actual tree relation
    for site, v in dm.items():
        for x in v["deep"]:
            ts = tree.get(by_label[x], set())
            assert any(t == r or t.startswith(r + ".")
                       for t in ts for r in v["roots"]), (
                f"{site} claims {x!r}, which sits at {sorted(ts)} and under "
                f"none of {v['roots']}")


def test_the_deep_column_declares_where_the_tree_merges_this_pages_sites():
    """`Head and Neck Neoplasms` subsumes oesophagus and thyroid in MeSH."""
    d, md = _doc(), MD.read_text()
    m = _mod()
    dm = m.deep_map()
    want = {a: sorted(b for b in m.SITES
                      if b != a and m.SITES[b] <= dm[a]["deep"])
            for a in m.SITES
            if any(b != a and m.SITES[b] <= dm[a]["deep"] for b in m.SITES)}
    assert d.get("deep_site_overlaps") == want, (
        f"the artifact records {d.get('deep_site_overlaps')} overlaps, the "
        f"tree gives {want}")
    for a, bs in want.items():
        assert f"`{a}` subsumes" in md, (
            f"the deep list for `{a}` contains {bs}, which this page lists as "
            "separate sites, and the report does not say so")
    if want:
        assert "double-counts across the page's own list" in md


def test_the_depth_is_measured_rather_than_traded_away_in_prose():
    """"a deep subtree walk would raise assignability, at the cost of a mapping
    nobody can audit" was asserted in four places and computed in none.
    """
    d, md, src = _doc(), MD.read_text(), SCRIPT.read_text()
    var = d.get("variants") or {}
    for k in ("shallow", "deep"):
        assert var.get(k, {}).get("assigned"), f"the {k} variant was not run"
    assert var["deep"]["assigned"] > var["shallow"]["assigned"], (
        "the deep list assigns no more than the shallow one, so the page's "
        "measured gain needs rewriting rather than restating")
    assert var["shallow"]["assigned"] == d["assigned"]
    for k in ("deep",):
        assert f"{100*var[k]['assigned']/d['census']:.1f}%" in md, (
            f"the {k} assignability is measured and not rendered")
    assert "nobody can audit" in md and "NEITHER HALF HAD BEEN MEASURED" in md, (
        "the unmeasured tradeoff is neither retracted nor replaced")
    # the per-site ratios must be rendered, since they are the finding
    for s, c in _mod()._pairs(d["sites"]):
        cd = dict(_mod()._pairs(var["deep"]["sites"])).get(s, 0)
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


def test_the_ranking_survives_a_json_round_trip():
    """`sort_keys=True` reordered the stored rankings alphabetically, so
    `--render-only` -- the script's own second documented invocation --
    regenerated a DIFFERENT report: every rank cell became a self-identity and
    the rank-change clause vanished. #730 fixed exactly this once already.
    """
    m, d = _mod(), _doc()
    assert isinstance(d["sites"], list), (
        "the per-site ranking is stored as a dict, which json.dumps("
        "sort_keys=True) reorders alphabetically")
    for k, v in (d.get("variants") or {}).items():
        assert isinstance(v["sites"], list), f"{k}'s ranking is a dict"
    # the committed report must be what the renderer produces from the
    # committed artifact, which is what --render-only does
    assert m.render(d) == MD.read_text(), (
        "the committed report is not what --render-only reproduces")
    # and rendering must not depend on key order at all
    shuffled = json.loads(json.dumps(d, sort_keys=True))
    assert m.render(shuffled) == MD.read_text(), (
        "the report changes when the artifact's keys are sorted, so the "
        "ranking is being read from something order-dependent")


def test_the_map_root_column_actually_places_its_descriptor():
    """`write_map`'s middle column is the whole auditability claim, and no
    guard read it: falsifying 95 of 233 roots left every test green.
    """
    m = _mod()
    tsv = REPO_ROOT / "analysis" / "site-descriptor-map.tsv"
    tree = m.c04_tree()
    lab = m.c04_labels()
    by_label = {v: k for k, v in lab.items()}
    n = 0
    for line in tsv.read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        site, root, desc = line.split("\t")
        ts = tree.get(by_label.get(desc, ""), set())
        assert any(t == root or t.startswith(root + ".") for t in ts), (
            f"{site}: {desc!r} sits at {sorted(ts)}, none of which is under "
            f"the root {root!r} the map says placed it")
        n += 1
    assert n >= 200


def test_the_tree_file_is_pinned():
    """Nothing pinned the 31KB committed input the whole walk reads. Dropping
    one row and regenerating consistently left every guard green.
    """
    m = _mod()
    tree, lab = m.c04_tree(), m.c04_labels()
    assert set(tree) == set(lab), (
        f"the tree file covers {len(tree)} descriptors and the definition "
        f"file {len(lab)}; a partial fetch reads as a smaller tree")
    pairs = sum(len(v) for v in tree.values())
    assert pairs >= len(lab), "fewer tree numbers than descriptors"
    # ancestor-closed: every internal node of a tree number must itself be a
    # tree number somewhere, or the file was truncated by depth
    allt = {t for v in tree.values() for t in v}
    for t in sorted(allt):
        parts = t.split(".")
        for i in range(2, len(parts)):
            anc = ".".join(parts[:i])
            assert anc in allt or not anc.startswith("C04."), (
                f"{t} has no ancestor {anc} in the file, so it is truncated")
    for t in allt:
        assert t.startswith("C04"), f"{t} is not a C04 tree number"
    # THE FILE DECLARES ITS OWN SIZE, like the descriptor file beside it, so a
    # partial fetch fails against the file rather than being invisible. It is
    # written by `atlas_baseline.py --refresh-mesh` now; its header used to
    # name that command while nothing in the repo wrote it.
    import re as _re
    head = (REPO_ROOT / "corpus" / "atlas" / "mesh"
            / "c04-tree-numbers.tsv").read_text().splitlines()
    dec = next((l for l in head if l.startswith("# descriptors:")), None)
    assert dec, "the tree file does not declare its own size"
    n_d, n_p = (int(x) for x in _re.findall(r"(\d+)", dec))
    assert n_d == len(tree), f"header says {n_d} descriptors, file has {len(tree)}"
    assert n_p == pairs, f"header says {n_p} pairs, file has {pairs}"
    src = (REPO_ROOT / "scripts" / "atlas_baseline.py").read_text()
    assert "def fetch_c04_tree_numbers(" in src, (
        "nothing in the repo writes this committed input, and its header "
        "names a command that does not produce it")


def test_the_rank_change_finding_cannot_be_deleted_silently():
    """`moved = []` deleted the "13 of 18 sites change rank" clause and every
    guard stayed green, including the meta-guard written for that shape.
    """
    m, d, md = _mod(), _doc(), MD.read_text()
    var = d["variants"]
    sh = m._pairs(var["shallow"]["sites"])
    dp = m._pairs(var["deep"]["sites"])
    r_sh = {s: i + 1 for i, (s, _c) in enumerate(sh)}
    r_dp = {s: i + 1 for i, (s, _c) in enumerate(dp)}
    merged = set(d.get("deep_site_overlaps") or {})
    moved = [s for s in r_sh
             if s not in merged and r_dp.get(s) and r_sh[s] != r_dp[s]]
    if moved:
        assert f"{len(moved)} of {len(r_sh)} sites change rank" in md, (
            f"{len(moved)} sites change rank between the two lists and the "
            "report does not say so")
    # a site whose subtree contains another of these sites gets no rank at all
    for s in merged:
        assert f"n/a (subsumes" in md
        assert f"| {s} |" in md


def test_the_scope_cost_names_more_than_benign_and_precursor():
    """9.9% of placements are experimental models, veterinary disease or named
    syndromes, and the page's only stated scope cost was benign/precursor.
    """
    d, md = _doc(), MD.read_text()
    nh = d.get("non_human_disease_placements") or {}
    if not nh:
        return
    tot = sum(len(v) for v in nh.values())
    assert f"{tot} placements are experimental" in md, (
        f"{tot} placements are not human disease at that site and the report "
        "does not say so")
    biggest = max(nh, key=lambda k: len(nh[k]))
    for x in nh[biggest][:4]:
        assert f"`{x}`" in md


def test_the_residue_is_characterised_by_its_own_descriptors():
    """The prose pointed at a descriptor list accumulated over EVERY
    unassigned record, including the two buckets it had just excluded.
    """
    d, md = _doc(), MD.read_text()
    u = d["unassigned"]
    res = u.get("residue_top_descriptors")
    assert res, "the residue has no descriptor profile of its own"
    allu = dict(u.get("top_descriptors") or [])
    resd = dict(res)
    assert "neoplasms" not in resd, (
        "the residue profile counts generic-`Neoplasms` records, which the "
        "residue excludes by definition")
    for k, v in resd.items():
        assert v <= allu.get(k, v), (
            f"{k} is commoner in the residue ({v:,}) than in every unassigned "
            f"record ({allu.get(k)}), which cannot be")
    assert "cancer at a site this list does not cover" in md.lower()


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
