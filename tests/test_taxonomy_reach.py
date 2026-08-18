"""Guards for the mechanism taxonomy's measured field of view (#730).

WHAT THIS ANALYSIS CLAIMS
-------------------------
That the keyword tagger can label roughly 6% of the cancer literature, that the
documented 0.20%-41.86% per-mechanism capture spread therefore describes
variation inside that fraction, and that the unlabelled remainder is NOT simply
literature a mechanism taxonomy should ignore -- over a tenth of it carries
explicit antineoplastic-therapy descriptors.

THE WAYS THIS GOES WRONG
------------------------
1. MEASURING THE WRONG INSTRUMENT. There are two: the keyword tagger, which
   actually labels the corpus, and the MeSH leaf map, which is PRECISION-FIRST
   by design and deliberately drops umbrella descriptors. Reporting the map's
   low reach as the taxonomy's field of view would criticise a precision
   instrument for not being a coverage one. Both are measured; the report must
   keep them distinct.

2. REIMPLEMENTING THE MATCHER. The reach number is only about the production
   tagger if it uses the production matcher. A local copy of the substring rule
   would measure the copy.

3. AN UNINFORMATIVE PROFILE PASSING AS ANALYSIS. The first version's
   top-descriptor table read `humans 92%, female 45.6%, male 37.4%` -- all true,
   all MeSH check-tags, and silent on what the literature is about. The
   check-tag exclusion is what makes that section mean anything.

4. A PARSER THAT SILENTLY UNDERCOUNTS. A hand-rolled line parser missed five of
   the thirty-one map descriptors because some entries are multi-line lists,
   which would have understated the map's reach in the exact place the report
   makes a claim about it.
"""

import json
import re
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPT = REPO_ROOT / "scripts" / "atlas_taxonomy_reach.py"
MD = REPO_ROOT / "analysis" / "atlas-taxonomy-reach.md"
JSON_OUT = REPO_ROOT / "analysis" / "atlas-taxonomy-reach.json"


def _doc():
    return json.loads(JSON_OUT.read_text())


def _mod():
    import importlib.util
    spec = importlib.util.spec_from_file_location("reach", SCRIPT)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def test_it_uses_the_production_matcher_not_a_copy():
    """The claim is about the production tagger, so it must run it."""
    src = SCRIPT.read_text()
    assert "from tag_articles import" in src, (
        "the reach scan no longer imports the production matcher; a local copy "
        "would measure the copy rather than the tagger")
    assert "def text_matches_keyword" not in src, (
        "the matcher has been reimplemented inside the scan")
    assert "import config" in src and "config.MECHANISM_KEYWORDS" in src, (
        "the scan no longer reads the production keyword vocabulary")


def test_the_two_instruments_are_reported_separately():
    """Conflating them would criticise a precision tool for low coverage."""
    d, md = _doc(), MD.read_text()
    assert d["keyword_hits"] > 0 and d["mesh_leaf_hits"] > 0
    kw_rate = 100 * d["keyword_hits"] / d["sampled"]
    mesh_rate = 100 * d["mesh_leaf_hits"] / d["census_total"]
    assert f"{kw_rate:.2f}%" in md, (
        f"the report does not state the keyword reach {kw_rate:.2f}%")
    assert f"{mesh_rate:.2f}%" in md, (
        f"the report does not state the MeSH-map reach {mesh_rate:.2f}%")
    assert "precision instrument" in md, (
        "the report no longer says the MeSH map is precision-first, so its low "
        "reach reads as a coverage failure")
    # and in the renderer, since the committed report does not change when the
    # generator does until someone re-runs it
    assert "not** a coverage failure" in SCRIPT.read_text(), (
        "the renderer no longer emits the precision-instrument caveat")


def test_the_sampled_rate_carries_an_interval():
    """Sampled figures carry intervals; the full-census figure must not.

    There are TWO sampled reaches now -- the raw keyword loop and the
    production matcher -- and both legitimately carry one. The property that
    matters is unchanged: the MeSH-map figure is a full-census count and must
    not be dressed with a sampling interval it does not have.
    """
    d, md = _doc(), MD.read_text()
    n_ci = md.count("95% CI")
    assert n_ci >= 1, "no sampled figure carries an interval"
    sampled_rates = [d["keyword_hits"]]
    if d.get("production_hits"):
        sampled_rates.append(d["production_hits"])
    assert n_ci == len(sampled_rates), (
        f"{n_ci} intervals are reported for {len(sampled_rates)} sampled "
        "figures; the full-census MeSH figure has no sampling error and must "
        "not carry one")
    mesh_pct = 100 * d["mesh_leaf_hits"] / d["census_total"]
    i = md.find(f"{mesh_pct:.2f}%")
    if i >= 0:
        assert "95% CI" not in md[i:i + 60], (
            "the full-census MeSH figure is shown with a sampling interval")


def test_the_descriptor_profile_excludes_check_tags():
    """`humans 92%` is true and says nothing; the exclusion is the analysis."""
    m, d = _mod(), _doc()
    assert "humans" in m.DEMOGRAPHIC_TAGS and "female" in m.DEMOGRAPHIC_TAGS
    # the artifact stores ORDERED PAIRS now, so take the keys
    top = [k for k, _v in m._pairs(d["untagged_top_mesh"])]
    assert top, "no descriptor profile was produced"
    leaked = [x for x in top if x in m.DEMOGRAPHIC_TAGS]
    assert not leaked, (
        f"demographic check-tags {leaked} are in the top-descriptor table, so "
        "it is reporting demographics rather than subject matter")
    # STUDY-DESIGN descriptors are excluded too, but they are REPORTED
    # separately rather than silently dropped: the section's own question is
    # what the unlabelled remainder is about, and `cell line, tumor` and
    # `retrospective studies` answer it.
    assert d.get("untagged_top_study_design"), (
        "study-design descriptors are excluded and no longer reported")
    # AND the code must still exclude them. The assertion above reads the
    # committed artifact, which does not move when the scan is edited, so on
    # its own it would pass for a script that had stopped filtering.
    src = SCRIPT.read_text()
    assert "if m in DEMOGRAPHIC_TAGS:" in src and "continue" in src, (
        "the scan no longer excludes demographic check-tags")
    assert any("neoplasm" in t or "antineoplastic" in t for t in top), (
        "no oncology descriptor survived the exclusion, which suggests the "
        "filter is removing too much")


def test_the_map_is_parsed_with_a_real_yaml_loader():
    """A line parser missed 5 of 31 descriptors and would understate reach."""
    m = _mod()
    d = m.mesh_leaf_descriptors()
    assert len(d) >= 30, (
        f"only {len(d)} descriptors parsed from mesh-mechanism-map.yaml; the "
        "hand-rolled parser this replaced found 26 of 31")
    import yaml
    doc = yaml.safe_load(m.MESH_MAP.read_text())
    expected = {x.strip().lower()
                for b in doc["mechanisms"].values()
                for x in (b.get("descriptors") or [])}
    assert d == expected, (
        f"parsed descriptors differ from the yaml by "
        f"{sorted(d ^ expected)[:5]}")


def test_an_empty_result_refuses_to_render():
    """Zero matches is what a broken import looks like, not a finding."""
    src = SCRIPT.read_text()
    # The CONDITION, not the message beside it. Asserting only the prose let a
    # mutation replacing the test with `if False:` pass untouched.
    assert 'if d["keyword_hits"] == 0:' in src, (
        "the zero-result check no longer tests the keyword hit count, so a "
        "scan that matched nothing would write 0% reach as a measurement")
    assert "is not a finding" in src and "raise SystemExit" in src, (
        "the refusal message or the raise is gone")


def test_the_verdict_is_derived_not_asserted():
    """The falsifier's answer must come from the counts."""
    d, md = _doc(), MD.read_text()
    untagged = d["sampled"] - d["keyword_hits"]
    # THE UNION the report publishes, not a re-summed set of overlapping
    # rows. Summing them double-counts every article carrying both, which is
    # the defect this analysis was corrected for.
    ther = d.get("untagged_therapy_union")
    assert ther is not None, "the therapy union is no longer measured"
    assert ther > 0, (
        "no therapy descriptor appears in the unlabelled profile, which is the "
        "evidence the verdict section rests on")
    assert f"{100*ther/untagged:.1f}%" in md, (
        "the report's therapy-descriptor share is not the one its own JSON "
        "supports")
    assert "partly true and not sufficient" in md, (
        "the verdict no longer states that the comfortable reading is "
        "insufficient, which is the finding")


def test_render_only_works_without_the_census():
    """The prose must be rebuildable from the artifact alone."""
    # Run against a COPY of the tree's artifacts. Invoking the generator here
    # rewrites the committed report, which makes the suite non-idempotent and --
    # measured the hard way -- lets a mutation sweep corrupt a committed file:
    # a mutated renderer ran through this test and its text landed on disk.
    import shutil, tempfile
    with tempfile.TemporaryDirectory() as td:
        shutil.copy2(MD, Path(td) / MD.name)
        try:
            res = subprocess.run([sys.executable, str(SCRIPT), "--render-only"],
                                 cwd=REPO_ROOT, capture_output=True, text=True)
        finally:
            shutil.copy2(Path(td) / MD.name, MD)
    assert res.returncode == 0, (
        f"--render-only failed:\n{res.stdout}\n{res.stderr}")


def test_render_is_invariant_to_json_key_ordering():
    """The bug, stated as a property.

    `json.dumps(..., sort_keys=True)` reordered the stored dicts, so
    `--render-only` produced an ALPHABETICAL document while the documented
    command produced a RANKED one. The committed .md was a second-generation
    render and the default command would have silently rewritten it.

    Round-tripping the artifact through sorted keys must not change a single
    byte of the rendered page.
    """
    m = _mod()
    d = _doc()

    # (1) THE STORAGE SHAPE is what makes ordering safe. A round-trip of an
    # ordered LIST through sort_keys is identity, so asserting only the
    # invariant is vacuous once the fix is in -- it would pass again the
    # moment someone stored dicts.
    for key in ("per_mechanism", "untagged_top_mesh", "untagged_pubtypes"):
        assert isinstance(d[key], list), (
            f"{key} is stored as a {type(d[key]).__name__}; "
            "`json.dumps(sort_keys=True)` reorders a dict, which is what made "
            "`--render-only` produce an alphabetical document while the "
            "documented command produced a ranked one")
        assert all(isinstance(x, list) and len(x) == 2 for x in d[key]), (
            f"{key} is a list but not of [key, value] pairs")

    # (2) AND the reader must re-sort a dict, so a stale artifact in the old
    # shape cannot reintroduce the bug.
    got = m._pairs({"a": 1, "b": 9, "c": 5})
    assert [k for k, _v in got] == ["b", "c", "a"], (
        f"_pairs does not re-sort a dict by value (got {got}), so an artifact "
        "in the old shape would render in key order again")

    # (3) the invariant itself, which must hold in both shapes
    assert m.render(d) == m.render(json.loads(json.dumps(d, sort_keys=True))), (
        "rendering is sensitive to JSON key order")
    as_dict = {k: (dict(v) if k in ("per_mechanism", "untagged_top_mesh",
                                    "untagged_pubtypes") else v)
               for k, v in d.items()}
    assert m.render(as_dict) == m.render(d), (
        "rendering an artifact in the OLD dict shape does not reproduce the "
        "ranked document")


def test_both_ranked_tables_are_ranked():
    """An alphabetical slice omitted the largest descriptor and the "
    second-largest mechanism, under headings that promise rank."""
    m, d, md = _mod(), _doc(), MD.read_text()
    for key, heading in (("untagged_top_mesh", "descriptor"),
                         ("per_mechanism", "mechanism")):
        pairs = m._pairs(d[key])
        vals = [v for _k, v in pairs]
        assert vals == sorted(vals, reverse=True), (
            f"{key} is not in descending order in the artifact")
        top = pairs[0][0]
        assert f"| {top} |" in md, (
            f"the largest entry of {key} (`{top}`) is missing from the "
            f"rendered {heading} table, which is headed as if ranked")


def test_the_publication_type_table_partitions():
    """Presented as a breakdown of the remainder, it summed past 100%."""
    m, d = _mod(), _doc()
    s, kw = d["sampled"], d["keyword_hits"]
    untagged = s - kw
    tot = sum(v for _k, v in m._pairs(d["untagged_pubtypes"]))
    assert tot <= untagged, (
        f"the publication-type buckets sum to {tot:,} against {untagged:,} "
        f"unlabelled articles ({100*tot/untagged:.1f}%); a breakdown cannot "
        "exceed the thing it breaks down")
    src = SCRIPT.read_text()
    assert "break" in src.split("for bucket, needles in PUBTYPE_BUCKETS")[1][:400], (
        "the bucket loop no longer stops at the first match, so a record can "
        "land in several and the table stops being a partition")


def test_overlapping_shares_are_published_as_unions():
    """Two overlapping buckets summed is not a share of anything."""
    d = _doc()
    for key in ("untagged_therapy_union", "untagged_soft_union"):
        assert d.get(key) is not None, (
            f"{key} is gone; the report would fall back to summing overlapping "
            "rows, which double-counts every article carrying both")
    m = _mod()
    ther_sum = sum(v for k, v in m._pairs(d["untagged_top_mesh"])
                   if "antineoplastic" in k or "chemotherapy" in k
                   or "radiotherapy" in k or "drug therapy" in k)
    assert d["untagged_therapy_union"] <= ther_sum, (
        "the union exceeds the sum of its parts, which is impossible")


def test_the_production_reach_is_reported_and_distinguished():
    """The published figure was a raw keyword loop, not the production matcher."""
    d, md = _doc(), MD.read_text()
    assert d.get("production_hits"), "the production matcher is no longer measured"
    s = d["sampled"]
    assert f"{100*d['production_hits']/s:.2f}%" in md, (
        "the production reach is measured and not rendered")
    assert "PRODUCTION matcher" in md
    src = SCRIPT.read_text()
    assert "match_mechanisms" in src, (
        "the production entry point is no longer called, so the reported "
        "production reach is a reimplementation of it")


def test_the_two_instruments_are_not_conflated_in_the_headline():
    """The capture spread is computed by the MeSH map, not the keyword tagger."""
    md = MD.read_text()
    for mt in re.finditer(r"0\.20%-41\.86%", md):
        window = md[max(0, mt.start() - 300):mt.end() + 300]
        assert "NOT" in window or "not variation inside" in window, (
            "the capture spread is attached to the keyword reach without the "
            "correction; it is computed entirely by the MeSH instrument, "
            "which is the conflation this script's docstring forbids")


def test_the_excluded_tags_are_split_and_disclosed():
    """'demographic check-tags excluded' also removed study-design descriptors."""
    m = _mod()
    assert m.DEMOGRAPHIC_TAGS and m.STUDY_DESIGN_TAGS
    assert not (m.DEMOGRAPHIC_TAGS & m.STUDY_DESIGN_TAGS), (
        "a tag is in both sets, so the split does not partition")
    for probe in ("humans", "female", "mice"):
        assert probe in m.DEMOGRAPHIC_TAGS
    for probe in ("retrospective studies", "prognosis", "treatment outcome",
                  "cell line, tumor"):
        assert probe in m.STUDY_DESIGN_TAGS, (
            f"`{probe}` is a study-design descriptor and must not be excluded "
            "as a demographic check-tag")
    d = _doc()
    assert d.get("untagged_top_study_design") is not None, (
        "the study-design descriptors are excluded and no longer reported, so "
        "the section still silently drops what answers its own question")


def test_scan_returns_ordered_pairs_not_dicts():
    """Checks the SCAN's output shape, not the committed artifact.

    Reverting the storage to dicts is a scan-level change, so it is invisible
    to any guard that reads the committed JSON or runs `--render-only` -- and
    the storage shape is the whole reason the ordering is safe. Run over a
    handful of shards so this costs seconds rather than the full scan.
    """
    m = _mod()
    shards = sorted((m.ATLAS / "records").glob("*.jsonl.gz"))
    assert shards, "no census shards available"

    real_glob = type(m.ATLAS).glob
    subset = shards[::400][:4] or shards[:2]

    class _Stub:
        def __truediv__(self, other):
            return self

        def glob(self, _pat):
            return iter(subset)

    saved = m.ATLAS
    try:
        m.ATLAS = _Stub()
        out = m.scan(sample_every=1)
    finally:
        m.ATLAS = saved

    for key in ("per_mechanism", "untagged_top_mesh", "untagged_pubtypes"):
        assert isinstance(out[key], list), (
            f"scan() returns {key} as a {type(out[key]).__name__}; a dict is "
            "reordered by `json.dumps(sort_keys=True)`, which is what made the "
            "committed artifact alphabetical")
        assert all(isinstance(x, list) and len(x) == 2 for x in out[key])
        vals = [v for _k, v in out[key]]
        assert vals == sorted(vals, reverse=True), (
            f"scan() returns {key} out of descending order")
