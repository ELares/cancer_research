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
    assert "from tag_articles import text_matches_keyword" in src, (
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
    """A sampled rate reported as a point invites over-reading."""
    m, d = _mod(), _doc()
    lo, hi = m.wilson(d["keyword_hits"], d["sampled"])
    assert lo < d["keyword_hits"] / d["sampled"] < hi
    md = MD.read_text()
    assert "95% CI" in md, "the keyword reach is reported without its interval"
    # the full-census figure must NOT be given a sampling interval
    assert md.count("95% CI") == 1, (
        "more than one interval is reported; the MeSH-map figure is a full "
        "census count and has no sampling error")


def test_the_descriptor_profile_excludes_check_tags():
    """`humans 92%` is true and says nothing; the exclusion is the analysis."""
    m, d = _mod(), _doc()
    assert "humans" in m.CHECK_TAGS and "female" in m.CHECK_TAGS
    top = list(d["untagged_top_mesh"])
    assert top, "no descriptor profile was produced"
    leaked = [t for t in top if t in m.CHECK_TAGS]
    assert not leaked, (
        f"check-tags {leaked} are in the top-descriptor table, so it is "
        "reporting demographics rather than subject matter")
    # AND the code must still exclude them. The assertion above reads the
    # committed artifact, which does not move when the scan is edited, so on
    # its own it would pass for a script that had stopped filtering.
    assert "if m not in CHECK_TAGS:" in SCRIPT.read_text(), (
        "the scan no longer excludes check-tags when building the profile")
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
    ther = sum(v for k, v in d["untagged_top_mesh"].items()
               if "antineoplastic" in k or "chemotherapy" in k
               or "radiotherapy" in k or "drug therapy" in k)
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
