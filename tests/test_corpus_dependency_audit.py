"""Guards for the corpus-dependency audit.

This artifact could be used to justify DELETING an analysis, so its two failure
modes are not symmetric and both are pinned. Both actually happened while it
was being written, one after the other, in opposite directions:

  1. A substring scan counted files whose only mention of the corpus is prose
     ABOUT it -- including `atlas_baseline.py`, whose sole mention is a
     docstring promising it leaves the corpus alone, and the audit itself,
     which names every marker as a constant.
  2. Narrowing the scan to code created the mirror error and the dangerous
     one: most of the pipeline imports its paths from `config` and never
     writes a literal, so `tag_articles.py`, `fetch_articles.py` and five more
     -- the scripts that BUILD the corpus -- were dropped from the table
     entirely. A false negative here reads as "nothing uses this".

So the guards check a KNOWN consumer stays in and a KNOWN non-consumer stays
out, rather than checking a count.
"""
import json
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
JSON = REPO / "analysis/corpus-dependency-audit.json"
MD = REPO / "analysis/corpus-dependency-audit.md"

# Reads the corpus only through paths imported from config -- the case the
# code-only scan dropped. If this leaves the table, the second failure mode
# has returned.
KNOWN_CONSUMERS_VIA_CONFIG = {"tag_articles.py", "fetch_articles.py"}
# Mentions a corpus path exclusively in a docstring stating it does NOT touch
# the corpus. If this enters the table, the first failure mode has returned.
KNOWN_NON_CONSUMER = "atlas_baseline.py"


@pytest.fixture(scope="module")
def d():
    return json.loads(JSON.read_text())


def test_a_consumer_reached_only_through_config_paths_is_counted(d):
    names = {r["script"] for r in d["consumers"]}
    missing = KNOWN_CONSUMERS_VIA_CONFIG - names
    assert not missing, (
        f"{sorted(missing)} read the corpus through paths imported from "
        "config and are absent from the audit; a code-only string scan drops "
        "exactly these, and their absence would read as nothing using them")
    for r in d["consumers"]:
        if r["script"] in KNOWN_CONSUMERS_VIA_CONFIG:
            assert "config path symbol" in r["route"]


def test_prose_about_the_corpus_is_not_counted_as_use(d):
    names = {r["script"] for r in d["consumers"]}
    assert KNOWN_NON_CONSUMER not in names, (
        f"{KNOWN_NON_CONSUMER}'s only corpus mention is a docstring saying it "
        "leaves the corpus alone; counting it means the detector is reading "
        "prose as use")
    assert KNOWN_NON_CONSUMER in d["excluded_mentions_only"]
    assert "corpus_dependency_audit.py" not in names, (
        "the audit classified itself, which it can only do by matching its own "
        "marker constants")


def test_the_audit_names_both_failure_modes_it_hit(d):
    """A correction that leaves no trace invites the same fix twice."""
    md = MD.read_text()
    assert "Named the corpus without reading it" in md
    src = (REPO / "scripts/corpus_dependency_audit.py").read_text()
    assert "created a false-negative one" in src, (
        "the second failure mode is not recorded, so a future narrowing of the "
        "scan will look like an improvement")


def test_classes_are_derived_from_the_fields_each_script_reads(d):
    """The class must follow from the recorded evidence, not stand beside it."""
    for r in d["consumers"]:
        if r["annotations"]:
            assert r["class"] == "needs project annotations", r["script"]
        elif r["uses_fulltext"]:
            assert r["class"] == "full text, available at census scale", r["script"]
        elif r["census_fields_read"]:
            assert r["class"] == "census-supplied fields only", r["script"]
        else:
            assert r["class"] == "no record fields resolved", r["script"]
    counts = {}
    for r in d["consumers"]:
        counts[r["class"]] = counts.get(r["class"], 0) + 1
    assert counts == d["by_class"]
    assert d["n_consumers"] == len(d["consumers"])


def test_the_corpus_bound_set_is_exactly_the_annotation_readers(d):
    assert sorted(d["corpus_bound"]) == sorted(
        r["script"] for r in d["consumers"] if r["annotations"])


def test_it_issues_no_keep_or_cut_verdict():
    """A dependency measurement that also decided what to delete would be
    making an editorial call under cover of an arithmetic one."""
    md = MD.read_text()
    assert "What this does not decide" in md
    for verdict in ("should be deleted", "safe to remove", "can be cut",
                    "recommend removing"):
        assert verdict not in md.lower(), (
            f"the audit issues a verdict ({verdict!r}); it reports dependencies")


def test_the_annotation_list_is_not_silently_narrowed(d):
    """Dropping a field from PROJECT_ANNOTATIONS would move consumers out of
    the corpus-bound class without any measurement changing."""
    for required in ("mechanisms", "cancer_types", "evidence_level",
                     "pathway_targets", "diagnostic_therapy_links"):
        assert required in d["project_annotations"], (
            f"`{required}` is a tag only this project produces and is no "
            "longer counted as a corpus dependency")
