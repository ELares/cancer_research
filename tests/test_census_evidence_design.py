"""Guards for the census design instrument that replaces the evidence ladder.

WHAT MAKES THIS EASY TO GET WRONG
---------------------------------
1. QUOTING THE WRONG DENOMINATOR. 44.5% of the census carries nothing that
   discriminates design. A distribution over the classifiable remainder,
   presented as a distribution over the literature, is the defect this repo has
   produced repeatedly. Both denominators must ship.

2. A CLASS NAMED FOR A TIER IT DOES NOT MEASURE. `animal-model` and
   `cell-culture` are MeSH check-tag inferences -- what organism or material
   appears in the study, not what kind of study it is. If the page ever calls
   them preclinical tiers it is claiming more than the label carries.

3. PRECEDENCE DRIFT. A trial that used a mouse model is a trial. If the class
   order changed silently, every share moves and nothing else would notice.

4. AN ACCURACY CLAIM. There is no tagger of ours here, so no accuracy figure is
   owed. If one appears, something has started scoring NLM against itself.
"""

import gzip
import importlib.util
import json
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPT = REPO_ROOT / "scripts" / "census_evidence_design.py"
MD = REPO_ROOT / "analysis" / "census-evidence-design.md"
JSON_OUT = REPO_ROOT / "analysis" / "census-evidence-design.json"


def _doc():
    return json.loads(JSON_OUT.read_text())


def _mod():
    spec = importlib.util.spec_from_file_location("ced", SCRIPT)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def test_the_classes_partition_the_census():
    d = _doc()
    assert sum(d["classes"].values()) == d["census"], (
        "the classes do not sum to the census, so a record is counted twice "
        "or dropped")
    assert d["classifiable"] == d["census"] - d["classes"]["undetermined"]
    assert d["census"] > 4_000_000


def test_both_denominators_are_published():
    """A share of the classifiable remainder is not a share of the literature."""
    d, md = _doc(), MD.read_text()
    und = d["classes"]["undetermined"]
    assert f"{und:,} of {d['census']:,} ({100*und/d['census']:.1f}%)" in md, (
        "the undetermined share is not stated against the census")
    assert "Coverage first" in md
    for k, v in d["classes"].items():
        assert f"{100*v/d['census']:.1f}%" in md, (
            f"`{k}`'s share of the CENSUS is not rendered")
        if k != "undetermined":
            assert f"{100*v/d['classifiable']:.1f}%" in md, (
                f"`{k}`'s share of the CLASSIFIABLE set is not rendered")


def test_the_classifier_is_a_pure_function_of_its_two_inputs():
    """Precedence is the judgement this file makes, so it is pinned."""
    m = _mod()
    # a trial that used a mouse model is a trial
    assert m.classify(["Clinical Trial, Phase III"],
                      ["Mice, Nude", "Xenograft Model Antitumor Assays"]) == "trial"
    # a review that discusses trials is not a trial
    assert m.classify(["Review", "Journal Article"], ["Humans"]) == "non-primary"
    # bare Journal Article with no discriminating descriptor
    assert m.classify(["Journal Article"], ["Humans", "Female"]) == "undetermined"
    # check-tag inferences, in order
    assert m.classify(["Journal Article"], ["Disease Models, Animal"]) == "animal-model"
    assert m.classify(["Journal Article"], ["Cell Line, Tumor"]) == "cell-culture"
    assert m.classify(["Journal Article"], ["Animals", "Rats"]) == "animal-other"
    # an animal study that also cultured cells is filed by the model, not the dish
    assert m.classify(["Journal Article"],
                      ["Mice, Nude", "Cell Line, Tumor"]) == "animal-model"
    # a human study with a cell line is not an animal study
    assert m.classify(["Journal Article"], ["Humans", "Animals"]) != "animal-other"
    # case-insensitive on descriptors, exact on publication types
    assert m.classify(["Journal Article"], ["CELL LINE, TUMOR"]) == "cell-culture"


def test_no_accuracy_claim_is_made_for_labels_we_did_not_assign():
    """There is no tagger of ours to score."""
    md = MD.read_text()
    assert "no accuracy figure is owed" in md
    import re
    for m_ in re.finditer(r"\b(precision|recall|F1|accuracy)\b", md, re.I):
        w = md[max(0, m_.start() - 300):m_.end() + 300]
        assert ("no accuracy figure is owed" in w or "46% exact-label" in w), (
            f"the page quotes {m_.group(0)!r} for labels NLM assigned; the "
            "only accuracy figure allowed here is the superseded ladder's, "
            "cited as what this replaces")


def test_the_check_tag_classes_are_named_for_what_they_measure():
    """`animal-model` and `cell-culture` are descriptor inferences, not tiers."""
    m, md = _mod(), MD.read_text()
    for k in ("animal-model", "cell-culture"):
        assert "descriptor" in m.WHAT_IT_MEASURES[k] or \
            "check tag" in m.WHAT_IT_MEASURES[k] or \
            "carries a" in m.WHAT_IT_MEASURES[k], (
            f"`{k}`'s stated meaning does not say it is a descriptor inference")
    assert "not what kind of study it is" in md
    assert "It is not an evidence hierarchy" in md
    # the withdrawn framing must not come back
    for bad in ("preclinical in vivo tier", "evidence tier", "seven-tier"):
        assert bad not in md.lower().replace("seven-tier evidence ladder", ""), (
            f"the page reintroduces {bad!r}, which is the taxonomy this "
            "instrument replaces rather than reproduces")


def test_the_excluded_stream_is_named():
    """Name the excluded stream without inventing its size for a fresh scan."""
    d, md = _doc(), MD.read_text()
    assert "Text-recovered census records are excluded" in md
    assert f"{d['census']:,} indexed records read by this analysis" in md, (
        "the page does not say which census stream its denominator is")
    assert "stream is not counted here" in md
    assert "783,271" not in md


def test_the_phase_column_is_real_and_derived():
    d, md = _doc(), MD.read_text()
    ph = d.get("phased_trials") or {}
    assert ph, "no phased trials were counted"
    assert sum(ph.values()) <= d["classes"]["trial"] * 2, (
        "more phase labels than trials could plausibly carry")
    for k, v in ph.items():
        assert f"| {k} | {v:,} |" in md
    assert d["classes"]["trial"] > 0


def test_the_committed_report_is_what_the_generator_produces():
    m = _mod()
    assert m.render(_doc()) == MD.read_text(), (
        "analysis/census-evidence-design.md is not what the renderer produces "
        "from the committed JSON -- re-run with --render-only")


def _scan_records(rows, tmp_path, monkeypatch):
    m = _mod()
    records = tmp_path / "records"
    records.mkdir()
    # Multiple shards exercise the default full scan rather than a prefix.
    for i, row in enumerate(rows):
        with gzip.open(records / f"part-{i:03d}.jsonl.gz", "wt",
                       encoding="utf-8") as stream:
            stream.write(json.dumps(row) + "\n")
    monkeypatch.setattr(m, "RECORDS", records)
    return m, m.scan()


@pytest.mark.parametrize("row,expected", [
    ({"pub_types": ["Review"]}, "non-primary"),
    ({"pub_types": ["Journal Article"], "mesh": ["Humans"]}, "undetermined"),
])
def test_readable_input_without_trials_renders(row, expected, tmp_path, monkeypatch):
    m, d = _scan_records([row], tmp_path, monkeypatch)

    assert d["census"] == 1
    assert d["classes"][expected] == 1
    assert d["classes"]["trial"] == 0
    assert d["phased_trials"] == {}
    md = m.render(d)
    assert "The trial column is the strong one" not in md
    assert "783,272" not in md
    if expected == "undetermined":
        assert d["classifiable"] == 0
        assert "percentages of the classifiable set are unavailable (N/A)" in md
        assert "1 of 1 (100.0%)" in md
        assert md.count("| N/A |") == 6
    else:
        assert d["classifiable"] == 1
        assert "N/A" not in md


def test_the_scan_counts_every_record_and_preserves_phase_and_year_rules(
        tmp_path, monkeypatch):
    """Run the real reader and scan on every design class, including unknowns."""
    _, d = _scan_records([
        {"pub_types": ["Clinical Trial, Phase I", "Clinical Trial, Phase II",
                       "Clinical Trial, Phase I", "Review"],
         "mesh": ["Mice, Nude"], "year": 1975},
        {"pub_types": ["Case Reports"], "year": 2026},
        {"mesh": ["Disease Models, Animal", "Cell Line, Tumor"], "year": 2000},
        {"mesh": ["Cell Line, Tumor"], "year": 1974},
        {"mesh": ["Rats"], "year": 2027},
        {"pub_types": ["Review"], "year": "2025"},
        {"pub_types": ["Journal Article"], "mesh": ["Humans"]},
        {"pub_types": [], "year": 2000},
    ], tmp_path, monkeypatch)

    assert d["census"] == 8
    assert d["classifiable"] == 6
    assert d["classes"] == {
        "trial": 1, "clinical-other": 1, "animal-model": 1,
        "cell-culture": 1, "animal-other": 1, "non-primary": 1,
        "undetermined": 2,
    }
    assert d["phased_trials"] == {
        "Clinical Trial, Phase I": 1, "Clinical Trial, Phase II": 1,
    }
    assert d["trial_share_by_year"] == {
        "1975": [1, 1], "2000": [0, 2], "2026": [0, 1],
    }


def test_bare_count_preserves_its_original_definition_and_describes_it(
        tmp_path, monkeypatch):
    m, d = _scan_records([
        {},
        {"pub_types": ["Journal Article", "English Abstract", "Multicenter Study",
                       "Comparative Study", "Evaluation Study"]},
        {"pub_types": ["Journal Article", "Research Support, Non-U.S. Gov't"]},
    ], tmp_path, monkeypatch)

    assert d["bare_or_uninformative_pub_types"] == 2
    md = m.render(d)
    assert "2 (66.7%) carry no publication types or only these nonspecific types" in md
    assert "funding" not in md
    for pub_type in m.UNINFORMATIVE:
        assert f"`{pub_type}`" in md
