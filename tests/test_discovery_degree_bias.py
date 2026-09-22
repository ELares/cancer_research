"""What the discovery rankers select, pinned against what they score.

`atlas-discovery-eval.md` observed that the harder a ranking corrects for
candidate degree, the worse it does -- ordering the methods by a hand-written
column describing what each formula does. `atlas-discovery-degree-bias.md`
measures that instead: `L` is the mean degree of a ranker's top-k over the mean
degree of the pool it drew them from, so `L = 1` is degree-neutral and
popularity is the ceiling by construction.

The measured ordering confirms the hand-written one and adds a magnitude. These
guards pin the relationship, not the numbers, and pin the limits with it --
because the interesting reading (the metric rewards not correcting) is one step
from an overclaim (hub-selection is wrong), and that step is not identifiable
from this corpus.
"""

import json
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
RAW = REPO_ROOT / "analysis" / "atlas-discovery-degree-bias.json"
DOC = REPO_ROOT / "analysis" / "atlas-discovery-degree-bias.md"
EVAL = REPO_ROOT / "analysis" / "atlas-discovery-eval.json"


def _raw():
    return json.loads(RAW.read_text())


def test_the_measure_is_calibrated_on_its_own_control():
    """Pin the retained random control's observed proximity to degree-neutral.

    Uniform sampling has expected mean L = 1 within each pool. The median over
    finite samples need not be 1; this tolerance describes the retained study,
    rather than imposing a universal requirement on alternative snapshots.
    """
    L = _raw()["median_L"]
    assert "random" in L, "the control ranking is missing"
    assert 0.8 <= L["random"] <= 1.25, (
        f"random sits at {L['random']:.2f}x, outside the retained control's range; "
        "revisit the snapshot and its interpretation")


def test_popularity_is_the_ceiling_and_jaccard_is_below_neutral():
    """Pin popularity's structural ceiling and Jaccard's observed location."""
    L = _raw()["median_L"]
    assert L["popularity"] == max(L.values()), (
        "popularity is degree, so nothing can be more degree-selective")
    assert L["jaccard"] < 1.0, (
        "Jaccard is no longer below neutral in the retained study; revisit "
        "the snapshot and its interpretation")


def test_selectivity_predicts_precision():
    """The relationship the document is about."""
    d = _raw()
    rho = d["spearman_L_vs_precision_excluding_random"]
    assert rho > 0.8, (
        f"L no longer predicts precision (rho={rho:.2f}); the document's central "
        "claim does not hold and its reading must be revisited")


def test_the_two_documents_describe_the_same_build():
    """L and precision are quoted side by side, so they must share a split.

    Quoting a selectivity computed on one graph beside a precision computed on
    another is the error that let a sibling analysis in this repo invert its own
    finding. Cheap to check, so checked.
    """
    d = _raw()
    ev = json.loads(EVAL.read_text())["headline"]
    assert d["split"] == ev["split_year"], (
        f"degree bias is computed at split {d['split']} but the precisions it "
        f"quotes come from {ev['split_year']}")
    assert d["top"] == ev["top_k"], "different k"
    # Split year and k are NOT a build fingerprint. Both matched across two
    # graphs that differed by 11% in edge count -- the corrections file this
    # analysis consumes had grown -- and this guard passed while the document
    # paired selectivity from one build with precision from the other.
    assert d.get("pairs_before") == ev["pairs_before"], (
        f"different graph builds: degree bias saw {d.get('pairs_before')} pairs "
        f"before the split, the evaluation saw {ev['pairs_before']}. Every L "
        "quoted beside a precision describes a different world.")
    for m, p in d["precision_at_k"].items():
        assert abs(p - ev["precision"][m]) < 1e-9, (
            f"the precision quoted for {m} is not the evaluation's")


def _flat() -> str:
    """The document with whitespace collapsed.

    The generator wraps prose, so a phrase can arrive split across a newline --
    `"not identifiable"` is stored as `"not\nidentifiable"`. Matching raw text
    makes a guard fail for a reason that has nothing to do with what it checks,
    and (worse) makes it PASS if a later rewrap happens to join the words.
    """
    return " ".join(DOC.read_text().split()).lower()


def test_the_document_does_not_claim_hub_selection_is_wrong():
    """One step past the finding is an unidentifiable claim."""
    txt = _flat()
    assert "does not show hub-selection is wrong" in txt, (
        "the document draws the identifiability limit nowhere")
    assert "not identifiable" in txt, (
        "the document does not say the separation is unidentifiable, so a "
        "reader may take the correlation as evidence that correcting for "
        "degree is the right thing to do")


def test_the_document_credits_the_prior_observation():
    """The direction was already in the evaluation; only the magnitude is new."""
    txt = DOC.read_text()
    assert "already suspected this" in txt, (
        "the document presents a measurement of an existing observation as a "
        "new observation")


def test_degree_snapshot_renders_offline_without_replacing_json(monkeypatch, tmp_path):
    """Prose repair must not silently rebuild a different graph or data snapshot."""
    import importlib
    import sys
    sys.path.insert(0, str(REPO_ROOT / "scripts"))
    mod = importlib.import_module("atlas_discovery_degree_bias")
    source = RAW.read_bytes()
    raw = tmp_path / "snapshot.json"
    md = tmp_path / "report.md"
    raw.write_bytes(source)
    md.write_text("old report")
    monkeypatch.setattr(mod, "OUT_JSON", raw)
    monkeypatch.setattr(mod, "OUT_MD", md)

    def unavailable(*args, **kwargs):
        raise AssertionError("offline replay accessed graph inputs")

    for name in ("atlas_root", "load_index", "pmid_years", "pair_first_year", "load_corrections"):
        monkeypatch.setattr(mod, name, unavailable)
    assert mod.main(["--render-only"]) == 0
    assert raw.read_bytes() == source
    assert md.read_text() == mod.render(json.loads(source)) == DOC.read_text()
    assert "candidate-selection control" in md.read_text()


def test_degree_render_failure_preserves_both_reports(monkeypatch, tmp_path):
    import importlib
    import sys
    import pytest
    sys.path.insert(0, str(REPO_ROOT / "scripts"))
    mod = importlib.import_module("atlas_discovery_degree_bias")
    source = RAW.read_bytes()
    raw, md = tmp_path / "snapshot.json", tmp_path / "report.md"
    raw.write_bytes(source)
    md.write_text("old report")
    monkeypatch.setattr(mod, "OUT_JSON", raw)
    monkeypatch.setattr(mod, "OUT_MD", md)

    def broken(_):
        raise ValueError("render failure")

    monkeypatch.setattr(mod, "render", broken)
    with pytest.raises(ValueError, match="render failure"):
        mod.main(["--render-only"])
    assert raw.read_bytes() == source
    assert md.read_text() == "old report"


def test_degree_snapshot_rejects_invalid_metrics_before_writing(monkeypatch, tmp_path):
    import importlib
    import sys
    import pytest
    sys.path.insert(0, str(REPO_ROOT / "scripts"))
    mod = importlib.import_module("atlas_discovery_degree_bias")
    raw, md = tmp_path / "snapshot.json", tmp_path / "report.md"
    monkeypatch.setattr(mod, "OUT_JSON", raw)
    monkeypatch.setattr(mod, "OUT_MD", md)
    for medians in ({}, {**_raw()["median_L"], "abc": float("nan")},
                    {**_raw()["median_L"], "abc": -1}):
        data = _raw()
        data["median_L"] = medians
        raw.write_text(json.dumps(data))
        source = raw.read_bytes()
        md.write_text("old report")
        with pytest.raises(ValueError):
            mod.main(["--render-only"])
        assert raw.read_bytes() == source
        assert md.read_text() == "old report"


def _render_alternative(medians, precision=None):
    """Render internally consistent alternative summaries without graph access."""
    import importlib
    import sys
    sys.path.insert(0, str(REPO_ROOT / "scripts"))
    mod = importlib.import_module("atlas_discovery_degree_bias")
    data = _raw()
    data["median_L"] = medians
    data["share_of_popularity"] = {
        m: (v - 1) / (medians["popularity"] - 1)
        for m, v in medians.items() if medians["popularity"] > 1
    }
    if precision is not None:
        data["precision_at_k"] = precision
    paired = [m for m in medians if m in data["precision_at_k"]]
    for key, methods in (
        ("spearman_L_vs_precision", paired),
        ("spearman_L_vs_precision_excluding_random", [m for m in paired if m != "random"]),
    ):
        data[key] = mod._spearman(
            [medians[m] for m in methods], [data["precision_at_k"][m] for m in methods]
        )
    return " ".join(mod.render(data).split())


@pytest.mark.parametrize("jaccard, relation", [(0.5, "below"), (1, "equal to"), (1.5, "above")])
def test_degree_interpretation_follows_observed_medians(jaccard, relation):
    medians = {**_raw()["median_L"], "abc": 7.356947036860076, "jaccard": jaccard}
    text = _render_alternative(medians)
    assert "ABC's median L is 7.36x, equal to popularity's 7.36x" in text
    assert "ABC retains 100%" in text
    assert f"Jaccard's median L is {jaccard:.2f}x, {relation} the degree-neutral value" in text
    assert "correction is real and it is large" not in text


def test_degree_neutral_popularity_has_no_relative_selectivity():
    text = _render_alternative(dict.fromkeys(_raw()["median_L"], 1.0))
    assert "share of popularity's excess selectivity is unavailable" in text
    assert "ABC retains 0%" not in text
    assert "| abc | **1.0x** | -- |" in text
    assert "unavailable (fewer than three methods or constant values)" in text


@pytest.mark.parametrize("direction", [1, -1])
def test_degree_correlation_does_not_determine_random_selectivity(direction):
    medians = _raw()["median_L"]
    # Both perfect positive and negative orderings can have random's median != 1.
    precision = {m: (v / 10 if direction == 1 else 1 - v / 10) for m, v in medians.items()}
    text = _render_alternative(medians, precision)
    assert f"Rank correlation between L and precision is **{direction:.2f}**" in text
    assert "over 7 methods with both measurements" in text
    assert "over 6 methods excluding random" in text
    assert "observed median L for `random` is 0.91x" in text
    assert "sits exactly" not in text
    assert "near-perfect ordering" not in text
    assert "expected mean L = 1" in text


@pytest.mark.parametrize("precision", [{}, {"abc": 0.1, "popularity": 0.2}])
def test_degree_missing_precision_is_not_reported_as_zero_correlation(precision):
    text = _render_alternative(_raw()["median_L"], precision)
    assert "Rank correlation between L and precision is unavailable" in text
    assert f"over {len(precision)} methods with both measurements" in text
    assert "Rank correlation between L and precision is **0.00**" not in text
