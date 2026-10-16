"""Guards for "is that a statement about MeSH, or about our map?" (#952).

EVERY GUARD HERE DRIVES `assemble()` ON A DOCTORED FIXTURE. None reads the
committed artifact to check a derived field, because a guard that reads only a
generated artifact CANNOT FAIL when its generator changes -- the artifact is
committed and does not move until somebody re-runs it. That defect has been
found in this repository seven times across three pull requests, four of them
in one sweep, so the rule here is absolute: the fixture is the input and
`assemble` is the thing under test.

THE FIXTURE IS A VOCABULARY, which is what makes the arms drivable. Both arms
are rules over descriptor NAMES applied to the census's descriptor vocabulary,
so a doctored vocabulary can make a descriptor appear, disappear, shrink below
the pool floor, or match inside a token -- and each of those is a way one of the
three verdicts can be wrong.
"""
import ast
import json
import math
import re
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts"))

import descriptor_absence_audit as daa  # noqa: E402
from config import MECHANISM_KEYWORDS  # noqa: E402

# Captured BEFORE the autouse fixture below can replace it, so the one guard
# that is about the repository rather than about a fixture can get the real
# reader back. Without this that guard asserts against the stub and is
# vacuous -- which is exactly what it was, until a read of this file caught it.
REAL_SITE_PRESENT = daa._site_present

SCRIPT = REPO / "scripts" / "descriptor_absence_audit.py"
SRC = SCRIPT.read_text()


# --------------------------------------------------------------------------
# Fixture: a small vocabulary and a map over it.
# --------------------------------------------------------------------------
def _vocab(rows, scanned=4_403_994):
    return {"scanned": scanned, "rows": [dict(r) for r in rows]}


def _row(ui, name, records, major=0):
    return {"ui": ui, "name": name, "records": records, "major": major}


# A descriptor whose name no mechanism keyword reaches, so a fixture can hold a
# vocabulary without accidentally satisfying an arm.
NEUTRAL = _row("D000001", "Neoplasms", 900_000, 400_000)

BASE_ROWS = [
    NEUTRAL,
    _row("D004599", "Electric Stimulation Therapy", 945, 434),
    _row("D058626", "Plasma Gases", 474, 182),
    _row("D003300", "Copper", 6_526, 1_333),
    _row("D000098560", "Disulfidptosis", 99, 1),
    _row("D000094243", "Proteolysis Targeting Chimera", 241, 49),
    _row("D059748", "Proteolysis", 5_288, 814),
    _row("D004572", "Electrolysis", 85, 27),
    _row("D053672", "Electrochemotherapy", 744, 360),
    _row("D019275", "Radiopharmaceuticals", 31_784, 8_217),
    _row("D016499", "Radioimmunotherapy", 3_240, 1_131),
    _row("D008564", "Membrane Potentials", 3_587, 163),
    _row("D007473", "Ion Channels", 1_624, 125),
    _row("D004574", "Electromagnetic Fields", 2_132, 694),
    _row("D055909", "Magnetic Field Therapy", 214, 71),
]


def _entry(verdict, covering=(), adjacent=(), leaves=None, reason="because."):
    out = {"verdict": verdict, "covering": list(covering),
           "adjacent": list(adjacent), "reason": reason}
    if leaves is not None:
        out["leaves"] = list(leaves)
    return out


def _map():
    """The real map's shape, written out so a guard can doctor one entry."""
    return {"unmeasurable": {
        "ttfields": _entry("absent", adjacent=["Electric Stimulation Therapy"]),
        "cold-atmospheric-plasma": _entry("over-broad", covering=["Plasma Gases"]),
        "bioelectric": _entry("absent",
                              adjacent=["Membrane Potentials", "Ion Channels"]),
        "frequency-therapy": _entry(
            "absent", adjacent=["Electromagnetic Fields", "Magnetic Field Therapy"]),
        "cuproptosis": _entry("over-broad", covering=["Copper"]),
        "disulfidptosis": _entry("name-exact-leaky", leaves=["Disulfidptosis"]),
        "electrolysis": _entry("name-exact-leaky", adjacent=["Electrochemotherapy"],
                               leaves=["Electrolysis"]),
        "radioligand-therapy": _entry("over-broad",
                                      covering=["Radiopharmaceuticals"],
                                      adjacent=["Radioimmunotherapy"]),
        "targeted-protein-degradation": _entry(
            "name-exact-leaky", covering=["Proteolysis"],
            leaves=["Proteolysis Targeting Chimera"]),
    }}


def _detail(rows, mechanisms, overlap=None):
    """Raw pass-2 counts for every descriptor the fixture holds."""
    return {
        "scanned": 4_403_994,
        "descriptor_records": {r["ui"]: r["records"] for r in rows},
        "text_arm_records": {m: 500 for m in mechanisms},
        "both": dict(overlap or {}),
        "by_year": {r["ui"]: {"2024": r["records"]} for r in rows},
        "pub_types": {r["ui"]: {"Journal Article": r["records"]} for r in rows},
    }


def _raw(rows=None, mapping=None, overlap=None):
    rows = BASE_ROWS if rows is None else rows
    mapping = _map() if mapping is None else mapping
    return {"vocabulary": _vocab(rows),
            "detail": _detail(rows, sorted(mapping["unmeasurable"]), overlap),
            "map": mapping}


@pytest.fixture(autouse=True)
def _sites_all_present(monkeypatch):
    """Registry phrases are a property of the repository, not of the fixture.

    Every guard below is about the descriptor measurement, so the site reader is
    stubbed to the state the repository is in after this pull request: the
    replacement wording present, the retired wording gone. The guards that are
    ABOUT the registry override this deliberately.
    """
    def present(path, phrase):
        retired = {r for _, _, r, _ in daa.carrying_sites() if r}
        return phrase not in retired
    monkeypatch.setattr(daa, "_site_present", present)


def _assemble(**kw):
    return daa.assemble(_raw(**kw))


# --------------------------------------------------------------------------
# The arms are re-derived, not echoed.
# --------------------------------------------------------------------------
def test_the_round_trip_rereads_its_inputs_rather_than_echoing_them():
    """`_split_stored` must call `load`, or every guard here is inert.

    A value assertion cannot tell an echo from a recompute when the two agree,
    which is why this is structural -- the lesson `acsl4_slope_bound.py` paid
    for with a sweep that survived four behavioural guards.
    """
    tree = ast.parse(SRC)
    fn = next(n for n in tree.body
              if isinstance(n, ast.FunctionDef) and n.name == "_split_stored")
    called = {n.func.id for n in ast.walk(fn)
              if isinstance(n, ast.Call) and isinstance(n.func, ast.Name)}
    assert "load" in called, "_split_stored must re-read the committed inputs"


def test_assemble_rederives_both_arms_from_the_vocabulary():
    """`assemble` must call `candidate_uis`, which is where the arms live."""
    tree = ast.parse(SRC)
    fn = next(n for n in tree.body
              if isinstance(n, ast.FunctionDef) and n.name == "assemble")
    called = {n.func.id for n in ast.walk(fn)
              if isinstance(n, ast.Call) and isinstance(n.func, ast.Name)}
    assert "candidate_uis" in called, (
        "assemble must re-run the arms over the committed vocabulary; reading a "
        "stored verdict would make every guard below inert")


def test_a_descriptor_appearing_in_the_vocabulary_refutes_an_absent_verdict():
    """The whole finding, driven: a leaf refutes `absent`."""
    clean = _assemble()
    assert clean["mechanisms"]["ttfields"]["checks"][
        "no_name_exact_where_absence_claimed"]
    assert clean["mechanisms"]["ttfields"]["n_arm_a"] == 0

    doctored = BASE_ROWS + [_row("D999001", "Tumor Treating Fields Therapy", 300)]
    out = _assemble(rows=doctored)
    mech = out["mechanisms"]["ttfields"]
    assert mech["n_arm_a"] == 1, "arm A must reach a name-exact descriptor"
    assert not mech["checks"]["no_name_exact_where_absence_claimed"]
    assert not out["all_checks_pass"]
    assert "Tumor Treating Fields Therapy" in daa.render(out)


def test_a_descriptor_appearing_in_the_vocabulary_refutes_an_over_broad_verdict():
    """`over-broad` is refused too, so breadth cannot stand in for absence."""
    doctored = BASE_ROWS + [_row("D999002", "Cold Atmospheric Plasma", 60)]
    out = _assemble(rows=doctored)
    mech = out["mechanisms"]["cold-atmospheric-plasma"]
    assert mech["n_arm_a"] == 1
    assert not mech["checks"]["no_name_exact_where_absence_claimed"]
    assert not out["all_checks_pass"]


def test_a_named_descriptor_that_matches_nothing_is_caught():
    """`atlas_thesis_rank.py` had 5 of 17 of its own stems in this state."""
    rows = [r for r in BASE_ROWS if r["name"] != "Plasma Gases"]
    out = _assemble(rows=rows)
    mech = out["mechanisms"]["cold-atmospheric-plasma"]
    assert not mech["checks"]["every_named_descriptor_exists"]
    assert mech["named"]["Plasma Gases"]["exists"] is False
    assert "**NO**" in daa.render(out)


def test_absent_may_not_name_a_descriptor_that_contains_the_mechanism():
    """The difference between `absent` and `over-broad` is exactly this field."""
    mapping = _map()
    mapping["unmeasurable"]["ttfields"] = _entry(
        "absent", covering=["Electric Stimulation Therapy"])
    out = daa.assemble(_raw(mapping=mapping))
    assert not out["mechanisms"]["ttfields"]["checks"]["absent_names_no_covering"]
    assert not out["all_checks_pass"]


def test_over_broad_needs_a_cover_that_is_materially_wider():
    """"Over-broad" with a narrow cover is a verdict about nothing."""
    rows = [dict(r, records=40) if r["name"] == "Copper" else r for r in BASE_ROWS]
    out = _assemble(rows=rows)
    mech = out["mechanisms"]["cuproptosis"]
    assert not mech["checks"]["over_broad_has_a_wide_cover"]
    assert out["pool_floor"] * out["over_broad_multiple"] > 40


def test_a_leaf_the_scan_finds_must_be_declared_in_the_map():
    """A descriptor found here and unrecorded there leaves the next reader stuck."""
    mapping = _map()
    mapping["unmeasurable"]["disulfidptosis"] = _entry("name-exact-leaky",
                                                       leaves=[])
    out = daa.assemble(_raw(mapping=mapping))
    mech = out["mechanisms"]["disulfidptosis"]
    assert mech["n_arm_a"] == 1, "the leaf is still found"
    assert not mech["checks"]["leaves_declared"]


def test_name_exact_leaky_requires_the_scan_to_find_a_leaf():
    """The claim cannot be made where nothing names the mechanism."""
    mapping = _map()
    mapping["unmeasurable"]["ttfields"] = _entry("name-exact-leaky",
                                                 leaves=["Nothing At All"])
    out = daa.assemble(_raw(mapping=mapping))
    mech = out["mechanisms"]["ttfields"]
    assert not mech["checks"]["name_exact_has_a_leaf"]
    assert not mech["checks"]["every_named_descriptor_exists"]


def test_the_verdict_vocabulary_is_closed():
    mapping = _map()
    mapping["unmeasurable"]["ttfields"] = _entry("unmeasurable-ish")
    out = daa.assemble(_raw(mapping=mapping))
    assert not out["mechanisms"]["ttfields"]["checks"]["verdict_is_known"]


# --------------------------------------------------------------------------
# Volume-measurable is not recall-measurable, and a small pool is neither.
# --------------------------------------------------------------------------
def test_a_leaf_below_the_pool_floor_is_found_but_not_volume_measurable():
    rows = [dict(r, records=daa.POOL_FLOOR - 1) if r["name"] == "Disulfidptosis"
            else r for r in BASE_ROWS]
    out = _assemble(rows=rows)
    mech = out["mechanisms"]["disulfidptosis"]
    assert mech["n_arm_a"] == 1
    assert mech["volume_measurable"] is False
    assert "cannot" in daa.render(out)


def test_no_mechanism_here_is_recall_measurable():
    """The registered restraint: a leaky leaf cannot be a reference label."""
    out = _assemble()
    for mech, v in out["mechanisms"].items():
        assert v["recall_measurable"] is False, mech
        for row in v["arm_a"]:
            keywords = MECHANISM_KEYWORDS.get(mech) or []
            assert daa.arm_a(row["name"], keywords), (
                f"{row['name']} must trip {mech}'s own keywords, which is what "
                "makes it LEAKY by mechanism_recall's classifier")


def test_the_pool_floor_is_the_one_mechanism_recall_uses():
    """Two copies of a floor disagree; this one is imported."""
    from mechanism_recall import N_MIN
    assert daa.POOL_FLOOR == N_MIN
    assert "from mechanism_recall import N_MIN as POOL_FLOOR" in SRC, (
        "the floor must be imported, not restated")


# --------------------------------------------------------------------------
# Absence is a third state, not a zero.
# --------------------------------------------------------------------------
def test_a_missing_text_overlap_is_withheld_rather_than_published_as_zero():
    """A 0.0 under a column headed agreement reads as measured disagreement."""
    out = _assemble(overlap={})
    row = out["mechanisms"]["disulfidptosis"]["arm_a"][0]
    assert row["overlap"] is None
    assert row["descriptor_found_by_text"] is None
    assert row["text_found_by_descriptor"] is None
    assert "0.0%" not in daa.render(out)


def test_a_real_overlap_is_reported_both_ways():
    out = _assemble(overlap={"disulfidptosis|D000098560": 50})
    row = out["mechanisms"]["disulfidptosis"]["arm_a"][0]
    assert row["overlap"] == 50
    assert row["descriptor_found_by_text"] == pytest.approx(50 / 99)
    assert row["text_found_by_descriptor"] == pytest.approx(50 / 500)


_NOT_A_NUMBER = re.compile(r"\b(?:nan|inf|-inf)\b", re.I)


def test_the_artifact_never_carries_a_nan_or_an_infinity():
    """`json.dumps` writes NaN as the bare token NaN, which is not JSON."""
    out = _assemble(overlap={})

    def walk(node):
        if isinstance(node, float):
            assert math.isfinite(node), "non-finite value in the artifact"
        elif isinstance(node, dict):
            for v in node.values():
                walk(v)
        elif isinstance(node, list):
            for v in node:
                walk(v)
    walk(out)
    hit = _NOT_A_NUMBER.search(daa.render(out))
    assert hit is None, f"rendered page carries {hit.group(0)!r}"


def test_the_page_does_not_depend_on_the_stored_dicts_key_order():
    rows = list(reversed(BASE_ROWS))
    mapping = _map()
    mapping["unmeasurable"] = dict(reversed(list(mapping["unmeasurable"].items())))
    shuffled = daa.assemble(_raw(rows=rows, mapping=mapping))
    assert daa.render(shuffled) == daa.render(_assemble())


# --------------------------------------------------------------------------
# The two arms, and the limit of the wider one.
# --------------------------------------------------------------------------
def test_arm_a_is_the_production_matcher_and_not_the_loose_one():
    """A sweep caught this: pointing arm A at arm B survived every guard.

    The guard that was supposed to catch it asserted only that each found row
    trips `arm_a` -- which is trivially true once `arm_a` IS the rule that
    found them. What discriminates is a case the two rules DISAGREE on, and
    `mechanism_recall.py`'s leak classifier is the reason it has to be arm A:
    if the two rules differ, "LEAKY there" and "name-exact here" stop meaning
    the same thing.
    """
    assert not daa.arm_a("Ion Channels", ["ion channel cancer"]), (
        "the production matcher is a substring rule and must not match a "
        "token-subset case")
    assert daa.arm_b("Ion Channels", ["ion channel cancer"]), (
        "the looser arm must, or the pair is not a disagreement"
    )
    assert daa.arm_a("Disulfidptosis", ["disulfidptosis"])
    assert "text_matches_keyword(lowered, k)" in SRC, (
        "arm A must call the live production matcher")


def test_arm_b_is_looser_than_arm_a_on_the_fixture():
    out = _assemble()
    for mech, v in out["mechanisms"].items():
        assert v["n_arm_b"] >= v["n_arm_a"], mech


def test_arm_a_can_reach_inside_a_token_where_arm_b_cannot_and_the_page_says_so():
    """Arm A is a SUBSTRING rule, so the convenient subset direction is not free.

    Without this the page could claim "an absence under both arms" while arm A
    carried hits arm B never saw.
    """
    rows = BASE_ROWS + [_row("D999003", "Protacanthopterygii", 70)]
    out = _assemble(rows=rows)
    mech = out["mechanisms"]["targeted-protein-degradation"]
    assert "D999003" in mech["arm_a_not_in_arm_b"]
    assert "D999003" in out["arm_a_outside_arm_b"]
    assert "arm B does not" in daa.render(out)
    assert daa.arm_a("Protacanthopterygii", ["protac"])
    assert not daa.arm_b("Protacanthopterygii", ["protac"])


def test_arm_b_drops_generic_words_from_the_keyword_side():
    """What the generic list is FOR: a keyword's distinctive tokens decide.

    AND ITS ASYMMETRY IS INERT, measured rather than assumed. A sweep that also
    subtracted the generics from the DESCRIPTOR side survived every guard, and
    the reason is structural: the required tokens come from the keyword's
    distinctive set, which already excludes generics, so removing them from the
    name side can only make the subset test HARDER and never manufactures a
    match. The comment in the source claims exactly that, so this records the
    asymmetry as a no-op rather than guarding a defect that cannot occur.
    """
    assert daa.arm_b("Radionuclide Imaging", ["radionuclide therapy"])
    assert not daa.arm_b("Imaging", ["radionuclide therapy"])
    assert "therapy" in daa.ARM_B_GENERIC

    def arm_b_both_sides(name, keywords):
        tokens = set(daa._normalise(name)) - daa.ARM_B_GENERIC
        for keyword in keywords:
            distinctive = [t for t in daa._normalise(keyword)
                           if t not in daa.ARM_B_GENERIC]
            if distinctive and set(distinctive) <= tokens:
                return True
        return False

    vocab = daa.read_vocabulary() if daa.VOCAB_TSV_GZ.exists() else None
    if vocab is None:
        pytest.skip("descriptor vocabulary not committed")
    differ = []
    for mech in sorted(daa.unusable(daa._yaml())):
        keywords = MECHANISM_KEYWORDS.get(mech) or []
        for row in vocab["rows"]:
            if bool(daa.arm_b(row["name"], keywords)) != bool(
                    arm_b_both_sides(row["name"], keywords)):
                differ.append((mech, row["name"]))
    assert not differ, (
        "subtracting the generics on the descriptor side changed a verdict, so "
        "the asymmetry is load-bearing after all and needs its own guard: "
        f"{differ[:5]}")


# --------------------------------------------------------------------------
# The registry, in both directions.
# --------------------------------------------------------------------------
def test_every_registered_replacement_is_present_and_every_retired_wording_gone(
        monkeypatch):
    """Drives the real repository, which is what the registry is about.

    AND IT HAS TO UNDO THE AUTOUSE FIXTURE to do that. The fixture stubs the
    site reader to the state this pull request leaves the tree in, which is
    right for every guard about the descriptor measurement and makes THIS one
    assert against its own stub. Found by reading the file rather than by a
    sweep, because a vacuous guard passes.
    """
    monkeypatch.setattr(daa, "_site_present", REAL_SITE_PRESENT)
    if not (daa.VOCAB_TSV_GZ.exists() and daa.DETAIL_JSON.exists()):
        pytest.skip("census derivatives not committed")
    out = daa.assemble(daa.load())
    assert out["sites_missing"] == [], out["sites_missing"]
    assert out["sites_retired_still_present"] == [], (
        out["sites_retired_still_present"])
    assert out["sites_misclassified"] == [], out["sites_misclassified"]


def test_a_retired_wording_that_comes_back_is_caught(monkeypatch):
    monkeypatch.setattr(daa, "_site_present", lambda path, phrase: True)
    out = _assemble()
    assert out["sites_retired_still_present"], (
        "the flattened wording returning must be reported")
    assert "**NO**" in daa.render(out)


def test_a_replacement_wording_that_moves_is_caught(monkeypatch):
    monkeypatch.setattr(daa, "_site_present", lambda path, phrase: False)
    out = _assemble()
    assert out["sites_missing"], "a registry entry pointing at moved text"


def test_a_control_site_naming_a_descriptor_bearing_mechanism_is_misclassified(
        monkeypatch):
    """A control is only a control while every mechanism it names is absent."""
    rows = daa.carrying_sites()
    monkeypatch.setattr(daa, "carrying_sites", lambda: tuple(
        (p, ph, r, ("disulfidptosis",) if r is None else c)
        for p, ph, r, c in rows))
    out = _assemble()
    assert out["sites_misclassified"], (
        "a control naming a mechanism with a descriptor is still shipping the "
        "flattening")


def test_a_correction_with_nothing_to_correct_is_misclassified(monkeypatch):
    rows = daa.carrying_sites()
    monkeypatch.setattr(daa, "carrying_sites", lambda: tuple(
        (p, ph, r, ("ttfields",) if r is not None else c)
        for p, ph, r, c in rows))
    out = _assemble()
    assert out["sites_misclassified"], (
        "a site called corrected while every mechanism it names is absent")


def test_the_registry_holds_both_corrected_sites_and_controls():
    out = _assemble()
    assert out["n_sites_corrected"] >= 1
    assert out["n_sites_control"] >= 1, (
        "a registry whose every entry is a defect has not shown where the line "
        "is")


def test_the_registry_names_the_manuscript_and_the_generated_pages():
    paths = {p for p, _, _, _ in daa.carrying_sites()}
    assert "article/drafts/v1.md" in paths
    assert any(p.startswith("analysis/") for p in paths)
    for path, _, _, concerns in daa.carrying_sites():
        assert (REPO / path).exists(), path
        assert concerns, path


# --------------------------------------------------------------------------
# The derived note the consumers import.
# --------------------------------------------------------------------------
def test_the_shared_note_fails_open_when_the_artifact_is_absent(monkeypatch, tmp_path):
    monkeypatch.setattr(daa, "OUT_JSON", tmp_path / "absent.json")
    assert daa.verdicts() == {}
    assert daa.absence_note(["ttfields"]) == []


def test_the_shared_note_splits_a_caller_set_by_verdict(monkeypatch, tmp_path):
    artifact = tmp_path / "audit.json"
    artifact.write_text(json.dumps({"mechanisms": {
        "ttfields": {"verdict": "absent"},
        "cuproptosis": {"verdict": "over-broad"},
        "disulfidptosis": {"verdict": "name-exact-leaky"}}}))
    monkeypatch.setattr(daa, "OUT_JSON", artifact)
    note = " ".join(daa.absence_note(
        ["ttfields", "cuproptosis", "disulfidptosis"]))
    assert "no descriptor names or contains `ttfields`" in note
    assert "CONTAINS `cuproptosis`" in note
    assert "NAME-EXACT descriptor for `disulfidptosis`" in note
    # A caller asking about one mechanism must not be handed the others'.
    only = " ".join(daa.absence_note(["ttfields"]))
    assert "cuproptosis" not in only and "disulfidptosis" not in only


# --------------------------------------------------------------------------
# The registered falsifier and the registered restraints.
# --------------------------------------------------------------------------
def test_the_mechanism_recall_headline_is_unchanged():
    """THE FALSIFIER (#952). A leaky leaf cannot move this, and nothing here may.

    Every descriptor this work found trips its own mechanism's keywords, so
    `mechanism_recall.py`'s classifier calls it LEAKY and its headline pool is
    untouched. If this number moves, the work is wrong.
    """
    recall = json.loads((REPO / "analysis" / "mechanism-recall.json").read_text())
    agg = recall["aggregates"]
    assert agg["n_measurable_mechanisms"] == 9
    assert agg["volume_weighted_recall_leakage_free"] == pytest.approx(0.906, abs=5e-4)
    assert agg["macro_recall_leakage_free"] == pytest.approx(0.892, abs=5e-4)


def test_nothing_here_constrains_the_simulation():
    report = json.loads(
        (REPO / "analysis" / "identifiability-report.json").read_text())
    assert report["data_constrained_in_production"] == 0


def test_the_map_still_carries_a_reason_for_every_unmeasurable_mechanism():
    """Read through `unusable`, which spans BOTH blocks of the map.

    `frequency-therapy` lives in `mechanisms:` with an empty descriptor list
    and the other eight in `unmeasurable:`, so a guard reading one block would
    silently cover eight of nine -- and the one it would miss is the only
    mechanism whose absence verdict a sibling page already publishes.
    """
    mapping = daa.unusable(daa._yaml())
    assert len(mapping) == 9
    assert "frequency-therapy" in mapping, (
        "the set must span both blocks, not the one a mechanism happens to be "
        "in")
    for mech, entry in mapping.items():
        assert entry["verdict"] in daa.VERDICTS, mech
        assert entry["reason"].strip(), mech
        assert mech in MECHANISM_KEYWORDS, mech
        if entry["verdict"] == "name-exact-leaky":
            assert entry.get("leaves"), mech


# --------------------------------------------------------------------------
# A descriptor can be older than the therapy it now names.
# --------------------------------------------------------------------------
def test_a_descriptor_whose_pool_predates_the_modality_is_flagged():
    """`census-unnamed-modalities.md` measures this trap on the keyword arm.

    A name-exact leaf is not exempt from it, so the page must say so rather
    than publishing the pool as a measurement of the modality.
    """
    rows = BASE_ROWS
    detail = _detail(rows, sorted(_map()["unmeasurable"]))
    detail["by_year"]["D004572"] = {"1960": 40, "2020": 45}
    out = daa.assemble({"vocabulary": _vocab(rows), "detail": detail,
                        "map": _map()})
    row = out["mechanisms"]["electrolysis"]["arm_a"][0]
    assert row["vintage_share"] == pytest.approx(40 / 85)
    assert row["vintage_suspect"] is True
    assert "electrolysis" in out["vintage_suspect"]
    md = daa.render(out)
    assert "older than the therapy it is being read as" in md
    assert "upper bound" in md


def test_a_descriptor_entirely_newer_than_the_modality_is_not_flagged():
    rows = BASE_ROWS
    detail = _detail(rows, sorted(_map()["unmeasurable"]))
    detail["by_year"]["D004572"] = {"2020": 85}
    out = daa.assemble({"vocabulary": _vocab(rows), "detail": detail,
                        "map": _map()})
    row = out["mechanisms"]["electrolysis"]["arm_a"][0]
    assert row["vintage_share"] == 0.0
    assert row["vintage_suspect"] is False
    assert "older than the therapy" not in daa.render(out)


def test_the_vintage_year_is_the_one_the_sibling_page_uses():
    """Two copies of a threshold disagree; the import would cycle, so pin them."""
    import census_unnamed_modalities as cum
    assert daa.VINTAGE_YEAR == cum.VINTAGE_YEAR
    assert daa.VINTAGE_SUSPECT_SHARE == pytest.approx(0.10), (
        "the sibling marks a row suspect at a tenth of its pool")


def test_an_undated_descriptor_withholds_its_vintage_rather_than_reading_zero():
    rows = BASE_ROWS
    detail = _detail(rows, sorted(_map()["unmeasurable"]))
    detail["by_year"]["D004572"] = {}
    out = daa.assemble({"vocabulary": _vocab(rows), "detail": detail,
                        "map": _map()})
    row = out["mechanisms"]["electrolysis"]["arm_a"][0]
    assert row["vintage_share"] is None
    assert row["vintage_suspect"] is None


# --------------------------------------------------------------------------
# The phrase rule, which is inert on the current text and drivable only here.
# --------------------------------------------------------------------------
def test_a_registered_phrase_is_found_across_a_line_wrap():
    assert daa._contains("the rest carry no\ndescriptor this map can use",
                         "no descriptor this map can use")


def test_an_absent_phrase_is_still_reported_absent():
    assert not daa._contains("the rest carry a descriptor this map can use",
                             "no descriptor this map can use")


def test_a_phrase_whose_words_join_only_if_whitespace_is_deleted_stays_absent():
    """Collapsing runs is not deleting them, and the two accept different strings."""
    assert not daa._contains("over broad", "overbroad")
    assert daa._contains("over   broad", "over broad")


def test_a_mechanism_with_no_usable_descriptor_and_no_verdict_is_refused():
    """A new unmeasurable mechanism must carry a verdict or the page raises.

    Failing open here would quietly drop a mechanism out of the audit while the
    pages that flatten its verdict carried on saying MeSH has no descriptor for
    it, which is the whole defect.
    """
    mapping = {"mechanisms": {"newcomer": {"descriptors": [], "note": "x"}},
               "unmeasurable": {}}
    with pytest.raises(SystemExit) as excinfo:
        daa.unusable(mapping)
    assert "newcomer" in str(excinfo.value)


def test_a_mechanism_with_a_usable_descriptor_is_not_in_the_unusable_set():
    mapping = {"mechanisms": {"measured": {"descriptors": ["D1"], "note": "x"}},
               "unmeasurable": {}}
    assert daa.unusable(mapping) == {}


def test_the_site_registry_is_data_rather_than_generator_source():
    """A quoted claim carries its author's figures; a generator may not.

    So the registry of quoted sentences lives in a committed CSV, which is the
    arrangement that makes the no-handwritten-figures rule and the quote-it-
    verbatim rule compatible. Checked structurally, because a copy pasted back
    into the source would still render the same page.
    """
    assert daa.SITES_CSV.exists()
    assert "CARRYING_SITES = (" not in SRC, (
        "the registry must not move back into generator source")
    tree = ast.parse(SRC)
    fn = next(n for n in tree.body
              if isinstance(n, ast.FunctionDef) and n.name == "carrying_sites")
    assert "SITES_CSV" in ast.dump(fn), "the loader must read the committed CSV"


def test_an_empty_or_missing_site_registry_raises_rather_than_passing(
        monkeypatch, tmp_path):
    """An empty registry would publish "no sites carry it", which is a finding."""
    monkeypatch.setattr(daa, "SITES_CSV", tmp_path / "gone.csv")
    with pytest.raises(SystemExit):
        daa.carrying_sites()
    empty = tmp_path / "empty.csv"
    empty.write_text("path,phrase,retired,concerns\n")
    monkeypatch.setattr(daa, "SITES_CSV", empty)
    with pytest.raises(SystemExit):
        daa.carrying_sites()


def test_a_mechanism_with_no_keywords_cannot_earn_an_absence_verdict(monkeypatch):
    """Both arms read the mechanism's own vocabulary, so an empty one accepts
    nothing and would hand `absent` to any mechanism that simply has no terms."""
    patched = dict(MECHANISM_KEYWORDS)
    patched["ttfields"] = []
    monkeypatch.setattr(daa, "MECHANISM_KEYWORDS", patched)
    out = _assemble()
    mech = out["mechanisms"]["ttfields"]
    assert mech["n_arm_a"] == 0, "nothing can match an empty vocabulary"
    assert not mech["checks"]["mechanism_has_keywords"]
    assert not out["all_checks_pass"]
    assert "no keyword vocabulary" in daa.render(out)


def test_every_real_unmeasurable_mechanism_has_a_keyword_vocabulary():
    out = _assemble()
    for mech, v in out["mechanisms"].items():
        assert v["checks"]["mechanism_has_keywords"], mech


def test_the_committed_map_passes_every_check_on_the_real_vocabulary():
    """The most important assertion here, and a sweep found it MISSING.

    Changing a declared verdict to one the census refutes left every guard
    green, because each check was driven only through doctored fixtures and
    nothing asserted that the COMMITTED map's own verdicts survive the real
    descriptor vocabulary. That is the state the page exists to publish.
    """
    if not (daa.VOCAB_TSV_GZ.exists() and daa.DETAIL_JSON.exists()):
        pytest.skip("census derivatives not committed")
    out = daa.assemble(daa.load())
    failed = {mech: sorted(k for k, ok in v["checks"].items() if not ok)
              for mech, v in out["mechanisms"].items() if not v["checks_pass"]}
    assert not failed, failed
    assert out["all_checks_pass"] is True
    assert out["detail_missing"] == [], (
        "a candidate with no committed detail means the vocabulary moved and "
        "the second census pass was not re-run")


def test_the_committed_verdicts_are_the_measured_ones():
    """A verdict the census contradicts must not be publishable as `absent`."""
    if not (daa.VOCAB_TSV_GZ.exists() and daa.DETAIL_JSON.exists()):
        pytest.skip("census derivatives not committed")
    out = daa.assemble(daa.load())
    for mech, v in out["mechanisms"].items():
        if v["verdict"] == "name-exact-leaky":
            assert v["n_arm_a"] >= 1, mech
        else:
            assert v["n_arm_a"] == 0, (
                f"{mech} is declared {v['verdict']} while the census carries a "
                "descriptor naming it")
    assert set(out["refuted_absence"]) == {
        m for m, v in out["mechanisms"].items()
        if v["verdict"] == "name-exact-leaky"}


def test_an_inert_threshold_is_declared_and_its_value_bounded():
    """An inert floor can only be caught by a bound on its committed VALUE.

    Deleting a filter that nothing on the real data reaches changes no output
    at all, so a behavioural guard cannot see it -- the shape #930 and #932
    both paid for. Both thresholds here are inert on the current vocabulary, so
    the page declares them and this bounds what they may be set to.
    """
    out = _assemble()
    names = {name for name, _, _ in out["inert_thresholds"]}
    assert names, "an inert threshold must be declared, not hidden"
    assert "POOL_FLOOR" in names
    assert 10 <= daa.POOL_FLOOR <= 100, (
        "the pool floor must stay a real floor; at zero every leaf is "
        "volume-measurable whatever its pool")
    assert 2 <= daa.OVER_BROAD_MULTIPLE <= 20, (
        "at or below one, `over-broad` is satisfied by a cover the size of the "
        "mechanism itself and the verdict stops meaning anything")
    assert "INERT on the current data" in daa.render(out)


def test_a_threshold_that_binds_is_not_declared_inert():
    """The declaration must be measured, not a standing caveat."""
    rows = [dict(r, records=daa.POOL_FLOOR - 1) if r["name"] == "Disulfidptosis"
            else r for r in BASE_ROWS]
    out = _assemble(rows=rows)
    names = {name for name, _, _ in out["inert_thresholds"]}
    assert "POOL_FLOOR" not in names, (
        "a leaf below the floor makes the floor load-bearing, so it must stop "
        "being declared inert")


def test_the_repository_guard_reads_the_repository_and_not_the_fixture_stub(
        monkeypatch):
    """Non-vacuity for the guard above: the real reader must be able to fail.

    Pointing the real reader at a tree where nothing matches has to produce a
    non-empty `sites_missing`; if it does not, that guard is checking a stub
    whatever it says in its name.
    """
    monkeypatch.setattr(daa, "REPO", Path("/nonexistent-tree"), raising=False)
    monkeypatch.setattr(daa, "_site_present", REAL_SITE_PRESENT)
    out = _assemble()
    assert out["sites_missing"], (
        "with no files to read, every registered phrase must be reported "
        "missing -- otherwise the reader is not reading")
