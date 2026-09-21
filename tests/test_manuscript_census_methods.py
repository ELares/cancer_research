"""Keep the active census methods tied to their committed inputs.

Scope is Sections 3.4 and 3.5, the adjacent matrix comparison in 3.6,
and the interpretation of mechanism-pair ordering in 3.13:
the historical keyword instrument is still described elsewhere, and a correct
number elsewhere cannot rescue a stale methods claim. These checks bind named
quantities and consequential method distinctions, without snapshotting the
prose or needing the raw census.
"""

import importlib.util
import itertools
import json
import re
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[1]


def _section(number):
    text = (ROOT / "article/drafts/v1.md").read_text()
    match = re.search(rf"^### {re.escape(number)} [^\n]+\n(.*?)(?=^##[# ]|\Z)",
                      text, re.M | re.S)
    assert match, f"Missing manuscript section {number}"
    return match.group(1)


def _paragraph(number, label):
    match = re.search(rf"^\*\*{re.escape(label)}\.\*\* (.*?)(?=\n\n|\Z)",
                      _section(number), re.M | re.S)
    assert match, f"Section {number} no longer identifies {label}"
    return " ".join(match.group(1).split())


def _data(name):
    return json.loads((ROOT / "analysis" / f"{name}.json").read_text())


def _module(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _claim(pattern, text):
    match = re.search(pattern, text, re.I)
    assert match, f"Missing methods claim matching {pattern!r}"
    return match


def _integer(value):
    return int(value.replace(",", ""))


def test_descriptor_group_counts_and_matching_follow_the_maps():
    groups = yaml.safe_load((ROOT / "analysis/mesh-mechanism-map.yaml").read_text())["mechanisms"]
    mapped = {name for name, group in groups.items() if group["descriptors"]}
    sites = {line.split("\t")[0] for line in
             (ROOT / "analysis/site-descriptor-map.tsv").read_text().splitlines()
             if line.strip() and not line.startswith("#")}
    mechanism = _paragraph("3.4", "Mechanism groups")
    site = _paragraph("3.4", "Cancer-site groups")
    assert int(_claim(r"(\d+) nonempty mechanism descriptor groups", mechanism)[1]) == len(mapped)
    assert int(_claim(r"(\d+) site groups", site)[1]) == len(sites)
    assert set(_data("census-mechanism-cancer-matrix")["mechanism_totals"]) == mapped
    assert set(_data("census-mechanism-cancer-matrix")["site_totals"]) == sites
    for descriptor in groups["sonodynamic"]["descriptors"]:
        assert f"`{descriptor}`" in mechanism
    _claim(r"exact and case-insensitive", mechanism)
    _claim(r"once to each group", mechanism)
    _claim(r"unavailable, not zero", mechanism)
    _claim(r"once to each matching site", site)
    _claim(r"assignments can exceed.*unique articles", site)


def test_design_inventory_and_claimed_precedence_match_the_classifier():
    text = _paragraph("3.4", "Study-design classes")
    order = _claim(r"precedence order: ([^.]+)\.", text)[1].split(", ")
    assert len(order) == len(set(order))
    assert set(order) == set(_data("census-evidence-design")["classes"])
    classifier = _module("census_evidence_design")
    examples = {
        "trial": ([sorted(classifier.TRIAL)[0]], []),
        "non-primary": ([sorted(classifier.NON_PRIMARY)[0]], []),
        "clinical-other": ([sorted(classifier.CLINICAL_OTHER)[0]], []),
        "animal-model": ([], [sorted(classifier.IN_VIVO)[0]]),
        "cell-culture": ([], [sorted(classifier.IN_VITRO)[0]]),
        "animal-other": ([], [sorted(classifier.ANIMAL_ANY)[0]]),
        "undetermined": ([], []),
    }
    # The display ORDER constant is not precedence. Exercise every ordering
    # the manuscript declares with records carrying both kinds of metadata.
    for earlier, later in itertools.combinations(order, 2):
        p1, m1 = examples[earlier]
        p2, m2 = examples[later]
        assert classifier.classify(p1 + p2, m1 + m2) == earlier, (earlier, later)
    _claim(r"animal-other requiring.*no `Humans` tag", text)
    _claim(r"design classes, not evidence-quality tiers", text)


def test_design_coverage_names_both_denominators_and_keeps_undetermined():
    data = _data("census-evidence-design")
    text = _paragraph("3.4", "Study-design classes")
    census = _integer(_claim(r"covers ([\d,]+) records", text)[1])
    classified = _integer(_claim(r"([\d,]+) classifiable", text)[1])
    undetermined = _claim(r"([\d,]+) undetermined \(([\d.]+)%\)", text)
    assert census == data["census"]
    assert classified == data["classifiable"] == census - data["classes"]["undetermined"]
    assert _integer(undetermined[1]) == data["classes"]["undetermined"]
    assert float(undetermined[2]) == round(100 * data["classes"]["undetermined"] / census, 1)
    _claim(r"full-census and classifiable denominators", text)
    _claim(r"undetermined class is retained", text)


def test_matrix_dimensions_and_marginals_use_the_both_label_population():
    data = _data("census-mechanism-cancer-matrix")
    text = _paragraph("3.5", "Mechanism-site matrix")
    dimensions = _claim(r"(\d+) mechanism groups by (\d+) site groups.*?(\d+) cells", text)
    m, s, cells = map(int, dimensions.groups())
    assert m == len(data["mechanism_totals"])
    assert s == len(data["site_totals"])
    assert cells == m * s == data["n_cells"]
    assert _integer(_claim(r"universe is ([\d,]+) articles", text)[1]) == data["universe"]
    _claim(r"both at least one mapped mechanism and at least one mapped site", text)
    _claim(r"mechanism marginal multiplied by the site marginal and divided by this universe", text)
    _claim(r"both marginals use the same eligible articles", text)


def test_matrix_floor_applies_to_expectation_and_ratios_are_not_opportunities():
    text = _paragraph("3.5", "Mechanism-site matrix")
    threshold = _claim(r"expected count of at least ([\d.]+)", text)
    assert float(threshold[1]) == _data("census-mechanism-cancer-matrix")["min_expected"]
    _claim(r"not an observed-count cutoff", text)
    _claim(r"descriptive ratios, not significance tests", text)
    _claim(r"(?:low ratios|empty cells).*do not establish.*research opportunity", text)


def test_cotags_count_articles_without_claiming_tested_combinations():
    text = _paragraph("3.5", "Mechanism co-tagging")
    _claim(r"unordered pair once per article", text)
    _claim(r"deduplicates.*per-mechanism lists", text)
    _claim(r"partner ordering ranks raw co-occurrence counts", text)
    _claim(r"without adjustment for partner prevalence or descriptor breadth", text)
    _claim(r"include reviews and incidental coverage", text)
    _claim(r"do not establish that a combination was tested", text)
    _claim(r"not a field-wide rate", text)


def test_profile_summaries_remain_conditional_on_descriptor_selection():
    text = _paragraph("3.5", "Profile interpretation")
    _claim(r"trial share, site enrichment and partner ordering.*articles selected", text)
    _claim(r"all three can change with descriptor coverage or indexing", text)
    _claim(r"normalization does not establish invariance or comparability", text)


def test_pair_ordering_requires_a_paired_comparison_not_an_aggregate_rate():
    text = " ".join(_section("3.13").split())
    _claim(r"pair counts and ordering.*describe co-tagging under the current MeSH descriptor map", text)
    _claim(r"aggregate co-occurrence-rate comparison does not establish ranking stability", text)
    _claim(r"paired comparison of partner rankings on the same articles using both labelling methods", text)
    _claim(r"no such comparison is reported", text)
    _claim(r"historical retrieved corpus", text)
    assert "stable in ordering under both instruments" not in text


def test_profile_site_enrichment_keeps_its_assignment_denominator_and_observed_floor():
    text = _paragraph("3.5", "Profile site enrichment")
    _claim(r"site assignments within a mechanism divided by.*site assignments across the indexed census", text)
    _claim(r"including articles without a mapped mechanism", text)
    _claim(r"denominators count assignments, not unique articles", text)
    _claim(r"must not be substituted for the matrix's expected-count floor", text)
    threshold = int(_claim(r"at least (\d+) observed articles", text)[1])
    data = _data("census-mechanism-profile")
    probe = {
        **data,
        "count": {"probe": threshold * 2},
        "site_totals": {"below": threshold * 3, "at": threshold * 4},
        "by_site": {"probe": {"below": threshold - 1, "at": threshold}},
    }
    row = _module("census_mechanism_profile").assemble(probe)["rows"][0]
    assert {site["site"] for site in row["top_sites"]} == {"at"}


def test_trial_share_uses_all_mechanism_articles_and_trial_specific_types():
    text = _paragraph("3.5", "Trial-publication representation")
    _claim(r"share of all articles in each mechanism group", text)
    _claim(r"denominator includes non-primary and undetermined", text)
    _claim(r"not restricted to classifiable designs", text)
    excluded = re.findall(r"`([^`]+)`", text)
    assert excluded, "The methods must name the non-trial publication types"
    trial_types = _module("census_mechanism_profile").TRIAL_TYPES
    assert not set(excluded) & trial_types
    assert set(excluded) == {"Multicenter Study", "Comparative Study", "Clinical Study"}
    _claim(r"alone do not qualify as trials", text)
    _claim(r"not unique trials.*highest evidence level", text)


def test_growth_windows_and_starting_count_floors_follow_the_analysis():
    data = _data("census-mechanism-growth")
    text = _paragraph("3.5", "Matched publication growth")
    windows = [tuple(map(int, pair)) for pair in re.findall(r"(\d{4})[–-](\d{4})", text)]
    assert windows == [(data["start_year"], data["end_year"]),
                       (data["recent_start_year"], data["end_year"])]
    decade = int(_claim(r"decade ratios require at least (\d+) articles", text)[1])
    recent = int(_claim(r"recent percentage changes require at least (\d+)", text)[1])
    assert decade == data["min_base"]
    # The recent floor is not stored in JSON. Check the manuscript's boundary
    # against the real assembler, without loading or scanning raw records.
    analysis = _module("census_mechanism_growth")
    probe = {**data, "mechanism_by_year": {
        "below": {str(data["recent_start_year"]): recent - 1},
        "at": {str(data["recent_start_year"]): recent},
    }}
    rows = {row["mechanism"]: row for row in analysis.assemble(probe)["rows"]}
    assert rows["below"]["recent_pct"] is None
    assert rows["at"]["recent_pct"] is not None
    _claim(r"field series includes every dated article", text)
    _claim(r"union series counts an article once", text)
    _claim(r"union's endpoint growth factor by the field's factor using unrounded endpoint counts", text)
    _claim(r"not clinical progress or efficacy", text)


def test_keyword_ladder_and_old_matrix_are_explicitly_historical():
    history = _paragraph("3.4", "Historical comparison")
    _claim(r"frozen retrieved corpus", history)
    _claim(r"instrument-comparison arm", history)
    _claim(r"do not define the current census classifications or their denominators", history)
    active = _section("3.4").split("**Historical comparison.**", 1)[0] + _section("3.5")
    # Permit the matrix's explicit historical comparison, but never let the
    # old matrix or evidence ladder become the active method again.
    active = re.sub(r"The historical [^.]+\.", "", active)
    assert not re.search(r"\b19\s*(?:x|×|-by-)\s*22\b", active)
    assert not re.search(r"\b(?:seven-tier|seven tiers|19 (?:keyword )?(?:mechanism )?tags|22 cancer-type tags)\b", active)
    assert "Classification was based on keyword detection" not in active
    limitations = " ".join(_section("3.6").split())
    _claim(r"both the population and the taxonomy differ from the historical matrix", limitations)
    _claim(r"comparison cannot isolate the effect of retrieval", limitations)
