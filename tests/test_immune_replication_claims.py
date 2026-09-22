"""Keep published replication claims tied to the complete archived study."""

import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize("path", ["README.md", "article/drafts/v1.md", "article/drafts/v1.tex", "docs/RESEARCH_NEXT_STEPS.md",
                                  "docs/IMMUNE_2D_REPLICATION_RESULTS.md"])
def test_primary_claim_surfaces_match_the_declared_endpoint_and_interval(path):
    data = json.loads((ROOT / "analysis/immune-2d-replication.json").read_text())
    primary = data["primary"]
    text = " ".join((ROOT / path).read_text().replace("**", "").split())
    assert f"{primary['mean_paired_difference']:.6f}" in text
    low, high = primary["percentile_interval"]
    assert f"[{low:.6f}, {high:.6f}]" in text
    assert "accumulated activation steps per initial tumor cell" in text
    assert "whole-block" in text and "validation" in text
    assert all(not ids for arm in data["threshold_positive_blocks"].values() for ids in arm.values()), (
        "Threshold-positive observations now exist: update every no-crossing claim")


def test_results_table_is_transcribed_from_all_twenty_blocks():
    data = json.loads((ROOT / "analysis/immune-2d-replication.json").read_text())
    rows = [line.split("|")[1:-1] for line in
            (ROOT / "docs/IMMUNE_2D_REPLICATION_RESULTS.md").read_text().splitlines() if line.startswith("|")][2:]
    keys = ["activation_sum_per_initial_tumor_cell", "immune_kills", "eligible_cell_steps", "mean_activation",
            "max_eligible_damp", "completed_releases", "censored_deaths", "release_lp_completed_cohort"]
    assert len(rows) == len(keys)
    for row, key in zip(rows, keys):
        for cell, arm in zip(row[1:], ("RSL3", "SDT")):
            values = [float(n) for n in cell.replace(",", "").replace("[", " ").replace("]", "").split()]
            expected = data["arm_distributions"][arm][key]
            assert expected["n_defined"] == 20
            assert values == pytest.approx([expected[k] for k in ("median", "minimum", "maximum")], rel=5e-6)


def test_secondary_claims_preserve_small_denominator_and_empty_population_limits():
    data = json.loads((ROOT / "analysis/immune-2d-replication.json").read_text())
    ratio = data["sdt_rsl3_immune_kill_ratio"]
    assert ratio["n_defined"] == 20 and ratio["undefined_blocks"] == []
    assert len(data["arm_distributions"]["Control"]["mean_activation"]["undefined_blocks"]) == 12
    for path in ("article/drafts/v1.md", "docs/IMMUNE_2D_REPLICATION_RESULTS.md"):
        text = " ".join((ROOT / path).read_text().replace("**", "").split())
        assert f"median {ratio['median']:.1f}" in text
        assert f"{ratio['minimum']:.0f}–{ratio['maximum']:.0f}" in text
        assert "Twelve Control blocks" in text and "undefined" in text
        assert "denominator" in text and "causal" in text


def test_pages_claim_of_a_contrast_in_every_block_has_all_twenty_observations():
    data = json.loads((ROOT / "analysis/immune-2d-replication.json").read_text())
    assert [b["block"] for b in data["blocks"]] == list(range(1, 21))
    assert all(b["arms"]["SDT"]["activation_sum_per_initial_tumor_cell"]
               > b["arms"]["RSL3"]["activation_sum_per_initial_tumor_cell"] for b in data["blocks"])
    text = " ".join((ROOT / "docs/index.html").read_text().split())
    assert "All 20 new seed blocks retain the activation contrast" in text
    assert "independent validation remains pending" in text
