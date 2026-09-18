"""Keep the combination claim tied to the actual uncertainty artifact.

The manuscript and P1 protocol called a lower prior-predictive bound rounded
to 1.000 robust supra-additivity, despite a sampled minimum of 0.953. They also
called an observed/expected ratio an excess and treated failure of a synergy
threshold as proof that the two biological pathways were not independent.

These checks read the committed report rather than repeating its numbers in
the tests. They guard the interpretation and decision-rule distinction in the
reader-facing sections, without modifying the historical preregistration.
"""

import json
import re
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
MANUSCRIPT = ROOT / "article/drafts/v1.md"
PROTOCOL = ROOT / "analysis/p1-wetlab-protocol.md"


def _section(text, heading):
    match = re.search(rf"^{re.escape(heading)}\n(.*?)(?=^###? |\Z)", text,
                      re.M | re.S)
    assert match, f"Missing section: {heading}"
    return match.group(1)


def _report_interval(label):
    report = (ROOT / "analysis/headline-uncertainty-report.md").read_text()
    match = re.search(rf"^\| {re.escape(label)} \| \[([\d.]+), ([\d.]+)\]",
                      report, re.M)
    assert match, f"Missing uncertainty quantity: {label}"
    return tuple(float(value) for value in match.groups())


@pytest.mark.parametrize("path", [MANUSCRIPT, PROTOCOL])
def test_combination_section_discloses_both_interval_and_sampled_range(path):
    text = path.read_text()
    if path == MANUSCRIPT:
        text = _section(text, "### 7.3 Mechanistic Combination Modeling: Why Drug Pairs Synergize")
    else:
        text = _section(text, "## Model prediction being tested")
    for label in ("95% prior-predictive interval", "full range (min, max)"):
        lo, hi = _report_interval(label)
        assert f"[{lo:.3f}, {hi:.3f}]" in text, (
            f"{path.name} no longer reports {label} from the committed artifact")
    assert "additive null" in text and "sub-additive" in text, (
        "The reader must see why the interval is not guaranteed synergy")
    assert "ratio" in text and "excess" in text, (
        "The ratio and excess must be distinguished at the point of use")


def test_uncertainty_summary_does_not_reinstate_robust_supra_additivity():
    lo, _ = _report_interval("95% prior-predictive interval")
    minimum, _ = _report_interval("full range (min, max)")
    if lo > 1 and minimum > 1:
        pytest.fail("The artifact changed to uniformly supra-additive draws; review the claim")
    text = MANUSCRIPT.read_text()
    abstract = _section(text, "## Abstract")
    design = _section(text, "### 5.2 Simulation Design")
    for body in (abstract, design):
        assert "supra-additive direction is robust" not in body
        assert f"{minimum:.3f}" in body, (
            "The sub-additive sampled tail is missing from an uncertainty summary")


def test_protocol_keeps_the_registered_threshold_and_limits_its_interpretation():
    registration = (ROOT / "PREREGISTRATION.md").read_text()
    p1 = registration.split("**P1.", 1)[1].split("**P2.", 1)[0]
    threshold = re.search(r"combination index greater than ([\d.]+)", p1).group(1)
    protocol = " ".join(PROTOCOL.read_text().replace("**", "").split())
    assert f"greater than {threshold}" in protocol
    assert "at or below the Bliss-independence prediction within assay error" in protocol
    assert "does not amend `PREREGISTRATION.md`" in protocol
    interpretation = _section(PROTOCOL.read_text(), "## Expected result and what a refutation means")
    assert "does not by itself show that the two defenses are not independent" in " ".join(interpretation.split())


def test_brequinar_is_a_separate_comparison_with_primary_source_attribution():
    for text in (MANUSCRIPT.read_text(), PROTOCOL.read_text()):
        assert "10.1038/s41586-021-03539-7" in text
        assert "10.1038/s41586-023-06269-0" in text
        assert "10.1038/s41586-023-06270-7" in text
        assert "iFSP1 or DHODH inhibitor brequinar" not in text
        assert "target engagement" in text or "target-engagement" in text


def test_roadmap_combination_reference_resolves_to_the_combination_section():
    text = MANUSCRIPT.read_text()
    sections = dict(re.findall(r"^### (\d+\.\d+) (.+)$", text, re.M))
    match = re.search(r"Dual-pathway combination response \(Section (\d+\.\d+)\)", text)
    assert match, "The roadmap no longer distinguishes the combination response"
    assert "Combination Modeling" in sections[match.group(1)]


def test_model_card_size_comparison_uses_the_committed_aggregate_measurement():
    data = json.loads((ROOT / "analysis/modality-module-depth.json").read_text())
    card = " ".join((ROOT / "MODEL_CARD.md").read_text().split())
    assert f"{data['dedicated_code_lines']:,} production lines" in card
    assert f"{data['line_ratio_narrow']} to {data['line_ratio_wide']}" in card
    assert "aggregate code-size comparison" in card


def test_current_arm_inventory_matches_the_treatment_enum():
    source = (ROOT / "simulations/ferroptosis-core/src/cell.rs").read_text()
    body = re.search(r"pub enum Treatment\s*\{(.*?)\n\}", source, re.S).group(1)
    variants = set(re.findall(r"^\s+([A-Z]\w*),$", body, re.M))
    assert "Control" in variants
    claim = f"{len(variants - {'Control'})} selectable treatment arms plus an untreated control"
    assert claim in MANUSCRIPT.read_text()
    assert claim in (ROOT / "MODEL_CARD.md").read_text()


def test_immune_claim_does_not_infer_a_per_dead_cell_ratio_from_population_means():
    source = (ROOT / "simulations/ferroptosis-core/src/params.rs").read_text()
    default = source.split("impl Default for Params", 1)[1].split("impl Params", 1)[0]
    threshold = re.search(r"death_threshold: ([\d.]+)", default).group(1)
    body = _section(MANUSCRIPT.read_text(), "### 8.2 Immune Coupling: DAMP-Mediated T Cell Activation")
    assert f"default threshold of {threshold}" in body
    assert "not LP conditional on death in this matched spatial run" in body
    assert "No treatment-specific per-dead-cell LP or DAMP ratio is established" in body
    for unsupported in ("~7.8", "~2.6× more DAMPs", "LP reaches ~20"):
        assert unsupported not in body
    caption_source = (ROOT / "scripts/generate_latex.py").read_text()
    caption = next(line for line in caption_source.splitlines()
                   if "'16': ('fig19_immune_coupling_flow'" in line)
    assert "total immune-kill ratio" in caption
    assert "not a per-dead-cell DAMP or immunogenicity ratio" in caption
