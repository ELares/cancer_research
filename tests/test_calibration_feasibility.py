"""Guards for the #616 calibration-feasibility scan.

`analysis/atlas-model-gaps.md` named four high-attention unmodelled ferroptosis
mechanisms. The layer-freeze policy requires a named calibration target before any
of them becomes a layer, and `scripts/calibration_feasibility.py` found none:
within-cohort mRNA z-scores recover the normal distribution for every gene tested,
so the route that partially anchored ACSL4 (#462) carries no gene-specific signal.

These pin that negative result, which is the load-bearing one -- a negative that
rots back into a positive is how an unanchored layer gets written.
"""

import json
import statistics as st
from math import erf, sqrt
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
CAL = REPO_ROOT / "analysis" / "calibration"
RAW = CAL / "calibration-feasibility.json"
DOC = CAL / "calibration-feasibility.md"
ACSL4_DOC = CAL / "acsl4-prevalence-calibration.md"
STATUS = REPO_ROOT / "simulations" / "calibration" / "CALIBRATION_STATUS.md"

GAPS = ["HMOX1", "TP53", "TFRC", "KEAP1"]
EXPECTED_LOW = 0.5 * (1 + erf(-1 / sqrt(2)))       # 15.87%
EXPECTED_VERYLOW = 0.5 * (1 + erf(-2 / sqrt(2)))   # 2.28%


def _raw():
    return json.loads(RAW.read_text())


def test_artifacts_exist_and_cover_every_gene():
    d = _raw()
    for g in ["ACSL4", "GPX4", "SLC7A11"] + GAPS:
        assert g in d["genes"], f"{g} missing from the feasibility scan"
        assert d["genes"][g]["n_types"] >= 30, g
    assert DOC.exists()


def test_the_zscore_route_recovers_the_normal_for_every_gene():
    """The finding: no gene's low-expression fraction is gene-specific.

    Every gene sits within 2.5 points of P(z < -1) = 15.87%. If a gene ever lands
    outside that band this test should fail and be re-read, not widened -- a real
    deviation would mean the route has become informative for that gene.
    """
    d = _raw()
    for g, s in d["genes"].items():
        assert abs(s["median_low"] - EXPECTED_LOW) < 0.025, (
            f"{g} median_low {s['median_low']:.3f} deviates from the normal "
            f"expectation {EXPECTED_LOW:.4f}; the route may now carry signal")


def test_no_gap_gene_has_a_calibration_target():
    """None of the four is anchored, so no layer may be written for them."""
    d = _raw()
    assert set(d["gaps"]) == set(GAPS)
    for g in GAPS:
        assert abs(d["genes"][g]["dev_from_normal"]) < 0.025, g


def test_tp53_deep_tail_is_the_one_surviving_signal():
    """TP53's z < -2 tail and inter-cancer spread do deviate; the others' do not.

    This is the honest counterweight to the negative above, and the reason TP53 is
    recorded as 'weak anchor available' rather than flatly data-blocked.
    """
    d = _raw()["genes"]
    tp53 = d["TP53"]
    others = [s for g, s in d.items() if g != "TP53"]
    assert tp53["median_verylow"] > EXPECTED_VERYLOW, tp53["median_verylow"]
    assert tp53["median_verylow"] > max(s["median_verylow"] for s in others)
    assert tp53["iqr_low"] > max(s["iqr_low"] for s in others)


def test_status_doc_records_every_gap_gene_as_not_built():
    """The #616 acceptance criterion: a row per gene, with the reason."""
    txt = STATUS.read_text()
    assert "Layers proposed and NOT built" in txt
    for g in GAPS:
        assert g in txt, f"{g} has no CALIBRATION_STATUS row"
    assert "data-blocked" in txt


def test_acsl4_doc_carries_the_correction_and_not_the_old_framing():
    """#462 called the figure 'the real, committed population prior'.

    Six control genes return the same number, so it describes the z-score, not
    ACSL4. The old phrasing must not come back without the control alongside it.
    """
    txt = ACSL4_DOC.read_text()
    assert "Correction (#616)" in txt
    assert "not an ACSL4-specific prior" in txt
    assert "the real, committed population prior" not in txt.split("Correction (#616)")[0]


def test_control_genes_are_measured_not_asserted():
    """The correction rests on measurements in the committed CSV, not on prose."""
    import csv
    rows = list(csv.DictReader((CAL / "acsl4_prevalence_tcga.csv").open()))
    for g in GAPS:
        vals = [float(r[f"{g}_frac_low"]) for r in rows if r.get(f"{g}_frac_low")]
        assert len(vals) >= 30, f"{g} not measured across the cohort"
        assert abs(st.median(vals) - _raw()["genes"][g]["median_low"]) < 1e-9, g
