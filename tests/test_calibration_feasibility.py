"""Guards for the #616 calibration-feasibility scan.

`analysis/atlas-model-gaps.md` named four high-attention unmodelled ferroptosis
mechanisms. The layer-freeze policy requires a named calibration target before
any becomes a layer, and `scripts/calibration_feasibility.py` found the answer
differs by threshold: the z<-1 cut recovers the normal for every gene, while the
z<-2 cut separates TP53 and ACSL4 from the rest.

BOTH halves need guarding, and an earlier version of this file guarded only the
first. That let the central negative be inverted without a single failure, and it
let a prose fix pass by pinning one literal string. Every assertion here is
therefore tied to a recomputed quantity or to a behaviour, and the deep cut is
covered for every gene rather than for TP53 alone.
"""

import csv
import json
import statistics as st
import sys
from math import erf, sqrt
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "scripts"))

from calibration_feasibility import quartile_spread, sign_test  # noqa: E402

CAL = REPO_ROOT / "analysis" / "calibration"
RAW = CAL / "calibration-feasibility.json"
DOC = CAL / "calibration-feasibility.md"
CSV_SRC = CAL / "acsl4_prevalence_tcga.csv"
ACSL4_DOC = CAL / "acsl4-prevalence-calibration.md"
ACSL4_RS = REPO_ROOT / "simulations" / "ferroptosis-core" / "src" / "acsl4.rs"
STATUS = REPO_ROOT / "simulations" / "calibration" / "CALIBRATION_STATUS.md"

GAPS = ["HMOX1", "TP53", "TFRC", "KEAP1"]
ALL_GENES = ["ACSL4", "GPX4", "SLC7A11"] + GAPS
EXPECTED_LOW = 0.5 * (1 + erf(-1 / sqrt(2)))       # 15.87%
EXPECTED_VERYLOW = 0.5 * (1 + erf(-2 / sqrt(2)))   # 2.28%


def _raw():
    return json.loads(RAW.read_text())


def _csv():
    return list(csv.DictReader(CSV_SRC.open()))


def test_artifacts_exist_and_cover_every_gene():
    d = _raw()
    for g in ALL_GENES:
        assert g in d["genes"], f"{g} missing from the feasibility scan"
        assert d["genes"][g]["n_types"] >= 30, g
    assert DOC.exists()


def test_every_committed_statistic_is_recomputable_from_the_csv():
    """The whole finding must be re-derivable, for EVERY gene at BOTH cuts.

    Covers all seven genes rather than the four gaps, because the #462
    correction turns on ACSL4's numbers and an earlier version left them
    unguarded.
    """
    rows, d = _csv(), _raw()
    for g in ALL_GENES:
        lo = [float(r[f"{g}_frac_low"]) for r in rows if r[f"{g}_frac_low"]]
        vl = [float(r[f"{g}_frac_verylow"]) for r in rows if r[f"{g}_frac_verylow"]]
        s = d["genes"][g]
        assert abs(st.median(lo) - s["median_low"]) < 1e-12, g
        assert abs(st.median(vl) - s["median_verylow"]) < 1e-12, g
        assert abs(quartile_spread(lo) - s["iqr_low"]) < 1e-12, g
        assert sum(1 for x in lo if x > EXPECTED_LOW) == s["types_above_expected_low"], g
        assert sum(1 for x in vl if x > EXPECTED_VERYLOW) == s["types_above_expected_verylow"], g


def test_the_shallow_cut_recovers_the_normal_for_every_gene():
    """Half one of the finding: z<-1 carries no gene-specific information.

    If a gene ever lands outside the band this must fail and be re-read, not
    widened -- a real deviation would mean the cut has become informative.
    """
    d = _raw()
    assert d["separated_low"] == [], (
        f"the z<-1 cut now separates {d['separated_low']}; the finding has changed")
    for g, s in d["genes"].items():
        assert abs(s["median_low"] - EXPECTED_LOW) < 0.025, (
            f"{g} median_low {s['median_low']:.4f} vs normal {EXPECTED_LOW:.4f}")


def test_the_deep_cut_separates_exactly_tp53_and_acsl4():
    """Half two, which the earlier version of this file did not cover at all.

    Pinned as a SET, so inflating any gene's deep tail fails here -- including
    the controls, whose flatness is what makes the ACSL4 result meaningful.
    """
    d = _raw()
    assert set(d["separated_verylow"]) == {"TP53", "ACSL4"}, d["separated_verylow"]
    for g in ["GPX4", "SLC7A11", "HMOX1", "TFRC", "KEAP1"]:
        s = d["genes"][g]
        assert not s["separated_verylow"], f"{g} now separates at z<-2"
        assert s["median_verylow"] < 0.025, (
            f"{g} deep tail {s['median_verylow']:.4f} is no longer near expectation")
    tp53, acsl4 = d["genes"]["TP53"], d["genes"]["ACSL4"]
    assert tp53["median_verylow"] > acsl4["median_verylow"] > EXPECTED_VERYLOW
    assert tp53["types_above_expected_verylow"] >= 30
    assert acsl4["types_above_expected_verylow"] >= 20


def test_no_gap_gene_is_anchored_at_the_shallow_cut():
    """The four are unanchored where the route is definitionally uninformative."""
    d = _raw()
    assert set(d["gaps"]) == set(GAPS)
    for g in GAPS:
        assert not d["genes"][g]["separated_low"], g


def test_only_tp53_among_the_gaps_has_any_signal():
    """HMOX1/TFRC/KEAP1 are unanchored at BOTH cuts; TP53 is the exception."""
    d = _raw()
    for g in ["HMOX1", "TFRC", "KEAP1"]:
        assert not d["genes"][g]["separated_verylow"], (
            f"{g} now has a deep-cut signal and its data-blocked row is stale")
    assert d["genes"]["TP53"]["separated_verylow"]


def test_sign_test_and_quartile_spread_are_correct():
    """The two statistics the finding rests on, checked against known values."""
    assert sign_test(16, 32) == 1.0
    assert sign_test(31, 32) < 1e-6
    assert 0.001 < sign_test(25, 32) < 0.01
    assert sign_test(0, 0) == 1.0
    # A true IQR, not the order-statistic difference this once used.
    assert abs(quartile_spread([1, 2, 3, 4, 5]) - 2.0) < 1e-12
    assert quartile_spread([7]) == 0.0


def test_status_doc_records_every_gap_gene_with_its_verdict():
    """#616's acceptance criterion: a row per gene, with the reason."""
    txt = STATUS.read_text()
    assert "Layers proposed and NOT built" in txt
    section = txt[txt.index("Layers proposed and NOT built"):]
    section = section[:section.index("What would unblock")]
    for g in GAPS:
        assert g in section, f"{g} has no row in the not-built table"
    assert section.count("data-blocked") >= 3
    assert "weak anchor available" in section, "TP53's verdict is missing"


def test_tfrc_is_not_claimed_to_be_covered_by_ferritinophagy():
    """Refuted reason: ferritin release is not transferrin-receptor import.

    #616 states the engine models ferritinophagy but NOT TFRC import, and
    params.rs documents the field as NCOA4-driven release from intracellular
    stores. An earlier draft used 'already absorbed' as a data-blocked reason.
    """
    for path in (STATUS, DOC):
        txt = path.read_text()
        assert "absorbed into `ferritinophagy_release`" not in txt, path
        # Each doc must state the distinction itself; an `or` across the two
        # would let one of them silently revert.
        assert "transferrin-receptor import" in txt, f"{path.name} drops the distinction"
        assert "uncovered axis" in txt, f"{path.name} drops the TFRC verdict"


def test_the_correction_is_scoped_to_the_shallow_cut_everywhere():
    """The #462 correction must not overreach onto the z<-2 row, which survives.

    Checked in all three places the framing lives, including the shipped Rust
    doc comment -- a reader landing on the library API would otherwise get the
    retracted story with no signal.
    """
    doc = ACSL4_DOC.read_text()
    assert "Correction (#616)" in doc
    assert "z < -2` row survives" in doc or "z < -2` row SURVIVES" in doc, \
        "the correction does not say the deep cut survives"
    assert "the usable population prior" not in doc, \
        "the section heading still calls the shallow cut the usable prior"

    rs = ACSL4_RS.read_text()
    assert "#616" in rs, "acsl4.rs was not updated"
    assert "median 14%) is the committed population prior" not in rs, \
        "acsl4.rs still carries the retracted z<-1 framing"
    # The prior it points at must be the DEEP cut. Pinning the corrected
    # sentence, not merely the absence of the old one -- flipping the threshold
    # back to z<-1 otherwise passes, since the doc legitimately names both.
    assert "fraction of tumors with `z < -2`" in rs, \
        "acsl4.rs no longer points the prior at the deep cut"
    assert "Use the DEEP cut, not `z < -1`" in rs


def test_generated_docs_do_not_contradict_each_other_on_tp53():
    """census-findings once said 'cannot anchor anything' while the feasibility
    doc granted TP53 a weak anchor. Two generated docs, opposite conclusions."""
    findings = (REPO_ROOT / "analysis" / "census-findings.md").read_text()
    assert "cannot anchor anything" not in findings
    if "#616" in findings:
        assert "weak anchor" in findings or "TP53 does separate" in findings
