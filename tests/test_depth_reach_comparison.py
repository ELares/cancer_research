"""Guards for `analysis/depth-reach-comparison.{md,json}`.

The page makes one claim the manuscript has never been able to make: that
"physical modality" spans a range, from light dying in millimetres to
megavoltage photons that barely notice a tumour. It is only worth anything if
the two columns it puts side by side are genuinely independent -- delivered
ENERGY from each modality's own attenuation law, and observed KILL from the
binary -- because the finding is that they come apart.

So the guards are about that independence, and about the caveats being real
rather than decorative.
"""
import importlib.util
import json
import math
import re
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
MD = REPO / "analysis/depth-reach-comparison.md"
JSON_ = REPO / "analysis/depth-reach-comparison.json"
CORE = REPO / "simulations/ferroptosis-core/src"


def _load():
    spec = importlib.util.spec_from_file_location(
        "depth_reach_comparison", REPO / "scripts/depth_reach_comparison.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


DR = _load()


@pytest.fixture(scope="module")
def d():
    return json.loads(JSON_.read_text())


@pytest.fixture(scope="module")
def md():
    return MD.read_text()


def test_the_constants_are_read_from_the_crate_not_typed(d):
    """Every attenuation constant must still be the Rust one.

    Panel (b) of Figure 8 has its own drift guard for the same reason; this is
    the artifact's. A constant typed here and changed there is how a page
    describing a binary stops describing it.
    """
    live = DR._constants()
    assert d["constants"] == live
    src = (CORE / "radiation.rs").read_text()
    m = re.search(r"MU_6MV_SOFT_TISSUE_PER_CM: f64 = ([0-9.]+)", src)
    assert m and float(m.group(1)) == d["constants"]["radiation_mu_per_cm"]


def test_delivered_energy_is_computed_from_the_law_not_from_the_kill(d):
    """The two columns must be INDEPENDENT, which is the whole point.

    Recomputed here from the law rather than read back, so a generator that
    started deriving one from the other fails.
    """
    c = d["constants"]
    for r in d["rows"]:
        tx, z = r["treatment"], r["deep_mm"]
        if r["delivered_at_deep_pct"] is None:
            continue
        assert abs(DR.delivered(tx, z, c) * 100.0 - r["delivered_at_deep_pct"]) < 1e-9
    # And they really do come apart: RSL3 delivers everything and kills least,
    # PDT delivers least and kills a lot at the surface. If this ever stops
    # being true the page's central paragraph is wrong.
    by = {r["treatment"]: r for r in d["rows"]}
    assert by["RSL3"]["delivered_at_deep_pct"] == 100.0
    assert by["RSL3"]["overall_kill_pct"] < by["PDT"]["overall_kill_pct"], (
        "RSL3 now outkills PDT overall; the dissociation paragraph needs "
        "re-deriving")
    assert by["PDT"]["delivered_at_deep_pct"] < 20.0
    assert by["PDT"]["surface_kill_pct"] > 50.0


def test_radiation_is_the_flat_arm_and_that_is_measured(d, md):
    """The reason the arm was added. Asserted against BOTH columns, because
    'reaches deeper' and 'works deeper' are different claims and the page's own
    caveat says so."""
    rad = next(r for r in d["rows"] if r["treatment"] == "Radiation")
    assert rad["delivered_at_deep_pct"] > 95.0, (
        f"radiation delivers only {rad['delivered_at_deep_pct']}% at "
        f"{rad['deep_mm']} mm; it is no longer the flat arm")
    assert rad["kill_retained_pct"] > 85.0
    for other in ("PDT", "SDT"):
        o = next(r for r in d["rows"] if r["treatment"] == other)
        assert rad["delivered_at_deep_pct"] > o["delivered_at_deep_pct"], other
    assert "range the phrase" in md


def test_the_radiation_kill_matches_the_linear_quadratic_expectation(d):
    """Radiation's column is the one with a published parameterisation behind
    it, so it can be checked against the model rather than only against
    itself. A drift here means the binary and the crate have parted.
    """
    src = (CORE / "radiation.rs").read_text()
    alpha = float(re.search(
        r"ALPHA_GBM_PARAMETERISATION_PER_GY: f64 = ([0-9.]+)", src).group(1))
    ratio = float(re.search(
        r"ALPHA_BETA_TUMOUR_GY: f64 = ([0-9.]+)", src).group(1))
    beta = alpha / ratio
    dose = 2.0  # sim-spatial's default --radiation-dose-gy
    expected = (1.0 - math.exp(-(alpha * dose + beta * dose * dose))) * 100.0
    rad = next(r for r in d["rows"] if r["treatment"] == "Radiation")
    assert abs(rad["overall_kill_pct"] - expected) < 2.0, (
        f"the binary killed {rad['overall_kill_pct']:.1f}% where LQ at "
        f"{dose} Gy predicts {expected:.1f}%; the DNA channel and the crate "
        "constants have drifted apart")
    # Non-vacuous: the tolerance must not admit an unrelated arm.
    for other in ("PDT", "SDT", "RSL3"):
        o = next(r for r in d["rows"] if r["treatment"] == other)
        assert abs(o["overall_kill_pct"] - expected) > 2.0, (
            f"{other} also lands within tolerance of the LQ expectation, so "
            "this check does not identify radiation")


def test_the_caveats_are_the_ones_that_actually_bind(md):
    """Each names a limit that would change a number if lifted, so none is
    decorative."""
    for frag in ("ferroptosis channel is OFF",
                 "Single fraction",
                 "no oxygen field",
                 "not like for like",
                 "identifiability-report.md"):
        assert frag in md, f"the page no longer states: {frag}"


def test_the_dropped_bin_is_declared_and_small(d, md):
    """The deepest bin is noise and the page says so. An earlier figure
    annotation pointed AT it and reported 100% kill for two modalities that
    disagree everywhere else."""
    assert d["dropped_last_bins"] >= 1
    assert "poles hold a handful of cells" in md
    for tx, series in d["binned"].items():
        if len(series) < 3:
            continue
        tail = series[-1]["n_cells"]
        typical = sorted(b["n_cells"] for b in series)[len(series) // 2]
        assert tail < typical, (
            f"{tx}'s deepest bin has {tail} cells against a median of "
            f"{typical}; it is no longer sparse and dropping it needs a reason")
