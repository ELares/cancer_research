"""Drift guards for the remaining quantitative simulation figures (#295 item 5).

Figs 8 (depth-kill) and 24 (flagship) already have data-drift guards
(test_depth_kill_physics_constants.py, test_flagship_figure_data.py). This file
extends the same pattern to the other three quantitative sim figures whose
captions carry hardcoded headline numbers fed by a gitignored sim output:

  - Figure 21 — hypoxia kill-collapse        (fig24_hypoxia_killcurve, sim-tme)
  - Figure 22 — Bliss dual-pathway synergy   (fig25_bliss_synergy, sim-combo-mech)
  - Figure 23 — vulnerability window         (fig26_vulnerability_window, sim-window)

Each figure's generator reads a gitignored JSON (output/{tme,combo-mech,window}/)
that CI never regenerates, and the caption numbers live as literal prose in
generate_latex.py + v1.md. So a sim re-run could silently leave a stale figure
and stale caption. To be effective IN CI, each test reads a COMMITTED fixture
(the exact rows the generator uses, tests/fixtures/*.json) and pins the caption
headline numbers; a separate test cross-checks the fixture against the live JSON
when present (dev machine) so a stale fixture is caught at dev time.

Pinned values mirror each generator's own extraction (see scripts/generate_figures.py).
"""

import json
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
FIX = REPO_ROOT / "tests" / "fixtures"
SIM = REPO_ROOT / "simulations" / "output"


def _load(path):
    return json.loads(path.read_text())


def _fixture(name):
    p = FIX / name
    assert p.exists(), f"committed fixture missing: {p}"
    return _load(p)


# ── Figure 21: hypoxia kill-collapse (overall kill, immune off) ──────────────
LAMBDAS = (80, 100, 120, 150)


def _hypoxia_kills(conds):
    """Mirror fig24_hypoxia_killcurve: normoxic = uniform; hypoxic = mean over
    the gradient_{80,100,120,150}um conditions; per treatment, immune off."""
    def ov(t, o2):
        for c in conds:
            if c["treatment"] == t and c.get("o2_condition") == o2 and c.get("immune_mode") == "off":
                return c["overall_kill_rate"]
        raise AssertionError(f"missing {t}/{o2}")
    out = {}
    for t in ("RSL3", "SDT"):
        norm = ov(t, "uniform") * 100
        hyp = sum(ov(t, f"gradient_{l}um") for l in LAMBDAS) / len(LAMBDAS) * 100
        out[t] = (norm, hyp)
    return out


def test_fig21_hypoxia_caption_numbers():
    """Caption: RSL3 3.7%->~0.1%, SDT 91.9%->87.8%."""
    k = _hypoxia_kills(_fixture("hypoxia_killcurve_rows.json")["conditions"])
    assert round(k["RSL3"][0], 1) == 3.7, k["RSL3"]
    assert round(k["RSL3"][1], 1) == 0.1, k["RSL3"]
    assert round(k["SDT"][0], 1) == 91.9, k["SDT"]
    assert round(k["SDT"][1], 1) == 87.8, k["SDT"]


def test_fig21_rsl3_collapses_while_sdt_holds():
    """The figure's thesis: hypoxia roughly kills RSL3 (>30x drop) but SDT holds."""
    k = _hypoxia_kills(_fixture("hypoxia_killcurve_rows.json")["conditions"])
    assert k["RSL3"][1] < k["RSL3"][0] / 10.0, "RSL3 should collapse under hypoxia"
    assert k["SDT"][1] > 80.0, "SDT should hold >80% under hypoxia"


# ── Figure 22: Bliss dual-pathway synergy (RSL3 + FSP1i) ─────────────────────
def test_fig22_bliss_caption_numbers():
    """Caption: RSL3+FSP1i 84.1%, Bliss-expected 42.2%, 1.99x synergy."""
    rf = _fixture("bliss_synergy.json")["rsl3_fsp1i"]
    assert round(rf["rate_combo"] * 100, 1) == 84.1, rf["rate_combo"]
    assert round(rf["bliss_prediction"] * 100, 1) == 42.2, rf["bliss_prediction"]
    assert round(rf["synergy_score"], 2) == 1.99, rf["synergy_score"]


def test_fig22_combination_is_synergistic():
    """Observed combo must beat the Bliss-independent expectation (synergy > 1)."""
    rf = _fixture("bliss_synergy.json")["rsl3_fsp1i"]
    assert rf["rate_combo"] > rf["bliss_prediction"], "combo should exceed Bliss expectation"
    assert rf["synergy_score"] > 1.0, "synergy score should be > 1"


# ── Figure 23: vulnerability window (RSL3 closes, SDT stays open) ─────────────
def _window_death(rows, t, day):
    for r in rows:
        if r["treatment"] == t and r["timepoint_days"] == day:
            return r["death_rate"] * 100
    raise AssertionError(f"missing {t} @ day {day}")


def test_fig23_window_caption_numbers():
    """Caption: RSL3 42.4%->1.4% by day 3 and ~0 by day 7; SDT ~99.5% by day 28."""
    rows = _fixture("vulnerability_window.json")["rows"]
    assert round(_window_death(rows, "RSL3", 0.0), 1) == 42.4
    assert round(_window_death(rows, "RSL3", 3.0), 1) == 1.4
    assert _window_death(rows, "RSL3", 7.0) < 0.5, "RSL3 ~0 by day 7"
    assert round(_window_death(rows, "SDT", 28.0), 1) == 99.5


def test_fig23_rsl3_window_closes_sdt_stays_open():
    """RSL3 kill collapses across the window while SDT stays high to day 28."""
    rows = _fixture("vulnerability_window.json")["rows"]
    assert _window_death(rows, "RSL3", 0.0) > 10 * _window_death(rows, "RSL3", 7.0)
    assert _window_death(rows, "SDT", 28.0) > 95.0


# ── Dev-time freshness: committed fixtures must match the live sim ───────────
def test_fig21_fixture_matches_live():
    p = SIM / "tme" / "tme_summary.json"
    if not p.exists():
        pytest.skip("live tme_summary.json not present (CI / fresh checkout)")
    live = _load(p)
    live = live["conditions"] if isinstance(live, dict) and "conditions" in live else live
    assert _hypoxia_kills(_fixture("hypoxia_killcurve_rows.json")["conditions"]) == _hypoxia_kills(live), (
        "hypoxia fixture STALE vs live sim — re-run sim-tme, refresh the fixture + caption"
    )


def test_fig22_fixture_matches_live():
    p = SIM / "combo-mech" / "combo_summary.json"
    if not p.exists():
        pytest.skip("live combo_summary.json not present")
    live = _load(p)
    combos = live["combinations"] if isinstance(live, dict) and "combinations" in live else live
    rf_live = next(c for c in combos if {c["drug_a"], c["drug_b"]} == {"RSL3", "FSP1i"})
    rf_fix = _fixture("bliss_synergy.json")["rsl3_fsp1i"]
    for key in ("rate_combo", "bliss_prediction", "synergy_score"):
        assert rf_fix[key] == rf_live[key], f"Bliss fixture STALE vs live sim on {key}"


def test_fig23_fixture_matches_live():
    p = SIM / "window" / "vulnerability_window.json"
    if not p.exists():
        pytest.skip("live vulnerability_window.json not present")
    live = _load(p)
    rows = _fixture("vulnerability_window.json")["rows"]
    for r in rows:
        match = next(
            (x for x in live if x["treatment"] == r["treatment"] and x["timepoint_days"] == r["timepoint_days"]),
            None,
        )
        assert match is not None, f"live missing {r['treatment']} @ {r['timepoint_days']}"
        # death_rate drives panel (a); mean_gpx4 drives panel (b)'s twin axis
        # (the "RSL3 collapse tracks GPX4 recovery" mechanism), so guard both.
        for field in ("death_rate", "mean_gpx4"):
            assert r[field] == match[field], (
                f"window fixture STALE vs live sim on {field} at "
                f"{r['treatment']} day {r['timepoint_days']}"
            )
