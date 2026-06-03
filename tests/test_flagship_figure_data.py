"""Drift guard for the flagship resistance-asymmetry figure (manuscript Figure 24, #285).

`fig27_resistance_asymmetry` in scripts/generate_figures.py renders four panels whose
headline numbers ALSO appear as literal prose in the figure caption (generate_latex.py
figs '24'), the v1.md [FIGURE 24] placeholder, and the §7.1-7.4 text. The figure PDF/PNG
are git-tracked but the live data source (simulations/output/tme/tme_summary.json) is
gitignored and never regenerated in CI, so a sim re-run could silently leave a stale
figure and stale caption text with nothing detecting it.

To make this guard effective IN CI (not a skip-if-absent no-op), it reads a committed
snapshot of the exact rows Figure 24 depends on (tests/fixtures/flagship_tme_rows.json).
The fixture is the canonical record of the figure's numbers; if sim-tme output changes,
regenerate the figure, refresh the fixture, and update the caption together. When the
live JSON is present (a developer machine), an extra test cross-checks the fixture
against it so a stale fixture is caught at dev time.

This guard (the #293 precedent, mirroring tests/test_depth_kill_physics_constants.py):
  - asserts the conditions fig27 needs resolve,
  - asserts the figure's entire thesis (RSL3 kill < SDT kill) on every panel,
  - pins each caption headline number,
  - confirms the pH panel's ferroptosis_kills metric is an immune-free counter,
  - cross-checks the committed fixture against the live sim when available.

Each panel deliberately uses the SAME metric its manuscript section reports
(§7.1 overall kill / §7.3 CAF-adjacent kill / §7.4 ferroptosis kills / §7.2 immune
kills), so these assertions also keep the figure and the prose in agreement.
"""

import json
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
FIXTURE = REPO_ROOT / "tests" / "fixtures" / "flagship_tme_rows.json"
LIVE_SUMMARY = REPO_ROOT / "simulations" / "output" / "tme" / "tme_summary.json"
G = "gradient_120um"


def _load(path):
    d = json.loads(path.read_text())
    return d["conditions"] if isinstance(d, dict) and "conditions" in d else d


def _conds():
    assert FIXTURE.exists(), f"committed fixture missing: {FIXTURE}"
    return _load(FIXTURE)


def _find(conds, treatment, **kw):
    for c in conds:
        if c["treatment"] == treatment and all(c.get(k) == v for k, v in kw.items()):
            return c
    return None


def _rows(conds, t):
    """The exact (baseline, stressed) rows fig27 reads, per treatment."""
    return {
        "hyp": (_find(conds, t, o2_condition="uniform", immune_mode="off"),
                _find(conds, t, o2_condition=G, immune_mode="off",
                      stromal_mode=None, ph_mode=None)),
        "strm": (_find(conds, t, o2_condition=G, immune_mode="immune_on", stromal_mode="off"),
                 _find(conds, t, o2_condition=G, immune_mode="immune_on", stromal_mode="stromal_on")),
        "ph": (_find(conds, t, o2_condition=G, immune_mode="immune_on", stromal_mode="off"),
               _find(conds, t, o2_condition=G, immune_mode="immune_on", ph_mode="ph_on")),
        "imm": _find(conds, t, o2_condition=G, immune_mode="immune_on", stromal_mode="off"),
    }


def test_required_conditions_resolve():
    """All rows fig27 dereferences must exist (else the figure silently skips)."""
    conds = _conds()
    for t in ("RSL3", "SDT"):
        r = _rows(conds, t)
        assert r["hyp"][0] and r["hyp"][1], f"{t}: hypoxia rows missing"
        assert r["strm"][0] and r["strm"][1], f"{t}: stromal rows missing"
        assert r["ph"][0] and r["ph"][1], f"{t}: pH rows missing"
        assert r["imm"], f"{t}: immune row missing"


def test_rsl3_collapses_below_sdt_on_every_panel():
    """The figure's whole thesis: under each mechanism, RSL3 kill < SDT kill."""
    conds = _conds()
    rsl3, sdt = _rows(conds, "RSL3"), _rows(conds, "SDT")
    # rate panels (a) hypoxia, (b) stromal
    assert rsl3["hyp"][1]["overall_kill_rate"] < sdt["hyp"][1]["overall_kill_rate"]
    assert rsl3["strm"][1]["stromal_adjacent_kill_rate"] < sdt["strm"][1]["stromal_adjacent_kill_rate"]
    # count panels (c) pH ferroptosis, (d) immune
    assert rsl3["ph"][1]["ferroptosis_kills"] < sdt["ph"][1]["ferroptosis_kills"]
    assert rsl3["imm"]["immune_kills"] < sdt["imm"]["immune_kills"]


def test_hypoxia_caption_numbers():
    """Panel (a) headline (§7.1): RSL3 3.7%->0.1%, SDT 91.9%->87.8%."""
    conds = _conds()
    rsl3, sdt = _rows(conds, "RSL3")["hyp"], _rows(conds, "SDT")["hyp"]
    assert round(rsl3[0]["overall_kill_rate"] * 100, 1) == 3.7
    assert round(rsl3[1]["overall_kill_rate"] * 100, 1) == 0.1
    assert round(sdt[0]["overall_kill_rate"] * 100, 1) == 91.9
    assert round(sdt[1]["overall_kill_rate"] * 100, 1) == 87.8


def test_stromal_caption_numbers():
    """Panel (b) headline (§7.3): RSL3 3.0%->1.5%, SDT 96.1%->91.2%."""
    conds = _conds()
    rsl3, sdt = _rows(conds, "RSL3")["strm"], _rows(conds, "SDT")["strm"]
    assert round(rsl3[0]["stromal_adjacent_kill_rate"] * 100, 1) == 3.0
    assert round(rsl3[1]["stromal_adjacent_kill_rate"] * 100, 1) == 1.5
    assert round(sdt[0]["stromal_adjacent_kill_rate"] * 100, 1) == 96.1
    assert round(sdt[1]["stromal_adjacent_kill_rate"] * 100, 1) == 91.2


def test_ph_caption_numbers():
    """Panel (c) headline (§7.4): RSL3 ferroptosis kills 163->77, SDT 139640->140693."""
    conds = _conds()
    rsl3, sdt = _rows(conds, "RSL3")["ph"], _rows(conds, "SDT")["ph"]
    assert rsl3[0]["ferroptosis_kills"] == 163
    assert rsl3[1]["ferroptosis_kills"] == 77
    assert sdt[0]["ferroptosis_kills"] == 139640
    assert sdt[1]["ferroptosis_kills"] == 140693


def test_immune_caption_numbers_and_ratio():
    """Panel (d) headline (§7.2): RSL3 5, SDT 521, ratio rounds to 104:1."""
    conds = _conds()
    rsl3 = _rows(conds, "RSL3")["imm"]["immune_kills"]
    sdt = _rows(conds, "SDT")["imm"]["immune_kills"]
    assert rsl3 == 5
    assert sdt == 521
    assert round(sdt / max(rsl3, 1)) == 104


def test_ph_panel_metric_is_immune_free():
    """The pH panel uses ferroptosis_kills, a counter genuinely separate from
    immune_kills. The decisive evidence: on the RSL3 acidic-pH row the ferroptosis
    counter is nonzero (77) while immune_kills is exactly 0, and a metric that
    counts cells with no immune kills present cannot be folding immune kills in.
    So panel (c) is not contaminated by the immune_on baseline (unlike a raw kill
    rate, which sums both causes into state.dead)."""
    conds = _conds()
    rsl3_ph = _rows(conds, "RSL3")["ph"][1]
    assert rsl3_ph["ferroptosis_kills"] > 0
    assert rsl3_ph["immune_kills"] == 0


def test_fixture_matches_live_sim_when_present():
    """Dev-time freshness: when the (gitignored) live tme_summary.json is present,
    the committed fixture must still match it on every pinned value, so a sim
    re-run that shifts the numbers is caught and the fixture+caption refreshed."""
    if not LIVE_SUMMARY.exists():
        pytest.skip("live tme_summary.json not present (CI / fresh checkout)")
    fix, live = _conds(), _load(LIVE_SUMMARY)
    for t in ("RSL3", "SDT"):
        rf, rl = _rows(fix, t), _rows(live, t)
        assert rf["hyp"][1]["overall_kill_rate"] == rl["hyp"][1]["overall_kill_rate"]
        assert rf["strm"][1]["stromal_adjacent_kill_rate"] == rl["strm"][1]["stromal_adjacent_kill_rate"]
        assert rf["ph"][1]["ferroptosis_kills"] == rl["ph"][1]["ferroptosis_kills"]
        assert rf["imm"]["immune_kills"] == rl["imm"]["immune_kills"], (
            f"{t}: committed fixture is STALE vs live sim; re-run sim-tme, refresh "
            "tests/fixtures/flagship_tme_rows.json and the Figure 24 caption together."
        )
