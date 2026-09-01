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
    # The ~104:1 headline as a band rather than `round(...) == 104` (a tautology
    # of the exact 5/521 goldens above — 521/5 = 104.2 always rounds to 104). The
    # exact goldens still gate this, so on a 2D fixture re-baseline update them and
    # this band together; the band documents the intended ~104:1 claim explicitly.
    assert 100 <= sdt / max(rsl3, 1) <= 108


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


FIG27 = REPO_ROOT / "article/figures/fig27_resistance_asymmetry.pdf"


def test_fig27_draws_each_panels_numbers_in_that_panel():
    """This figure is the manuscript's flagship and NOTHING opened it (#793).

    Every other guard in this file validates the INPUT rows -- "asserts the
    conditions fig27 needs resolve" -- so nothing compared what the figure
    DRAWS against anything. And a check that a value is merely PRESENT in a
    four-panel figure cannot say which panel it landed in, while the panel IS
    the claim: (a) hypoxia, (b) stromal, (c) pH, (d) immune. A number in the
    wrong panel inverts the comparison the figure exists to make.

    THE PANEL IS THE DRAWN AXES RECTANGLE. Two earlier attempts split on
    geometry that only looks like a boundary: the midpoint between panel
    TITLES falls inside the left panel, because a title sits at its panel's
    left edge; and the widest horizontal gap is INSIDE panel (a), between two
    bar labels. matplotlib draws each panel as an axes rectangle, and that
    rectangle is the panel -- the same lesson this repo already learned about
    clip rectangles.

    `stromal_adjacent_kill_rate`, which #793 first called a missing field and
    then correctly retracted, is carried by the stromal rows of the fixture, so
    panel (b) is reachable with no `sim-tme` run. I re-made that same error
    while writing this, by checking row 0 instead of the union of all rows.
    """
    fitz = pytest.importorskip("fitz")
    doc = fitz.open(FIG27)
    page = doc[0]
    spans = [(round(s["bbox"][0], 1), round(s["bbox"][1], 1), s["text"].strip())
             for b in page.get_text("dict")["blocks"]
             for l in b.get("lines", []) for s in l["spans"] if s["text"].strip()]
    rects = [d["rect"] for d in page.get_drawings()
             if d["rect"].width > 150 and d["rect"].height > 100
             and d["rect"].width < page.rect.width * 0.8]
    doc.close()

    assert len(rects) == 4, (
        f"fig27 draws {len(rects)} panel rectangles; it is a 2x2 figure and a "
        "change in that count changes what it claims")
    rects.sort(key=lambda r: (round(r.y0), round(r.x0)))
    letter_of = {id(r): L for r, L in zip(rects, "abcd")}

    def panel_of(x, y):
        for r in rects:
            if r.x0 - 42 <= x <= r.x1 + 8 and r.y0 - 26 <= y <= r.y1 + 30:
                return letter_of[id(r)]
        return None

    conds = _conds()
    checks = []
    for tx in ("SDT", "RSL3"):
        base, stressed = _rows(conds, tx)["hyp"]
        for r in (base, stressed):
            checks.append(("a", f"{r['overall_kill_rate'] * 100:.1f}"))
    for c in conds:
        if "stromal_adjacent_kill_rate" in c:
            checks.append(("b", f"{c['stromal_adjacent_kill_rate'] * 100:.1f}"))
        if c.get("ph_mode") and c.get("ferroptosis_kills") is not None:
            checks.append(("c", str(int(c["ferroptosis_kills"]))))
        if c.get("immune_mode") == "on" and c.get("immune_kills"):
            checks.append(("d", str(int(c["immune_kills"]))))

    # A short value can legitimately belong to two panels -- "0.1" is panel
    # (a)'s RSL3 hypoxic kill AND a stromal rate -- so only values expected in
    # exactly ONE panel can be bound. Checking the ambiguous ones would report
    # a correct figure as wrong, which is the failure mode a guard must not
    # have.
    from collections import defaultdict
    expected = defaultdict(set)
    for letter, value in checks:
        expected[value].add(letter)
    unambiguous = [(l, v) for l, v in checks if len(expected[v]) == 1]

    misplaced, drawn = [], 0
    for letter, value in unambiguous:
        located = {p for x, y, t in spans if t.rstrip("%") == value
                   for p in [panel_of(x, y)] if p}
        if not located:
            continue
        drawn += 1
        if letter not in located:
            misplaced.append(
                f"{value} drawn in panel(s) {sorted(located)}, expected {letter}")
    assert not misplaced, (
        "fig27 draws values outside the panel whose data they are: "
        + "; ".join(misplaced[:5]))
    assert drawn >= 6, (
        f"only {drawn} of {len(unambiguous)} unambiguous fixture values appear "
        "in the figure at all, so the panel binding above is checking almost "
        "nothing")
