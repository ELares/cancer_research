"""The manuscript must describe the per-cell parameter draw the code performs.

Section 5.3 said each cell's parameters are "drawn from a +/-20% uniform random
variation around the population mean". `gen_cell` in
`ferroptosis-core/src/cell.rs` draws each parameter from `Normal::new(mean, sd)`
with per-parameter CV between 10% and 40%, floored at a physiological minimum.

Three ways that was wrong, and the third is why this file exists:

  family   uniform vs Gaussian;
  width    "+/-20%" against 95% spans of +/-20% to +/-78%;
  SUPPORT  bounded against unbounded.

The width line is itself a correction. It first read 12% to 25%, which is the
range over the Glycolytic block ALONE -- the first phenotype in the file, and so
the first one read. Across all four it is 10% to 40%, and the widest draw is
nearly twice what the narrow reading suggested. `test_the_stated_cv_range_matches
_the_code` computes the range from the source rather than restating it, which is
what caught this; the docstring beside it had to be fixed by hand.

A bounded draw makes a cell outside the box impossible, so a death rate reported
as zero would be exactly zero and no sample size could ever find an event. The
Gaussian draw leaves every threshold a positive-probability tail, so a reported
zero is an upper bound set by n. The distinction decides whether rare-event work
at large n is meaningful at all -- it was found during pre-flight for exactly
such a run, where the manuscript's description would have made the exercise
pointless.
"""

import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
CELL_RS = REPO_ROOT / "simulations" / "ferroptosis-core" / "src" / "cell.rs"
MS = REPO_ROOT / "article" / "drafts" / "v1.md"


def _draws():
    """(name, mean, sd) for every normal draw in gen_cell's phenotype blocks."""
    src = CELL_RS.read_text()
    body = src[src.index("pub fn gen_cell"):]
    return [(n, float(m), float(s))
            for n, m, s in re.findall(r"(\w+):\s*norm\(rng,\s*([\d.]+),\s*([\d.]+)\)", body)]


def test_the_code_still_draws_from_normals():
    """The premise. If this ever becomes uniform, the prose must change back."""
    assert "Normal::new" in CELL_RS.read_text(), (
        "gen_cell no longer uses a normal distribution; the manuscript's "
        "unbounded-support argument depends on it")
    assert _draws(), "no norm(rng, mean, sd) draws found in gen_cell"


def test_the_manuscript_does_not_call_the_draw_uniform():
    txt = MS.read_text()
    para = txt[txt.index("Each cell in the simulation receives its own parameter set"):]
    para = para[:para.index("\n\n")]
    assert "uniform random variation" not in para or "earlier version" in para, (
        "the manuscript describes the per-cell draw as uniform; the code draws "
        "from normals, and the difference decides whether a reported 0% death "
        "rate is exactly zero or merely below the sample's resolution")
    assert "normal" in para.lower()
    assert "UNBOUNDED" in para or "unbounded" in para, (
        "the consequential property -- unbounded support -- is not stated")


def test_the_stated_cv_range_matches_the_code():
    cvs = sorted({round(100 * sd / mean) for _, mean, sd in _draws() if mean})
    txt = MS.read_text()
    para = txt[txt.index("Each cell in the simulation receives its own parameter set"):]
    para = para[:para.index("\n\n")]
    assert f"{min(cvs)}% and {max(cvs)}%" in para, (
        f"the manuscript's CV range does not match the code's {min(cvs)}-{max(cvs)}%")


def test_a_reported_zero_is_described_as_an_upper_bound():
    """The operational consequence, which is what a reader needs."""
    txt = MS.read_text()
    para = txt[txt.index("Each cell in the simulation receives its own parameter set"):]
    para = para[:para.index("\n\n")]
    assert "upper bound" in para, (
        "the manuscript does not say that a zero death rate is an upper bound "
        "set by the sample size rather than a measurement of zero")
