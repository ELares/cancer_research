"""
Drift guard for the depth-kill figure (manuscript Figure 8, #285).

Panel (b) of `fig8_simulation_by_treatment` in scripts/generate_figures.py is a
Python re-implementation of the Rust penetration physics, with the attenuation
constants HARDCODED. Panel (a) of the same figure comes from the `sim-spatial`
binary, which reads those constants from `ferroptosis-core/src/params.rs`
(`SpatialParams::default`). If the two drift apart, the figure's two panels
silently disagree.

This test pins the figure's hardcoded constants to the Rust defaults: if anyone
retunes `pdt_mu_eff` / `sdt_alpha` / `sdt_freq_mhz` in params.rs (or edits the
figure's copies), it fails, signalling that BOTH must be updated together (and
the depth_kill_curves.csv regenerated with default sim-spatial flags).
"""

import re
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
PARAMS_RS = REPO_ROOT / "simulations" / "ferroptosis-core" / "src" / "params.rs"
GENERATE_FIGURES = REPO_ROOT / "scripts" / "generate_figures.py"

# (params.rs field name, figure constant name) -> both must hold this value.
COUPLED = [
    ("pdt_mu_eff", "PDT_MU_EFF_PER_MM"),
    ("sdt_alpha", "SDT_ALPHA_DB_CM_MHZ"),
    ("sdt_freq_mhz", "SDT_FREQ_MHZ"),
]


def _rust_default(field: str) -> float:
    text = PARAMS_RS.read_text()
    m = re.search(rf"\b{re.escape(field)}:\s*([0-9]+\.?[0-9]*)", text)
    assert m, f"{field} not found in {PARAMS_RS}"
    return float(m.group(1))


def _figure_const(name: str) -> float:
    text = GENERATE_FIGURES.read_text()
    m = re.search(rf"\b{re.escape(name)}\s*=\s*([0-9]+\.?[0-9]*)", text)
    assert m, f"{name} not found in {GENERATE_FIGURES} (fig8 panel b)"
    return float(m.group(1))


@pytest.mark.parametrize("rust_field,fig_const", COUPLED)
def test_figure_physics_matches_rust_defaults(rust_field, fig_const):
    """fig8 panel (b)'s hardcoded constant must equal the params.rs default."""
    rust_val = _rust_default(rust_field)
    fig_val = _figure_const(fig_const)
    assert rust_val == fig_val, (
        f"Depth-kill figure drift: params.rs {rust_field}={rust_val} but "
        f"generate_figures.py {fig_const}={fig_val}. Update fig8 panel (b) to "
        f"match params.rs (and regenerate depth_kill_curves.csv with default "
        f"sim-spatial flags), or vice versa."
    )
