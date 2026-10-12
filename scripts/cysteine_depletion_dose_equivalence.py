#!/usr/bin/env python3
"""What is a cyst(e)inase worth in erastin, and what will the engine not say?

WHY
---
`cysteine_depletion_feasibility.py` established that the engineered-bacteria arm
of PMID 39696900 needs no new axis: `erastin_xc_inhib` is a supply fraction, an
enzyme delivers a residual concentration, and the transporter's measured Km
(179.4 uM, PMID 34287957) against a formulation reference (207.7 uM) maps one to
the other. This asks the next question. The engine already carries a CALIBRATED
erastin dose-response -- Hill `dose^h/(dose^h+K)` fitted to 785 CTRPv2 curves --
so the two maps compose, and a depletion depth can be quoted as the erastin dose
it is worth.

WHAT SURVIVES WITHOUT A CASCADE, AND IT IS MOST OF THE PAGE
-----------------------------------------------------------
Three results need no simulation at all, which is why they are the headline:

1. DOSE EQUIVALENCE. Composing the two maps is arithmetic over published
   constants. The whole operating range of the enzyme -- a quarter of the
   cystine removed up to essentially all of it -- lands inside a ~2.5x erastin
   dose window, every point of it within the calibrated range. That narrowness
   is the Hill coefficient of 6, not a property of the enzyme.

2. A STRUCTURAL CEILING. The enzyme cannot cut supply below zero, so through
   this mechanism it cannot exceed a saturating erastin dose. That is topology,
   not a measurement, and it is the useful kind of negative: the advantage the
   paper claims -- "efficacy, persistence, and precision" -- has to come from
   targeting and persistence, because mechanistic depth is already spent.

3. THE ESCAPE IS AN EXACT NO-OP SWITCH. `xc_cystine_factor()` is
   `1 - erastin_xc_inhib * (1 - transsulfuration_floor)`, so at a floor of 1 the
   factor is EXACTLY 1 for every inhibition: a fully competent transsulfuration
   pathway makes the entire cystine arm a precise no-op, at any depletion depth
   and at any cascade. More generally the block delivers only `1 - floor` of its
   nominal cut, linearly. This is the quantitative reason the second bacterial
   paper (PMID 41704298) had to add H2S to inhibit the bypass, and it needs no
   calibration to state.

WHAT THIS REFUSES TO SAY, AND WHY THE REFUSAL IS THE RESULT
-----------------------------------------------------------
No absolute kill rate. Not "not yet measured" -- REFUSED, because both available
cascades fail, in opposite directions, and the repo already has the gate that
catches it:

  * At the IN-VITRO FITTED cascade (lp_propagation 1.0, lp_rate 0.4, the one that
    reproduces the real CTRPv2 erastin curve) the arm kills up to ~52% -- but the
    UNTREATED Glycolytic control dies at 8.9%, four times the manuscript's own
    2% baseline-viability constraint, and three of four phenotypes die at 89-100%
    with no treatment at all. `headline_at_fitted.py` already found and published
    this for every other headline; it holds here too.
  * At the ADMISSIBLE in-vivo default cascade every control is ~0% and a COMPLETE
    cystine block kills exactly 0.00% at every depth. The arm does nothing.

So the magnitude of this arm is not predictable by this engine today, for a
reason already on record rather than a new one. A first draft of this analysis
reported the fitted-cascade kill table as a result; running the untreated
controls is what stopped it, and that check is the only reason this page says
what it says.

OFFLINE CONTRACT
----------------
The arithmetic needs nothing but the committed constants and runs in CI. The
cascade measurements need the compiled `ferroptosis_core` extension, which CI
does not build, so they are taken with `--measure` locally and committed into
the JSON. Guards that read them skip where the extension is absent.

Usage:
    python scripts/cysteine_depletion_dose_equivalence.py
    python scripts/cysteine_depletion_dose_equivalence.py --measure
    python scripts/cysteine_depletion_dose_equivalence.py --render-only
"""

import argparse
import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from config import PROJECT_ROOT  # noqa: E402

# IMPORTED, not restated. The depletion->supply map belongs to #909 and this
# page composes on top of it; a second copy beside the real one is how the two
# drift apart, which is the precedent `atlas_landscape.PHARMACOLOGICAL` was
# hoisted to module level for. They agreed across 99 depletions when this page
# carried its own copy -- agreement today is exactly what makes a duplicate
# comfortable to leave in place.
from cysteine_depletion_feasibility import (  # noqa: E402
    inhibition_from_depletion,
    saturation,
)

CAL = PROJECT_ROOT / "analysis" / "calibration"
CONST_SRC = CAL / "cysteine_depletion_constants.csv"
ERASTIN_SRC = CAL / "erastin-calibration.json"
# The repo's own baseline-viability record. A module constant rather than an
# inline path so a guard can point it elsewhere and prove the constraint is
# READ: a value-level check cannot tell a read from a hardcoded literal that
# happens to be right today, and a mutation planting the correct 0.02 survived
# exactly that way.
ADMISSIBILITY_SRC = PROJECT_ROOT / "analysis" / "headline-at-fitted-cascade.json"
OUT_MD = CAL / "cysteine-depletion-dose-equivalence.md"
OUT_JSON = CAL / "cysteine-depletion-dose-equivalence.json"

DEPLETIONS = [0.25, 0.50, 0.75, 0.90, 0.95, 0.99]
FLOORS = [0.0, 0.25, 0.50, 0.75, 1.0]
# The repo's own baseline-viability constraint, read from the committed
# admissibility record rather than typed here.
PHENOTYPES = ["Glycolytic", "OXPHOS", "Persister", "PersisterNrf2"]
# The order the refusal ARGUES in: the responsive-but-inadmissible cascade
# first, then the admissible-but-unresponsive one. Established here rather
# than inherited from dict iteration -- `json.dumps(sort_keys=True)` reorders
# a stored dict, so a table walking `.items()` publishes whatever
# serialisation happened to produce. That is the third instance of this
# defect in this repo and the second of mine today.
CASCADE_ORDER = ["in_vitro_fitted", "in_vivo_default"]


class MissingInputError(RuntimeError):
    """Raised rather than rendering a page from absent or partial inputs."""


def _load(path: Path) -> dict:
    if not path.exists():
        raise MissingInputError(f"missing committed input: {path}")
    return json.loads(path.read_text())


def load_inputs() -> dict:
    if not CONST_SRC.exists():
        raise MissingInputError(f"missing committed constants: {CONST_SRC}")
    rows = {r["constant"]: r for r in csv.DictReader(CONST_SRC.open())}
    for k in ("xct_cystine_km", "rpmi1640_cystine_2hcl_mass", "cystine_2hcl_molar_mass"):
        if k not in rows:
            raise MissingInputError(f"{CONST_SRC} lacks {k}")
    er = _load(ERASTIN_SRC)
    fit = er.get("calibrated_params") or {}
    for k in ("k_erastin", "hill"):
        if k not in fit:
            raise MissingInputError(f"{ERASTIN_SRC} lacks calibrated_params.{k}")
    # The admissibility constraint is the repo's, not this page's.
    adm = _load(ADMISSIBILITY_SRC)
    constraint = (adm.get("ctrpv2_point", {}).get("admissibility", {}) or {}).get("constraint")
    if constraint is None:
        raise MissingInputError(
            "headline-at-fitted-cascade.json carries no admissibility constraint; "
            "this page's refusal rests on it and must not invent one")
    return {
        "km_uM": float(rows["xct_cystine_km"]["value"]),
        "reference_uM": (float(rows["rpmi1640_cystine_2hcl_mass"]["value"])
                         / float(rows["cystine_2hcl_molar_mass"]["value"]) * 1000.0),
        "k_erastin_uM": float(fit["k_erastin"]),
        "hill": float(fit["hill"]),
        "erastin_fit_rmse": er.get("fit_rmse"),
        # The committed saturating-erastin effect, so the structural ceiling can
        # be CHECKED rather than only argued. At the top of the committed grid
        # the Hill map is already at ~99.99% inhibition, so that increment IS a
        # complete-block increment measured independently of this page.
        "committed_baseline_death": er.get("control_baseline_death_at_fitted_cascade"),
        "committed_top_dose_increment": er.get("erastin_increment_top_dose"),
        "erastin_fit_phenotype": er.get("phenotype"),
        "erastin_dose_grid_uM": er.get("dose_grid_um"),
        "erastin_cascade": {k: v for k, v in fit.items()
                            if k in ("lp_propagation", "lp_rate")},
        "viability_constraint": float(constraint),
    }


def erastin_dose_for_inhibition(i: float, k: float, h: float):
    """Invert the production Hill fit: `dose^h / (dose^h + k^h) = i`.

    `k` is RAISED TO THE POWER h, so `k` is the half-maximal dose and the
    inverse is `k * (i/(1-i))^(1/h)`.

    AN EARLIER VERSION INVERTED `dose^h/(dose^h + k)` -- k used raw -- giving
    `(k*i/(1-i))^(1/h)`, which is low by a factor of `k^(1-1/h)`: 4.451 at the
    committed k=6, h=6. It put the half-maximal dose at 1.35 uM when the
    half-maximal dose is k=6 uM BY DEFINITION, which is the check that should
    have been run first. The form was taken from a docstring summary instead of
    from `calibrate_erastin.dose_to_inhib`, and the round-trip guard could not
    catch it because it inverted this function with this function -- the
    guard-computes-its-own-expectation shape. The guard now round-trips through
    the PRODUCTION function.

    Returns None at the endpoints rather than raising: a zero block needs no
    dose and a total block needs an infinite one, and reporting either as a
    number would be a lie about an asymptote.
    """
    if i <= 0.0 or i >= 1.0:
        return None
    return k * (i / (1.0 - i)) ** (1.0 / h)


def escape_delivery(floor: float) -> float:
    """The fraction of its nominal cut the block delivers at a given floor.

    Exact and cascade-free: `xc_cystine_factor = 1 - inhib*(1-floor)`, so the
    supply actually removed is `inhib*(1-floor)` and the delivered fraction is
    `1-floor`, independent of `inhib` entirely.
    """
    return 1.0 - min(max(floor, 0.0), 1.0)


def cystine_factor(inhib: float, floor: float) -> float:
    """The engine's own expression, restated so the identity can be checked."""
    return 1.0 - min(max(inhib, 0.0), 1.0) * (1.0 - min(max(floor, 0.0), 1.0))


def measure_cascades(inp: dict) -> dict:
    """Run both cascades. Needs the compiled extension; `--measure` only."""
    import ferroptosis_core as fc  # noqa: PLC0415

    fitted = dict(inp["erastin_cascade"])
    out = {}
    for name, kw in (("in_vitro_fitted", fitted), ("in_vivo_default", {})):
        controls = {}
        for ph in PHENOTYPES:  # declared order, not dict-iteration order
            r = fc.sim_batch(ph, "Control", 4000, 11, context="2d",
                             erastin_xc_inhib=0.0, transsulfuration_floor=0.0, **kw)
            controls[ph] = r["death_rate"]
        complete = fc.sim_batch("Glycolytic", "Control", 4000, 11, context="2d",
                                erastin_xc_inhib=1.0, transsulfuration_floor=0.0, **kw)
        noop = fc.sim_batch("Glycolytic", "Control", 4000, 11, context="2d",
                            erastin_xc_inhib=1.0, transsulfuration_floor=1.0, **kw)
        out[name] = {
            "untreated_controls": controls,
            "complete_block_death": complete["death_rate"],
            "complete_block_at_full_floor_death": noop["death_rate"],
            "untreated_glycolytic": controls["Glycolytic"],
            # The identity, measured rather than asserted.
            "full_floor_is_exact_noop":
                noop["death_rate"] == controls["Glycolytic"],
            "admissible": max(controls.values()) <= inp["viability_constraint"],
            "n_cells": 4000,
            "seed": 11,
        }
    return out


def assemble(raw: dict) -> dict:
    inp = raw["inputs"]
    km, s0 = inp["km_uM"], inp["reference_uM"]
    k, h = inp["k_erastin_uM"], inp["hill"]
    rows = []
    for d in DEPLETIONS:
        i = inhibition_from_depletion(d, s0, km)
        rows.append({"depletion": d, "residual_uM": s0 * (1 - d), "inhibition": i,
                     "erastin_equiv_uM": erastin_dose_for_inhibition(i, k, h)})
    doses = [r["erastin_equiv_uM"] for r in rows if r["erastin_equiv_uM"]]
    escape = [{"floor": f, "delivered_fraction": escape_delivery(f),
               "factor_at_full_block": cystine_factor(1.0, f)} for f in FLOORS]
    # Cross-check the ceiling against the committed calibration. Two things it
    # buys: the ceiling claim stops being pure topology, and it confirms this
    # page drives the engine the same way `calibrate_erastin.py` did -- a
    # harness that disagreed with the committed numbers would be measuring
    # something else and every figure above would be suspect.
    cross = None
    cascades = raw.get("cascades") or {}
    fitted = cascades.get("in_vitro_fitted")
    if fitted and inp.get("committed_top_dose_increment") is not None:
        mine = fitted["complete_block_death"] - fitted["untreated_glycolytic"]
        theirs = float(inp["committed_top_dose_increment"])
        top = inp["erastin_dose_grid_uM"][-1]
        cross = {
            "my_complete_block_increment": mine,
            "committed_top_dose_increment": theirs,
            "absolute_disagreement": abs(mine - theirs),
            "inhibition_at_committed_top_dose":
                top ** inp["hill"] / (top ** inp["hill"] + inp["k_erastin_uM"] ** inp["hill"]),
            "top_dose_uM": top,
        }

    out = {k2: raw[k2] for k2 in ("inputs", "cascades") if k2 in raw}
    out.update({
        "equivalence_rows": rows,
        "erastin_window_uM": [min(doses), max(doses)] if doses else None,
        "erastin_window_ratio": (max(doses) / min(doses)) if doses else None,
        "escape_rows": escape,
        "ceiling_cross_check": cross,
        "full_floor_factor_is_exactly_one":
            all(cystine_factor(i, 1.0) == 1.0
                for i in (0.0, 0.25, 0.5, 0.75, 0.9, 1.0)),
    })
    return out


def pct(x):
    return f"{100 * x:.1f}%"


def _ceiling_cross_check_prose(d: dict) -> str:
    """The ceiling, checked against the committed calibration rather than argued."""
    x = d.get("ceiling_cross_check")
    if not x:
        return ("The ceiling is not cross-checked here yet; run `--measure` with the "
                "compiled extension available.")
    return (
        f"**And it is not only an argument.** At the top of the committed dose grid "
        f"({x['top_dose_uM']:.0f} uM) the Hill map already sits at "
        f"{100 * x['inhibition_at_committed_top_dose']:.3f}% inhibition -- written to "
        f"three decimals because `pct` rounds it to a flat 100% and it is not 100% -- "
        f"so the committed "
        f"saturating-erastin increment IS a complete-block increment measured "
        f"independently of this page. It is {x['committed_top_dose_increment']:.4f}; a "
        f"complete cystine block here gives {x['my_complete_block_increment']:.4f}, "
        f"agreeing to {x['absolute_disagreement']:.4f} absolute. That does double duty: "
        f"the ceiling is confirmed, and this page is confirmed to drive the engine the "
        f"same way `calibrate_erastin.py` did -- a harness disagreeing with the committed "
        f"numbers would be measuring something else and every figure above would be "
        f"suspect.")


def render(d: dict) -> str:
    inp = d["inputs"]
    c = d.get("cascades") or {}
    L = [
        "# What a cyst(e)inase is worth in erastin, and what the engine will not say",
        "",
        "Generated by `scripts/cysteine_depletion_dose_equivalence.py`. No layer is "
        "added here; every parameter used is one the engine already carries.",
        "",
        "## Dose equivalence, which needs no simulation",
        "",
        f"The #909 map turns a residual cystine concentration into a supply fraction "
        f"(Km {inp['km_uM']} uM, reference {inp['reference_uM']:.1f} uM). The committed "
        f"erastin fit turns a dose into the same quantity "
        f"(`dose^{inp['hill']:.0f}/(dose^{inp['hill']:.0f} + "
        f"{inp['k_erastin_uM']:.0f}^{inp['hill']:.0f})`, so "
        f"{inp['k_erastin_uM']:.0f} uM is the half-maximal dose, "
        f"RMSE {inp['erastin_fit_rmse']}, {inp['erastin_fit_phenotype']}). Composing "
        "them is arithmetic over published constants.",
        "",
        "| cystine removed | residual (uM) | supply cut | erastin equivalent |",
        "| --: | --: | --: | --: |",
    ]
    for r in d["equivalence_rows"]:
        dose = f"{r['erastin_equiv_uM']:.2f} uM" if r["erastin_equiv_uM"] else "--"
        L.append(f"| {pct(r['depletion'])} | {r['residual_uM']:.1f} | "
                 f"{pct(r['inhibition'])} | {dose} |")
    lo, hi = d["erastin_window_uM"]
    L += [
        "",
        f"The enzyme's whole operating range lands in **{lo:.2f}-{hi:.2f} uM** erastin, a "
        f"**{d['erastin_window_ratio']:.1f}x** window, every point inside the calibrated "
        f"dose grid ({inp['erastin_dose_grid_uM'][0]}-{inp['erastin_dose_grid_uM'][-1]} uM). "
        f"That narrowness is the Hill coefficient of {inp['hill']:.0f}, not a property of "
        "the enzyme: a steep switch compresses a wide range of supply cuts into a short "
        "stretch of dose.",
        "",
        "## A structural ceiling",
        "",
        "The enzyme cannot cut supply below zero, so **through this mechanism it cannot "
        "exceed a saturating erastin dose**. That is the engine's topology rather than a "
        "measurement -- both routes enter through the same supply factor -- and it is the "
        "useful kind of negative. The advantage PMID 39696900 claims, \"efficacy, "
        "persistence, and precision\", has to come from targeting and persistence, because "
        "mechanistic depth is already spent by a drug that exists.",
        "",
        _ceiling_cross_check_prose(d),
        "",
        "## The transsulfuration escape is an exact no-op switch",
        "",
        "`xc_cystine_factor()` is `1 - erastin_xc_inhib * (1 - transsulfuration_floor)`, so "
        "the supply actually removed is `inhib * (1 - floor)` and the fraction delivered is "
        "`1 - floor` -- independent of the inhibition entirely.",
        "",
        "| transsulfuration floor | nominal cut delivered | factor at a COMPLETE block |",
        "| --: | --: | --: |",
    ]
    for r in d["escape_rows"]:
        L.append(f"| {pct(r['floor'])} | {pct(r['delivered_fraction'])} "
                 f"| {r['factor_at_full_block']:.3f} |")
    L += [
        "",
        f"At a floor of 1 the factor is **exactly 1 for every inhibition** "
        f"(`{d['full_floor_factor_is_exactly_one']}`): a fully competent transsulfuration "
        "pathway makes the entire cystine arm a precise no-op, at any depletion depth and "
        "at any cascade. This is the quantitative reason the second bacterial paper "
        "(PMID 41704298) had to generate H2S to inhibit the bypass, and it needs no "
        "calibration to state.",
        "",
        "## What this refuses to say",
        "",
        "**No absolute kill rate.** Not \"not yet measured\" -- refused, because both "
        "available cascades fail in opposite directions and the repo already has the gate "
        f"that catches it (baseline viability {pct(inp['viability_constraint'])}, from "
        "`headline-at-fitted-cascade.json`).",
        "",
    ]
    if not c:
        L += ["The cascade measurements are not committed yet; run `--measure` with the "
              "compiled extension available.", ""]
    else:
        L += ["| cascade | untreated Glycolytic | worst untreated phenotype | admissible "
              "| complete block |",
              "| --- | --: | --: | --: | --: |"]
        unknown = sorted(set(c) - set(CASCADE_ORDER))
        if unknown:
            raise MissingInputError(
                f"cascades not in CASCADE_ORDER, so their row order would be "
                f"serialisation-dependent: {unknown}")
        for name in [n for n in CASCADE_ORDER if n in c]:
            m = c[name]
            worst = max(m["untreated_controls"].values())
            L.append(f"| `{name}` | {pct(m['untreated_glycolytic'])} | {pct(worst)} "
                     f"| {m['admissible']} | {pct(m['complete_block_death'])} |")
        fitted = c.get("in_vitro_fitted", {})
        default = c.get("in_vivo_default", {})
        L += [
            "",
            f"The in-vitro fitted cascade reproduces the real CTRPv2 erastin curve and the "
            f"arm reaches {pct(fitted.get('complete_block_death', 0))} at a complete block "
            f"-- but its untreated Glycolytic control dies at "
            f"{pct(fitted.get('untreated_glycolytic', 0))}, over the "
            f"{pct(inp['viability_constraint'])} constraint, and its worst untreated "
            f"phenotype at {pct(max(fitted.get('untreated_controls', {0: 0}).values()))}. "
            f"The admissible in-vivo default has every control at or near zero and a "
            f"complete block kills "
            f"{pct(default.get('complete_block_death', 0))} -- the arm does nothing.",
            "",
            "So the magnitude of this arm is not predictable by this engine today, for a "
            "reason already on record rather than a new one: the same in-vitro/in-vivo "
            "cascade disjunction `headline_at_fitted.py` published for every other "
            "headline. A FIRST DRAFT OF THIS PAGE REPORTED THE FITTED-CASCADE KILL TABLE "
            "AS A RESULT. Running the untreated controls is what stopped it, and that "
            "check is the only reason this section exists.",
            "",
            "The no-op identity is also confirmed by measurement rather than left as "
            "algebra: at a full floor a complete block reproduces the untreated rate "
            "EXACTLY at both cascades "
            f"({fitted.get('full_floor_is_exact_noop')}, "
            f"{default.get('full_floor_is_exact_noop')}).",
            "",
        ]
    L += [
        "## What would unblock the magnitude",
        "",
        "A cascade that is admissible AND reproduces an inducer dose-response. That is the "
        "standing `#500` problem, not a cysteine-specific one. The cysteine-specific piece "
        "is cheaper and named in `cysteine-depletion-feasibility.md`: the transporter Km "
        "measured in a pancreatic line, since the only kinetic analysis in the cancer "
        "literature is a single lung-cell-line paper.",
        "",
    ]
    return "\n".join(L)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--measure", action="store_true",
                    help="run both cascades (needs the compiled extension) and commit them")
    ap.add_argument("--render-only", action="store_true",
                    help="rebuild the page from the committed JSON")
    a = ap.parse_args()
    try:
        if a.render_only:
            if not OUT_JSON.exists():
                print(f"missing {OUT_JSON}", file=sys.stderr)
                return 1
            data = assemble(json.loads(OUT_JSON.read_text()))
        else:
            raw = {"inputs": load_inputs()}
            if OUT_JSON.exists():
                prior = json.loads(OUT_JSON.read_text())
                if "cascades" in prior:
                    raw["cascades"] = prior["cascades"]
            if a.measure:
                raw["cascades"] = measure_cascades(raw["inputs"])
            data = assemble(raw)
            OUT_JSON.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
    except MissingInputError as e:
        print(f"refusing to render: {e}", file=sys.stderr)
        return 1
    OUT_MD.write_text(render(data))
    print(f"wrote {OUT_MD}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
