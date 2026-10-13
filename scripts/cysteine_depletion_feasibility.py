#!/usr/bin/env python3
"""Can the engine express bacterial cysteine depletion, and what would anchor it?

WHY
---
`corpus_retrieval_search.py` turned up PMID 39696900 (Adv. Mater. 2025): engineered
bacteria that colonise a pancreatic tumour and, under the stimulus of tumour
hypoxia, secrete cyst(e)inase (CGL) to destroy the cysteine the cell needs for
glutathione, collapsing the ferroptosis defence. The engine already carries
`erastin_xc_inhib` and `transsulfuration_floor` (#502), so the question the
layer-freeze policy forces is not "add a bacterial layer" but: is this arm
already expressible, and does a calibration target exist for the one piece that
is genuinely missing?

THE FIRST THING TO CORRECT IS MY OWN FRAMING
--------------------------------------------
The search that found the paper reported a pool of exactly ONE record across
7,727,217 -- `ferroptosis AND "engineered bacteria"` -- and I read that as an
unexplored direction. It is not. That number measures a PHRASE PAIR. The
mechanism it stands for, cysteine/cystine supply against ferroptosis, returns
3,404 records; `cyst(e)inase` alone returns 20, and they include Cramer 2017
(Nat. Med., PMID 27869804), Badgley 2020 (Science, PMID 32241947, the same
mechanism in the same tumour type by both Slc7a11 deletion AND cyst(e)inase),
and a SECOND bacterial paper, PMID 41704298. So the pool of one was a fact about
vocabulary, and the inference drawn from it was an overclaim. Both numbers are
reported here, from the same index, so neither can be quoted alone.

WHAT IS ALREADY EXPRESSIBLE, AND THE UNIT PROBLEM THAT IS NOT
-------------------------------------------------------------
`xc_cystine_factor() = 1 - erastin_xc_inhib * (1 - transsulfuration_floor)` is a
SUPPLY multiplier on glutathione resynthesis, not a transporter model. Erastin
shuts the door; cyst(e)inase destroys what is outside it. Badgley 2020 reaches
ferroptosis in PDAC by BOTH routes, consistent with the collapse onto one supply
factor the engine already makes -- though two routes to one endpoint do not
establish that they are quantitatively interchangeable, and nothing here tests
that. What IS certain is a unit mismatch: `erastin_xc_inhib` is a FRACTION OF
SUPPLY and an enzyme delivers a CONCENTRATION. The map between them is the
transporter's own saturation curve.

THE CALIBRATION TARGET, AND IT IS A FAILABLE ONE
------------------------------------------------
Cystine uptake is saturable and someone has measured it: Km = 179.4 +/- 26.7 uM
in A549 (PMID 34287957). The in-vitro reference concentration is not a
measurement at all but a formulation constant -- RPMI-1640 specifies L-cystine
2HCl at 65.0 mg/L, which is 207.7 uM. Those two numbers fix the map with nothing
left to choose, and they make it FAIL if the shape is wrong.

THE HEADLINE IS THAT DEPLETION BUYS LESS THAN IT LOOKS LIKE, AND THE GAP IS BIG
-------------------------------------------------------------------------------
Culture medium sits just past half-saturation (207.7 against a Km of 179.4), so
the first substrate removed barely moves the flux: the marginal return on
depletion starts at Km/(Km+S0), well under one. It does not merely approach one --
it CROSSES it and ends at (Km+S0)/Km, so the last percent of substrate removed is
worth several times the first. Taking HALF the cystine out of the medium cuts
supply by about a third, and a 90% supply cut needs about 95% depletion. The
derivative RISES with depth, the opposite of a saturating drug -- so what matters
is holding the concentration near zero, not the average depletion. That is a mechanistic argument for a continuously-producing
living source over a bolus enzyme, and it is the quantitative content of the
paper's own complaint that existing methods "fail to efficiently block cysteine
metabolism".

THE PIECE THAT IS GENUINELY MISSING IS WIRING, NOT DATA
-------------------------------------------------------
PMID 41422256 reports that hypoxia activates the transsulfuration pathway and so
confers ferroptosis RESISTANCE under cystine limitation. The bacteria are
hypoxia-TRIGGERED. So the trigger co-locates with the escape, and the second
bacterial paper (41704298) independently added H2S to block transsulfuration --
which is the modification the engine's own parameter structure predicts you would
need. The engine CANNOT state that today, and the reason is checked against the
source rather than asserted: `transsulfuration_floor` is a `Params` field with no
per-cell counterpart, while the hypoxia field is per-cell. Raising `cell.nrf2`
does not substitute, because resynthesis is MULTIPLIED by the cystine factor, so
a complete block is zero at every NRF2 level. A hypoxia-graded escape needs the
floor moved to `Cell` or a per-cell params path; that is a wiring change with a
byte-identity obligation, not a missing dataset.

NOTHING IS BUILT HERE. Per CONTRIBUTING.md's layer policy this names the target
and the gap; it adds no axis.

OFFLINE CONTRACT
----------------
Everything this script reads is committed: the published constants, the measured
pool counts, the landmark pool membership, and the engine source. It needs no
corpus, no network and no index, so CI runs it in seconds.

The pool counts and the landmark membership WERE measured against the local
21 GB `corpus/atlas/retrieval.sqlite`, which is gitignored, so they cannot be
re-derived from a fresh clone. The queries are committed beside the counts for
exactly that reason -- re-measuring is
`SELECT COUNT(*) FROM doc WHERE doc MATCH ?` per row of
`cysteine_depletion_pools.csv`, and for the landmarks the same queries with
`SELECT pmid` intersected against the named papers. A count that moves because
the corpus grew is expected; the `measured_against_index_records` column records
which build each was taken on so a stale figure is visible rather than silent.

Usage:
    python scripts/cysteine_depletion_feasibility.py
    python scripts/cysteine_depletion_feasibility.py --render-only
"""

import argparse
import csv
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from config import PROJECT_ROOT  # noqa: E402

CONST_SRC = PROJECT_ROOT / "analysis" / "calibration" / "cysteine_depletion_constants.csv"
POOL_SRC = PROJECT_ROOT / "analysis" / "calibration" / "cysteine_depletion_pools.csv"
LANDMARK_SRC = PROJECT_ROOT / "analysis" / "calibration" / "cysteine_depletion_landmarks.csv"
PARAMS_RS = PROJECT_ROOT / "simulations" / "ferroptosis-core" / "src" / "params.rs"
BIOCHEM_RS = PROJECT_ROOT / "simulations" / "ferroptosis-core" / "src" / "biochem.rs"
CELL_RS = PROJECT_ROOT / "simulations" / "ferroptosis-core" / "src" / "cell.rs"
# The one spatial binary the #502 arm can be driven in (#918) and the only
# consumer of the per-cell floor (#920). Scanned, not assumed: the claim this
# page used to publish -- that `erastin_xc_inhib` is reachable in no simulation
# binary -- was TRUE when written and false the moment #918 landed, and nothing
# failed, because the sentence was prose rather than a measurement.
TME3D_RS = PROJECT_ROOT / "simulations" / "sim-tme-3d" / "src" / "main.rs"
OUT_MD = PROJECT_ROOT / "analysis" / "calibration" / "cysteine-depletion-feasibility.md"
OUT_JSON = PROJECT_ROOT / "analysis" / "calibration" / "cysteine-depletion-feasibility.json"

# Depletion fractions the table reports. Chosen to straddle the knee, not to
# flatter the curve: 0.5 is the "half the substrate gone" intuition the gap
# refutes, and 0.99 is where a living source is doing something a bolus cannot.
DEPLETIONS = [0.25, 0.5, 0.75, 0.9, 0.95, 0.99]
# Supply cuts the inverse table reports, i.e. the `erastin_xc_inhib` values a
# consumer would actually set.
SUPPLY_CUTS = [0.25, 0.5, 0.75, 0.9]
# The cut the Km-sensitivity table and the regime comparison are both quoted at.
# Named rather than repeated, so the column header interpolates it instead of
# spelling "90%" beside a derived value.
SENSITIVITY_CUT = 0.9
# The depletion fraction the regime contrast is quoted at.
REGIME_DEPLETION = 0.5


class MissingInputError(RuntimeError):
    """Raised rather than rendering a page from absent or partial inputs."""


def load_constants(path: Path) -> dict:
    if not path.exists():
        raise MissingInputError(f"missing committed constants: {path}")
    rows = {r["constant"]: r for r in csv.DictReader(path.open())}
    need = ["xct_cystine_km", "rpmi1640_cystine_2hcl_mass", "cystine_2hcl_molar_mass"]
    absent = [k for k in need if k not in rows]
    if absent:
        raise MissingInputError(f"{path} lacks required constants: {', '.join(absent)}")
    return rows


def load_pools(path: Path) -> dict:
    if not path.exists():
        raise MissingInputError(f"missing committed pool counts: {path}")
    # dict preserves insertion order, so this is the CSV's own row order.
    rows = {r["key"]: r for r in csv.DictReader(path.open())}
    for k in ["phrase_pair_ferro_engbact", "cystine_import_ferro"]:
        if k not in rows:
            raise MissingInputError(f"{path} lacks the {k} row, which the correction needs")
    return rows


def load_landmarks(path: Path) -> list:
    """Which pool each named paper was actually found in, MEASURED.

    Written as data because the first draft of this page asserted in prose that
    PMID 41704298 was among the 20 `cyst(e)inase` records. It is not -- it sits
    in the cysteine-depletion, bacteria and transsulfuration pools. A sentence
    naming a paper's pool is a claim about a measurement, so it is one here.
    """
    if not path.exists():
        raise MissingInputError(f"missing committed landmark membership: {path}")
    rows = [r for r in csv.DictReader(path.open())]
    if not rows:
        raise MissingInputError(f"{path} holds no landmark rows")
    for r in rows:
        if not r["found_in_pools"].strip():
            raise MissingInputError(
                f"landmark {r['pmid']} matched no pool; the prose cannot cite it")
    return [{"pmid": r["pmid"],
             "pools": r["found_in_pools"].split("|"),
             "note": r["note"]} for r in rows]


def saturation(s: float, km: float) -> float:
    """Michaelis-Menten fractional occupancy `s/(km+s)`. Zero at zero substrate."""
    if s <= 0:
        return 0.0
    return s / (km + s)


def relative_supply(s: float, s0: float, km: float) -> float:
    """Uptake at `s` as a fraction of uptake at the reference `s0`."""
    ref = saturation(s0, km)
    if ref <= 0:
        raise MissingInputError("reference concentration gives zero uptake")
    return saturation(s, km) / ref


def inhibition_from_depletion(d: float, s0: float, km: float) -> float:
    """The `erastin_xc_inhib` a depletion fraction `d` is worth."""
    return 1.0 - relative_supply(s0 * (1.0 - d), s0, km)


def depletion_for_inhibition(target: float, s0: float, km: float) -> float:
    """Invert the map: the depletion fraction a target supply cut requires.

    Closed form, not a search. `rel = sat(s)/sat(s0)` and `sat` is invertible:
    `sat(s) = y` gives `s = km*y/(1-y)`. A search here would make the test that
    checks the round trip a test of the search's tolerance.
    """
    rel = 1.0 - target
    y = rel * saturation(s0, km)
    if y >= 1.0:
        raise MissingInputError("target implies occupancy at or above one")
    s = km * y / (1.0 - y)
    return 1.0 - s / s0


def marginal_return(d: float, s0: float, km: float) -> float:
    """d(inhibition)/d(depletion), the extra supply cut one more unit buys.

    Analytic: with `r = s0/km`, inhibition is `1 - (1-d)(1+r)/(1+r(1-d))`, whose
    derivative is `(1+r)/(1+r-r*d)^2`. At `d=0` this is `km/(km+s0)` -- under one
    whenever there is any substrate at all -- and at `d=1` it is `1+r`, so it
    crosses one rather than approaching it. An earlier draft of this file said it
    "only reaches 1 asymptotically", which the table it rendered contradicted on
    the same page.
    """
    r = s0 / km
    return (1.0 + r) / (1.0 + r - r * d) ** 2


def marginal_crossover(s0: float, km: float) -> float:
    """The depletion at which one more unit removed starts buying more than its share.

    Solves `(1+r)/(1+r-r*d)^2 = 1` for `d`, i.e. `d = (1+r-sqrt(1+r))/r`. Below it
    depletion under-delivers, above it over-delivers; it exists for any `r > 0`.
    """
    r = s0 / km
    return (1.0 + r - (1.0 + r) ** 0.5) / r


def read_structural_facts() -> dict:
    """Re-derive the engine-structure claims from the Rust source.

    Written this way because a guard that reads only this script's own JSON
    cannot fail when the engine changes underneath it -- the defect
    `tests/test_generator_guards_pin_source.py` exists to block. Every claim
    below is a grep over committed source at generation time.
    """
    for p in (PARAMS_RS, BIOCHEM_RS, CELL_RS):
        if not p.exists():
            raise MissingInputError(f"missing engine source: {p}")
    params_src = PARAMS_RS.read_text()
    biochem_src = BIOCHEM_RS.read_text()
    cell_src = CELL_RS.read_text()
    # Production code only -- everything before the crate's test module.
    tme3d_src = TME3D_RS.read_text()
    cut = tme3d_src.find("\n#[cfg(test)]")
    tme3d_prod = tme3d_src if cut < 0 else tme3d_src[:cut]

    # A `Params` field declaration, not a mention: `pub <name>: f64,`.
    def declares(src: str, name: str) -> bool:
        return re.search(rf"^\s*pub\s+{re.escape(name)}\s*:", src, re.M) is not None

    # The resynthesis STATEMENTS, not lines. A line-based scan broke the moment
    # `cargo fmt` reflowed these across four lines: `gsh +=` landed on the first
    # and `xc_cystine_factor` on the fourth, so the gated count read 0 and this
    # page was about to publish that NRF2 CAN rescue a complete cystine block --
    # a false structural claim produced by formatting, not by a code change.
    # Statements end at the first `;`, which is what the claim is actually about.
    resynth = [
        m.group(0) for m in re.finditer(r"\bgsh\s*\+=.*?;", biochem_src, re.S)
    ]
    gated = [ln for ln in resynth if "xc_cystine_factor" in ln]
    nrf2_scaled = [ln for ln in resynth if "nrf2" in ln]
    return {
        "floor_is_params_field": declares(params_src, "transsulfuration_floor"),
        "floor_is_cell_field": declares(cell_src, "transsulfuration_floor"),
        "inhib_is_params_field": declares(params_src, "erastin_xc_inhib"),
        "nrf2_is_cell_field": declares(cell_src, "nrf2"),
        "gsh_resynthesis_sites": len(resynth),
        "resynthesis_sites_gated_by_cystine": len(gated),
        "resynthesis_sites_scaled_by_nrf2": len(nrf2_scaled),
        # The load-bearing one: if EVERY resynthesis site is multiplied by the
        # cystine factor, a complete block is zero regardless of NRF2.
        "nrf2_cannot_rescue_complete_block": bool(resynth) and len(gated) == len(resynth),
        "escape_can_vary_per_cell": declares(cell_src, "transsulfuration_floor"),
        # Is the #502 arm driveable in a spatial binary at all, and does anything
        # write the per-cell floor? Both measured on PRODUCTION code only: the
        # `#[cfg(test)]` module is cut first, because an A/B test setting a
        # parameter is not a consumer, and this repo has already shipped a
        # detector whose only match was inside a test block.
        "arm_is_spatially_reachable": re.search(
            r"apply_param_overrides\(\s*&mut\s+params", tme3d_prod
        )
        is not None,
        "floor_has_spatial_consumer": re.search(
            r"\.transsulfuration_floor\s*=\s*Some\(", tme3d_prod
        )
        is not None,
    }


def read_inputs() -> dict:
    """The RAW inputs, with nothing derived: committed constants + pool counts.

    Kept separate from `assemble` so that `assemble` can be handed the stored
    artifact and re-derive every reported quantity from it. That is what makes
    `--render-only` recompute rather than echo, and it is what
    `tests/test_artifact_freshness.py` enforces: a generator whose render-only
    branch re-renders stored derived fields makes every guard reading one of
    those fields inert.
    """
    consts = load_constants(CONST_SRC)
    pools = load_pools(POOL_SRC)
    return {
        "landmarks": load_landmarks(LANDMARK_SRC),
        # Whole rows, `sd` included. An earlier version kept only a subset here
        # and derived `km_sd` from the dropped column, so a second pass through
        # `assemble` could not reproduce the first -- the idempotency the
        # render-only contract rests on.
        "constants": {k: dict(v) for k, v in consts.items()},
        # ORDERED PAIRS, not a dict: `json.dumps(sort_keys=True)` reorders a
        # dict, so `--render-only` would render the table alphabetically while a
        # fresh run rendered it in query order, producing two different documents
        # from one artifact. That exact defect is on record in this repo from
        # `atlas_taxonomy_reach.py`; the render-only guard caught it again here.
        "pool_rows": [{"key": k, "query": v["query"], "pool": int(v["pool"])}
                      for k, v in pools.items()],
        "index_records": int(
            pools["phrase_pair_ferro_engbact"]["measured_against_index_records"]),
    }


RAW_KEYS = ("constants", "pool_rows", "index_records", "landmarks")


def assemble(raw: dict) -> dict:
    """Derive every reported quantity from the raw inputs. Idempotent.

    Handed its own output, it recomputes the same derived fields from the same
    raw fields and returns an equal dict -- so `--render-only` can be fed the
    committed JSON and still do the arithmetic.
    """
    absent = [k for k in RAW_KEYS if k not in raw]
    if absent:
        raise MissingInputError(f"stored artifact lacks raw inputs: {', '.join(absent)}")
    consts = raw["constants"]
    km = float(consts["xct_cystine_km"]["value"])
    km_sd = float(consts["xct_cystine_km"]["sd"])
    mass = float(consts["rpmi1640_cystine_2hcl_mass"]["value"])
    mw = float(consts["cystine_2hcl_molar_mass"]["value"])
    # mg/L / (g/mol) = mmol/L; x1000 -> uM.
    s0 = mass / mw * 1000.0

    rows = [{
        "depletion": d,
        "residual_uM": s0 * (1.0 - d),
        "inhibition": inhibition_from_depletion(d, s0, km),
        "marginal_return": marginal_return(d, s0, km),
    } for d in DEPLETIONS]
    inverse = []
    for t in SUPPLY_CUTS:
        need = depletion_for_inhibition(t, s0, km)
        inverse.append({"supply_cut": t, "depletion_required": need,
                        "residual_uM": s0 * (1.0 - need)})

    # Km sensitivity. The Km is ONE cell line in ONE paper and the target disease
    # is pancreatic, so the scale is the weak part of the anchor and the report
    # must show how much rides on it.
    sens = [{
        "label": label,
        "km_uM": k,
        "saturation_at_reference": saturation(s0, k),
        "initial_marginal_return": marginal_return(0.0, s0, k),
        "depletion_for_cut": depletion_for_inhibition(SENSITIVITY_CUT, s0, k),
    } for label, k in [("km_low", km - km_sd), ("km_point", km), ("km_high", km + km_sd),
                       ("km_half", km / 2), ("km_double", km * 2)]]

    # THE HEADLINE IS REGIME-SPECIFIC, and saying so is the point of this block.
    # The gap between depletion and supply cut exists because culture medium sits
    # ABOVE the Km. Well BELOW the Km the saturation curve is almost linear and
    # depletion maps nearly one-to-one, so the finding is about cell culture and
    # does not transfer unexamined to a tumour. The low reference used here is the
    # committed free-CYSTEINE figure, which is NOT free cystine, so this is an
    # ILLUSTRATION of the regime effect and not a claim about the in-vivo map.
    low_ref = float(consts["plasma_free_cysteine_baseline"]["value"])
    # The two comparison points are STORED so the prose interpolates them instead
    # of naming "50%" and "90%" beside derived values -- the shape
    # `tests/test_no_handwritten_figures_in_generators.py` exists to stop.
    d_cmp, cut_cmp = REGIME_DEPLETION, SENSITIVITY_CUT
    regime = {
        "depletion_compared": d_cmp,
        "supply_cut_compared": cut_cmp,
        "depletion_for_cut_reference": depletion_for_inhibition(cut_cmp, s0, km),
        "low_reference_uM": low_ref,
        "low_reference_is_free_cystine": False,
        "saturation_at_low_reference": saturation(low_ref, km),
        "inhibition_at_half_depletion_low": inhibition_from_depletion(d_cmp, low_ref, km),
        "inhibition_at_half_depletion_reference": inhibition_from_depletion(d_cmp, s0, km),
        "depletion_for_cut_low": depletion_for_inhibition(cut_cmp, low_ref, km),
    }

    out = {k: raw[k] for k in RAW_KEYS}
    out.update({
        "reference_uM": s0,
        "km_uM": km,
        "km_sd_uM": km_sd,
        "saturation_at_reference": saturation(s0, km),
        "initial_marginal_return": marginal_return(0.0, s0, km),
        "final_marginal_return": marginal_return(1.0, s0, km),
        "marginal_crossover_depletion": marginal_crossover(s0, km),
        "depletion_rows": rows,
        "inverse_rows": inverse,
        "km_sensitivity": sens,
        "regime": regime,
        "structure": read_structural_facts(),
    })
    return out


def pct(x: float) -> str:
    return f"{100 * x:.1f}%"


def render(d: dict) -> str:
    st = d["structure"]
    pool = {r["key"]: r["pool"] for r in d["pool_rows"]}
    query = {r["key"]: r["query"] for r in d["pool_rows"]}
    # Every figure and PMID in the prose comes from the stored constants, not
    # from a literal typed beside them. "Derived number, handwritten prose" is a
    # recorded defect class in this repo: the figure updates and the sentence
    # next to it does not.
    c = d["constants"]
    g = d["regime"]
    half = next(r for r in d["depletion_rows"] if r["depletion"] == 0.5)
    deep = next(r for r in d["inverse_rows"] if r["supply_cut"] == 0.9)
    L = [
        "# Bacterial cysteine depletion: is it expressible, and what would anchor it?",
        "",
        "Generated by `scripts/cysteine_depletion_feasibility.py`. No layer is added here.",
        "",
        "## The pool of one measured vocabulary, not the field",
        "",
        f"`{query['phrase_pair_ferro_engbact']}` returns "
        f"**{pool['phrase_pair_ferro_engbact']}** record across "
        f"{d['index_records']:,}, which is how PMID 39696900 was found and which I read as an "
        "unexplored direction. That reading was wrong. The same index, asked about the "
        f"mechanism rather than the phrase, returns **{pool['cystine_import_ferro']:,}** "
        f"for `{query['cystine_import_ferro']}`. Both are reported here so neither "
        "can be quoted alone, and deliberately NOT as a ratio: a multiplier over a pool of one "
        "would invite that one to be read as a baseline.",
        "",
        "| pool | query | records |",
        "| --- | --- | --: |",
    ]
    for r in d["pool_rows"]:
        L.append(f"| `{r['key']}` | `{r['query']}` | {r['pool']:,} |")
    L += [
        "",
        "The papers this page leans on, with the pools each was actually found in -- measured, "
        "because an earlier draft asserted in prose that PMID 41704298 was among the "
        f"{pool['cysteinase_any']} `cyst(e)inase` records and it is not:",
        "",
        "| PMID | found in | why it matters |",
        "| --- | --- | --- |",
    ]
    for lm in d["landmarks"]:
        L.append(f"| {lm['pmid']} | {', '.join(f'`{x}`' for x in lm['pools'])} "
                 f"| {lm['note']} |")
    two_bact = [lm["pmid"] for lm in d["landmarks"]
                if "bacteria_cysteine_ferro" in lm["pools"]]
    L += [
        "",
        f"The `bacteria_cysteine_ferro` pool of {pool['bacteria_cysteine_ferro']} holds "
        f"{len(two_bact)} bacterial papers on this mechanism ({', '.join(two_bact)}), not one. "
        "This is a literature, not a gap.",
        "",
        "## The target: a measured Km and a formulation constant",
        "",
        f"- cystine-uptake Km **{d['km_uM']} +/- {d['km_sd_uM']} uM** "
        f"(PMID {c['xct_cystine_km']['source_pmid']}, A549)",
        f"- RPMI-1640 reference **{d['reference_uM']:.1f} uM** "
        f"({c['rpmi1640_cystine_2hcl_mass']['value']} mg/L L-cystine 2HCl / "
        f"{c['cystine_2hcl_molar_mass']['value']} g/mol), a formulation constant, "
        "not a measurement",
        "",
        f"Medium therefore sits at **{pct(d['saturation_at_reference'])}** of saturation -- just "
        "past the knee. Nothing in the map is chosen, so it can fail.",
        "",
        f"**One cell line is not a shortcut, it is the literature.** Of "
        f"{pool['cystine_uptake_any']} records naming cystine uptake, only "
        f"{pool['cystine_uptake_kinetic']} report kinetic constants, and one of those two "
        "studies cultured fibroblasts in cystinuria and Hartnup disease (PMID 5054467, dated "
        "in the landmark table rather than here) -- not a tumour. So in cancer cell lines "
        "there is exactly ONE kinetic analysis across "
        f"{d['index_records']:,} records, which is what PMID "
        f"{c['xct_cystine_km']['source_pmid']} says of itself: very few reports exist. "
        f"{pool['cystine_uptake_pancreatic']} records pair cystine uptake with pancreas and "
        "NONE of them is kinetic (they measure expression, drug response or isotope tracing). "
        "That makes measuring this Km in a PDAC line a genuinely novel experiment rather than "
        "a tidying-up exercise, and it is the single cheapest thing that would tighten this "
        "page.",
        "",
        "What the map assumes, stated because these are where it would break before the "
        "arithmetic does (no count in this sentence: a spelled one goes stale the moment "
        "another assumption is added):",
        "",
        "0. **The map sets `erastin_xc_inhib` alone.** `xc_cystine_factor()` is "
        "`1 - erastin_xc_inhib * (1 - transsulfuration_floor)`, so the figures below equal "
        "the engine's supply factor exactly only at a zero floor; a non-zero floor adds "
        "supply back on top, which is correct because the floor is a SEPARATE source "
        "(methionine, not import) and not part of the transporter curve. Set the two "
        "independently and do not count the floor twice.",
        "1. **Glutathione resynthesis is proportional to cystine uptake.** This is the "
        "engine's own assumption, not one added here -- #502 made the NRF2-driven "
        "resynthesis term *be* the cystine supply -- but it is the step that lets a "
        "transporter rate stand in for a supply fraction, so it is the first thing to "
        "doubt if the map disagrees with an assay.",
        "2. **The transporter is the rate limit.** Macropinocytosis is a parallel route to "
        "cysteine that an enzyme in the medium does not touch (PMID 35287706, where it "
        "mitigates sorafenib-induced ferroptosis), and the map is blind to it. A cell line "
        "using it would be depleted on paper and supplied in fact.",
        "3. **An initial-rate constant is being used as a steady-state one.** The Km comes "
        f"from a {c['xct_assay_linear_window']['value']}-"
        f"{c['xct_assay_linear_window']['unit']} linear-uptake window; the engine applies "
        "the factor every step of a long run. Nothing here measures whether the constant "
        "holds once the intracellular pool has drawn down.",
        "",
        "## Depletion buys less supply reduction than its percentage suggests",
        "",
        "| cystine removed | residual (uM) | supply cut it is worth | marginal return |",
        "| --: | --: | --: | --: |",
    ]
    for r in d["depletion_rows"]:
        L.append(f"| {pct(r['depletion'])} | {r['residual_uM']:.1f} | {pct(r['inhibition'])} "
                 f"| {r['marginal_return']:.2f} |")
    L += [
        "",
        f"Removing **{pct(half['depletion'])}** of the cystine is worth an "
        f"`erastin_xc_inhib` of only **{pct(half['inhibition'])}**. The marginal return starts "
        f"at **{d['initial_marginal_return']:.2f}** -- `Km/(Km+S0)`, necessarily under one "
        f"wherever there is substrate -- crosses one at "
        f"**{pct(d['marginal_crossover_depletion'])}** depletion and ends at "
        f"**{d['final_marginal_return']:.2f}**, so the last substrate removed is worth several "
        "times the first. That is the opposite of a saturating drug. The quantity that matters "
        "is therefore how close to zero the concentration is HELD, not the average depletion, "
        "which is a mechanistic argument for a continuously producing living source over a "
        "bolus enzyme.",
        "",
        "| supply cut wanted | depletion required | residual (uM) |",
        "| --: | --: | --: |",
    ]
    for r in d["inverse_rows"]:
        L.append(f"| {pct(r['supply_cut'])} | {pct(r['depletion_required'])} "
                 f"| {r['residual_uM']:.1f} |")
    L += [
        "",
        f"A {pct(deep['supply_cut'])} supply cut needs {pct(deep['depletion_required'])} "
        f"depletion, down to {deep['residual_uM']:.1f} uM. That is the quantitative content of "
        "PMID 39696900's own complaint that existing methods fail to block cysteine "
        "metabolism efficiently.",
        "",
        "### The gap is a CULTURE phenomenon, and that qualifies the headline",
        "",
        f"Everything above is measured against a medium reference ABOVE the Km "
        f"({d['reference_uM']:.1f} against {d['km_uM']} uM), which is why the curve is flat "
        "where the first substrate is removed. Well BELOW the Km the same curve is almost "
        "linear and depletion maps nearly one-to-one. Taking the committed free-cysteine "
        f"figure of {g['low_reference_uM']} uM as a low reference, "
        f"{pct(g['depletion_compared'])} depletion is worth "
        f"**{pct(g['inhibition_at_half_depletion_low'])}** instead of "
        f"**{pct(g['inhibition_at_half_depletion_reference'])}**, and a "
        f"{pct(g['supply_cut_compared'])} supply cut needs "
        f"**{pct(g['depletion_for_cut_low'])}** depletion instead of "
        f"**{pct(g['depletion_for_cut_reference'])}**.",
        "",
        "That figure is an ILLUSTRATION of the regime effect, not a claim about the in-vivo "
        "map: it is free CYSTEINE, not free cystine, and the transporter imports cystine. The "
        "point that survives is structural -- the in-vitro regime the engine's erastin fit "
        "lives in and the in-vivo regime an enzyme therapy acts in sit on OPPOSITE sides of "
        "the transporter's half-saturation point, so a fraction fitted in culture does not "
        "carry across unexamined. Establishing the free-cystine concentration in tumour "
        "interstitium is what would close this, and nothing here does.",
        "",
        "## How much rides on the Km",
        "",
        "One cell line (A549, lung) in one paper, against a pancreatic target, so the scale is "
        "the weak half of the anchor. The SHAPE survives; the depletion a given cut requires "
        "does move:",
        "",
        "| Km case | Km (uM) | saturation at reference | initial marginal return "
        f"| depletion for a {pct(SENSITIVITY_CUT)} cut |",
        "| --- | --: | --: | --: | --: |",
    ]
    for s in d["km_sensitivity"]:
        L.append(f"| `{s['label']}` | {s['km_uM']:.1f} | {pct(s['saturation_at_reference'])} "
                 f"| {s['initial_marginal_return']:.2f} | {pct(s['depletion_for_cut'])} |")
    lo = min(s["depletion_for_cut"] for s in d["km_sensitivity"])
    hi = max(s["depletion_for_cut"] for s in d["km_sensitivity"])
    L += [
        "",
        f"Across a doubling and a halving of the Km the required depletion spans "
        f"{pct(lo)}-{pct(hi)} -- so the conclusion that deep depletion is needed is robust to "
        "the constant, while any specific figure is not. Measuring the Km in a pancreatic line "
        "is the cheap next step and would need no new model code.",
        "",
        ("## The spatial escape, now wired end to end"
         if st["floor_has_spatial_consumer"]
         else "## The spatial escape, now expressible in the library"
         if st["escape_can_vary_per_cell"]
         else "## What is NOT expressible, and it is wiring rather than data"),
        "",
        f"- `transsulfuration_floor` is a `Params` field: {st['floor_is_params_field']}",
        f"- it has a per-cell `Cell` counterpart: {st['floor_is_cell_field']}",
        f"- `cell.nrf2` is a per-cell field: {st['nrf2_is_cell_field']}",
        f"- the #502 arm is driveable in sim-tme-3d: {st['arm_is_spatially_reachable']}",
        f"- something writes the per-cell floor in sim-tme-3d: "
        f"{st['floor_has_spatial_consumer']}",
        f"- glutathione resynthesis sites in `biochem.rs`: {st['gsh_resynthesis_sites']}, "
        f"of which {st['resynthesis_sites_gated_by_cystine']} are multiplied by "
        f"`xc_cystine_factor()` and {st['resynthesis_sites_scaled_by_nrf2']} by `cell.nrf2`",
        "",
        "Because EVERY resynthesis site is multiplied by the cystine factor, a complete block "
        f"is zero at every NRF2 level (`nrf2_cannot_rescue_complete_block`: "
        f"{st['nrf2_cannot_rescue_complete_block']}). That is the right topology -- "
        "transsulfuration is the only cysteine source once import is gone -- and it is also "
        "why raising `cell.nrf2` is not a substitute for a per-cell escape.",
        "",
        "This matters because of a tension in the paper's own design. PMID 41422256 reports "
        "that hypoxia activates the transsulfuration pathway and so confers ferroptosis "
        "resistance under cystine limitation, and these bacteria are hypoxia-TRIGGERED: the "
        "trigger co-locates with the escape. The second bacterial paper, PMID 41704298, "
        "independently added H2S to inhibit the transsulfuration bypass -- which is the one "
        "escape this parameter structure leaves open. Read that as a CORRESPONDENCE, not a "
        "prediction: both papers were in hand before the structure was consulted, and nothing "
        "here was forecast in advance."
        + (
            " **The engine can now be asked** (#916 library, #918 reachability, #920 "
            "consumer): sim-tme-3d writes each tumour cell's floor from its own O2 "
            "supply, and the #502 arm can be driven there, so the escape is largest "
            "exactly where the trigger is. WHAT IT MEASURED IS DELIBERATELY NOT "
            "QUOTED HERE -- it comes from a spatial run this page cannot "
            "recompute, and a figure copied in from elsewhere is a figure nothing on "
            "this page can keep honest. It lives with the consumer: sim-tme-3d's own "
            "A/B tests and the layer's row in "
            "`simulations/calibration/CALIBRATION_STATUS.md`. This paragraph has now "
            "been rewritten TWICE by the engine changing under it, which is the whole "
            "reason the facts above are re-parsed rather than asserted."
            if st["floor_has_spatial_consumer"]
            else " **The library can now express it** (#916): "
            "`Cell::transsulfuration_floor` overrides the global field per cell, so "
            "the escape can be made largest exactly where the trigger is. What is "
            "still missing is a SPATIAL CONSUMER"
            + (
                ", though the #502 arm is now driveable in sim-tme-3d (#918), so the "
                "floor would no longer multiply a no-op."
                if st["arm_is_spatially_reachable"]
                else ", and the blocker is not this axis -- `erastin_xc_inhib` is "
                "reachable in no simulation binary at all, and sim-tme-3d applies no "
                "parameter overrides, so a per-cell floor there would multiply a "
                "cystine factor that no erastin dose has moved off its no-op value."
            )
            if st["escape_can_vary_per_cell"]
            else " The engine cannot state it today, because the hypoxia field is "
            "per-cell and the escape is global. Moving the floor to `Cell` is a wiring "
            "change carrying a byte-identity obligation, not a missing dataset."
        ),
        "",
        "## Verdict",
        "",
        "- the supply arm is **already expressible**, with no new axis, via `erastin_xc_inhib` "
        "plus `transsulfuration_floor`",
        f"- the concentration-to-fraction map has a **named, failable target** "
        f"(Km {d['km_uM']} uM, PMID {c['xct_cystine_km']['source_pmid']}; reference "
        f"{d['reference_uM']:.1f} uM by formulation)",
        ("- a hypoxia-graded escape is **wired end to end** (#916 library, #918 "
         "reachability, #920 spatial consumer), so it is now a question the engine can be "
         "asked rather than a gap; what it answers is reported with the consumer, not here"
         if st["floor_has_spatial_consumer"]
         else "- a hypoxia-graded escape is **expressible in the library** (#916) but has no "
         "spatial consumer"
         if st["escape_can_vary_per_cell"]
         else "- a hypoxia-graded escape is **not expressible**: a per-cell floor is missing"),
        "- the depletion-buys-less headline is **in-vitro only**: the gap exists because "
        "medium sits above the Km, and below the Km the map is nearly linear, so the "
        "in-vitro fit and the in-vivo therapy sit on opposite sides of half-saturation",
        "- the dose magnitude stays **unanchored** -- nothing here ties enzyme activity to a "
        "residual concentration in vivo, and the plasma pools in the constants table are TOTAL "
        "cyst(e)ine and free cysteine, neither of which is free cystine, so they cannot serve "
        "as the in-vivo reference",
        "",
        "No layer was added. Per `CONTRIBUTING.md`'s layer policy this names the target and "
        "the gap.",
        "",
    ]
    return "\n".join(L)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--render-only", action="store_true",
                    help="rebuild the page from the committed JSON")
    a = ap.parse_args()
    try:
        if a.render_only:
            if not OUT_JSON.exists():
                print(f"missing {OUT_JSON}", file=sys.stderr)
                return 1
            # Re-derive from the stored RAW inputs rather than echoing the stored
            # derived fields, so the guards reading those fields stay live.
            data = assemble(json.loads(OUT_JSON.read_text()))
        else:
            data = assemble(read_inputs())
            OUT_JSON.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
    except MissingInputError as e:
        print(f"refusing to render: {e}", file=sys.stderr)
        return 1
    OUT_MD.write_text(render(data))
    print(f"wrote {OUT_MD}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
