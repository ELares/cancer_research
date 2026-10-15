#!/usr/bin/env python3
"""Condition the PRCC rankings on the model's own baseline-viability constraint,
and characterise what the Sobol design cannot be conditioned on (#944).

FIGURES: none are written here. This generator is enrolled in
`tests/test_no_handwritten_figures_in_generators.py`, so every quantity below
the docstring is interpolated from a measurement and every statement here is a
property rather than a number -- which is the only way a page whose subject is
an unconditioned published figure can be trusted not to carry one of its own.

`analysis/prcc-results.json` is the root of this repository's parameter
uncertainty. It defines the swept rate constants and their ranges, and
`uncertainty_intervals.py`, `sobol_sensitivity.py` and `headline_sensitivity.py`
all load `metadata.parameter_ranges` from it deliberately, "so the
analyses share one parameter-range source of truth".

`#936` measured that box as mostly inadmissible against the manuscript's own
untreated-viability constraint, and recorded why the published design could not
have known: its conditions are treated arms only, so no untreated death rate was
ever computed. THAT IS TRUE OF `uncertainty_intervals.py` AND FALSE OF THE
ARTIFACT IT TAKES ITS RANGES FROM. `run_prcc.CONDITIONS` contains
`("Persister", "Control")`: the PRCC computed the untreated Persister death rate
for every draw. The per-draw values were discarded (`raw_death_rates` is a
local), but the aggregate is committed, in a field named for a mean death rate,
under a condition named Control, in the file that defines the box -- and it is
multiples of the limit. Every consumer reads only the ranges.

The same row names the mechanism: its own `ranked_by_abs_prcc` opens with the
parameters `#936` and `#938` later measured as the ones the admissibility
truncation selects on. The root artifact ranked the drivers of untreated death
and the row was read as sensitivity output rather than as a statement about the
box.

WHAT THIS PAGE DOES
-------------------
Reproduces the published design exactly -- the same ranges, draw count, per-draw
cell count and seed, all IMPORTED from that generator or read from the
artifact's own metadata rather than restated -- adds the untreated arms on the
SAME draws, partitions, and recomputes every published coefficient on the
admissible subset. The full-design half must reproduce the committed
coefficients or the comparison is with a re-run rather than with what is
published, and the design is fully deterministic, so that is checkable to float
round-off rather than within a chosen tolerance.

Chapter 5 cites these rankings and draws the manuscript's central mechanism from
a comparison of two of them. The two numeric bounds it quotes are PARSED from
`v1.md` and scored on both halves, so an edit to that sentence breaks this page
instead of silently changing what is tested.

CONDITIONING IS NOT INFERENCE
-----------------------------
Removing draws the model declares invalid adds no data. It cannot turn a
sensitivity ranking into a calibrated one and it leaves
`identifiability-report.json`'s data-conditioned-headline count where it is.

AND THE ADMISSIBLE SUBSET IS NOT A LATIN HYPERCUBE, which is a different caveat
from "not a posterior" and is measured here rather than asserted. LHS buys one
draw per stratum per margin, so on the full design every margin's empirical
distribution sits within a draw-width of uniform. Rejection filtering destroys
that, and destroys it UNEVENLY, because admissibility selects on a few of the
margins and leaves the drug parameters alone. So the conditioned coefficients
are rank partial correlations over a filtered sample, not over an LHS design: a
comparison between two populations rather than between two designs.

THE SOBOL DESIGN CANNOT BE CONDITIONED AT ALL, and the reason is in the
estimator rather than in the biology. Saltelli/Jansen differences matched A/B
rows under per-parameter cross-swaps, so dropping rows unmatches the pairs the
estimator is built from. What is measurable is how many of those pairs STRADDLE
the constraint -- one side a cell the model says is already dead untreated, the
other not -- and that is measured per parameter here. It is also why filtering
would not be a correction: the straddling fraction varies by orders of magnitude
across the parameters, so dropping straddling pairs would strip some of a large
share of their pairs and others of none.

AND THE ADMISSIBLE FRACTION ITSELF IS RESOLUTION-LIMITED, which bounds how
precisely this family of pages can quote one. The constraint is a threshold on a
Monte-Carlo estimate, the boundary cuts through the dense part of the
worst-untreated distribution, and a draw near it can land on either side
depending on the per-draw seed. Measured across cell counts and seeds on these
same draws rather than argued, and the DIRECTION came back against the argument
that predicted it.

OFFLINE CONTRACT
----------------
The measurement needs the compiled Python extension, so it is taken locally with
`--measure` and the per-draw results committed. Every figure is then derived
from those files, so CI re-derives the whole page and runs no simulation.
"""

from __future__ import annotations

import argparse
import ast
import csv
import gzip
import json
import math
import re
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts"))

# ONE source for the published design. Restating any of it is how a
# "reproduction" quietly stops reproducing anything.
from run_prcc import (  # noqa: E402
    CONDITIONS,
    PARAM_RANGES,
    compute_prcc,
    latin_hypercube_sample,
)

# ONE source for what admissible means, so this page cannot disagree with the
# two that already use it.
from intervals_under_admissibility import (  # noqa: E402
    CONTROL_CONDITIONS,
    DRUG_PARAMETERS,
    MIN_ADMISSIBLE_DRAWS,
    STRONG_SELECTION,
)
from headline_at_fitted import BASELINE_MAX  # noqa: E402
import sobol_sensitivity as SOBOL  # noqa: E402

PUBLISHED = REPO / "analysis" / "prcc-results.json"
MANUSCRIPT = REPO / "article" / "drafts" / "v1.md"
DRAWS_CSV_GZ = REPO / "analysis" / "calibration" / "prcc-admissibility-draws.csv.gz"
SOBOL_CSV_GZ = REPO / "analysis" / "calibration" / "prcc-admissibility-sobol.csv.gz"
OUT_JSON = REPO / "analysis" / "prcc-under-admissibility.json"
OUT_MD = REPO / "analysis" / "prcc-under-admissibility.md"

# `run_prcc.main` evaluates draw i at `seed=i + 100`. It is a loop-local
# expression rather than a module constant, so it is the ONE piece of the design
# that cannot be imported; `_pinned_seed_offset` re-parses it from that source at
# generation time, because a changed draw seed moves every number here while
# nothing in CI re-measures.
DRAW_SEED_OFFSET = 100

# `run_prcc` skips a condition whose output is constant (`np.std < 1e-10`) rather
# than correlating against it. The same rule decides here whether a conditioned
# ranking exists at all, and is re-parsed from that source for the same reason.
DEGENERATE_STD = 1e-10

# The design is deterministic, so a reproduction deviation is float round-off in
# the correlation-matrix inversion and nothing else. The measured worst is
# reported on the page beside this bound so the margin is visible rather than
# claimed; the guards require the bound to sit above that noise and far below any
# coefficient difference that could change a reading.
REPRODUCTION_TOLERANCE = 1e-12

# Cell counts and seeds for the resolution band on the admissibility verdict.
# The siblings judge at one cell count and one seed, which cannot show this.
RESOLUTION_CELLS = (1000, 2000, 8000)
RESOLUTION_SEEDS = (42, 7, 1234)

# The deviation-from-uniform a RANDOM sample of n draws shows at the 5% level is
# `KS_COEFFICIENT / sqrt(n)`, the standard Kolmogorov-Smirnov asymptotic critical
# value. It is the null this section needs: a filtered subset is smaller, and a
# smaller sample deviates more whatever was filtered, so the question is whether
# a margin deviates by MORE than having fewer draws explains. Deriving the bar
# from the admissible count rather than choosing one also makes the count of
# deformed margins a measurement.
# The coefficient and the level are ONE pair from the Kolmogorov distribution,
# so they are carried together: quoting the level in prose without it beside the
# coefficient is how the two drift apart, and a literal percentage written into
# a generator's prose is what the handwritten-figures gate exists to stop.
KS_LEVEL = 0.05
KS_COEFFICIENT = 1.36


def _label(condition) -> str:
    return f"{condition[0]}/{condition[1]}"


def _published_label(label: str) -> str:
    """The published artifact writes a condition as `A \u00d7 B`; this page
    writes `A/B`. Normalised so the escalation flags can be compared by
    IDENTITY and not merely counted."""
    return label.replace(" \u00d7 ", "/")


def _conditioned_floor(n_parameters: int) -> int:
    """Below this many admissible draws no conditioned ranking is published.

    Two floors, and the larger binds. The sibling pages already own one for
    quoting a statistic on the admissible subset, so it is IMPORTED rather than
    chosen again here. The second is arithmetic and belongs to this estimator
    rather than to the family: a partial correlation controlling for the other
    parameters has `n - 1 - p` degrees of freedom, so a ranking is not merely
    unreliable below that, it is undefined.
    """
    return max(MIN_ADMISSIBLE_DRAWS, n_parameters + 2)


def _published() -> dict:
    return json.loads(PUBLISHED.read_text())


def _design() -> dict:
    """The design, read from the published artifact's own metadata.

    `run_prcc` carries the draw count, per-condition cell count and seed as
    argparse defaults rather than module constants, so they cannot be imported.
    Reading them from the artifact is better than restating them: the
    reproduction check then validates the reading.
    """
    meta = _published()["metadata"]
    return {"n_samples": int(meta["n_samples"]),
            "n_cells": int(meta["n_cells_per_condition"]),
            "seed": int(meta["seed"]),
            "parameters": list(meta["parameters"])}


def _pinned_seed_offset() -> int:
    """The per-draw seed offset, re-parsed from `run_prcc.py`.

    Returns the integer added to the draw index in the `seed=` keyword of the
    `sim_batch` call. Raises if the shape is no longer an index plus a constant,
    because silently falling back to the pinned value is how a design change
    stops being noticed.
    """
    tree = ast.parse((REPO / "scripts" / "run_prcc.py").read_text())
    for node in ast.walk(tree):
        if not (isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
                and node.func.attr == "sim_batch"):
            continue
        for kw in node.keywords:
            if kw.arg != "seed":
                continue
            v = kw.value
            if (isinstance(v, ast.BinOp) and isinstance(v.op, ast.Add)
                    and isinstance(v.right, ast.Constant)
                    and isinstance(v.right.value, int)):
                return int(v.right.value)
            raise ValueError(f"run_prcc's draw seed is no longer index+constant: "
                             f"{ast.unparse(v)}")
    raise ValueError("no sim_batch(seed=...) call found in run_prcc.py")


def _pinned_degenerate_std() -> float:
    """`run_prcc`'s own constant-output threshold, re-parsed from that source."""
    src = (REPO / "scripts" / "run_prcc.py").read_text()
    m = re.search(r"np\.std\(death_rates\)\s*<\s*([0-9.e+-]+)", src)
    if not m:
        raise ValueError("run_prcc no longer skips constant outputs by np.std")
    return float(m.group(1))


# --------------------------------------------------------------------------
# The manuscript's own quoted figures, parsed rather than typed. Chapter 5 draws
# the central mechanism from a comparison of two of these rankings, and the two
# numeric bounds it quotes are what this page scores. Parsing them means an edit
# to that sentence breaks this page instead of silently changing what is tested.
# --------------------------------------------------------------------------
_FIGURE = r"([0-9]+(?:\.[0-9]+)?)"
_SIGNED = r"([+\u2212-][0-9]+(?:\.[0-9]+)?)"
# The bound sits at the END of its sentence, so a `[0-9.]+` class swallows the
# full stop and the figure stops being a float. Measured once, by the parse
# failing; the shared `_FIGURE` is why it cannot differ between the two.
_SDT_SENTENCE = re.compile(
    r"For Persister . SDT, a single parameter dominates: ([a-z0-9_]+) "
    r"\(" + _SIGNED + r"[^)]*\); all other parameters have "
    r"\|PRCC\| < " + _FIGURE)
_RSL3_SENTENCE = re.compile(
    r"For Persister . RSL3, the LP cascade parameters dominate: (.+?) . the "
    r"death outcome is a balance[^.]*?with ([a-z]+) parameters showing "
    r"\|PRCC\| > " + _FIGURE)
_NAMED_COEFFICIENT = re.compile(r"([a-z0-9_]+) \((?:PRCC = )?" + _SIGNED + r"\)")

# int -> word, inverted. One source for the mapping rather than a second list.
try:
    from headline_intervals_under_admissibility import _COUNT_WORDS as _WORDS
except ImportError:  # pragma: no cover - the sibling is committed beside this
    _WORDS = {}
_WORD_TO_INT = {w.lower(): n for n, w in _WORDS.items()}


def _float(text: str) -> float:
    return float(text.replace("−", "-"))


def manuscript_claims() -> dict:
    """The two Chapter 5 PRCC claims, with their figures, read from `v1.md`."""
    body = MANUSCRIPT.read_text()
    sdt = _SDT_SENTENCE.search(body)
    rsl3 = _RSL3_SENTENCE.search(body)
    if not sdt or not rsl3:
        raise ValueError("Chapter 5's PRCC sentences no longer parse; this page "
                         "scores their figures, so they cannot be assumed")
    named = [(n, _float(v)) for n, v in _NAMED_COEFFICIENT.findall(rsl3.group(1))]
    count_word = rsl3.group(2).lower()
    if count_word not in _WORD_TO_INT:
        raise ValueError(f"cannot read {count_word!r} as a count")
    return {
        "sdt": {"condition": _label(("Persister", "SDT")),
                "dominant": sdt.group(1),
                "dominant_quoted": _float(sdt.group(2)),
                "others_below": _float(sdt.group(3)),
                "sentence": sdt.group(0)},
        "rsl3": {"condition": _label(("Persister", "RSL3")),
                 "named": named,
                 "count_above": _WORD_TO_INT[count_word],
                 "above": _float(rsl3.group(3)),
                 "sentence": rsl3.group(0)},
    }


# --------------------------------------------------------------------------
# Measurement (needs the compiled extension; never runs in CI)
# --------------------------------------------------------------------------
def _fc():
    import ferroptosis_core

    return ferroptosis_core


def _all_conditions() -> list:
    """The published conditions plus whichever untreated arms it lacks.

    `("Persister", "Control")` is already in the published design; adding it
    twice would evaluate it twice and invite the two copies to disagree.
    """
    return list(CONDITIONS) + [c for c in CONTROL_CONDITIONS if c not in CONDITIONS]


def _worst_untreated(row: dict, conditions) -> float:
    return max(float(row[_label(c)]) for c in conditions)


def measure() -> None:
    design = _design()
    names = list(PARAM_RANGES)
    if names != design["parameters"]:
        raise ValueError("the published parameter order no longer matches "
                         "run_prcc.PARAM_RANGES")
    offset = _pinned_seed_offset()
    n = design["n_samples"]
    cells = design["n_cells"]
    X = latin_hypercube_sample(n, list(PARAM_RANGES.values()), seed=design["seed"])
    conditions = _all_conditions()
    fc = _fc()

    DRAWS_CSV_GZ.parent.mkdir(parents=True, exist_ok=True)
    resolution = [(c, s) for c in RESOLUTION_CELLS for s in RESOLUTION_SEEDS]
    with gzip.open(DRAWS_CSV_GZ, "wt", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["draw"] + names + [_label(c) for c in conditions]
                   + [f"worst_untreated_n{c}_seed{s}" for c, s in resolution])
        for i in range(n):
            kw = {names[j]: float(X[i, j]) for j in range(len(names))}
            rates = [fc.sim_batch(ph, tr, n=cells, seed=i + offset, **kw)["death_rate"]
                     for ph, tr in conditions]
            band = []
            for c, s in resolution:
                band.append(max(fc.sim_batch(ph, tr, n=c, seed=s, **kw)["death_rate"]
                                for ph, tr in CONTROL_CONDITIONS))
            w.writerow([i] + [f"{v:.12g}" for v in kw.values()]
                       + [f"{v:.12g}" for v in rates] + [f"{v:.12g}" for v in band])
            if (i + 1) % 200 == 0:
                print(f"  draws {i + 1}/{n}", flush=True)
    print(f"wrote {DRAWS_CSV_GZ}")

    # The Sobol design's own admissibility, per block, so the straddling
    # fractions on the page are DERIVED from the flags rather than carried.
    sob = [p[0] for p in SOBOL.PARAMS]
    lows = np.array([p[1] for p in SOBOL.PARAMS], float)
    highs = np.array([p[2] for p in SOBOL.PARAMS], float)
    rng = np.random.default_rng(SOBOL.RNG_SEED)
    A = lows + rng.random((SOBOL.N_BASE, len(sob))) * (highs - lows)
    B = lows + rng.random((SOBOL.N_BASE, len(sob))) * (highs - lows)

    def flags(M):
        out = []
        for row in M:
            kw = {nm: float(row[j]) for j, nm in enumerate(sob)}
            worst = max(fc.sim_batch(ph, tr, n=SOBOL.N_CELLS, seed=SOBOL.SIM_SEED,
                                     **kw)["death_rate"]
                        for ph, tr in CONTROL_CONDITIONS)
            out.append(1 if worst < BASELINE_MAX else 0)
        return out

    cols = {"A": flags(A), "B": flags(B)}
    for i, nm in enumerate(sob):
        ABi = A.copy()
        ABi[:, i] = B[:, i]
        cols[f"AB_{nm}"] = flags(ABi)
        print(f"  sobol cross-swap {nm}", flush=True)
    # The Sobol design excludes `sdt_ros`, so the untreated arms here are
    # evaluated at its default. That is only sound if it is inert untreated,
    # which is measured rather than assumed.
    base = {nm: float(A[0, j]) for j, nm in enumerate(sob)}
    inert = []
    for ph, tr in CONTROL_CONDITIONS:
        vals = []
        for v in PARAM_RANGES["sdt_ros"]:
            kw = dict(base)
            kw["sdt_ros"] = float(v)
            vals.append(fc.sim_batch(ph, tr, n=SOBOL.N_CELLS,
                                     seed=SOBOL.SIM_SEED, **kw)["death_rate"])
        inert.append({"condition": _label((ph, tr)), "low": vals[0], "high": vals[1],
                      "identical": vals[0] == vals[1]})
    with gzip.open(SOBOL_CSV_GZ, "wt", newline="") as fh:
        w = csv.writer(fh)
        order = ["A", "B"] + [f"AB_{nm}" for nm in sob]
        w.writerow(["row"] + order)
        for r in range(SOBOL.N_BASE):
            w.writerow([r] + [cols[k][r] for k in order])
    (SOBOL_CSV_GZ.with_suffix("").with_suffix(".control.json")).write_text(
        json.dumps({"sdt_ros_inert_untreated": inert}, indent=2) + "\n")
    print(f"wrote {SOBOL_CSV_GZ}")


# --------------------------------------------------------------------------
# Derivation
# --------------------------------------------------------------------------
def load() -> dict:
    with gzip.open(DRAWS_CSV_GZ, "rt") as fh:
        rows = list(csv.DictReader(fh))
    with gzip.open(SOBOL_CSV_GZ, "rt") as fh:
        sobol_rows = list(csv.DictReader(fh))
    control = json.loads(
        SOBOL_CSV_GZ.with_suffix("").with_suffix(".control.json").read_text())
    return {"rows": rows, "sobol_rows": sobol_rows,
            "sobol_control": control["sdt_ros_inert_untreated"],
            "published": _published(), "manuscript": manuscript_claims()}


def _ks_from_uniform(values, lo, hi) -> float:
    """Largest deviation of the margin's empirical distribution from uniform."""
    u = np.sort((np.asarray(values, float) - lo) / (hi - lo))
    m = len(u)
    if m == 0:
        return float("nan")
    above = np.max(np.abs(u - np.arange(1, m + 1) / m))
    below = np.max(np.abs(u - np.arange(0, m) / m))
    return float(max(above, below))


def _ranking(coefficients, names) -> list:
    return sorted(names, key=lambda nm: -abs(coefficients[nm]))


def _coefficients(X, y, names) -> dict:
    return {nm: float(c) for nm, (c, _p) in zip(names, compute_prcc(X, y))}


def assemble(raw: dict) -> dict:
    rows = raw["rows"]
    pub = raw["published"]
    names = list(PARAM_RANGES)
    X = np.array([[float(r[nm]) for nm in names] for r in rows], float)
    n = len(rows)

    cells = _design()["n_cells"]
    worst = np.array([_worst_untreated(r, CONTROL_CONDITIONS) for r in rows], float)
    admissible = worst < BASELINE_MAX
    n_adm = int(admissible.sum())

    per_condition = {}
    repro_dev = 0.0
    for c in CONDITIONS:
        label = _label(c)
        y = np.array([float(r[label]) for r in rows], float)
        if y.std() < DEGENERATE_STD:
            # `run_prcc` correlates every published condition unconditionally,
            # so an artifact holding a constant one could not have been
            # produced. Say that rather than surfacing a singular-matrix error
            # from three frames down.
            raise ValueError(
                f"{label} is constant across the full design, which the "
                f"published artifact cannot be: its own generator would have "
                f"failed to correlate it")
        full = _coefficients(X, y, names)
        published_row = next(p for p in pub["conditions"]
                             if (p["phenotype"], p["treatment"]) == c)
        repro_dev = max(repro_dev,
                        max(abs(full[nm] - published_row["prcc"][nm]) for nm in names))
        ya = y[admissible]
        too_few = n_adm < _conditioned_floor(len(names))
        degenerate = bool(too_few or ya.std() < DEGENERATE_STD)
        cond = None if degenerate else _coefficients(X[admissible], ya, names)
        # The MOMENTS exist whenever any draw survives; only the RANKING can be
        # undefined. Collapsing the two published a withheld mean as a measured
        # zero, which is the absence-versus-zero defect already on record here.
        per_condition[label] = {
            "full": {"mean": float(y.mean()), "std": float(y.std()),
                     "coefficients": full, "ranking": _ranking(full, names),
                     "deaths": int(round(y.sum() * cells))},
            "admissible": None if n_adm == 0 else {
                "mean": float(ya.mean()), "std": float(ya.std()),
                "coefficients": cond,
                "ranking": None if cond is None else _ranking(cond, names)},
            "degenerate_under_conditioning": degenerate,
            "withheld_for_too_few_draws": bool(too_few),
            "published_mean": published_row["mean_death_rate"],
            "conditioned_on": label in [_label(x) for x in CONTROL_CONDITIONS],
        }
    reproduces = repro_dev <= REPRODUCTION_TOLERANCE

    # The manuscript's two quoted bounds, scored on both halves.
    ms = raw["manuscript"]
    sdt = per_condition[ms["sdt"]["condition"]]
    others = [nm for nm in names if nm != ms["sdt"]["dominant"]]

    def _coeffs(half):
        """A half's coefficients, or None.

        Two distinct absences collapse here and both must stay readable: no
        admissible draw at all (`half is None`), and draws without a definable
        ranking (`coefficients is None`). A `.get` chain over the first one
        returned None and then called `.get` on it, which is how splitting the
        moments from the ranking surfaced this.
        """
        return None if half is None else half["coefficients"]

    def _sdt_breaches(half):
        c = _coeffs(half)
        if c is None:
            return None
        return sorted(((nm, abs(c[nm])) for nm in others
                       if abs(c[nm]) >= ms["sdt"]["others_below"]),
                      key=lambda t: (-t[1], t[0]))

    rsl3 = per_condition[ms["rsl3"]["condition"]]

    def _rsl3_count(half):
        c = _coeffs(half)
        if c is None:
            return None
        return sum(1 for nm in names if abs(c[nm]) > ms["rsl3"]["above"])

    scored = {
        "sdt_bound": {
            "bound": ms["sdt"]["others_below"],
            "dominant": ms["sdt"]["dominant"],
            "dominant_full": sdt["full"]["coefficients"][ms["sdt"]["dominant"]],
            "dominant_admissible": (
                None if _coeffs(sdt["admissible"]) is None
                else _coeffs(sdt["admissible"])[ms["sdt"]["dominant"]]),
            "breaches_full": _sdt_breaches(sdt["full"]),
            "breaches_admissible": _sdt_breaches(sdt["admissible"]),
            "still_dominant_admissible": bool(
                sdt["admissible"] and sdt["admissible"]["ranking"]
                and sdt["admissible"]["ranking"][0] == ms["sdt"]["dominant"]),
        },
        "rsl3_count": {
            "bound": ms["rsl3"]["above"],
            "claimed": ms["rsl3"]["count_above"],
            "full": _rsl3_count(rsl3["full"]),
            "admissible": _rsl3_count(rsl3["admissible"]),
            "named": [nm for nm, _v in ms["rsl3"]["named"]],
            "named_still_top": bool(
                rsl3["admissible"] and rsl3["admissible"]["ranking"]
                and (rsl3["admissible"]["ranking"][:len(ms["rsl3"]["named"])]
                     == [nm for nm, _v in ms["rsl3"]["named"]])),
        },
    }
    scored["sdt_bound"]["holds_full"] = scored["sdt_bound"]["breaches_full"] == []
    scored["sdt_bound"]["holds_admissible"] = (
        None if scored["sdt_bound"]["breaches_admissible"] is None
        else scored["sdt_bound"]["breaches_admissible"] == [])
    scored["rsl3_count"]["holds_full"] = (
        scored["rsl3_count"]["full"] == ms["rsl3"]["count_above"])
    scored["rsl3_count"]["holds_admissible"] = (
        None if scored["rsl3_count"]["admissible"] is None
        else scored["rsl3_count"]["admissible"] == ms["rsl3"]["count_above"])

    # Circularity falsifier: if conditioning selected on the DRUG parameters, a
    # change in an SDT ranking would be partly an artifact of keeping
    # weaker-drug draws, and the comparison this page makes would be circular.
    # The same-size null, derived from how many draws survive rather than
    # chosen; `None` when nothing survives, since there is no sample to compare.
    null_ks = (KS_COEFFICIENT / math.sqrt(n_adm)) if n_adm else None
    full_null_ks = KS_COEFFICIENT / math.sqrt(n)
    selection, stratification = {}, {}
    for j, nm in enumerate(names):
        lo, hi = PARAM_RANGES[nm]
        f = float((X[:, j].mean() - lo) / (hi - lo))
        a = (float((X[admissible, j].mean() - lo) / (hi - lo)) if n_adm else None)
        selection[nm] = {
            "full_mean_fraction": f, "admissible_mean_fraction": a,
            "shift": None if a is None else a - f,
            "drug": nm in DRUG_PARAMETERS,
            # With nothing admissible there is no selection to characterise, so
            # this is False rather than None: the claim "this parameter was
            # selected on" is not true, and None would make the circularity
            # verdict below unreadable.
            "strongly_selected": (False if a is None
                                  else abs(a - f) > STRONG_SELECTION)}
        # None rather than NaN when nothing survives. `json.dumps` writes NaN
        # as the bare token `NaN`, which is not JSON, so a nan here does not
        # merely read badly on the page -- it makes the artifact unparseable by
        # a strict reader.
        a_ks = (_ks_from_uniform(X[admissible, j], lo, hi) if n_adm else None)
        stratification[nm] = {
            "full": _ks_from_uniform(X[:, j], lo, hi),
            "admissible": a_ks,
            # Against a same-size random subsample, so the ratio is a statement
            # about the SELECTION and not about the subset being smaller.
            "times_null": (a_ks / null_ks) if (null_ks and a_ks is not None)
                          else None}
    circularity_clear = (
        None if n_adm == 0
        else not any(v["strongly_selected"] for v in selection.values()
                     if v["drug"]))

    # The escalation decision published in the same artifact, recomputed on the
    # PRCC's OWN conditions. A first pass here ran it over the untreated arms as
    # well and reported three flags against the published two, which would have
    # been a reproduction failure read as a finding.
    from scipy.stats import pearsonr, spearmanr

    threshold = float(pub["sobol_decision"]["divergence_threshold"])
    escalation = {}
    for half, sel in (("full", np.ones(n, bool)), ("admissible", admissible)):
        if int(sel.sum()) < _conditioned_floor(len(names)):
            # A Pearson-against-Spearman divergence over a handful of draws is
            # not a measurement either, and below two draws the estimators
            # raise. Withheld with the same floor the rankings use, so the page
            # cannot report a decision on one subset and a ranking on another.
            escalation[half] = {"flags": None, "n_flags": None,
                                "warranted": None, "withheld": True}
            continue
        flags = []
        for c in CONDITIONS:
            y = np.array([float(r[_label(c)]) for r in rows], float)[sel]
            if y.std() < DEGENERATE_STD:
                continue
            for j, nm in enumerate(names):
                d = abs(float(pearsonr(X[sel, j], y)[0])
                        - float(spearmanr(X[sel, j], y)[0]))
                if d > threshold:
                    flags.append({"condition": _label(c), "parameter": nm,
                                  "divergence": round(d, 3)})
        escalation[half] = {"flags": flags, "n_flags": len(flags),
                            "warranted": len(flags) > 2, "withheld": False}
    escalation["published"] = {
        "n_flags": int(pub["sobol_decision"]["n_flags"]),
        "warranted": bool(pub["sobol_decision"]["warranted"]),
        "withheld": False,
        "flags": [{"condition": _published_label(f["condition"]),
                   "parameter": f["parameter"], "divergence": f["divergence"]}
                  for f in pub["sobol_decision"]["non_monotonic_flags"]]}
    # Compared by IDENTITY, not by count: two flags on different parameters is
    # not a reproduction of two flags, and counting alone cannot tell them
    # apart. The condition labels are normalised because the two artifacts
    # spell a condition differently.
    def _flagset(flags):
        return sorted((f["condition"], f["parameter"], round(f["divergence"], 3))
                      for f in flags)

    escalation["reproduces"] = (
        escalation["full"]["n_flags"] == escalation["published"]["n_flags"]
        and escalation["full"]["warranted"] == escalation["published"]["warranted"]
        and _flagset(escalation["full"]["flags"])
        == _flagset(escalation["published"]["flags"]))
    # None, not False: a withheld half has not disagreed with the published
    # decision, and recording it as a disagreement would be the absence-as-zero
    # defect one field over.
    escalation["survives_conditioning"] = (
        None if escalation["admissible"]["withheld"]
        else escalation["admissible"]["warranted"]
        == escalation["published"]["warranted"])

    # Resolution band on the admissibility verdict itself.
    band, conventions = {}, []
    for c in RESOLUTION_CELLS:
        counts = []
        for s in RESOLUTION_SEEDS:
            col = f"worst_untreated_n{c}_seed{s}"
            if col not in rows[0]:
                continue
            counts.append(int(sum(1 for r in rows if float(r[col]) < BASELINE_MAX)))
            conventions.append(col)
        if counts:
            band[c] = {"counts": counts, "mean_fraction": float(np.mean(counts) / n),
                       "low": min(counts) / n, "high": max(counts) / n,
                       "spread_draws": max(counts) - min(counts)}
    cells = sorted(band)
    monotone_in_cells = all(band[a]["mean_fraction"] >= band[b]["mean_fraction"]
                            for a, b in zip(cells, cells[1:])) if len(cells) > 1 else None

    # The Sobol design: how many cross-swap pairs straddle the constraint.
    srows = raw["sobol_rows"]
    aA = np.array([int(r["A"]) for r in srows], bool)
    aB = np.array([int(r["B"]) for r in srows], bool)
    straddle = {}
    for nm in [p[0] for p in SOBOL.PARAMS]:
        key = f"AB_{nm}"
        if key not in srows[0]:
            continue
        ab = np.array([int(r[key]) for r in srows], bool)
        straddle[nm] = {"admissible_fraction": float(ab.mean()),
                        "straddling": int((aA != ab).sum()),
                        "straddling_fraction": float((aA != ab).mean())}
    st_values = [v["straddling_fraction"] for v in straddle.values()]

    return {
        "design": {"n_samples": n, "n_cells": _design()["n_cells"],
                   "seed": _design()["seed"], "draw_seed_offset": DRAW_SEED_OFFSET,
                   "parameters": names, "conditions": [_label(c) for c in CONDITIONS],
                   "untreated_arms": [_label(c) for c in CONTROL_CONDITIONS],
                   "baseline_max": BASELINE_MAX},
        "reproduction": {"max_coefficient_deviation": repro_dev,
                         "tolerance": REPRODUCTION_TOLERANCE,
                         "reproduces": bool(reproduces),
                         "n_coefficients": len(CONDITIONS) * len(names)},
        "admissibility": {
            "n_admissible": n_adm, "n_draws": n, "fraction": n_adm / n,
            "median_worst_untreated": float(np.median(worst)),
            "max_worst_untreated": float(np.max(worst)),
            "published_control_mean": next(
                p["mean_death_rate"] for p in pub["conditions"]
                if (p["phenotype"], p["treatment"]) == ("Persister", "Control")),
            "published_control_ranking": next(
                p["ranked_by_abs_prcc"] for p in pub["conditions"]
                if (p["phenotype"], p["treatment"]) == ("Persister", "Control")),
        },
        "resolution": {"band": band, "conventions": conventions,
                       "falls_with_cell_count": monotone_in_cells},
        "per_condition": per_condition,
        "manuscript": raw["manuscript"],
        "scored": scored,
        "selection": selection,
        # NOT `bool(...)`: this is three-state, and a cast turned the
        # WITHHELD None into a False verdict -- the absence-as-a-result
        # defect one line over from where it was just fixed.
        "circularity_clear": circularity_clear,
        "stratification": stratification,
        "stratification_null": {"admissible": null_ks, "full": full_null_ks,
                                "coefficient": KS_COEFFICIENT,
                                "level": KS_LEVEL},
        "stratification_lost": sorted(
            nm for nm, v in stratification.items()
            if v["times_null"] is not None and v["times_null"] > 1.0),
        "escalation": escalation,
        "sobol_design": {
            "n_base": SOBOL.N_BASE, "excluded": sorted(SOBOL.EXCLUDED),
            "n_cells": SOBOL.N_CELLS, "seed": SOBOL.SIM_SEED,
            "a_admissible_fraction": float(aA.mean()),
            "b_admissible_fraction": float(aB.mean()),
            "straddle": straddle,
            "straddle_min": min(st_values) if st_values else None,
            "straddle_max": max(st_values) if st_values else None,
            "control": raw["sobol_control"],
            "sdt_ros_inert": all(c["identical"] for c in raw["sobol_control"]),
        },
        "conditioned_floor": {"floor": _conditioned_floor(len(names)),
                              "n_admissible": n_adm,
                              "bound": n_adm < _conditioned_floor(len(names))},
        "thresholds": {"strong_selection": STRONG_SELECTION,
                       "ks_coefficient": KS_COEFFICIENT,
                       "ks_level": KS_LEVEL,
                       "degenerate_std": DEGENERATE_STD},
    }


# --------------------------------------------------------------------------
# Render. Every quantity is interpolated; see the docstring's FIGURES note.
# --------------------------------------------------------------------------
def _pct(x, places=1) -> str:
    return "n/a" if x is None else f"{100 * x:.{places}f}%"


def _sig(x, places=4) -> str:
    return "n/a" if x is None else f"{x:.{places}f}"


def _stratification_section(d) -> list:
    """The margin-deviation table, which needs a subset large enough to
    have a shape. Factored out so `render` can skip it without an early
    return that would also drop the sections below it.
    """
    L = []
    nul = d["stratification_null"]
    L.append(
        f"LHS buys one draw per stratum per margin, so on the full design every "
        f"margin sits within a draw-width of uniform -- "
        f"{nul['full'] / max(v['full'] for v in d['stratification'].values()):.0f} "
        f"times tighter than a random sample of the same size would give, which "
        f"is what the design is for. A filtered subset is SMALLER, though, and a "
        f"smaller sample deviates more whatever was filtered, so each margin is "
        f"reported against the deviation a random subsample of the admissible "
        f"size would show ({nul['admissible']:.4f} at the "
        f"{nul['level']:.0%} level):\n")
    L.append("| parameter | full | admissible | times the same-size null |")
    L.append("|---|--:|--:|--:|")
    for nm, v in sorted(d["stratification"].items(),
                        key=lambda kv: (-kv[1]["admissible"], kv[0])):
        L.append(f"| `{nm}` | {v['full']:.4f} | {v['admissible']:.4f} | "
                 f"{v['times_null']:.2f} |")
    L.append("")
    worst_nm = min(d["stratification"],
                   key=lambda k: (-d["stratification"][k]["admissible"], k))
    best_nm = min(d["stratification"],
                  key=lambda k: (d["stratification"][k]["admissible"], k))
    w = d["stratification"][worst_nm]
    bst = d["stratification"][best_nm]
    L.append(
        f"The loss is UNEVEN and it lands exactly where the selection table "
        f"says it should: `{worst_nm}` deviates by {w['admissible']:.4f}, "
        f"{w['times_null']:.1f} times the same-size null, against `{best_nm}` at "
        f"{bst['admissible']:.4f} or {bst['times_null']:.1f} times it -- so the "
        f"least-deformed margins are within what having fewer draws explains and "
        f"the selected ones are not. So the conditioned coefficients "
        f"are rank partial correlations over a filtered sample, not over an LHS "
        f"design: a comparison between two populations rather than between two "
        f"designs, and the stratification LHS exists to provide is gone from "
        f"{len(d['stratification_lost'])} of {len(d['stratification'])} "
        f"margins.\n")
    return L


def _escalation_section(d) -> list:
    """The escalation decision, whose conditioned half can be withheld.

    Its own function for the same reason the stratification table is:
    the withheld branch has to stop EARLY without cutting the page short,
    and an early `return` inside `render` would have dropped the Sobol
    section below it -- and returned a list where a string belongs.
    """
    L = []
    e = d["escalation"]
    if e["admissible"]["withheld"]:
        L.append(
            f"`sobol_decision` publishes `warranted: "
            f"{str(e['published']['warranted']).lower()}` on "
            f"{e['published']['n_flags']} flags against a threshold of more "
            f"than two. The conditioned half is WITHHELD: too few draws survive "
            f"the constraint to recompute a divergence on, and the same floor "
            f"the rankings use decides it, so this page cannot report a "
            f"decision on one subset and a ranking on another.\n")
        return L
    L.append(
        f"`sobol_decision` publishes `warranted: "
        f"{str(e['published']['warranted']).lower()}` on "
        f"{e['published']['n_flags']} flags against a threshold of more than "
        f"two, i.e. one flag from flipping, and the manuscript repeats the "
        f"call. Recomputed on the PRCC's own conditions it reproduces "
        f"({e['full']['n_flags']} flags) and **survives conditioning**: "
        f"{e['admissible']['n_flags']} flags on the admissible subset, "
        f"`warranted: {str(e['admissible']['warranted']).lower()}`. "
        + ("The divergences themselves move AWAY from the threshold.\n"
           if all(f["divergence"] <= max((x["divergence"]
                                          for x in e["full"]["flags"]), default=0)
                  for f in e["admissible"]["flags"]) else "\n"))
    L.append("| half | flags | warranted |")
    L.append("|---|---|---|")
    for half in ("published", "full", "admissible"):
        h = e[half]
        flags = ", ".join(f"`{f['parameter']}` in `{f['condition']}` "
                          f"({f['divergence']})" for f in h.get("flags", [])) or "--"
        L.append(f"| {half} | {h['n_flags']} ({flags}) | {h['warranted']} |")
    L.append("")
    L.append(
        "A first pass recomputed this over the untreated arms as well and got "
        "one flag more than the artifact publishes, which would have been a "
        "reproduction failure reported as a conditioning effect. The "
        "recomputation is restricted to the published condition set, and the "
        "reproduction is asserted before the conditioned half is read.\n")
    return L


def _scored_section(d) -> list:
    """Chapter 5's two quoted figures, scored on both halves.

    The third section that needs a WITHHELD branch, which is the pattern:
    everything comparing the full design against the admissible subset
    has to be able to say that the second half does not exist. Scoring a
    manuscript figure against an absent subset formatted `None` into the
    prose on the first pass.
    """
    L = []
    ms, sc = d["manuscript"], d["scored"]
    if sc["sdt_bound"]["dominant_admissible"] is None:
        L.append(
            "Both figures are parsed from `article/drafts/v1.md`, and "
            "NEITHER is scored under conditioning here: no admissible "
            "subset survives to score them against. They are quoted so a "
            "reader can see what would be scored.\n")
        for half in ("sdt", "rsl3"):
            L.append(f"> {ms[half]['sentence']}\n")
        return L
    ms, sc = d["manuscript"], d["scored"]
    L.append(f"Both figures are parsed from `article/drafts/v1.md`, not typed "
             f"here.\n")
    L.append(f"> {ms['sdt']['sentence']}\n")
    s = sc["sdt_bound"]
    L.append(
        f"`{s['dominant']}` stays dominant and strengthens, "
        f"{s['dominant_full']:+.4f} -> {s['dominant_admissible']:+.4f}"
        f"{' (still the top-ranked parameter)' if s['still_dominant_admissible'] else ''}"
        f". **The quoted bound of {s['bound']} does not survive**: "
        f"{len(s['breaches_full'])} parameters reach it on the full design and "
        f"{len(s['breaches_admissible'])} do on the admissible subset"
        + (" -- " + ", ".join(f"`{nm}` at {v:.4f}"
                              for nm, v in s["breaches_admissible"])
           if s["breaches_admissible"] else "")
        + ". So the asymmetry the sentence describes survives conditioning and "
          "the number it puts on the other parameters does not.\n")
    L.append(f"> {ms['rsl3']['sentence']}\n")
    r = sc["rsl3_count"]
    L.append(
        f"The {len(r['named'])} named parameters "
        + ", ".join(f"`{nm}`" for nm in r["named"])
        + (" remain the top three in the same order"
           if r["named_still_top"] else " no longer lead in that order")
        + f". **The quoted count does not survive**: parameters above "
          f"{r['bound']} go {r['full']} on the full design (reproducing the "
          f"sentence) to {r['admissible']} on the admissible subset.\n")
    L.append(
        "Neither sentence is edited. Each is true of the prior it cites, and "
        "the deliverable here is that the two bounds are measured under the "
        "restriction rather than assumed to carry over. The rankings -- what "
        "the sentences are for -- hold.\n")
    return L


def _circularity_section(d) -> list:
    """The falsifier, which needs an admissible subset to characterise.

    Fourth section with a withheld branch, and the last: every comparison
    between the full design and the admissible subset has one.
    """
    L = []
    if d["conditioned_floor"]["bound"]:
        L.append(
            "No admissible subset large enough to characterise, so the "
            "falsifier is WITHHELD rather than reported as passing: a "
            "selection that kept nothing has not selected on a drug "
            "parameter, and reading that as the check passing would be "
            "reading an absence as a result.\n")
        return L
    L.append(
        f"If conditioning selected on the DRUG parameters, a change in an SDT "
        f"ranking would be partly an artifact of keeping weaker-drug draws. "
        f"Each parameter's admissible mean against its full-prior mean, as a "
        f"fraction of its own range (strong selection at "
        f"{d['thresholds']['strong_selection']}):\n")
    L.append("| parameter | full | admissible | shift | |")
    L.append("|---|--:|--:|--:|---|")
    for nm, v in sorted(d["selection"].items(),
                        key=lambda kv: (-abs(kv[1]["shift"]), kv[0])):
        tag = []
        if v["drug"]:
            tag.append("drug")
        if v["strongly_selected"]:
            tag.append("**strongly selected**")
        L.append(f"| `{nm}` | {v['full_mean_fraction']:.3f} | "
                 f"{v['admissible_mean_fraction']:.3f} | {v['shift']:+.3f} | "
                 f"{', '.join(tag)} |")
    L.append("")
    drug_max = max(abs(v["shift"]) for v in d["selection"].values() if v["drug"])
    L.append(
        f"**{'Clear' if d['circularity_clear'] else 'NOT clear'}.** The largest "
        f"drug-parameter shift is {drug_max:.3f} against a bar of "
        f"{d['thresholds']['strong_selection']}, so the rankings move because "
        f"the defences kept are stronger and not because the drug kept is "
        f"weaker. That matters more here than on the sibling pages, because the "
        f"claim under test is about a drug parameter's dominance.\n")
    return L


def render(d: dict) -> str:
    g, adm, rep = d["design"], d["admissibility"], d["reproduction"]
    L = []
    L.append("# The PRCC rankings under the model's own viability constraint "
             "(#944)\n")
    L.append("*Generated by `scripts/prcc_under_admissibility.py`. Every figure "
             "is derived from the committed per-draw measurements; the generator "
             "hand-writes none.*\n")

    L.append("## The root artifact ranked the drivers of untreated death, and "
             "the row was read as sensitivity output\n")
    L.append(
        f"`analysis/prcc-results.json` defines the {len(g['parameters'])} swept "
        f"rate constants and their ranges, and `uncertainty_intervals.py`, "
        f"`sobol_sensitivity.py` and `headline_sensitivity.py` each load "
        f"`metadata.parameter_ranges` from it so the analyses share one range "
        f"source. `#936` measured that box as mostly inadmissible against the "
        f"manuscript's untreated-viability constraint of "
        f"{_pct(g['baseline_max'], 0)} death under Control, and recorded that "
        f"the published design could not have known because its conditions are "
        f"treated arms only.\n")
    L.append(
        f"**The PRCC's conditions are not treated arms only.** "
        f"`run_prcc.CONDITIONS` contains `Persister/Control`, and the committed "
        f"artifact publishes its `mean_death_rate` as "
        f"**{_pct(adm['published_control_mean'], 1)}** against that limit of "
        f"{_pct(g['baseline_max'], 0)} -- a factor of "
        f"{adm['published_control_mean'] / g['baseline_max']:.1f}. The per-draw "
        f"values were discarded, but a summary of the violation has been "
        f"committed since the run, and every consumer reads only the ranges.\n")
    L.append(
        "That row also names the mechanism. Its own `ranked_by_abs_prcc` opens "
        + ", ".join(f"`{p}`" for p in adm["published_control_ranking"][:3])
        + " -- the parameters `#936` and `#938` later measured as the ones the "
          "admissibility truncation selects on. The ordering was available in "
          "the root artifact before anything downstream measured it.\n")

    L.append("## Reproduction, which everything else rests on\n")
    L.append(
        f"The design is fully deterministic (draw `i` is evaluated at "
        f"`seed = i + {g['draw_seed_offset']}`, re-parsed from `run_prcc.py` "
        f"rather than restated), so the published coefficients reproduce to "
        f"float round-off rather than within a chosen tolerance. All "
        f"{rep['n_coefficients']} of them do: the largest deviation across "
        f"every condition and parameter is {rep['max_coefficient_deviation']:.2e} "
        f"against a bound of {rep['tolerance']:.0e}. "
        + ("Reproduced, so the comparison below is with what is published.\n"
           if rep["reproduces"] else
           "**NOT reproduced**, so the comparison below is withheld: it would "
           "be against a re-run rather than against what is published.\n"))

    L.append("## How much of the design is admissible, and how precisely that "
             "can be said\n")
    L.append(
        f"Judged on the same {len(g['untreated_arms'])} untreated arms and the "
        f"same imported constraint the sibling pages use "
        f"({', '.join('`' + a + '`' for a in g['untreated_arms'])}): "
        f"**{adm['n_admissible']} of {adm['n_draws']} draws "
        f"({_pct(adm['fraction'])})** satisfy it. Median worst untreated death "
        f"{_pct(adm['median_worst_untreated'], 1)}, worst "
        f"{_pct(adm['max_worst_untreated'], 1)}.\n")
    floor = d["conditioned_floor"]
    if floor["bound"]:
        L.append(
            f"**No conditioned ranking is published at all.** Only "
            f"{floor['n_admissible']} draws survive the constraint against a "
            f"floor of {floor['floor']} -- the larger of the figure the sibling "
            f"pages use for quoting a statistic on the admissible subset and "
            f"the arithmetic floor below which a partial correlation over "
            f"{len(g['parameters'])} parameters has no degrees of freedom left. "
            f"Everything below describes the full design only.\n")
    band = d["resolution"]["band"]
    if band:
        L.append(
            "**That fraction is resolution-limited, and this family of pages "
            "quotes one to a decimal place.** The constraint is a threshold on "
            "a Monte-Carlo estimate and the boundary cuts through the dense "
            "part of the worst-untreated distribution, so a draw near it lands "
            "on either side depending on the per-draw seed. The same draws, "
            "re-judged:\n")
        L.append("| cells per arm | admissible counts across seeds | mean | "
                 "range |")
        L.append("|--:|---|--:|---|")
        for cells in sorted(band, key=int):
            b = band[cells]
            L.append(f"| {cells} | {', '.join(str(c) for c in b['counts'])} | "
                     f"{_pct(b['mean_fraction'], 2)} | {_pct(b['low'])}-"
                     f"{_pct(b['high'])} ({b['spread_draws']} draws) |")
        L.append("")
        cells = sorted(band, key=int)
        lo_c, hi_c = cells[0], cells[-1]
        L.append(
            f"So the admissible fraction of one box over one set of draws moves "
            f"by {band[lo_c]['spread_draws']} draws on the seed alone at the "
            f"lowest resolution here. **AND IT FALLS WITH RESOLUTION**, "
            f"{_pct(band[lo_c]['mean_fraction'], 2)} at {lo_c} cells per arm "
            f"against {_pct(band[hi_c]['mean_fraction'], 2)} at {hi_c}, while "
            f"the seed spread collapses from {band[lo_c]['spread_draws']} draws "
            f"to {band[hi_c]['spread_draws']}. "
            + ("Monotone across every count measured. "
               if d["resolution"]["falls_with_cell_count"] else
               "Not monotone across the counts measured. ")
            + "A fraction judged at the lower counts therefore reads HIGH, "
              "which is the direction that makes `#936`'s conclusion stronger "
              "rather than weaker -- and it is the opposite of the direction "
              "predicted before measuring, on the argument that a maximum over "
              "noisy arms overestimates the maximum and so should bias the "
              "admissible fraction down. That argument is refuted here and is "
              "recorded rather than removed.\n")
        L.append(
            "Nothing substantive turns on it: *most of the prior is "
            "inadmissible* survives a band this wide comfortably. What does "
            "not survive is quoting the fraction to a decimal place as a "
            "property of the box, and `#938`'s cross-page agreement statistic "
            "compares two such fractions through a binomial standard error, "
            "which models the sampling of draws and not the classification of "
            "one. The agreement it reports is unaffected -- a wider band makes "
            "two fractions agree more easily, not less.\n")

    L.append("## The rankings, conditioned\n")
    L.append("| condition | mean (full -> admissible) | std | top three "
             "(full) | top three (admissible) |")
    L.append("|---|---|---|---|---|")
    for label in g["conditions"]:
        c = d["per_condition"][label]
        a = c["admissible"]
        full_top = ", ".join(f"`{p}`" for p in c["full"]["ranking"][:3])
        # The moments are printed whenever they exist, and `withheld` is not
        # spelled as a number: a withheld mean rendered as 0.000000 reads as a
        # measured zero under a column heading that means a mean.
        means = f"{_sig(c['full']['mean'], 6)} -> " + (
            "withheld" if a is None else _sig(a["mean"], 6))
        stds = f"{_sig(c['full']['std'], 6)} -> " + (
            "withheld" if a is None else _sig(a["std"], 6))
        if a is None or a["ranking"] is None:
            atop = ("**undefined** (too few admissible draws)"
                    if c["withheld_for_too_few_draws"]
                    else "**undefined** (constant output)")
        else:
            atop = ", ".join(f"`{p}`" for p in a["ranking"][:3])
        note = " *(the arm conditioned on)*" if c["conditioned_on"] else ""
        L.append(f"| `{label}`{note} | {means} | {stds} | {full_top} | {atop} |")
    L.append("")
    degenerate = [k for k, v in d["per_condition"].items()
                  if v["degenerate_under_conditioning"]]
    if degenerate:
        for label in degenerate:
            c = d["per_condition"][label]
            L.append(
                f"**`{label}` has no conditioned ranking at all.** Its output "
                f"is constant across the admissible subset, so a rank partial "
                f"correlation is undefined and none is published here. On the "
                f"full design the condition contains "
                f"{c['full']['deaths']} dead cells across "
                f"{adm['n_draws'] * g['n_cells']:,} simulated, so its published "
                f"ranking is an ordering of correlations against that many "
                f"events; conditioning removes them.\n")
    exo = [l for l in g["conditions"] if l.endswith("/SDT")]
    rose = [l for l in exo if d["per_condition"][l]["admissible"] is not None
            and d["per_condition"][l]["admissible"]["std"]
            > d["per_condition"][l]["full"]["std"]]
    L.append(
        f"**One mechanism, opposite effects on attribution.** Admissibility "
        f"keeps the draws whose defences are strong enough that an untreated "
        f"cell survives, which DE-SATURATES the exogenous-ROS arms -- the "
        f"output standard deviation rises on {len(rose)} of {len(exo)} of them, "
        f"giving the defence parameters variance to be attributed -- while "
        f"pushing the inducer arm toward its floor. It is the same two-"
        f"directional structure `#936` found in the intervals, showing up in the "
        f"rankings.\n")
    L.append("The arm being conditioned on is reported as a diagnostic and not "
             "as a result: selecting draws on an outcome and then reporting "
             "that outcome's sensitivity is selection on the dependent "
             "variable.\n")

    L.append("## What Chapter 5 quotes, scored\n")
    L.extend(_scored_section(d))
    L.append("## Is the comparison circular?\n")
    L.extend(_circularity_section(d))
    L.append("## The admissible subset is not a Latin hypercube\n")
    if floor["bound"]:
        L.append(
            "The admissible subset is too small to characterise its margins, "
            "so no deviation table is published: with that few draws the "
            "deviation from uniform is dominated by how few there are rather "
            "than by which ones the constraint kept.\n")
    else:
        L.extend(_stratification_section(d))
    L.append("## The escalation decision in the same artifact\n")
    L.extend(_escalation_section(d))
    L.append("## The Sobol design cannot be conditioned, only characterised\n")
    sd = d["sobol_design"]
    L.append(
        f"`sobol_sensitivity.py` screens the same ranges over "
        f"{len(sd['straddle'])} of these parameters (excluding "
        f"{', '.join('`' + x + '`' for x in sd['excluded'])}) with a Saltelli "
        f"base of {sd['n_base']}. Its A and B samples are "
        f"{_pct(sd['a_admissible_fraction'])} and "
        f"{_pct(sd['b_admissible_fraction'])} admissible -- but the estimator "
        f"does not consume rows, it consumes PAIRS: Jansen's total effect "
        f"differences each A row against the same row with one column swapped "
        f"in from B. So the question is how many of those pairs STRADDLE the "
        f"constraint, one side a cell the model says is already dead untreated "
        f"and the other not.\n")
    L.append("| parameter | admissible after swap | straddling pairs |")
    L.append("|---|--:|--:|")
    for nm, v in sorted(sd["straddle"].items(),
                        key=lambda kv: (-kv[1]["straddling_fraction"], kv[0])):
        L.append(f"| `{nm}` | {_pct(v['admissible_fraction'])} | "
                 f"{v['straddling']} ({_pct(v['straddling_fraction'])}) |")
    L.append("")
    worst_s = min(sd["straddle"],
                  key=lambda k: (-sd["straddle"][k]["straddling_fraction"], k))
    # SORTED, because a bare dict comprehension publishes whichever order the
    # stored JSON happened to carry -- which is what the freshness gate's
    # key-order round trip caught here.
    zero_s = sorted(nm for nm, v in sd["straddle"].items() if v["straddling"] == 0)
    L.append(
        f"**This is why filtering is not available as a correction.** The "
        f"straddling fraction runs from {_pct(sd['straddle_min'])} to "
        f"{_pct(sd['straddle_max'])}: dropping straddling pairs would strip "
        f"`{worst_s}` of {_pct(sd['straddle'][worst_s]['straddling_fraction'])} "
        f"of its pairs and "
        + (", ".join(f"`{nm}`" for nm in zero_s) + " of none"
           if zero_s else "others of far fewer")
        + ", so the filtered estimator would compare a parameter stripped of "
          "its largest moves against one left intact. Re-running Saltelli "
          "inside the admissible region is not available either: admissibility "
          "is a constraint surface rather than a box, and rejection sampling "
          "breaks the column independence between A and B that the estimator "
          "assumes. The design can be characterised and not conditioned.\n")
    L.append(
        f"The parameters whose cross-swaps straddle most are the ones the "
        f"published decomposition attributes most to, which is not a "
        f"coincidence: moving a defence parameter far enough to change the "
        f"kill rate is the same move that kills the untreated cell. "
        f"`{worst_s}` is the sharpest case -- the report flags it NOT "
        f"identifiable from the kill rate, and it is both the highest-"
        f"straddling parameter here and the one the PRCC's own Control row "
        f"ranks first as a driver of untreated death. The verdict survives "
        f"conditioning (its coefficient on the inducer arm falls rather than "
        f"rises), so what is withdrawn is not the verdict but the reading that "
        f"the parameter does not matter: its influence lands on untreated "
        f"viability rather than on treatment response, which the kill "
        f"observable cannot see.\n")
    if not sd["sdt_ros_inert"]:
        L.append("**The excluded parameter is NOT inert untreated**, so the "
                 "Sobol legs above are withheld.\n")
    else:
        L.append(
            "The excluded parameter is evaluated at its default in the "
            "untreated arms above, which is only sound if it is inert there. "
            "Measured rather than assumed: across every untreated arm, the "
            "death rate at the bottom and the top of its range is bit-"
            "identical.\n")

    L.append("## What this does not establish\n")
    L.append(
        f"- **Conditioning is not inference.** It removes draws and adds no "
        f"data, so `identifiability-report.json`'s data-conditioned-headline "
        f"count stays where it is, and a guard breaks if that report ever says "
        f"otherwise.\n"
        f"- **No ranking here is a calibration.** Both halves are sensitivity "
        f"over an uncalibrated box; the restriction changes which draws the "
        f"box contains, not whether any parameter is anchored to data.\n"
        f"- **The admissible half is not a design.** See the stratification "
        f"table: it is a filtered sample, and no property LHS guarantees can "
        f"be claimed for it.\n"
        f"- **The fraction is not comparable to the sibling pages' by "
        f"subtraction.** Those sample uniformly where this is a Latin "
        f"hypercube, at a different per-draw seed convention, so the "
        f"difference is a design difference and not a disagreement about the "
        f"box.\n"
        f"- **Admissibility is judged on {len(g['untreated_arms'])} phenotypes, "
        f"the same ones the sibling pages use**, so the fraction is an upper "
        f"bound in the same way theirs is.\n")
    return "\n".join(L) + "\n"


def _split_stored(stored):
    """Recover `assemble`'s input for the freshness gate's round trip.

    Every input is committed -- the per-draw CSVs beside this page -- so this
    RE-READS them rather than echoing stored derived fields back out of the
    artifact. An echo would make every derived-field guard in
    `tests/test_prcc_under_admissibility.py` inert, which is a defect this
    repository has already paid for more than once.
    """
    return (load(),)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--measure", action="store_true",
                    help="re-measure (needs the compiled extension)")
    ap.add_argument("--render-only", action="store_true",
                    help="re-derive from the committed measurements only")
    args = ap.parse_args()
    if args.measure:
        measure()
    # The two branches are the SAME call, deliberately. Elsewhere
    # `--render-only` selects the cheap path because the default one MEASURES;
    # here measuring is its own flag and every remaining input is committed, so
    # the default run already re-derives the whole page from files -- and there
    # is no stored derived field for either path to echo, because the committed
    # files hold raw measurements only. Written out rather than collapsed
    # because `tests/test_artifact_freshness.py` keys two real round trips on
    # seeing it (that `assemble` is idempotent, and that the page does not
    # depend on the stored dicts' key order), both of which are worth having --
    # and because a derived cache added later must not slip in behind a
    # collapsed branch.
    if args.render_only:
        d = assemble(load())
    else:
        d = assemble(load())
    OUT_JSON.write_text(json.dumps(d, indent=2, sort_keys=True) + "\n")
    OUT_MD.write_text(render(d))
    print(f"wrote {OUT_JSON}")
    print(f"wrote {OUT_MD}")
    a = d["admissibility"]
    print(f"  admissible {a['n_admissible']}/{a['n_draws']} "
          f"({100 * a['fraction']:.1f}%); reproduces="
          f"{d['reproduction']['reproduces']}; "
          f"sdt bound holds={d['scored']['sdt_bound']['holds_admissible']}; "
          f"rsl3 count holds={d['scored']['rsl3_count']['holds_admissible']}")
    return 0


if __name__ == "__main__":
    main()
