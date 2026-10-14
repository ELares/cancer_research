#!/usr/bin/env python3
"""Condition the spatial and combination headline intervals on the model's own
baseline-viability constraint (#938).

`#936` conditioned the SINGLE-CELL prior-predictive intervals on the constraint
and found most of the published prior inadmissible. It deliberately left the
other headline families alone and said so on its own page, because leaving
them silently beside a corrected one would imply they had been checked. This is
that correction.

The families are `headline_uncertainty.py`'s `--headline
{bliss,sim-tme,penetration}`. They sample the SAME joint uniform prior over the
SAME PRCC rate constants, loaded from the same `analysis/prcc-results.json`
ranges the single-cell design uses -- verified here rather than assumed, so a
draw's admissibility is well defined and the two pages cannot disagree about what
admissible means.

WHAT MAKES THIS CHEAPER THAN IT LOOKS
-------------------------------------
Admissibility is a property of the PARAMETER VECTOR, judged by the untreated
arms through the Python extension in milliseconds. The expensive binaries do not
enter it at all: they supply the observable, admissibility supplies the mask, and
the two are measured on the same draws. The sim-tme family draws a PREFIX of what the
other families draw (a property of `sample_prior`, guarded), so ONE
admissibility pass covers every family and the partitions are nested rather
than independent samples.

ONE SAMPLE PARTITIONED, NOT TWO SAMPLED. Every draw is run and then split; a
second sample would confound the conditioning with sampling noise.

THE REPRODUCTION GATE
---------------------
The full-prior half must reproduce each committed report's own published figures
within a tolerance DERIVED from that report's printing precision. A family that
does not reproduce has its comparison WITHHELD, because a page that compares
against a re-run is not correcting what is published.

OFFLINE CONTRACT
----------------
The measurement needs the compiled binaries and the Python extension, so it is
taken locally with `--measure` and the per-draw results committed. Every figure
is then derived from that one file, so CI re-derives the whole page and runs
nothing.

NOTE FOR WHOEVER MEASURES NEXT: `headline_sensitivity._default_binary` looks only
in `simulations/target/{release,debug}`, so a machine whose cargo config
redirects `target-dir` has no discoverable binary and `--binary-dir` must be
passed.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import os
import re
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts"))

# ONE source for the published headline design: draws, counts, seed and the
# committed report paths. Restating any of them is how a reproduction quietly
# stops reproducing anything.
from headline_uncertainty import (  # noqa: E402
    DEFAULT_PEN_SAMPLES,
    DEFAULT_SAMPLES,
    DEFAULT_TME_SAMPLES,
    PEN_REPORT,
    REPORT as BLISS_REPORT,
    SEED,
    TME_REPORT,
    sample_prior,
)
# ONE source for the observables, so this page cannot measure a different
# quantity from the one the published intervals describe.
from headline_sensitivity import (  # noqa: E402
    PARAM_NAMES,
    PENETRATION_TISSUES,
    run_bliss_row,
    run_sim_tme_observables,
    run_tissue_pk_observables,
)
# ONE source for what "admissible" means, and for the circularity falsifier.
from intervals_under_admissibility import (  # noqa: E402
    CONTROL_CONDITIONS,
    MIN_ADMISSIBLE_DRAWS,
    STRONG_SELECTION,
    _label,
    _selection,
)
# ONE source for the constraint, whose own guard parses it from the manuscript.
from headline_at_fitted import BASELINE_MAX  # noqa: E402
# ONE source for the per-draw cell count and simulation seed.
from uncertainty_intervals import N_CELLS, QUANTILES, SIM_SEED  # noqa: E402

OUT_MD = REPO / "analysis" / "headline-intervals-under-admissibility.md"
OUT_JSON = REPO / "analysis" / "headline-intervals-under-admissibility.json"
DRAWS_CSV_GZ = REPO / "analysis" / "calibration" / "headline-admissibility-draws.csv.gz"

# The widest published draw count; the other families are prefixes of it.
N_ROWS = max(DEFAULT_SAMPLES, DEFAULT_TME_SAMPLES, DEFAULT_PEN_SAMPLES)

CONTROL_COLUMNS = tuple(_label(c) for c in CONTROL_CONDITIONS)
BLISS_COLUMN = "bliss_synergy"
# `synergy_score` is a RATIO, so a change in it has two possible sources and the
# page must be able to say which. These are its parts, straight from the row the
# score comes from: the two single-agent rates, the combination rate (numerator)
# and the Bliss independence prediction (denominator).
BLISS_PART_COLUMNS = ("bliss_rate_rsl3", "bliss_rate_fsp1i",
                      "bliss_rate_combo", "bliss_prediction")
TME_COLUMNS = ("tme_hypoxia_gap", "tme_immune")
PEN_COLUMNS = tuple(f"pen_{key}" for key, _ in PENETRATION_TISSUES)

# A leg is (columns it fills, how many draws the published family ran, binary).
LEGS = {
    "admissibility": (CONTROL_COLUMNS, N_ROWS, None),
    "bliss": ((BLISS_COLUMN,) + BLISS_PART_COLUMNS, DEFAULT_SAMPLES,
               "sim-combo-mech"),
    "tme": (TME_COLUMNS, DEFAULT_TME_SAMPLES, "sim-tme"),
    "penetration": (PEN_COLUMNS, DEFAULT_PEN_SAMPLES, "sim-tissue-pk"),
}

COLUMNS = ["draw"] + list(PARAM_NAMES) + list(CONTROL_COLUMNS) \
    + [BLISS_COLUMN] + list(BLISS_PART_COLUMNS) \
    + list(TME_COLUMNS) + list(PEN_COLUMNS)

# How often a long leg flushes its partial results to disk.
CHECKPOINT_EVERY = 10

# A subprocess that errored, kept distinct from "not measured" (empty) and from
# an undefined-but-successful observable ("nan", which the published generator
# counts separately rather than dropping silently).
FAILED = "fail"

# Published reports print percentages to one decimal and the Bliss ratio to
# three, so parsing a figure back loses at most half of its last digit. The
# tolerance is DERIVED from that rather than picked: #936's first attempt used a
# hand-picked round number that was 22x looser than the largest real deviation,
# which would accept a genuine discrepancy while reporting a reproduction.
ROUNDING_SLACK = 2.0


def _tolerance(decimals: int, scale: float = 1.0) -> float:
    """Half of the last printed digit, times the slack, on the value's own scale."""
    return ROUNDING_SLACK * 0.5 * 10 ** (-decimals) * scale


# A quantile moving by less than this FRACTION OF THE PUBLISHED INTERVAL WIDTH is
# not a correction worth reporting as one. Expressed as a fraction because the
# three families are on three different scales (a synergy ratio, a signed kill
# gap, kill fractions) and one absolute number there would be three thresholds
# wearing one.
MATERIAL_FRACTION = 0.10

# Sites that draw a conclusion from the SUB-ADDITIVE draws. Every one of those
# draws turns out to be inadmissible, so the conclusion rests on a corner the
# model rejects -- and these are handwritten interpretation rather than
# generated prose, so this page NAMES them instead of editing them. Checked, so
# the note cannot outlive the sentences it points at.
# ONE ENTRY PER SENTENCE, not per file. This registry's job is stated in the
# paragraph it renders -- "the naming is checked against their current text so it
# cannot outlive the sentences it points at" -- and it points at SENTENCES while
# it was keyed on FILES, one representative phrase each. Measured (#942): the
# manuscript carries FIVE sentences resting on this corner and the registry
# verified one of them, the unverified set including `v1.md`'s opening summary,
# which is the most prominent statement of the claim anywhere here.
#
# (#938) named two files and (#940) found a third by a repo-wide search, so both
# of those rounds were about WHICH FILES and neither noticed that the question
# was the wrong one. The reason the choice went unnoticed is mechanical and is
# fixed below: `load` keyed the result dict on the PATH, so several entries for
# one file collapsed into one and a per-sentence registry could not have worked.
#
# EVERY sentence carrying the concept is registered, with no classifier deciding
# which of them "draws a conclusion" -- a keyword rule cannot measure that, and
# this repository has that measurement on record. Two of the seven below report
# the range rather than concluding from it; that reading is a JUDGEMENT stated
# here and is deliberately not computed.
SUBADDITIVE_CITATIONS = (
    ("README.md",
     "some draws are sub-additive, so synergy is not guaranteed"),
    ("analysis/p1-wetlab-protocol.md",
     "some draws are sub-additive, so supra-additivity is not guaranteed"),
    # the book's opening summary
    ("article/drafts/v1.md",
     "includes additivity and some sub-additive draws"),
    # Chapter 5's restatement of the identifiability verdict
    ("article/drafts/v1.md",
     "includes sub-additive draws, so its supra-additive direction"),
    # Section 7.3's own statement of what the interval does not establish
    ("article/drafts/v1.md",
     "the full range includes sub-additive draws"),
    # reports the range, with its minimum (judged descriptive)
    ("article/drafts/v1.md",
     "and some sub-additive draws (minimum"),
    ("article/drafts/v1.md",
     "The lower bound rounds to additivity and some draws are sub-additive"),
)

# Spelled out so the prose naming these sites cannot carry a stale count. The
# first version wrote "Two sites" by hand beside a registry of two, which is the
# derived-number-handwritten-prose defect. (#940) then derived the count and
# derived it from the WRONG QUANTITY: it counted the FILES and the sentence it
# renders says sites.
_COUNT_WORDS = {1: "One", 2: "Two", 3: "Three", 4: "Four", 5: "Five",
                6: "Six", 7: "Seven", 8: "Eight"}


def _counted(n: int, noun: str = "sentence", cap: bool = False) -> str:
    """"<word> <noun>(s)", with the plural computed.

    The number word and the noun are kept APART deliberately: a literal
    pairing a spelled count with the noun is what
    `tests/test_no_handwritten_figures_in_generators.py` flags, and it is right
    to -- such a phrase in a generator is a sentence that can outlive the
    registry it describes, which is exactly how the first version shipped.
    """
    word = _COUNT_WORDS.get(n, str(n))
    if not cap:
        word = word.lower() if word[0].isalpha() else word
    return f"{word} {noun}{'' if n == 1 else 's'}"

# A Bliss ratio equals this when the combination kill exactly matches the
# independence prediction, because the engine computes `rate_combo /
# bliss_prediction`. Definitional for a ratio rather than a measured figure --
# and it is the line the published report's own direction claim rests on.
BLISS_NO_SYNERGY = 1.0

# Registered in the issue before any observable existed: #926 established the
# penetration ordering is robust BY CONSTRUCTION (one shared monotone
# dose-response read at three decreasing fixed concentrations), so admissibility
# cannot break it. Stating the expectation here is what stops "the ordering
# survived conditioning" being reported as new evidence.
ORDERING_IS_STRUCTURAL = True


def draw_matrix() -> np.ndarray:
    """The published families' own draws, regenerated from their shared seed."""
    return sample_prior(N_ROWS, seed=SEED)


def _fmt(v) -> str:
    if v is None:
        return FAILED
    if isinstance(v, float) and not math.isfinite(v):
        return "nan"
    return f"{v:.6f}"


def _blank_table(draws) -> list:
    rows = []
    for i, row in enumerate(draws):
        rec = {c: "" for c in COLUMNS}
        rec["draw"] = str(i)
        for j, name in enumerate(PARAM_NAMES):
            rec[name] = f"{row[j]:.10g}"
        rows.append(rec)
    return rows


def _read_table(draws) -> list:
    """The committed per-draw file if it exists, else a blank table. Measuring a
    leg must not discard another leg's hours of binary runs."""
    if not DRAWS_CSV_GZ.exists():
        return _blank_table(draws)
    with gzip.open(DRAWS_CSV_GZ, "rt") as fh:
        rows = [dict(r) for r in csv.DictReader(fh)]
    for rec in rows:                       # tolerate a file written before a
        for c in COLUMNS:                  # column existed
            rec.setdefault(c, "")
    if len(rows) != len(draws):
        raise SystemExit(
            f"{DRAWS_CSV_GZ.name} holds {len(rows)} rows against {len(draws)} draws; "
            "the design changed, so re-measure every leg rather than merging")
    return rows


def _write_table(rows) -> None:
    DRAWS_CSV_GZ.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(DRAWS_CSV_GZ, "wt", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=COLUMNS)
        writer.writeheader()
        writer.writerows(rows)


def _binary(binary_dir, name):
    cand = Path(binary_dir) / name
    if not cand.exists():
        raise SystemExit(f"no {name} binary at {cand}; build it and pass --binary-dir")
    return cand


def measure(leg, binary_dir, workers) -> None:
    """Fill one leg's columns on the shared draws, leaving every other leg alone."""
    draws = draw_matrix()
    rows = _read_table(draws)
    columns, n, binary_name = LEGS[leg]

    if leg == "admissibility":
        import uncertainty_intervals as U
        for i in range(n):
            for (ph, tr), col in zip(CONTROL_CONDITIONS, columns):
                overrides = {name: float(draws[i][j])
                             for j, name in enumerate(PARAM_NAMES)}
                rate = U._fc().sim_batch(ph, tr, n=N_CELLS, seed=SIM_SEED,
                                         **overrides)["death_rate"]
                rows[i][col] = _fmt(rate)
            if (i + 1) % 50 == 0:
                print(f"  admissibility {i + 1} of {n}", flush=True)
        _write_table(rows)
        print(f"wrote {DRAWS_CSV_GZ.relative_to(REPO)} ({leg})")
        return

    binary = _binary(binary_dir, binary_name)
    runner = {"bliss": _bliss_observables, "tme": run_sim_tme_observables,
              "penetration": run_tissue_pk_observables}[leg]

    done = [0]

    def _one(i):
        try:
            out = runner(draws[i], binary)
        except Exception as exc:                      # subprocess / parse failure
            print(f"  draw {i} failed: {exc}", flush=True)
            out = None
        done[0] += 1
        if done[0] % 10 == 0:
            print(f"  {leg} {done[0]} of {n}", flush=True)
        return i, out

    # sim-tme costs minutes per draw, so a leg is an hour of work that a dropped
    # connection or a SIGTERM would otherwise discard whole. The table is 300
    # rows, so checkpointing it is free beside one run of the binary.
    def _checkpoint(k):
        if k % CHECKPOINT_EVERY == 0:
            _write_table(rows)

    # Which key of the runner's dict lands in which column. Every runner here
    # returns a dict keyed to match its leg's columns in order.
    keys = {"bliss": (BLISS_COLUMN,) + BLISS_PART_COLUMNS,
            "tme": ("hypoxia", "immune"),
            "penetration": tuple(key for key, _ in PENETRATION_TISSUES)}[leg]

    with ThreadPoolExecutor(max_workers=workers) as ex:  # subprocesses release the GIL
        for k, (i, out) in enumerate(ex.map(_one, range(n)), start=1):
            if out is None:
                for col in columns:
                    rows[i][col] = FAILED
            else:
                for col, key in zip(columns, keys):
                    rows[i][col] = _fmt(out[key])
            _checkpoint(k)
    _write_table(rows)
    print(f"wrote {DRAWS_CSV_GZ.relative_to(REPO)} ({leg})")


# The published generator's own tissue labels, imported so the parser keys on
# what that generator PRINTS and cannot drift if a tissue is renamed.
from headline_uncertainty import _TISSUE_LABELS, _pctiles  # noqa: E402

# Published figures are printed to this many decimals in all three reports, so
# parsing one back loses at most half of its last digit.
PUBLISHED_DECIMALS = 3

# The sibling page this one extends, read only to check the two agree about the
# prior. A disagreement means the two definitions of admissible have diverged.
SINGLE_CELL_JSON = REPO / "analysis" / "intervals-under-admissibility.json"

# The standing refusal: this page must not be able to move that report's
# data-conditioned-headline count, and breaks if that count ever changes.
IDENTIFIABILITY_JSON = REPO / "analysis" / "identifiability-report.json"

# Where the engine decides a Bliss ratio is undefined.
COMBO_MECH_RS = REPO / "simulations" / "sim-combo-mech" / "src" / "main.rs"

# Beyond this many standard errors the two pages are not describing one prior.
AGREEMENT_Z = 3.0

FAMILY_SAMPLES = {"bliss": DEFAULT_SAMPLES,
                  "tme": DEFAULT_TME_SAMPLES,
                  "penetration": DEFAULT_PEN_SAMPLES}

# (column, family, label as this page prints it, key in the parsed report)
#
# The labels carry NO COMMAS, deliberately: several sections list the affected
# observables inline, and a label with a comma in it makes such a list
# unreadable -- there is no telling where one entry ends.
OBSERVABLES = (
    (BLISS_COLUMN, "bliss", "Bliss synergy (RSL3 + FSP1i)", "synergy"),
    ("tme_hypoxia_gap", "tme", "Hypoxia kill-collapse gap (SDT minus RSL3)",
     "hypoxia"),
    ("tme_immune", "tme", "SDT immune kill fraction (final counts)", "immune"),
    ("pen_well_vascularized", "penetration",
     "Vessel-wall RSL3 kill: well-vascularized", "well_vascularized"),
    ("pen_poorly_vascularized", "penetration",
     "Vessel-wall RSL3 kill: poorly-vascularized", "poorly_vascularized"),
    ("pen_cns_bbb", "penetration", "Vessel-wall RSL3 kill: CNS / BBB", "cns_bbb"),
)


def _bliss_from_row(row) -> dict:
    """One sim-combo-mech RSL3 + FSP1i row, keyed to this leg's columns.

    The single-agent rates are resolved BY DRUG NAME. The shared reader matches
    `{drug_a, drug_b}` as a SET, so which column holds RSL3 is not fixed by the
    format, and taking `rate_a` positionally would silently swap the two agents
    on any run where the binary emits them the other way round. Pure, so that
    property is unit-tested without the binary."""
    rate = {row["drug_a"]: float(row["rate_a"]),
            row["drug_b"]: float(row["rate_b"])}
    return {BLISS_COLUMN: float(row["synergy_score"]),
            "bliss_rate_rsl3": rate["RSL3"],
            "bliss_rate_fsp1i": rate["FSP1i"],
            "bliss_rate_combo": float(row["rate_combo"]),
            "bliss_prediction": float(row["bliss_prediction"])}


def _bliss_observables(params_row, binary):
    """The Bliss score AND the parts it is a ratio of, from one run."""
    return _bliss_from_row(run_bliss_row(params_row, binary))


def undefined_bliss_floor(rust_src) -> float:
    """The denominator below which the engine refuses the ratio, PARSED from the
    binary's own source rather than typed. It is a property of the engine, and a
    constant here would be a second copy free to disagree with it."""
    m = re.search(r"if\s+bliss\s*>\s*([0-9.]+)", rust_src)
    if not m:
        raise SystemExit("cannot find the undefined-Bliss floor in sim-combo-mech")
    return float(m.group(1))


def _cells(line):
    if not line.strip().startswith("|"):
        return None
    return [c.strip() for c in line.strip().strip("|").split("|")]


def _numbers(text):
    return [float(x) for x in re.findall(r"[-+]?\d*\.\d+|\d+", text)]


def _quantity_table(section) -> dict:
    """The `| quantity | value |` shape the Bliss and sim-tme reports both use."""
    out = {}
    for line in section.split("\n"):
        cells = _cells(line)
        if not cells or len(cells) != 2:
            continue
        label, value = cells[0].lower(), cells[1]
        nums = _numbers(value)
        if label.startswith("default point estimate") and len(nums) == 1:
            out["point"] = nums[0]
        elif label.startswith("prior-predictive median") and len(nums) == 1:
            out["median"] = nums[0]
        elif "prior-predictive interval" in label and len(nums) == 2:
            out["low"], out["high"] = nums
        elif label.startswith("full range") and len(nums) == 2:
            out["min"], out["max"] = nums
    return out


def parse_published(texts) -> dict:
    """Each committed report's OWN published figures, so the comparison is with
    what is published rather than with a re-run of it."""
    out = {"synergy": _quantity_table(texts["bliss"])}

    for section in texts["tme"].split("\n## "):
        if section.lower().startswith("hypoxia kill-collapse"):
            out["hypoxia"] = _quantity_table(section)
        elif section.lower().startswith("immune kill fraction"):
            out["immune"] = _quantity_table(section)

    label_to_key = {v: k for k, v in _TISSUE_LABELS.items()}
    for line in texts["penetration"].split("\n"):
        cells = _cells(line)
        if not cells or len(cells) != 5 or cells[0] not in label_to_key:
            continue
        interval, full = _numbers(cells[3]), _numbers(cells[4])
        if len(interval) != 2 or len(full) != 2:
            continue
        out[label_to_key[cells[0]]] = {
            "point": _numbers(cells[1])[0], "median": _numbers(cells[2])[0],
            "low": interval[0], "high": interval[1],
            "min": full[0], "max": full[1]}

    ordering = re.search(r"ordering preserved[^:]*:\s*\*\*([\d.]+)%\*\*\s*of\s*(\d+)",
                         texts["penetration"])
    if ordering:
        out["ordering"] = {"fraction": float(ordering.group(1)) / 100.0,
                           "n": int(ordering.group(2))}
    return out


def _col(rows, column, n):
    """(index, value) for every draw of a family whose column holds a number,
    keeping the index so the admissibility mask still lines up after drops."""
    got = []
    for i in range(n):
        raw = rows[i].get(column, "")
        if raw in ("", FAILED):
            continue
        got.append((i, float(raw)))
    return got


def _refuse_incomplete(rows) -> None:
    """A leg measured under an older column set, or half-measured, would publish
    an interval over whatever happened to be present. Refuse instead."""
    missing = []
    for leg, (columns, n, _) in LEGS.items():
        for col in columns:
            absent = sum(1 for i in range(n) if rows[i].get(col, "") == "")
            if absent:
                missing.append(f"{leg}:{col} ({absent} of {n} draws unmeasured)")
    if missing:
        raise SystemExit("the per-draw file is incomplete, so no interval can be "
                         "published from it:\n  " + "\n  ".join(missing)
                         + f"\n\nRe-measure with --measure, then regenerate.")


def _identity_error(rows, n) -> float:
    """How far the stored Bliss parts stray from the engine's own
    `rate_a + rate_b - rate_a*rate_b`. A column mix-up is otherwise invisible,
    because a swapped pair still looks like a plausible ratio. Rows whose run
    failed are skipped rather than crashing the page."""
    worst = 0.0
    for i in range(n):
        vals = {}
        for col in ("bliss_prediction", "bliss_rate_rsl3", "bliss_rate_fsp1i"):
            raw = rows[i].get(col, "")
            if raw in ("", FAILED):
                break
            vals[col] = float(raw)
        if len(vals) != 3:
            continue
        a, b = vals["bliss_rate_rsl3"], vals["bliss_rate_fsp1i"]
        worst = max(worst, abs(vals["bliss_prediction"] - (a + b - a * b)))
    return worst


def assemble(raw) -> dict:
    rows = raw["rows"]
    if len(rows) < N_ROWS:
        raise SystemExit(
            f"the per-draw file holds {len(rows)} rows against {N_ROWS} draws in "
            "the widest published family, so a family's partition would run off "
            "the end of it")
    _refuse_incomplete(rows)
    published = parse_published(raw["published"])
    floor = undefined_bliss_floor(raw["rust"])
    tol = _tolerance(PUBLISHED_DECIMALS)

    worst = [max(float(r[c]) for c in CONTROL_COLUMNS) for r in rows]
    admissible = [w < BASELINE_MAX for w in worst]

    # The prior is one box, so the three families' fractions are three samples
    # of one quantity and the sibling page measured it on a fourth. Checked
    # rather than asserted: a real disagreement would mean the two pages no
    # longer judge admissibility the same way.
    sibling = raw["sibling"]["admissibility"]
    families = {}
    for fam, n in FAMILY_SAMPLES.items():
        k = sum(admissible[:n])
        p = sibling["fraction"]
        se = math.sqrt(p * (1.0 - p) / n)
        families[fam] = {
            "draws": n, "n_admissible": k, "fraction": k / n,
            "usable": k >= MIN_ADMISSIBLE_DRAWS,
            "z_against_sibling": (k / n - p) / se if se > 0 else None,
        }
    agreement = {
        "sibling_fraction": sibling["fraction"],
        "sibling_draws": sibling["n_draws"],
        "max_abs_z": max(abs(f["z_against_sibling"]) for f in families.values()),
        "bar": AGREEMENT_Z,
    }
    agreement["agrees"] = agreement["max_abs_z"] < AGREEMENT_Z

    observables = []
    for column, fam, label, key in OBSERVABLES:
        n = FAMILY_SAMPLES[fam]
        got = _col(rows, column, n)
        failed = sum(1 for i in range(n) if rows[i].get(column, "") == FAILED)
        undefined = [i for i, v in got if not math.isfinite(v)]
        full = [v for _, v in got if math.isfinite(v)]
        adm = [v for i, v in got if math.isfinite(v) and admissible[i]]
        pub = published.get(key)

        f = _pctiles(np.asarray(full, float)) if full else None
        a = _pctiles(np.asarray(adm, float)) if adm else None

        # The full-prior half must reproduce the published figures, or the page
        # is correcting a re-run rather than what is published.
        checks, repro = {}, None
        if pub and f:
            pairs = (("median", "median"), ("low", "p2_5"), ("high", "p97_5"),
                     ("min", "min"), ("max", "max"))
            for pk, fk in pairs:
                if pk in pub:
                    checks[pk] = {"published": pub[pk], "recomputed": f[fk],
                                  "deviation": abs(f[fk] - pub[pk]),
                                  "reproduces": abs(f[fk] - pub[pk]) <= tol}
            repro = all(c["reproduces"] for c in checks.values()) and len(checks) > 0

        width = (pub["high"] - pub["low"]) if pub else None
        # A conditioned TAIL needs enough admissible draws to rest on; a
        # conditioned MEDIAN needs far fewer. The floor is the sibling page's,
        # imported -- and the costliest family's published design draws no more
        # than the floor in total, so it can only meet the floor if the
        # constraint removes nothing at all from it. That is reported, not
        # worked around: the tails are withheld and the median is quoted.
        quotable_tail = families[fam]["usable"]
        shifts = {}
        if pub and f and a:
            for pk, fk in (("median", "median"), ("low", "p2_5"), ("high", "p97_5")):
                d = a[fk] - f[fk]
                quotable = True if pk == "median" else quotable_tail
                shifts[pk] = {"full": f[fk], "admissible": a[fk], "delta": d,
                              "quotable": quotable,
                              "material": ((abs(d) > MATERIAL_FRACTION * width)
                                           if width and width > 0 else None)}

        # Where the manuscript's own unperturbed point estimate falls once the
        # invalid draws are gone. #936's sharpest single number was a median
        # crossing from one side of it to the other.
        point = None
        if pub and "point" in pub and f and a:
            point = {
                "value": pub["point"],
                "full_median_above_point": f["median"] > pub["point"],
                "admissible_median_above_point": a["median"] > pub["point"],
                "crossed": (f["median"] > pub["point"]) != (a["median"] > pub["point"]),
                # Only where the tails are quotable. Placing the point inside an
                # interval this page has just refused to publish would be using
                # a bound it says it does not have.
                "inside_admissible_interval": (
                    (a["p2_5"] <= pub["point"] <= a["p97_5"])
                    if quotable_tail else None),
            }

        observables.append({
            "column": column, "family": fam, "label": label, "published_key": key,
            "n_full": len(full), "n_admissible": len(adm),
            "n_failed": failed, "n_undefined": len(undefined),
            "n_undefined_admissible": sum(1 for i in undefined if admissible[i]),
            "published": pub, "published_width": width,
            "full": f, "admissible_only": a,
            "reproduction": checks, "reproduces": repro,
            "shifts": shifts, "point": point,
            "interval_quotable": quotable_tail,
            "n_admissible_floor": MIN_ADMISSIBLE_DRAWS,
            # The identifiability report's per-headline verdicts are DIRECTION
            # claims, so what would move them is a sign change, not a magnitude
            # change. Measured per observable rather than argued.
            "sign_preserved": (None if not (f and a) else
                               (f["median"] >= 0) == (a["median"] >= 0)),
            # Something can be said about this observable: it reproduces what is
            # published and has admissible draws to condition on. Whether its
            # TAILS can be quoted is `interval_quotable`, separately.
            "reportable": bool(repro) and a is not None,
        })

    # What the truncation selects, and the circularity falsifier with it.
    selection = _selection(rows, admissible)

    # The Bliss score is a ratio, so its move has two possible sources.
    parts = {}
    for col in BLISS_PART_COLUMNS:
        got = _col(rows, col, FAMILY_SAMPLES["bliss"])
        full = [v for _, v in got if math.isfinite(v)]
        adm = [v for i, v in got if math.isfinite(v) and admissible[i]]
        parts[col] = {"full": _pctiles(np.asarray(full, float)) if full else None,
                      "admissible_only": _pctiles(np.asarray(adm, float)) if adm else None}
    den = parts["bliss_prediction"]
    num = parts["bliss_rate_combo"]
    den_adm, num_adm = den["admissible_only"], num["admissible_only"]
    decomposition = {
        "parts": parts,
        "floor": floor,
        "admissible_min_denominator": den_adm["min"] if den_adm else None,
        "denominator_clears_floor": (den_adm["min"] > floor) if den_adm else None,
        # Both move; the ratio follows whichever moves more. Stating which is
        # the difference between a mechanism and an assertion.
        "numerator_ratio": (num_adm["median"] / num["full"]["median"]
                            if num_adm and num["full"]["median"] else None),
        "denominator_ratio": (den_adm["median"] / den["full"]["median"]
                              if den_adm and den["full"]["median"] else None),
        # The identity the engine computes, checked on the stored parts: a
        # column mix-up would otherwise be invisible.
        "identity_max_error": _identity_error(rows, FAMILY_SAMPLES["bliss"]),
    }
    # The published report's conclusion is about the lower BOUND, not the
    # median: it reads the supra-additive direction as robust because the
    # interval's lower end does not go below the no-synergy line. Which draws
    # sit at or below that line, and whether the model accepts any of them, is
    # therefore the thing to measure rather than describe.
    bliss_vals = [(i, v) for i, v in _col(rows, BLISS_COLUMN,
                                          FAMILY_SAMPLES["bliss"])
                  if math.isfinite(v)]
    at_or_below = [i for i, v in bliss_vals if v <= BLISS_NO_SYNERGY]
    on_line = [i for i, v in bliss_vals if v == BLISS_NO_SYNERGY]
    saturated = [float(rows[i]["bliss_rate_rsl3"]) for i in on_line
                 if rows[i]["bliss_rate_rsl3"] not in ("", FAILED)]
    decomposition["no_synergy_line"] = BLISS_NO_SYNERGY
    decomposition["n_at_or_below_line"] = len(at_or_below)
    decomposition["n_at_or_below_line_admissible"] = sum(
        1 for i in at_or_below if admissible[i])
    decomposition["n_on_line"] = len(on_line)
    # Reported as the weakest single-arm kill among the on-line draws rather
    # than classified against a saturation threshold nobody measured.
    decomposition["min_single_arm_kill_on_line"] = min(saturated) if saturated else None

    # EXPLAINING TWO OF SIX READS AS EXPLAINING ALL SIX. The saturation note
    # above covers only the draws sitting exactly on the line, so the corner's
    # own mechanism is measured across EVERY draw in it: this ratio can reach
    # the line only when the independence prediction is already high, which
    # needs both single arms killing hard, which in this model is the same
    # strong-peroxidation regime that kills the untreated control. Reported as
    # MINIMA over the corner, so each figure is a floor every corner draw
    # clears rather than a mean one draw could carry.
    def _corner_min(cols, combine):
        vals = []
        for i in at_or_below:
            try:
                vals.append(combine([float(rows[i][c]) for c in cols]))
            except (ValueError, KeyError):
                continue
        return min(vals) if vals else None

    decomposition["corner_min_worst_untreated"] = _corner_min(CONTROL_COLUMNS, max)
    decomposition["corner_min_larger_single_arm"] = _corner_min(
        ("bliss_rate_rsl3", "bliss_rate_fsp1i"), max)
    decomposition["corner_min_independence_prediction"] = _corner_min(
        ("bliss_prediction",), max)
    decomposition["sub_additive_corner_is_wholly_inadmissible"] = (
        len(at_or_below) > 0
        and decomposition["n_at_or_below_line_admissible"] == 0)

    decomposition["ratio_follows_denominator"] = (
        decomposition["denominator_ratio"] is not None
        and decomposition["numerator_ratio"] is not None
        and decomposition["denominator_ratio"] < decomposition["numerator_ratio"])

    # Registered in the issue BEFORE any observable existed: by #926 this
    # ordering is robust by construction, so it cannot break here and
    # "it survived conditioning" is not evidence. It is a harness check.
    def _ordering(idx):
        kept = [i for i in idx
                if all(rows[i].get(c, "") not in ("", FAILED) for c in PEN_COLUMNS)]
        ok = sum(1 for i in kept
                 if float(rows[i]["pen_well_vascularized"])
                 >= float(rows[i]["pen_poorly_vascularized"])
                 >= float(rows[i]["pen_cns_bbb"]))
        return {"n": len(kept), "preserved": ok,
                "fraction": ok / len(kept) if kept else None}

    n_pen = FAMILY_SAMPLES["penetration"]
    ordering = {
        "full": _ordering(range(n_pen)),
        "admissible_only": _ordering([i for i in range(n_pen) if admissible[i]]),
        "structural": ORDERING_IS_STRUCTURAL,
        "published": published.get("ordering"),
    }
    ordering["reproduces_published"] = (
        ordering["published"] is not None
        and ordering["full"]["fraction"] is not None
        and abs(ordering["full"]["fraction"] - ordering["published"]["fraction"]) <= tol)

    # The restraint, guarded: conditioning is not inference.
    data_conditioned = raw["identifiability"]["data_constrained_in_production"]

    reportable = [o for o in observables if o["reportable"]]
    withheld_tails = [o["label"] for o in observables if not o["interval_quotable"]]
    return {
        "constraint": {"baseline_max": BASELINE_MAX, "n_cells": N_CELLS,
                       "sim_seed": SIM_SEED, "seed": SEED,
                       "controls": list(CONTROL_COLUMNS)},
        "worst_untreated": {"median": float(np.median(worst)),
                            "max": float(max(worst)),
                            "n_over_half": int(sum(1 for w in worst if w > 0.5))},
        "families": families,
        "agreement": agreement,
        "observables": observables,
        "selection": selection,
        "bliss_decomposition": decomposition,
        "ordering": ordering,
        "tolerance": tol,
        "material_fraction": MATERIAL_FRACTION,
        "data_conditioned_headlines": data_conditioned,
        "n_reportable": len(reportable),
        "withheld_tails": withheld_tails,
        "subadditive_citations": raw["subadditive_citations"],
        "sign_flips": [o["label"] for o in observables
                       if o["sign_preserved"] is False],
        "n_material": sum(1 for o in reportable
                          if any(s["material"] for s in o["shifts"].values()
                                 if s["quotable"])),
        "all_reproduce": all(o["reproduces"] for o in observables),
        "verdict": (len(reportable) == len(observables)
                    and all(o["reproduces"] for o in observables)
                    and selection["drug_parameters_untouched"] is True
                    and agreement["agrees"]
                    and data_conditioned == 0),
    }


def _flat(text: str) -> str:
    """Collapse runs of whitespace, so a wrapped sentence still matches."""
    return re.sub(r"\s+", " ", text)


def _phrase_present(path, phrase) -> bool:
    """Does `path` still carry `phrase`?

    Its own function so a guard can drive BOTH answers. Inlined in `load`, the
    read is unreachable from any test that doctors the loaded inputs, and
    hard-coding it to True survived a mutation sweep for exactly that reason.
    """
    # WHITESPACE-INSENSITIVE, because the documents are wrapped and a sentence
    # does not stop being present when a line break falls inside it. The raw
    # comparison is why (#940)'s fragments had to stop mid-clause, which is how a
    # registry ends up keyed on a phrase short enough to be ambiguous.
    return _flat(phrase) in _flat((REPO / path).read_text())


def load() -> dict:
    """Read the committed per-draw results and the published reports.
    No extension, no binary, no network."""
    with gzip.open(DRAWS_CSV_GZ, "rt") as fh:
        rows = list(csv.DictReader(fh))
    return {
        "rows": rows,
        "rust": COMBO_MECH_RS.read_text(),
        # Both of these gate the verdict, so they belong in the INPUTS rather
        # than being read inside `assemble`: a check whose input a guard cannot
        # doctor is a check nothing can drive.
        # KEYED ON THE SITE, not the path. Keyed on the path, several entries
        # for one file collapsed into one -- which is why the registry was
        # file-granular without the choice ever being made (#942).
        "subadditive_citations": {f"{path} :: {phrase}":
                                  _phrase_present(path, phrase)
                                  for path, phrase in SUBADDITIVE_CITATIONS},
        "sibling": json.loads(SINGLE_CELL_JSON.read_text()),
        "identifiability": json.loads(IDENTIFIABILITY_JSON.read_text()),
        "published": {
            "bliss": BLISS_REPORT.read_text(),
            "tme": TME_REPORT.read_text(),
            "penetration": PEN_REPORT.read_text(),
        },
    }


def _pct(x, places=1):
    return "n/a" if x is None else f"{100.0 * x:.{places}f}%"


def _num(x, places=3):
    return "n/a" if x is None else f"{x:.{places}f}"


def _signed(x, places=3):
    return "n/a" if x is None else f"{x:+.{places}f}"


def _yesno(v):
    return "n/a" if v is None else ("yes" if v else "no")


def _confidence_label():
    """The interval's own confidence, derived from the published quantiles rather
    than typed: a figure in prose here could outlive the design it describes."""
    low, _, high = QUANTILES
    return f"{(high - low) * 100.0:.0f}%"


def render(d) -> str:
    L = []
    cons, fams = d["constraint"], d["families"]
    sel, dec, order = d["selection"], d["bliss_decomposition"], d["ordering"]
    conf = _confidence_label()

    L.append("# The spatial and combination headline intervals, conditioned on "
             "the model's own viability constraint\n")
    L.append("Generated by `scripts/headline_intervals_under_admissibility.py`. "
             "`scripts/headline_uncertainty.py` publishes prior-predictive "
             f"{conf} intervals for several headline families -- the RSL3 + FSP1i "
             "Bliss synergy, the two spatial `sim-tme` headlines, and the "
             "per-tissue vessel-wall kill -- by sweeping a joint uniform prior "
             "over the same PRCC rate constants the single-cell design uses. "
             "`intervals-under-admissibility.md` (#936) conditioned the "
             "single-cell intervals on the manuscript's own baseline-viability "
             "constraint and left these explicitly outstanding. This "
             "applies the same conditioning to them.\n")

    L.append("## The constraint, and what it leaves of each family's prior\n")
    L.append(f"A draw is INADMISSIBLE when its worst untreated phenotype reaches "
             f"the manuscript's stated ceiling of {_pct(cons['baseline_max'], 0)} "
             f"death under Control -- the same ceiling, imported from the same "
             f"generator, that `headline_at_fitted.py` uses to reject the "
             f"in-vitro-fitted cascade, judged on "
             f"{', '.join('`' + c + '`' for c in cons['controls'])} at the "
             f"published per-draw cell count. The untreated arms do not appear "
             f"anywhere in the published design, so a draw's admissibility was "
             f"not knowable from any committed artifact.\n")
    L.append("| family | published draws | admissible | share |")
    L.append("|---|---|---|---|")
    for fam in FAMILY_SAMPLES:          # the published families, in one order
        f = fams[fam]
        L.append(f"| {fam} | {f['draws']} | {f['n_admissible']} | "
                 f"{_pct(f['fraction'])} |")
    L.append("")
    w = d["worst_untreated"]
    L.append(f"Over the widest draw set the median worst-untreated death rate is "
             f"{_pct(w['median'])} and {w['n_over_half']} draws kill more than "
             f"half the untreated population before any treatment is applied.\n")

    # The rows above are NOT independent estimates and the page must not read
    # as though they were: families drawing the same count draw the SAME rows,
    # and a smaller family's draws are a prefix of a larger one's.
    counts = sorted({f["draws"] for f in fams.values()})
    shared = {}
    for name in FAMILY_SAMPLES:         # the design's order, not the file's
        shared.setdefault(fams[name]["draws"], []).append(name)
    dupes = [names for count, names in sorted(shared.items()) if len(names) > 1]
    L.append(f"Those rows are NESTED, not independent: "
             + ("; ".join(f"`{'` and `'.join(names)}` draw the same rows"
                          for names in dupes) + ", and " if dupes else "")
             + f"the smaller sample is a PREFIX of the larger, so there "
               f"{'are' if len(counts) != 1 else 'is'} {len(counts)} distinct "
               f"draw set{'s' if len(counts) != 1 else ''} here rather than "
               f"{len(fams)}. That is deliberate -- one admissibility pass then "
               f"covers every family, and the partitions cannot disagree "
               f"because of sampling.\n")

    ag = d["agreement"]
    L.append(f"The prior is ONE box, so these are samples of one quantity that "
             f"the sibling page measured on {ag['sibling_draws']} draws "
             f"({_pct(ag['sibling_fraction'])}). Checked rather than asserted: "
             f"the largest deviation here is {ag['max_abs_z']:.2f} standard "
             f"errors against a bar of {ag['bar']:.0f}, so the two pages "
             f"**{'agree' if ag['agrees'] else 'DISAGREE'}** about the prior. A "
             f"real disagreement would mean the two definitions of admissible "
             f"had drifted apart, which is the failure this check exists for.\n")

    L.append("## The published intervals, reproduced and then conditioned\n")
    L.append(f"The full-prior half must reproduce each committed report's own "
             f"figures to within {_num(d['tolerance'], 4)}, a tolerance DERIVED "
             f"from the reports' own printed precision rather than chosen. "
             f"Without that, this page would be correcting a re-run rather than "
             f"what is published; a family that fails it has its comparison "
             f"withheld.\n")

    L.append(f"| headline | published median | published {conf} | admissible "
             f"median | admissible {conf} | reproduces |")
    L.append("|---|---|---|---|---|---|")
    for o in d["observables"]:
        pub, a = o["published"], o["admissible_only"]
        pub_i = (f"[{_num(pub['low'])}, {_num(pub['high'])}]" if pub else "n/a")
        if a is None:
            adm_i = "n/a"
        elif not o["interval_quotable"]:
            adm_i = f"withheld (n={a['n']})"
        else:
            adm_i = f"[{_num(a['p2_5'])}, {_num(a['p97_5'])}]"
        L.append(f"| {o['label']} | {_num(pub['median']) if pub else 'n/a'} | "
                 f"{pub_i} | {_num(a['median']) if a else 'n/a'} | {adm_i} | "
                 f"{_yesno(o['reproduces'])} |")
    L.append("")
    L.append(f"Reproduction: {sum(1 for o in d['observables'] if o['reproduces'])} "
             f"of {len(d['observables'])} headline observables, each on every "
             f"published quantile and both range endpoints.\n")

    if d["withheld_tails"]:
        thin = min(d["families"], key=lambda k: d["families"][k]["draws"])
        fam = d["families"][thin]
        floor = next(o["n_admissible_floor"] for o in d["observables"]
                     if not o["interval_quotable"])
        L.append(f"**Those tails are withheld, and the reason is structural.** "
                 f"Quoting a conditioned tail needs at least {floor} admissible "
                 f"draws (the sibling page's floor, imported rather than chosen "
                 f"here), and the `{thin}` family's published design draws "
                 f"{fam['draws']} in TOTAL because each costs minutes. So it "
                 f"could only ever meet that floor if the constraint removed "
                 f"NOTHING from it: in any world where conditioning does "
                 f"anything at all, this family cannot quote a conditioned "
                 f"tail. Here it leaves {fam['n_admissible']}, and a "
                 f"{_confidence_label()} tail resting on that many draws is not "
                 f"a bound. Its median is quoted and its tails are not. The "
                 f"published report already says to read its median rather than "
                 f"its bounds; this says the conditioned version cannot offer "
                 f"bounds at all at that draw count.\n")

    L.append("### How far each quantile moves\n")
    L.append("A shift counts as MATERIAL when it exceeds "
             f"{_pct(d['material_fraction'], 0)} of that observable's own "
             "published interval width. These observables are on different "
             "scales -- a ratio, a signed difference, and kill fractions -- so "
             "a single absolute bar would mean a different effective bar for "
             "each of them.\n")
    L.append("| headline | quantile | full prior | admissible | change | material |")
    L.append("|---|---|---|---|---|---|")
    for o in d["observables"]:
        for qname in ("low", "median", "high"):
            sh = o["shifts"].get(qname)
            if not sh:
                continue
            verdict = (_yesno(sh["material"]) if sh["quotable"]
                       else "not quoted (too few admissible draws)")
            L.append(f"| {o['label']} | {qname} | {_num(sh['full'])} | "
                     f"{_num(sh['admissible'])} | {_signed(sh['delta'])} | "
                     f"{verdict} |")
    L.append("")

    L.append("### The correction does not run one way\n")
    rising = [o for o in d["observables"]
              if o["shifts"].get("median", {}).get("delta", 0) > 0]
    falling = [o for o in d["observables"]
               if o["shifts"].get("median", {}).get("delta", 0) < 0]
    L.append(f"Admissibility keeps the draws whose defences are strong enough "
             f"that an untreated cell survives, which is a single mechanism -- "
             f"and across these headlines it moves "
             f"{len(falling)} medians DOWN and {len(rising)} UP. "
             + (f"Down: {', '.join(o['label'] for o in falling)}. " if falling else "")
             + (f"Up: {', '.join(o['label'] for o in rising)}. " if rising else "")
             + "Reporting whichever was looked at first would have flattered the "
               "page in whichever direction that happened to be.\n")

    crossed = [o for o in d["observables"] if (o["point"] or {}).get("crossed")]
    if crossed:
        L.append("The sharpest consequence is where this leaves the "
                 "manuscript's own unperturbed point estimate. For "
                 + ", ".join(f"**{o['label']}** (point {_num(o['point']['value'])}, "
                             f"full-prior median {_num(o['full']['median'])}, "
                             f"admissible median "
                             f"{_num(o['admissible_only']['median'])})"
                             for o in crossed)
                 + " the median crosses from one side of the published point to "
                   "the other, so the conditioned prior disagrees with the "
                   "full prior about whether the headline figure is optimistic "
                   "or pessimistic.\n")
    # Three states, not two: an observable whose tails are withheld cannot be
    # said to contain the point either way, and calling that "outside" would
    # manufacture a finding out of a missing bound.
    inside = [o for o in d["observables"]
              if o["point"] and o["point"]["inside_admissible_interval"] is True]
    outside = [o for o in d["observables"]
               if o["point"] and o["point"]["inside_admissible_interval"] is False]
    unassessed = [o for o in d["observables"]
                  if o["point"] and o["point"]["inside_admissible_interval"] is None]
    L.append(f"The published point estimate still falls inside the admissible "
             f"interval for {len(inside)} of {len(d['observables'])} observables"
             + (f", falls OUTSIDE it for {', '.join(o['label'] for o in outside)}"
                if outside else "")
             + (f", and cannot be placed either way for "
                f"{', '.join(o['label'] for o in unassessed)}, whose tails are "
                f"withheld above" if unassessed else "")
             + ".\n")

    L.append("## Is the Bliss move real, or is its denominator collapsing?\n")
    L.append("`synergy_score` is a RATIO of the combination kill to the Bliss "
             "independence prediction, and conditioning lowers kill, so a rise "
             "in the ratio could be the numerator holding or the denominator "
             "falling away. The parts come from the same row as the score, so "
             "the question is answerable rather than arguable.\n")
    L.append("| part | full-prior median | admissible median | ratio |")
    L.append("|---|---|---|---|")
    for col in BLISS_PART_COLUMNS:
        pr = dec["parts"][col]
        adm = pr["admissible_only"]
        r = (adm["median"] / pr["full"]["median"]
             if adm and pr["full"]["median"] else None)
        L.append(f"| `{col}` | {_num(pr['full']['median'])} | "
                 f"{_num(adm['median']) if adm else 'n/a'} | "
                 f"{_num(r) if r is not None else 'n/a'} |")
    L.append("")
    L.append(f"Both move down. The denominator falls to "
             f"{_num(dec['denominator_ratio'])} of its full-prior median and the "
             f"numerator to {_num(dec['numerator_ratio'])}, so the ratio rises "
             f"**because the Bliss baseline falls faster than the combination "
             f"kill does** -- "
             f"{'the denominator is the larger mover' if dec['ratio_follows_denominator'] else 'the numerator is the larger mover'}"
             f". That is a mechanism, not a measurement of more synergy, and it "
             f"is the honest reading of the rise.\n")
    L.append(f"The engine refuses the ratio entirely below a Bliss baseline of "
             f"{_num(dec['floor'], 4)} (parsed from `sim-combo-mech`, not typed "
             f"here). The smallest admissible denominator is "
             f"{_num(dec['admissible_min_denominator'], 4)}, which "
             f"{'clears' if dec['denominator_clears_floor'] else 'does NOT clear'} "
             f"it, so the conditioned interval is not being computed next to the "
             f"undefined corner. The engine's own identity for that baseline "
             f"reproduces on the stored parts to "
             f"{_num(dec['identity_max_error'], 9)}, so the columns are the "
             f"parts they are labelled as.\n")

    L.append("### And the lower bound rises for a different reason\n")
    L.append(f"The published report's conclusion is about the lower BOUND "
             f"rather than the median: it reads the supra-additive direction as "
             f"robust because the interval does not reach below the no-synergy "
             f"line at {_num(dec['no_synergy_line'], 1)}, where the combination "
             f"kill exactly equals the independence prediction. "
             f"{dec['n_at_or_below_line']} of the full prior's draws sit at or "
             f"below that line, and "
             + ("**not one of them is admissible**"
                if dec["sub_additive_corner_is_wholly_inadmissible"] else
                f"{dec['n_at_or_below_line_admissible']} of them are admissible")
             + f". So the sub-additive corner of this prior is "
               f"{'entirely' if dec['sub_additive_corner_is_wholly_inadmissible'] else 'partly'} "
               f"made of draws the model itself rejects, which is why the "
               f"conditioned lower end clears the line.\n")
    if dec["n_on_line"] and dec["min_single_arm_kill_on_line"] is not None:
        L.append(f"And the draws sitting EXACTLY on the line are there for a "
                 f"mechanical reason rather than a biological one: every one of "
                 f"them has a single-agent RSL3 kill of at least "
                 f"{_num(dec['min_single_arm_kill_on_line'])}, so both arms "
                 f"saturate, the independence prediction is driven to its own "
                 f"ceiling, and the ratio is "
                 f"{_num(dec['no_synergy_line'], 1)} by arithmetic. A saturated "
                 f"arm is also the easiest way to fail the untreated "
                 f"constraint, which is why these are the first draws "
                 f"conditioning removes.\n")
    if (dec["corner_min_worst_untreated"] is not None
            and dec["corner_min_larger_single_arm"] is not None
            and dec["corner_min_independence_prediction"] is not None):
        L.append(f"THAT EXPLAINS THE DRAWS ON THE LINE AND NOT THE WHOLE "
                 f"CORNER, so the corner is measured across every draw in it. "
                 f"MEASURED, NOT ENTAILED, and the distinction matters because "
                 f"a low-kill pair admits a sub-additive ratio perfectly "
                 f"well: where both arms are weak the independence prediction "
                 f"is small too, and any combination falling below it undercuts "
                 f"the line, so nothing forces the corner to sit where it does. "
                 f"In THIS "
                 f"prior it does: every draw at or below the line has an arm "
                 f"killing at least "
                 f"{_num(dec['corner_min_larger_single_arm'])} against an "
                 f"independence prediction of at least "
                 f"{_num(dec['corner_min_independence_prediction'])}, and every "
                 f"one of them fails the untreated constraint by a wide margin "
                 f"rather than at its boundary -- the LEAST inadmissible draw "
                 f"in the corner still kills "
                 f"{_pct(dec['corner_min_worst_untreated'], 2)} of its "
                 f"untreated population, "
                 f"{dec['corner_min_worst_untreated'] / d['constraint']['baseline_max']:.1f}x "
                 f"the limit. That the two co-occur is CONSISTENT with one "
                 f"strong-peroxidation regime driving both, and the selection "
                 f"table above is the evidence for that reading rather than "
                 f"this paragraph: the constraint selects on the peroxidation "
                 f"and defence parameters and leaves the drug parameters alone. "
                 f"Each figure is a MINIMUM over the corner, so it is a floor "
                 f"every draw in it clears rather than a mean one draw could "
                 f"carry.\n")
    L.append(f"The two Bliss movements are therefore NOT one finding: the "
             f"MEDIAN rises because the denominator falls faster than the "
             f"numerator, and the LOWER BOUND rises because the sub-additive "
             f"corner is inadmissible. Only the second strengthens the "
             f"published conclusion, and conflating them would read as this "
             f"page having found more synergy.\n")

    still = [k for k, present in sorted(d["subadditive_citations"].items())
             if present]
    if still and dec["sub_additive_corner_is_wholly_inadmissible"]:
        # SENTENCES AND FILES ARE DIFFERENT COUNTS and the earlier version of
        # this paragraph reported the second under the first's name, which is
        # how four of the manuscript's five such sentences went unnamed (#942).
        files = sorted({k.split(" :: ", 1)[0] for k in still})
        per_file = {f: sum(1 for k in still if k.startswith(f + " :: "))
                    for f in files}
        L.append(f"**{_counted(len(still), cap=True)} across "
                 f"{_counted(len(files), 'file')} draw a conclusion from those "
                 f"sub-additive draws, and they are NOT edited here**: "
                 + ", ".join(
                     f"`{f}`" + (f" ({_counted(per_file[f])})"
                                 if per_file[f] > 1 else "")
                     for f in files)
                 + ". Each says the interval reaches the additive null and that "
                   "some draws fall below it, so synergy is not guaranteed. "
                   "That is TRUE of the full prior they cite and this page does "
                   "not contradict it -- but every one of those draws is "
                   "inadmissible, so the conclusion rests on a corner the model "
                   "rejects. They are handwritten interpretation rather than "
                   "generated prose, so they are named here rather than "
                   "rewritten, and the naming is checked against their current "
                   "text SENTENCE BY SENTENCE so it cannot outlive the "
                   "sentences it points at -- which is what this registry says "
                   "it does and what, keyed on files, it did not do.\n")

    L.append("## The ordering check, registered in advance as NOT evidence\n")
    L.append(f"The penetration report's load-bearing claim is a within-draw "
             f"ordering, which the penetration-access analysis (#926) showed to "
             f"be robust BY CONSTRUCTION: "
             f"one shared monotone dose-response read at decreasing fixed "
             f"concentrations, so conditioning CANNOT break it, and that was "
             f"registered on the issue before any observable here existed. "
             f"Measured: {_pct(order['full']['fraction'])} of "
             f"{order['full']['n']} full-prior draws and "
             f"{_pct(order['admissible_only']['fraction'])} of "
             f"{order['admissible_only']['n']} admissible ones. The full-prior "
             f"figure {'reproduces' if order['reproduces_published'] else 'does NOT reproduce'} "
             f"the published one. Read it as a check that this harness evaluates "
             f"the same observable, not as the ordering surviving a test.\n")

    L.append("## The circularity falsifier\n")
    L.append("Admissibility keeps strong-defence draws, which is a claim about "
             "the PARAMETERS, so it is measured as one: each parameter's mean in "
             "the admissible subset against its full-prior mean, as a fraction "
             "of its own range.\n")
    L.append("| parameter | full-prior mean | admissible mean | shift | strongly selected |")
    L.append("|---|---|---|---|---|")
    for name, row in sorted(sel["per_parameter"].items(),
                            key=lambda kv: -abs(kv[1]["shift"])):
        mark = " (drug)" if row["is_drug_parameter"] else ""
        L.append(f"| `{name}`{mark} | {_num(row['full_mean_fraction'])} | "
                 f"{_num(row['admissible_mean_fraction'])} | "
                 f"{_signed(row['shift'])} | {_yesno(row['strongly_selected'])} |")
    L.append("")
    L.append(f"The falsifier is the DRUG parameters "
             f"({', '.join('`' + n + '`' for n in sel['drug_parameters'])}): had "
             f"conditioning pushed them down, the kill reductions above would be "
             f"partly an artifact of keeping weaker-drug draws and the whole "
             f"comparison would be circular. Largest drug-parameter shift "
             f"{_num(sel['max_drug_shift'])} against a bar of "
             f"{_num(STRONG_SELECTION, 2)}: "
             f"**{'untouched' if sel['drug_parameters_untouched'] else 'SELECTED, so the comparison is circular'}**. "
             f"The kill falls because the defences kept are stronger, not "
             f"because the drug kept is weaker.\n")

    L.append("## Does this move any per-headline verdict?\n")
    flips = d["sign_flips"]
    L.append(f"`identifiability-report.md` reads these intervals for its "
             f"per-headline verdicts, and those are DIRECTION claims rather "
             f"than magnitudes -- so what would move one is a SIGN change, not "
             f"a shifted quantile. Measured: "
             + (f"the conditioned median changes sign for "
                f"{', '.join(flips)}, so those verdicts have to be re-examined."
                if flips else
                "no conditioned median changes sign, so no verdict in that "
                "report is contradicted by this page. The magnitudes it already "
                "declines to quote are the ones that move.")
             + "\n")

    L.append("## What this does NOT establish\n")
    L.append("**Conditioning is not inference.** Removing draws the model itself "
             "declares invalid adds no data. It cannot turn a prior-predictive "
             "interval into a data-conditioned one, and "
             "`identifiability-report.json` still reports "
             f"{d['data_conditioned_headlines']} data-conditioned headline "
             "outputs -- unchanged by everything on this page, and guarded, so "
             "this refusal breaks if that report ever says otherwise.\n")
    L.append("**The admissible subset is not a posterior.** It is the same "
             "uniform box with a region removed, so its median is a property of "
             "what remains rather than an estimate of anything.\n")
    L.append("**The constraint is a modelling choice**, stated in the "
             "manuscript's prose and imported from the generator whose own guard "
             "parses it from there. It is not an experimental observation.\n")
    # Which family is thinnest is a property of the published design, so it is
    # derived: naming one here would be a sentence free to outlive the design.
    thin = min(fams, key=lambda k: fams[k]["draws"])
    L.append(f"**The `{thin}` family cannot be given a conditioned interval "
             f"at its published draw count.** It draws {fams[thin]['draws']} "
             f"because each costs minutes, and conditioning leaves "
             f"{fams[thin]['n_admissible']}. Its median is reported and its "
             f"tails are withheld, which is a statement about the published "
             f"design's resolution rather than about the biology.\n")
    L.append("**Nothing here recalibrates a parameter.** The prior's ranges are "
             "unchanged; only draws the model rejects are removed.\n")

    L.append("## Verdict\n")
    L.append(f"{d['n_reportable']} of {len(d['observables'])} headline "
             f"observables are reportable -- each reproduces what is published "
             f"and has admissible draws to condition on -- of which "
             f"{len(d['observables']) - len(d['withheld_tails'])} can also have "
             f"their tails quoted, and {d['n_material']} carry at least one "
             f"material shift among the quantiles this page will quote. "
             f"Overall: **{'reportable' if d['verdict'] else 'WITHHELD'}**.\n")
    return "\n".join(L)




def _split_stored(stored):
    """Recover `assemble`'s inputs for the freshness gate's round trip.

    Every input is committed, so this re-reads them rather than echoing stored
    derived fields back: an echo would make every derived-field guard inert,
    which is a defect this repo has already paid for once.
    """
    return (load(),)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--measure", choices=sorted(LEGS),
                    help="take one leg's measurements (needs the binaries and the "
                         "Python extension); legs share one file, so run them "
                         "ONE AT A TIME")
    ap.add_argument("--binary-dir", default=str(REPO / "simulations" / "target" / "release"))
    ap.add_argument("--workers", type=int, default=min(8, (os.cpu_count() or 4)))
    ap.add_argument("--render-only", action="store_true",
                    help="accepted for consistency with the other generators, "
                         "and a no-op here: see the note in main()")
    a = ap.parse_args()
    if a.measure:
        measure(a.measure, a.binary_dir, a.workers)
        # Re-render, so a measurement never leaves the page describing the
        # previous one.
    # The two branches below are the SAME call, and that is deliberate rather
    # than sloppy: elsewhere `--render-only` selects the cheap path because the
    # default one MEASURES, while here measuring is its own flag and every
    # remaining input is committed, so the default run already re-derives the
    # whole page from files. The branch is written out instead of collapsed
    # because `tests/test_artifact_freshness.py` keys two real round trips on
    # seeing it -- that `assemble` is idempotent, and that the rendered page
    # does not depend on the stored dicts' key order -- and both are worth
    # having. `_split_stored` below is what those round trips call.
    if a.render_only:
        d = assemble(load())
    else:
        d = assemble(load())
    OUT_JSON.write_text(json.dumps(d, indent=2, sort_keys=True) + "\n")
    OUT_MD.write_text(render(d))
    print(f"wrote {OUT_MD.relative_to(REPO)} and {OUT_JSON.relative_to(REPO)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
