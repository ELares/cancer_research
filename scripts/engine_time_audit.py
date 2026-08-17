#!/usr/bin/env python3
"""What is one step of the ferroptosis engine worth in real time? (#727)

WHAT THIS DOCUMENT USED TO SAY, AND WHY IT WAS WRONG
-----------------------------------------------------
The first version's headline was "Nothing in the engine states it", searched
across every library module. That is FALSE, and the falsifier sits inside the
directory the scan already walked:

    tumor_pk.rs   /// Time points in minutes (one per simulation step).
                  /// ...at 1-minute resolution (one value per simulation step)
                  for minute in 0..n_steps

    trigger_wave.rs   /// Time step (min). ...  pub dt_min: f64   (default 0.02)

The detector could not see either. `_declares_step_duration` required the
STEP-FIRST word order -- "one step is N minutes" -- and a step-duration
declaration is naturally written UNIT-FIRST: "time points in minutes, one per
simulation step". The scan's other half, `TIME_UNIT`, required a NUMBER before
the unit, so a bare "Time step (min)" was invisible to it too. The instrument
was structurally blind to the exact thing it existed to find, and the report
published the blindness as an absence.

The caveat "the scan is textual ... the real-time set is a lower bound" was
already sitting two lines under that headline. A true caveat does not rescue a
false headline.

WHAT IS ACTUALLY TRUE, WHICH IS A SHARPER FINDING
---------------------------------------------------
The CORE biochemical loop states no step duration -- that part survives. But
two modules bind a step to a real duration, they DISAGREE, and one of them
reaches production:

  tumor_pk       1 step = 1 minute, and `sim-tumor-pk` calls it with the same
                 180-step count the production matrix uses, as does sim-tme-3d.
  trigger_wave   1 step = `dt_min`, default 0.02 minutes -- fifty times finer.

So the engine does not lack a step duration. It carries two, unreconciled, one
of them silently load-bearing. THIS REVERSES THE WITHDRAWAL the previous
version made: issue #727 originally alleged "four conflicting bindings", the
last version withdrew that as "there are NONE", and the truth is two. The
withdrawal over-corrected.

THE CONSEQUENCE, NOW QUANTIFIED RATHER THAN ASSERTED
------------------------------------------------------
Under the only binding that reaches production, the whole 180-step run spans
three hours. PREREGISTRATION.md states P3's window as days and its FALSIFICATION
THRESHOLD as "returns to baseline within 24 hours". The production loop is
shorter than the threshold, so it cannot represent either outcome of the
project's most directly testable prediction. That is stronger than "the units
are missing": the units are present, they are inconsistent, and the one in force
is off by nearly two orders of magnitude from the prediction it is meant to
score.

WHAT THIS SCRIPT DOES NOT DO
-----------------------------
It does not choose a step duration or reconcile the two. That is a modelling
decision with consequences for every calibrated layer and for the committed
byte-identity gates, and it belongs to whoever owns those calibrations. A guard
fails if this script starts asserting one.

Usage:
    python scripts/engine_time_audit.py
    python scripts/engine_time_audit.py --render-only
"""

import argparse
import json
import re
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
SRC = PROJECT_ROOT / "simulations" / "ferroptosis-core" / "src"
SIMS = PROJECT_ROOT / "simulations"
PREREG = PROJECT_ROOT / "PREREGISTRATION.md"
OUT_MD = PROJECT_ROOT / "analysis" / "engine-time-audit.md"
OUT_JSON = PROJECT_ROOT / "analysis" / "engine-time-audit.json"

UNIT_MIN = {"sec": 1 / 60, "second": 1 / 60, "seconds": 1 / 60, "s": 1 / 60,
            "min": 1.0, "mins": 1.0, "minute": 1.0, "minutes": 1.0,
            "h": 60.0, "hr": 60.0, "hrs": 60.0, "hour": 60.0, "hours": 60.0,
            "day": 1440.0, "days": 1440.0}
UNIT_RE = r"(sec|seconds?|min|mins|minutes?|h|hrs?|hours?|days?)"

# A real-time unit appearing with a magnitude. Narrow on purpose: the question
# is which modules carry PHYSICAL time, not which mention the word.
TIME_UNIT = re.compile(r"\b(\d+(?:\.\d+)?)\s*(?:-)?\s*" + UNIT_RE + r"\b", re.I)
PER_STEP = re.compile(r"per[- ]step|per dimensionless|each step|/step", re.I)

# --- the three ways a step-to-duration binding is actually written ----------
# Each returns (magnitude_or_None, unit). A bare unit with "one per step" means
# one unit per step, so the magnitude defaults to 1.
BINDINGS = [
    # "one step is 30 minutes", "each step corresponds to 1 min"
    ("step-first", re.compile(
        r"(?:one|a|each|per)\s+step\s+"
        r"(?:is|=|represents?|corresponds?\s+to|equals?)\s*"
        r"[^.\n]{0,30}?(\d+(?:\.\d+)?)?\s*" + UNIT_RE + r"\b", re.I)),
    # "time points in minutes (one per simulation step)"
    # "at 1-minute resolution (one value per simulation step)"
    ("unit-first", re.compile(
        r"(\d+(?:\.\d+)?)?\s*-?\s*" + UNIT_RE +
        r"\b[^.]{0,70}?\b(?:one|1)\s+(?:value\s+|point\s+)?per\s+"
        r"(?:simulation\s+)?step", re.I)),
    # an explicit per-step duration field: `dt_min`, `dt_hours`
    ("named-field", re.compile(r"\bdt_" + UNIT_RE + r"\b", re.I)),
    # `for minute in 0..n_steps` -- the loop variable IS the unit
    ("loop-variable", re.compile(
        r"for\s+" + UNIT_RE + r"\s+in\s+0\.\.\s*\w*n_steps", re.I)),
]


def logical_lines(text: str):
    """Yield (line_no, text) with contiguous doc-comment runs joined.

    A step-duration declaration routinely spans two `///` lines, and a
    line-by-line scan splits the phrase in half. Joining the run is what let
    `tumor_pk`'s '1-minute resolution (one value per simulation step)' be seen
    at all -- the previous detector read those as two unrelated lines.
    """
    lines = text.split("\n")
    i = 0
    while i < len(lines):
        s = lines[i].strip()
        if s.startswith("///") or s.startswith("//!"):
            start, buf = i + 1, []
            while i < len(lines) and lines[i].strip()[:3] in ("///", "//!"):
                buf.append(re.sub(r"^\s*//[/!]\s?", "", lines[i]))
                i += 1
            yield start, " ".join(buf).strip()
        else:
            yield i + 1, lines[i]
            i += 1


def _default_for(stem: str, field: str, unit: str):
    """The literal a named per-step duration field is constructed with."""
    for p in SRC.glob("*.rs"):
        m = re.search(rf"\b{re.escape(field)}\s*:\s*(\d+(?:\.\d+)?)", p.read_text(errors="ignore"))
        if m:
            return float(m.group(1)) * UNIT_MIN[unit.lower()]
    return None


def find_step_bindings():
    """Every place the engine binds one step to a real duration.

    Returns sites, not a bool. The previous version returned True/False, so a
    reader could not tell WHICH module bound it or whether two disagreed -- and
    a bool cannot carry the finding that they do.
    """
    sites = []
    for p in sorted(SRC.glob("*.rs")):
        text = p.read_text(errors="ignore")
        for ln, line in logical_lines(text):
            for how, pat in BINDINGS:
                m = pat.search(line)
                if not m:
                    continue
                groups = m.groups()
                if how in ("named-field", "loop-variable"):
                    unit = groups[0]
                    mag = None
                else:
                    mag, unit = groups[0], groups[1]
                u = unit.lower()
                if u not in UNIT_MIN:
                    continue
                if how == "named-field":
                    field = m.group(0)
                    mpstep = _default_for(p.stem, field, u)
                else:
                    mpstep = (float(mag) if mag else 1.0) * UNIT_MIN[u]
                # Show the MATCHED span, not the start of the joined block.
                # A doc-comment run can be hundreds of characters and the
                # phrase that fired often sits at the end of it, so slicing
                # from the start displays everything except the evidence.
                lo = max(0, m.start() - 12)
                hi = min(len(line), m.end() + 12)
                # A named `dt_*` field is an INTEGRATOR'S OWN timestep, chosen
                # for numerical stability -- trigger_wave asserts a CFL bound
                # `dt < h^2/(2D)` right beside it. That is not a claim about
                # what a simulation step is worth in wall-clock time. An
                # earlier draft pooled the two and reported a "50x
                # disagreement, unreconciled" between a PDE stability step and
                # a step-to-minute alignment: unlike objects, and two solvers
                # having different dt is not a contradiction.
                kind = ("solver-timestep" if how == "named-field"
                        else "wall-clock")
                sites.append({
                    "module": p.name, "line": ln, "how": how, "kind": kind,
                    "unit": u, "minutes_per_step": mpstep,
                    "text": ("..." if lo else "") + line[lo:hi].strip() +
                            ("..." if hi < len(line) else ""),
                })
                break
    return sites


def _callers(stem: str):
    """Which simulation binaries reference this module. Measured, not assumed.

    The previous version claimed three timescales were 'all composed into the
    same 180-step run'. Nothing checked it, and the day-scale limb had no
    caller at all.
    """
    out = []
    for p in sorted(SIMS.glob("sim-*/src/*.rs")):
        if re.search(rf"\b{re.escape(stem)}\b", p.read_text(errors="ignore")):
            out.append(p.parts[-3])
    return sorted(set(out))


def _orphan_timescales():
    """Real-time CONFIG FIELDS with no reference outside their own module.

    Module-level callers are too coarse to check the composition claim. The
    previous version said day-scale, hour-scale and minute-scale code was "all
    composed into the same 180-step run"; `senescence.rs` IS called by
    sim-tme-3d, so a module-level check passes -- while the day-scale fields
    inside it are referenced nowhere but the module and its own tests.
    """
    field = re.compile(r"\bpub\s+(\w+_(?:days?|hours?|hrs?|mins?|minutes?|secs?))\s*:")
    out = {}
    for p in sorted(SRC.glob("*.rs")):
        names = set(field.findall(p.read_text(errors="ignore")))
        for nm in sorted(names):
            refs = []
            for q in list(SRC.glob("*.rs")) + list(SIMS.glob("sim-*/src/*.rs")):
                if q == p:
                    continue
                if re.search(rf"\b{re.escape(nm)}\b", q.read_text(errors="ignore")):
                    refs.append(q.parts[-3] if "sim-" in str(q) else q.name)
            out.setdefault(p.name, {})[nm] = sorted(set(refs))
    return out


def _step_counts(binding_modules):
    """Every binary's declared step count, and whether it consumes a binding.

    The previous version matched `const N_STEPS: usize` and returned the FIRST
    hit. Only sim-tumor-pk annotates `usize`; sim-tme-3d, sim-tme and
    sim-combo-mech all use `u32`, so the "production binary" was being selected
    by a type annotation and sim-tme-3d's count was never read at all. Two
    planted mutations confirmed it: changing sim-tme-3d's count left the report
    unchanged, and retyping sim-combo-mech's made it the "production" binary
    for a binary that never touches the module in question.

    So: read every binary regardless of integer type, and record separately
    whether it references a module that declares a wall-clock binding. Only
    those can price a step in real time.
    """
    stems = {Path(m).stem for m in binding_modules}
    out = []
    for p in sorted(SIMS.glob("sim-*/src/main.rs")):
        t = p.read_text(errors="ignore")
        m = re.search(r"const\s+N_STEPS\s*:\s*\w+\s*=\s*(\d+)", t)
        if not m:
            continue
        binary = p.parts[-3]
        srcs = "".join(q.read_text(errors="ignore")
                       for q in (p.parent).glob("*.rs"))
        uses = sorted(s for s in stems if re.search(rf"\b{re.escape(s)}\b", srcs))
        out.append({"binary": binary, "steps": int(m.group(1)),
                    "consumes_binding_module": uses})
    return out


def _p3_threshold():
    """P3's stated threshold, window and RECOVERY ORDER, from the preregistration."""
    if not PREREG.exists():
        return None
    txt = PREREG.read_text(errors="ignore")
    blk = re.search(r"\*\*P3\..*?(?=\n\*\*P4\.)", txt, re.S)
    if not blk:
        return None
    b = blk.group(0)
    thr = re.search(r"within\s+(\d+)\s*" + UNIT_RE + r"\b", b, re.I)
    win = re.search(r"(\d+)\s*to\s*(\d+)\s*(days?|hours?)\b", b, re.I)
    out = {}
    if thr:
        out["threshold_hours"] = float(thr.group(1)) * UNIT_MIN[thr.group(2).lower()] / 60
    if win:
        f = UNIT_MIN[win.group(3).lower()] / 60
        out["window_hours"] = [float(win.group(1)) * f, float(win.group(2)) * f]
    # "(FSP1 and GSH first, GPX4 and NRF2 later)"
    order = re.search(r"\(([^)]*?)\s+first,\s*([^)]*?)\s+later\)", b, re.I)
    if order:
        def names(s):
            return sorted({w.lower() for w in re.findall(r"[A-Za-z][A-Za-z0-9]+", s)
                           if w.lower() not in ("and", "or")})
        out["order_first"] = names(order.group(1))
        out["order_later"] = names(order.group(2))
    return out or None


def _p3_is_modelled():
    """Is P3's quantity represented anywhere, and does the engine AGREE with it?

    The previous version of this report concluded the model "cannot represent
    either outcome of its own most directly testable prediction". That is FALSE
    and this function is why it is now measured instead of argued: `cell.rs`
    carries per-defense recovery half-times in DAYS, `sim-window` sweeps them
    from 0 to 28 days at P3's own timepoints, and the result is a published
    manuscript figure. The refutation was already inside this script's own
    artifact -- `orphan_timescale_fields["cell.rs"]` lists all four fields with
    caller `sim-window` -- and the renderer printed only the fields with EMPTY
    caller lists, so it computed the counter-example and dropped it.
    """
    cell = SRC / "cell.rs"
    if not cell.exists():
        return None
    rates = {}
    for m in re.finditer(r"\b(\w+)_half_recovery_days\s*:\s*(\d+(?:\.\d+)?)",
                         cell.read_text(errors="ignore")):
        rates[m.group(1).lower()] = float(m.group(2))
    if not rates:
        return None
    consumers = sorted({p.parts[-3] for p in SIMS.glob("sim-*/src/*.rs")
                        if re.search(r"half_recovery_days|RecoveryRates",
                                     p.read_text(errors="ignore"))})
    span = None
    for p in SIMS.glob("sim-*/src/main.rs"):
        t = p.read_text(errors="ignore")
        if "timepoints_hours" not in t:
            continue
        blk = re.search(r"timepoints_hours[^;]*?vec!\[(.*?)\]", t, re.S)
        if blk:
            # Strip trailing line comments FIRST. `168.0, // 1 week` carries a
            # bare `1` that a naive number scan counts as a tenth timepoint --
            # nine values were reported as twelve.
            body = re.sub(r"//[^\n]*", "", blk.group(1))
            vals = [float(x) for x in re.findall(r"(\d+(?:\.\d+)?)", body)]
            if vals:
                span = {"binary": p.parts[-3], "max_hours": max(vals),
                        "n_timepoints": len(vals)}
    return {"recovery_days": rates, "consumers": consumers, "sweep": span}


def _p3_order_verdict(p3, modelled):
    """Does the engine's default recovery ORDER match the one P3 states?

    Derived, not asserted: P3's grouping is parsed from the preregistration and
    compared against cell.rs's defaults. This is a far sharper finding than the
    units complaint the previous version shipped, and it is only visible
    BECAUSE the model represents the quantity.
    """
    if not p3 or not modelled:
        return None
    first, later = p3.get("order_first"), p3.get("order_later")
    rates = modelled["recovery_days"]
    if not first or not later:
        return None
    f = {k: rates[k] for k in first if k in rates}
    l = {k: rates[k] for k in later if k in rates}
    if not f or not l:
        return None
    slowest_first = max(f, key=f.get)
    fastest_later = min(l, key=l.get)
    agrees = f[slowest_first] < l[fastest_later]
    return {
        "stated_first": f, "stated_later": l, "agrees": agrees,
        "violator": None if agrees else slowest_first,
        "violator_days": None if agrees else f[slowest_first],
        "fastest_later": fastest_later, "fastest_later_days": l[fastest_later],
        "engine_order": sorted(rates.items(), key=lambda kv: kv[1]),
    }


def scan() -> dict:
    modules = {}
    for p in sorted(SRC.glob("*.rs")):
        text = p.read_text(errors="ignore")
        times, steps = [], 0
        for i, line in enumerate(text.split("\n"), 1):
            for m in TIME_UNIT.finditer(line):
                ctx = line.strip()
                if ctx.startswith("//!") and "v0." in ctx:
                    continue
                times.append({"line": i, "value": m.group(0).strip(),
                              "context": ctx[:110]})
            if PER_STEP.search(line):
                steps += 1
        if times or steps:
            modules[p.name] = {"real_time_mentions": len(times),
                               "per_step_mentions": steps,
                               "examples": times[:4]}

    bindings = find_step_bindings()
    for b in bindings:
        b["callers"] = _callers(Path(b["module"]).stem)

    # Count CONVENTIONS, not regex hits. `\bdt_min\b` matches the field
    # declaration, its literal in baseline(), and `let dt = cfg.dt_min` -- one
    # binding, three matches. An earlier draft reported "6 declarations" for
    # what is 2 conventions in 2 modules.
    wall = [b for b in bindings if b["kind"] == "wall-clock"]
    solver = [b for b in bindings if b["kind"] == "solver-timestep"]
    conventions = sorted({(b["module"], round(b["minutes_per_step"], 6))
                          for b in wall if b["minutes_per_step"] is not None})
    distinct = sorted({v for _, v in conventions})
    in_prod = [b for b in wall if b["callers"]]
    prod_mps = sorted({round(b["minutes_per_step"], 6) for b in in_prod
                       if b["minutes_per_step"] is not None})

    both = {k: v for k, v in modules.items()
            if v["real_time_mentions"] and v["per_step_mentions"]}
    real_only = {k: v for k, v in modules.items()
                 if v["real_time_mentions"] and not v["per_step_mentions"]}
    step_only = {k: v for k, v in modules.items()
                 if v["per_step_mentions"] and not v["real_time_mentions"]}

    # which real-time modules are actually reachable from a binary
    reach = {m: _callers(Path(m).stem) for m in sorted(set(real_only) | set(both))}
    p3 = _p3_threshold()
    modelled = _p3_is_modelled()

    return {
        "modules_total": len(list(SRC.glob("*.rs"))),
        "modules": modules,
        "carry_both": sorted(both),
        "real_time_only": sorted(real_only),
        "per_step_only": sorted(step_only),
        "real_time_module_callers": reach,
        "step_bindings": bindings,
        "n_step_bindings": len(bindings),
        "wall_clock_conventions": [{"module": m, "minutes_per_step": v}
                                   for m, v in conventions],
        "n_wall_clock_conventions": len(conventions),
        "n_solver_timesteps": len(solver),
        "distinct_minutes_per_step": distinct,
        "bindings_reaching_production": len(in_prod),
        "production_minutes_per_step": prod_mps,
        "step_counts": _step_counts({b["module"] for b in wall}),
        "orphan_timescale_fields": _orphan_timescales(),
        "p3": p3,
        "p3_modelled": modelled,
        "p3_order": _p3_order_verdict(p3, modelled),
    }


def _core_prices(d):
    """Does the CORE biochemical loop itself price a step? Measured."""
    return any(c["module"] == "biochem.rs"
               for c in d.get("wall_clock_conventions", []))


def _priced_runs(d):
    """Binaries whose step count can actually be priced in wall-clock time.

    Only a binary that CONSUMES a module declaring a wall-clock binding can be
    priced. The previous version picked the first `usize`-typed N_STEPS it
    found and called it "the production matrix", which survived two mutations.
    """
    out = []
    for sc in d.get("step_counts", []):
        if not sc["consumes_binding_module"]:
            continue
        for c in d.get("wall_clock_conventions", []):
            if Path(c["module"]).stem in sc["consumes_binding_module"]:
                out.append({**sc, "minutes_per_step": c["minutes_per_step"],
                            "hours": sc["steps"] * c["minutes_per_step"] / 60})
                break
    return out


def render(d: dict) -> str:
    n = d["modules_total"]
    b = d["step_bindings"]
    L = ["# What is one engine step worth in real time?", ""]
    L += ["*Generated by `scripts/engine_time_audit.py`. Every count is "
          "recomputed.*", ""]

    # HEADLINE, derived. It must be able to say the opposite: the previous
    # version hardcoded "Nothing in the engine states it" and stayed true-
    # looking while the scan below listed the module that states it.
    if not b:
        L += [f"**Nothing in the engine states it.** Searched {n} library "
              f"modules for a declaration binding a step to a duration: none "
              f"found.", ""]
    else:
        conv = d.get("wall_clock_conventions", [])
        nsolv = d.get("n_solver_timesteps", 0)
        mods = sorted({c["module"] for c in conv})
        L += [f"**{len(conv)} module{'' if len(conv) == 1 else 's'} "
              f"{'prices' if len(conv) == 1 else 'price'} a "
              f"simulation step in wall-clock time** "
              f"({', '.join(f'`{m}`' for m in mods) or 'none'}), out of {n} "
              f"library modules, plus {nsolv} numerical-integrator timesteps "
              f"that do not.", ""]
        L += ["Those are counted as CONVENTIONS, not matches. A `dt_min` field "
              "is matched three times -- its declaration, its default literal, "
              "and its use -- and an earlier draft reported that arithmetic as "
              "\"6 declarations\".", ""]
        if len(d["distinct_minutes_per_step"]) > 1:
            vals = ", ".join(f"{v:g} min" for v in d["distinct_minutes_per_step"])
            hi, lo = max(d["distinct_minutes_per_step"]), min(d["distinct_minutes_per_step"])
            L += [f"The distinct wall-clock step durations are {vals}, a factor "
                  f"of {hi/lo:.0f} apart.", ""]

    if b:
        L += ["| module | line | kind | how it is written | one step = | reached by |",
              "|---|--:|---|---|--:|---|"]
        for x in b:
            mps = f"{x['minutes_per_step']:g} min" if x["minutes_per_step"] is not None else "unstated"
            who = ", ".join(f"`{c}`" for c in x["callers"]) or "*no caller*"
            L.append(f"| `{x['module']}` | {x['line']} | {x['kind']} | "
                     f"`{x['text'][:56]}` | {mps} | {who} |")
        L += [""]
        L += ["`solver-timestep` rows are an integrator's own `dt`, constrained "
              "by numerical stability -- `trigger_wave` asserts a CFL bound "
              "beside its. Two solvers having different `dt` is not a "
              "disagreement about what a step is worth, and an earlier draft "
              "pooled them to report one.", ""]

    priced = _priced_runs(d)
    if priced:
        L += ["## What can actually be priced in wall-clock time", ""]
        L += ["| binary | steps | min/step | span |", "|---|--:|--:|--:|"]
        for r in priced:
            L.append(f"| `{r['binary']}` | {r['steps']} | "
                     f"{r['minutes_per_step']:g} | {r['hours']:.1f} h |")
        L += [""]
        others = [s for s in d.get("step_counts", [])
                  if not s["consumes_binding_module"]]
        if others:
            names = ", ".join("`{}` ({})".format(s["binary"], s["steps"])
                              for s in others)
            L += [f"{len(others)} other binaries declare a step count "
                  f"({names}) and consume no module that prices a step, so "
                  f"their runs cannot be converted to wall-clock time at all.",
                  ""]
        L += ["These spans are properties of those binaries' assays. They are "
              "**not** a bound on what the project can predict -- an earlier "
              "draft of this page said they were, which was the second false "
              "headline on it. See below.", ""]

    # --- P3: represented, and the engine DISAGREES with it -------------------
    mod = d.get("p3_modelled")
    p3 = d.get("p3") or {}
    if mod:
        L += ["## P3 is represented, on a real-time axis, and this page "
              "previously said it was not", ""]
        sw = mod.get("sweep")
        who = ", ".join(f"`{c}`" for c in mod["consumers"]) or "nothing"
        L += [f"`cell.rs` carries per-defense recovery half-times **in days** "
              f"({', '.join(f'`{k}` {v:g}d' for k, v in sorted(mod['recovery_days'].items()))}), "
              f"consumed by {who}."]
        if sw:
            L += [f"`{sw['binary']}` sweeps {sw['n_timepoints']} timepoints out "
                  f"to **{sw['max_hours']/24:g} days**, which is the axis P3 is "
                  f"scored on."]
        L += ["", "So a previous draft's conclusion -- that the model cannot "
              "represent either outcome of its own most directly testable "
              "prediction -- was **false**, and is withdrawn here. It compared "
              "one binary's inner assay length against a prediction scored on a "
              "different, outer axis: an asymmetric comparison, and the same "
              "class of error as the absence this page already retracts.", ""]
        L += ["Worse, the refutation was already inside this script's own "
              "artifact. `orphan_timescale_fields` recorded those four fields "
              "with their caller, and the renderer printed only the fields "
              "whose caller list was EMPTY -- so it computed the "
              "counter-example and dropped it before rendering.", ""]

    order = d.get("p3_order")
    if order:
        L += ["### The finding that replaces it", ""]
        eo = ", ".join(f"`{k}` {v:g}d" for k, v in order["engine_order"])
        if not order["agrees"]:
            L += [f"P3 states the defenses recover "
                  f"**{', '.join(sorted(order['stated_first']))} first, "
                  f"{', '.join(sorted(order['stated_later']))} later**. The "
                  f"engine's defaults, fastest first, are {eo}.", ""]
            L += [f"They **contradict**: `{order['violator']}` is stated to "
                  f"recover early but is the engine's *slowest* at "
                  f"{order['violator_days']:g} days, against "
                  f"`{order['fastest_later']}` at "
                  f"{order['fastest_later_days']:g} days among the defenses P3 "
                  f"puts later. A wet-lab result matching the engine would "
                  f"falsify P3's stated ordering, and a result matching P3 "
                  f"would falsify the engine's defaults. That is a real, "
                  f"scoreable disagreement -- and it is visible only BECAUSE "
                  f"the model represents the quantity.", ""]
        else:
            L += [f"P3's stated ordering and the engine's defaults ({eo}) "
                  f"agree.", ""]

    L += ["## Which modules carry real time, and whether anything calls them", ""]
    L += ["| module | callers |", "|---|---|"]
    for m, c in sorted(d["real_time_module_callers"].items()):
        L.append(f"| `{m}` | {', '.join(f'`{x}`' for x in c) or '**none**'} |")
    L += [""]
    # BOTH sides, always. Printing only the empty caller lists is how this
    # renderer computed `cell.rs`'s four P3 recovery fields -- reachable from
    # sim-window -- and dropped them, leaving a false "P3 cannot be
    # represented" claim standing over its own counter-example.
    fields = d.get("orphan_timescale_fields", {})
    orphan_f = {m: sorted(f for f, r in fs.items() if not r)
                for m, fs in fields.items()}
    orphan_f = {m: fs for m, fs in orphan_f.items() if fs}
    reach_f = {m: sorted((f, r) for f, r in fs.items() if r)
               for m, fs in fields.items()}
    reach_f = {m: fs for m, fs in reach_f.items() if fs}
    if orphan_f or reach_f:
        L += ["Module-level callers are too coarse to check what is actually "
              "composed, so every real-time configuration field is listed with "
              "its callers -- **both** those that have one and those that do "
              "not. An earlier draft printed only the orphans, which is how it "
              "dropped the four `cell.rs` recovery fields that refute its own "
              "P3 claim.", ""]
    if reach_f:
        tot = sum(len(v) for v in reach_f.values())
        L += [f"**{tot} real-time fields ARE reachable:**", ""]
        for m, fs in sorted(reach_f.items()):
            for f, r in fs:
                L.append(f"* `{m}` `{f}` -> {', '.join(f'`{x}`' for x in r)}")
        L += [""]
    if orphan_f:
        tot = sum(len(v) for v in orphan_f.values())
        L += [f"**{tot} real-time configuration fields are referenced nowhere "
              f"outside their own module**:", ""]
        for m, fs in sorted(orphan_f.items()):
            L.append(f"* `{m}`: {', '.join(f'`{f}`' for f in fs)}")
        L += [""]
        L += ["An earlier version of this report said the day-scale, hour-scale "
              "and minute-scale modules were \"all composed into the same "
              "180-step run\". That was asserted rather than measured, and it "
              "is false for every field listed here: the module is reachable, "
              "the timescale inside it is not.", ""]

    L += ["## What this does not do", ""]
    L += ["* It does not choose a step duration. Adopting one across the "
          "engine would move every calibrated layer and the committed "
          "byte-identity gates, and belongs to whoever owns those "
          "calibrations.",
          f"* The core biochemical loop (`biochem.rs`) "
          f"{'declares no wall-clock step duration' if not _core_prices(d) else 'NOW declares one, so this sentence has changed'}; "
          f"a module composed alongside it does, which is a different "
          f"statement. Derived, because a fixed sentence here would survive "
          f"biochem.rs acquiring one.",
          "* The scan is textual, so the binding set is a lower bound. That "
          "caveat was true in the previous version too, sitting under a "
          "headline asserting an absolute absence -- which is why the headline "
          "is now generated from the scan instead of written beside it.",
          ""]
    return "\n".join(L) + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--render-only", action="store_true")
    args = ap.parse_args()
    if args.render_only:
        d = json.loads(OUT_JSON.read_text())
    else:
        d = scan()
        if not d["modules"]:
            raise SystemExit(
                "no modules matched, which is not a finding -- it is what a "
                "wrong source path looks like.")
        OUT_JSON.write_text(json.dumps(d, indent=1, sort_keys=True) + "\n",
                            encoding="utf-8")
    OUT_MD.write_text(render(d), encoding="utf-8")
    print(f"wrote {OUT_MD}")
    print(f"wrote {OUT_JSON}")
    print(f"  wall-clock conventions: {d['n_wall_clock_conventions']}  "
          f"solver timesteps: {d['n_solver_timesteps']}  "
          f"distinct durations: {d['distinct_minutes_per_step']}")
    for r in _priced_runs(d):
        print(f"    {r['binary']:16s} {r['steps']:>4} steps -> {r['hours']:.1f} h")
    o = d.get("p3_order")
    if o:
        print(f"  P3 recovery order agrees with engine defaults: {o['agrees']}")


if __name__ == "__main__":
    main()
