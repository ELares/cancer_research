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
                sites.append({
                    "module": p.name, "line": ln, "how": how,
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


def _production_steps():
    """The step count a production binary actually runs, read from source."""
    for p in sorted(SIMS.glob("sim-*/src/main.rs")):
        m = re.search(r"const\s+N_STEPS\s*:\s*usize\s*=\s*(\d+)", p.read_text(errors="ignore"))
        if m:
            return {"binary": p.parts[-3], "steps": int(m.group(1))}
    return None


def _p3_threshold():
    """P3's falsification threshold, parsed from the preregistration."""
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
    return out or None


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

    distinct = sorted({round(b["minutes_per_step"], 6) for b in bindings
                       if b["minutes_per_step"] is not None})
    in_prod = [b for b in bindings if b["callers"]]
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

    return {
        "modules_total": len(list(SRC.glob("*.rs"))),
        "modules": modules,
        "carry_both": sorted(both),
        "real_time_only": sorted(real_only),
        "per_step_only": sorted(step_only),
        "real_time_module_callers": reach,
        "step_bindings": bindings,
        "n_step_bindings": len(bindings),
        "distinct_minutes_per_step": distinct,
        "bindings_reaching_production": len(in_prod),
        "production_minutes_per_step": prod_mps,
        "production": _production_steps(),
        "orphan_timescale_fields": _orphan_timescales(),
        "p3": _p3_threshold(),
    }


def _span_hours(d):
    """Hours the production loop spans, if a production binding exists."""
    prod, mps = d.get("production"), d.get("production_minutes_per_step")
    if not prod or not mps:
        return None
    return prod["steps"] * min(mps) / 60


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
        mods = sorted({x["module"] for x in b})
        L += [f"**{len(b)} declarations bind a step to a real duration**, in "
              f"{', '.join(f'`{m}`' for m in mods)}, out of {n} library "
              f"modules. The core biochemical loop states none; these do.", ""]
        if len(d["distinct_minutes_per_step"]) > 1:
            vals = ", ".join(f"{v:g} min" for v in d["distinct_minutes_per_step"])
            hi, lo = max(d["distinct_minutes_per_step"]), min(d["distinct_minutes_per_step"])
            L += [f"**They disagree.** The distinct step durations declared are "
                  f"{vals} -- a factor of {hi/lo:.0f} apart, unreconciled.", ""]

    if b:
        L += ["| module | line | how it is written | one step = | reached by |",
              "|---|--:|---|--:|---|"]
        for x in b:
            mps = f"{x['minutes_per_step']:g} min" if x["minutes_per_step"] is not None else "unstated"
            who = ", ".join(f"`{c}`" for c in x["callers"]) or "*no caller*"
            L.append(f"| `{x['module']}` | {x['line']} | `{x['text'][:64]}` | "
                     f"{mps} | {who} |")
        L += [""]

    span = _span_hours(d)
    prod = d.get("production")
    if span is not None and prod:
        L += ["## The binding that is load-bearing", ""]
        L += [f"`{prod['binary']}` runs **{prod['steps']} steps**, the same "
              f"count as the production matrix, against a module documented at "
              f"{min(d['production_minutes_per_step']):g} minute per step. So "
              f"the whole production run spans **{span:.1f} hours** -- a fact "
              f"stated nowhere, and the only step duration that reaches a "
              f"binary.", ""]
        p3 = d.get("p3") or {}
        if "threshold_hours" in p3:
            t = p3["threshold_hours"]
            L += [f"`PREREGISTRATION.md` sets P3's falsification threshold at "
                  f"**{t:g} hours**. The production loop spans {span:.1f} "
                  f"hours, which is **{t/span:.0f}x shorter than the threshold "
                  f"it would be scored against**.", ""]
            if "window_hours" in p3:
                lo, hi = p3["window_hours"]
                L += [f"P3's predicted window is {lo/24:g} to {hi/24:g} days "
                      f"({lo:g}-{hi:g} h), i.e. {lo/span:.0f}x to {hi/span:.0f}x "
                      f"the entire run. The model cannot represent either "
                      f"outcome of its own most directly testable prediction -- "
                      f"not because the units are missing, but because the one "
                      f"in force is inconsistent with the prediction's.", ""]

    L += ["## Which modules carry real time, and whether anything calls them", ""]
    L += ["| module | callers |", "|---|---|"]
    for m, c in sorted(d["real_time_module_callers"].items()):
        L.append(f"| `{m}` | {', '.join(f'`{x}`' for x in c) or '**none**'} |")
    L += [""]
    orphan_f = {m: [f for f, r in fs.items() if not r]
                for m, fs in d.get("orphan_timescale_fields", {}).items()}
    orphan_f = {m: fs for m, fs in orphan_f.items() if fs}
    if orphan_f:
        tot = sum(len(v) for v in orphan_f.values())
        L += [f"Module-level callers are too coarse to check what is actually "
              f"composed. **{tot} real-time configuration fields are referenced "
              f"nowhere outside their own module**:", ""]
        for m, fs in sorted(orphan_f.items()):
            L.append(f"* `{m}`: {', '.join(f'`{f}`' for f in fs)}")
        L += [""]
        L += ["An earlier version of this report said the day-scale, hour-scale "
              "and minute-scale modules were \"all composed into the same "
              "180-step run\". That was asserted rather than measured, and it "
              "is false for every field listed here: the module is reachable, "
              "the timescale inside it is not.", ""]

    L += ["## What this does not do", ""]
    L += ["* It does not choose a step duration, and does not reconcile the "
          "declarations that disagree. That moves every calibrated layer and "
          "the committed byte-identity gates.",
          "* It does not claim the core loop declares one. It does not; the "
          "modules composed alongside it do, which is a different statement "
          "and the one that took three versions to get right.",
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
    print(f"  step bindings found: {d['n_step_bindings']}  "
          f"distinct durations: {d['distinct_minutes_per_step']}")
    span = _span_hours(d)
    if span is not None:
        print(f"  production run spans {span:.1f} h")


if __name__ == "__main__":
    main()
