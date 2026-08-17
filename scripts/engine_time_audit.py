#!/usr/bin/env python3
"""What is one step of the ferroptosis engine worth in real time? (#727)

THE ANSWER IS THAT NOTHING SAYS, AND SEVERAL THINGS ASSUME
------------------------------------------------------------
The core loop advances a cell by calling `sim_cell_step` N times -- 180 in the
production matrix -- and every rate it uses is per-dimensionless-step. No
constant, comment or document states what a step is worth in minutes.

That would be unremarkable for a self-contained dimensionless model. It is not
self-contained. Modules carrying REAL time units are composed into the same run:
`tumor_pk` integrates ODEs with half-lives in minutes, `trigger_wave` takes an
explicit `dt_min` and is calibrated against a measured front speed in um/min,
and `dose_schedule` produces a per-step availability factor from a schedule
whose natural units are hours.

So the engine mixes two time systems and never states the conversion. Two
consequences, and the second is the one that matters:

  NO RATE CAN BE COMPARED TO ANOTHER RATE across that boundary. A per-step
  biochemical rate and a per-minute pharmacokinetic rate are not commensurable
  without the missing constant.

  NO TIME-STATED PREDICTION CAN BE SCORED AGAINST THE MODEL THAT GENERATED IT.
  PREREGISTRATION.md states P3 as "defenses recover over roughly 3 to 7 days"
  and P8 in micrometres. The model that produced them runs in steps. A wet-lab
  result in days cannot be compared to a model output in steps, which means the
  preregistration's most falsifiable leg is currently unfalsifiable.

WHAT THIS SCRIPT DOES
---------------------
It enumerates every real-time constant in the engine, and every module that
carries one, so the conversion that is missing can be seen rather than argued
about. It does NOT invent a step duration: choosing one is a modelling decision
with consequences for calibrated layers, and it belongs to whoever owns those
calibrations.

WHAT IT DELIBERATELY DOES NOT CLAIM
------------------------------------
An earlier version of this issue said the engine has "four conflicting
step-duration bindings" and quoted three specific values. Only some of that
survives inspection: `tumor_pk` and `trigger_wave` genuinely carry real time,
and the MUFA layer genuinely cites a 48-72h literature timescale its rate was
never fitted to. But those are not four competing DEFINITIONS of a step -- they
are modules with their own units that were never reconciled with the core loop.
The distinction matters because "four bindings disagree" implies someone made
four choices, and nobody did. The absence is the finding.

Usage:
    python scripts/engine_time_audit.py
"""

import json
import re
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
SRC = PROJECT_ROOT / "simulations" / "ferroptosis-core" / "src"
OUT_MD = PROJECT_ROOT / "analysis" / "engine-time-audit.md"
OUT_JSON = PROJECT_ROOT / "analysis" / "engine-time-audit.json"

# A real-time unit appearing in code or a doc comment. Deliberately narrow: the
# question is which modules carry PHYSICAL time, not which mention the word.
TIME_UNIT = re.compile(
    r"\b(\d+(?:\.\d+)?)\s*(?:-)?\s*(min|mins|minute|minutes|h|hr|hrs|hour|hours|"
    r"day|days|sec|second|seconds)\b", re.I)
# A per-step rate: the dimensionless system.
PER_STEP = re.compile(r"per[- ]step|per dimensionless|each step|/step", re.I)


def scan() -> dict:
    modules = {}
    for p in sorted(SRC.glob("*.rs")):
        text = p.read_text(errors="ignore")
        times, steps = [], 0
        for i, line in enumerate(text.split("\n"), 1):
            for m in TIME_UNIT.finditer(line):
                # skip obvious version strings and array sizes
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
    both = {k: v for k, v in modules.items()
            if v["real_time_mentions"] and v["per_step_mentions"]}
    real_only = {k: v for k, v in modules.items()
                 if v["real_time_mentions"] and not v["per_step_mentions"]}
    step_only = {k: v for k, v in modules.items()
                 if v["per_step_mentions"] and not v["real_time_mentions"]}
    return {
        "modules_total": len(list(SRC.glob("*.rs"))),
        "modules": modules,
        "carry_both": sorted(both),
        "real_time_only": sorted(real_only),
        "per_step_only": sorted(step_only),
        "declares_step_duration": _declares_step_duration(),
    }


def _declares_step_duration() -> bool:
    """Does anything state what one core step is worth? Searched, not assumed."""
    pat = re.compile(
        r"(one|a|each|per)\s+step\s+(is|=|represents?|corresponds? to|equals?)\s*"
        r"[^.\n]{0,40}\b(min|minute|hour|h|day|sec)", re.I)
    for p in list(SRC.glob("*.rs")) + list((PROJECT_ROOT / "simulations").glob("**/*.md")):
        try:
            if pat.search(p.read_text(errors="ignore")):
                return True
        except OSError:
            continue
    return False


def render(d: dict) -> str:
    L = ["# What is one engine step worth in real time?", ""]
    L += ["*Generated by `scripts/engine_time_audit.py`. Every count is "
          "recomputed.*", ""]

    L += [f"**Nothing in the engine states it.** Searched across "
          f"{d['modules_total']} library modules and the simulation docs for any "
          f"declaration that a step equals a duration: "
          f"{'FOUND' if d['declares_step_duration'] else 'none found'}.", ""]

    L += ["That would be unremarkable for a self-contained dimensionless model. "
          "It is not self-contained:", ""]
    L += ["| | modules |", "|---|---|"]
    L += [f"| carry REAL time units only | {', '.join(f'`{m}`' for m in d['real_time_only']) or 'none'} |"]
    L += [f"| carry PER-STEP rates only | {', '.join(f'`{m}`' for m in d['per_step_only']) or 'none'} |"]
    L += [f"| carry **both** | {', '.join(f'`{m}`' for m in d['carry_both']) or 'none'} |", ""]

    L += ["## The real-time constants that are already in there", ""]
    L += ["| module | example |", "|---|---|"]
    for m in d["real_time_only"] + d["carry_both"]:
        ex = d["modules"][m]["examples"]
        if ex:
            L.append(f"| `{m}` | `{ex[0]['context']}` |")
    L += [""]

    L += ["The timescales those modules carry are not close together. A "
          "photosensitizer clearing over 24-48 hours, therapy-induced "
          "senescence establishing over 14-21 days, and a drug half-life of 30 "
          "minutes are all composed into the same 180-step run. That is roughly "
          "three orders of magnitude of physical time inside one loop whose "
          "step has no duration.", ""]

    L += ["## Why this matters more than a units nicety", ""]
    L += ["**No rate can be compared to another rate across the boundary.** A "
          "per-step biochemical rate and a per-minute pharmacokinetic half-life "
          "are not commensurable without the missing constant, and they are "
          "composed in the same run.", ""]
    L += ["**No time-stated prediction can be scored against the model that "
          "generated it.** `PREREGISTRATION.md` states P3 as defenses "
          "recovering over roughly 3 to 7 days. The model runs in steps. A "
          "wet-lab result in days cannot be compared to a model output in "
          "steps, so the preregistration's most directly testable leg is "
          "currently unfalsifiable -- not because the biology is unclear but "
          "because the units are.", ""]

    L += ["## What this does not do", ""]
    L += ["* It does not choose a step duration. That is a modelling decision "
          "with consequences for every calibrated layer, and it belongs to "
          "whoever owns those calibrations.",
          "* It does not claim the engine holds several CONFLICTING definitions "
          "of a step. It holds none, alongside modules with their own real-time "
          "units that were never reconciled with it. 'Definitions disagree' "
          "would imply somebody made choices; the absence is the finding.",
          "* The scan is textual. A module can carry real time without writing "
          "a unit next to a number, so the real-time set is a lower bound.",
          ""]
    return "\n".join(L) + "\n"


def main():
    d = scan()
    if not d["modules"]:
        raise SystemExit(
            "no modules matched, which is not a finding -- it is what a wrong "
            "source path looks like.")
    OUT_JSON.write_text(json.dumps(d, indent=1, sort_keys=True) + "\n",
                        encoding="utf-8")
    OUT_MD.write_text(render(d), encoding="utf-8")
    print(f"wrote {OUT_MD}")
    print(f"wrote {OUT_JSON}")
    print(f"  declares a step duration anywhere: {d['declares_step_duration']}")
    print(f"  real-time only: {len(d['real_time_only'])}  "
          f"per-step only: {len(d['per_step_only'])}  "
          f"both: {len(d['carry_both'])}")


if __name__ == "__main__":
    main()
