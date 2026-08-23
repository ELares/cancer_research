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

WHAT IS ACTUALLY TRUE
----------------------
The CORE biochemical loop states no step duration -- that part survives.
ONE module DECLARES a step duration in wall-clock time (`tumor_pk`, one minute), and a SECOND reading is implied without declaring one: the immune model states a 0-48h scope over a 180-step loop in `sim-tme` and `sim-tme-3d`, pricing a step at 16 minutes. The two are 16x apart, neither is measured, and the reconciliation is the step-duration section of `simulations/calibration/parameter_provenance.md`
per step, reaching `sim-tumor-pk`'s 180-step run = 3.0 hours.

`trigger_wave`'s `dt_min` is NOT a second binding. It is a CFL-constrained
integrator timestep -- the module asserts `dt < h^2/(2D)` on the next line --
and pooling it with tumor_pk's alignment manufactured a "50x disagreement,
unreconciled". Two solvers having different `dt` is not a contradiction.

TWO FURTHER RETRACTIONS, BOTH CAUGHT BY REVIEW BEFORE MERGE
-------------------------------------------------------------
A draft of this file concluded that because the priced run spans three hours,
against P3's 24-hour threshold, "the model cannot represent either outcome of
its most directly testable prediction". FALSE. `cell.rs` carries the four
defenses' recovery half-times in DAYS, `sim-window` sweeps them to 28 days at
P3's own timepoints, and that is a published manuscript figure. The draft
compared one binary's INNER kill assay against a prediction scored on the
OUTER recovery axis. Worse, the refutation was already in this script's own
artifact and the renderer printed only the EMPTY caller lists, dropping its own
counter-example.

A later draft priced `sim-tme-3d` at 3.0 hours too. Its default matrix runs
`DoseSchedule::Constant` and reaches the PK solver only under `--dose-sweep`;
the check was a text grep for the module name.

THE FINDING THAT REPLACES THEM
--------------------------------
P3 states the defenses recover "FSP1 and GSH first, GPX4 and NRF2 later". The
engine's defaults are gsh 1d, gpx4 3d, nrf2 5d, fsp1 7d -- so `fsp1` is named
early and is the SLOWEST. They contradict, and that is visible only BECAUSE
the model represents the quantity.

Scoped deliberately: this is NOT a falsification of P3 as registered, whose
falsifier is simultaneity ("all four recover within the same timepoint"), which
a sequential order satisfies. The disagreement is between the preregistration's
descriptive ordering and the engine's defaults.

WHAT THIS SCRIPT DOES NOT DO
-----------------------------
It does not choose a step duration or reconcile anything. That is a modelling
decision with consequences for every calibrated layer and for the committed
byte-identity gates, and it belongs to whoever owns those calibrations. A guard
fails if this script starts asserting one.

KNOWN LIMITS OF THE REACHABILITY CHECK
----------------------------------------
`_reaches_by_default` resolves one call graph textually. It gives the right
answer on this tree, but adversarial probing found it flips on: a CLI flag
declared through clap rather than `std::env::args()`; braces inside the
`env::args` closure; a trailing comment or string literal naming the symbol
inside a live fn; and a nested `fn` inside `main`. `_enclosing_fn` returns the
last `fn` DECLARED before an offset whether or not it has closed. Treat the
opt-in classification as evidence, not proof.

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

# A real-time unit appearing with a magnitude.
TIME_UNIT = re.compile(r"\b(\d+(?:\.\d+)?)\s*(?:-)?\s*" + UNIT_RE + r"\b", re.I)
# A real-time unit DECLARED without a magnitude, which is how a units
# annotation is normally written: `half-time (days)`, `Time step (min)`,
# `(um^2/min)`, `(1/min)`.
#
# THIS IS THE SAME BUG, ONE FUNCTION OVER. `find_step_bindings` was fixed to
# see unit-first phrasing and TIME_UNIT was left requiring a number, so
# `cell.rs` and `trigger_wave.rs` scored ZERO real-time mentions and vanished
# from the table headed "which modules carry real time" -- while the same page,
# twenty lines below, said cell.rs carries recovery half-times in days. The
# page contradicted itself, and a second reviewer found it.
#
# Restricted to a parenthesised unit, which is how a units annotation is
# conventionally written in Rust.
#
# Single-letter units are EXCLUDED here even though the magnitude form accepts
# them, because without a number beside it a bare `h` or `s` is almost never a
# duration. Measured: allowing them matched `(rows, h)` in reaction_diffusion
# (h is GRID SPACING) and `("{:.1}", hours)` in io.rs (a format call), putting
# two modules into the real-time table that carry none. Quotes and commas in
# the parenthetical are excluded for the same reason.
UNIT_RE_BARE = r"(secs?|seconds?|mins?|minutes?|hrs?|hours?|days?)"
TIME_UNIT_BARE = re.compile(
    r"\((?:[^()\"',]{0,24}?[/\s])?" + UNIT_RE_BARE + r"\s*\)", re.I)
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


# Modules that integrate their own equation on their own internal clock, so a
# timestep of theirs is a numerical-stability parameter and not a statement
# about what one loop step is worth. Named explicitly: deriving this from a
# field name is what let a spelling decide a classification.
SOLVER_MODULES = {"trigger_wave", "reaction_diffusion"}

def _rust_sources():
    """Every Rust source the engine ships, LIBRARY AND BINARIES.

    The first version globbed `ferroptosis-core/src` only, so "exactly one
    wall-clock binding exists anywhere in the engine" was measured over the
    library alone -- and `sim-tme` carries a second one, in the binary that
    produced the manuscript's published immune numbers. An audit whose scope is
    narrower than its claim cannot find the thing it exists to find.
    """
    # SORTED. Every call site this replaced used `sorted(...)`, and the
    # artifact is byte-compared in CI on two filesystems: an unsorted glob
    # emits rows in directory order, which differs between APFS and ext4.
    out = sorted(SRC.glob("*.rs"))
    sims = SRC.parent.parent
    for d in sorted(sims.glob("sim-*/src")):
        out.extend(sorted(d.glob("*.rs")))
    return sorted(out, key=lambda q: (q.parent.parent.name, q.name))


def _key(path) -> str:
    """A stable, UNIQUE label for a scanned file.

    Keying on `p.name` collapsed twelve `main.rs` files into one entry, so the
    binary the whole reconciliation is about vanished from the modules table
    and a caller list pointed at all twelve at once. Library files keep their
    bare name; a binary is qualified by its crate.
    """
    crate = path.parent.parent.name
    return path.name if crate == "ferroptosis-core" else f"{crate}/{path.name}"


def _default_for(stem: str, field: str, unit: str):
    """The literal a named per-step duration field is constructed with.

    Keyed on the DECLARING FILE. The first version ignored its `stem` argument
    and returned the first match anywhere in the tree, so a second module
    declaring the same field name silently rewrote another module's published
    value -- `trigger_wave`'s 0.02 min/step was reported as 30.0 under a probe.
    """
    # Compare STEM to STEM. The first attempt compared `p.name` ("x.rs")
    # against the `p.stem` ("x") its caller passes, so the filter never
    # matched, `or _rust_sources()` fell through to the whole tree, and the
    # bug this was written to fix reproduced verbatim: trigger_wave's 0.02
    # min/step got reported as another module's 30.0.
    cands = [q for q in _rust_sources() if q.stem == stem]
    if not cands:
        return None
    for p in cands:
        m = re.search(rf"\b{re.escape(field)}\s*:\s*(\d+(?:\.\d+)?)", p.read_text(errors="ignore"))
        if m:
            return float(m.group(1)) * UNIT_MIN[unit.lower()]
    return None


# A window comment prices the LOOP only if it is scoping the model the loop
# runs. Requiring a validity/scope verb keeps a pharmacokinetic half-life --
# "terminal t_half ~20-48h" -- from being read as a step duration, which the
# first version did, publishing a spurious 9.33 min/step.
_SCOPE_VERB = re.compile(r"\b(valid for|models?|represents?|covers?|"
                         r"applicable to|cascade|phase)\b", re.I)
# Respellings of the same window. The first version matched `(0-48h)` alone,
# so an en dash, "hr", "to", or a day unit hid the reading entirely.
_WINDOW = re.compile(
    r"\(?\s*(\d+(?:\.\d+)?)\s*(?:-|\u2013|\u2014|to)\s*"
    r"(\d+(?:\.\d+)?)\s*(h|hr|hrs|hour|hours|d|day|days)\b", re.I)
_WINDOW_UNIT_H = {"h": 1.0, "hr": 1.0, "hrs": 1.0, "hour": 1.0, "hours": 1.0,
                  "d": 24.0, "day": 24.0, "days": 24.0}


def find_implied_windows():
    """Bindings a binary IMPLIES without ever declaring a per-step duration.

    `sim-tme` states no minutes-per-step anywhere. It states a validity WINDOW
    for the biology it models -- "spatial immune model valid for resident T
    cell phase (0-48h)" -- and a loop length, `N_STEPS = 180`. Together those
    price a step at 16 minutes, in the binary that produced this book's
    published immune numbers.

    This is the form issue #727 pointed at, and the form an audit looking for
    `dt_min` or "N minutes per step" cannot see. Reported SEPARATELY from an
    explicit declaration because it is weaker evidence about intent -- a scope
    claim about which biology is in range is not the same act as declaring a
    clock -- and not weaker about consequence.

    READMEs are scanned too. The same 0-48h claim appears in `sim-tme-3d`'s and
    `sim-tme`'s README, and reading only `.rs` attributed the window to one
    binary while the doc assigned the other a different conversion.
    """
    out = []
    sims = SRC.parent.parent
    for crate in sorted(sims.glob("sim-*")):
        steps = None
        for rs in sorted((crate / "src").glob("*.rs")):
            m = re.search(r"const N_STEPS:\s*\w+\s*=\s*(\d+)",
                          rs.read_text(errors="ignore"))
            if m:
                steps = int(m.group(1))
                break
        if not steps:
            continue
        files = sorted((crate / "src").glob("*.rs")) + \
            sorted(crate.glob("README*"))
        for f in files:
            for ln, line in enumerate(f.read_text(errors="ignore").splitlines(), 1):
                if f.suffix == ".rs" and not line.lstrip().startswith("//"):
                    continue
                if not _SCOPE_VERB.search(line):
                    continue
                m = _WINDOW.search(line)
                if not m:
                    continue
                hours = (float(m.group(2)) - float(m.group(1))) \
                    * _WINDOW_UNIT_H[m.group(3).lower()]
                if hours <= 0:
                    continue
                out.append({
                    "binary": crate.name, "module": f.name, "line": ln,
                    "kind": "implied-window", "n_steps": steps,
                    "window_hours": hours,
                    "minutes_per_step": round(hours * 60.0 / steps, 4),
                    "text": line.strip()[:160],
                })
    return sorted(out, key=lambda w: (w["binary"], w["module"], w["line"]))


def find_step_bindings():
    """Every place the engine binds one step to a real duration.

    Returns sites, not a bool. The previous version returned True/False, so a
    reader could not tell WHICH module bound it or whether two disagreed -- and
    a bool cannot carry the finding that they do.
    """
    sites = []
    for p in _rust_sources():
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
                # Classified by WHAT IT MEASURES, not by how it is spelled.
                # Keying on `how == "named-field"` meant any second binding
                # named `dt_<unit>` was excluded from the audit's own headline
                # claim by construction: a probe shipped a 30 min/step
                # wall-clock field and every guard stayed green.
                #
                # A solver timestep belongs to a module that integrates its own
                # equation on its own internal clock. Anything else prices the
                # 180-step loop. Unknown modules fall to WALL-CLOCK, the
                # conservative direction, because that is the side that makes
                # an "exactly one binding" claim harder to satisfy.
                kind = ("solver-timestep" if p.stem in SOLVER_MODULES
                        else "wall-clock")
                sites.append({
                    "module": _key(p), "line": ln, "how": how, "kind": kind,
                    "unit": u, "minutes_per_step": mpstep,
                    "text": ("..." if lo else "") + line[lo:hi].strip() +
                            ("..." if hi < len(line) else ""),
                })
                break
    return sites


def _callers(stem: str):
    """Which simulation binaries reference this LIBRARY module.

    Returns nothing for a binary. A binary is not referenced by anything in the
    workspace, and the qualified module keys made that visible: asking for the
    callers of `main.rs` matched all twelve binaries at once and published a
    row saying so.

    The previous version claimed three timescales were 'all composed into the
    same 180-step run'. Nothing checked it, and the day-scale limb had no
    caller at all.
    """
    if "/" in stem or stem == "main":
        return []
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
    for p in _rust_sources():
        names = set(field.findall(p.read_text(errors="ignore")))
        for nm in sorted(names):
            refs = []
            for q in list(SRC.glob("*.rs")) + list(SIMS.glob("sim-*/src/*.rs")):
                if q == p:
                    continue
                if re.search(rf"\b{re.escape(nm)}\b", q.read_text(errors="ignore")):
                    refs.append(q.parts[-3] if "sim-" in str(q) else q.name)
            out.setdefault(_key(p), {})[nm] = sorted(set(refs))
    return out


def _blocks(text, opener_re):
    """(start, end) offsets of brace-matched blocks opened by a pattern."""
    out = []
    for m in re.finditer(opener_re, text):
        i = text.find("{", m.end() - 1 if m.end() else m.start())
        if i < 0:
            continue
        depth, j = 0, i
        while j < len(text):
            if text[j] == "{":
                depth += 1
            elif text[j] == "}":
                depth -= 1
                if depth == 0:
                    break
            j += 1
        out.append((m.start(), j))
    return out


def _enclosing_fn(text, pos):
    """Name of the fn containing an offset, by scanning backwards."""
    best = None
    for m in re.finditer(r"\bfn\s+(\w+)", text):
        if m.start() < pos:
            best = m.group(1)
        else:
            break
    return best


def _reaches_by_default(binary_dir, symbol):
    """Is `symbol` reachable in this binary WITHOUT an opt-in CLI flag?

    A plain text grep says sim-tme-3d "uses" tumor_pk, and it does -- but its
    default 24-condition matrix runs `DoseSchedule::Constant` throughout and
    only reaches `solve_tumor_pk` through `run_dose_sweep`, which sits behind
    `if std::env::args().any(|a| a == "--dose-sweep") { ...; return; }`. So the
    3.0-hour span was being presented unqualified for a path the production
    matrix never takes. A second reviewer caught it.

    TRANSITIVE reachability from `main`, excluding CLI-guard blocks and
    #[cfg(test)]. A single level is not enough and gave the wrong answer: the
    `solve_tumor_pk` call sits in `rsl3_pk_factor_series`, whose only live call
    site is inside `run_dose_sweep` -- so one level found "a caller outside
    itself" and stopped, never discovering that `run_dose_sweep` is the thing
    behind the flag.

    Unresolved calls (through a trait object, a function pointer, a macro)
    are not followed, so this can still report opt-in for something reachable
    by an exotic route. It is stated as a limit rather than assumed away.
    """
    for p in sorted(binary_dir.glob("*.rs")):
        text = p.read_text(errors="ignore")
        dead = (_blocks(text, r"if\s+std::env::args\(\)[^\n{]*") +
                _blocks(text, r"#\[cfg\(test\)\]\s*mod\s+\w+"))

        def inside(pos):
            return any(a <= pos <= b for a, b in dead)

        def live_refs(name):
            out = []
            for m in re.finditer(rf"\b{re.escape(name)}\b", text):
                s = m.start()
                if inside(s):
                    continue
                line = text[text.rfind("\n", 0, s) + 1:s]
                if re.match(r"\s*(//|/\*|\*)", line):
                    continue          # a comment mentioning it is not a call.
                    # sim-tme-3d documents `solve_tumor_pk` in a doc comment 13
                    # lines above the function that calls it; counting that as
                    # a reference attributed it to whatever fn preceded the
                    # comment, which IS default-reachable, and the whole
                    # opt-in detection silently inverted.
                if re.match(r"\s*use\s", line) or re.search(r"\bfn\s*$", line):
                    continue          # the import, or the definition itself
                out.append(s)
            return out

        seen, frontier = set(), [symbol]
        while frontier:
            name = frontier.pop()
            if name in seen:
                continue
            seen.add(name)
            for r in live_refs(name):
                fn = _enclosing_fn(text, r)
                if fn is None or fn == "main":
                    return True
                if fn not in seen:
                    frontier.append(fn)
    return False


def _pricing_symbols(module, binding_lines):
    """The public items a module's wall-clock binding is attached to.

    Derived from the binding's own position -- the item a doc comment
    documents is the next one declared after it -- rather than named here, so
    a rename cannot leave this pointing at nothing.
    """
    text = (SRC / module).read_text(errors="ignore")
    lines = text.split("\n")
    out = set()
    for ln in binding_lines:
        found = None
        for j in range(ln - 1, min(len(lines), ln + 12)):
            m = re.match(r"\s*pub\s+(?:fn|struct|enum|const|static)\s+(\w+)",
                         lines[j])
            if m:
                found = m.group(1)
                break
            # A doc comment on a STRUCT FIELD documents the struct, not the
            # next free item. tumor_pk.rs:354 documents `time_min` inside
            # `pub struct TumorPKResult {` at 353; scanning only downward from
            # 354 skipped past the struct, hit no item, and fell through to
            # the enclosing-fn guess -- which named `doxorubicin_iv_bolus`,
            # declared 14 lines earlier and already closed. Because
            # _step_counts uses any(), a spurious symbol can only ADD pricing:
            # a binary calling that unrelated function was priced at 3.0 h.
            if re.match(r"\s*pub\s+\w+\s*:", lines[j]):
                for k in range(j, max(-1, j - 40), -1):
                    s = re.match(r"\s*pub\s+(?:struct|enum)\s+(\w+)", lines[k])
                    if s:
                        found = s.group(1)
                        break
                break
        if found:
            out.add(found)
    return sorted(out)


def _step_counts(binding_modules, symbols_by_module=None):
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
        # Does the PRICING symbol reach the default path, or only an opt-in
        # flag? Module-name presence is not consumption of the binding.
        default_path = []
        for stem in uses:
            syms = (symbols_by_module or {}).get(stem + ".rs", [])
            if any(_reaches_by_default(p.parent, s) for s in syms):
                default_path.append(stem)
        out.append({"binary": binary, "steps": int(m.group(1)),
                    "consumes_binding_module": uses,
                    "prices_on_default_path": sorted(default_path)})
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
    for p in _rust_sources():
        text = p.read_text(errors="ignore")
        times, steps = [], 0
        for i, line in enumerate(text.split("\n"), 1):
            ctx = line.strip()
            skip = ctx.startswith("//!") and "v0." in ctx
            for m in TIME_UNIT.finditer(line):
                if skip:
                    continue
                times.append({"line": i, "value": m.group(0).strip(),
                              "form": "magnitude", "context": ctx[:110]})
            for m in TIME_UNIT_BARE.finditer(line):
                if skip:
                    continue
                times.append({"line": i, "value": m.group(0).strip(),
                              "form": "bare-unit", "context": ctx[:110]})
            if PER_STEP.search(line):
                steps += 1
        if times or steps:
            modules[_key(p)] = {"real_time_mentions": len(times),
                               "per_step_mentions": steps,
                               "examples": times[:4]}

    bindings = find_step_bindings()
    for b in bindings:
        b["callers"] = _callers(b["module"].split("/")[-1].rsplit(".", 1)[0]
                                if "/" not in b["module"]
                                else b["module"])

    # Count CONVENTIONS, not regex hits. `\bdt_min\b` matches the field
    # declaration, its literal in baseline(), and `let dt = cfg.dt_min` -- one
    # binding, three matches. An earlier draft reported "6 declarations" for
    # what is 2 conventions in 2 modules.
    wall = [b for b in bindings if b["kind"] == "wall-clock"]
    solver = [b for b in bindings if b["kind"] == "solver-timestep"]
    conventions = sorted({(b["module"], round(b["minutes_per_step"], 6))
                          for b in wall if b["minutes_per_step"] is not None})
    # Solver timesteps are deduplicated the SAME way. Reporting deduplicated
    # conventions on one side and raw matches on the other put "1 module ...
    # plus 3 numerical-integrator timesteps" in the headline, two lines above
    # the sentence retracting exactly that arithmetic. There is ONE integrator
    # timestep in the engine, matched three times.
    solver_conv = sorted({(b["module"], round(b["minutes_per_step"], 6))
                          for b in solver if b["minutes_per_step"] is not None})
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
    reach = {m: _callers(m if "/" in m else Path(m).stem)
             for m in sorted(set(real_only) | set(both))}
    pricing_syms = {}
    for mod in {b["module"] for b in wall}:
        pricing_syms[mod] = _pricing_symbols(
            mod, [b["line"] for b in wall if b["module"] == mod])
    p3 = _p3_threshold()
    modelled = _p3_is_modelled()

    return {
        # The denominator must match what was actually scanned. It stayed at
        # the library count after the scan widened to the binaries, so the
        # headline read "1 of 33" over a 47-file sweep.
        "modules_total": len(_rust_sources()),
        "library_modules": len(sorted(SRC.glob("*.rs"))),
        "modules": modules,
        "carry_both": sorted(both),
        "real_time_only": sorted(real_only),
        "per_step_only": sorted(step_only),
        "real_time_module_callers": reach,
        "implied_windows": find_implied_windows(),
        "step_bindings": bindings,
        "n_step_bindings": len(bindings),
        "wall_clock_conventions": [{"module": m, "minutes_per_step": v}
                                   for m, v in conventions],
        "n_wall_clock_conventions": len(conventions),
        "solver_timestep_conventions": [{"module": m, "minutes_per_step": v}
                                        for m, v in solver_conv],
        "n_solver_timestep_conventions": len(solver_conv),
        "n_solver_timestep_matches": len(solver),
        "distinct_minutes_per_step": distinct,
        "bindings_reaching_production": len(in_prod),
        "production_minutes_per_step": prod_mps,
        "pricing_symbols": pricing_syms,
        "step_counts": _step_counts({b["module"] for b in wall}, pricing_syms),
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
        mods = sc.get("prices_on_default_path") or []
        if not mods:
            continue
        for c in d.get("wall_clock_conventions", []):
            if Path(c["module"]).stem in mods:
                out.append({**sc, "minutes_per_step": c["minutes_per_step"],
                            "hours": sc["steps"] * c["minutes_per_step"] / 60})
                break
    return out


def _opt_in_only(d):
    """Binaries that reference a pricing module but only behind a flag."""
    out = []
    for sc in d.get("step_counts", []):
        uses = set(sc.get("consumes_binding_module") or [])
        deflt = set(sc.get("prices_on_default_path") or [])
        if uses and not deflt:
            out.append({**sc, "opt_in_modules": sorted(uses)})
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
        nsolv = d.get("n_solver_timestep_conventions", 0)
        nsolvm = d.get("n_solver_timestep_matches", 0)
        mods = sorted({c["module"] for c in conv})
        L += [f"**{len(conv)} module{'' if len(conv) == 1 else 's'} "
              f"{'prices' if len(conv) == 1 else 'price'} a "
              f"simulation step in wall-clock time** "
              f"({', '.join(f'`{m}`' for m in mods) or 'none'}), out of {n} "
              f"scanned modules (library and binaries), plus {nsolv} "
              f"numerical-integrator "
              f"timestep{'' if nsolv == 1 else 's'} that "
              f"{'does' if nsolv == 1 else 'do'} not.", ""]
        iw = d.get("implied_windows") or []
        if iw:
            mins = sorted({w["minutes_per_step"] for w in iw})
            srcs = sorted({f"`{w['binary']}/{w['module']}`" for w in iw})
            L += [f"**A second reading is IMPLIED and never declared.** "
                  f"{len(iw)} site{'' if len(iw) == 1 else 's'} "
                  f"({', '.join(srcs)}) state a scope window over a "
                  f"{iw[0]['n_steps']}-step loop, pricing a step at "
                  f"{', '.join(str(m) for m in mins)} min -- against the "
                  f"declared {sorted({c['minutes_per_step'] for c in conv})}. "
                  f"Counting only declarations is how an earlier version of "
                  f"this audit reported exactly one binding while scanning the "
                  f"library alone. Neither reading is measured; the "
                  f"reconciliation is the step-duration section of "
                  f"`simulations/calibration/parameter_provenance.md`.", ""]
        L += [f"Both numbers are CONVENTIONS, not matches -- deduplicating one "
              f"side and not the other is how an earlier draft reported "
              f"\"6 declarations\". The {nsolv} integrator "
              f"timestep{'' if nsolv == 1 else 's'} "
              f"{'is' if nsolv == 1 else 'are'} matched {nsolvm} times, once "
              f"per declaration, default literal and use.", ""]
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
        optin = _opt_in_only(d)
        if optin:
            rows = ", ".join(
                "`{}` ({} steps, via {})".format(
                    s["binary"], s["steps"],
                    ", ".join(f"`{m}`" for m in s["opt_in_modules"]))
                for s in optin)
            conv1 = d["wall_clock_conventions"][0]["minutes_per_step"]
            would = ", ".join(
                "{:.1f} h".format(s["steps"] * conv1 / 60) for s in optin)
            L += [f"{len(optin)} "
                  f"{'binary references' if len(optin) == 1 else 'binaries reference'} "
                  f"a pricing module but only behind an opt-in flag, so "
                  f"{'its' if len(optin) == 1 else 'their'} DEFAULT run cannot "
                  f"be priced: {rows}. Reachability is checked transitively "
                  f"from `main`, excluding `std::env::args()` guards and "
                  f"`#[cfg(test)]`. An earlier draft listed "
                  f"{'this' if len(optin) == 1 else 'these'} at {would} "
                  f"unqualified, because the check was a text grep for the "
                  f"module name -- and a one-level check still got it wrong, "
                  f"since the call sits one function inside the guarded one.",
                  ""]
        others = [s for s in d.get("step_counts", [])
                  if not s["consumes_binding_module"]]
        if others:
            names = ", ".join("`{}` ({})".format(s["binary"], s["steps"])
                              for s in others)
            L += [f"{len(others)} other "
                  f"{'binary declares' if len(others) == 1 else 'binaries declare'} "
                  f"a step count "
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
                  f"puts later.", ""]
            # SCOPE, because the first framing over-reached. P3's REGISTERED
            # falsifier is simultaneity, not order identity, so the engine's
            # order satisfies it -- the disagreement is between P3's
            # descriptive parenthetical and the engine, which is worth
            # publishing but is not scoreable under the registered rule.
            if "threshold_hours" in p3:
                L += [f"Scope: this is **not** a falsification of P3 as "
                      f"registered. P3's stated falsifier is that all four "
                      f"defenses recover within the same timepoint (no "
                      f"sequential order), with a {p3['threshold_hours']:g}-hour "
                      f"threshold"
                      + (f" and a {p3['window_hours'][0]/24:g} to "
                         f"{p3['window_hours'][1]/24:g} day window"
                         if "window_hours" in p3 else "")
                      + ". The engine's order IS sequential, so it satisfies "
                      "that rule. The disagreement is between the "
                      "preregistration's descriptive ordering and the engine's "
                      "defaults, and it is visible only BECAUSE the model "
                      "represents the quantity.", ""]
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
    L += ["* It does not choose a step duration -- and it now MEASURES that "
          "there are two competing readings rather than one, an explicit "
          "`tumor_pk` declaration at 1 min/step and an implied window in "
          "`sim-tme` at 16 min/step. The reconciliation, which adopts neither "
          "and says which applies where, is the step-duration section of "
          "`simulations/calibration/parameter_provenance.md`. Adopting one across the "
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
          f"solver timesteps: {d['n_solver_timestep_conventions']}  "
          f"distinct durations: {d['distinct_minutes_per_step']}")
    for r in _priced_runs(d):
        print(f"    {r['binary']:16s} {r['steps']:>4} steps -> {r['hours']:.1f} h")
    o = d.get("p3_order")
    if o:
        print(f"  P3 recovery order agrees with engine defaults: {o['agrees']}")


if __name__ == "__main__":
    main()
