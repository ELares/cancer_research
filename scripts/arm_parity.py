#!/usr/bin/env python3
"""How far each treatment arm is from the ferroptosis arm, axis by axis.

WHY THIS EXISTS
---------------
The standing criticism of this repository is that it is a ferroptosis project
wearing a cancer-research title. `analysis/scope-audit.md` measures what the
committed WORK is about, `analysis/modality-coverage.md` measures what the
engine can be ASKED, and `analysis/modality-module-depth.md` measures the size
of the gap IN AGGREGATE -- one ratio, roughly eleven to thirteen times, for all
the arms at once.

None of them answers the question a campaign to close the gap has to ask every
week: *which arm, and how far?* An aggregate ratio hides the difference between
an arm with a calibrated leg and no prose and an arm with three pages of prose
and ninety lines of code. This joins the existing artifacts per arm and reports
the distance on each axis separately, because the axes do not move together and
a single "parity score" would let a cheap axis pay for an expensive one.

WHAT PARITY MEANS HERE, AND WHAT IT CANNOT MEAN
-----------------------------------------------
Six axes, each measured somewhere else and joined here:

  engine      production code lines and public functions the arm owns
  tests       test functions that exercise the arm's own module
  calibration the verdict from `analysis/modality-calibration.md`
  book        manuscript words under headings that name the arm
  figures     FIGURES.yaml entries whose generator draws the arm
  predictions registered, falsifiable predictions in PREREGISTRATION.md

LINES OF CODE ARE THE WEAKEST AXIS AND ARE REPORTED FIRST ONLY BECAUSE THE
CRITICISM WAS ABOUT SIZE. A module can be large and wrong; an arm can reach
parity on every column here and still be a worse piece of science than the arm
it is measured against. What the table can support is the negative: an arm at
one twentieth of the ferroptosis engine on every axis is not "represented" in
any sense a reader would accept, and that is a claim about this project rather
than about the biology.

THE ATTRIBUTION RULE, WHICH IS THE PART THAT CAN BE WRONG
---------------------------------------------------------
Engine code is attributed by MODULE OWNERSHIP, reusing the DEDICATED/SHARED
split `modality_module_depth.py` already publishes rather than inventing a
second one -- two documents counting the same crate two ways is a defect this
repository has already had to fix. Shared machinery is credited to NO arm, so
every arm's engine column is a lower bound, and the shared total is printed
beside the table rather than divided up.

Book words are attributed by HEADING, not by keyword frequency: a section
counts for an arm when its own heading names the arm. That rule is narrow and
it is stated because it under-counts -- Chapter 6 discusses every arm under
headings that name none of them, and that whole chapter therefore lands in the
unattributed remainder, which is printed. A keyword rule over body text was
rejected: `analysis/scope-audit.md` records what happened the last time this
repository counted subject by vocabulary, which was that the count moved on a
rename and an empty file named for a therapy was filed as a therapy.

Census volume is NOT a parity axis. It is the size of the literature, printed
for context, and two of the arms here -- radiation and chemotherapy -- have no
taxonomy row at all, so their census cell reads `no row` rather than zero. A
zero there would say the field does not exist.
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
CORE = REPO / "simulations" / "ferroptosis-core" / "src"
SIMS = REPO / "simulations"
ANALYSIS = REPO / "analysis"
MANUSCRIPT = REPO / "article" / "drafts" / "v1.md"
PREREG = REPO / "PREREGISTRATION.md"
FIGURES = REPO / "FIGURES.yaml"
OUT_MD = ANALYSIS / "arm-parity.md"
OUT_JSON = ANALYSIS / "arm-parity.json"

# The arms, as the engine defines them, with everything needed to find each one
# in six different artifacts. `mechanism` is the taxonomy row, or None where the
# taxonomy has no row for the modality -- which is itself a finding and is
# rendered as such.
ARMS = [
    {"arm": "RSL3", "label": "Ferroptosis induction", "mechanism": None,
     "modules": "ENGINE", "heading": ["ferroptosis engine"],
     "topic": ["ferroptosis", "rsl3", "gpx4", "lipid perox"],
     "calibration": None,
     "note": "the comparator: every module the other arms are measured against"},
    {"arm": "SDT", "label": "Sonodynamic therapy", "mechanism": "sonodynamic",
     "modules": [], "heading": ["sonodynamic"], "topic": ["sonodynamic", "sdt"],
     "calibration": None,
     "note": "shares the exogenous-ROS path and the depth physics with PDT"},
    {"arm": "PDT", "label": "Photodynamic therapy", "mechanism": None,
     "modules": ["photosensitizer_pk"], "heading": ["photodynamic"],
     "topic": ["photodynamic", "pdt", "photosensitiser", "photosensitizer"],
     "calibration": None,
     "note": "the one non-ferroptosis arm with a module of its own from the start"},
    {"arm": "Radiation", "label": "Ionizing radiation", "mechanism": None,
     "modules": ["radiation"], "heading": ["radiation", "radiotherapy"],
     "topic": ["radiation", "radiotherapy", "linear-quadratic"],
     "calibration": "Radiation (DNA channel)",
     "note": "no taxonomy row; the largest literature of any arm here"},
    {"arm": "Immunotherapy", "label": "Checkpoint blockade",
     "mechanism": "immunotherapy", "modules": [],
     "heading": ["checkpoint", "immune coupling"],
     "topic": ["checkpoint", "pd-1", "pd-l1", "immune coupling", "immunotherapy"],
     "calibration": "Checkpoint blockade",
     "note": "reaches only through shared immune machinery: owns no module"},
    {"arm": "AdoptiveCell", "label": "Adoptive cell therapy (CAR-T)",
     "mechanism": "car-t", "modules": ["adoptive"], "heading": ["car-t", "adoptive"],
     "topic": ["car-t", "adoptive", "effector"],
     "calibration": "CAR-T (adoptive transfer)", "note": ""},
    {"arm": "OncolyticVirus", "label": "Oncolytic virus",
     "mechanism": "oncolytic-virus", "modules": ["oncolytic"],
     "heading": ["oncolytic"], "topic": ["oncolytic", "virus"],
     "calibration": "Oncolytic virus", "note": ""},
    {"arm": "Ablation", "label": "Thermal and electrical ablation",
     "mechanism": "hifu", "modules": ["ablation"],
     "heading": ["ablation", "hifu"], "topic": ["ablation", "hifu", "electroporation"],
     "calibration": "Ablation (thermal)",
     "note": "two mechanisms in the taxonomy (hifu, electrochemical-therapy) "
             "share one arm"},
    {"arm": "Chemotherapy", "label": "Cytotoxic chemotherapy",
     "mechanism": None, "modules": ["chemo"],
     "heading": ["chemotherapy", "cytotoxic"],
     "topic": ["chemotherapy", "cytotoxic", "cell cycle", "dose density"],
     "calibration": "Chemotherapy (cell-cycle)",
     "note": "no taxonomy row, and the only arm whose dose-response target is "
             "unreachable rather than merely unfitted"},
    {"arm": "AntibodyDrugConjugate", "label": "Antibody-drug conjugate",
     "mechanism": "antibody-drug-conjugate", "modules": ["adc"],
     "heading": ["antibody-drug conjugate", "adc"],
     "topic": ["antibody-drug", "adc", "bystander"],
     "calibration": "ADC bystander effect", "note": ""},
]

# Arms the engine does NOT have, listed because a parity table that shows only
# what exists reports the campaign's progress and hides its scope.
#
# Cytotoxic chemotherapy was the sharpest entry here and has left the list: it
# now has a module, a `Treatment` variant and a row in the panel. What it does
# NOT have is a fitted dose-response, and that is recorded as a calibration
# verdict in the table above rather than as an absence here -- the two are
# different claims and collapsing them would let a built-but-uncalibrated arm
# read as a finished one.
ABSENT_ARMS = [
    ("Targeted small-molecule therapy",
     "PARP synthetic lethality is inside `radiation`; there is no arm for "
     "kinase inhibition, and the taxonomy's synthetic-lethality row is served "
     "by a boost on radiation's alpha"),
    ("Hormone therapy",
     "no arm and no taxonomy row, against a literature the census can measure "
     "only through its drug descriptors"),
    ("Radioligand therapy",
     "a taxonomy row and a diagnostic-therapy chain, but no engine arm: the "
     "radiation module models external beam only"),
]


def _json(name):
    return json.loads((ANALYSIS / name).read_text())


def _rust_tests_for(module: str) -> int:
    """Test functions inside a module's own file.

    Counted in the file rather than by name, because a test named for an arm
    can live anywhere and a test in the arm's file is exercising the arm by
    construction."""
    path = CORE / f"{module}.rs"
    if not path.exists():
        return 0
    return len(re.findall(r"^\s*#\[test\]", path.read_text(), flags=re.MULTILINE))


def _engine_rust_tests() -> int:
    depth = _json("modality-module-depth.json")
    owned = {m["module"] for m in depth["dedicated"]} | {
        m["module"] for m in depth["shared"]}
    return sum(_rust_tests_for(p.stem) for p in CORE.glob("*.rs")
               if p.stem not in owned and p.stem != "lib")


def _sections(md: str):
    """(chapter, heading, words) per section, and (chapter, words) per chapter.

    Both, because a chapter titled for an arm is attributable whole. Reading
    sections alone credited the COMPARATOR with zero words -- Chapter 5 is
    called "The Ferroptosis Engine" and none of its section headings say
    ferroptosis -- which is an error in the direction that flatters every
    other arm, and therefore the direction to look for first.
    """
    sections, chapter_words = [], {}
    chapter, heading, buf = None, None, []

    def flush():
        if heading is not None:
            n = len(" ".join(buf).split())
            sections.append((chapter, heading, n))
            chapter_words[chapter] = chapter_words.get(chapter, 0) + n

    for line in md.split("\n"):
        ch = re.match(r"^## (Chapter \d+|Appendix [A-Z]): (.+)$", line)
        sec = re.match(r"^### (.+)$", line)
        if ch or sec:
            flush()
            buf = []
            if ch:
                chapter, heading = f"{ch.group(1)}: {ch.group(2)}", None
            else:
                heading = sec.group(1)
        elif heading is not None:
            buf.append(line)
    flush()
    return sections, chapter_words


def _book_words(units, needles) -> tuple[int, list[str]]:
    """Whole chapters first, then sections in chapters not already counted."""
    sections, chapter_words = units
    chapters = [c for c in chapter_words
                if c and any(k in c.lower() for k in needles)]
    words = sum(chapter_words[c] for c in chapters)
    hits = list(chapters)
    for c, h, n in sections:
        if c in chapters:
            continue
        if any(k in h.lower() for k in needles):
            words += n
            hits.append(h)
    return words, hits


def _figures_for(entries, needles) -> list[str]:
    out = []
    for e in entries:
        hay = " ".join(str(e.get(k, "")) for k in ("filename", "note")).lower()
        gen = str(e.get("generator", "")).lower()
        if any(k in hay or k in gen for k in needles):
            out.append(e["filename"])
    return sorted(set(out))


def _predictions_for(prereg: str, needles) -> list[str]:
    out = []
    for m in re.finditer(r"^\*\*(P\d+)\. (.+?)\*\*", prereg, flags=re.MULTILINE):
        if any(k in m.group(2).lower() for k in needles):
            out.append(m.group(1))
    return out


def scan() -> dict:
    depth = _json("modality-module-depth.json")
    calib = {a["arm"]: a["verdict"] for a in _json("modality-calibration.json")["arms"]}
    profile = {r["mechanism"]: r for r in _json("census-mechanism-profile.json")["rows"]}
    per_module = {m["module"]: m for m in depth["dedicated"] + depth["shared"]}
    units = _sections(MANUSCRIPT.read_text())
    prereg = PREREG.read_text()
    try:
        import yaml
        figure_entries = yaml.safe_load(FIGURES.read_text())["figures"]
    except Exception:                                     # pragma: no cover
        figure_entries = []

    rows = []
    for spec in ARMS:
        if spec["modules"] == "ENGINE":
            lines, fns = depth["engine_code_lines"], depth["engine_pub_fns"]
            modules = depth["ferroptosis_engine_modules"]
            tests = _engine_rust_tests()
        else:
            names = spec["modules"]
            lines = sum(per_module[m]["code_lines"] for m in names
                        if m in per_module)
            fns = sum(per_module[m]["pub_fns"] for m in names if m in per_module)
            # A module the depth report does not classify (photosensitizer_pk
            # is engine, not a modality file) is still the arm's own code, so
            # it is measured directly rather than dropped.
            for m in names:
                if m not in per_module and (CORE / f"{m}.rs").exists():
                    import importlib.util
                    spec_mc = importlib.util.spec_from_file_location(
                        "modality_coverage", REPO / "scripts" / "modality_coverage.py")
                    mc = importlib.util.module_from_spec(spec_mc)
                    spec_mc.loader.exec_module(mc)
                    body = mc.strip_rust_comments(
                        mc.strip_test_blocks((CORE / f"{m}.rs").read_text()))
                    lines += len([l for l in body.split("\n") if l.strip()])
                    fns += len(re.findall(r"^\s*pub fn ", body, flags=re.MULTILINE))
            modules = len(names)
            tests = sum(_rust_tests_for(m) for m in names)

        words, heads = _book_words(units, spec["heading"])
        rows.append({
            "arm": spec["arm"],
            "label": spec["label"],
            "note": spec["note"],
            "mechanism": spec["mechanism"],
            "census": profile[spec["mechanism"]]["census"]
                      if spec["mechanism"] in profile else None,
            "trials": profile[spec["mechanism"]]["trials"]
                      if spec["mechanism"] in profile else None,
            "modules": modules,
            "code_lines": lines,
            "pub_fns": fns,
            "rust_tests": tests,
            "calibration": calib.get(spec["calibration"]) if spec["calibration"] else None,
            "book_words": words,
            "book_sections": heads,
            "figures": _figures_for(figure_entries, spec["topic"]),
            "predictions": _predictions_for(prereg, spec["topic"]),
        })

    total_words = sum(n for _, _, n in units[0])
    attributed = sum(r["book_words"] for r in rows)
    return {
        "arms": rows,
        "shared_modules": depth["shared_modules"],
        "shared_code_lines": depth["shared_code_lines"],
        "shared_pub_fns": depth["shared_pub_fns"],
        "manuscript_words_in_sections": total_words,
        "manuscript_words_attributed": attributed,
        "absent_arms": [{"arm": a, "why": w} for a, w in ABSENT_ARMS],
    }


def assemble(raw: dict) -> dict:
    base = next(r for r in raw["arms"] if r["arm"] == "RSL3")
    for r in raw["arms"]:
        for axis, key in (("lines", "code_lines"), ("fns", "pub_fns"),
                          ("tests", "rust_tests"), ("words", "book_words")):
            ref = base[key]
            r[f"parity_{axis}"] = round(r[key] / ref, 3) if ref else None
    # The distance to close, on the axis the criticism was about. Reported as a
    # SHORTFALL rather than a ratio because "0.02x" reads as a rounding error
    # and "3,893 lines" reads as the work it is.
    for r in raw["arms"]:
        r["lines_short_of_parity"] = max(0, base["code_lines"] - r["code_lines"])
    raw["comparator"] = base["arm"]
    raw["arms_at_parity"] = sum(1 for r in raw["arms"]
                                if r["arm"] != base["arm"]
                                and (r["parity_lines"] or 0) >= 1.0)
    raw["n_arms"] = len(raw["arms"])
    return raw


def _cell(v, none="--"):
    return none if v is None else (f"{v:,}" if isinstance(v, int) else str(v))


def render(d: dict) -> str:
    base = next(r for r in d["arms"] if r["arm"] == d["comparator"])
    L = ["# How far each arm is from the ferroptosis arm", "",
         "*Generated by `scripts/arm_parity.py --render-only`. Offline; joins "
         "`modality-module-depth.json`, `modality-calibration.json`, "
         "`census-mechanism-profile.json`, `FIGURES.yaml`, `PREREGISTRATION.md` "
         "and the manuscript.*", "",
         "The aggregate gap is already published: the modality arms are roughly "
         "an order of magnitude smaller than the ferroptosis engine by line "
         "count. That single ratio cannot say WHICH arm or HOW FAR, which is "
         "the question a campaign to close it has to answer. This is the same "
         "measurement, per arm, on the axes the work actually has to move.", "",
         "| arm | engine lines | pub fns | rust tests | calibration | book words | figures | predictions |",
         "|---|--:|--:|--:|---|--:|--:|--:|"]
    for r in d["arms"]:
        L.append(
            f"| {r['label']} | {_cell(r['code_lines'])} | {_cell(r['pub_fns'])} "
            f"| {_cell(r['rust_tests'])} | {r['calibration'] or '--'} "
            f"| {_cell(r['book_words'])} | {len(r['figures'])} "
            f"| {len(r['predictions'])} |")
    L += ["",
          f"**{base['code_lines']:,} lines of production code carry the "
          f"ferroptosis arm.** No other arm reaches a tenth of it, and the "
          f"table's own worst row is the one to read first: an arm with no "
          f"module of its own has no engine column to shrink.", ""]

    L += ["## What each arm is short, on the axis the criticism was about", "",
          "| arm | lines | share of ferroptosis | lines short of parity |",
          "|---|--:|--:|--:|"]
    for r in d["arms"]:
        if r["arm"] == d["comparator"]:
            continue
        L.append(f"| {r['label']} | {r['code_lines']:,} | "
                 f"{r['parity_lines']:.3f}x | {r['lines_short_of_parity']:,} |")
    L += ["",
          f"Shared machinery -- {d['shared_modules']} modules, "
          f"{d['shared_code_lines']:,} lines -- is credited to NO arm, so every "
          "row above is a lower bound. It is not divided up because dividing it "
          "would be a judgement about how much of the immune model belongs to "
          "checkpoint blockade rather than to CAR-T, and no such split is "
          "measurable from the code.", ""]

    L += ["## Where the arms are in the literature", "",
          "Volume is context, not a parity axis, and it is not comparable "
          "across rows: it counts how NLM has labelled articles, so an arm "
          "with no taxonomy row reads `no row` rather than zero.", "",
          "| arm | census articles | trials | taxonomy row |",
          "|---|--:|--:|---|"]
    for r in d["arms"]:
        L.append(f"| {r['label']} | {_cell(r['census'], 'no row')} | "
                 f"{_cell(r['trials'], 'no row')} | "
                 f"{r['mechanism'] or 'none'} |")

    L += ["", "## Arms this engine does not have at all", "",
          "A parity table listing only what exists measures the campaign's "
          "progress and hides its scope.", ""]
    for a in d["absent_arms"]:
        L.append(f"- **{a['arm']}** -- {a['why']}")

    L += ["", "## What this does not measure", "",
          "Quality, correctness, or whether any of it is used. A module can be "
          "large and wrong, and an arm could reach every number in the first "
          "table and still be worse science than the arm it is measured "
          "against. What the table supports is the negative: an arm at a "
          "fraction of the comparator on every axis is not represented in any "
          "sense a reader would accept.", "",
          "The book column is attributed by HEADING -- a section counts for an "
          "arm when its own heading names it. That under-counts deliberately: "
          f"{d['manuscript_words_attributed']:,} of "
          f"{d['manuscript_words_in_sections']:,} words in numbered sections "
          "are attributed, and the whole multi-modality chapter sits in the "
          "remainder because its headings name no arm. A keyword rule over body "
          "text was rejected; `analysis/scope-audit.md` records what happened "
          "the last time subject was counted by vocabulary.", ""]
    return "\n".join(L) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--render-only", action="store_true")
    a = ap.parse_args()
    d = assemble(json.loads(OUT_JSON.read_text()) if a.render_only else scan())
    OUT_JSON.write_text(json.dumps(d, indent=1) + "\n")
    OUT_MD.write_text(render(d))
    print(f"wrote {OUT_MD}")
    print(f"wrote {OUT_JSON}")
    base = next(r for r in d["arms"] if r["arm"] == d["comparator"])
    print(f"  comparator {base['label']}: {base['code_lines']:,} lines; "
          f"{d['arms_at_parity']} of {d['n_arms'] - 1} arms at parity")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
