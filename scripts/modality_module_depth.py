"""How much engine the non-ferroptosis arms actually have.

WHY THIS EXISTS
---------------
Chapter 6 carried a sentence since it was written, claiming each modality arm
was one function and a configuration struct beside a ferroptosis engine of many
modules. That was true, it was the honest thing to say, and it stopped being
true without anyone noticing -- which is precisely the defect this repository
keeps finding in its own prose.

The sentence is a MEASUREMENT and was written as an assertion. This measures
it, so the chapter can state where the work actually stands and cannot flatter
itself as modules land or, equally, claim a parity it has not reached.

WHAT IS COUNTED, AND THE CHOICE THAT MATTERS
--------------------------------------------
Public functions, types and constants in PRODUCTION code, per module -- test
blocks and comments stripped, using the same scanner
`analysis/modality-coverage.md` uses, because two documents counting the same
crate two ways is a defect this campaign has already had to fix once.

Modules are split into DEDICATED (a modality's own file) and SHARED (machinery
several arms reach through). THAT SPLIT IS A JUDGEMENT AND IT IS LOAD-BEARING,
which an earlier version of this script hid: the two lists were hardcoded, and
a reviewer moved ONE module between them, moved every headline figure on the
page, and every guard stayed green. Worse, the SHARED bucket is subtracted from
the ENGINE denominator, so every module placed there shrinks the comparator --
an error in the direction that flatters the arms, which is the direction to go
looking first.

Two things follow, and both are in the output rather than in a caveat. The
engine is reported under BOTH allocations (shared counted as engine, and
shared counted as neither), so the ratio ships as a RANGE and no reader
inherits one judgement as a fact. And the bucket membership is pinned by a
test, so moving a module is a decision somebody has to make explicitly rather
than a number that quietly moves.

How many arms actually reach a shared module is MEASURED (`_reach`), not
asserted: an earlier draft of the rendered page hand-wrote a count of arms
reaching through `immune.rs` beside a table row that listed a different one,
and neither figure came from the code.

WHAT THIS DOES NOT MEASURE
--------------------------
Quality, calibration, or whether any of it is used. A module can be large and
wrong. `analysis/modality-calibration.md` carries what is fitted and
`CALIBRATION_STATUS.md` carries what feeds a reported number -- for these arms,
nothing does. Lines of code are the weakest possible evidence of depth and are
reported only because the sentence they replace was about size too.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import re
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
CORE = REPO / "simulations" / "ferroptosis-core" / "src"

OUT_MD = REPO / "analysis" / "modality-module-depth.md"
OUT_JSON = REPO / "analysis" / "modality-module-depth.json"

# A modality's OWN file. Everything else is either shared machinery or the
# ferroptosis engine proper.
DEDICATED = {
    "checkpoint": "Checkpoint blockade: occupancy, the PD-L1 brake, "
                  "mutational-burden antigenicity and the resistance modes",
    "chemo": "Cytotoxic chemotherapy: the cell cycle, log kill and scheduling",
    "radiation": "Radiation + synthetic lethality (PARP)",
    "ablation": "HIFU + irreversible electroporation",
    "oncolytic": "Oncolytic virus spread",
    "adc": "Antibody-drug conjugate bystander effect",
    "adoptive": "CAR-T trafficking, infiltration and activation barriers",
}
# Machinery more than one arm reaches through. Credited to none of them.
SHARED = {
    "immune": "checkpoint blockade, CAR-T, bispecifics, oncolytic ICD, microbiome, mRNA vaccine",
    "immune_spatial": "the spatial immune model",
    "drug_transport": "nanocarrier + ADC delivery profiles",
    "nutrient": "metabolic targeting",
    "cell": "CRISPR knockouts",
}


def _mc():
    spec = importlib.util.spec_from_file_location(
        "modality_coverage", REPO / "scripts" / "modality_coverage.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _measure(mc, stem: str) -> dict:
    path = CORE / f"{stem}.rs"
    if not path.exists():
        return {}
    code = mc.strip_test_blocks(mc.strip_rust_comments(path.read_text()))
    return {
        "module": stem,
        "pub_fns": len(re.findall(r"\bpub fn ", code)),
        "pub_types": len(re.findall(r"\bpub (?:struct|enum) ", code)),
        "pub_consts": len(re.findall(r"\bpub const ", code)),
        "code_lines": len([l for l in code.splitlines() if l.strip()]),
    }


def _reach(mc, stem: str) -> int:
    """How many OTHER crate files name this module's path, in production code.

    The rendered page used to hand-write how many arms reached through
    `immune.rs`, beside a table row listing a different count, and neither
    number was read from anything. This counts files that write `<stem>::` --
    every module in the core crate plus every binary in the workspace, the
    module's own file excluded, comments and test blocks stripped first so a
    doc comment naming the symbol is not a caller.
    """
    pat = re.compile(rf"\b{re.escape(stem)}::")
    # DEDUPLICATED BY PATH. `simulations/*/src/*.rs` fully CONTAINS
    # `ferroptosis-core/src/*.rs`, so the first version counted every core
    # module twice and published "35 other files call `cell.rs`" against a
    # true 25 -- a measured-but-wrong number in a bolded sentence, on the page
    # written to retire a hand-written one. The sibling generators glob
    # `sim-*/src/*.rs` for exactly this reason.
    seen = {p.resolve() for p in CORE.glob("*.rs")}
    seen |= {p.resolve() for p in (REPO / "simulations").glob("sim-*/src/*.rs")}
    n = 0
    for f in sorted(seen):
        if f.stem == stem:
            continue
        code = mc.strip_test_blocks(mc.strip_rust_comments(f.read_text()))
        if pat.search(code):
            n += 1
    return n


def scan() -> dict:
    mc = _mc()
    all_stems = sorted(p.stem for p in CORE.glob("*.rs") if p.stem != "lib")
    dedicated = [dict(_measure(mc, s), serves=DEDICATED[s], reach=_reach(mc, s))
                 for s in sorted(DEDICATED) if (CORE / f"{s}.rs").exists()]
    shared = [dict(_measure(mc, s), serves=SHARED[s], reach=_reach(mc, s))
              for s in sorted(SHARED) if (CORE / f"{s}.rs").exists()]
    ded_names = set(DEDICATED) | set(SHARED)
    engine = [_measure(mc, s) for s in all_stems if s not in ded_names]
    return {"dedicated": dedicated, "shared": shared,
            "engine_module_names": [e["module"] for e in engine],
            "ferroptosis_engine_modules": len(engine),
            "total_modules": len(all_stems),
            "engine_pub_fns": sum(e["pub_fns"] for e in engine),
            "engine_code_lines": sum(e["code_lines"] for e in engine)}


def assemble(raw: dict) -> dict:
    """Totals under BOTH allocations of the shared bucket.

    `engine_*` counts shared machinery as NEITHER side (the narrow engine);
    `engine_*_with_shared` counts it as engine (the wide one). Both are
    defensible and they differ substantially, so the page publishes the ratio
    as a range instead of picking the judgement that makes the arms look
    larger.
    """
    d, sh = raw["dedicated"], raw["shared"]
    ded_lines = sum(x["code_lines"] for x in d)
    ded_fns = sum(x["pub_fns"] for x in d)
    wide_lines = raw["engine_code_lines"] + sum(x["code_lines"] for x in sh)
    wide_fns = raw["engine_pub_fns"] + sum(x["pub_fns"] for x in sh)
    return dict(raw,
                dedicated_modules=len(d),
                dedicated_pub_fns=ded_fns,
                dedicated_code_lines=ded_lines,
                shared_modules=len(sh),
                shared_pub_fns=sum(x["pub_fns"] for x in sh),
                shared_code_lines=sum(x["code_lines"] for x in sh),
                engine_modules_with_shared=raw["ferroptosis_engine_modules"] + len(sh),
                engine_pub_fns_with_shared=wide_fns,
                engine_code_lines_with_shared=wide_lines,
                line_ratio_narrow=round(raw["engine_code_lines"] / ded_lines, 1),
                line_ratio_wide=round(wide_lines / ded_lines, 1),
                fn_ratio_narrow=round(raw["engine_pub_fns"] / ded_fns, 1),
                fn_ratio_wide=round(wide_fns / ded_fns, 1))


def render(d: dict) -> str:
    L = ["# How much engine the non-ferroptosis arms actually have", "",
         "*Generated by `scripts/modality_module_depth.py --render-only`. "
         "Offline; counts production code with test blocks and comments "
         "stripped, using the same scanner `modality-coverage.md` uses.*", "",
         "Chapter 6 carried a sentence since it was written, calling each "
         "modality arm *\"one function and a configuration struct\"* beside "
         "a ferroptosis engine of many modules — which was true, was the "
         "honest thing to say, and stopped being true without anyone "
         "noticing. That is exactly the defect this repository keeps finding "
         "in its own prose, so the sentence is measured here rather than "
         "asserted there.", "",
         "## Modules a modality owns", "",
         "| module | serves | pub fns | types | consts | lines | callers |",
         "|---|---|--:|--:|--:|--:|--:|"]
    for m in d["dedicated"]:
        L.append(f"| `{m['module']}` | {m['serves']} | {m['pub_fns']} | "
                 f"{m['pub_types']} | {m['pub_consts']} | {m['code_lines']} | {m['reach']} |")
    unreached = [m["module"] for m in d["dedicated"] if m["reach"] == 0]
    L += ["",
          "The last column is why it is there. **" +
          (f"`{'`, `'.join(unreached)}` " +
           ("has" if len(unreached) == 1 else "have") +
           " no production caller at all**, so " +
           ("its" if len(unreached) == 1 else "their") +
           " functions are counted below and reached by nothing — a layer "
           "without a caller, which is the defect this campaign keeps finding "
           "in its own work. The count that follows includes them, because "
           "hiding them would be the more flattering error."
           if unreached else
           "Every module here has a production caller**, which was not true "
           "when this column was added: two of them had none, and the table "
           "printed their functions with no way for a reader to tell."), "",
          f"**{d['dedicated_modules']} dedicated modules, "
          f"{d['dedicated_pub_fns']} public functions, "
          f"{d['dedicated_code_lines']:,} lines of production code.** Whatever "
          "else is true, it is not one function and a configuration struct.", "",
          "## Machinery several arms reach through", "",
          "| module | serves | pub fns | lines | files that call it |",
          "|---|---|--:|--:|--:|"]
    for m in d["shared"]:
        L.append(f"| `{m['module']}` | {m['serves']} | {m['pub_fns']} | "
                 f"{m['code_lines']} | {m['reach']} |")
    deepest = max(d["shared"], key=lambda m: m["reach"])
    L += ["",
          "Credited to no single arm, deliberately. "
          f"`{deepest['module']}.rs` is the deepest of them and "
          f"**{deepest['reach']} other files in the workspace call it**, so "
          "assigning its weight to any one arm would overstate that arm and "
          "assigning it to none would understate the engine. Neither number "
          "alone is the answer, so both are reported. (That reach is counted "
          "from the code. An earlier version of this page said \"four arms\" "
          "beside a table row naming six, and neither figure came from "
          "anything.)", "",
          "## Against the ferroptosis engine", "",
          "**And here the judgement above becomes load-bearing, so the answer "
          "is a range rather than a number.** The shared bucket is subtracted "
          "from the engine as well as from the arms, so every module placed "
          "in it shrinks the comparator — an error in the direction that "
          "flatters the arms, which is the direction to check first. Both "
          "allocations are defensible and both are reported:", "",
          "| shared machinery counted as | engine modules | pub fns | lines | "
          "lines per line of modality code |", "|---|--:|--:|--:|--:|",
          f"| neither side | {d['ferroptosis_engine_modules']} | "
          f"{d['engine_pub_fns']} | {d['engine_code_lines']:,} | "
          f"{d['line_ratio_narrow']}x |",
          f"| the engine | {d['engine_modules_with_shared']} | "
          f"{d['engine_pub_fns_with_shared']} | "
          f"{d['engine_code_lines_with_shared']:,} | {d['line_ratio_wide']}x |",
          "",
          "So the modality arms are somewhere between "
          f"**1/{d['line_ratio_narrow']} and 1/{d['line_ratio_wide']} of the "
          f"engine by line count**, and 1/{d['fn_ratio_narrow']} to "
          f"1/{d['fn_ratio_wide']} by public function. On either reading the "
          "ferroptosis engine is still the larger body of work, and it "
          "carries something none of the new modules do: legs fitted against "
          "independent published data, and numbers the manuscript actually "
          "reports. (That clause read \"years of calibration\" until a "
          "reviewer checked the first commit date against it. The repository "
          "is months old, so the claim was false — and it was the one clause "
          "in a paragraph built around a measured interval that had been "
          "exempted from measurement.)", "",
          "## What this does NOT measure", "",
          "**Quality, calibration, or use.** A module can be large and wrong. "
          "`analysis/modality-calibration.md` carries what is fitted — and "
          "records one arm as inadmissible and one as having no fittable "
          "target at all — while `CALIBRATION_STATUS.md` carries what feeds a "
          "reported number. For every arm counted above, that is `N`. **That "
          "is not the same as invisible, and this page said it was:** it "
          "claimed none of this appears in a figure or a claim the manuscript "
          "makes, while `FIGURES.yaml` feeds `analysis/modality-panel.json` "
          "into fig31 and Chapter 6 cites that figure and quotes its numbers. "
          "Three of the modules counted above supply rows in it. What `N` "
          "actually means is narrower and worth stating exactly: none of this "
          "is FITTED to an independent dataset, and nothing here feeds a "
          "number in the manuscript's quantitative chapters. It does appear, "
          "as a described comparison carrying its own uncalibrated status.", "",
          "**Lines of code are the weakest possible evidence of depth** and "
          "are reported only because the sentence they replace was about size "
          "too. A count going up is not progress on its own, and this page "
          "should not be read as saying it is.", ""]
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
    print(f"  dedicated: {d['dedicated_modules']} modules, "
          f"{d['dedicated_pub_fns']} fns, {d['dedicated_code_lines']:,} lines")
    print(f"  ferroptosis engine: {d['ferroptosis_engine_modules']} modules, "
          f"{d['engine_pub_fns']} fns, {d['engine_code_lines']:,} lines")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
