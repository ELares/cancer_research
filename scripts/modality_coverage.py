#!/usr/bin/env python3
"""Which modalities does the ENGINE model, against what the field publishes?

WHY
---
`scripts/scope_audit.py` counts what this project's work is ABOUT and answers
"overwhelmingly ferroptosis": 14 ferroptosis analyses against 1 for any other
therapy, and 8 of 8 preregistered predictions. That is the subject measured,
and it is a fair criticism of a repository whose README invites the whole
cancer community.

What it does not say is what to do instead. This does: it puts the simulation
engine's capability next to the census's attention, per mechanism, so the
ordering of what to build next is a measurement rather than a preference.

WHAT "THE ENGINE MODELS IT" MEANS, and why the strictest reading is used
-----------------------------------------------------------------------
Three tiers, because they are genuinely different and collapsing them flatters
the engine:

  TREATMENT   a `Treatment` variant exists, so the modality can be the thing a
              run APPLIES. Today: Control, RSL3, SDT, PDT -- four, all
              ferroptosis or physical-ROS.
  MODIFIER    a module or parameter changes how an existing treatment lands,
              but cannot be applied on its own. The whole immune stack is here:
              checkpoint panels, Treg/MDSC fields, DC subsets, exhaustion --
              rich machinery reached only through ferroptotic death.
  PROSE-ONLY  the term appears only inside a comment. Reported in its own
              column and never counted as coverage: a docstring naming
              MR-guided focused ultrasound as context is not a model of HIFU,
              and crediting it made HIFU look modelled when nothing runs.

A modality with a checkpoint-blockade parameter but no way to give checkpoint
blockade AS the treatment is MODIFIER, not TREATMENT. That distinction is the
finding, so it is not softened.

WHAT THIS IS NOT
----------------
Not a claim that census volume is importance. A large literature can be large
because it is easy to publish in, and this project's own thesis deliberately
sits on a thin one (the sonodynamic leg is ~30 papers). Volume plus TRIAL SHARE
plus GROWTH is offered so a reader can weigh those separately, and the gap
column is arithmetic on them, not a recommendation.

Not a survey of the field's mechanisms either. The 16 rows are this repo's own
taxonomy, so a modality the taxonomy cannot name is invisible here exactly as
it is everywhere else -- see `analysis/atlas-untagged-partner.md`, where
radiotherapy has no lane and is absorbed into immunotherapy's count.

AN IRONY WORTH RECORDING, and it is a measurement rather than a joke.
`scope_audit.py` classifies this document into its FERROPTOSIS column. That is
correct behaviour: the audit's body route matches how often an analysis MENTIONS
a subject, and a document arguing the repository is ferroptosis-centred names
ferroptosis on nearly every line. The audit says so itself -- "what it measures
is how many analyses MENTION a therapy, which is a different question". Nothing
is wrong here, but the next reader who sees the ferroptosis count go UP when a
non-ferroptosis analysis lands should know why. It is deliberately not stated in
the rendered document: `scope_audit` reads `analysis/*.md`, so a sentence here
about that classification would have to be regenerated after the audit that
reads it, and a two-artifact regeneration cycle is how stale pairs get shipped.

Offline: reads committed JSON and the crate source. No corpus, no census scan.
"""
from __future__ import annotations

import argparse
import ast
import json
import re
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
PROFILE = REPO / "analysis" / "census-mechanism-profile.json"
CORE = REPO / "simulations" / "ferroptosis-core" / "src"
BINS = REPO / "simulations"

OUT_MD = REPO / "analysis" / "modality-coverage.md"
OUT_JSON = REPO / "analysis" / "modality-coverage.json"

# Terms that indicate a module is ABOUT a mechanism, keyed to the census
# taxonomy's names. Deliberately generous on the engine's side: a false
# positive here understates the gap, which is the safe direction for a document
# arguing that the gap is large.
#
# Every term is matched on WORD BOUNDARIES, and a trailing `*` means "stem, any
# suffix". Unbounded substrings do not work here and the failure is silent: the
# first draft matched `t cell`, which is inside `mut cell`, and credited 16 of
# 33 modules -- including `physics`, `stats` and `grid` -- with modelling
# immunotherapy. tests/test_modality_coverage.py pins a must-match and a
# must-not-match string for every term, read out of this table rather than
# restated.
ENGINE_TERMS = {
    # `immun*` not `immune*`: the stem is what precedes the suffix, so
    # `immune*` reaches `immunesuppression` and misses `immunity`,
    # `immunotherapy` and `immunogenic` -- three words and not the fourth.
    "immunotherapy": ("immun*", "checkpoint*", "pd-1", "pd-l1", "ctla-4",
                      "t cell", "t-cell*", "cd8", "dendritic", "icd",
                      "treg*", "mdsc*"),
    "epigenetic": ("hdac*", "histone*", "methylation", "epigenetic*", "dnmt*",
                   "chromatin"),
    "nanoparticle": ("nanoparticle*", "liposom*", "micelle*", "nanocarrier*"),
    "car-t": ("car-t", "car t", "chimeric antigen"),
    "metabolic-targeting": ("glycolysis", "glycolytic", "oxphos", "metabolic*",
                            "glutamin*", "warburg", "lactate"),
    "antibody-drug-conjugate": ("antibody-drug", "adc", "payload*"),
    "synthetic-lethality": ("parp*", "synthetic lethal*", "brca*",
                            "homologous recombination"),
    "oncolytic-virus": ("oncolytic", "virotherapy", "adenovir*"),
    "crispr": ("crispr", "cas9", "guide rna", "sgrna*"),
    "bispecific-antibody": ("bispecific*", "bite", "t-cell engager*"),
    "electrochemical-therapy": ("electroporation", "electrochemical",
                                "irreversible electro*"),
    "sonodynamic": ("sonodynamic", "sdt", "ultrasound", "sonosensitiz*"),
    "hifu": ("hifu", "focused ultrasound", "thermal ablation"),
    "phagocytosis-checkpoint": ("cd47", "sirp*", "phagocytos*"),
    "microbiome": ("microbiome", "microbiota", "bacteri*"),
    "mrna-vaccine": ("mrna vaccine", "neoantigen*", "lipid nanoparticle vaccine"),
}

# The engine's own subject, included so the comparison has its own baseline.
# `ferropto*` rather than `ferroptos*`, which misses `ferroptotic` -- the
# commonest form in this crate. A stem must stop before the point the forms
# diverge, and one measured form is not enough to find that point.
FERROPTOSIS_TERMS = ("ferropto*", "gpx4", "lipid perox*", "fsp1", "acsl4",
                     "slc7a11", "rsl3", "erastin*")


def term_pattern(term: str) -> str:
    """Regex for one term: left-bounded always, right-bounded unless a stem.

    A trailing `*` marks a stem (`immune*` reaches `immunity`, `immunogenic`).
    Without it both ends are bounded, so `t cell` cannot fire inside
    `mut cell` -- the miss that made 16 modules look immunological.
    """
    stem = term.endswith("*")
    body = re.escape(term[:-1] if stem else term)
    return r"\b" + body + ("" if stem else r"\b")


def _matches(text: str, terms) -> bool:
    return any(re.search(term_pattern(t), text) for t in terms)


def _strip_comments(src: str) -> str:
    """Rust source with `//`, `///`, `//!` and `/* */` removed.

    Prose about a mechanism is not a model of it -- see the PROSE-ONLY tier.
    """
    src = re.sub(r"/\*.*?\*/", " ", src, flags=re.S)
    return "\n".join(re.sub(r"//.*$", "", ln) for ln in src.splitlines())


def _module_text() -> tuple[dict, dict]:
    """(code-only text, full text) per module, both lowercased."""
    raw = {p.stem: p.read_text() for p in sorted(CORE.glob("*.rs"))}
    return ({k: _strip_comments(v).lower() for k, v in raw.items()},
            {k: v.lower() for k, v in raw.items()})


def _treatment_variants() -> list:
    """The `Treatment` enum's variants."""
    src = (CORE / "cell.rs").read_text()
    m = re.search(r"pub enum Treatment\s*\{(.*?)\}", src, re.S)
    if not m:
        return []
    return [v.strip().rstrip(",") for v in m.group(1).splitlines()
            if v.strip() and not v.strip().startswith("//")]


def scan() -> dict:
    """Raw counts only: no tier, no total, nothing derived.

    Everything here is a direct reading -- census rows copied from the
    committed profile, module names whose CODE matches a term, module names
    whose comments do and whose code does not. `assemble` turns those into
    tiers and totals, so `--render-only` can rebuild every derived field from
    the artifact instead of trusting it.
    """
    profile = json.loads(PROFILE.read_text())
    mods, full = _module_text()
    variants = _treatment_variants()

    mechanisms = []
    for r in profile["rows"]:
        terms = ENGINE_TERMS.get(r["mechanism"], ())
        code = sorted(n for n, t in mods.items() if _matches(t, terms))
        prose = sorted(n for n, t in full.items()
                       if n not in code and _matches(t, terms))
        mechanisms.append({
            "mechanism": r["mechanism"],
            "census": r["census"],
            "trials": r["trials"],
            "trial_share": r["trial_share"],
            "growth": r.get("growth"),
            "code_modules": code,
            "prose_only_modules": prose,
        })
    return {
        "census": profile["census"],
        "treatment_variants": variants,
        "ferroptosis_modules": sorted(n for n, t in mods.items()
                                      if _matches(t, FERROPTOSIS_TERMS)),
        "module_count": len(mods),
        "binaries": sorted(p.name for p in BINS.glob("sim-*") if p.is_dir()),
        "mechanisms": mechanisms,
    }


def _tier(mech: str, code: list, variants: list) -> str:
    terms = ENGINE_TERMS.get(mech, ())
    for v in variants:
        flat = v.lower()
        if any(re.search(term_pattern(t.replace(" ", "")), flat) for t in terms):
            return "treatment"
    return "modifier" if code else "absent"


def assemble(raw: dict) -> dict:
    """Derive every tier and total from the raw reading."""
    variants = raw["treatment_variants"]
    rows = []
    for m in raw["mechanisms"]:
        rows.append(dict(m, engine_tier=_tier(m["mechanism"],
                                              m["code_modules"], variants)))
    rows.sort(key=lambda x: -x["census"])
    absent = [r for r in rows if r["engine_tier"] == "absent"]
    return dict(raw, rows=rows, absent_count=len(absent),
                absent_census=sum(r["census"] for r in absent),
                absent_trials=sum(r["trials"] for r in absent),
                total_census=sum(r["census"] for r in rows))


def render(d: dict) -> str:
    rows = d["rows"]
    variants = d["treatment_variants"]
    ferro = d["ferroptosis_modules"]
    by_tier: dict = {}
    for r in rows:
        by_tier.setdefault(r["engine_tier"], []).append(r)
    modifier = by_tier.get("modifier", [])
    treatment = by_tier.get("treatment", [])

    L = ["# What the engine models, against what the field publishes", "",
         "*Generated by `scripts/modality_coverage.py --render-only`. Offline; "
         "reads committed JSON and the crate source.*", ""]
    L += [f"The engine expresses **{len(variants)} treatments** — "
          f"`{'`, `'.join(variants)}` — every one of them ferroptosis or "
          f"physical-ROS, and **{len(ferro)} of its {d['module_count']} "
          "modules** name ferroptosis chemistry in code. Against this repo's "
          f"own {len(rows)}-mechanism taxonomy over {d['census']:,} census "
          "articles:", ""]
    L += ["| mechanism | census | trials | trial % | growth | engine | prose only |",
          "|---|--:|--:|--:|--:|---|---|"]
    for r in rows:
        g = f"{r['growth']}×" if r.get("growth") is not None else "—"
        tier = {"treatment": "**treatment**", "modifier": "modifier",
                "absent": "—"}[r["engine_tier"]]
        prose = ", ".join(f"`{m}`" for m in r["prose_only_modules"]) or "—"
        L.append(f"| {r['mechanism']} | {r['census']:,} | {r['trials']:,} | "
                 f"{r['trial_share']:.1f}% | {g} | {tier} | {prose} |")
    L += ["",
          f"**{d['absent_count']} of {len(rows)} mechanisms have no engine "
          f"representation at all**, and they carry {d['absent_census']:,} "
          "census articles "
          f"({d['absent_census'] / d['total_census'] * 100:.0f}% of the "
          f"taxonomy's volume) and {d['absent_trials']:,} registered trials "
          "between them.", ""]

    if modifier:
        names = ", ".join(f"`{r['mechanism']}`" for r in modifier)
        verb = "is a MODIFIER" if len(modifier) == 1 else "are MODIFIERS"
        L += [f"**{names} {verb}, not treatments** — the distinction this "
              "table exists to draw. A modifier changes how an existing "
              "treatment lands and cannot be the thing under test.", ""]
        imm = next((r for r in modifier if r["mechanism"] == "immunotherapy"),
                   None)
        if imm:
            mods_ = ", ".join(f"`{m}`" for m in imm["code_modules"])
            L += ["Immunotherapy is the sharpest case, and the mechanism is "
                  f"explicit in the code: {mods_} reach T-cell killing through "
                  "`immune_kill_probability(activation, rate, brake)` = "
                  "`activation · rate · (1 − brake)`, whose `activation` is "
                  "`dc_activation(local_damp, kd)` = `damp/(damp + kd)`. Every "
                  "writer of the DAMP field in both consumers is the same "
                  "line — `damp_field[idx] += cell.lp_at_grace_end · "
                  "damp_per_lp` — so DAMPs are literally proportional to "
                  "lipid peroxidation at death. There is no "
                  "ferroptosis-independent source of antigen anywhere in the "
                  "engine. With no ferroptotic kill, `activation` is 0 and the "
                  "product is 0 **at every checkpoint-blockade setting**: "
                  "anti-PD-1 enters only through `brake`, a multiplier on a "
                  "term that is already zero. The four-axis checkpoint panel, "
                  "the Treg/MDSC suppressor field, exhaustion and the DC "
                  "subsets are all real and all downstream. Blockade can be "
                  "observed as a coefficient on ferroptosis; it cannot be "
                  "given as a treatment arm. That is a modifier carrying "
                  f"{imm['census']:,} census articles and {imm['trials']:,} "
                  "trials — the largest trial count in the table.", ""]
    if treatment:
        names = ", ".join(f"`{r['mechanism']}`" for r in treatment)
        L += [f"Modelled as a treatment: {names}. `PDT` has no row because "
              "photodynamic therapy is not one of the taxonomy's 16 "
              "mechanisms — a limit of the table, not of the engine.", ""]

    L += ["## What this does not say", "",
          "**Volume is NOT comparable across mechanisms**, and this table is "
          "in volume order only for legibility. `analysis/"
          "census-mechanism-profile.md` states the reason and it applies "
          "unchanged here: descriptor breadth varies enormously, so ordering "
          "by census count is substantially an ordering of how broad each "
          "descriptor is. Read the trial and growth columns beside every "
          "count, and treat none of the three as a ranking of importance — "
          "this project's own thesis deliberately sits on a thin literature, "
          "the sonodynamic leg being ~30 papers. The totals above are "
          "arithmetic, not a recommendation.", "",
          "**A mechanism with no MeSH descriptor is not in the table at all**, "
          "which is a stronger absence than a zero. Tumour-treating fields "
          "have FDA approval in two indications and no descriptor, so they "
          "appear neither as a row nor as a gap; bioelectric modulation is the "
          "same. The engine does not model either, and this document cannot "
          "say so.", "",
          "**The taxonomy bounds the table.** These 16 rows are this repo's "
          "own mechanism list, so a modality it cannot name is invisible here "
          "exactly as it is everywhere else. Radiotherapy is the worked "
          "example: it has no lane, and `analysis/atlas-untagged-partner.md` "
          "measures where it goes instead — of its 88 frozen-corpus articles, "
          "45 are recorded as `immunotherapy`, 11 as `ttfields`, 7 as "
          "`nanoparticle` and 6 as `bispecific-antibody`. So a row above can "
          "be inflated by a modality that has nowhere else to sit. "
          "Chemotherapy and surgery have the same problem and are larger "
          "still.", "",
          "**The engine side errs generous, within code.** A module counts as "
          "covering a mechanism if its *code* names any one of that "
          "mechanism's terms, which over-credits — the safe direction for a "
          "document arguing the gap is large. Comments are stripped first, "
          "because prose about a thing is not a model of it: the only HIFU "
          "term anywhere in the crate is a docstring citing MR-guided focused "
          "ultrasound as context, and counting it made HIFU read as modelled "
          "while nothing runs. Those mentions are kept in the prose-only "
          "column rather than dropped, so both errors stay visible.", ""]
    return "\n".join(L) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--render-only", action="store_true")
    a = ap.parse_args()
    # RE-ASSEMBLE rather than re-render the stored derived fields: the JSON
    # carries `mechanisms`, which is the raw reading, so every tier and total
    # is recomputed here and the stored ones are checkable rather than trusted.
    d = assemble(json.loads(OUT_JSON.read_text()) if a.render_only else scan())
    OUT_JSON.write_text(json.dumps(d, indent=1) + "\n")
    OUT_MD.write_text(render(d))
    print(f"wrote {OUT_MD}")
    print(f"wrote {OUT_JSON}")
    print(f"  {d['absent_count']} of {len(d['rows'])} mechanisms absent from "
          f"the engine; {d['absent_census']:,} census articles")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
