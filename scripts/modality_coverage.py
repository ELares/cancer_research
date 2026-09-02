#!/usr/bin/env python3
"""Which modalities does the ENGINE model, against what the field publishes?

WHY
---
`scripts/scope_audit.py` counts what this project's work is ABOUT and answers
"overwhelmingly ferroptosis": 15 committed ferroptosis analyses against 1 for
any other therapy, and 8 of 13 preregistered predictions. That is the subject
measured, and it is a fair criticism of a repository whose README invites the
whole cancer community.

What it does not say is what the SIMULATION can express. This does. It puts the
engine's capability next to the census's attention per mechanism, so a reader
can see which modalities the engine cannot be asked about at all.

WHAT "THE ENGINE MODELS IT" MEANS, and why the strictest reading is used
-----------------------------------------------------------------------
Three tiers, because they are genuinely different and collapsing them flatters
the engine:

  TREATMENT   a `Treatment` variant exists, so the modality can be the thing a
              run APPLIES.
  MODIFIER    a module's code names the mechanism, so it changes how an
              existing treatment lands, but it cannot be applied on its own.
              The whole immune stack is here.
  PROSE-ONLY  the term appears only inside a comment or a `#[cfg(test)]` block.
              Reported in its own column and never counted as coverage.

THE PROSE-ONLY TIER EXISTS BECAUSE IT WAS MEASURED, TWICE. The only HIFU term
anywhere in the crate is a docstring citing MR-guided focused ultrasound as
context, and counting it made HIFU read as modelled while nothing runs. And an
`assert!` MESSAGE inside a `#[cfg(test)]` block -- "sasp-only config has an
immune effect" -- credited `senescence.rs` with modelling immunotherapy. Both
are prose about a thing rather than a model of it, and `scope_audit.py` had
already been bitten by the second one (`analysis/scope-audit.md`: "four modules
pass the in-code check only via their `#[cfg(test)]` block").

WHAT THIS IS NOT, and the refusal is load-bearing
-------------------------------------------------
**It is NOT a ranking, and no ordering of work is drawn from the volume
column.** `analysis/census-mechanism-profile.md` refuses to draw a
cross-mechanism ranking from volume, because descriptor breadth varies
enormously and a volume ordering is substantially an ordering of how broad each
descriptor is. That refusal is inherited here in full, including the part that
bites: the rows are SORTED by census count for legibility only, and the
actionable content of this document is the ENGINE column, which is binary and
has nothing to do with volume. "This modality cannot be expressed as a
treatment at all" is a fact about our own code. "This modality has the most
articles" is a fact about MeSH.

Not a survey of the field's mechanisms either. The rows are the subset of this
repo's taxonomy that carries a discriminative MeSH descriptor, so a modality
without one is invisible here -- and this document reports how many are in that
position rather than leaving the reader to assume the taxonomy is the table.

AN IRONY WORTH RECORDING, and it is a measurement rather than a joke.
`scope_audit.py` classifies this document into its FERROPTOSIS column. That is
correct behaviour: the audit's body route matches how often an analysis
MENTIONS a subject, and a document arguing the repository is ferroptosis-centred
names ferroptosis on nearly every line. The audit says so itself -- "what it
measures is how many analyses MENTION a therapy, which is a different
question". It is deliberately not stated in the rendered document:
`scope_audit` reads `analysis/*.md`, so a sentence here about that
classification would have to be regenerated after the audit that reads it, and
a two-artifact regeneration cycle is how stale pairs get shipped.

WHAT IT DOES NOT AUTHORISE
--------------------------
`CONTRIBUTING.md`'s layer-freeze policy (calibrate-or-cut, #501) requires a
named calibration target before any new simulation axis lands, and prefers
calibrating an existing layer to adding one. Nothing here bypasses that: this
document says which modalities the engine cannot express, not that any of them
should be built. `analysis/atlas-model-gaps.md`, the nearest sibling
"what is missing" page, carries the same sentence and this one inherits it.

RELATED ISSUES, because three of them already propose parts of this:
#723 (radiotherapy has no lane anywhere -- the worked example below), #724 (the
physical-modality class excludes radiotherapy), #726 (add ionizing radiation
through the exogenous-ROS path, whose oxygen-enhancement-ratio groundwork has
already shipped), #728 (normal-tissue selectivity).

Offline: reads committed JSON/YAML and the crate source. No corpus, no census
scan.
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import yaml

REPO = Path(__file__).resolve().parent.parent
PROFILE = REPO / "analysis" / "census-mechanism-profile.json"
MECH_MAP = REPO / "analysis" / "mesh-mechanism-map.yaml"
CORE = REPO / "simulations" / "ferroptosis-core" / "src"
BINS = REPO / "simulations"

OUT_MD = REPO / "analysis" / "modality-coverage.md"
OUT_JSON = REPO / "analysis" / "modality-coverage.json"

# `lib.rs` is the crate ROOT, not a module: its entire comment-stripped body is
# `pub mod <name>;` lines, so `pub mod immune;` matched `immun*` and credited it
# with modelling immunotherapy. `scope_audit.py` excludes it for the same
# reason and calls it "the crate root rather than a module"; excluding it here
# also makes the two documents' denominators agree, which they did not.
NOT_A_MODULE = {"lib"}

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
#
# WORD BOUNDARIES ARE NOT ENOUGH ON THEIR OWN, and the second draft proved it.
# `glycolytic` and `oxphos` are properly bounded, are real words, and are also
# the names of two `Phenotype` ENUM VARIANTS -- baseline cell states assigned at
# generation. They credited eight modules with modelling METABOLIC TARGETING, a
# therapy, when what the engine models is tumour metabolism as a substrate. The
# therapy-side terms below fire nowhere in the crate, which is the correct
# answer. A must-not-match case pins it.
ENGINE_TERMS = {
    # `immun*` not `immune*`: the stem is what PRECEDES the suffix, so
    # `immune*` reaches neither `immunity` nor `immunotherapy` nor
    # `immunogenic` -- three of the four words that matter.
    "immunotherapy": ("immun*", "checkpoint*", "pd-1", "pd-l1", "ctla-4",
                      "t cell", "t-cell*", "cd8", "dendritic", "icd",
                      "treg*", "mdsc*"),
    "epigenetic": ("hdac*", "histone*", "methylation", "epigenetic*", "dnmt*",
                   "chromatin"),
    "nanoparticle": ("nanoparticle*", "liposom*", "micelle*", "nanocarrier*"),
    # `cart` looks like a false-positive waiting to happen and is not: it is
    # how RUST spells it. The crate's variant is `EffectorSource::CarT`, which
    # lowercases to `cart`, and neither `car-t` nor `car t` can reach a
    # CamelCase identifier. This is the same lesson as the underscore
    # boundary above -- the boundary rule has to know the naming conventions
    # of the language being scanned, and Rust uses CamelCase for types and
    # snake_case for everything else. Bounded on the right, `cart` cannot fire
    # inside `cartilage` or `cartesian`; the case table pins both.
    "car-t": ("car-t", "car t", "cart", "chimeric antigen"),
    "metabolic-targeting": ("glycolysis", "metabolic*", "glutamin*", "warburg",
                            "lactate", "2-deoxyglucose", "dichloroacetate"),
    "antibody-drug-conjugate": ("antibody-drug", "adc", "payload*"),
    "synthetic-lethality": ("parp*", "synthetic lethal*", "brca*",
                            "homologous recombination"),
    "oncolytic-virus": ("oncolytic", "virotherapy", "adenovir*"),
    "crispr": ("crispr", "cas9", "guide rna", "sgrna*"),
    # `bite` was here as the BiTE acronym and matched the ENGLISH VERB: a
    # comment-stripped line in `adoptive.rs` reading "it must actually bite
    # over a run" credited that module with modelling bispecific antibodies.
    # Same class as `t cell` inside `mut cell` and `glycolytic` as a Phenotype
    # variant -- the term was right and the referent was wrong, which no
    # boundary rule can fix because the string genuinely appears. The acronym
    # is unreachable case-insensitively and is dropped; the two remaining
    # terms are unambiguous, and `blinatumomab` is added because a named drug
    # cannot collide with English.
    "bispecific-antibody": ("bispecific*", "t-cell engager*", "blinatumomab"),
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

# The untreated arm. Counted separately from the treatments, because calling it
# one made the headline read "4 treatments, every one of them ferroptosis or
# physical-ROS" when the fourth applies nothing at all
# (`physics.rs`: `Treatment::Control => 0.0`).
CONTROL_VARIANT = "Control"

# What each arm's DOMINANT lethal channel is. Not decoration: the headline
# sentence used to assert "every one of them ferroptosis or physical-ROS", and
# the first arm that was neither -- `Radiation`, whose dominant lesion is the
# DNA double-strand break and which does not pass through `CellState` at all --
# made it false the moment it landed. A guard requires every variant to be
# classified here, so the next arm cannot slip in under a sentence that no
# longer describes it.
TREATMENT_KIND = {
    "RSL3": "ferroptosis",
    "SDT": "physical-ROS",
    "PDT": "physical-ROS",
    "Radiation": "DNA damage, with a separate ferroptosis channel",
    "Immunotherapy": "immune cascade, no ferroptotic death required",
    "AdoptiveCell": "redirected effectors, bypassing DC priming",
    "OncolyticVirus": "lysis into the shared ICD chain",
    "Ablation": "threshold destruction, not a dose-response",
    "Chemotherapy": "log kill weighted by the cell cycle, not a redox state",
    "AntibodyDrugConjugate": "ferroptosis payload, delivery-limited",
}


# `\b` is the wrong boundary for Rust, and the failure is silent and one-sided.
# `_` is a word character, so `\bsdt\b` does NOT match `sdt_ros`, `\bacsl4\b`
# does not match `acsl4_multiplier`, and `\bcd47\b` does not match
# `cd47_blockade` -- the crate's whole naming convention is invisible. Measured:
# `acsl4.rs`, the ACSL4 module, was not counted as naming ferroptosis chemistry
# in code. The error runs AGAINST this document's own argument (it understates
# what the engine expresses), which is the direction that keeps it unnoticed.
#
# So the boundary is "not a letter or a digit", which makes `_` a separator
# while still refusing `mut cell` for `t cell` and `cd80` for `cd8`.
_LEFT = r"(?<![A-Za-z0-9])"
_RIGHT = r"(?![A-Za-z0-9])"


def term_pattern(term: str) -> str:
    """Regex for one term: left-bounded always, right-bounded unless a stem.

    A trailing `*` marks a stem (`immun*` reaches `immunity`, `immunogenic`).
    Without it both ends are bounded, so `t cell` cannot fire inside
    `mut cell` -- the miss that made 16 modules look immunological -- while
    `sdt_ros` and `acsl4_multiplier` still count, which `\b` refused.
    """
    stem = term.endswith("*")
    body = re.escape(term[:-1] if stem else term)
    return _LEFT + body + ("" if stem else _RIGHT)


def _matches(text: str, terms) -> bool:
    return any(re.search(term_pattern(t), text) for t in terms)


def strip_rust_comments(src: str) -> str:
    """Rust source with comments removed, string literals preserved.

    A regex cannot do this correctly and the failures are the kind that credit
    prose as code:

    * Rust block comments NEST. `/\\*.*?\\*/` closes at the first inner `*/`
      and leaves the outer comment's tail behind as "code".
    * `//` inside a string literal (`"http://x"`) truncates real code.
    * `/*` inside a string literal swallows everything to the next `*/`.

    None of the three is live in this crate today -- there is no `/*` anywhere
    in it -- so this is a scanner written against latent defects rather than
    measured ones, and its tests plant each construct deliberately.
    """
    out = []
    i, n = 0, len(src)
    while i < n:
        c = src[i]
        # Raw string: r"..." / r#"..."# / r##"..."##
        if c == "r" and i + 1 < n and src[i + 1] in '#"':
            j = i + 1
            hashes = 0
            while j < n and src[j] == "#":
                hashes += 1
                j += 1
            if j < n and src[j] == '"':
                close = '"' + "#" * hashes
                end = src.find(close, j + 1)
                end = n if end == -1 else end + len(close)
                out.append(src[i:end])
                i = end
                continue
        if c == '"':
            j = i + 1
            while j < n:
                if src[j] == "\\":
                    j += 2
                    continue
                if src[j] == '"':
                    j += 1
                    break
                j += 1
            out.append(src[i:j])
            i = j
            continue
        if c == "'":
            # A char literal is `'x'` or `'\n'`; anything else starting with a
            # quote is a lifetime (`&'a str`) and must not open a literal.
            m = re.match(r"'(\\.|[^\\'])'", src[i:])
            if m:
                out.append(m.group(0))
                i += m.end()
                continue
            out.append(c)
            i += 1
            continue
        if c == "/" and i + 1 < n and src[i + 1] == "/":
            j = src.find("\n", i)
            i = n if j == -1 else j
            continue
        if c == "/" and i + 1 < n and src[i + 1] == "*":
            depth, j = 1, i + 2
            while j < n and depth:
                if src.startswith("/*", j):
                    depth += 1
                    j += 2
                elif src.startswith("*/", j):
                    depth -= 1
                    j += 2
                else:
                    j += 1
            # LINE-PRESERVING: keep one newline per line removed. Every
            # stripper here reports line numbers into the real file, and the
            # first version counted them over the STRIPPED text -- two of the
            # seven DAMP writers in the committed JSON pointed at the wrong
            # line, including the 3D binary's principal one, because a
            # mid-file `#[cfg(test)]` item shifted everything after it by six.
            out.append("\n" * src.count("\n", i, j))
            i = j
            continue
        out.append(c)
        i += 1
    return "".join(out)


# `#[cfg(test)]`, but also `#[cfg(all(test, feature = "x"))]` and `#[cfg( test )]`.
# The first version matched only the bare spelling, so a perfectly ordinary
# `#[cfg(all(test, ...))]` left an entire test module in the scanned code --
# the PROSE-ONLY false positive the tier exists to prevent.
_CFG_TEST = re.compile(r"#\[\s*cfg\s*\((?:[^\[\]])*?\)\s*\]")


def _cfg_test_matches(text: str) -> bool:
    """Does this `#[cfg(...)]` attribute gate ON `test` (not on NOT-test)?

    `re.search(r"\\btest\\b", ...)` also matches `#[cfg(not(test))]`, whose
    item is PRODUCTION code -- it compiles when tests are off. Stripping it is
    a false negative in the direction that argues for building code that
    already exists, which is the direction this whole file is most careful
    about. Latent here (the crate has none), planted in a test.
    """
    if not re.search(r"\btest\b", text):
        return False
    # Strip every `not(...)` group, innermost first, then re-check: if the only
    # occurrences of `test` were inside a negation, this does not gate on test.
    stripped = text
    while True:
        reduced = re.sub(r"not\s*\([^()]*\)", "", stripped)
        if reduced == stripped:
            break
        stripped = reduced
    return re.search(r"\btest\b", stripped) is not None


def strip_test_blocks(src: str) -> str:
    """Remove every `test`-gated item, preserving line numbers.

    Run AFTER comment stripping so a brace in a comment cannot unbalance the
    match. Four things the first version got wrong, none of them live in this
    crate today and all of them legal Rust -- the same class its sibling
    `strip_rust_comments` was hardened against and this one was not:

    * `#[cfg(all(test, feature = "x"))]` and `#[cfg( test )]` did not match at
      all, so the whole block survived as "engine code".
    * The item after the attribute is not always a `mod`: `use`, `const`, `fn`
      and `impl` are all legal. Brace-matching from the next `{` anywhere in
      the file swallowed the following PRODUCTION item, or the rest of the
      file.
    * A `'{'` or `'"'` char literal inside the block unbalanced the depth
      counter and dropped everything after it.
    * Raw strings (`r#"{...}"#`) were not recognised, and three already exist
      inside test blocks in this crate. They are harmless today only because
      those blocks run to EOF, where an overshoot and a correct match are
      indistinguishable.

    Every one of those is planted in a test rather than waited for.
    """
    out, i, n = [], 0, len(src)
    while i < n:
        m = _CFG_TEST.search(src, i)
        while m and not _cfg_test_matches(m.group(0)):
            m = _CFG_TEST.search(src, m.end())
        if not m:
            out.append(src[i:])
            break
        out.append(src[i:m.start()])
        j, end = m.end(), None
        # A brace-bodied item (`mod`, `fn`, `impl`, `struct`) ends at its
        # matching `}`; a statement item (`use`, `const`, `type`) ends at the
        # first `;`. Whichever terminator comes first decides which it is.
        brace = _scan_for(src, j, "{")
        semi = _scan_for(src, j, ";")
        if brace is not None and (semi is None or brace < semi):
            end = _matching_brace(src, brace)
        elif semi is not None:
            end = semi + 1
        if end is None:
            end = n
        out.append("\n" * src.count("\n", m.start(), end))
        i = end
    return "".join(out)


def _scan_for(src: str, i: int, ch: str):
    """Index of the next `ch` at depth zero, skipping strings and chars."""
    n = len(src)
    while i < n:
        i = _skip_literal(src, i)
        if i >= n:
            return None
        if src[i] == ch:
            return i
        i += 1
    return None


def _matching_brace(src: str, open_idx: int) -> int:
    """Index just past the `}` matching the `{` at `open_idx`."""
    depth, i, n = 1, open_idx + 1, len(src)
    while i < n and depth:
        i = _skip_literal(src, i)
        if i >= n:
            break
        if src[i] == "{":
            depth += 1
        elif src[i] == "}":
            depth -= 1
        i += 1
    return i


def _skip_literal(src: str, i: int) -> int:
    """If a literal starts at `i`, return the index just past it.

    Handles raw strings, byte strings, ordinary strings and char literals, so
    a brace or a quote inside one cannot move the brace depth. A `'` that does
    not open a char literal is a lifetime and is left alone.
    """
    n = len(src)
    if i >= n:
        return i
    m = re.match(r"(?:b?r(#*)\"|b?\")", src[i:])
    if m:
        if m.group(1) is not None and "r" in m.group(0):
            close = '"' + "#" * len(m.group(1))
            end = src.find(close, i + m.end())
            return n if end == -1 else end + len(close)
        j = i + m.end()
        while j < n:
            if src[j] == "\\":
                j += 2
                continue
            if src[j] == '"':
                return j + 1
            j += 1
        return n
    m = re.match(r"'(?:\\.|[^\\'])'", src[i:])
    if m:
        return i + m.end()
    return i


def _module_text() -> tuple[dict, dict]:
    """(engine-code text, full text) per module, both lowercased.

    Engine code = comments stripped, `#[cfg(test)]` items stripped, `lib.rs`
    excluded. Full text is everything, and the difference between them is the
    PROSE-ONLY column.
    """
    raw = {p.stem: p.read_text() for p in sorted(CORE.glob("*.rs"))
           if p.stem not in NOT_A_MODULE}
    code = {k: strip_test_blocks(strip_rust_comments(v)).lower()
            for k, v in raw.items()}
    return code, {k: v.lower() for k, v in raw.items()}


def _treatment_variants() -> list:
    """The `Treatment` enum's variants.

    Attribute lines are skipped: `#[default]` above a variant is legal Rust and
    the first version counted it as one, so a derive change would have moved
    the headline count.
    """
    src = (CORE / "cell.rs").read_text()
    m = re.search(r"pub enum Treatment\s*\{(.*?)\}", src, re.S)
    if not m:
        return []
    out = []
    for line in m.group(1).splitlines():
        v = line.strip().rstrip(",")
        if not v or v.startswith("//") or v.startswith("#["):
            continue
        out.append(v)
    return out


def _baseline_antigenicity() -> dict:
    """Is there a ferroptosis-INDEPENDENT antigen source, and is it on?

    Derived, because this file's central paragraph used to say both immune
    models are "gated on ferroptotic death by construction" and the
    immunotherapy arm (#728) made that false the moment it landed. A sentence
    describing the engine has to be read off the engine.
    """
    params = strip_test_blocks(strip_rust_comments((CORE / "params.rs").read_text()))
    immune = strip_test_blocks(strip_rust_comments((CORE / "immune.rs").read_text()))
    spatial = strip_test_blocks(
        strip_rust_comments((CORE / "immune_spatial.rs").read_text()))
    field = "baseline_antigenicity"
    m = re.search(rf"{field}:\s*([0-9.]+)", params)
    return {
        "field": field,
        "exists": field in params,
        "default": float(m.group(1)) if m else None,
        "consumed_by": sorted(
            n for n, t in (("immune", immune), ("immune_spatial", spatial))
            if field in t),
    }


def _immune_models() -> list:
    """The engine's immune kill paths and every file that calls them.

    Derived, because the first version NAMED five modules as reaching T-cell
    killing and four of them did not -- one on the strength of two `assert!`
    messages. There are TWO immune models, not one, and the second has four
    consumers the first version did not mention at all.
    """
    files = sorted(CORE.glob("*.rs")) + sorted(BINS.glob("sim-*/src/*.rs"))
    models = [
        {"kill_fn": "immune_kill_probability", "module": "immune_spatial",
         "kind": "spatial DAMP field"},
        {"kill_fn": "immune_cascade", "module": "immune",
         "kind": "well-mixed cascade"},
    ]
    for mo in models:
        callers = []
        for p in files:
            body = strip_test_blocks(strip_rust_comments(p.read_text()))
            calls = re.findall(rf"\b{mo['kill_fn']}\s*\(", body)
            if not calls:
                continue
            name = p.stem if p.parent == CORE else p.parent.parent.name
            if name != mo["module"]:
                callers.append(name)
        mo["callers"] = sorted(set(callers))
    return models


def _rust_files() -> list:
    return sorted(CORE.glob("*.rs")) + sorted(BINS.glob("sim-*/src/*.rs"))


def scan() -> dict:
    """Raw counts only: no tier, no total, nothing derived."""
    profile = json.loads(PROFILE.read_text())
    mech_map = yaml.safe_load(MECH_MAP.read_text())
    code, full = _module_text()
    variants = _treatment_variants()

    named = sorted(set(mech_map.get("mechanisms", {}))
                   | set(mech_map.get("unmeasurable", {})))
    rows_present = {r["mechanism"] for r in profile["rows"]}

    mechanisms = []
    for r in profile["rows"]:
        terms = ENGINE_TERMS.get(r["mechanism"], ())
        hits = sorted(n for n, t in code.items() if _matches(t, terms))
        prose = sorted(n for n, t in full.items()
                       if n not in hits and _matches(t, terms))
        mechanisms.append({
            "mechanism": r["mechanism"],
            "census": r["census"],
            "trials": r["trials"],
            "trial_share": r["trial_share"],
            "growth": r.get("growth"),
            "code_modules": hits,
            "prose_only_modules": prose,
        })
    return {
        "census": profile["census"],
        "treatment_variants": variants,
        "treatment_kinds": {v: TREATMENT_KIND.get(v) for v in variants
                            if v != CONTROL_VARIANT},
        "ferroptosis_modules": sorted(n for n, t in code.items()
                                      if _matches(t, FERROPTOSIS_TERMS)),
        "module_count": len(code),
        "taxonomy_named": named,
        "taxonomy_without_a_row": sorted(set(named) - rows_present),
        "immune_models": _immune_models(),
        "baseline_antigenicity": _baseline_antigenicity(),
        # The two detector defects the report narrates, RECOMPUTED under the
        # rules the report currently describes rather than quoted from the
        # version that had them. Both were typed, and both had gone stale:
        # "eight modules" was measured before `lib.rs` was excluded and test
        # blocks stripped, and "16 of 33" used a denominator that counted the
        # crate root.
        "phenotype_term_credits": sorted(
            n for n, t in code.items() if _matches(t, ("glycolytic", "oxphos"))),
        "unbounded_t_cell_credits": sorted(
            n for n, t in code.items() if "t cell" in t),
        "mechanisms": mechanisms,
    }


# WHICH `Treatment` VARIANT, IF ANY, LETS A RUN SELECT THIS MECHANISM.
#
# This was a REGEX AGAINST THE VARIANT'S SPELLING, and it was a spelling
# accident rather than a capability measurement. `oncolytic` is right-bounded
# so it cannot match inside `OncolyticVirus`; `adc` cannot match inside
# `AntibodyDrugConjugate`; `cart` cannot match `AdoptiveCell`; `hifu` cannot
# match `Ablation`. Four mechanisms whose arms exist and are RUN by
# `sim-modality-panel` were reported as modifiers, so the page told a reader
# they could not be selected while the panel selected them.
#
# No regex can express this, because the enum is named after what the arm DOES
# and the taxonomy is named after what the literature CALLS it. So it is an
# explicit mapping -- a judgement, recorded where it can be disputed, and
# pinned by `tests/test_modality_coverage.py` so a new variant cannot land
# without someone deciding which mechanism (if any) it makes selectable.
VARIANT_FOR = {
    "sonodynamic": ["SDT"],
    "immunotherapy": ["Immunotherapy"],
    "car-t": ["AdoptiveCell"],
    "bispecific-antibody": ["AdoptiveCell"],
    "oncolytic-virus": ["OncolyticVirus"],
    "antibody-drug-conjugate": ["AntibodyDrugConjugate"],
    "hifu": ["Ablation"],
    "electrochemical-therapy": ["Ablation"],
}
# Variants that make NO taxonomy mechanism selectable, each with the reason.
# `RSL3`, `PDT` and `Radiation` are real arms whose mechanisms this taxonomy
# does not name -- GPX4 inhibitors, photodynamic therapy and radiotherapy have
# no row in the 16, which is a limit of the taxonomy and not of the engine.
UNMAPPED_VARIANTS = {
    "RSL3": "ferroptosis induction has no taxonomy row",
    "PDT": "photodynamic therapy has no taxonomy row",
    "Radiation": "radiotherapy has no taxonomy row",
    "Chemotherapy": "cytotoxic chemotherapy has no taxonomy row -- the arm "
                    "that reaches the most patients is one this project's "
                    "mechanism vocabulary cannot name",
}


def _tier(mech: str, code: list, variants: list) -> str:
    live = {v for v in variants if v != CONTROL_VARIANT}
    if set(VARIANT_FOR.get(mech, ())) & live:
        return "treatment"
    return "modifier" if code else "absent"


def assemble(raw: dict) -> dict:
    """Derive every tier and total from the raw reading."""
    variants = raw["treatment_variants"]
    rows = [dict(m, engine_tier=_tier(m["mechanism"], m["code_modules"], variants))
            for m in raw["mechanisms"]]
    rows.sort(key=lambda x: -x["census"])
    absent = [r for r in rows if r["engine_tier"] == "absent"]
    return dict(raw, rows=rows,
                active_treatments=[v for v in variants if v != CONTROL_VARIANT],
                absent_count=len(absent),
                absent_census=sum(r["census"] for r in absent),
                absent_trials=sum(r["trials"] for r in absent),
                total_census=sum(r["census"] for r in rows))


def _antigen_paragraph(d: dict) -> list:
    """The gating sentence, derived from the crate rather than asserted.

    It used to read "both are gated on ferroptotic death by construction",
    which the immunotherapy arm (#728) falsified the day it landed -- one
    model now has a ferroptosis-independent path. Which model, and whether it
    is on, are read off `params.rs` and the two modules.
    """
    b = d.get("baseline_antigenicity") or {}
    if not b.get("exists"):
        return ["Both paths are gated on ferroptotic death by construction: "
                "there is no ferroptosis-independent antigen source anywhere "
                "in the engine.", ""]
    on = b.get("default") not in (0.0, None)
    consumers = ", ".join(f"`{c}`" for c in b.get("consumed_by", [])) or "nothing"
    ungated = [m["module"] for m in d["immune_models"]
               if m["module"] in b.get("consumed_by", [])]
    gated = [m["module"] for m in d["immune_models"] if m["module"] not in ungated]
    return [
        f"**One of them is no longer gated on ferroptotic death.** "
        f"`ImmuneParams::{b['field']}` (default `{b['default']}`) adds a "
        "ferroptosis-independent presenting fraction, which is the mechanism "
        "real tumours use and the reason checkpoint blockade is a "
        f"monotherapy at all. It is consumed by {consumers}"
        + (f", so `{'`, `'.join(gated)}` remains gated and "
           f"`{'`, `'.join(ungated)}` does not" if gated and ungated else "")
        + ". "
        + ("It is ON by default, so every immune number in this repository "
           "now includes it."
           if on else
           "It is **OFF by default**, so every committed number is unmoved and "
           "the engine can be ASKED the question rather than answering it: "
           "nothing is fit to the published anti-PD-1 monotherapy response "
           "band yet, and `CALIBRATION_STATUS.md` says so."),
        "",
        ("That is why this row read MODIFIER for so long. The tier is about "
         "whether a run can APPLY the modality, which needs a `Treatment` "
         "variant; the field made the question reachable, and the variant "
         "made it askable. Before either, anti-PD-1 multiplied a term that "
         "was structurally zero."
         if any(r["mechanism"] == "immunotherapy" and r["engine_tier"] == "treatment"
                for r in d["rows"]) else
         "That is why this row still reads MODIFIER rather than **treatment**. "
         "The tier is about whether a run can APPLY the modality, and that "
         "needs a `Treatment` variant and a binary that uses it. What changed "
         "is that the question is now reachable at all — before it, anti-PD-1 "
         "multiplied a term that was structurally zero."),
        "",
    ]


def render(d: dict) -> str:
    rows = d["rows"]
    active = d["active_treatments"]
    ferro = d["ferroptosis_modules"]
    by_tier: dict = {}
    for r in rows:
        by_tier.setdefault(r["engine_tier"], []).append(r)
    modifier = by_tier.get("modifier", [])
    treatment = by_tier.get("treatment", [])
    n_named = len(d["taxonomy_named"])
    n_norow = len(d["taxonomy_without_a_row"])

    L = ["# What the engine models, against what the field publishes", "",
         "*Generated by `scripts/modality_coverage.py --render-only`. Offline; "
         "reads committed JSON/YAML and the crate source.*", ""]
    kinds: dict = {}
    for v in active:
        kinds.setdefault(TREATMENT_KIND.get(v, "unclassified"), []).append(v)
    kind_text = "; ".join(
        f"{', '.join('`' + v + '`' for v in vs)} — {k}"
        for k, vs in kinds.items())
    L += [f"The engine expresses **{len(active)} treatments** ({kind_text}), "
          f"plus an untreated `{CONTROL_VARIANT}` arm that applies nothing. "
          f"**{len(ferro)} of its {d['module_count']} modules** name "
          "ferroptosis chemistry in code, with comments and `#[cfg(test)]` "
          "blocks stripped first. Against the "
          f"{len(rows)} mechanisms of this repo's taxonomy that carry a "
          f"discriminative MeSH descriptor, over {d['census']:,} census "
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
    n_treat = len(by_tier.get("treatment", []))
    n_mod = len(by_tier.get("modifier", []))
    if d["absent_count"]:
        L += ["",
              f"**{d['absent_count']} of {len(rows)} mechanisms have no engine "
              f"representation at all**, and they carry {d['absent_census']:,} "
              "census articles "
              f"({d['absent_census'] / d['total_census'] * 100:.0f}% of the "
              f"table's volume) and {d['absent_trials']:,} registered trials "
              "between them. That is the only actionable content here, and it "
              "is binary: it says which questions the engine cannot be asked, "
              "not which are worth asking.", ""]
    else:
        L += ["",
              f"**Every one of the {len(rows)} mechanisms now has some engine "
              "representation.** That was not true when this document was "
              "written — it opened at thirteen absent, carrying 90,019 census "
              "articles — and the column that mattered then is now empty.", "",
              "**Which makes the remaining distinction the whole content, and "
              "it is a harder one.** Presence is not applicability. "
              f"{n_treat} of {len(rows)} can be APPLIED as a treatment — a "
              "`Treatment` variant a run can select — and the other "
              f"{n_mod} are MODIFIERS: real code, reachable, and only ever a "
              "coefficient on something else. An engine where every mechanism "
              "is present and one is applicable answers a narrower set of "
              "questions than the first count suggests, and this table's job "
              "is now to say so rather than to count absences.", "",
              "**And presence says nothing about calibration.** "
              "`simulations/calibration/CALIBRATION_STATUS.md` carries a row "
              "per layer with its named target and its "
              "used-in-any-reported-number status; most of these arms are "
              "`N`. A mechanism the engine can name, cannot apply, and has not "
              "fitted is a long way from one it can answer a question with.", ""]

    if modifier:
        names = ", ".join(f"`{r['mechanism']}`" for r in modifier)
        verb = "is a MODIFIER" if len(modifier) == 1 else "are MODIFIERS"
        L += [f"**{names} {verb}, not treatments** — the distinction this "
              "table exists to draw. A modifier changes how an existing "
              "treatment lands and cannot be the thing under test.", ""]
        imm = next((r for r in rows if r["mechanism"] == "immunotherapy"), None)
        if imm:
            models = d["immune_models"]
            desc = "; ".join(
                f"`{m['module']}::{m['kill_fn']}` ({m['kind']}), called from "
                + (", ".join(f"`{c}`" for c in m["callers"]) or "nowhere")
                for m in models)
            tier_now = imm["engine_tier"]
            opener = (
                "Immunotherapy was the sharpest case in this table, and it is "
                "the one that moved. It is now **applicable as a treatment** — "
                "`Treatment::Immunotherapy` — where for most of this "
                "document's life it was a MODIFIER, and the reason it was one "
                "is worth keeping rather than deleting, because it is the "
                "clearest example of what that tier means."
                if tier_now == "treatment" else
                "Immunotherapy is the sharpest case.")
            L += [f"{opener} The engine has "
                  f"**{len(models)} immune kill paths** — {desc}.", "",
                  "The two gate on ferroptotic death in DIFFERENT ways, and "
                  "conflating them was wrong in an earlier draft of this "
                  "paragraph. The spatial one gates on DAMPs: "
                  "`activation · rate · (1 − brake)` with "
                  "`activation = dc_activation(local_damp, kd)` = "
                  "`damp/(damp + kd)`, and DAMPs are written as "
                  "`lp_at_grace_end · damp_per_lp`. The well-mixed one gates "
                  "on the COUNT: `mature_dcs = activation · maturation_rate · "
                  "n_dead`, so zero ferroptotic deaths give zero kills even if "
                  "DAMPs are somehow non-zero — verified by mutation, where "
                  "injecting a constant `damp_per_dead_cell` produces no kills "
                  "and injecting a constant `mature_dcs` does. Either way "
                  "anti-PD-1 enters only through `brake`, a multiplier on a "
                  "term that is already zero. The four-axis checkpoint panel, "
                  "the Treg/MDSC suppressor field, exhaustion and the DC "
                  "subsets are all real and all downstream of it. **At the "
                  "committed defaults that is still the whole story** — see "
                  "the next paragraph for the one path that now has an "
                  "alternative, and note it is off. "
                  + ("Turning that path on is what made the arm applicable: "
                     "an unconfigured run still behaves exactly as it did."
                     if tier_now == "treatment" else
                     "That is a modifier.")
                  + f" It carries {imm['census']:,} census articles and "
                  f"{imm['trials']:,} trials, the largest trial count in the "
                  "table.", "",
                  ] + _antigen_paragraph(d) + [
                  "**What is proved, and what is not.** The well-mixed chain "
                  "is proved END TO END in Rust: "
                  "`no_ferroptotic_death_end_to_end_means_no_kills_at_any_"
                  "blockade` generates cells, runs the ferroptosis engine with "
                  "death made impossible, collects lipid peroxidation at "
                  "death and feeds the immune cascade, asserting zero kills "
                  "across three treatments, three anti-PD-1 efficacies and "
                  "both blockade arms — and pairs that with a POSITIVE case, "
                  "because a test that only asserts zero is satisfied by an "
                  "engine that does nothing. Both kill formulas are pinned "
                  "at the function level. **What no test here covers is a "
                  "ferroptosis-independent antigen source added inside "
                  "`sim-tme` or `sim-tme-3d`'s own loops**, which compose the "
                  "spatial model outside the crate.", "",
                  "That limit is stated rather than papered over, because "
                  "three successive attempts to cover it with a source scan "
                  "were each defeated by an ordinary Rust idiom that put the "
                  "mutation and the field name on different lines — a helper "
                  "function, an `iter_mut().for_each`, a plain `for` loop. A "
                  "guard that records a property its scan cannot DECIDE is "
                  "worse than no guard, so that scan was removed rather than "
                  "widened a fourth time.", ""]

    if treatment:
        names = ", ".join(f"`{r['mechanism']}`" for r in treatment)
        L += [f"Modelled as a treatment: {names}. `PDT` and `RSL3` have no row "
              "because photodynamic therapy and GPX4 inhibitors are not among "
              "the mechanisms this taxonomy names — a limit of the table, not "
              "of the engine.", ""]

    L += ["## What this does not say", "",
          "**This is not a ranking, and no ordering of work is drawn from the "
          "volume column.** `analysis/census-mechanism-profile.md` refuses a "
          "cross-mechanism ranking because descriptor breadth varies "
          "enormously, so a volume ordering is substantially an ordering of "
          "how broad each descriptor is. That refusal is inherited here in "
          "full. The rows are sorted by census count for legibility only; the "
          "content is the ENGINE column, which is binary and independent of "
          "volume. `Cannot be expressed as a treatment` is a fact about this "
          "repository's code. `Has the most articles` is a fact about MeSH, "
          "and this project's own thesis deliberately sits on a thin "
          "literature.", "",
          f"**The taxonomy names {n_named} mechanisms and only {len(rows)} "
          "have a row.** The rest carry no discriminative MeSH descriptor and "
          "so cannot be counted at census scale at all: "
          + ", ".join(f"`{m}`" for m in d["taxonomy_without_a_row"])
          + f" ({n_norow}). Their absence from this table is a stronger "
          "statement than a zero, and it is not a statement about the engine — "
          "tumour-treating fields have FDA approval in two indications, are "
          "counted elsewhere in this repository from the frozen corpus, and "
          "are simply unmeasurable through MeSH. The criterion is the "
          "descriptor, not the mechanism's importance and not always its "
          "existence: `analysis/mesh-mechanism-map.yaml` records that "
          "cold-atmospheric-plasma HAS a descriptor and is excluded for "
          "breadth.", "",
          "**And a modality the taxonomy never names is invisible even from "
          "that count.** Radiotherapy is the worked example. Of the 88 "
          "articles in the frozen corpus whose TITLES are about it, 75 carry "
          "some other mechanism's tag — `immunotherapy` 45, `ttfields` 11, "
          "`nanoparticle` 7, `bispecific-antibody` 6, and the tags overlap, so "
          "24 of the 88 carry two or more "
          "(`analysis/atlas-untagged-partner.md`). Half the corpus mentions "
          "radiotherapy somewhere (2,436 articles, 50.4%). So a row above can "
          "be inflated by a modality that has nowhere else to sit. "
          "Chemotherapy and surgery have the same problem and are larger "
          "still.", "",
          "**The engine side errs generous, within production code.** A module "
          "counts as covering a mechanism if its code names any one of that "
          "mechanism's terms, which over-credits — the safe direction for a "
          "document arguing the gap is large. Three things are stripped first, "
          "each because prose about a thing is not a model of it and each "
          "was measured crediting prose as code: comments "
          "(the only HIFU term in the crate is a docstring citing MR-guided "
          "focused ultrasound as context), `#[cfg(test)]` blocks (an "
          "`assert!` MESSAGE — \"sasp-only config has an immune effect\" — "
          "was `senescence`'s ONLY immunological match, so it was credited "
          "for a string in a test; it is credited again below, but now for "
          "`sasp_immune_mult`, a production field), and `lib.rs`, whose whole "
          "body is `pub mod` declarations. Those mentions are kept in the "
          "prose-only column rather than dropped, so both directions of error "
          "stay visible.", "",
          "**A companion audit counts the same crate and gets a different "
          "number, and the difference is the point.** "
          "`analysis/scope-audit.md` asks how many modules name *ferroptosis "
          "or the physical-ROS modalities* in production code; this asks how "
          "many name ferroptosis *chemistry*. Same "
          f"{d['module_count']} modules, same `lib.rs` "
          "and `#[cfg(test)]` exclusions, different question, so the two "
          "counts differ by exactly the modules their term sets disagree "
          "about: that audit's list carries `photodynamic`, `sonodynamic`, "
          "`pdt`, `sdt` and `photosensitiz`, which this one does not, and "
          "this one carries `acsl4`, which that one does not. Neither is "
          "wrong; a reader comparing them without this sentence would "
          "reasonably think one of them was.", "",
          "**Word boundaries are not enough on their own, and they are also "
          "the wrong boundary.** Two opposite failures, both measured, both "
          "invisible in the output:", "",
          "*Too wide, unbounded.* The first version matched `t cell`, which "
          "sits inside `mut cell`. Under the rules this document now applies "
          f"it would credit {len(d['unbounded_t_cell_credits'])} of "
          f"{d['module_count']} modules — `"
          + "`, `".join(d["unbounded_t_cell_credits"][:4])
          + "`, and so on — with modelling immunotherapy. Recomputed here "
          "rather than quoted from the version that had the defect, where it "
          "read 16 of 33 over a denominator that still counted the crate "
          "root.", "",
          "*Too wide, bounded.* `glycolytic` and `oxphos` are properly "
          "bounded, are "
          "real words, and are the names of two `Phenotype` enum variants — "
          "baseline cell states, not a therapy. Under the rules this document "
          f"now applies they would credit {len(d['phenotype_term_credits'])} "
          "of its "
          f"{d['module_count']} modules with modelling metabolic-targeting, "
          "when what the engine models is tumour metabolism as a substrate. "
          "The therapy-side terms fire nowhere, which is the correct answer. "
          "(That count is recomputed here rather than quoted from the version "
          "that had the defect, where it was eight over a denominator that "
          "still counted the crate root.)", "",
          "*Too narrow.* `_` is a word character, so `\\b` cannot match a Rust "
          "identifier: `\\bsdt\\b` misses `sdt_ros`, `\\bacsl4\\b` misses "
          "`acsl4_strength`, `\\bcd47\\b` misses `cd47_blockade`. The crate's "
          "entire naming convention was invisible, and `acsl4.rs` — the ACSL4 "
          "module — did not count as naming ferroptosis chemistry. The "
          "boundary is now \"not a letter or a digit\", which still refuses "
          "`mut cell` for `t cell` and `cd80` for `cd8`. This error ran "
          "AGAINST this document's own argument, which is the direction that "
          "keeps one unnoticed.", ""]
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
