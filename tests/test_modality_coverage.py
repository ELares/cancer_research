"""Guards for `analysis/modality-coverage.{md,json}`.

The document argues that this engine models one mechanism and that twelve of
the taxonomy's sixteen have no representation at all. An argument of that shape
is only as good as its detector, and a coverage detector fails in two opposite
directions that a count cannot tell apart:

* **Too wide** and absent mechanisms look modelled. The first draft matched the
  unbounded substring `t cell`, which sits inside `mut cell`, and credited 16
  of 33 modules -- `physics`, `stats`, `grid` -- with modelling immunotherapy.
  Nothing in the output looked wrong; the table just looked better.
* **Too narrow** and modelled mechanisms look absent, which in a document that
  proposes what to BUILD is an argument for writing code that already exists.

So every term carries a string that must match and a string that must not, and
the case table is keyed on the terms read out of the generator rather than
restated -- a table that lists its own terms cannot notice a new one.

The second half pins the prose. The report's central claim is not a count: it
is that checkpoint blockade cannot be a treatment arm because DAMPs are
proportional to lipid peroxidation at death, so `activation` is zero without
ferroptosis and every blockade setting multiplies zero. That is a claim about
Rust the report cannot check by reading itself, and it is exactly the kind of
sentence this repo has shipped wrong before, so it is measured against the
crate here.
"""
import importlib.util
import json
import re
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
MD = REPO / "analysis/modality-coverage.md"
JSON_ = REPO / "analysis/modality-coverage.json"
CORE = REPO / "simulations/ferroptosis-core/src"
SIMS = REPO / "simulations"


def _load():
    spec = importlib.util.spec_from_file_location(
        "modality_coverage", REPO / "scripts/modality_coverage.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


MC = _load()


@pytest.fixture(scope="module")
def d():
    return json.loads(JSON_.read_text())


@pytest.fixture(scope="module")
def md():
    return MD.read_text()


# ---------------------------------------------------------------- the detector

# For every term: (a string it MUST match, a string it must NOT match). The
# negatives are the point -- each is a real substring collision, most of them
# taken from Rust the crate actually contains.
CASES = {
    "immun*": ("immunity falls", "the immediate neighbour"),
    "checkpoint*": ("checkpoint brake", "check the point"),
    "pd-1": ("anti pd-1 efficacy", "pd-12 axis"),
    "pd-l1": ("dc pd-l1 protection", "pd-l10"),
    "ctla-4": ("ctla-4 brake", "ctla-40"),
    "t cell": ("one t cell kills", "let mut cell = grid"),
    "t-cell*": ("t-cells prime", "at-cell"),
    "cd8": ("cd8 infiltration", "cd80 costimulation"),
    "dendritic": ("dendritic uptake", "dendritics"[:9] + "x"),
    "icd": ("icd signal", "predicdt"),
    "treg*": ("treg suppression", "retreg"),
    "mdsc*": ("mdscs arrive", "amdsc"),
    "hdac*": ("hdac inhibitor", "shdac"),
    "histone*": ("histone marks", "ahistone"),
    "methylation": ("dna methylation", "methylations"),
    "epigenetic*": ("epigenetically locked", "aepigenetic"),
    "dnmt*": ("dnmt1 loss", "xdnmt"),
    "chromatin": ("chromatin state", "chromatins"),
    "nanoparticle*": ("nanoparticles deliver", "ananoparticle"),
    "liposom*": ("liposomal payload", "aliposome"),
    "micelle*": ("micelles form", "amicelle"),
    "nanocarrier*": ("nanocarrier load", "ananocarrier"),
    "car-t": ("car-t infusion", "car-total"),
    "car t": ("car t product", "scar t"),
    "chimeric antigen": ("chimeric antigen receptor", "chimeric antigens"),
    "glycolysis": ("aerobic glycolysis", "glycolysises"),
    "glycolytic": ("glycolytic tumors", "glycolytics"),
    "oxphos": ("oxphos cells", "oxphosx"),
    "metabolic*": ("metabolic reprogramming", "ametabolic"),
    "glutamin*": ("glutaminolysis", "aglutamine"),
    "warburg": ("warburg effect", "warburgs"),
    "lactate": ("lactate export", "lactates"),
    "antibody-drug": ("antibody-drug conjugate", "antibody-drugs"),
    "adc": ("adc payload", "adcs"),
    "payload*": ("payloads released", "apayload"),
    "parp*": ("parp inhibitor", "aparp"),
    "synthetic lethal*": ("synthetic lethality", "asynthetic lethal"),
    "brca*": ("brca1 mutant", "abrca"),
    "homologous recombination": ("homologous recombination repair",
                                 "homologous recombinations"),
    "oncolytic": ("oncolytic virus", "oncolytics"),
    "virotherapy": ("virotherapy arm", "virotherapies"),
    "adenovir*": ("adenovirus vector", "aadenovirus"),
    "crispr": ("crispr screen", "crisprs"),
    "cas9": ("cas9 nuclease", "cas90"),
    "guide rna": ("guide rna pool", "guide rnas"),
    "sgrna*": ("sgrnas designed", "asgrna"),
    "bispecific*": ("bispecific antibody", "abispecific"),
    "bite": ("bite construct", "bites"),
    "t-cell engager*": ("t-cell engagers", "at-cell engager"),
    "electroporation": ("irreversible electroporation", "electroporations"),
    "electrochemical": ("electrochemical therapy", "electrochemicals"),
    "irreversible electro*": ("irreversible electroporation",
                              "airreversible electro"),
    "sonodynamic": ("sonodynamic therapy", "sonodynamics"),
    "sdt": ("sdt arm", "sdts"),
    "ultrasound": ("ultrasound pulse", "ultrasounds"),
    "sonosensitiz*": ("sonosensitizer dose", "asonosensitizer"),
    "hifu": ("hifu ablation", "hifus"),
    "focused ultrasound": ("focused ultrasound beam", "focused ultrasounds"),
    "thermal ablation": ("thermal ablation zone", "thermal ablations"),
    "cd47": ("cd47 blockade", "cd470"),
    "sirp*": ("sirpalpha", "asirp"),
    "phagocytos*": ("phagocytosis checkpoint", "aphagocytosis"),
    "microbiome": ("gut microbiome", "microbiomes"),
    "microbiota": ("microbiota shift", "microbiotas"),
    "bacteri*": ("bacterial vector", "abacteria"),
    "mrna vaccine": ("mrna vaccine arm", "mrna vaccines"),
    "neoantigen*": ("neoantigens presented", "aneoantigen"),
    "lipid nanoparticle vaccine": ("lipid nanoparticle vaccine dose",
                                   "lipid nanoparticle vaccines"),
    "ferropto*": ("ferroptotic death", "aferroptosis"),
    "gpx4": ("gpx4 inhibition", "gpx40"),
    "lipid perox*": ("lipid peroxidation", "alipid perox"),
    "fsp1": ("fsp1 rescue", "fsp10"),
    "acsl4": ("acsl4 high", "acsl40"),
    "slc7a11": ("slc7a11 export", "slc7a110"),
    "rsl3": ("rsl3 dose", "rsl30"),
    "erastin*": ("erastin arm", "aerastin"),
}

ALL_TERMS = sorted({t for ts in MC.ENGINE_TERMS.values() for t in ts}
                   | set(MC.FERROPTOSIS_TERMS))


def test_every_term_the_generator_uses_has_a_case_in_both_directions():
    """A case table that lists its own terms cannot notice a new one, and an
    entry emptied on one side passes a coverage check that compares keys."""
    assert set(CASES) == set(ALL_TERMS), (
        f"cases missing for {sorted(set(ALL_TERMS) - set(CASES))}; "
        f"cases for terms that no longer exist: "
        f"{sorted(set(CASES) - set(ALL_TERMS))}")
    for term, (yes, no) in CASES.items():
        assert yes.strip() and no.strip(), f"{term} has an empty case"


@pytest.mark.parametrize("term", ALL_TERMS)
def test_each_term_matches_its_positive_and_misses_its_negative(term):
    yes, no = CASES[term]
    pat = MC.term_pattern(term)
    assert re.search(pat, yes), f"{term!r} ({pat}) failed to match {yes!r}"
    assert not re.search(pat, no), f"{term!r} ({pat}) wrongly matched {no!r}"


def test_the_motivating_collision_stays_fixed():
    """`t cell` inside `mut cell` credited 16 of 33 modules with immunotherapy.

    Pinned on the real Rust rather than on the abstract rule, and pinned in
    both directions: the fix must still find an actual T cell.
    """
    assert not MC._matches("let mut cell = grid.cells[idx];", ("t cell",))
    assert MC._matches("each t cell rolls once", ("t cell",))


def test_a_stem_and_a_bounded_term_differ():
    """`*` is the only thing separating the two rules, so it must do work."""
    assert MC._matches("immunotherapy arm", ("immun*",))
    assert not MC._matches("immunotherapy arm", ("immun",))


def test_a_stem_must_stop_before_the_forms_diverge():
    """Two stems shipped one letter too long and each missed the commonest
    form in the crate: `ferroptos*` does not reach `ferroptotic`, and
    `immune*` reaches neither `immunity` nor `immunotherapy`. One measured
    form is not enough to place a stem's boundary."""
    for form in ("ferroptosis", "ferroptotic", "ferroptotically"):
        assert MC._matches(form, ("ferropto*",)), form
        assert MC._matches(form, MC.FERROPTOSIS_TERMS), form
    assert not MC._matches("ferroptotic", ("ferroptos*",))
    for form in ("immunity", "immunotherapy", "immunogenic", "immune"):
        assert MC._matches(form, MC.ENGINE_TERMS["immunotherapy"]), form
    assert not MC._matches("immunity", ("immune*",))


# ------------------------------------------------------- comments are not code

def test_comments_are_stripped_and_code_is_not():
    """Both directions. Stripping too much would empty the detector silently,
    which reads as `the engine models nothing` -- the direction that argues
    for building things twice."""
    src = "\n".join([
        "/// HIFU is mentioned here as context only.",
        "//! module doc naming crispr",
        "/* block comment naming parp */",
        "pub struct Checkpoint { pub brake: f64 }",
        "let x = 1; // trailing comment naming cd47",
    ])
    out = MC._strip_comments(src).lower()
    for gone in ("hifu", "crispr", "parp", "cd47"):
        assert gone not in out, f"{gone} survived comment stripping"
    assert "pub struct checkpoint" in out
    assert "let x = 1;" in out


def test_hifu_is_prose_only_and_the_report_says_why(d, md):
    """The measured case that motivated the tier. If a HIFU model ever lands,
    this fails and the paragraph explaining the tier must be rewritten with
    it -- which is the intent."""
    row = next(r for r in d["rows"] if r["mechanism"] == "hifu")
    assert row["engine_tier"] == "absent", (
        "hifu is no longer prose-only; the report's worked example for the "
        "PROSE-ONLY tier is stale and must be replaced")
    assert row["prose_only_modules"], (
        "hifu has no prose mention either, so the report's example is false")
    assert "prose about a thing is not a model of it" in md


# ----------------------------------------------- the claims the prose makes

def test_the_treatment_variants_are_read_not_asserted(d):
    """The headline count comes from `cell.rs`, so a new arm must move it."""
    src = (CORE / "cell.rs").read_text()
    body = re.search(r"pub enum Treatment\s*\{(.*?)\}", src, re.S).group(1)
    named = [v.strip().rstrip(",") for v in body.splitlines()
             if v.strip() and not v.strip().startswith("//")]
    assert d["treatment_variants"] == named
    assert f"**{len(named)} treatments**" in MD.read_text()


def test_no_treatment_variant_is_a_non_ros_modality(d):
    """The report's headline -- every arm is ferroptosis or physical-ROS -- is
    a claim about the enum, so it is checked against the enum."""
    tiers = {r["mechanism"]: r["engine_tier"] for r in d["rows"]}
    modelled = {m for m, t in tiers.items() if t == "treatment"}
    assert modelled <= {"sonodynamic"}, (
        f"{sorted(modelled)} are now treatment arms; the report's "
        "'every one of them ferroptosis or physical-ROS' sentence and its "
        "immunotherapy paragraph both need re-deriving")


def _damp_writes():
    """Every `damp_field[...] +=` in the simulation binaries."""
    out = []
    for p in sorted(SIMS.glob("sim-*/src/*.rs")):
        for i, line in enumerate(p.read_text().splitlines(), 1):
            stripped = line.split("//")[0]
            if re.search(r"damp_field\s*\[[^\]]+\]\s*\+=", stripped):
                out.append((p.name, i, stripped.strip()))
    return out


def test_every_damp_source_is_lipid_peroxidation_at_death():
    """The report's central claim, and the reason immunotherapy is a MODIFIER.

    If any writer ever adds a ferroptosis-independent antigen source, the
    paragraph asserting blockade multiplies zero becomes false -- and it is
    the paragraph the whole 'what to build next' argument rests on.
    """
    writes = _damp_writes()
    assert writes, "no DAMP writer found; the scan is broken, not the engine"
    for name, line_no, text in writes:
        assert "lp_at_grace_end" in text, (
            f"{name}:{line_no} adds DAMPs without lipid peroxidation at "
            f"death: {text!r}. The report claims there is no "
            "ferroptosis-independent antigen source in the engine.")


def test_activation_is_zero_without_damp_and_gates_the_kill():
    """`damp/(damp+kd)` and `activation · rate · (1 − brake)`, both quoted in
    the report, both read out of the crate."""
    src = (CORE / "immune_spatial.rs").read_text()
    act = re.search(r"pub fn dc_activation\([^)]*\)\s*->\s*f64\s*\{(.*?)\n\}",
                    src, re.S).group(1)
    assert "local_damp / (local_damp + kd)" in act, (
        f"dc_activation is no longer damp/(damp+kd): {act.strip()!r}")
    kill = re.search(
        r"pub fn immune_kill_probability\([^)]*\)\s*->\s*f64\s*\{(.*?)\n\}",
        src, re.S).group(1)
    assert re.search(r"activation\s*\*\s*kill_rate\s*\*\s*\(1\.0\s*-\s*effective_brake\)",
                     kill), (
        f"immune_kill_probability no longer multiplies activation: {kill.strip()!r}")


# ----------------------------------------------------------- the refusals

def test_the_report_refuses_a_volume_ranking(md):
    """Inherited verbatim from `census-mechanism-profile.md`, because the
    reason is the same one and a table in volume order invites the reading it
    refuses."""
    assert "Volume is NOT comparable across mechanisms" in md
    assert "how broad each descriptor is" in md


def test_the_report_names_what_the_taxonomy_cannot_see(md):
    """A mechanism with no MeSH descriptor is absent rather than zero, and
    TTFields is FDA-approved in two indications."""
    assert "Tumour-treating fields" in md or "tumour-treating fields" in md
    assert "no descriptor" in md
    assert "Radiotherapy" in md


def test_the_totals_are_the_rows(d):
    """Every headline number is arithmetic on the table beside it."""
    absent = [r for r in d["rows"] if r["engine_tier"] == "absent"]
    assert d["absent_count"] == len(absent)
    assert d["absent_census"] == sum(r["census"] for r in absent)
    assert d["absent_trials"] == sum(r["trials"] for r in absent)
    assert d["total_census"] == sum(r["census"] for r in d["rows"])
    assert f"**{d['absent_count']} of {len(d['rows'])} mechanisms" in MD.read_text()


def test_the_rows_match_the_census_profile(d):
    """The counts are copied, not recomputed, so they must still be copies."""
    prof = json.loads((REPO / "analysis/census-mechanism-profile.json").read_text())
    src = {r["mechanism"]: r for r in prof["rows"]}
    assert d["census"] == prof["census"]
    assert {r["mechanism"] for r in d["rows"]} == set(src)
    for r in d["rows"]:
        for k in ("census", "trials", "trial_share"):
            assert r[k] == src[r["mechanism"]][k], f"{r['mechanism']}.{k} drifted"
