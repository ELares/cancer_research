"""Guards for `analysis/modality-coverage.{md,json}`.

The document argues that this engine models one mechanism and that thirteen of
the sixteen it can count have no representation at all. An argument of that
shape is only as good as its detector, and a coverage detector fails in two
opposite directions that a count cannot tell apart:

* **Too wide** and absent mechanisms look modelled. Three measured instances,
  each caught by a different reviewer: the unbounded substring `t cell` sits
  inside `mut cell` and credited 16 of 33 modules -- `physics`, `stats`,
  `grid` -- with modelling immunotherapy; `glycolytic` and `oxphos` are
  properly bounded real words that are also `Phenotype` ENUM VARIANTS, and
  credited eight modules with modelling metabolic targeting; and an `assert!`
  MESSAGE inside a `#[cfg(test)]` block credited `senescence` with the same.
* **Too narrow** and modelled mechanisms look absent, which in a document that
  describes what the engine CANNOT express is an argument for writing code
  that already exists.

So every term carries a string that must match and a string that must not, the
case table is keyed on the terms read out of the generator rather than
restated, and the tier assertions are EQUALITIES rather than bounds -- the
first version asserted `modelled <= {"sonodynamic"}`, which the empty set
satisfies, so deleting the entire TREATMENT tier left all ninety guards green
and shipped a document calling sonodynamic a modifier.

The second half pins the prose. The report's central claim is not a count: it
is that checkpoint blockade cannot be a treatment arm because every DAMP source
is proportional to lipid peroxidation at death, so the activation term is zero
without ferroptosis and every blockade setting multiplies zero. That is a claim
about Rust the report cannot check by reading itself.
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
# negatives are the point -- each is a real substring collision, and each
# CONTAINS the term's text, so a matcher that dropped its boundaries entirely
# would fail here. Three of the first version's negatives did not contain the
# term at all and so proved nothing.
CASES = {
    "immun*": ("immunity falls", "a preimmunisation step"),
    "checkpoint*": ("checkpoint brake", "the recheckpoint pass"),
    "pd-1": ("anti pd-1 efficacy", "pd-1x axis"),
    "pd-l1": ("dc pd-l1 protection", "pd-l1x"),
    "ctla-4": ("ctla-4 brake", "ctla-4x"),
    "t cell": ("one t cell kills", "let mut cell = grid"),
    "t-cell*": ("t-cells prime", "at-cell"),
    "cd8": ("cd8 infiltration", "cd80 costimulation"),
    "dendritic": ("dendritic uptake", "adendritic"),
    "icd": ("icd signal", "aicd"),
    "treg*": ("treg suppression", "atreg"),
    "mdsc*": ("mdscs arrive", "amdsc"),
    "hdac*": ("hdac inhibitor", "shdac"),
    "histone*": ("histone marks", "ahistone"),
    "methylation": ("dna methylation", "amethylation"),
    "epigenetic*": ("epigenetically locked", "aepigenetic"),
    "dnmt*": ("dnmt1 loss", "xdnmt"),
    "chromatin": ("chromatin state", "achromatin"),
    "nanoparticle*": ("nanoparticles deliver", "ananoparticle"),
    "liposom*": ("liposomal payload", "aliposome"),
    "micelle*": ("micelles form", "amicelle"),
    "nanocarrier*": ("nanocarrier load", "ananocarrier"),
    "car-t": ("car-t infusion", "car-tx"),
    "car t": ("car t product", "scar t"),
    "chimeric antigen": ("chimeric antigen receptor", "chimeric antigenx"),
    "glycolysis": ("aerobic glycolysis", "aglycolysis"),
    "metabolic*": ("metabolic reprogramming", "ametabolic"),
    "glutamin*": ("glutaminolysis", "aglutamine"),
    "warburg": ("warburg effect", "awarburg"),
    "lactate": ("lactate export", "alactate"),
    "2-deoxyglucose": ("2-deoxyglucose dose", "x2-deoxyglucose"),
    "dichloroacetate": ("dichloroacetate arm", "adichloroacetate"),
    "antibody-drug": ("antibody-drug conjugate", "antibody-drugx"),
    "adc": ("adc payload", "adcx"),
    "payload*": ("payloads released", "apayload"),
    "parp*": ("parp inhibitor", "aparp"),
    "synthetic lethal*": ("synthetic lethality", "asynthetic lethal"),
    "brca*": ("brca1 mutant", "abrca"),
    "homologous recombination": ("homologous recombination repair",
                                 "ahomologous recombination"),
    "oncolytic": ("oncolytic virus", "aoncolytic"),
    "virotherapy": ("virotherapy arm", "avirotherapy"),
    "adenovir*": ("adenovirus vector", "aadenovirus"),
    "crispr": ("crispr screen", "acrispr"),
    "cas9": ("cas9 nuclease", "cas9x"),
    "guide rna": ("guide rna pool", "guide rnax"),
    "sgrna*": ("sgrnas designed", "asgrna"),
    "bispecific*": ("bispecific antibody", "abispecific"),
    "bite": ("bite construct", "bitex"),
    "t-cell engager*": ("t-cell engagers", "at-cell engager"),
    "electroporation": ("irreversible electroporation", "aelectroporation"),
    "electrochemical": ("electrochemical therapy", "aelectrochemical"),
    "irreversible electro*": ("irreversible electroporation",
                              "airreversible electro"),
    "sonodynamic": ("sonodynamic therapy", "asonodynamic"),
    "sdt": ("sdt arm", "sdtx"),
    "ultrasound": ("ultrasound pulse", "aultrasound"),
    "sonosensitiz*": ("sonosensitizer dose", "asonosensitizer"),
    "hifu": ("hifu ablation", "hifux"),
    "focused ultrasound": ("focused ultrasound beam", "afocused ultrasound"),
    "thermal ablation": ("thermal ablation zone", "athermal ablation"),
    "cd47": ("cd47 blockade", "cd47x"),
    "sirp*": ("sirpalpha", "asirp"),
    "phagocytos*": ("phagocytosis checkpoint", "aphagocytosis"),
    "microbiome": ("gut microbiome", "amicrobiome"),
    "microbiota": ("microbiota shift", "amicrobiota"),
    "bacteri*": ("bacterial vector", "abacteria"),
    "mrna vaccine": ("mrna vaccine arm", "amrna vaccine"),
    "neoantigen*": ("neoantigens presented", "aneoantigen"),
    "lipid nanoparticle vaccine": ("lipid nanoparticle vaccine dose",
                                   "alipid nanoparticle vaccine"),
    "ferropto*": ("ferroptotic death", "aferroptosis"),
    "gpx4": ("gpx4 inhibition", "gpx4x"),
    "lipid perox*": ("lipid peroxidation", "alipid perox"),
    "fsp1": ("fsp1 rescue", "fsp1x"),
    "acsl4": ("acsl4 high", "acsl4x"),
    "slc7a11": ("slc7a11 export", "slc7a11x"),
    "rsl3": ("rsl3 dose", "rsl3x"),
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


def test_every_negative_case_actually_contains_the_term():
    """The negatives must probe the BOUNDARY, not the word.

    Three of the first version's -- `immun*`/"the immediate neighbour",
    `checkpoint*`/"check the point", `virotherapy`/"virotherapies" -- did not
    contain the term's characters at all, so a matcher with its boundaries
    stripped out entirely would still have passed them.
    """
    for term, (_, no) in CASES.items():
        body = term[:-1] if term.endswith("*") else term
        assert body in no, (
            f"{term!r}'s negative {no!r} does not contain {body!r}, so it "
            "cannot detect a matcher that lost its boundaries")


@pytest.mark.parametrize("term", ALL_TERMS)
def test_each_term_matches_its_positive_and_misses_its_negative(term):
    yes, no = CASES[term]
    pat = MC.term_pattern(term)
    assert re.search(pat, yes), f"{term!r} ({pat}) failed to match {yes!r}"
    assert not re.search(pat, no), f"{term!r} ({pat}) wrongly matched {no!r}"


def test_the_motivating_collision_stays_fixed():
    """`t cell` inside `mut cell` credited 16 of 33 modules with immunotherapy.

    Pinned on the real Rust rather than on the abstract rule, and in both
    directions: the fix must still find an actual T cell.
    """
    assert not MC._matches("let mut cell = grid.cells[idx];", ("t cell",))
    assert MC._matches("each t cell rolls once", ("t cell",))


def test_a_stem_must_stop_before_the_forms_diverge():
    """Two stems shipped one letter too long. `ferroptos*` does not reach
    `ferroptotic`, the commonest form in this crate, and `immune*` reaches
    neither `immunity` nor `immunotherapy` nor `immunogenic`. One measured form
    is not enough to place a stem's boundary."""
    for form in ("ferroptosis", "ferroptotic", "ferroptotically"):
        assert MC._matches(form, ("ferropto*",)), form
        assert MC._matches(form, MC.FERROPTOSIS_TERMS), form
    assert not MC._matches("ferroptotic", ("ferroptos*",))
    for form in ("immunity", "immunotherapy", "immunogenic", "immune"):
        assert MC._matches(form, MC.ENGINE_TERMS["immunotherapy"]), form
    assert not MC._matches("immunity", ("immune*",))


def test_a_phenotype_is_not_a_therapy():
    """`glycolytic` and `oxphos` are bounded, real, and the wrong thing.

    They name two `Phenotype` enum variants -- baseline cell states -- and
    credited eight modules with modelling metabolic TARGETING. Word boundaries
    do not help against a real word used for something else, so this is pinned
    on the enum itself.
    """
    src = MC.strip_rust_comments((CORE / "cell.rs").read_text()).lower()
    body = re.search(r"pub enum phenotype\s*\{(.*?)\}", src, re.S).group(1)
    assert "glycolytic" in body and "oxphos" in body, (
        "the Phenotype enum no longer names these; re-check whether the "
        "exclusion still makes sense")
    terms = MC.ENGINE_TERMS["metabolic-targeting"]
    for state in ("glycolytic", "oxphos"):
        assert not MC._matches(state, terms), (
            f"{state} is a cell phenotype and is scoring as a therapy again")
    # And the therapy side must still work.
    assert MC._matches("2-deoxyglucose blocks glycolysis", terms)


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
    out = MC.strip_rust_comments(src).lower()
    for gone in ("hifu", "crispr", "parp", "cd47"):
        assert gone not in out, f"{gone} survived comment stripping"
    assert "pub struct checkpoint" in out
    assert "let x = 1;" in out


def test_the_stripper_handles_the_three_constructs_a_regex_cannot():
    """All three are legal Rust, all three are latent in this crate today, and
    the first version got all three wrong. Planted deliberately, because
    waiting for one to appear means shipping the wrong count first."""
    # Nested block comments: a non-greedy regex closes at the INNER `*/` and
    # leaves the outer comment's tail behind as code.
    nested = "/* outer /* inner */ hdac_inhibitor in prose */ let a = 1;"
    out = MC.strip_rust_comments(nested)
    assert "hdac_inhibitor" not in out, f"nested comment leaked: {out!r}"
    assert "let a = 1;" in out
    # `//` inside a string literal must not truncate the rest of the line.
    s = 'let doc = "see http://x"; let p = Checkpoint::Pd1;'
    out = MC.strip_rust_comments(s)
    assert "Checkpoint::Pd1" in out, f"string-literal // truncated code: {out!r}"
    # `/*` inside a string literal must not open a comment.
    s = 'let a = "/*"; let panel = CheckpointPanel::new(); let b = "*/";'
    out = MC.strip_rust_comments(s)
    assert "CheckpointPanel::new()" in out, f"string-literal /* ate code: {out!r}"
    # A lifetime is not a char literal.
    s = "fn f<'a>(x: &'a str) -> Checkpoint { todo!() }"
    assert "Checkpoint" in MC.strip_rust_comments(s)
    # A raw string survives intact.
    s = 'let j = r#"{"sdt_ros": 5.0}"#; let q = 1;'
    assert "let q = 1;" in MC.strip_rust_comments(s)


def test_test_blocks_are_stripped_and_the_measured_case_stays_fixed():
    """An `assert!` MESSAGE credited `senescence` with modelling immunotherapy.

    Pinned twice: on a synthetic block, and on the real file, so the exclusion
    cannot quietly stop applying.
    """
    src = "\n".join([
        "pub fn real() -> f64 { 1.0 }",
        "#[cfg(test)]",
        "mod tests {",
        "    #[test]",
        "    fn t() { assert!(x, \"has an immune effect {y}\"); }",
        "}",
        "pub fn after() -> f64 { 2.0 }",
    ])
    out = MC.strip_test_blocks(src)
    assert "immune effect" not in out, f"test block leaked: {out!r}"
    assert "pub fn real()" in out and "pub fn after()" in out, (
        "stripping ate code outside the test block")
    # The real file. `senescence` IS credited with immunotherapy, and after
    # this fix it is credited for the right reason: `sasp_immune_mult` and
    # `sasp_immune_multiplier` are production identifiers for a SASP field that
    # modulates immune killing. What must NOT be there is the assertion message
    # that credited it before, so both halves are pinned -- dropping the
    # test-block strip would leave this test green on the first assertion
    # alone.
    code, _ = MC._module_text()
    sen = code["senescence"]
    assert "sasp_immune_mult" in sen, (
        "the production identifier is gone; re-derive why senescence is "
        "credited with immunotherapy before trusting the tier")
    assert "sasp-only config has an immune effect" not in sen, (
        "the #[cfg(test)] assertion message is back in the scanned code")
    assert "#[test]" not in sen and "assert!(" not in sen, (
        "test code survived stripping in a real module")


def test_the_crate_root_is_not_a_module(d):
    """`lib.rs`'s whole body is `pub mod <name>;`, so `pub mod immune;`
    credited it with modelling immunotherapy. `scope_audit.py` excludes it for
    the same reason, and the two denominators must agree."""
    assert "lib" in MC.NOT_A_MODULE
    assert "lib" not in d["ferroptosis_modules"]
    for r in d["rows"]:
        assert "lib" not in r["code_modules"], r["mechanism"]
    live = len([p for p in CORE.glob("*.rs") if p.stem != "lib"])
    assert d["module_count"] == live == 33, (
        f"module count {d['module_count']} against {live} on disk")


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
             if v.strip() and not v.strip().startswith(("//", "#["))]
    assert d["treatment_variants"] == named
    assert d["active_treatments"] == [v for v in named if v != "Control"]
    assert f"**{len(d['active_treatments'])} treatments**" in MD.read_text()
    # Every arm must be CLASSIFIED, so the next one cannot land under a
    # headline sentence that no longer describes it. `Radiation` made the old
    # blanket "every one of them ferroptosis or physical-ROS" false the moment
    # it landed, and nothing caught it but this file's own regeneration.
    for v in d["active_treatments"]:
        assert MC.TREATMENT_KIND.get(v), (
            f"{v} has no entry in TREATMENT_KIND, so the headline sentence "
            "says nothing about it")
        assert d["treatment_kinds"][v] == MC.TREATMENT_KIND[v]
        assert f"`{v}`" in MD.read_text()
    assert set(MC.TREATMENT_KIND) == set(d["active_treatments"]), (
        "TREATMENT_KIND classifies arms that do not exist, or misses one")


def test_the_control_arm_is_not_counted_as_a_treatment(d, md):
    """`Control` applies nothing -- `physics.rs` returns 0.0 for it and
    `biochem.rs` gives it no exogenous ROS -- so calling it one made the
    headline read `4 treatments, every one of them ferroptosis or
    physical-ROS` when the fourth is neither."""
    assert "Control" in d["treatment_variants"]
    assert "Control" not in d["active_treatments"]
    phys = MC.strip_rust_comments((CORE / "physics.rs").read_text())
    assert re.search(r"Treatment::Control\s*=>\s*0\.0", phys), (
        "Control now does something; the headline sentence needs re-deriving")
    assert "untreated `Control` arm that applies nothing" in md


def test_an_attribute_line_is_not_a_variant():
    """`#[default]` above a variant is legal Rust and was being counted, so a
    derive change would have moved the headline count."""
    assert MC._treatment_variants(), "the enum parse returned nothing"
    for v in MC._treatment_variants():
        assert not v.startswith("#["), v


def test_no_treatment_variant_is_a_non_ros_modality(d):
    """An EQUALITY, not a bound. The first version asserted
    `modelled <= {"sonodynamic"}`, which the empty set satisfies: deleting the
    whole TREATMENT tier left every guard green and shipped a document calling
    sonodynamic a modifier."""
    tiers = {r["mechanism"]: r["engine_tier"] for r in d["rows"]}
    modelled = {m for m, t in tiers.items() if t == "treatment"}
    assert modelled == {"sonodynamic"}, (
        f"treatment tier is {sorted(modelled)}, expected exactly "
        "{'sonodynamic'}. If an arm was added, the report's 'every one of them "
        "ferroptosis or physical-ROS' sentence and its immunotherapy paragraph "
        "both need re-deriving; if the tier emptied, the detector is broken.")


def test_every_taxonomy_mechanism_has_terms():
    """A mechanism with no entry scores `absent` through an `any()` over an
    empty tuple, silently inflating the headline."""
    prof = json.loads((REPO / "analysis/census-mechanism-profile.json").read_text())
    for r in prof["rows"]:
        assert MC.ENGINE_TERMS.get(r["mechanism"]), (
            f"{r['mechanism']} has no ENGINE_TERMS entry, so it scores absent "
            "without anything being measured")


def test_the_immune_claim_is_proved_in_rust_and_its_limit_is_stated(md):
    """No Python guard here decides the central claim, and the report says so.

    Three successive scans tried and each was defeated by an ordinary Rust
    idiom that put the mutation and the field name on different lines: a
    helper function, an `iter_mut().for_each`, a plain `for` loop. A guard
    recording a property its scan cannot DECIDE is worse than no guard, so the
    scan was removed rather than widened a fourth time.

    What replaces it is a behaviour test in the crate, and this asserts the
    report is honest about its scope rather than re-implementing the scan.
    """
    assert "proved END TO END in Rust" in md
    assert "no test here covers" in md
    assert "sim-tme-3d`'s own loops" in md
    assert "worse than no guard" in md
    # And the scan really is gone -- not merely unused.
    src = (REPO / "scripts/modality_coverage.py").read_text()
    for gone in ("_damp_surface", "SURFACE_DIGEST", "_FRESH", "damp_writers"):
        assert gone not in src, f"{gone} is still in the generator"
    d = json.loads(JSON_.read_text())
    assert "damp_surface" not in d and "damp_writers" not in d


def test_the_end_to_end_rust_proof_exists_and_carries_a_positive_case(md):
    """The proof must assert BOTH directions.

    A test that only asserts zero kills is satisfied by an engine that does
    nothing, so the Rust test pairs the no-death case with a case where the
    same chain DOES kill and blockade DOES move the number. Both halves are
    required here, keyed on the assertions rather than on the function
    existing -- an earlier version of this guard checked only that a `for`
    loop was present, and passed against a body that asserted nothing.
    """
    src = (CORE / "immune.rs").read_text()
    body = re.search(
        r"fn no_ferroptotic_death_end_to_end_means_no_kills_at_any_blockade\(\)"
        r"\s*\{(.*?)\n    }", src, re.S)
    assert body, "the end-to-end proof is gone"
    b = body.group(1)
    assert "death_threshold: f64::INFINITY" in b, (
        "the negative case no longer makes death impossible, so it proves "
        "nothing about ferroptosis")
    assert 'assert_eq!(\n                        kills, 0.0' in b or "kills, 0.0" in b
    assert "dead > 0" in b, "the positive case no longer requires real deaths"
    assert "kills_on > kills_off" in b, (
        "the positive case no longer requires blockade to move the number, so "
        "the negative case could pass on a disconnected chain")
    assert "kills_on < survivors" in b, (
        "the survivor cap is no longer excluded, and at the default rates "
        "both arms saturate at the same number")


def test_both_immune_models_are_named_and_each_gate_is_the_right_one(d, md):
    """The first version named five modules as reaching T-cell killing and four
    of them did not -- `lib.rs` on a `pub mod` line, `senescence` on two test
    assertion messages, `params.rs` on a struct, and `immune.rs`, which is a
    SECOND immune model with its own DC activation and PD-1 brake rather than a
    caller of the first."""
    models = {m["kill_fn"]: m for m in d["immune_models"]}
    assert set(models) == {"immune_kill_probability", "immune_cascade"}
    assert models["immune_kill_probability"]["callers"] == ["sim-tme", "sim-tme-3d"]
    assert models["immune_cascade"]["callers"] == ["sim-combo", "sim-icd"]
    assert "2 immune kill paths" in md
    # The two gates are DIFFERENT and an earlier draft conflated them. The
    # well-mixed model multiplies by `n_dead`, so it is gated on the COUNT of
    # ferroptotic deaths, not on DAMPs -- verified by mutation: injecting a
    # constant `damp_per_dead_cell` produces no kills, injecting a constant
    # `mature_dcs` does.
    assert "gate on ferroptotic death in DIFFERENT ways" in md

    spatial = MC.strip_rust_comments((CORE / "immune_spatial.rs").read_text())
    assert "local_damp / (local_damp + kd)" in spatial, (
        "dc_activation is no longer damp/(damp+kd)")
    assert re.search(
        r"activation\s*\*\s*kill_rate\s*\*\s*\(1\.0\s*-\s*effective_brake\)",
        spatial), "immune_kill_probability no longer multiplies activation"

    wellmixed = MC.strip_rust_comments((CORE / "immune.rs").read_text())
    assert re.search(
        r"damp_per_dead_cell\s*/\s*\(damp_per_dead_cell\s*\+\s*params\.dc_activation_kd\)",
        wellmixed), "immune_cascade's DC activation changed shape"
    assert re.search(r"primed_tcells\s*\*\s*params\.tcell_kill_rate\s*\*\s*kill_efficiency",
                     wellmixed), "immune_cascade no longer gates kills on priming"
    assert re.search(r"mature_dcs\s*=\s*dc_activation_fraction\s*\*\s*"
                     r"params\.dc_maturation_rate\s*\*\s*n_dead as f64",
                     wellmixed), (
        "the well-mixed model no longer multiplies by n_dead, so the report's "
        "'gated on the COUNT' sentence needs re-deriving")
    # And priming is gated on dead cells in BOTH: an empty death list must
    # produce no DAMPs at all.
    assert "if dead_cell_lps.is_empty()" in wellmixed


# ----------------------------------------------------------- the refusals

def test_the_report_refuses_a_volume_ranking(md):
    """Inherited from `census-mechanism-profile.md`, INCLUDING the clause that
    bites -- the profile refuses a cross-mechanism RANKING, and a table sorted
    by census invites exactly that reading."""
    assert "not a ranking" in md
    assert "how broad each descriptor is" in md
    assert "for legibility only" in md
    assert "the content is the ENGINE column" in md
    # The source's own refusal must still be there to inherit.
    profile_md = (REPO / "analysis/census-mechanism-profile.md").read_text()
    assert "how broad each descriptor is" in profile_md


def test_the_report_counts_what_the_taxonomy_cannot_measure(d, md):
    """A mechanism with no MeSH descriptor is absent rather than zero, and the
    first version called the sixteen rows `this repo's own mechanism list`
    when the taxonomy names twenty-five."""
    import yaml
    mp = yaml.safe_load((REPO / "analysis/mesh-mechanism-map.yaml").read_text())
    named = set(mp.get("mechanisms", {})) | set(mp.get("unmeasurable", {}))
    assert set(d["taxonomy_named"]) == named
    assert len(named) > len(d["rows"]), (
        "every named mechanism now has a row; the paragraph about the "
        "unmeasurable ones is stale")
    assert set(d["taxonomy_without_a_row"]) == named - {r["mechanism"] for r in d["rows"]}
    assert f"taxonomy names {len(named)} mechanisms" in md
    for m in ("ttfields", "frequency-therapy"):
        assert m in d["taxonomy_without_a_row"], m
        assert f"`{m}`" in md
    assert "cold-atmospheric-plasma HAS a descriptor" in md
    assert "counted elsewhere in this repository" in md


def test_the_radiotherapy_figures_name_their_denominator(md):
    """88 is the count TITLED about radiotherapy, not its corpus presence, and
    the tag counts overlap. Both were stated loosely, and both are checked here
    against the artifact they come from."""
    assert "whose TITLES are about it" in md
    assert "the tags overlap" in md
    src = (REPO / "analysis/atlas-untagged-partner.md").read_text()
    for frag in ("| radiotherapy | 88 | 13 (14.8%) | **75** (85.2%) | 24 (27.3%) |",
                 "`immunotherapy` 45, `ttfields` 11, `nanoparticle` 7, "
                 "`bispecific-antibody` 6"):
        assert frag in src, f"the source no longer says {frag!r}"
    assert "| radiotherapy | 88 (1.8%) | 907 (18.8%) | 2,436 (50.4%) |" in src


def test_the_totals_are_the_rows(d):
    """Every headline number is arithmetic on the table beside it."""
    absent = [r for r in d["rows"] if r["engine_tier"] == "absent"]
    assert d["absent_count"] == len(absent)
    assert d["absent_census"] == sum(r["census"] for r in absent)
    assert d["absent_trials"] == sum(r["trials"] for r in absent)
    assert d["total_census"] == sum(r["census"] for r in d["rows"])
    assert f"**{d['absent_count']} of {len(d['rows'])} mechanisms" in MD.read_text()


def test_the_ferroptosis_baseline_is_guarded(d, md):
    """The headline's other half. Nothing checked it before, and four of the
    fifteen it then claimed were `lib.rs` or test-only matches."""
    code, _ = MC._module_text()
    live = sorted(n for n, t in code.items() if MC._matches(t, MC.FERROPTOSIS_TERMS))
    assert d["ferroptosis_modules"] == live
    assert f"**{len(live)} of its {d['module_count']} modules**" in md
    # `io` and `slab` matched only inside `#[cfg(test)]` -- a `"rsl3"` CSV
    # fixture and a `gpx4` field name in a byte-identity assert -- and must
    # stay out. `lib.rs` is the crate root and is excluded upstream. This is
    # the same exclusion `scope_audit.py` documents for itself.
    for excluded in ("io", "slab", "lib"):
        assert excluded not in live, (
            f"{excluded} is credited with ferroptosis code again; check "
            "whether the match is production code or a test fixture")
    # And the underscore-boundary fix must still be working: these four are
    # credited ONLY through snake_case identifiers, which `\b` refused.
    for snake in ("acsl4", "copper", "ifngamma", "trigger_wave"):
        assert snake in live, (
            f"{snake} is uncredited again -- the term boundary has stopped "
            "matching Rust identifiers like `acsl4_strength` / `gpx4_defense`")


def test_the_two_module_counts_are_reconciled_not_left_to_collide(d, md):
    """Two front-door artifacts count the same crate and publish different
    numbers, and until this test neither mentioned the other.

    They are genuinely comparable -- same 32 modules, same `lib.rs` and
    `#[cfg(test)]` exclusions -- so the gap has to be the term sets and
    nothing else. Recomputed live from `scope_audit`'s own regex rather than
    read from either artifact, so a change to that regex fails here instead of
    silently widening the gap.
    """
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "scope_audit", REPO / "scripts/scope_audit.py")
    sa = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(sa)

    theirs = []
    for f in sorted(CORE.glob("*.rs")):
        if f.stem == "lib":
            continue
        t = f.read_text()
        i = t.find("#[cfg(test)]")
        prod = t if i < 0 else t[:i]
        code = "\n".join(ln for ln in prod.split("\n")
                         if not ln.strip().startswith(("//", "*", "/*")))
        if sa.FERRO.search(code):
            theirs.append(f.stem)
    mine = set(d["ferroptosis_modules"])
    assert set(theirs) - mine == {"photosensitizer_pk"}, (
        f"scope_audit now counts {sorted(set(theirs) - mine)} that this does "
        "not; the reconciliation paragraph names only photosensitizer_pk")
    assert mine - set(theirs) == {"acsl4", "ifngamma"}, (
        f"this now counts {sorted(mine - set(theirs))} that scope_audit does "
        "not; the reconciliation paragraph attributes both to the single term "
        "`acsl4`, which that audit's regex lacks")
    # The stated CAUSE, not just the stated effect: their list has the
    # physical-ROS terms and no acsl4.
    for physical in ("photodynamic", "sonodynamic", "photosensitiz"):
        assert physical in sa.FERRO.pattern, physical
    assert "acsl4" not in sa.FERRO.pattern, (
        "scope_audit's term list now has acsl4, so the reconciliation "
        "paragraph is stale -- and the two counts should have converged")
    assert "a different number, and the difference is the point" in md


def test_the_narrated_defect_counts_are_recomputed_not_quoted(d, md):
    """The paragraphs explaining the detector's own defects had HAND-TYPED
    figures from the version that had them, and both had already gone stale:
    "eight modules" was measured before `lib.rs` was excluded and test blocks
    stripped, and "16 of 33" used a denominator counting the crate root. My
    own dominant defect class, inside the paragraph explaining a defect."""
    code, _ = MC._module_text()
    pheno = sorted(n for n, t in code.items()
                   if MC._matches(t, ("glycolytic", "oxphos")))
    tcell = sorted(n for n, t in code.items() if "t cell" in t)
    assert d["phenotype_term_credits"] == pheno
    assert d["unbounded_t_cell_credits"] == tcell
    assert f"would credit {len(pheno)} of its {d['module_count']} modules" in md
    assert f"credit {len(tcell)} of {d['module_count']} modules" in md
    # Non-vacuous: both must be a real, non-empty subset, or the sentences
    # they support say nothing.
    assert 0 < len(pheno) < d["module_count"], pheno
    assert 0 < len(tcell) < d["module_count"], tcell
    # And the terms must still be the ones the prose names.
    assert "glycolytic" not in str(MC.ENGINE_TERMS["metabolic-targeting"])
    assert "t cell" in MC.ENGINE_TERMS["immunotherapy"]


def test_the_test_block_stripper_handles_the_forms_it_used_to_miss():
    """Four legal-Rust failures, none live in this crate today, all planted --
    the same treatment `strip_rust_comments` got and this one did not."""
    # A gated attribute that is not the bare spelling.
    for attr in ("#[cfg(all(test, feature = \"x\"))]", "#[cfg( test )]",
                 "#[cfg(any(test, doctest))]"):
        src = f"pub fn a() {{1}}\n{attr}\nmod t {{ fn z() {{ let s = \"crispr\"; }} }}\npub fn b() {{2}}"
        out = MC.strip_test_blocks(src)
        assert "crispr" not in out, f"{attr} left the test block behind: {out!r}"
        assert "pub fn a()" in out and "pub fn b()" in out, out
    # A non-`mod` item ends at its semicolon, not at some later brace.
    src = "#[cfg(test)]\nuse std::collections::HashMap;\npub fn keep() {9}"
    out = MC.strip_test_blocks(src)
    assert "HashMap" not in out and "pub fn keep()" in out, out
    src = "#[cfg(test)]\nconst K: usize = 3;\npub fn keep() {9}"
    out = MC.strip_test_blocks(src)
    assert "const K" not in out and "pub fn keep()" in out, out
    # Braces and quotes inside literals must not move the depth.
    src = "#[cfg(test)]\nmod t { fn z() { let c = '{'; let q = '\"'; } }\npub fn keep() {9}"
    out = MC.strip_test_blocks(src)
    assert "pub fn keep()" in out, out
    src = '#[cfg(test)]\nmod t { fn z() { let j = r#"{"a": 1}"#; } }\npub fn keep() {9}'
    out = MC.strip_test_blocks(src)
    assert "pub fn keep()" in out, out
    # Deeper nesting: `all(any(test, doctest), ...)` is two levels, and a
    # one-level pattern left the whole module behind as engine code.
    src = ('pub fn a() {1}\n#[cfg(all(any(test, doctest), feature = "x"))]\n'
           'mod t { fn z() { let s = "crispr"; } }\npub fn b() {2}')
    out = MC.strip_test_blocks(src)
    assert "crispr" not in out, f"two-level cfg nesting left the block: {out!r}"
    assert "pub fn a()" in out and "pub fn b()" in out, out

    # A cfg that does NOT gate on test must survive untouched -- and the
    # control has to vary the NEGATION, not just the key. `\btest\b` matches
    # inside `not(test)`, whose item is PRODUCTION code, so stripping it is a
    # false negative in the direction that argues for rebuilding what exists.
    for surviving in ('#[cfg(feature = "tls")]',
                      "#[cfg(not(test))]",
                      "#[cfg(all(not(test), unix))]"):
        src = (f'{surviving}\nmod t {{ fn crispr_thing() {{}} }}\n'
               'pub fn keep() {9}')
        out = MC.strip_test_blocks(src)
        assert "crispr_thing" in out, f"{surviving} was wrongly stripped: {out!r}"
    # Both directions of the predicate itself, so a rule that returns a
    # constant fails here.
    for attr, gated in (("#[cfg(test)]", True),
                        ("#[cfg(all(test, feature = \"x\"))]", True),
                        ("#[cfg(all(any(test, doctest), unix))]", True),
                        ("#[cfg(not(test))]", False),
                        ("#[cfg(all(not(test), unix))]", False),
                        ("#[cfg(feature = \"tls\")]", False)):
        assert MC._cfg_test_matches(attr) is gated, attr
    # Line-preserving, so line numbers stay true to the source.
    src = "a\n#[cfg(test)]\nmod t {\n  fn z() {}\n}\nb"
    assert MC.strip_test_blocks(src).count("\n") == src.count("\n")


def test_the_gating_sentence_is_derived_from_the_crate(d, md):
    """The paragraph used to say both immune models are "gated on ferroptotic
    death by construction". The immunotherapy arm (#728) made that false the
    day it landed, and nothing but this file's regeneration noticed.

    So the sentence is read off `params.rs` now, and this pins both halves:
    the field's existence and its DEFAULT. A nonzero default would move every
    committed immune number in the repository.
    """
    b = d["baseline_antigenicity"]
    params = MC.strip_test_blocks(
        MC.strip_rust_comments((CORE / "params.rs").read_text()))
    assert b["exists"] == (b["field"] in params)
    if not b["exists"]:
        assert "no ferroptosis-independent antigen source" in md
        return
    assert b["default"] == 0.0, (
        f"{b['field']} defaults to {b['default']}, not 0 — every committed "
        "immune number has moved and the report's 'OFF by default' sentence "
        "is false")
    assert "OFF by default" in md
    # Which modules consume it must be MEASURED, not listed: the sentence
    # names one as still gated and one as not.
    live = sorted(
        n for n in ("immune", "immune_spatial")
        if b["field"] in MC.strip_test_blocks(
            MC.strip_rust_comments((CORE / f"{n}.rs").read_text())))
    assert b["consumed_by"] == live
    assert live, "nothing consumes the field, so the arm is inert"
    for name in live:
        assert f"`{name}`" in md
    # And the tier must NOT have moved on the strength of a knob: applying a
    # modality needs a Treatment variant, which this does not add.
    imm = next(r for r in d["rows"] if r["mechanism"] == "immunotherapy")
    assert imm["engine_tier"] == "modifier", (
        "immunotherapy is no longer a MODIFIER; if a `Treatment` variant "
        "landed, the report's explanation of why it was one needs rewriting")
    assert "still reads MODIFIER" in md


def test_the_rows_match_the_census_profile(d):
    """The counts are copied, not recomputed, so they must still be copies."""
    prof = json.loads((REPO / "analysis/census-mechanism-profile.json").read_text())
    src = {r["mechanism"]: r for r in prof["rows"]}
    assert d["census"] == prof["census"]
    assert {r["mechanism"] for r in d["rows"]} == set(src)
    for r in d["rows"]:
        for k in ("census", "trials", "trial_share"):
            assert r[k] == src[r["mechanism"]][k], f"{r['mechanism']}.{k} drifted"
