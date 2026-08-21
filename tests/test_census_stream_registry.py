"""Every consumer of a census stream must declare which streams it reads.

WHY THIS EXISTS
---------------
The cancer census is stored as parallel streams of gzipped JSONL under
corpus/atlas/, and consumers enumerate them BY NAME:

    records/            the MeSH-indexed census, ~4.4M articles
    records_c04only/    the same, restricted to MeSH tree C04 (a sensitivity
                        control, a SUBSET of records/ rather than more articles)
    records_unindexed/  ~783k text-recovered articles MeSH has not indexed,
                        disjoint from records/
    records_updates/    the daily update stream: new articles plus revisions

A stream gets added, an existing consumer is never updated, and it keeps
answering for the world as it was. This happened FOUR times in one day, three
shipped and one introduced while fixing another, and every instance failed
silently -- no error, no empty result, just a quietly wrong number:

  1. `census_pmids` read records/ alone, so 20,345 articles already held in
     records_unindexed/ were reported as NEW -- a 24% overstatement.
  2. `load_pmcid_map` read two of three, so a full-text recovery matched zero
     articles while reporting success.
  3. `load_cancer_pmids` read two of three, so the newest 65,966 articles
     carried no typed relations, which reads as "the extractor has nothing for
     them" rather than "they were never offered to the filter".
  4. Re-running one source of a four-source ingest left the other three
     filtering against the pre-update census, so those articles had relations
     but zero entity annotations.

WHAT THIS GUARD DOES, AND WHAT IT DELIBERATELY DOES NOT
--------------------------------------------------------
It does NOT try to decide which streams a consumer should read. That judgement
is real work and often the answer is records/ alone -- the frozen census is
what every committed manuscript figure was computed on, and widening it there
would silently change published numbers. A guard that demanded all streams
everywhere would be wrong more often than the bug it replaced.

What it does is remove the SILENCE. Every module naming a stream must appear in
the registry below with the streams it reads and the reason. A new consumer,
or an existing one that starts naming a different set, fails here until someone
writes down what it is for. The decision still belongs to a person; this only
guarantees the decision gets made.
"""

import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPTS = REPO_ROOT / "scripts"

STREAMS = ("records", "records_c04only", "records_unindexed", "records_updates")

# Matches a stream name used as a path component: in a quoted literal, a glob,
# or a path expression. Deliberately broad -- a false positive costs one
# registry line, a false negative costs a silent wrong number.
_PATTERN = re.compile(
    r"""["'](records(?:_c04only|_unindexed|_updates)?)["']"""
    r"""|(?:/|\\)(records(?:_c04only|_unindexed|_updates)?)(?:/|["'])"""
)

# module -> (streams it reads, why that set is the right one)
REGISTRY = {
    "census_evidence_design.py": (
        {"records"},
        "Study design is read from NLM publication types and MeSH check tags, "
        "which only the INDEXED stream carries. records_unindexed has neither, "
        "so there is no label of this kind to compute for it -- and the report "
        "states that denominator rather than quietly using the smaller one."),
    "census_oa_bias.py": (
        {"records"},
        "Compares mechanism ranking with and without a PMC identifier on "
        "identical MeSH descriptors. Both arms must come from the same stream "
        "or the comparison confounds availability with indexing; "
        "records_unindexed carries no descriptors and so cannot supply either "
        "arm."),
    "census_diagnostic_chains.py": (
        {"records"},
        "Chain membership is matched over title, MeSH and abstract. The "
        "MeSH channel is the one the published corpus figure used and is "
        "carried only by the indexed stream, so folding in records_unindexed "
        "would change the instrument mid-measurement -- exactly the thing this "
        "analysis exists to hold fixed between its two arms."),
    "census_mechanism_sites.py": (
        {"records"},
        "Both the site assignment and the mechanism labels are MeSH "
        "descriptors, and records_unindexed has no descriptors at all, so it "
        "can supply neither axis. Every record it holds would land in the "
        "unassigned bucket and depress the assignability rate without "
        "carrying any information about site."),
    "census_thesis_direction.py": (
        {"records"},
        "Both selectors are MeSH descriptors -- `Ferroptosis` and `Drug "
        "Resistance, Neoplasm` -- so records_unindexed cannot supply the leg "
        "at all. Its text is read only to CLASSIFY articles the descriptors "
        "already selected, which is a different role from selection and does "
        "not extend to a stream where selection is impossible."),
    "census_modality_comparison.py": (
        {"records"},
        "Trial status comes from NLM publication types and one arm from MeSH "
        "descriptors, neither of which records_unindexed carries. Its text arm "
        "COULD read that stream, and deliberately does not: the two arms are "
        "compared against each other as a descriptor-validity test, so they "
        "must read the same population or a divergence would measure the "
        "streams rather than the descriptor."),
    "census_mechanism_cancer_matrix.py": (
        {"records"},
        "Both axes are MeSH descriptors -- mechanism and anatomical site -- so "
        "records_unindexed can supply neither. Worse than useless here: its "
        "records would enter neither the universe nor the marginals, but "
        "counting them in the census total would make the universe look like a "
        "smaller share of the literature than it is."),
    "census_translation_lag.py": (
        {"records"},
        "Dates a mechanism's literature and its first trial. Trial status "
        "comes from NLM publication types, which only the indexed stream "
        "carries, so the trial end of every lag is undefined for "
        "records_unindexed -- and its text end would then start clocks the "
        "trial end could never stop, turning every mechanism censored."),
    "census_external_check.py": (
        {"records"},
        "Compares the census against PubMed's own index, and the comparison "
        "only holds if both sides admit records the same way. `neoplasms[mh]` "
        "returns C04-indexed records, so the census side must be the C04 CORE "
        "of the indexed stream -- records_unindexed carries no MeSH at all and "
        "every record in it would be a census-side record PubMed's query "
        "cannot return, turning the whole check into a measurement of the "
        "streams' definitions."),
    "census_fulltext_ceiling.py": (
        {"records"},
        "Measures how much of the INDEXED stream's design-label gap open-access "
        "full text could reach, so both the numerator and the denominator must "
        "come from that stream. records_unindexed carries no publication types, "
        "so every record in it is undetermined BY CONSTRUCTION -- folding it in "
        "would inflate the gap with records that were never a labelling failure "
        "and depress the ceiling with them."),
    "corpus_dependency_audit.py": (
        {"records"},
        "Reads the census only to learn its record SCHEMA -- which fields a "
        "consumer could be pointed at -- and one shard is as informative as "
        "all of them. records_unindexed would add nothing, since a "
        "text-recovered record carries a SUBSET of the same fields, and a "
        "consumer needing a field it lacks is already accounted for by the "
        "indexed stream's schema."),
    "census_mechanism_profile.py": (
        {"records"},
        "Joins mechanism descriptors, C04 site descriptors and NLM publication "
        "types, all three of which only the indexed stream carries. Adding "
        "records_unindexed would contribute records to no mechanism, no site "
        "and no trial column while inflating nothing but the pass count."),
    "census_mechanism_growth.py": (
        {"records"},
        "Mechanism labels are MeSH descriptors, so the numerator can only come "
        "from the indexed stream, and the field denominator must come from the "
        "same stream or the growth ratio compares a MeSH-labelled numerator "
        "against a denominator that includes records MeSH never reached -- "
        "which would bias the ratio DOWN in exactly the recent years the claim "
        "is about."),
    "atlas_ingest_sensitivity.py": (
        {"records"},
        "The qualifier measurement itself re-parses the raw baseline XML, which "
        "is not a census stream at all -- the committed records carry only "
        "DescriptorName, which is the very gap this analysis measures. It reads "
        "records/ for ONE thing: checking its own sample's era coverage against "
        "the census, since the 8 sampled shards have committed counterparts "
        "with a year. records/ alone is right there because the question is "
        "what the MeSH-indexed census contains; the text-recovered stream has "
        "no MeSH and so nothing to compare against."),
    # --- reads every stream, because it is asking what the project HOLDS ---
    "atlas_baseline.py": (
        {"records", "records_unindexed", "records_updates"},
        "census_pmids answers 'do we already have this article', so it must see "
        "both census streams plus anything an interrupted run already wrote. "
        "Reading records/ alone called 20,345 held articles new."),
    "atlas_fulltext.py": (
        {"records", "records_unindexed", "records_updates"},
        "load_pmcid_map decides which PMC bulk packages can match, so an article "
        "in any stream is fair game. Missing records_updates/ left the PMC13 "
        "recovery cliff in place while the run reported success."),
    "atlas_relations.py": (
        {"records", "records_unindexed", "records_updates"},
        "the PubTator bulk files are filtered against the census PMID set; any "
        "held article should get its relations and entity annotations."),
    "atlas_recent_window.py": (
        {"records", "records_unindexed", "records_updates"},
        "its whole subject is the difference between the streams."),

    # --- reads the frozen census alone, ON PURPOSE ---
    "atlas_coverage.py": (
        {"records"},
        "measures the frozen corpus against the frozen census; both sides must "
        "be the snapshot the manuscript quotes or the ratio is not comparable."),
    "atlas_landscape.py": (
        {"records"},
        "recomputes the manuscript's central corpus claim, so the denominator "
        "has to be the census the manuscript reports."),
    "manuscript_vs_census.py": (
        {"records"},
        "re-tests published manuscript claims; widening the denominator would "
        "test a different claim than the one in print."),
    "atlas_thesis_position.py": (
        {"records"},
        "positions the thesis against the census figure the repo publishes."),
    "atlas_prediction_position.py": (
        {"records"},
        "same frozen denominator as the other position analyses, so the legs "
        "are comparable with each other."),
    "atlas_model_gaps.py": (
        {"records"},
        "ranks field attention via PubTator gene annotations, which cover the "
        "frozen census; adding articles the annotation layer has not processed "
        "would depress every rank without evidence."),
    "atlas_evidence_check.py": (
        {"records"},
        "validates the evidence tagger against NLM publication types on the "
        "frozen records the tagger was measured on."),
    "atlas_retraction_exposure.py": (
        {"records", "records_updates"},
        "joins against the relation graph, which was filtered over ALL census "
        "streams, so records/ alone left 393,972 of its PMIDs unreachable. "
        "records_unindexed/ is excluded on a MEASURED ground rather than an "
        "assumed one: none of its 783,271 records carries a pub_types field, "
        "so it cannot carry a retraction flag."),
    "atlas_site_coverage.py": (
        {"records"},
        "gates a burden-weighted analysis on the same frozen census every other "
        "coverage figure divides by, so the assignability rate it reports is "
        "comparable with them."),
    "atlas_descriptor_recall.py": (
        {"records"},
        "measures how completely two MeSH descriptors recall their concepts, "
        "so it must read the stream those descriptors are assigned in. "
        "records_unindexed/ carries no MeSH at all, so including it would put "
        "articles in the denominator that CANNOT contribute to the numerator "
        "and depress both recalls by construction -- and unequally, since the "
        "un-indexed share differs by topic. The comparison is between two "
        "descriptors on the same articles, which is what makes it symmetric."),
    "atlas_thesis_rank.py": (
        {"records"},
        "ranks the ferroptosis-modality intersection over the same frozen "
        "census every other position analysis uses, so the legs stay comparable "
        "with the thesis-position figures they are meant to sit beside."),
    "atlas_modality_ratio.py": (
        {"records"},
        "recomputes the manuscript's pharmacological-versus-physical ratio, so "
        "it must divide by the same frozen census the manuscript's own figure "
        "was computed on or the comparison is not the one being tested."),
    "atlas_taxonomy_reach.py": (
        {"records"},
        "measures the mechanism taxonomy's field of view against the census the "
        "rest of the repo quotes, so it must use the same 4,403,994 denominator "
        "every capture figure is compared to."),
    "atlas_unindexed.py": (
        {"records_unindexed"},
        "it WRITES that stream; it is the recovery layer itself."),

    # --- dating tables, which must cover every stream the graph contains ---
    #
    # BOTH ENTRIES HERE WERE WRONG WHEN FIRST WRITTEN, and wrong in the way
    # this file exists to prevent. They read {records, records_c04only} and the
    # reason given was that the module "reports emergence both ways so a result
    # cannot rest on the nine adjacent descriptors alone". That was a
    # rationalisation of a set nobody chose: c04only is a strict SUBSET of
    # records, both directories merge into ONE dict, and one table is emitted,
    # so there was no second arm at all. A plausible reason for an accidental
    # set is worse than no reason, because it stops the next reader looking.
    "atlas_emergence.py": (
        {"records", "records_c04only", "records_unindexed", "records_updates"},
        "a PMID->year table joined against the relation graph, not a "
        "denominator. The graph covers every stream, so a narrower year map "
        "silently DROPS pairs it cannot date -- and the dropped ones are "
        "overwhelmingly recent, which is exactly what this layer measures."),
    "atlas_discovery_eval.py": (
        {"records", "records_c04only", "records_unindexed", "records_updates"},
        "the same dating table, and the one with the widest blast radius: four "
        "other analyses import it, so a narrow year map propagated everywhere."),
}


def _streams_named(text: str) -> set:
    found = set()
    for m in _PATTERN.finditer(text):
        name = m.group(1) or m.group(2)
        if name in STREAMS:
            found.add(name)
    return found


def _consumers() -> dict:
    out = {}
    for p in sorted(SCRIPTS.glob("*.py")):
        named = _streams_named(p.read_text(encoding="utf-8"))
        if named:
            out[p.name] = named
    return out


def test_the_detector_finds_planted_samples():
    """A scan returning nothing because it is broken looks like a clean repo.

    This repository has that exact lesson written down, so the detector is
    exercised against samples rather than trusted.
    """
    assert _streams_named('root / "records_updates"') == {"records_updates"}
    assert _streams_named('for d in ("records", "records_unindexed"):') == {
        "records", "records_unindexed"}
    assert _streams_named('glob(str(root / "records" / "*.jsonl.gz"))') == {"records"}
    assert _streams_named("path = base / 'records_c04only' / name") == {"records_c04only"}
    assert _streams_named("corpus/atlas/records_unindexed/x.gz") == {"records_unindexed"}
    # and it must not fire on unrelated words
    assert _streams_named('"medical_records"') == set()
    assert _streams_named("# the records were parsed") == set()


def test_every_census_consumer_is_registered():
    """The whole point: a new consumer cannot be silent about what it reads."""
    found = _consumers()
    missing = sorted(set(found) - set(REGISTRY))
    assert not missing, (
        "these scripts name a census stream but are not in REGISTRY:\n  "
        + "\n  ".join(f"{m} (names {', '.join(sorted(found[m]))})" for m in missing)
        + "\n\nAdd each to tests/test_census_stream_registry.py with the streams "
          "it reads AND why that set is right. Reading records/ alone is often "
          "correct -- the frozen census is the comparable denominator -- but it "
          "has to be a decision someone made, not one nobody noticed.")


def test_the_registry_has_no_entries_for_scripts_that_stopped_reading():
    """A stale registry entry is a claim about code that no longer exists."""
    found = _consumers()
    ghosts = sorted(set(REGISTRY) - set(found))
    assert not ghosts, (
        f"REGISTRY describes {ghosts}, which no longer name any census stream. "
        "Remove the entries, or the registry is documenting a world that moved.")


def test_each_registered_set_matches_what_the_module_names():
    """The registry must describe the code, not the intention.

    Detected names are a superset in principle -- a module may MENTION a stream
    in a docstring without reading it -- so this checks the registry never
    claims a stream the module does not name at all, which is the direction
    that hides a missing read.
    """
    found = _consumers()
    wrong = []
    for mod, (declared, _why) in REGISTRY.items():
        named = found.get(mod, set())
        phantom = declared - named
        if phantom:
            wrong.append(f"{mod}: registry claims {sorted(phantom)} which the "
                         f"module never names (it names {sorted(named)})")
    assert not wrong, "\n  ".join([""] + wrong)


def test_every_registry_entry_gives_a_reason():
    """A registry of bare sets would record the bug rather than prevent it."""
    thin = [m for m, (_s, why) in REGISTRY.items() if len(why.split()) < 8]
    assert not thin, (
        f"these registry entries give no real reason: {thin}. The reason is the "
        "part that survives; the set alone is what four separate consumers "
        "already got wrong.")


def test_the_four_known_instances_stay_fixed():
    """Pin the specific reads whose absence caused a measured wrong number."""
    for mod in ("atlas_baseline.py", "atlas_fulltext.py", "atlas_relations.py"):
        declared = REGISTRY[mod][0]
        assert "records_updates" in declared, (
            f"{mod} no longer reads records_updates/; that omission is what "
            "made 65,966 articles invisible to three different layers")
    src = (SCRIPTS / "atlas_fulltext.py").read_text()
    assert '("records", "records_unindexed", "records_updates")' in src, (
        "load_pmcid_map's stream tuple changed; the PMC13 recovery depends on it")
