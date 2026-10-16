"""Is "MeSH has no descriptor for it" a statement about MeSH, or about our map?

THE CLAIM AND WHERE IT COMES FROM
---------------------------------
Several mechanisms are reported throughout this project as *not measurable at
census scale*, and the reason given -- in generated pages and in the manuscript
-- is that **MeSH has no descriptor for them**.

Not one of those sites measures that. Every consumer derives it the same way::

    # census_unnamed_modalities._arms
    mesh = {k.lower() for k, v in mp.items() if v["descriptors"]}
    invisible = {k: v for k, v in text.items() if k not in mesh}
    # census_translation_lag / census_modality_comparison / census_external_check
    if not r["mesh_measurable"]: ...

`mesh_measurable` is a property of `analysis/mesh-mechanism-map.yaml`, a
handwritten file whose own header records the corpus its counts were verified
against -- the frozen corpus, orders of magnitude smaller than the census. So
what is measured is *this file names no usable descriptor for X*. What is
carried is *MeSH has no descriptor for X, so the census cannot see it*. One
word of scope -- the map, not MeSH -- and a sample far smaller than the
population the claim is about.

This repository has already paid for that class once on this exact axis.
`atlas_thesis_rank.py` records it: "a hand-written no-descriptor list was
refuted by the census it sits beside" -- `Plasma Gases` exists and carries
census records, against a note saying no descriptor did. The fix there was to
name a candidate descriptor per modality and check it, which narrows the check
to the candidate the author thought of and is the same shape one level down.
`census_unnamed_modalities.py` came closest and stopped one level short: its
docstring says "not measurable BY DESCRIPTOR, which is a statement about MeSH"
-- correcting one level and leaving the second, since it is a statement about
the map -- while its rendered first line still read "because MeSH has no
descriptor for them. That is the right call."

WHAT THIS PAGE DOES INSTEAD
---------------------------
It makes the check as wide as the claim, over the census's OWN descriptor
vocabulary -- every distinct descriptor NLM applied to any census record, with
its record count -- rather than over a hand list.

Two arms, because a single rule that licenses a negative claim is the thing
this page exists to criticise:

  ARM A, the production rule. A descriptor is a candidate when the mechanism's
  own `MECHANISM_KEYWORDS` accept its NAME under `tag_articles`'
  `text_matches_keyword`. Keyed on the live matcher so it cannot drift from the
  tagger, and it is the same rule `mechanism_recall.py` uses to decide whether
  a descriptor is LEAKY -- so an arm-A hit is by construction a descriptor that
  names the mechanism in the mechanism's own words.

  ARM B, a looser control. Token-subset matching with light normalisation
  (plural stripping, punctuation folding, and the generic words a therapy
  descriptor shares with every other one removed from the KEYWORD side only, so
  they can never manufacture a match the descriptor's name does not support).
  Wider than arm A, so an absence verdict is an absence under both -- and arm
  A's own blind spot is measured rather than assumed away, since arm A is a
  SUBSTRING rule and can match inside a token where arm B cannot. Arm B is a
  recall instrument with poor precision by design (`cap therapy` reaches `RNA
  Caps`) and its loose-only hits are published with their counts so a reader
  judges them rather than taking a verdict.

NEITHER ARM IS SEMANTIC, and that is the honest limit. Both are rules over
descriptor NAMES, so a descriptor that means the mechanism under different
words is invisible to both, exactly as it is invisible to a hand list. What
changes is the scope: from the candidates an author recalled to every one of
the census's descriptors that the mechanism's own production vocabulary
accepts. Wider is not complete.

THE VERDICT STAYS HANDWRITTEN. Whether a broad descriptor CONTAINS a mechanism
(so a count of it over-counts) or merely RESEMBLES it (so a count of it is a
different thing) is a judgement about MeSH semantics that no name rule can
make, and this project's own rule is that judgements live in the handwritten
map. So the map declares, per mechanism, a `verdict` plus the descriptors it
considers `covering` and `adjacent`, and THIS page checks the checkable half:

  * every descriptor the map names must EXIST in the census vocabulary -- the
    defect `atlas_thesis_rank.py` found in its own list, where several stems
    matched nothing at all;
  * `absent` is refused where the map names a covering descriptor, or where
    either arm finds a name-exact one;
  * `over-broad` is refused where either arm finds a name-exact descriptor, so
    breadth can never be reported where absence is not even the question;
  * `name-exact-leaky` is refused unless arm A finds one, and the map must name
    every leaf arm A finds, so a descriptor discovered here cannot sit
    unrecorded in the file the consumers read.

VOLUME-MEASURABILITY AND RECALL-MEASURABILITY ARE DIFFERENT CLAIMS and
collapsing them is half of why this went unnoticed. A descriptor whose name
trips the mechanism's keywords is LEAKY by `mechanism_recall.py`'s own
classifier, so it CANNOT serve as an independent reference label and that
page's headline is right to exclude it. It can perfectly well carry a census
VOLUME count, which is what `census_mechanism_profile.py` reports for every
measurable mechanism and what the sentences corrected here deny is possible.
The recall verdict survives; the absence claim does not.

AND A NAME-EXACT DESCRIPTOR IS NOT AUTOMATICALLY A CLEAN COUNT. A descriptor
can be older than the therapy it now names -- the trap
`census-unnamed-modalities.md` documents for the keyword arm, where a copper
ionophore term reaches pre-cuproptosis chemistry -- so the share of each pool
predating the modality is measured and a suspect pool is reported as an upper
bound rather than a measurement.

AND THE SECOND ARM BECOMES AVAILABLE for the mechanisms with a name-exact
leaf, which changes a sibling page's stated justification rather than merely
adding to it. `census-unnamed-modalities.md` rests on the descriptor-vs-text
agreement check being unavailable for these mechanisms "BY CONSTRUCTION --
there is no second arm, because the absence of a second arm is what puts them
on this list". Where a leaf exists that is no longer true, so this page
measures the agreement both ways and that page is corrected at its generator.

Offline contract: the census bulk is never committed. The committed
derivatives are the descriptor vocabulary (every descriptor, its record count
and its major-topic count), the per-candidate detail (year distribution,
publication types, and the two-way text-arm agreement) and the site registry.
Both arms are RE-DERIVED from the committed vocabulary on every render through
the live matcher, so no verdict in this page is echoed back out of its own
artifact, and CI runs no census scan.
"""
import argparse
import csv
import gzip
import io
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from census_input import iter_census_records  # noqa: E402
from config import MECHANISM_KEYWORDS  # noqa: E402
from tag_articles import text_matches_keyword  # noqa: E402

REPO = Path(__file__).resolve().parent.parent
MECH_MAP = REPO / "analysis" / "mesh-mechanism-map.yaml"
VOCAB_TSV_GZ = REPO / "analysis" / "census-descriptor-vocabulary.tsv.gz"
DETAIL_JSON = REPO / "analysis" / "census-descriptor-candidates.json"
OUT_JSON = REPO / "analysis" / "descriptor-absence-audit.json"
OUT_MD = REPO / "analysis" / "descriptor-absence-audit.md"

# The three verdicts, each licensing a DIFFERENT claim. They are not severity
# tiers and must not be collapsed: `CLAUDE.md` records that reporting breadth
# as absence is the error this page exists for.
VERDICTS = ("absent", "over-broad", "name-exact-leaky")

# Generic words a therapy descriptor shares with every other therapy
# descriptor. Dropped from arm B's token sets so a keyword's DISTINCTIVE tokens
# decide the match -- otherwise every keyword ending in "therapy" reaches every
# descriptor ending in "Therapy". They are dropped from the KEYWORD side only,
# so they can never manufacture a match the descriptor's own name does not
# support.
ARM_B_GENERIC = frozenset({
    "the", "of", "and", "a", "an", "to", "for", "in", "with",
    "therapy", "treatment", "cancer", "tumor", "tumour", "therapeutic",
})

# A pool this small cannot support a census volume count whatever its label
# says, so a name-exact leaf below it is reported as found-but-unusable rather
# than as a correction. Declared, not measured: it is the same order as
# `mechanism_recall.N_MIN`, which is imported rather than restated so the two
# floors cannot drift.
from mechanism_recall import N_MIN as POOL_FLOOR  # noqa: E402

# A DESCRIPTOR CAN BE OLDER THAN THE THERAPY IT NOW NAMES, which is the trap
# `census-unnamed-modalities.md` documents for the keyword arm -- `copper
# ionophore` applied to disulfiram in 1982, four decades before cuproptosis was
# named -- and it applies to a NAME-EXACT descriptor just as hard. Restated
# here rather than imported because that module imports this one, so the import
# would cycle; a guard pins the two equal.
VINTAGE_YEAR = 1990
VINTAGE_SUSPECT_SHARE = 0.10

# A covering descriptor has to be materially wider than the mechanism for
# "over-broad" to mean anything. Expressed as a multiple of the pool floor
# rather than an absolute count, so it moves with the floor it is relative to.
OVER_BROAD_MULTIPLE = 3


_MAP_CORPUS = re.compile(r"mesh_terms \((\d+) records")


def map_corpus_size() -> int:
    """The record count the map's own header says its notes were verified on.

    PARSED rather than typed, because the whole finding is that a figure
    measured on that corpus was carried as a reason about the census. A number
    written here could drift from the file it describes, which is the same
    defect one level up.
    """
    text = MECH_MAP.read_text(encoding="utf-8")
    hit = _MAP_CORPUS.search(text)
    if not hit:
        raise SystemExit(
            f"{MECH_MAP.name} no longer states the corpus its notes were "
            "verified against; this page quotes it")
    return int(hit.group(1))


def _yaml():
    import yaml

    return yaml.safe_load(MECH_MAP.read_text(encoding="utf-8"))


def unusable(mapping: dict) -> dict:
    """Every mechanism the map gives no usable descriptor, from BOTH blocks.

    A mechanism with no usable leaf can sit in `unmeasurable:` or in
    `mechanisms:` with an empty `descriptors` list, and `frequency-therapy` has
    always been the second because `census_mechanism_profile.py` and
    `census_external_check.py` key on it being there. Deriving the set from
    both means a verdict does not depend on which block a mechanism happens to
    occupy -- and listing a mechanism in both is refused by
    `tests/test_mechanism_recall.py`, so the two are disjoint by construction.
    `reason` falls back to `note`, which is what the `mechanisms:` block calls
    the same field.
    """
    out = {}
    for mech, entry in (mapping.get("unmeasurable") or {}).items():
        out[mech] = entry
    for mech, entry in (mapping.get("mechanisms") or {}).items():
        if entry.get("descriptors"):
            continue
        if "verdict" not in entry:
            raise SystemExit(
                f"{mech} has no usable descriptor and no verdict field; add "
                "verdict/covering/adjacent to its entry in "
                f"{MECH_MAP.name}")
        out[mech] = dict(entry, reason=entry.get("reason") or entry.get("note", ""))
    return out


def _normalise(text: str) -> list:
    """Arm B's token rule: fold punctuation, strip a plural s."""
    out = []
    for token in re.split(r"[^a-z0-9]+", text.lower()):
        if not token:
            continue
        if len(token) > 3 and token.endswith("s") and not token.endswith("ss"):
            token = token[:-1]
        out.append(token)
    return out


def arm_a(name: str, keywords) -> list:
    """Descriptor names the mechanism in the mechanism's own words.

    The LIVE production matcher, so this arm and `mechanism_recall.py`'s leak
    classifier can never disagree about what "the descriptor's own string would
    trip the tagger" means.
    """
    lowered = name.lower()
    return [k for k in keywords if text_matches_keyword(lowered, k)]


def arm_b(name: str, keywords) -> list:
    """Looser control: every distinctive keyword token present in the name."""
    tokens = set(_normalise(name))
    hits = []
    for keyword in keywords:
        distinctive = [t for t in _normalise(keyword) if t not in ARM_B_GENERIC]
        if distinctive and set(distinctive) <= tokens:
            hits.append(keyword)
    return hits


def _mechanism_patterns(mechanisms):
    """Word-boundary text-arm patterns, the rule census_unnamed_modalities uses."""
    out = {}
    for mech in mechanisms:
        terms = MECHANISM_KEYWORDS.get(mech) or []
        if terms:
            out[mech] = re.compile(
                "|".join(rf"\b{re.escape(t.lower())}\b" for t in terms))
    return out


def measure_vocabulary(stride: int = 1) -> dict:
    """Pass 1: every descriptor the census carries, with its record counts."""
    names, records, major = {}, Counter(), Counter()
    scanned = 0
    for record in iter_census_records(stride=stride):
        scanned += 1
        uis = record.get("mesh_ui") or []
        labels = record.get("mesh") or []
        for ui, label in zip(uis, labels):
            names[ui] = label
            records[ui] += 1
        for ui in (record.get("mesh_major") or []):
            major[ui] += 1
    return {"scanned": scanned, "names": names,
            "records": records, "major": major}


def write_vocabulary(vocab: dict) -> None:
    rows = sorted(vocab["names"],
                  key=lambda u: (-vocab["records"][u], vocab["names"][u]))
    # mtime=0 so an unchanged census regenerates a BYTE-IDENTICAL file, the
    # property `atlas_variant_drug_map.py` had to learn the hard way.
    with gzip.GzipFile(VOCAB_TSV_GZ, "wb", mtime=0) as raw:
        fh = io.TextIOWrapper(raw, encoding="utf-8", newline="")
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["ui", "name", "records", "major"])
        writer.writerow(["#scanned", str(vocab["scanned"]), "", ""])
        for ui in rows:
            writer.writerow([ui, vocab["names"][ui],
                             vocab["records"][ui], vocab["major"].get(ui, 0)])
        fh.flush()
        fh.detach()
    print(f"wrote {VOCAB_TSV_GZ} ({len(rows):,} descriptors, "
          f"{vocab['scanned']:,} records)")


def read_vocabulary() -> dict:
    if not VOCAB_TSV_GZ.exists():
        raise SystemExit(
            f"{VOCAB_TSV_GZ} is missing; run with --measure (needs the census)")
    scanned = None
    rows = []
    with gzip.open(VOCAB_TSV_GZ, "rt", encoding="utf-8") as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if not row or row[0] == "ui":
                continue
            if row[0] == "#scanned":
                scanned = int(row[1])
                continue
            rows.append({"ui": row[0], "name": row[1],
                         "records": int(row[2]), "major": int(row[3])})
    if scanned is None:
        raise SystemExit(f"{VOCAB_TSV_GZ} carries no scanned-record total")
    return {"scanned": scanned, "rows": rows}


def candidate_uis(vocab_rows, unmeasurable) -> dict:
    """Which descriptors each arm reaches, and which the map names.

    Pure: no census, no stored derived field. Called by `measure` to choose what
    to look at in pass 2 AND by `assemble` to re-derive the arms on every
    render, so the two can never describe different candidate sets.
    """
    by_name = {}
    for row in vocab_rows:
        by_name.setdefault(row["name"].lower(), row)
    out = {}
    for mech, entry in sorted(unmeasurable.items()):
        keywords = MECHANISM_KEYWORDS.get(mech) or []
        a, b = {}, {}
        for row in vocab_rows:
            hit_a = arm_a(row["name"], keywords)
            if hit_a:
                a[row["ui"]] = hit_a
            hit_b = arm_b(row["name"], keywords)
            if hit_b:
                b[row["ui"]] = hit_b
        named = {}
        for field in ("covering", "adjacent", "leaves"):
            for name in (entry.get(field) or []):
                row = by_name.get(name.lower())
                named[name] = {"field": field,
                               "ui": row["ui"] if row else None,
                               "records": row["records"] if row else None,
                               "major": row["major"] if row else None}
        out[mech] = {"arm_a": a, "arm_b": b, "named": named,
                     "keywords": list(keywords)}
    return out


def measure_detail(uis, mechanisms, stride: int = 1) -> dict:
    """Pass 2: year, publication type and two-way text agreement per candidate.

    The text agreement is the second arm `census-unnamed-modalities.md` reports
    as unavailable BY CONSTRUCTION for these mechanisms. Both directions are
    counted because they answer different questions: how much of the
    descriptor's pool the project's own vocabulary finds, and how much of the
    project's text-arm pool NLM filed under the descriptor.
    """
    uis = set(uis)
    patterns = _mechanism_patterns(mechanisms)
    prefilter = (re.compile("|".join(p.pattern for p in patterns.values()))
                 if patterns else None)
    years = defaultdict(Counter)
    types = defaultdict(Counter)
    desc_total = Counter()
    text_total = Counter()
    both = Counter()
    scanned = 0
    for record in iter_census_records(stride=stride):
        scanned += 1
        record_uis = set(record.get("mesh_ui") or []) & uis
        blob = f"{record.get('title') or ''} {record.get('abstract') or ''}".lower()
        text_hit = set()
        if blob and prefilter is not None and prefilter.search(blob):
            text_hit = {m for m, p in patterns.items() if p.search(blob)}
        for mech in text_hit:
            text_total[mech] += 1
        if not record_uis:
            continue
        year = record.get("year")
        pub_types = record.get("pub_types") or []
        for ui in record_uis:
            desc_total[ui] += 1
            if isinstance(year, int):
                years[ui][year] += 1
            for pub_type in pub_types:
                types[ui][pub_type] += 1
        for mech in text_hit:
            for ui in record_uis:
                both[f"{mech}|{ui}"] += 1
    return {"scanned": scanned,
            "descriptor_records": dict(desc_total),
            "text_arm_records": dict(text_total),
            "both": dict(both),
            "by_year": {u: dict(sorted(c.items())) for u, c in years.items()},
            "pub_types": {u: dict(c) for u, c in types.items()}}


def measure(stride: int = 1) -> None:
    vocab = measure_vocabulary(stride)
    write_vocabulary(vocab)
    stored = read_vocabulary()
    unmeasurable = unusable(_yaml())
    cands = candidate_uis(stored["rows"], unmeasurable)
    wanted = set()
    for mech, found in cands.items():
        wanted |= set(found["arm_a"]) | set(found["arm_b"])
        wanted |= {v["ui"] for v in found["named"].values() if v["ui"]}
    detail = measure_detail(wanted, sorted(unmeasurable), stride)
    DETAIL_JSON.write_text(json.dumps(detail, indent=1, sort_keys=True) + "\n")
    print(f"wrote {DETAIL_JSON} ({len(wanted)} candidate descriptors)")


def load() -> dict:
    if not DETAIL_JSON.exists():
        raise SystemExit(
            f"{DETAIL_JSON} is missing; run with --measure (needs the census)")
    return {"vocabulary": read_vocabulary(),
            "detail": json.loads(DETAIL_JSON.read_text()),
            "map": _yaml()}


# THE REGISTRY LIVES IN DATA, not here, and the reason is a rule this project
# already enforces: a quoted claim legitimately contains its author's figures
# and a generator may not, so a registry of quoted sentences cannot sit in
# generator source without hand-writing quantities into it. Keeping the
# quotations in `analysis/descriptor-absence-sites.csv` is what makes the two
# rules compatible -- the same arrangement `intervals_under_admissibility.py`
# uses for `published-interval-claims.csv`.
#
# FOUR COLUMNS, and the third is what makes this more than a reading list.
# `phrase` must still be in the file; `retired` is the flattened wording this
# work removed and must now be ABSENT. A correction verified only by its
# replacement being present cannot tell a fix from a fix shipped beside a
# re-added original, which is the "verify the old text is gone" shape already
# on record here. `concerns` names the mechanisms the sentence is about, so the
# scoring checks the sentence against THEIR verdicts rather than against the
# set as a whole; whether a wording over-claims is a JUDGEMENT, encoded by
# whether a `retired` phrase is given at all, because a keyword classifier
# cannot measure that and this project has the measurement on record.
SITES_CSV = REPO / "analysis" / "descriptor-absence-sites.csv"


def carrying_sites() -> tuple:
    if not SITES_CSV.exists():
        raise SystemExit(f"{SITES_CSV} is missing; it is the site registry")
    rows = []
    with SITES_CSV.open(encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            rows.append((row["path"], row["phrase"],
                         row["retired"] or None,
                         tuple(row["concerns"].split())))
    if not rows:
        raise SystemExit(f"{SITES_CSV} holds no rows")
    return tuple(rows)


def _squash(text: str) -> str:
    """Whitespace-insensitive comparison.

    The raw version has already been paid for here: a registry keyed on a
    phrase that must not cross a line wrap ends up keyed on a fragment short
    enough to be ambiguous.
    """
    return re.sub(r"\s+", " ", text)


def _contains(haystack: str, needle: str) -> bool:
    """Whitespace-insensitive containment, COLLAPSING runs rather than deleting.

    Deleting whitespace would match a phrase whose words only join if the gaps
    vanish, which is a different string. Extracted from `_site_present` so the
    rule is drivable: on the current text it is INERT -- no registered phrase's
    match depends on it, measured -- so a fixture is the only thing that can
    exercise it, and an inert check reachable no other way is a recorded defect
    class here.
    """
    return _squash(needle) in _squash(haystack)


def _site_present(path: str, phrase: str) -> bool:
    target = REPO / path
    if not target.exists():
        return False
    return _contains(target.read_text(encoding="utf-8"), phrase)


def _first_last_year(by_year) -> tuple:
    years = sorted(int(y) for y in by_year)
    return (years[0], years[-1]) if years else (None, None)


TRIAL_TYPES = frozenset({
    "Clinical Trial", "Randomized Controlled Trial", "Controlled Clinical Trial",
    "Clinical Trial, Phase I", "Clinical Trial, Phase II",
    "Clinical Trial, Phase III", "Clinical Trial, Phase IV",
    "Pragmatic Clinical Trial", "Adaptive Clinical Trial",
})


def _descriptor_row(ui, name, records, major, detail, mech) -> dict:
    by_year = (detail.get("by_year") or {}).get(ui) or {}
    pub_types = (detail.get("pub_types") or {}).get(ui) or {}
    first, last = _first_last_year(by_year)
    trials = sum(v for k, v in pub_types.items() if k in TRIAL_TYPES)
    dated = sum(by_year.values())
    pre = sum(v for k, v in by_year.items() if int(k) < VINTAGE_YEAR)
    vintage_share = None if not dated else pre / dated
    text_arm = (detail.get("text_arm_records") or {}).get(mech)
    overlap = (detail.get("both") or {}).get(f"{mech}|{ui}")
    return {
        "ui": ui, "name": name, "records": records, "major": major,
        "first_year": first, "last_year": last,
        "trials": trials,
        "late_phase": sum(v for k, v in pub_types.items()
                          if k in ("Clinical Trial, Phase III",
                                   "Clinical Trial, Phase IV")),
        "by_year": {str(k): v for k, v in sorted(by_year.items())},
        "vintage_share": vintage_share,
        "vintage_suspect": (None if vintage_share is None
                            else vintage_share >= VINTAGE_SUSPECT_SHARE),
        # The second arm `census-unnamed-modalities.md` reports as unavailable
        # BY CONSTRUCTION for these mechanisms. Both directions, because they
        # answer different questions and neither is a precision figure: the
        # descriptor side is how much of NLM's pool our own words find, the text
        # side how much of our pool NLM filed under the descriptor.
        "text_arm_records": text_arm,
        "overlap": overlap,
        "descriptor_found_by_text": (
            None if not records or overlap is None else overlap / records),
        "text_found_by_descriptor": (
            None if not text_arm or overlap is None else overlap / text_arm),
    }


def assemble(raw: dict) -> dict:
    vocab = raw["vocabulary"]
    detail = raw["detail"]
    unmeasurable = unusable(raw["map"])
    # BOTH ARMS ARE RE-DERIVED HERE, from the committed vocabulary through the
    # live matcher. Nothing about which descriptors exist is read back out of
    # this page's own artifact, so every guard below drives a real derivation.
    found = candidate_uis(vocab["rows"], unmeasurable)
    by_ui = {r["ui"]: r for r in vocab["rows"]}

    mechanisms = {}
    for mech, entry in sorted(unmeasurable.items()):
        hits = found[mech]
        arm_a_rows = [
            _descriptor_row(ui, by_ui[ui]["name"], by_ui[ui]["records"],
                            by_ui[ui]["major"], detail, mech)
            for ui in sorted(hits["arm_a"],
                             key=lambda u: (-by_ui[u]["records"], u))]
        loose_only = [
            _descriptor_row(ui, by_ui[ui]["name"], by_ui[ui]["records"],
                            by_ui[ui]["major"], detail, mech)
            for ui in sorted(set(hits["arm_b"]) - set(hits["arm_a"]),
                             key=lambda u: (-by_ui[u]["records"], u))]
        named = {}
        for name, info in sorted(hits["named"].items()):
            row = (None if info["ui"] is None else
                   _descriptor_row(info["ui"], name, info["records"],
                                   info["major"], detail, mech))
            named[name] = {"field": info["field"], "exists": info["ui"] is not None,
                           "descriptor": row}

        verdict = entry["verdict"]
        usable = [r for r in arm_a_rows if r["records"] >= POOL_FLOOR]
        covering = [n for n, v in named.items() if v["field"] == "covering"]
        # THE CHECKABLE HALF OF THE HANDWRITTEN VERDICT. Each condition exists
        # because its absence is a way a verdict can be wrong, and the
        # `name_exact` ones are what refute an absence claim.
        checks = {
            # A name the map invented matches nothing -- `atlas_thesis_rank.py`
            # found 5 of 17 of its own stems in that state.
            "every_named_descriptor_exists": all(
                v["exists"] for v in named.values()),
            "verdict_is_known": verdict in VERDICTS,
            # A MECHANISM WITH NO KEYWORDS GETS `absent` FOR FREE, because both
            # arms are rules over the mechanism's own vocabulary and an empty
            # vocabulary accepts nothing. Without this the strongest verdict on
            # the page would be available to any mechanism that simply has no
            # terms, which is a vacuity rather than a measurement.
            "mechanism_has_keywords": bool(hits["keywords"]),
            # `absent` and `over-broad` BOTH require that nothing names the
            # mechanism, so a found leaf refutes either.
            "no_name_exact_where_absence_claimed": (
                verdict == "name-exact-leaky" or not arm_a_rows),
            # `absent` additionally requires that nothing CONTAINS it.
            "absent_names_no_covering": (verdict != "absent" or not covering),
            # `over-broad` is meaningless without something materially wider.
            "over_broad_has_a_wide_cover": (
                verdict != "over-broad"
                or any(named[n]["descriptor"]
                       and named[n]["descriptor"]["records"]
                       >= POOL_FLOOR * OVER_BROAD_MULTIPLE
                       for n in covering)),
            # A leaf found by the scan must be recorded in the file the
            # consumers read, or the next reader is back where this started.
            "leaves_declared": (
                verdict != "name-exact-leaky"
                or ({r["name"] for r in arm_a_rows}
                    <= {n for n, v in named.items() if v["field"] == "leaves"})),
            "name_exact_has_a_leaf": (
                verdict != "name-exact-leaky" or bool(arm_a_rows)),
        }
        mechanisms[mech] = {
            "verdict": verdict,
            "reason": entry["reason"].strip(),
            "keywords": hits["keywords"],
            "named": named,
            "arm_a": arm_a_rows,
            "arm_b_loose_only": loose_only,
            "n_arm_a": len(arm_a_rows),
            "n_arm_b": len(hits["arm_b"]),
            # ARM A IS NOT A SUBSET OF ARM B IN PRINCIPLE, and the page says so
            # rather than assuming the convenient direction: arm A is a
            # substring rule, so a keyword can match INSIDE a token
            # (`protac` within a species name) where arm B's token-subset rule
            # cannot. Measured per mechanism instead of asserted, because "an
            # absence under both arms" is only as strong as the overlap.
            "arm_a_not_in_arm_b": sorted(
                set(hits["arm_a"]) - set(hits["arm_b"])),
            # Volume-measurable is NOT recall-measurable, and the page reports
            # them apart because collapsing them is half of this defect.
            "volume_measurable": bool(usable),
            "volume_pool": max((r["records"] for r in usable), default=0),
            "recall_measurable": False,
            "checks": checks,
            "checks_pass": all(checks.values()),
        }

    refuted = sorted(m for m, v in mechanisms.items()
                     if v["verdict"] == "name-exact-leaky")
    absent = sorted(m for m, v in mechanisms.items() if v["verdict"] == "absent")
    over_broad = sorted(m for m, v in mechanisms.items()
                        if v["verdict"] == "over-broad")

    sites = []
    for path, phrase, retired, concerns in carrying_sites():
        present = _site_present(path, phrase)
        gone = retired is None or not _site_present(path, retired)
        named_verdicts = sorted({mechanisms[m]["verdict"] for m in concerns
                                 if m in mechanisms})
        over_claimed = any(mechanisms[m]["verdict"] != "absent"
                           for m in concerns if m in mechanisms)
        sites.append({
            "path": path, "phrase": phrase, "retired": retired,
            "corrected": retired is not None,
            "concerns": list(concerns), "present": present,
            "retired_absent": gone,
            "verdicts_named": named_verdicts,
            # MISCLASSIFIED when the judgement disagrees with the measurement:
            # corrected while every mechanism it names is `absent` (so there
            # was nothing to correct), or left as a control while one of them
            # has a descriptor (so the flattening is still shipped).
            "misclassified": (retired is not None) != over_claimed,
        })

    return {
        "census_records": vocab["scanned"],
        "map_corpus_records": map_corpus_size(),
        "census_descriptors": len(vocab["rows"]),
        "detail_records": detail["scanned"],
        "pool_floor": POOL_FLOOR,
        "over_broad_multiple": OVER_BROAD_MULTIPLE,
        # DECLARED INERT rather than left looking load-bearing. Measured, not
        # assumed: every leaf found clears the pool floor and every covering
        # descriptor clears the breadth bar by orders of magnitude, so deleting
        # either filter would change no output on this data.
        "inert_thresholds": [
            [name, value, why] for name, value, why, is_inert in (
                ("POOL_FLOOR", POOL_FLOOR,
                 "every leaf found clears it", all(
                     r["records"] >= POOL_FLOOR
                     for v in mechanisms.values() for r in v["arm_a"])),
                ("OVER_BROAD_MULTIPLE x POOL_FLOOR",
                 OVER_BROAD_MULTIPLE * POOL_FLOOR,
                 "every covering descriptor clears it", all(
                     named[n]["descriptor"]
                     and named[n]["descriptor"]["records"]
                     >= POOL_FLOOR * OVER_BROAD_MULTIPLE
                     for v in mechanisms.values()
                     for named in [v["named"]]
                     for n in named
                     if named[n]["field"] == "covering")),
            ) if is_inert],
        "arm_b_generic": sorted(ARM_B_GENERIC),
        "mechanisms": mechanisms,
        "n_unmeasurable": len(mechanisms),
        "absent": absent,
        "over_broad": over_broad,
        "refuted_absence": refuted,
        "all_checks_pass": all(v["checks_pass"] for v in mechanisms.values()),
        "sites": sites,
        "n_sites_corrected": sum(1 for s in sites if s["corrected"]),
        "n_sites_control": sum(1 for s in sites if not s["corrected"]),
        "sites_missing": [s["path"] for s in sites if not s["present"]],
        "sites_retired_still_present": [s["path"] for s in sites
                                        if not s["retired_absent"]],
        "sites_misclassified": [s["path"] for s in sites if s["misclassified"]],
        # A candidate with no committed detail means the vocabulary moved and
        # pass 2 was not re-run. A state to report, not a crash.
        "vintage_suspect": sorted(
            m for m, v in mechanisms.items()
            if any(r["vintage_suspect"] for r in v["arm_a"])),
        "arm_a_outside_arm_b": sorted(
            ui for v in mechanisms.values() for ui in v["arm_a_not_in_arm_b"]),
        "detail_missing": sorted(
            r["ui"] for v in mechanisms.values()
            for r in v["arm_a"] + v["arm_b_loose_only"]
            if r["ui"] not in (detail.get("descriptor_records") or {})),
    }


def _n(x) -> str:
    return "n/a" if x is None else f"{x:,}"


def _pct(x) -> str:
    return "n/a" if x is None else f"{100 * x:.1f}%"


def _descriptor_table(rows) -> list:
    L = ["| descriptor | UI | census records | major topic | years | trials | "
         "late phase |", "|---|---|--:|--:|---|--:|--:|"]
    for r in rows:
        span = ("n/a" if r["first_year"] is None
                else (str(r["first_year"]) if r["first_year"] == r["last_year"]
                      else f"{r['first_year']}-{r['last_year']}"))
        L.append(f"| `{r['name']}` | {r['ui']} | {_n(r['records'])} | "
                 f"{_n(r['major'])} | {span} | {_n(r['trials'])} | "
                 f"{_n(r['late_phase'])} |")
    return L


def _mechanism_section(mech, v) -> list:
    L = [f"### `{mech}` -- {v['verdict']}", ""]
    L.append(v["reason"])
    L.append("")
    if v["arm_a"]:
        n = len(v["arm_a"])
        L.append(f"**Arm A finds {n} name-exact "
                 f"descriptor{'' if n == 1 else 's'}.** "
                 f"The production matcher accepts "
                 f"{'this descriptor name' if n == 1 else 'these descriptor names'}"
                 f" from this mechanism's own keywords, so "
                 f"{'it names' if n == 1 else 'each names'} the mechanism in "
                 f"the mechanism's own words.")
        L.append("")
        L += _descriptor_table(v["arm_a"])
        L.append("")
        pool = v["volume_pool"]
        if v["volume_measurable"]:
            L.append(f"So a census VOLUME count is available at a pool of "
                     f"{_n(pool)} records against a floor of "
                     f"{POOL_FLOOR}. A recall measurement is not: every leaf "
                     f"here is LEAKY by `mechanism_recall.py`'s own "
                     f"classifier, since the descriptor's name is the keyword, "
                     f"so agreement with it would be near-tautological.")
            suspect = [r for r in v["arm_a"] if r["vintage_suspect"]]
            if suspect:
                L.append("")
                L.append(
                    "**But the descriptor is older than the therapy it is "
                    "being read as**, and a name-exact leaf is not exempt from "
                    "that: "
                    + ", ".join(
                        f"`{r['name']}` carries {_pct(r['vintage_share'])} of "
                        f"its pool before {VINTAGE_YEAR}"
                        for r in suspect)
                    + f", above the {_pct(VINTAGE_SUSPECT_SHARE)} at which "
                      f"`census-unnamed-modalities.md` marks its keyword rows "
                      f"suspect. None of these modalities existed then, so that "
                      f"share is the descriptor matching an older use of its "
                      f"own words and the volume count is an upper bound rather "
                      f"than a measurement of the modality.")
        else:
            L.append(f"The pool is {_n(pool)} records against a floor of "
                     f"{POOL_FLOOR}, so the descriptor is found but cannot "
                     f"carry a volume count either.")
        L.append("")
    else:
        L.append("**Arm A finds nothing**: no descriptor in the census "
                 "vocabulary carries a name this mechanism's own keywords "
                 "accept.")
        L.append("")
    named = v["named"]
    if named:
        L.append("Descriptors the map names, measured rather than trusted:")
        L.append("")
        L.append("| descriptor | map calls it | exists in census | records |")
        L.append("|---|---|---|--:|")
        for name, info in named.items():
            row = info["descriptor"]
            L.append(f"| `{name}` | {info['field']} | "
                     f"{'yes' if info['exists'] else '**NO**'} | "
                     f"{_n(row['records']) if row else 'n/a'} |")
        L.append("")
    loose = v["arm_b_loose_only"]
    if loose:
        L.append(f"Arm B additionally reaches {len(loose)} descriptor"
                 f"{'' if len(loose) == 1 else 's'} arm A "
                 f"does not. Arm B is a recall instrument with poor precision "
                 f"by design, so these are shown with their counts to be "
                 f"judged rather than taken as candidates:")
        L.append("")
        L.append("| descriptor | records |")
        L.append("|---|--:|")
        for r in loose[:8]:
            L.append(f"| `{r['name']}` | {_n(r['records'])} |")
        L.append("")
    elif v["arm_a"]:
        L.append("Arm B reaches nothing arm A does not, so the two rules agree "
                 "on exactly which descriptors name this mechanism.")
        L.append("")
    else:
        L.append("Arm B reaches nothing either, so the absence is an absence "
                 "under both rules.")
        L.append("")
    if v["arm_a_not_in_arm_b"]:
        L.append("Arm A reaches "
                 + ", ".join(f"`{u}`" for u in v["arm_a_not_in_arm_b"])
                 + " and arm B does not, which is possible because arm A is a "
                   "substring rule and can match inside a token. Those hits "
                   "rest on one arm.")
        L.append("")
    failed = sorted(k for k, ok in v["checks"].items() if not ok)
    if not v["checks"]["mechanism_has_keywords"]:
        L.append("**This mechanism has no keyword vocabulary**, so both arms "
                 "are rules with nothing to apply and an absence here would be "
                 "a property of the empty set rather than a measurement.")
        L.append("")
    if failed:
        L.append("**The declared verdict fails its own checks: "
                 + ", ".join(f"`{k}`" for k in failed) + ".**")
        L.append("")
    return L


def render(d: dict) -> str:
    L = [f"# Is \"MeSH has no descriptor\" a statement about MeSH?\n"]
    L.append(
        f"Generated by `scripts/descriptor_absence_audit.py` over the census's "
        f"own descriptor vocabulary: **{d['census_descriptors']:,} distinct "
        f"MeSH descriptors** across {d['census_records']:,} indexed cancer "
        f"records. {d['n_unmeasurable']} mechanisms are reported throughout "
        f"this project as *not measurable at census scale*, and "
        f"{d['n_sites_corrected'] + d['n_sites_control']} sites across the "
        f"generated pages and the manuscript give the reason as MeSH having no "
        f"descriptor for them. Measured below, "
        f"{d['n_sites_corrected']} of those over-claimed and are corrected "
        f"here; {d['n_sites_control']} were about mechanisms whose verdict "
        f"really is `absent` and are left exactly as written.\n")
    L.append(
        f"**That is not what any of those sites measures.** Every consumer "
        f"derives it from `analysis/mesh-mechanism-map.yaml` listing no usable "
        f"descriptor -- a handwritten file whose own header records that its "
        f"counts were verified against a corpus of "
        f"{d['map_corpus_records']:,} records. What is measured is *this file "
        f"names none*; what is carried is *MeSH has none*, about a population "
        f"{d['census_records'] // d['map_corpus_records']:,} times larger.\n")
    L.append(
        f"So the check here is as wide as the claim. **Arm A** is the "
        f"production rule: a descriptor is a candidate when the mechanism's own "
        f"`MECHANISM_KEYWORDS` accept its NAME under `tag_articles`' "
        f"`text_matches_keyword`, the same live matcher "
        f"`mechanism_recall.py` uses to decide leakiness. **Arm B** is a "
        f"looser token-subset control with plural stripping and "
        f"{len(d['arm_b_generic'])} generic words dropped from the keyword "
        f"side, so an absence is an absence under both and arm A's own blind "
        f"spot -- it is a substring rule, so a plural keyword cannot reach a "
        f"singular descriptor -- is measured rather than assumed away. Neither "
        f"arm is semantic: both are rules over descriptor NAMES, so a "
        f"descriptor meaning the mechanism in different words is invisible to "
        f"them exactly as it is to a hand list. Wider is not complete.\n")

    L.append("## The verdict\n")
    refuted = d["refuted_absence"]
    if refuted:
        # DERIVED rather than asserted, and the distinction matters: the map's
        # own notes did not deny a descriptor for all three -- one said the
        # descriptor existed and gave a pool measured on the frozen corpus --
        # while the sentences the registry retired swept all three into an
        # absence. Counting the retired rows that NAME each one says exactly
        # how far the over-claim reached.
        named_by = {m: sum(1 for s in d["sites"]
                           if s["corrected"] and m in s["concerns"])
                    for m in refuted}
        L.append(
            f"**{len(refuted)} of {d['n_unmeasurable']} are refuted**: the "
            f"census carries a name-exact MeSH descriptor for "
            + ", ".join(f"`{m}` (swept into {named_by[m]} corrected "
                        f"sentence{'' if named_by[m] == 1 else 's'})"
                        for m in refuted)
            + ". Each is volume-measurable at census scale, and each was named "
              "in a sentence saying MeSH has no descriptor for it.\n")
    else:
        L.append(
            f"**No absence assertion is refuted.** Arm A finds a name-exact "
            f"descriptor for none of the {d['n_unmeasurable']}, so the "
            f"whole-vocabulary scan agrees with the hand list and what remains "
            f"of this page is the verdict classes and the stale pools.\n")
    L.append(
        f"The remaining {len(d['absent'])} + {len(d['over_broad'])} split into "
        f"two verdicts that were being reported as one word, and `CLAUDE.md` "
        f"already records that they must stay apart. **Absent** -- nothing "
        f"names the mechanism and nothing contains it: "
        + ", ".join(f"`{m}`" for m in d["absent"])
        + f". **Over-broad** -- a descriptor contains the mechanism but means "
          f"far more, so a count of it over-counts, which was never absence: "
        + ", ".join(f"`{m}`" for m in d["over_broad"]) + ".\n")
    L.append("| mechanism | verdict | arm A | arm B | volume-measurable | "
             "pool | recall-measurable |")
    L.append("|---|---|--:|--:|---|--:|---|")
    for mech, v in d["mechanisms"].items():
        L.append(f"| `{mech}` | {v['verdict']} | {v['n_arm_a']} | "
                 f"{v['n_arm_b']} | "
                 f"{'yes' if v['volume_measurable'] else 'no'} | "
                 f"{_n(v['volume_pool']) if v['volume_pool'] else '--'} | "
                 f"{'yes' if v['recall_measurable'] else 'no'} |")
    L.append("")
    suspect = sorted(m for m, v in d["mechanisms"].items()
                     if any(r["vintage_suspect"] for r in v["arm_a"]))
    if suspect:
        L.append(
            f"**A name-exact descriptor is not automatically a clean count**, "
            f"and one of the refuted rows shows it: "
            + ", ".join(f"`{m}`" for m in suspect)
            + f" carries a descriptor whose own pool reaches back past "
              f"{VINTAGE_YEAR}, when the modality did not exist, so its volume "
              f"is an upper bound. That is the same vintage trap "
              f"`census-unnamed-modalities.md` measures on the keyword arm, and "
              f"finding a descriptor does not escape it.\n")
    L.append(
        f"**Volume-measurable and recall-measurable are different claims, and "
        f"collapsing them is half of why this went unnoticed.** A descriptor "
        f"whose name trips the mechanism's keywords is LEAKY by "
        f"`mechanism_recall.py`'s own classifier, so it cannot serve as an "
        f"independent reference label and that page's headline is right to "
        f"exclude it. It can carry a census volume count perfectly well, which "
        f"is what `census_mechanism_profile.py` reports for every measurable "
        f"mechanism and what the sentences below deny is possible. The recall "
        f"verdict survives for all {d['n_unmeasurable']}; the absence claim "
        f"does not.\n")

    L.append("## The handwritten verdict, checked where it is checkable\n")
    L.append(
        "Whether a broad descriptor CONTAINS a mechanism or merely RESEMBLES "
        "it is a judgement about MeSH semantics that no rule over names can "
        "make, so the verdict stays in the handwritten map. What this page "
        "checks is the half that is checkable: every descriptor the map names "
        "must exist in the census vocabulary -- `atlas_thesis_rank.py` found "
        "several of its own stems matching nothing -- `absent` is refused where "
        "a covering descriptor is named or a name-exact one is found, "
        "`over-broad` is refused where a name-exact one is found, and "
        f"`name-exact-leaky` is refused unless arm A finds one and the map "
        f"lists it.\n")
    L.append(f"All checks pass: **{d['all_checks_pass']}**.\n")
    inert = d["inert_thresholds"]
    if inert:
        L.append(
            f"**{len(inert)} of those thresholds are INERT on the current data "
            f"and are declared rather than left to look load-bearing**: "
            + ", ".join(
                f"`{name}` at {value:,} ({why})"
                for name, value, why in inert)
            + ". An inert floor reachable only through a fixture is a recorded "
              "defect class here, so each is driven by a doctored vocabulary in "
              "`tests/test_descriptor_absence_audit.py` and its committed VALUE "
              "is bounded -- because a sweep that deletes a filter nothing "
              "reaches changes no output at all.\n")

    L.append("## Per mechanism\n")
    for mech, v in d["mechanisms"].items():
        L += _mechanism_section(mech, v)

    L.append("## The second arm, for the mechanisms that now have one\n")
    pairs = [(m, r) for m, v in d["mechanisms"].items() for r in v["arm_a"]]
    if pairs:
        L.append(
            "`analysis/census-unnamed-modalities.md` rests on the "
            "descriptor-versus-text agreement check being unavailable for "
            "these mechanisms *by construction* -- \"there is no second arm, "
            "because the absence of a second arm is what puts them on this "
            "list\". For the ones with a name-exact leaf that is no longer "
            "true, so the agreement is measured here in both directions. "
            "Neither figure is a precision measurement: the descriptor side is "
            "how much of NLM's pool this project's own words find, the text "
            "side how much of this project's pool NLM filed under the "
            "descriptor.\n")
        L.append("| mechanism | descriptor | descriptor pool | text-arm pool | "
                 "both | of descriptor pool | of text-arm pool |")
        L.append("|---|---|--:|--:|--:|--:|--:|")
        for mech, r in pairs:
            L.append(f"| `{mech}` | `{r['name']}` | {_n(r['records'])} | "
                     f"{_n(r['text_arm_records'])} | {_n(r['overlap'])} | "
                     f"{_pct(r['descriptor_found_by_text'])} | "
                     f"{_pct(r['text_found_by_descriptor'])} |")
        L.append("")
        scored = [(m, r) for m, r in pairs
                  if r["descriptor_found_by_text"] is not None]
        if scored:
            best = max(scored, key=lambda x: x[1]["descriptor_found_by_text"])
            L.append(
                f"The direction that carries information is the first column: "
                f"`{best[1]['name']}` is "
                f"{_pct(best[1]['descriptor_found_by_text'])} contained in "
                f"`{best[0]}`'s keyword pool, so where NLM applies the "
                f"descriptor this project's own words agree -- while the "
                f"keyword pool is {best[1]['text_arm_records'] / best[1]['records']:.1f}x "
                f"larger, which is the gap a single-arm count was reporting "
                f"without a check. The second column is NOT a recall figure for "
                f"the tagger: a keyword pool legitimately contains papers NLM "
                f"indexed under something else, and for a descriptor this "
                f"recent it also contains papers indexed before the descriptor "
                f"existed.\n")
    else:
        L.append(
            "No mechanism here has a name-exact leaf, so the second arm stays "
            "unavailable for all of them and that page's justification "
            "stands.\n")

    L.append("## Where the flattened claim is written\n")
    L.append(
        f"{d['n_sites_corrected']} sites stated it in a form at least one "
        f"mechanism refutes and are corrected here; {d['n_sites_control']} "
        f"state the same phrasing about mechanisms whose verdict really is "
        f"`absent` and are left exactly as written -- registered as controls, "
        f"because a registry whose every entry is a defect has not shown where "
        f"the line is, and because the corrected pages must still be allowed "
        f"to say this. Whether a wording over-claims is a JUDGEMENT encoded in "
        f"the registry, not a rule: a keyword classifier cannot measure that, "
        f"and this project has the measurement on record. Each corrected row "
        f"is checked in BOTH directions -- the replacement present AND the "
        f"flattened original gone -- because a correction verified only by its "
        f"replacement cannot tell a fix from a fix shipped beside a re-added "
        f"original.\n")
    L.append("| file | corrected | wording now | mechanisms | their verdicts | "
             "present | old wording gone |")
    L.append("|---|---|---|--:|---|---|---|")
    for site in d["sites"]:
        # Three rows share one path, so the wording has to appear or they are
        # indistinguishable in the table that exists to tell them apart.
        quoted = _squash(site["phrase"])
        if len(quoted) > 64:
            quoted = quoted[:64].rsplit(" ", 1)[0] + " ..."
        L.append(f"| `{site['path']}` | "
                 f"{'yes' if site['corrected'] else 'control'} | "
                 f"{quoted} | "
                 f"{len(site['concerns'])} | "
                 f"{', '.join(site['verdicts_named'])} | "
                 f"{'yes' if site['present'] else '**NO**'} | "
                 f"{'yes' if site['retired_absent'] else '**NO**'} |")
    L.append("")
    if d["sites_retired_still_present"]:
        L.append("**Flattened wording still present at: "
                 + ", ".join(f"`{x}`" for x in d["sites_retired_still_present"])
                 + ".** The correction was not applied, or was reverted.\n")
    if d["sites_missing"]:
        L.append("**Registered phrases no longer present: "
                 + ", ".join(f"`{p}`" for p in d["sites_missing"])
                 + ".** A registry entry pointing at text that has moved "
                   "verifies nothing.\n")
    if d["sites_misclassified"]:
        L.append("**Registry judgement disagrees with the measurement at: "
                 + ", ".join(f"`{p}`" for p in d["sites_misclassified"])
                 + ".**\n")

    L.append("## What this does not establish\n")
    L.append(
        f"**A descriptor existing is not a mechanism being measurable** in the "
        f"sense the corrected sentences need. The deliverable is that the "
        f"absence assertion is false, not that these mechanisms are now "
        f"countable alongside the rest: every leaf found is leaky, so none can "
        f"carry a concordance check, and a recent descriptor's census volume is "
        f"dominated by indexing vintage -- the lesson "
        f"`analysis/atlas-recent-window.md` records, where a third of "
        f"significantly-rising descriptors are flat against a recent "
        f"comparator.\n")
    L.append(
        "**Both arms are name rules, so neither bounds MeSH's semantic "
        "coverage.** A descriptor meaning a mechanism under words the project "
        "does not use is invisible to arm A, arm B and the hand list alike. "
        "What changed is the scope of the check, from the candidates an author "
        "recalled to every descriptor in the census whose name the mechanism's "
        "own production vocabulary accepts.\n")
    L.append(
        "**Nothing here touches the simulation.** No parameter is constrained "
        "and `analysis/identifiability-report.json`'s data-conditioned-headline "
        "count stays where it is.\n")
    L.append(
        f"**And the vocabulary is a snapshot.** It is the descriptor set of "
        f"one PubMed baseline, so an absence verdict here has the same expiry "
        f"the claims it corrects had: `Disulfidptosis` carried no census "
        f"records before 2024. The difference is that this page's verdicts are "
        f"re-derived from a committed vocabulary on every render, so refreshing "
        f"the census is what it takes to re-check them.\n")
    if d["detail_missing"]:
        L.append("**Candidates with no committed detail: "
                 + ", ".join(f"`{u}`" for u in d["detail_missing"])
                 + ".** The vocabulary moved and the second pass was not "
                   "re-run.\n")
    return "\n".join(L) + "\n"


def verdicts(mechanisms=None) -> dict:
    """The committed per-mechanism verdicts, for a consumer that states one.

    FAIL-OPEN: a consumer must not be blocked by this page's artifact being
    absent, and an empty answer is read by every caller as "say nothing extra"
    rather than as "nothing exists". The alternative -- failing closed -- would
    make every census generator that states a verdict unrunnable in a tree that has not run the
    census scan, which is every CI tree.
    """
    if not OUT_JSON.exists():
        return {}
    stored = json.loads(OUT_JSON.read_text()).get("mechanisms") or {}
    out = {m: v["verdict"] for m, v in stored.items()}
    if mechanisms is None:
        return out
    return {m: v for m, v in out.items() if m in set(mechanisms)}


def absence_note(mechanisms) -> list:
    """One derived sentence a consumer appends instead of saying "no descriptor".

    The sentence a consumer used to write -- "MeSH has no descriptor for these"
    -- is a claim about MeSH derived from the map, and for the mechanisms with a
    name-exact leaf it is false. This splits whatever set the caller has by the
    verdict the audit measured, so the caller states only what its own set
    supports and cannot drift from this page.
    """
    found = verdicts(mechanisms)
    if not found:
        return []
    groups = {}
    for mech, verdict in sorted(found.items()):
        groups.setdefault(verdict, []).append(mech)
    parts = []
    if groups.get("absent"):
        parts.append("no descriptor names or contains "
                     + ", ".join(f"`{m}`" for m in groups["absent"]))
    if groups.get("over-broad"):
        parts.append("a descriptor CONTAINS "
                     + ", ".join(f"`{m}`" for m in groups["over-broad"])
                     + " while meaning far more, so a count of it over-counts")
    if groups.get("name-exact-leaky"):
        parts.append("the census carries a NAME-EXACT descriptor for "
                     + ", ".join(f"`{m}`" for m in groups["name-exact-leaky"])
                     + ", so these are volume-measurable at census scale and "
                       "unusable only as a recall reference")
    if not parts:
        return []
    return ["Measured against the census's whole descriptor vocabulary in "
            "`analysis/descriptor-absence-audit.md`: " + "; ".join(parts)
            + f". Those {len(groups)} verdicts license different claims and "
              f"collapsing them into \"MeSH has no descriptor\" is what that "
              f"page corrects.\n"]

def _split_stored(stored):
    """Recover `assemble`'s input for the freshness gate's round trip.

    RE-READS the committed inputs rather than echoing stored derived fields out
    of the artifact. An echo would make every derived-field guard in
    `tests/test_descriptor_absence_audit.py` inert, which is a defect this
    repository has paid for repeatedly and now gates for.
    """
    return (load(),)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--measure", action="store_true",
                    help="re-scan the census (two passes over every shard)")
    ap.add_argument("--render-only", action="store_true",
                    help="re-derive from the committed vocabulary and detail")
    ap.add_argument("--stride", type=int, default=1,
                    help="sample every Nth shard (development only)")
    args = ap.parse_args()
    if args.measure:
        measure(args.stride)
    # The two branches are the SAME call, deliberately: every remaining input
    # is committed and holds raw counts only, so the default run already
    # re-derives the whole page -- both arms included -- and there is no stored
    # derived field for either path to echo. Written out rather than collapsed
    # because `tests/test_artifact_freshness.py` keys two real round trips on
    # seeing it, and because a derived cache added later must not slip in
    # behind a collapsed branch.
    if args.render_only:
        d = assemble(load())
    else:
        d = assemble(load())
    OUT_JSON.write_text(json.dumps(d, indent=2, sort_keys=True) + "\n")
    OUT_MD.write_text(render(d))
    print(f"wrote {OUT_JSON}")
    print(f"wrote {OUT_MD}")
    print(f"  {d['census_descriptors']:,} descriptors over "
          f"{d['census_records']:,} records")
    print(f"  absent={len(d['absent'])} over-broad={len(d['over_broad'])} "
          f"name-exact-leaky={len(d['refuted_absence'])} "
          f"checks_pass={d['all_checks_pass']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
