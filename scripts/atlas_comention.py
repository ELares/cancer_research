#!/usr/bin/env python3
"""Atlas: sentence-level entity co-mention over open-access full text (#ATLAS).

WHY
---
The relation graph is built from PubTator3's ABSTRACT-level extraction, and its
edge recall is too low for anything that reasons from absence. Measured directly:
literature-based discovery over that graph returned ERK, caspase-3 and cyclin D1
as top "undiscovered" GPX4 links, and 12 of 12 checkable candidates turned out to
be already co-mentioned in PubMed -- GPX4 and caspase-3 share 236 abstracts and
no graph edge. Absence in that graph is not informative.

This raises recall using material already on disk: 520,143 open-access full
texts. It does NOT need an LLM. The entity vocabulary is already normalized --
529,652 surface forms mapped to 168,385 NCBI Gene and MeSH identifiers, from
PubTator3's own annotations -- so the same dictionary can be applied to full text
directly, deterministically and offline.

WHAT IT PRODUCES, AND WHAT IT DOES NOT
--------------------------------------
Co-mention within a SENTENCE. Two entities named in one sentence are far more
likely related than two named in the same 30,000-word paper, but this is still
CO-MENTION, not a typed relation: there is no predicate, no direction, and no
polarity. A sentence saying "X does not inhibit Y" yields the same edge as one
saying it does.

So this layer complements the PubTator relations rather than replacing them:
  * PubTator gives PRECISION and predicates over abstracts;
  * this gives RECALL over full text, which is what absence-based reasoning needs.

DISAMBIGUATION
--------------
The alias map is the weak point. It contains short and ambiguous surface forms
that would produce enormous numbers of false mentions in running text -- `CAR`,
`AGE`, `ICE`, single letters, common words. Aliases are therefore filtered:
minimum length, at least one character that is not lowercase-alphabetic OR a
multi-word form, and an explicit stoplist of English words that happen to be gene
symbols. The filter is reported so its aggressiveness is visible, and it is
deliberately strict: a missing entity costs recall, a spurious one poisons every
downstream count.

That filter is about SHAPE, and shape says nothing about whether a form means one
thing. Length happens to catch some collisions -- `psa`, `p21`, `p62` and `er` are
all below the four-character minimum -- but that is luck, not disambiguation, and
75 forms `scripts/atlas_ambiguity.py` measured as sense collisions came through
it. Among them `cox-2` and `fsp1`, whose majority votes are WRONG: left alone this
layer counts full-text COX-2 co-mentions against mitochondrial cytochrome c
oxidase rather than PTGS2, and FSP1 against a spastic-paraplegia gene rather than
the ferroptosis suppressor.

So a second, SENSE filter now runs after it: a colliding form is redirected to its
measured cancer-domain sense where one exists (1 of the 75) and dropped otherwise
(the other 74). Both counts are reported at build time. Dropping costs recall,
which this module already treats as the cheaper loss.

AUDITABILITY
------------
The build also writes `comention/audit-sample.jsonl.gz`: a uniform reservoir
sample of matched sentences with their resolved entity names. This layer's
precision has never been measured, and the alias map is its acknowledged weak
point -- an unauditable weak point is indistinguishable from a sound one. The
sample is drawn across the whole run rather than from the first shard, because
shards are ordered by PMCID and a prefix would sample the oldest literature.

Usage:
    python scripts/atlas_comention.py --limit 1     # one shard, smoke test
    python scripts/atlas_comention.py               # all shards
    python scripts/atlas_comention.py --status
"""

import argparse
import collections
import gzip
import json
import pickle
import random
import re
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from atlas_baseline import atlas_root  # noqa: E402
from atlas_fulltext import fulltext_root  # noqa: E402
from atlas_graph import _ambiguity, load_index  # noqa: E402

MIN_ALIAS_LEN = 4
# A surface form must be attested this often across the census before it may be
# matched in running text. MEASURED, not guessed: on a 257-mention labelled
# sample, correct matches have a median census support of 1,067 and wrong ones
# 16. The graph's own MIN_MENTIONS of 3 is the proximate cause of the junk --
# it admits a generic phrase mis-annotated in a handful of abstracts, and the
# majority vote then assigns it an identifier.
MIN_ALIAS_SUPPORT = 50
# A surface form must also account for this share of its identifier's total
# mentions. Support alone does not catch the largest false-positive families:
# `tumor cells` has 312 census mentions and `t cell` 358, comfortably past the
# support bar, but they are 2.8% and 0.8% of Glucagonoma and Lymphoma T-Cell
# respectively -- a rare mis-annotation, not a name. Measured on 152 labelled
# false positives this kills 91.4% of them, with their median share at 0.33%,
# while real names sit at 46-98%.
#
# Both thresholds were chosen on that same labelled sample, so they are a FIT.
# A fresh audit sample after the next rebuild is what would confirm them.
MIN_ALIAS_SHARE = 0.05
MAX_NGRAM = 5
MIN_SENT_ENTITIES = 2
MAX_SENT_ENTITIES = 12   # a 40-entity sentence is a table or a gene list, not a claim

# English words that are also gene symbols or MeSH surface forms. Matching these
# in running text produces vastly more noise than signal.
STOPWORD_ALIASES = {
    "cancer", "tumor", "tumour", "cell", "cells", "gene", "genes", "protein",
    "proteins", "patients", "study", "control", "controls", "human", "mice",
    "mouse", "rat", "rats", "damage", "impact", "camp", "face", "lung", "skin",
    "liver", "brain", "breast", "colon", "blood", "water", "light", "heavy",
    "large", "small", "type", "types", "group", "groups", "level", "levels",
    "time", "times", "rate", "rates", "risk", "score", "stage", "grade",
    "mass", "band", "chip", "spin", "star", "arms", "beta", "alpha", "gamma",
    "delta", "kappa", "lambda", "sigma", "theta", "omega", "clock", "sleep",
    "rest", "size", "shape", "form", "line", "lines", "set", "sets", "map",
    "maps", "not", "and", "the", "for", "with", "from", "this", "that",
}

_SENT_SPLIT = re.compile(r"(?<=[.!?])\s+")
_TOKEN = re.compile(r"[A-Za-z0-9][A-Za-z0-9\-]*")


# Connectives that make an alias a sentence FRAGMENT rather than a name. PubTator's
# surface forms include things like "and sensory" (from "motor and sensory
# neuropathy"), which would match constantly in running text.
_EDGE_WORDS = {"and", "or", "of", "the", "a", "an", "in", "on", "for", "with",
               "to", "by", "at", "from", "as", "is", "are", "was", "were"}


def usable_alias(a: str) -> bool:
    """Keep only surface forms specific enough to match safely in running text."""
    if len(a) < MIN_ALIAS_LEN or a in STOPWORD_ALIASES:
        return False
    # A trailing or leading hyphen is a tokenizer stub, not a name. The token
    # pattern leaves `ifn-` behind from "IFN-gamma", and `ifn-` resolves to
    # IFNA1 -- so interferon gamma was being counted as interferon alpha.
    if a.startswith("-") or a.endswith("-"):
        return False
    if " " in a:
        words = a.split()
        # a fragment that starts or ends on a connective is not a name
        if words[0] in _EDGE_WORDS or words[-1] in _EDGE_WORDS:
            return False
        return True
    # A single token used to require a digit or hyphen, as a proxy for "looks
    # like a symbol rather than an English word". The proxy was expensive and
    # unnecessary: it dropped every purely alphabetic drug name -- `erastin`,
    # `cisplatin`, `pembrolizumab` -- while catching only 3 of 152 measured
    # false positives, because the real offenders are multi-word and were
    # exempt from it anyway.
    #
    # Specificity is now MEASURED downstream, by census support and by the
    # share of its identifier's mentions a form accounts for. Removing the
    # proxy kills the same 141 of 152 false positives and retains strictly more
    # true ones, so it is deleted rather than layered on top of.
    return True


def build_alias_map(idx: dict) -> tuple:
    """Usable surface form -> identifier, with sense collisions handled.

    `usable_alias` is a SHAPE filter -- length, connectives, digits -- and shape
    says nothing about whether a form means one thing. Its length rule happens to
    exclude `psa`, `p21`, `p62` and `er`, but that is luck rather than
    disambiguation: 75 forms measured as sense collisions by
    `scripts/atlas_ambiguity.py` pass it, including `cox-2` and `fsp1`, whose
    majority votes are wrong. Unfixed, this layer attributes full-text COX-2
    co-mentions to mitochondrial cytochrome c oxidase and FSP1 co-mentions to a
    spastic-paraplegia gene.

    So a blocklisted form is either redirected to its curated cancer-domain
    sense, where one has been MEASURED (`analysis/atlas-domain-sense-validation.md`
    puts those at 89.6%-100% of declaring papers), or dropped. Dropping costs
    recall; keeping a wrong identifier poisons every downstream count, and this
    module's own docstring already chooses recall as the cheaper loss.
    """
    blocked, domain = _ambiguity()
    support = idx.get("alias_support") or {}
    ident_tot = idx.get("ident_mentions") or {}
    # Per-(form, identifier) counts. Falls back to the cross-sense total on an
    # index built before this field existed, which reproduces the old (wrong)
    # behaviour rather than crashing -- the rebuild is what fixes it.
    ident_support = idx.get("alias_ident_support") or {}
    out, redirected, dropped, thin, minority = {}, 0, 0, 0, 0
    for a, i in idx["alias"].items():
        if not usable_alias(a):
            continue
        # The shape filter exempts EVERY multi-word form from the specificity
        # test it applies to single tokens, which is how `tumor cells`,
        # `overall survival` and `et al` reached the matcher. Support is the
        # test that catches them: 132 of 152 measured false positives were
        # multi-word, and they separate from true matches by two orders of
        # magnitude of census attestation.
        if support and support.get(a, 0) < MIN_ALIAS_SUPPORT:
            thin += 1
            continue
        if ident_tot:
            # The numerator must be this form's count FOR THIS IDENTIFIER, not
            # its total across every sense it carries. Using the cross-sense
            # total made the ratio exceed 1 for ambiguous forms and inverted the
            # filter's intent -- see the note in atlas_graph.build_index.
            share = ident_support.get(a, support.get(a, 0)) / max(1, ident_tot.get(i, 0))
            if share < MIN_ALIAS_SHARE:
                minority += 1
                continue
        if a in blocked:
            if a in domain:
                out[a] = domain[a]["id"]
                redirected += 1
            else:
                dropped += 1
            continue
        out[a] = i
    return out, {"redirected": redirected, "dropped_ambiguous": dropped,
                 "dropped_thin": thin, "dropped_minority": minority}


def sentence_entities(sentence: str, alias: dict) -> set:
    """Identifiers named in one sentence, by longest-match n-gram lookup."""
    toks = _TOKEN.findall(sentence.lower())
    found = set()
    n = len(toks)
    i = 0
    while i < n:
        hit = None
        for size in range(min(MAX_NGRAM, n - i), 0, -1):
            gram = " ".join(toks[i:i + size])
            ident = alias.get(gram)
            if ident:
                hit = (ident, size)
                break
        if hit:
            found.add(hit[0])
            i += hit[1]
        else:
            i += 1
    return found


# A uniform sample of matched sentences, written alongside the counts so this
# layer's precision can be AUDITED. Without it the build emits pair counts and
# nothing that would let anyone check whether a match is real -- the alias map
# is this module's acknowledged weak point, and an unauditable weak point is
# indistinguishable from a sound one.
AUDIT_SAMPLE = 400
_AUDIT_SEED = 20260803


def process_shard(path: Path, alias: dict, pairs: collections.Counter,
                  audit=None, rng=None) -> dict:
    docs = sents = mentions = kept_sents = 0
    with gzip.open(path, "rt", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            rec = json.loads(line)
            docs += 1
            for sentence in _SENT_SPLIT.split(rec.get("text", "")):
                if len(sentence) < 40 or len(sentence) > 1200:
                    continue
                sents += 1
                ents = sentence_entities(sentence, alias)
                if not (MIN_SENT_ENTITIES <= len(ents) <= MAX_SENT_ENTITIES):
                    continue
                kept_sents += 1
                mentions += len(ents)
                ordered = sorted(ents)
                for a in range(len(ordered)):
                    for b in range(a + 1, len(ordered)):
                        pairs[(ordered[a], ordered[b])] += 1
                if audit is not None:
                    # Algorithm R over every kept sentence in the whole run, so
                    # the sample is uniform across shards rather than a prefix
                    # of the first one.
                    audit["seen"] += 1
                    row = {"pmid": rec.get("pmid"), "sentence": sentence,
                           "entities": ordered}
                    if len(audit["rows"]) < AUDIT_SAMPLE:
                        audit["rows"].append(row)
                    else:
                        j = rng.randrange(audit["seen"])
                        if j < AUDIT_SAMPLE:
                            audit["rows"][j] = row
    return {"docs": docs, "sentences": sents, "kept_sentences": kept_sents,
            "mentions": mentions}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--status", action="store_true")
    ap.add_argument("--rebuild", action="store_true",
                    help="reprocess every shard from scratch. Clears the manifest "
                         "AND the pair table together, which is the only safe way: "
                         "the run merges its results into any existing table, so "
                         "reprocessing shards that are already counted there would "
                         "double every pair.")
    ap.add_argument("--recent-first", action="store_true",
                    help="process the highest PMCID shards first (modern literature)")
    args = ap.parse_args()

    root = atlas_root()
    out_dir = root / "comention"
    out_dir.mkdir(parents=True, exist_ok=True)
    man_path = out_dir / "manifest.json"
    man = json.loads(man_path.read_text()) if man_path.exists() else {"shards": {}}

    if args.status:
        done = [k for k, v in man["shards"].items() if v.get("done")]
        print(f"shards done   : {len(done)}")
        print(f"documents     : {sum(v.get('docs', 0) for v in man['shards'].values()):,}")
        print(f"sentences kept: {sum(v.get('kept_sentences', 0) for v in man['shards'].values()):,}")
        pf = out_dir / "pairs.tsv.gz"
        if pf.exists():
            print(f"pair file     : {pf} ({pf.stat().st_size/1e6:.0f} MB)")
        return

    print("loading the entity alias map ...", flush=True)
    idx = load_index(root)
    alias, amb = build_alias_map(idx)
    print(f"  {len(idx['alias']):,} aliases -> {len(alias):,} usable "
          f"({len(alias)/len(idx['alias']):.1%} kept after disambiguation filtering)")
    print(f"  sense collisions: {amb['redirected']} redirected to a measured "
          f"cancer-domain sense, {amb['dropped_ambiguous']} dropped as unresolvable")
    print(f"  thinly attested (< {MIN_ALIAS_SUPPORT} census mentions): "
          f"{amb['dropped_thin']:,} dropped")
    print(f"  minority form (< {100*MIN_ALIAS_SHARE:.0f}% of its identifier's "
          f"mentions): {amb['dropped_minority']:,} dropped")

    shards = sorted((fulltext_root() / "shards").glob("*.jsonl.gz"))
    if args.recent_first:
        shards = list(reversed(shards))
    if args.rebuild:
        # Both, or neither. Clearing the manifest alone re-counts every shard on
        # top of the existing table; clearing the table alone loses the shards
        # that will not be reprocessed.
        pf_old = out_dir / "pairs.tsv.gz"
        if pf_old.exists():
            pf_old.unlink()
        man = {"shards": {}}
        man_path.write_text(json.dumps(man, indent=1, sort_keys=True), encoding="utf-8")
        print("  --rebuild: manifest and pair table cleared", flush=True)

    todo = [s for s in shards if not man["shards"].get(s.name, {}).get("done")]
    if args.limit:
        todo = todo[:args.limit]
    print(f"shards: {len(shards)} total, {len(todo)} this run", flush=True)

    pairs: collections.Counter = collections.Counter()
    audit = {"rows": [], "seen": 0}
    arng = random.Random(_AUDIT_SEED)
    for i, s in enumerate(todo, 1):
        t0 = time.time()
        stats = process_shard(s, alias, pairs, audit=audit, rng=arng)
        stats["done"] = True
        stats["seconds"] = round(time.time() - t0, 1)
        man["shards"][s.name] = stats
        man_path.write_text(json.dumps(man, indent=1, sort_keys=True), encoding="utf-8")
        print(f"  [{i}/{len(todo)}] {s.name}: {stats['docs']:,} docs, "
              f"{stats['kept_sentences']:,} usable sentences, {len(pairs):,} pairs so far, "
              f"{stats['seconds']}s", flush=True)

    # Merge with any previous run. Safe only because a shard already recorded in
    # the manifest is skipped above, so its pairs are never added twice --
    # --rebuild clears both together for exactly this reason.
    pf = out_dir / "pairs.tsv.gz"
    if pf.exists():
        print("merging with the previous pair table ...", flush=True)
        with gzip.open(pf, "rt", encoding="utf-8") as fh:
            for line in fh:
                a, b, c = line.rstrip("\n").split("\t")
                pairs[(a, b)] += int(c)
    with gzip.open(pf, "wt", encoding="utf-8") as fh:
        for (a, b), c in pairs.most_common():
            fh.write(f"{a}\t{b}\t{c}\n")
    print(f"\nwrote {pf}: {len(pairs):,} distinct co-mentioned pairs")

    # The audit sample, written with the entity NAMES resolved so a reader can
    # judge a match without querying the index.
    if audit["rows"]:
        ap_path = out_dir / "audit-sample.jsonl.gz"
        with gzip.open(ap_path, "wt", encoding="utf-8") as fh:
            for row in audit["rows"]:
                row = dict(row)
                row["entity_names"] = [idx["canon"].get(e, e) for e in row["entities"]]
                fh.write(json.dumps(row) + "\n")
        print(f"wrote {ap_path}: {len(audit['rows'])} sentences uniformly sampled "
              f"from {audit['seen']:,} kept, for precision auditing")


if __name__ == "__main__":
    main()
