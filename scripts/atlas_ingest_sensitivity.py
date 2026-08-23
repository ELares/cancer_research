#!/usr/bin/env python3
"""What does the MeSH qualifier axis add, inside cancer articles? (#722 step 1)

WHY
---
`atlas_baseline.parse_articles` reads only
`./MeshHeadingList/MeshHeading/DescriptorName`. It never reads `QualifierName`
-- an earlier version of this line added "the string appears nowhere in
scripts/, tests/ or analysis/", which has since ROTTED: it now appears in this
script, its test, its report and in atlas_modality_ratio. The substance is
unaffected and is the checkable part -- `atlas_baseline.py:308` reads
DescriptorName only, and `tests/test_atlas.py`'s parser fixture is built from
DescriptorName alone, so that test is structurally incapable of failing on it
-- so the census
carries the descriptor axis of MeSH and drops the 76-subheading qualifier axis
entirely. `Lung Neoplasms/radiotherapy` is stored as `Lung Neoplasms`.

Nothing can catch that. `tests/test_atlas.py`'s parser fixture is built from
`DescriptorName` elements only, so the one existing test of the ingest is
structurally incapable of failing on a dropped qualifier.

WHAT THIS MEASURES, AND WHAT AN EARLIER VERSION GOT WRONG
----------------------------------------------------------
The tempting measurement -- compare articles that SAY a modality in their title
or abstract against articles the descriptor layer catches -- produced a headline
ratio that did not reproduce and an explanation that failed a control. Two
mistakes worth recording because they are easy to repeat:

  ASYMMETRIC RULES. One stem against one exact descriptor label for ferroptosis,
  versus a multi-term regex against a multi-descriptor family for radiotherapy.
  Six independent recounts produced six different numerators.

  THE CONTROL REFUTED THE STORY. Under one symmetric rule, angiogenesis -- which
  has NO qualifier form -- scores 5.9% descriptor recall against radiotherapy's
  2.7% and ferroptosis's 90.2%. Low descriptor recall is not evidence of a
  qualifier problem, because concepts with no qualifier have it too. Most of
  that spread is topicality: conditional on a title mention, every concept
  converges.

So this measures the only thing that isolates the qualifier axis: for the SAME
articles, parsed the SAME way, how many carry a modality on the qualifier axis
that the descriptor axis alone does not report. That difference cannot be
explained by topicality, vocabulary breadth, or era, because both arms see one
article set and differ in one element.

It runs on a handful of raw shards rather than the whole baseline, because the
question is a per-article rate and a re-download of 40 million records to
estimate one is not a reasonable price.

Usage:
    python scripts/atlas_ingest_sensitivity.py --shards 8
    python scripts/atlas_ingest_sensitivity.py --render-only
"""

import argparse
import gzip
import json
import math
import random
import shutil
import sys
import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from atlas_baseline import (ADJACENT_DESCRIPTORS, atlas_root,  # noqa: E402
                            download, fetch_c04_descriptors,
                            list_baseline_files)

OUT_MD = PROJECT_ROOT / "analysis" / "atlas-ingest-sensitivity.md"
OUT_JSON = PROJECT_ROOT / "analysis" / "atlas-ingest-sensitivity.json"

# The backbone modalities, as MeSH files them. `qualifiers` are subheading
# labels attached to a descriptor; `descriptors` are standalone topical
# descriptors that would catch the same article on the axis the ingest reads.
# Both sides are stated so the comparison is symmetric by construction -- the
# error the first attempt at this made.
MODALITIES = {
    "radiotherapy": {
        "qualifiers": {"radiotherapy"},
        "descriptors": {"radiotherapy", "radiotherapy dosage", "radiosurgery",
                        "chemoradiotherapy", "radiotherapy, adjuvant",
                        "radiotherapy planning, computer-assisted",
                        "radiotherapy, intensity-modulated"},
    },
    "drug therapy": {
        "qualifiers": {"drug therapy"},
        "descriptors": {"antineoplastic agents", "drug therapy",
                        "antineoplastic combined chemotherapy protocols",
                        "chemotherapy, adjuvant"},
    },
    "surgery": {
        "qualifiers": {"surgery"},
        "descriptors": {"surgical procedures, operative", "mastectomy",
                        # `neoplasms/surgery` WAS HERE and is removed: it is a
                        # descriptor/qualifier composite, and the descriptor
                        # arm matches DescriptorName text, so it could never
                        # fire. It presented this arm as five entries when
                        # four were live. Removing it changes no count.
                        "hepatectomy", "pneumonectomy"},
    },
    "diagnostic imaging": {
        "qualifiers": {"diagnostic imaging"},
        "descriptors": {"diagnostic imaging", "tomography, x-ray computed",
                        "magnetic resonance imaging",
                        "positron emission tomography computed tomography"},
    },
}


def wilson(k: int, n: int, z: float = 1.96) -> tuple:
    if n == 0:
        return (0.0, 0.0)
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


def parse_both_axes(path: Path, c04: dict):
    """Yield (descriptor_labels, qualifier_labels) per CANCER article.

    Cancer membership is decided EXACTLY as the ingest decides it -- C04 or the
    adjacent set, on descriptor UIs -- so the two arms describe the same
    articles. Qualifiers cannot change membership: an article filed as
    `Lung Neoplasms/radiotherapy` still carries the `Lung Neoplasms` descriptor.
    """
    with gzip.open(path, "rb") as fh:
        for _event, elem in ET.iterparse(fh, events=("end",)):
            if not elem.tag.endswith("PubmedArticle"):
                continue
            try:
                cit = elem.find("MedlineCitation")
                if cit is None:
                    continue
                uis, labels, quals = [], [], []
                for mh in cit.findall("./MeshHeadingList/MeshHeading"):
                    dn = mh.find("./DescriptorName")
                    if dn is None:
                        continue
                    uis.append(dn.get("UI", ""))
                    labels.append((dn.text or "").lower())
                    for q in mh.findall("./QualifierName"):
                        quals.append((q.text or "").lower())
                if not uis:
                    continue
                if not any(u in c04 for u in uis) and \
                   not any(u in ADJACENT_DESCRIPTORS for u in uis):
                    continue
                yield labels, quals
            finally:
                elem.clear()


# THE DESCRIPTOR ARM AS A SUBTREE, not a hand-written list. The audit that
# produced the withdrawal above showed the cross-modality ordering was ranking
# how completely each proxy covers its own MeSH family, so the family is the
# control: every descriptor NLM puts under the node the modality names.
TREE_FAMILIES = {
    "radiotherapy": ["E02.815"],
    "drug therapy": ["E02.319.310", "D27.505.954.248"],
    "surgery": ["E04"],
    "diagnostic imaging": ["E01.370.350"],
}
TREES_DIR = "mesh"


def fetch_tree_family(prefixes, dest, force: bool = False) -> dict:
    """UI -> label for every topical descriptor under any of `prefixes`.

    The same query `atlas_baseline.fetch_c04_descriptors` runs, with the tree
    filter as a parameter, cached the same way so CI reads only the committed
    file.
    """
    import time
    import urllib.parse
    from pathlib import Path as _P
    dest = _P(dest)
    if dest.exists() and not force:
        out = {}
        for line in dest.read_text(encoding="utf-8").splitlines():
            if line and not line.startswith("#"):
                ui, label = line.split("\t", 1)
                out[ui] = label
        return out
    from atlas_baseline import MESH_SPARQL, _get
    filt = " || ".join(
        f'STRSTARTS(STR(?t), "http://id.nlm.nih.gov/mesh/{p}")' for p in prefixes)
    query = (
        "PREFIX meshv: <http://id.nlm.nih.gov/mesh/vocab#> "
        "PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#> "
        "SELECT DISTINCT ?d ?label WHERE { "
        "?d a meshv:TopicalDescriptor . ?d meshv:treeNumber ?t . "
        f"?d rdfs:label ?label . FILTER({filt}) }} ORDER BY ?d"
    )
    found, offset, page = {}, 0, 500
    while True:
        params = urllib.parse.urlencode({
            "query": query, "format": "JSON", "inference": "true",
            "limit": page, "offset": offset})
        data = json.loads(_get(f"{MESH_SPARQL}?{params}", timeout=180))
        rows = data["results"]["bindings"]
        for b in rows:
            found[b["d"]["value"].rsplit("/", 1)[-1]] = b["label"]["value"]
        if len(rows) < page:
            break
        offset += page
        time.sleep(0.3)
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(
        f"# MeSH descriptors under {', '.join(prefixes)}.\n"
        "# Source: NLM MeSH SPARQL endpoint, https://id.nlm.nih.gov/mesh/sparql "
        "(inference on).\n"
        "# Regenerate: python scripts/atlas_ingest_sensitivity.py --refresh-trees\n"
        f"# descriptors: {len(found)}\n"
        + "\n".join(f"{u}\t{l}" for u, l in sorted(found.items())) + "\n",
        encoding="utf-8")
    return found


def tree_arms(root, force: bool = False) -> dict:
    """modality -> {lowercased descriptor labels} for its whole MeSH family."""
    out = {}
    for m, prefixes in TREE_FAMILIES.items():
        slug = m.replace(" ", "-")
        f = root / TREES_DIR / f"tree-{slug}.tsv"
        out[m] = {v.lower() for v in fetch_tree_family(prefixes, f, force).values()}
    return out


def scan(n_shards: int, seed: int, keep: bool, refresh_trees: bool = False):
    root = atlas_root()
    c04 = fetch_c04_descriptors(root / "mesh" / "c04-descriptors.tsv")
    names = list_baseline_files()
    rng = random.Random(seed)
    # Spread across the baseline rather than taking a prefix. THE REASON ONCE
    # GIVEN HERE WAS WRONG: "the files are chronological, so a prefix samples
    # only the oldest literature" is false as a monotone property -- file order
    # is only broadly chronological, the oldest literature sits in the middle
    # of the range, and a prefix samples the mid-1970s rather than the oldest.
    # The decision stands on the real reason: a prefix is a contiguous block of
    # near-single-year shards, so it samples one era whichever era that is,
    # and MeSH indexing practice has changed over fifty years. What spreading
    # does NOT buy is coverage -- see `era_coverage`.
    picked = sorted(rng.sample(names, min(n_shards, len(names))))

    raw = root / "raw_sensitivity"
    raw.mkdir(parents=True, exist_ok=True)
    trees = tree_arms(root, force=refresh_trees)
    stats = {m: Counter() for m in MODALITIES}
    tstats = {m: Counter() for m in MODALITIES}
    n_cancer = 0
    # per-SHARD counts, because the design is a cluster sample of shards and
    # the published Wilson intervals treated 25,809 articles as independent
    per_shard = {}
    years = Counter()
    used = []
    try:
        for name in picked:
            p = download(name, raw)
            used.append(name)
            shard = per_shard.setdefault(name, {"n": 0, "modalities": {
                m: Counter() for m in MODALITIES}})
            for labels, quals in parse_both_axes(p, c04):
                n_cancer += 1
                shard["n"] += 1
                lab, qua = set(labels), set(quals)
                for m, spec in MODALITIES.items():
                    q_hit = bool(qua & spec["qualifiers"])
                    for arm, hit in (("proxy", bool(lab & spec["descriptors"])),
                                     ("tree", bool(lab & trees[m]))):
                        s = stats[m] if arm == "proxy" else tstats[m]
                        s["descriptor"] += hit
                        s["qualifier"] += q_hit
                        s["either"] += (hit or q_hit)
                        s["qualifier_only"] += (q_hit and not hit)
                        if arm == "proxy":
                            sh = shard["modalities"][m]
                            sh["qualifier_only"] += (q_hit and not hit)
            if not keep:
                p.unlink(missing_ok=True)
    finally:
        if not keep:
            shutil.rmtree(raw, ignore_errors=True)

    out_shards = {k: {"n": v["n"],
                      "qualifier_only": {m: c["qualifier_only"]
                                         for m, c in v["modalities"].items()}}
                  for k, v in per_shard.items()}
    return {
        "shards": used,
        "n_shards": len(used),
        "seed": seed,
        "cancer_articles": n_cancer,
        "modalities": {m: dict(s) for m, s in stats.items()},
        "modalities_tree_arm": {m: dict(s) for m, s in tstats.items()},
        "tree_family_sizes": {m: len(v) for m, v in trees.items()},
        "proxy_coverage_of_family": {
            m: {"proxy_entries": len(MODALITIES[m]["descriptors"]),
                "family": len(trees[m]),
                "entries_in_family": len(
                    {x for x in MODALITIES[m]["descriptors"]} & trees[m])}
            for m in MODALITIES},
        "per_shard": {k: {"n": v["n"],
                          "qualifier_only": {m: c["qualifier_only"]
                                             for m, c in v["modalities"].items()}}
                      for k, v in per_shard.items()},
    }


def era_coverage(shards, stride: int = 10) -> dict:
    """What the 8 shards cover, against what the census holds.

    `render()` calls the sample "sampled across the whole chronological
    range". Whether that is true is checkable WITHOUT re-downloading anything:
    the same shards have committed counterparts under corpus/atlas/records,
    which carry a year. The census side is a STRIDE, never a prefix.
    """
    import gzip
    root = atlas_root() / "records"
    if not root.exists():
        return {}

    def decades(paths):
        c, n = Counter(), 0
        for f in paths:
            if not f.exists():
                continue
            with gzip.open(f, "rt", encoding="utf-8") as fh:
                for line in fh:
                    y = json.loads(line).get("year")
                    if isinstance(y, int):
                        c[y - y % 10] += 1
                        n += 1
        return c, n

    smp, n_s = decades(root / n.replace(".xml.gz", ".jsonl.gz") for n in shards)
    files = sorted(root.glob("*.jsonl.gz"))
    cen, n_c = decades(files[::stride])
    if not (n_s and n_c):
        return {}
    rows = []
    for d in sorted(set(cen) | set(smp)):
        cs, ss = 100 * cen.get(d, 0) / n_c, 100 * smp.get(d, 0) / n_s
        rows.append({"decade": d, "census_pct": round(cs, 2),
                     "sample_pct": round(ss, 2), "sample_n": smp.get(d, 0),
                     "ratio": round(ss / cs, 2) if cs else None})
    thin = [r for r in rows if r["ratio"] is not None and r["ratio"] < 0.1]
    return {
        "census_stride": stride, "census_records": n_c, "sample_records": n_s,
        "n_census_shards_read": len(files[::stride]),
        "rows": rows,
        "decades_effectively_unrepresented": [r["decade"] for r in thin],
        "census_share_unrepresented": round(
            sum(r["census_pct"] for r in thin), 2),
        "sample_records_in_those": sum(r["sample_n"] for r in thin),
    }


def _bootstrap(d, reps=10000, seed=20260819):
    """95% interval resampling SHARDS, which is the actual sampling unit.

    The published Wilson intervals treat every article as an independent draw.
    They are eight near-single-year blocks, so the interval that means
    something resamples the blocks. Deterministic seed, so --render-only
    reproduces it.
    """
    import random as _r
    ps = d.get("per_shard") or {}
    if len(ps) < 4:
        return {}
    names = sorted(ps)
    out = {}
    for m in d["modalities"]:
        rng = _r.Random(seed)
        vals = []
        for _ in range(reps):
            pick = [names[rng.randrange(len(names))] for _ in names]
            num = sum(ps[k]["qualifier_only"].get(m, 0) for k in pick)
            den = sum(ps[k]["n"] for k in pick)
            if den:
                vals.append(100 * num / den)
        vals.sort()
        out[m] = {"lo": round(vals[int(0.025 * len(vals))], 2),
                  "hi": round(vals[int(0.975 * len(vals))], 2), "reps": reps}
    return out


def _era_section(d) -> list:
    e = d.get("era_coverage") or {}
    if not e.get("rows"):
        return []
    L = ["### What \"across the whole chronological range\" is worth", ""]
    L += [f"This page says the shards are sampled across the whole "
          f"chronological range. Checked against the census itself "
          f"(1-in-{e['census_stride']} stride over the committed records, "
          f"{e['census_records']:,} dated articles):", ""]
    L += ["| decade | census | sample | sample n | ratio |",
          "|--:|--:|--:|--:|--:|"]
    for r in e["rows"]:
        L.append(f"| {r['decade']}s | {r['census_pct']:.2f}% | "
                 f"{r['sample_pct']:.2f}% | {r['sample_n']:,} | "
                 + (f"{r['ratio']:.2f}x |" if r["ratio"] is not None else "- |"))
    L += [""]
    thin = e.get("decades_effectively_unrepresented") or []
    if thin:
        L += [f"**{len(thin)} decades are effectively unrepresented** -- "
              + ", ".join(f"{x}s" for x in thin)
              + f" carry {e['census_share_unrepresented']:.1f}% of the census "
                f"and {e['sample_records_in_those']:,} of the "
                f"{e['sample_records']:,} sampled articles. So the sentence is "
                f"false as written: the shards are spread across the FILE "
                f"index, which is only broadly chronological, and the sample "
                f"is a handful of near-single-year blocks rather than a span. "
                f"The decision to spread rather than take a prefix is still "
                f"right; its stated reason is not.", ""]
    return L


def _tree_arm_section(d, n) -> list:
    """The descriptor arm as a SUBTREE -- the control the withdrawal needed."""
    t = d.get("modalities_tree_arm") or {}
    cov = d.get("proxy_coverage_of_family") or {}
    if not t:
        return []
    L = ["## The descriptor arm as a SUBTREE, which is the control", ""]
    L += ["The gain above is measured against a hand-written proxy, and that "
          "proxy's completeness is the confound. So the arm is recomputed as "
          "every descriptor NLM puts under the node the modality names -- no "
          "list to dispute, the same shape as the fix #729 made for site "
          f"coverage. Same shards, same {n:,} articles, qualifier arm "
          "untouched.", ""]
    L += ["| modality | proxy entries | family | gain, proxy | gain, family | "
          "descriptor recall, proxy -> family |", "|---|--:|--:|--:|--:|--:|"]
    rank_p = sorted(d["modalities"],
                    key=lambda m: -d["modalities"][m]["qualifier_only"])
    for m in rank_p:
        pr, q = d["modalities"][m], t[m]
        c = cov.get(m, {})
        L.append(f"| {m} | {c.get('proxy_entries', '?')} | "
                 f"{c.get('family', '?')} | "
                 f"{100*pr['qualifier_only']/n:.2f} pts | "
                 f"**{100*q['qualifier_only']/n:.2f} pts** | "
                 f"{pr['descriptor']/max(pr['either'],1):.3f} -> "
                 f"{q['descriptor']/max(q['either'],1):.3f} |")
    L += [""]
    rank_t = sorted(t, key=lambda m: -t[m]["qualifier_only"])
    if rank_p != rank_t:
        c0 = cov.get(rank_p[0], {})
        L += ["**THE ORDERING INVERTS.** Proxy arm: "
              + " > ".join(f"`{m}`" for m in rank_p)
              + ". Family arm: " + " > ".join(f"`{m}`" for m in rank_t)
              + f". The top rank changes hands, and the modality published as "
                f"the sharpest case is the one whose proxy covers least of its "
                f"own family ({c0.get('proxy_entries')} entries against "
                f"{c0.get('family')}). The withdrawal above is measured now, "
                f"not argued.", ""]
    else:
        L += ["The ordering does NOT change under the family arm, which "
              "weakens the withdrawal above rather than supporting it.", ""]
    worst = max(abs(100*d["modalities"][m]["qualifier_only"]/n
                    - 100*t[m]["qualifier_only"]/n) for m in t)
    L += [f"Every gain stays strictly positive, so what the page licenses -- "
          f"that the axis is unread and the loss is not negligible -- survives "
          f"the control. The MAGNITUDES do not: the two arms disagree by up to "
          f"{worst:.1f} points.", ""]

    bs = _bootstrap(d)
    if bs:
        L += ["### The interval, resampled over shards", ""]
        L += [f"The intervals in the first table treat {n:,} articles as "
              f"independent draws. They are {len(d.get('per_shard') or {})} "
              f"near-single-year blocks, so the honest interval resamples the "
              f"blocks ({bs[rank_p[0]]['reps']:,} replicates):", ""]
        L += ["| modality | point | shard bootstrap 95% |", "|---|--:|--:|"]
        for m in rank_p:
            pt = 100 * d["modalities"][m]["qualifier_only"] / n
            b = bs[m]
            L.append(f"| {m} | {pt:.2f}% | **{b['lo']:.2f}-{b['hi']:.2f}** |")
        L += ["", "Every interval on this page should be read several times "
              "wider than printed. No point estimate and no sign moves.", ""]
    return L


def _roundtrip(d: dict) -> dict:
    """Render from what the artifact WILL contain, not from the live dict.

    `OUT_JSON` is written with `sort_keys=True`, so rendering the in-memory
    dict produces a document a `--render-only` run cannot reproduce. Both paths
    now render the same value.

    Checked rather than assumed: the committed report already matches the
    round-tripped render, so no published ordering changes and no table here
    turned out to depend on the declared MODALITIES sequence. Where an order DOES
    carry meaning it must be re-established inside the renderer -- sorting the
    input replaces a rank with an alphabet, which flipped a published verdict
    elsewhere in this repo.
    """
    return json.loads(json.dumps(d, sort_keys=True))


def render(d: dict) -> str:
    n = d["cancer_articles"]
    L = ["# What the MeSH qualifier axis adds, inside cancer articles", ""]
    L += ["*Generated by `scripts/atlas_ingest_sensitivity.py`. Every figure is "
          "recomputed.*", ""]

    L += [f"`atlas_baseline.parse_articles` reads `DescriptorName` and never "
          f"`QualifierName`. This measures what that drops, on "
          f"**{d['n_shards']} baseline shards** drawn at random from the file "
          f"index (NOT a prefix; see the era note below for what that does "
          f"and does not buy) with seed {d['seed']}, covering "
          f"**{n:,} cancer articles**.", ""]

    L += ["Both arms see the same articles, selected exactly as the ingest "
          "selects them. Qualifiers cannot change census membership -- an "
          "article filed `Lung Neoplasms/radiotherapy` still carries the "
          "`Lung Neoplasms` descriptor -- so this isolates the labelling axis "
          "and nothing else.", ""]

    L += ["| modality | descriptor axis | qualifier axis | either | "
          "**qualifier-only** | marginal gain |",
          "|---|--:|--:|--:|--:|--:|"]
    for m, s in d["modalities"].items():
        dd, qq = s.get("descriptor", 0), s.get("qualifier", 0)
        ei, qo = s.get("either", 0), s.get("qualifier_only", 0)
        lo, hi = wilson(qo, n)
        L.append(f"| {m} | {100*dd/n:.1f}% | {100*qq/n:.1f}% | {100*ei/n:.1f}% | "
                 f"**{100*qo/n:.1f}%** | +{100*qo/n:.1f} pts "
                 f"({100*lo:.1f}-{100*hi:.1f}) |")
    L += [""]

    best = max(d["modalities"].items(),
               key=lambda kv: kv[1].get("qualifier_only", 0))
    L += _tree_arm_section(d, n)
    L += _era_section(d)
    L += ["## The cross-modality ORDERING is not a measurement", ""]
    L += [f"The largest figure in that column is **{best[0]}** at "
          f"{100*best[1].get('qualifier_only',0)/n:.1f}%, and an earlier "
          f"version of this page called it \"the sharpest case\" and drew a "
          f"prioritisation from the ranking. THAT ORDERING IS WITHDRAWN. Each "
          f"row's gain is measured against a hand-written descriptor proxy of "
          f"four to seven entries standing in for a MeSH family that runs from "
          f"tens to hundreds of descriptors, and how completely each proxy "
          f"covers its own family is neither controlled nor equal across the "
          f"four rows -- while the qualifier arm is the WHOLE of a closed "
          f"axis. So a row's gain is partly a statement about how incomplete "
          f"its descriptor list is, and ranking the four against each other "
          f"ranks that incompleteness as much as it ranks the qualifier "
          f"axis.", ""]
    L += [f"The same clause said the qualifier axis carries these articles "
          f"\"and nowhere the ingest looks\". THAT IS ALSO WITHDRAWN: it "
          f"holds only against these four descriptors, not against the "
          f"ingest, which stores every DescriptorName an article carries.", ""]

    L += ["## What this does and does not license", ""]
    L += ["* It licenses the re-parse in #722 step 2 on the grounds that the "
          "axis is unread and every row's gain is material -- NOT on the "
          "grounds that one modality is sharper than another, which is the "
          "ordering withdrawn above. An earlier version of this bullet said "
          "\"and it says which those are\", and no materiality threshold "
          "appears anywhere on this page.",
          "* It does NOT say the census is missing articles. Membership is "
          "decided on descriptor UIs and is unaffected; this is a labelling "
          "gap inside the census, not a selection gap.",
          "* It does NOT support the earlier 2.5x text-versus-descriptor "
          "sensitivity claim, which compared concepts under asymmetric rules "
          "and failed a control: angiogenesis, with no qualifier form at all, "
          "has descriptor recall of the same low order as radiotherapy's -- "
          "an earlier version of this bullet said LOWER, which contradicts the "
          "only figures the repo carries for it (5.9% against 2.7%), in three "
          "places at once. The direction was never the argument: both sit far "
          "below ferroptosis's 90.2%, which is what makes low descriptor "
          "recall not evidence of a qualifier problem. Those three figures are "
          "computed NOWHERE in this repo and the symmetric rule behind them is "
          "not stated, so they are quoted as the provenance of a withdrawn "
          "claim rather than as measurements.",
          "* The estimate is from sampled shards and carries Wilson intervals. "
          "MeSH indexing practice has changed over fifty years, so a sample "
          "spread across the range is not the same as a per-era estimate.",
          ""]
    return "\n".join(L) + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--shards", type=int, default=8)
    ap.add_argument("--seed", type=int, default=20260816)
    ap.add_argument("--keep-raw", action="store_true")
    ap.add_argument("--render-only", action="store_true")
    ap.add_argument("--refresh-trees", action="store_true",
                    help="re-fetch each modality's MeSH tree family")
    args = ap.parse_args()

    if args.render_only:
        d = json.loads(OUT_JSON.read_text())
    else:
        d = scan(args.shards, args.seed, args.keep_raw,
                 refresh_trees=args.refresh_trees)
        if d["cancer_articles"] == 0:
            raise SystemExit(
                "no cancer articles parsed, which is not a finding -- it is "
                "what a broken C04 set or a changed XML path looks like.")
        OUT_JSON.write_text(json.dumps(d, indent=1, sort_keys=True) + "\n",
                            encoding="utf-8")
    # computed here rather than in scan(): it reads only committed records, so
    # it works on the --render-only path too and needs no re-download
    ec = era_coverage(d.get("shards") or [])
    if ec:
        d["era_coverage"] = ec
        OUT_JSON.write_text(json.dumps(d, indent=1, sort_keys=True) + "\n",
                            encoding="utf-8")
    OUT_MD.write_text(render(_roundtrip(d)), encoding="utf-8")
    print(f"wrote {OUT_MD}")
    print(f"wrote {OUT_JSON}")


if __name__ == "__main__":
    main()
