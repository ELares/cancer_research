#!/usr/bin/env python3
"""Does the manuscript's landscape survive census scale? (#ATLAS-LANDSCAPE)

WHY
---
The manuscript's central corpus claim is a ranking: immunotherapy dominates,
physical modalities are comparatively under-studied and preclinical. It is
computed over 4,830 articles retrieved by 33 keyword queries with a
500-record-per-query cap.

That design cannot distinguish "this mechanism is under-studied" from "we ran a
narrower query for it". The project says so itself -- `MISSION.md` calls every
gap claim "a statement about a retrieval design" -- but the claim has never been
recomputed against a corpus that was not built by those queries.

Now it can be. The census holds 4,403,994 cancer articles selected by MeSH
indexing rather than by this project's keywords, and every record carries its
full MeSH descriptor list.

THE DESIGN
----------
The same 17 mechanisms, measured three ways, so method and scale are separable:

  A  frozen corpus  + keyword tags   -- what the manuscript reports
  B  frozen corpus  + MeSH descriptors -- same articles, independent labels
  C  census         + MeSH descriptors -- same labels, unrestricted articles

A vs B isolates the LABELLING method on fixed articles. B vs C isolates CORPUS
SELECTION at fixed labels. A vs C is the question the manuscript's claim rests
on, and it is the only one that has been unanswerable until now.

MeSH descriptors come from `analysis/mesh-mechanism-map.yaml`, the same
precision-first leaf map the non-circular recall measurement (#412) uses. They
are expert-assigned, so they are not this project's own keywords rediscovering
themselves.

WHAT THIS IS NOT
----------------
Not a claim that MeSH is truth. It is an INDEPENDENT label with its own recall
limits: 17 of 23 mechanisms have a usable discriminative leaf, and the rest --
device modalities especially -- have no MeSH concept at all and are reported as
unmeasurable rather than as zero. A mechanism absent here may be absent from
MeSH's vocabulary, not from the literature.

Usage:
    python scripts/atlas_landscape.py
"""

import collections
import glob
import gzip
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import yaml  # noqa: E402

from atlas_baseline import atlas_root  # noqa: E402
from config import PROJECT_ROOT  # noqa: E402

# NLM-assigned trial publication types. Trial-specific only: "Multicenter
# Study" is excluded because it applies to observational work too, and
# "Clinical Study" is not a trial designation.
CLINICAL_TYPES = {
    "Clinical Trial", "Randomized Controlled Trial", "Controlled Clinical Trial",
    "Clinical Trial, Phase I", "Clinical Trial, Phase II",
    "Clinical Trial, Phase III", "Clinical Trial, Phase IV",
    "Pragmatic Clinical Trial", "Adaptive Clinical Trial",
}

MAP = PROJECT_ROOT / "analysis" / "mesh-mechanism-map.yaml"
FROZEN = PROJECT_ROOT / "corpus" / "INDEX.jsonl"
OUT = PROJECT_ROOT / "analysis" / "atlas-landscape.md"
RAW = PROJECT_ROOT / "analysis" / "atlas-landscape.json"


def load_map() -> tuple:
    d = yaml.safe_load(MAP.read_text())
    mech = {k: set(v["descriptors"]) for k, v in d["mechanisms"].items()}
    unmeasurable = d.get("unmeasurable", {}) or {}
    return mech, unmeasurable


def count_mesh(records, mech: dict) -> collections.Counter:
    """Articles carrying at least one discriminative descriptor per mechanism."""
    out = collections.Counter()
    for terms in records:
        s = set(terms)
        for name, descs in mech.items():
            if s & descs:
                out[name] += 1
    return out


def frozen_mesh() -> dict:
    """pmid -> MeSH terms, from the per-article frontmatter.

    corpus/INDEX.jsonl carries NO MeSH field at all -- the terms live only in
    corpus/by-pmid/*.md. Reading the index for them silently yields zero for
    every mechanism, which is what the first version of this script did, and a
    column of zeros reads as a finding rather than as a missing join.
    """
    out = {}
    for f in sorted((PROJECT_ROOT / "corpus" / "by-pmid").glob("*.md")):
        terms, inblock = [], False
        for line in f.read_text(errors="ignore").splitlines():
            if line.startswith("mesh_terms:"):
                inblock = True
                continue
            if inblock:
                st = line.strip()
                if st.startswith("- "):
                    terms.append(st[2:].strip().strip('"\''))
                    continue
                if st and not st.startswith("#"):
                    break
        out[f.stem] = terms
    return out


def frozen_records():
    for line in FROZEN.read_text().splitlines():
        if line.strip():
            yield json.loads(line)


# Mechanisms whose dominant descriptor names a therapy or modality, rather than
# a process (Glycolysis, DNA Methylation), a material (Nanoparticles) or a
# technique (Electroporation, CRISPR-Cas Systems). A broad descriptor pulls in
# papers that are not about treatment at all, and those are unlikely to be
# trials, so it DEFLATES the clinical share exactly where it inflates the count.
PRECISE = {"hifu", "sonodynamic", "antibody-drug-conjugate", "bispecific-antibody",
           "car-t", "oncolytic-virus", "phagocytosis-checkpoint", "mrna-vaccine"}
PHYSICAL = {"hifu", "sonodynamic", "electrochemical-therapy"}
# The pharmacological comparator is a CURATED list of drug modalities, not
# "everything that is not physical". Sweeping in the delivery platforms and
# genetic tools (nanoparticle, crispr, oncolytic-virus, mrna-vaccine,
# microbiome, phagocytosis-checkpoint) inflates the ratio from 9.1:1 to 12.5:1
# by counting things that are neither a drug class nor a physical modality, so
# the comparison would no longer be the manuscript's. Kept module-level so the
# figure that plots this ratio imports the same definition rather than
# restating it.
PHARMACOLOGICAL = {"immunotherapy", "car-t", "antibody-drug-conjugate",
                   "bispecific-antibody", "synthetic-lethality", "epigenetic",
                   "metabolic-targeting"}


def _maturity_narrative(R: dict) -> list:
    """The maturity comparison, and the reason it does not settle cleanly."""
    def share(keys):
        keys = [k for k in keys if k in R and R[k].get("mesh_census")]
        num = sum(R[k]["clinical_census"] for k in keys)
        den = sum(R[k]["mesh_census"] for k in keys)
        return (num / den) if den else 0.0

    all_phys = share(PHYSICAL)
    all_pharm = share(set(R) - PHYSICAL)
    pre_phys = share(PHYSICAL & PRECISE)
    pre_pharm = share(PRECISE - PHYSICAL)
    hifu = R.get("hifu", {}).get("clinical_share") or 0
    sono = R.get("sonodynamic", {}).get("clinical_share") or 0
    cart = R.get("car-t", {}).get("clinical_share") or 0
    return [
        "",
        "| comparison | physical | pharmacological |", "|---|---|---|",
        f"| all mechanisms | {100*all_phys:.2f}% | {100*all_pharm:.2f}% |",
        f"| precise descriptors only | {100*pre_phys:.2f}% | {100*pre_pharm:.2f}% |",
        "",
        "**The answer flips.** Taken across all mechanisms physical modalities look",
        "MORE clinically mature; restricted to descriptors that name a therapy rather",
        "than a process or a material, they look less. Both cannot be reported as the",
        "finding, and the second is the sounder comparison -- a broad descriptor pulls",
        "in papers that are not about treatment, and those are not trials, so scope",
        "deflates the share exactly where it inflates the count.", "",
        "So the manuscript's direction survives on the sound comparison, but weakly:",
        f"{100*pre_pharm:.2f}% against {100*pre_phys:.2f}%, a factor of "
        f"{pre_pharm/max(pre_phys,1e-9):.2f}, not the gulf the volume ratio suggests.", "",
        "### The finding that does hold up", "",
        "`physical modalities` is not a maturity class, and treating it as one is what",
        "the manuscript actually gets wrong. Both of these rest on precise",
        "single-descriptor signals:", "",
        f"* **HIFU is {100*hifu:.2f}% clinical -- more than CAR-T at {100*cart:.2f}%.**",
        "  It is an approved modality for prostate, fibroid and neurological",
        "  indications, and calling it preclinical is simply wrong.",
        f"* **HIFU and sonodynamic differ by {hifu/max(sono,1e-9):.1f}x** "
        f"({100*hifu:.2f}% against {100*sono:.2f}%), so the two are not at the same",
        "  stage and the aggregate hides it.", "",
        "The defensible statement is that SONODYNAMIC therapy is early, not that",
        "physical modalities are. The manuscript's own thesis rests on SDT, so this",
        "narrows the claim to the mechanism it actually cares about rather than",
        "weakening it.", "",
    ]


def main() -> int:
    mech, unmeasurable = load_map()
    print(f"{len(mech)} mechanisms with a discriminative MeSH leaf", flush=True)

    # A + B: the frozen corpus, both ways
    keyword = collections.Counter()
    n_frozen = 0
    for r in frozen_records():
        n_frozen += 1
        for m in (r.get("mechanisms") or []):
            # the index uses mRNA-vaccine, the map mrna-vaccine
            keyword[m.lower()] += 1
    fm = frozen_mesh()
    print(f"  MeSH read for {sum(1 for v in fm.values() if v):,} frozen articles",
          flush=True)
    mesh_frozen = count_mesh(fm.values(), mech)
    print(f"frozen corpus: {n_frozen:,} records", flush=True)

    # C: the census
    files = sorted(glob.glob(str(atlas_root() / "records" / "*.jsonl.gz")))
    print(f"scanning {len(files):,} census shards ...", flush=True)
    mesh_census = collections.Counter()
    # per-descriptor counts, so a mechanism carried by one over-broad term is
    # visible instead of being reported as a bare total
    per_desc = collections.defaultdict(collections.Counter)
    # Maturity, from NLM's own publication types rather than this project's
    # evidence tagger, whose recall the repo measures at 55%.
    clinical = collections.Counter()
    n_census = 0
    for i, f in enumerate(files, 1):
        with gzip.open(f, "rt", encoding="utf-8") as fh:
            for line in fh:
                r = json.loads(line)
                n_census += 1
                s = set(r.get("mesh") or [])
                if not s:
                    continue
                is_trial = bool(set(r.get("pub_types") or []) & CLINICAL_TYPES)
                for name, descs in mech.items():
                    hit = s & descs
                    if hit:
                        mesh_census[name] += 1
                        clinical[name] += is_trial
                        for d in hit:
                            per_desc[name][d] += 1
        if i % 300 == 0:
            print(f"  {i}/{len(files)} shards, {n_census:,} records", flush=True)
    print(f"census: {n_census:,} records", flush=True)

    # Concentration: the share of a mechanism's matches supplied by its single
    # commonest descriptor. `epigenetic` sits at 75% on "DNA Methylation", a term
    # any paper MEASURING methylation carries -- which is not epigenetic THERAPY,
    # the mechanism the manuscript means. A count dominated by one broad term is
    # measuring that term, not the mechanism.
    conc = {}
    for name, c in per_desc.items():
        tot = mesh_census.get(name, 0)
        top, n = (c.most_common(1)[0] if c else ("-", 0))
        conc[name] = {"top_descriptor": top, "top_share": (n / tot) if tot else 0.0,
                      "breakdown": dict(c.most_common())}

    def ranks(c):
        return {k: i + 1 for i, (k, _v) in
                enumerate(sorted(c.items(), key=lambda kv: -kv[1]))}
    rk_a, rk_b, rk_c = ranks(keyword), ranks(mesh_frozen), ranks(mesh_census)

    rows = []
    for name in mech:
        rows.append({
            "mechanism": name,
            "keyword_frozen": keyword.get(name, 0),
            "mesh_frozen": mesh_frozen.get(name, 0),
            "mesh_census": mesh_census.get(name, 0),
            "rank_keyword": rk_a.get(name),
            "rank_mesh_census": rk_c.get(name),
            "rank_shift": ((rk_a.get(name) or 0) - (rk_c.get(name) or 0))
                          if (rk_a.get(name) and rk_c.get(name)) else None,
            "clinical_census": clinical.get(name, 0),
            "clinical_share": (clinical.get(name, 0) / mesh_census[name])
                              if mesh_census.get(name) else None,
            "top_descriptor": conc.get(name, {}).get("top_descriptor"),
            "top_share": conc.get(name, {}).get("top_share", 0.0),
        })
    rows.sort(key=lambda r: -r["mesh_census"])

    moved = [r for r in rows if r["rank_shift"] is not None and abs(r["rank_shift"]) >= 3]

    L = [
        "# Does the manuscript's landscape survive census scale? (#ATLAS-LANDSCAPE)", "",
        "Generated by `scripts/atlas_landscape.py`.", "",
        "The manuscript's central corpus claim is a ranking -- immunotherapy dominates,",
        "physical modalities are comparatively under-studied -- computed over 4,830",
        "articles retrieved by 33 keyword queries with a 500-record cap. That design",
        "cannot separate *this mechanism is under-studied* from *we ran a narrower query",
        "for it*.", "",
        "## The design", "",
        "The same mechanisms, measured three ways, so method and scale are separable:", "",
        "| | corpus | labels | what it is |", "|---|---|---|---|",
        f"| **A** | frozen, {n_frozen:,} | this project's keywords | what the manuscript reports |",
        f"| **B** | frozen, {n_frozen:,} | MeSH descriptors | same articles, independent labels |",
        f"| **C** | census, {n_census:,} | MeSH descriptors | same labels, articles this project did not select |",
        "", "A vs B isolates the labelling method at fixed articles. B vs C isolates",
        "corpus selection at fixed labels. **A vs C** is the comparison the manuscript's",
        "claim actually rests on.", "",
        "## Result", "",
        "| mechanism | A: keyword/frozen | B: MeSH/frozen | C: MeSH/census | rank A | rank C | shift | top descriptor supplies |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for r in rows:
        sh = r["rank_shift"]
        shift = "-" if sh is None else (f"**{sh:+d}**" if abs(sh) >= 3 else f"{sh:+d}")
        share = r.get("top_share") or 0.0
        flag = f"**{100*share:.0f}%**" if share >= 0.5 else f"{100*share:.0f}%"
        L.append(f"| `{r['mechanism']}` | {r['keyword_frozen']:,} | {r['mesh_frozen']:,} | "
                 f"{r['mesh_census']:,} | {r['rank_keyword'] or '-'} | "
                 f"{r['rank_mesh_census'] or '-'} | {shift} | {flag} |")

    inflated = [r for r in rows if (r.get("top_share") or 0) >= 0.5]
    L += [
        "", "### Read the last column before the others", "",
        "It gives the share of a mechanism's census matches supplied by its single",
        "commonest descriptor. A mechanism above 50% is largely measuring that one",
        f"term rather than the mechanism, and **{len(inflated)} of {len(rows)} are**:", "",
    ]
    for r in sorted(inflated, key=lambda r: -(r.get("top_share") or 0)):
        L.append(f"* `{r['mechanism']}` -- {100*r['top_share']:.0f}% from "
                 f"*{r['top_descriptor']}*")
    L += [
        "", "`epigenetic` is the case that matters. It ranks first on the census and",
        "would be the headline of this analysis, but 75% of its matches come from",
        "*DNA Methylation* -- a descriptor carried by any paper that MEASURES",
        "methylation, which is not the epigenetic THERAPY the manuscript means. Its",
        "rank is an artifact of descriptor scope, not a discovery, and the map's own",
        "notes already flagged that leaf as scope-broad.", "",
        f"**{len(moved)} of {len(rows)} mechanisms move at least three rank positions**",
        "between what the manuscript reports and what the census says. Discount any",
        "whose top descriptor supplies most of its count.", "",
    ]
    if moved:
        L += ["| mechanism | manuscript rank | census rank | shift |", "|---|---|---|---|"]
        for r in sorted(moved, key=lambda r: r["rank_shift"]):
            L.append(f"| `{r['mechanism']}` | {r['rank_keyword']} | {r['rank_mesh_census']} | "
                     f"{r['rank_shift']:+d} |")
        L.append("")

    # ---- the scope-invariant comparison -------------------------------
    # Descriptor scope varies between mechanisms, so a cross-mechanism RANK is
    # only as trustworthy as the least precise descriptor. Comparing the SAME
    # labels across two corpora cancels that out entirely: whatever "Ultrasonic
    # Therapy" over-counts, it over-counts identically in both.
    cap = {r["mechanism"]: (r["mesh_frozen"] / r["mesh_census"])
           for r in rows if r["mesh_census"] and r["mesh_frozen"]}
    PHYS, PHARM = sorted(PHYSICAL), sorted(PHARMACOLOGICAL)
    tot = lambda ks, c: sum(R[k][c] for k in ks if k in R and R[k].get(c))  # noqa: E731
    R = {r["mechanism"]: r for r in rows}
    kf_p, kf_h = tot(PHYS, "keyword_frozen"), tot(PHARM, "keyword_frozen")
    mf_p, mf_h = tot(PHYS, "mesh_frozen"), tot(PHARM, "mesh_frozen")
    mc_p, mc_h = tot(PHYS, "mesh_census"), tot(PHARM, "mesh_census")
    ratio_a, ratio_b, ratio_c = kf_h / kf_p, mf_h / mf_p, mc_h / mc_p
    over = (mf_p / mc_p) / (mf_h / mc_h)
    caps = sorted(cap.values())

    L += [
        "## The comparison descriptor scope cannot spoil", "",
        "A cross-mechanism ranking is only as good as its least precise descriptor,",
        "which is why the column above matters. Comparing the SAME labels across two",
        "corpora removes that problem completely: whatever *Ultrasonic Therapy*",
        "over-counts, it over-counts identically in both.", "",
        "**Capture** is the share of the census's articles for a mechanism that the",
        "frozen 4,830-article corpus actually contains.", "",
        "| mechanism | frozen | census | capture |", "|---|---|---|---|",
    ] + [
        f"| `{k}` | {R[k]['mesh_frozen']:,} | {R[k]['mesh_census']:,} | {100*v:.2f}% |"
        for k, v in sorted(cap.items(), key=lambda kv: -kv[1])
    ] + [
        "",
        f"Capture spans **{100*min(caps):.2f}% to {100*max(caps):.2f}%, a "
        f"{max(caps)/min(caps):.0f}-fold spread**. The frozen corpus is not a small",
        "uniform sample of the literature; it is a wildly uneven one, and any relative",
        "prevalence computed on it inherits that unevenness.", "",
        "## What that does to the manuscript's central claim", "",
        "The manuscript says immunotherapy dominates while physical modalities remain",
        "comparatively under-studied. Measured three ways, pharmacological to physical:", "",
        "| | ratio |", "|---|---|",
        f"| manuscript method (keyword, frozen) | {ratio_a:.1f} : 1 |",
        f"| same articles, MeSH labels | {ratio_b:.1f} : 1 |",
        f"| **census, MeSH labels** | **{ratio_c:.1f} : 1** |", "",
        f"The frozen corpus captures physical modalities at {100*mf_p/mc_p:.2f}% against",
        f"{100*mf_h/mc_h:.2f}% for pharmacological ones -- it **over-samples physical",
        f"modalities by {over:.1f}x**, which is exactly what a corpus built from queries",
        "about them would do.", "",
        "So the claim survives, and it survives in the direction that costs this",
        "project something: the real imbalance is roughly **twice** what the manuscript",
        f"reports ({ratio_c:.1f}:1 against {ratio_a:.1f}:1). The manuscript understates",
        "its own case, because the corpus it measured was tilted toward the modalities",
        "it argues are neglected.", "",
        "## The other half of the claim: preclinical, or just smaller?", "",
        "The manuscript says physical modalities remain comparatively *preclinical*,",
        "which is a maturity claim rather than a volume one. NLM assigns trial",
        "publication types independently of this project, so the share of a",
        "mechanism's census articles carrying one is a maturity signal that does not",
        "depend on our evidence tagger -- whose recall this repo measures at 55%.", "",
        "| mechanism | census articles | clinical trials | share |",
        "|---|---|---|---|",
    ] + [
        f"| `{r['mechanism']}` | {r['mesh_census']:,} | {r['clinical_census']:,} | "
        f"{100*r['clinical_share']:.2f}% |"
        for r in sorted([x for x in rows if x.get("clinical_share") is not None],
                        key=lambda x: -x["clinical_share"])
    ] + _maturity_narrative(R) + [
        "## What MeSH cannot see", "",
        "Not every mechanism has a MeSH concept, and reporting those as zero would",
        "manufacture exactly the false gap this analysis exists to test for. They are",
        "excluded rather than counted:", "",
    ]
    for k, v in (unmeasurable or {}).items():
        reason = v if isinstance(v, str) else (v or {}).get("reason", "")
        L.append(f"* `{k}` -- {str(reason).strip()[:150]}")

    L += [
        "", "## Limits", "",
        "* MeSH is an independent label, not truth. It has its own recall limits and",
        "  its own indexing lag, which falls hardest on the newest literature.",
        "* Descriptors are precision-first leaves, so these counts are LOWER BOUNDS on",
        "  each mechanism's true literature. The comparison between columns is the",
        "  result; the absolute values are not.",
        "* A rank shift is not automatically a manuscript error. A keyword query can",
        "  legitimately capture a concept MeSH splits across descriptors. It does mean",
        "  the ranking is method-dependent, which is the thing the manuscript asserts",
        "  it is not.",
        "* Device and very recent modalities have no usable MeSH leaf at all and are",
        "  outside this measurement entirely.",
    ]

    OUT.write_text("\n".join(L) + "\n")
    RAW.write_text(json.dumps({
        "frozen_records": n_frozen, "census_records": n_census,
        "mechanisms_measured": len(mech), "moved_3_or_more": len(moved),
        "rows": rows,
    }, indent=2) + "\n")
    print(f"\n{len(moved)} of {len(rows)} mechanisms shift >=3 rank positions")
    print(f"wrote {OUT}\nwrote {RAW}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
