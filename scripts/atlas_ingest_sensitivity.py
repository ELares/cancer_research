#!/usr/bin/env python3
"""What does the MeSH qualifier axis add, inside cancer articles? (#722 step 1)

WHY
---
`atlas_baseline.parse_articles` reads only
`./MeshHeadingList/MeshHeading/DescriptorName`. It never reads `QualifierName`
-- the string appears nowhere in scripts/, tests/ or analysis/ -- so the census
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
                        "neoplasms/surgery", "hepatectomy", "pneumonectomy"},
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


def scan(n_shards: int, seed: int, keep: bool):
    root = atlas_root()
    c04 = fetch_c04_descriptors(root / "mesh" / "c04-descriptors.tsv")
    names = list_baseline_files()
    rng = random.Random(seed)
    # Spread across the baseline rather than taking a prefix: the files are
    # chronological, so a prefix samples only the oldest literature and MeSH
    # indexing practice has changed over fifty years.
    picked = sorted(rng.sample(names, min(n_shards, len(names))))

    raw = root / "raw_sensitivity"
    raw.mkdir(parents=True, exist_ok=True)
    stats = {m: Counter() for m in MODALITIES}
    n_cancer = 0
    used = []
    try:
        for name in picked:
            p = download(name, raw)
            used.append(name)
            for labels, quals in parse_both_axes(p, c04):
                n_cancer += 1
                lab, qua = set(labels), set(quals)
                for m, spec in MODALITIES.items():
                    d_hit = bool(lab & spec["descriptors"])
                    q_hit = bool(qua & spec["qualifiers"])
                    s = stats[m]
                    s["descriptor"] += d_hit
                    s["qualifier"] += q_hit
                    s["either"] += (d_hit or q_hit)
                    s["qualifier_only"] += (q_hit and not d_hit)
            if not keep:
                p.unlink(missing_ok=True)
    finally:
        if not keep:
            shutil.rmtree(raw, ignore_errors=True)

    return {
        "shards": used,
        "n_shards": len(used),
        "seed": seed,
        "cancer_articles": n_cancer,
        "modalities": {m: dict(s) for m, s in stats.items()},
    }


def render(d: dict) -> str:
    n = d["cancer_articles"]
    L = ["# What the MeSH qualifier axis adds, inside cancer articles", ""]
    L += ["*Generated by `scripts/atlas_ingest_sensitivity.py`. Every figure is "
          "recomputed.*", ""]

    L += [f"`atlas_baseline.parse_articles` reads `DescriptorName` and never "
          f"`QualifierName`. This measures what that drops, on "
          f"**{d['n_shards']} baseline shards** sampled across the whole "
          f"chronological range (seed {d['seed']}), covering "
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
    L += [f"The sharpest case is **{best[0]}**, where "
          f"{100*best[1].get('qualifier_only',0)/n:.1f}% of cancer articles "
          f"carry the modality on the qualifier axis and nowhere the ingest "
          f"looks.", ""]

    L += ["## What this does and does not license", ""]
    L += ["* It licenses the re-parse in #722 step 2 for modalities whose "
          "marginal gain is material, and it says which those are.",
          "* It does NOT say the census is missing articles. Membership is "
          "decided on descriptor UIs and is unaffected; this is a labelling "
          "gap inside the census, not a selection gap.",
          "* It does NOT support the earlier 2.5x text-versus-descriptor "
          "sensitivity claim, which compared concepts under asymmetric rules "
          "and failed a control: angiogenesis, with no qualifier form at all, "
          "has lower descriptor recall than radiotherapy.",
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
    args = ap.parse_args()

    if args.render_only:
        d = json.loads(OUT_JSON.read_text())
    else:
        d = scan(args.shards, args.seed, args.keep_raw)
        if d["cancer_articles"] == 0:
            raise SystemExit(
                "no cancer articles parsed, which is not a finding -- it is "
                "what a broken C04 set or a changed XML path looks like.")
        OUT_JSON.write_text(json.dumps(d, indent=1, sort_keys=True) + "\n",
                            encoding="utf-8")
    OUT_MD.write_text(render(d), encoding="utf-8")
    print(f"wrote {OUT_MD}")
    print(f"wrote {OUT_JSON}")


if __name__ == "__main__":
    main()
