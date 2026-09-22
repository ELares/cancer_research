#!/usr/bin/env python3
"""Audit extracted direction overlap and ambiguity-associated flag rates.

The diagnostic cohort uses raw entity identifiers and distinct
(pair, PMID, direction) incidences. It differs from atlas_contradictions.py,
which uses corrected identifiers and all-predicate relation totals. These
structural diagnostics do not measure extraction accuracy or biological
agreement. Same-paper overlap can reflect different contexts within a paper.

The committed aggregate JSON is a historical snapshot. --render-only derives
counts and point estimates from it without a relation dump or ambiguity scan,
leaves that JSON untouched, and identifies its unreproducible historical
bootstrap interval. A fresh analysis preserves the original cohort definition
and orders pairs before seeded resampling.

Usage:
    python scripts/atlas_contradiction_quality.py --render-only
    python scripts/atlas_contradiction_quality.py
"""

import argparse
import collections
import gzip
import json
import math
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from atlas_baseline import atlas_root  # noqa: E402
from atlas_contradictions import MIN_TOTAL, MIN_WEAK  # noqa: E402
from config import PROJECT_ROOT  # noqa: E402

SCAN = PROJECT_ROOT / "analysis" / "atlas-ambiguity.json"
OUT_MD = PROJECT_ROOT / "analysis" / "atlas-contradiction-quality.md"
OUT_JSON = PROJECT_ROOT / "analysis" / "atlas-contradiction-quality.json"
BOOTSTRAP = 2000
BOOT_SEED = 20260803
MAX_BUCKET = 10
DIRECTIONS = {"positive_correlate", "negative_correlate"}
RELATIONS = DIRECTIONS | {
    "associate", "treat", "cause", "inhibit", "stimulate", "cotreat",
    "interact", "compare", "prevent", "drug_interact",
}


def load_directional(root: Path):
    """Load raw-ID pairs -> distinct PMIDs per direction; reject invalid input.

    A nonempty, valid file with only nondirectional relations is a valid zero
    result. Missing, empty, malformed and undecodable input is not a zero result.
    """
    path = root / "relations" / "relations.tsv.gz"
    pos, neg = collections.defaultdict(set), collections.defaultdict(set)
    count = 0
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        for count, line in enumerate(fh, 1):
            p = line.rstrip("\r\n").split("\t")
            if (len(p) < 4 or not p[0].isascii() or not p[0].isdigit()
                    or int(p[0]) <= 0 or p[1] not in RELATIONS):
                raise ValueError(f"{path}:{count}: malformed relation row")
            ids = []
            for entity in p[2:4]:
                parts = entity.split("|")
                if (len(parts) != 2 or any(not part or any(c.isspace() for c in part)
                                           for part in parts)):
                    raise ValueError(f"{path}:{count}: expected Type|ID entity")
                ids.append(parts[1])
            if p[1] in DIRECTIONS:
                key = tuple(sorted(ids))
                (pos if p[1] == "positive_correlate" else neg)[key].add(p[0])
    if not count:
        raise ValueError(f"{path}: empty relation input")
    return pos, neg


def load_contested(scan: dict) -> set:
    """Validate every required ambiguity group before selecting its IDs."""
    if not isinstance(scan, dict) or not isinstance(scan.get("by_type"), dict):
        raise ValueError("ambiguity scan requires a by_type object")
    contested = set()
    for kind in ("gene", "chemical", "disease"):
        group = scan["by_type"].get(kind)
        if not isinstance(group, dict) or not isinstance(group.get("sense_rows"), list):
            raise ValueError(f"ambiguity scan requires {kind}.sense_rows list")
        for row in group["sense_rows"]:
            for side in ("top", "runner_up"):
                entry = row.get(side) if isinstance(row, dict) else None
                ident = entry.get("id") if isinstance(entry, dict) else None
                if not isinstance(ident, str) or not ident or ident != ident.strip():
                    raise ValueError(f"ambiguity scan requires a nonempty {kind}.{side}.id")
                contested.add(ident)
    return contested


def _strata(rows):
    strata = collections.defaultdict(lambda: {"amb": [0, 0], "clean": [0, 0]})
    for n, is_ambiguous, flagged in rows:
        bucket = min(n.bit_length() - 1, MAX_BUCKET)
        cell = strata[bucket]["amb" if is_ambiguous else "clean"]
        cell[0] += 1
        cell[1] += flagged
    return strata


def _ratio(an, ac, cn, cc):
    return (ac / an) / (cc / cn) if an and cn and cc else None


def _mh(strata):
    numerator = denominator = 0.0
    for bucket in sorted(strata, key=int):
        s = strata[bucket]
        an, ac = s["amb"]
        cn, cc = s["clean"]
        if an and cn:
            numerator += ac * cn / (an + cn)
            denominator += cc * an / (an + cn)
    return numerator / denominator if denominator else None


def analyze(pos, neg, contested, *, bootstrap=BOOTSTRAP, seed=BOOT_SEED) -> dict:
    """Compute this diagnostic only; no file access or biological adjudication.

    PMID sets deduplicate repeated rows. A PMID may contribute to many pairs,
    and to both directions of the same pair. Bootstrap sampling units are pairs.
    """
    if isinstance(bootstrap, bool) or not isinstance(bootstrap, int) or bootstrap < 0:
        raise ValueError("bootstrap must be a nonnegative integer")
    _count(seed, "bootstrap seed")
    rows = []
    conflicting_pairs = overlap_pairs = overlap = total = 0
    papers, overlap_papers = set(), set()
    for key in sorted(set(pos) | set(neg)):
        positive, negative = pos.get(key, set()), neg.get(key, set())
        n = len(positive) + len(negative)
        if n < MIN_TOTAL:
            continue
        flagged = min(len(positive), len(negative)) >= MIN_WEAK
        rows.append((n, any(ident in contested for ident in key), flagged))
        if flagged:
            both = positive & negative
            conflicting_pairs += 1
            overlap_pairs += bool(both)
            overlap += len(both)
            total += n
            papers.update(positive | negative)
            overlap_papers.update(both)
    strata = _strata(rows)
    an = sum(s["amb"][0] for s in strata.values())
    ac = sum(s["amb"][1] for s in strata.values())
    cn = sum(s["clean"][0] for s in strata.values())
    cc = sum(s["clean"][1] for s in strata.values())
    rng = random.Random(seed)
    boots = []
    for _ in range(bootstrap):
        sampled = [rows[rng.randrange(len(rows))] for _ in range(len(rows))]
        value = _mh(_strata(sampled))
        if value is not None and math.isfinite(value):
            boots.append(value)
    boots.sort()
    ci = ([boots[int(0.025 * len(boots))], boots[int(0.975 * len(boots))]]
          if boots else None)
    return {
        "conflicting_pairs": conflicting_pairs,
        "pairs_with_self_contradiction": overlap_pairs,
        "self_contradicting_assertions": overlap,
        "total_assertions_in_conflicts": total,
        "unique_papers_in_conflicts": len(papers),
        "unique_papers_with_overlap": len(overlap_papers),
        "eligible_pairs": len(rows),
        "ambiguous": {"n": an, "conflicted": ac},
        "clean": {"n": cn, "conflicted": cc},
        "crude_risk_ratio": _ratio(an, ac, cn, cc),
        "mantel_haenszel": _mh(strata),
        "mh_ci95": ci,
        "bootstrap": {
            "resamples": bootstrap,
            "seed": seed,
            "finite_resamples": len(boots),
            "undefined_resamples": bootstrap - len(boots),
            "sampling_unit": "entity pair",
            "row_order": "sorted raw identifier pairs",
        },
        "strata": {str(b): dict(v) for b, v in sorted(strata.items())},
    }


def _count(value, name):
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return value


def _validate(raw):
    """Check the retained counts before computing any report or writing output."""
    if not isinstance(raw, dict):
        raise ValueError("analysis must be a JSON object")
    for key in ("conflicting_pairs", "pairs_with_self_contradiction",
                "self_contradicting_assertions", "total_assertions_in_conflicts",
                "eligible_pairs"):
        _count(raw.get(key), key)
    if not isinstance(raw.get("strata"), dict):
        raise ValueError("strata must be an object")
    totals = {"amb": [0, 0], "clean": [0, 0]}
    for bucket, stratum in raw["strata"].items():
        if (not isinstance(bucket, str) or not bucket.isdigit()
                or str(int(bucket)) != bucket
                or not MIN_TOTAL.bit_length() - 1 <= int(bucket) <= MAX_BUCKET
                or not isinstance(stratum, dict)):
            raise ValueError(f"invalid assertion-count stratum: {bucket}")
        for group in totals:
            cell = stratum.get(group)
            if not isinstance(cell, list) or len(cell) != 2:
                raise ValueError(f"stratum {bucket}.{group} must be [pairs, flagged]")
            n, flagged = (_count(v, f"stratum {bucket}.{group}") for v in cell)
            if flagged > n:
                raise ValueError(f"stratum {bucket}.{group}: flagged exceeds pairs")
            totals[group][0] += n
            totals[group][1] += flagged
    for stored, group in (("ambiguous", "amb"), ("clean", "clean")):
        counts = raw.get(stored)
        if not isinstance(counts, dict):
            raise ValueError(f"missing {stored} counts")
        cell = [_count(counts.get(k), f"{stored}.{k}") for k in ("n", "conflicted")]
        if cell != totals[group]:
            raise ValueError(f"{stored} counts disagree with strata")
    if sum(t[0] for t in totals.values()) != raw["eligible_pairs"]:
        raise ValueError("eligible pair count disagrees with strata")
    if sum(t[1] for t in totals.values()) != raw["conflicting_pairs"]:
        raise ValueError("conflicting pair count disagrees with strata")
    overlap_pairs = raw["pairs_with_self_contradiction"]
    overlap = raw["self_contradicting_assertions"]
    total = raw["total_assertions_in_conflicts"]
    minimum_support = maximum_support = 0
    unbounded_support = False
    for bucket, stratum in raw["strata"].items():
        flagged = stratum["amb"][1] + stratum["clean"][1]
        lower = 2 ** int(bucket)
        minimum_support += flagged * lower
        maximum_support += flagged * (2 * lower - 1)
        unbounded_support |= int(bucket) == MAX_BUCKET and flagged > 0
    if total < minimum_support or (not unbounded_support and total > maximum_support):
        raise ValueError("directional incidence total is outside the retained stratum bounds")
    if (not overlap_pairs <= raw["conflicting_pairs"] or overlap < overlap_pairs
            or bool(overlap) != bool(overlap_pairs) or total < 2 * overlap
            or total < MIN_TOTAL * raw["conflicting_pairs"]
            or (not raw["conflicting_pairs"] and total)):
        raise ValueError("inconsistent within-pair overlap counts")
    unique_keys = ("unique_papers_in_conflicts", "unique_papers_with_overlap")
    if any(k in raw for k in unique_keys):
        papers, both_papers = (_count(raw.get(k), k) for k in unique_keys)
        if (papers > total - overlap or both_papers > min(papers, overlap)
                or bool(papers) != bool(total) or bool(both_papers) != bool(overlap)):
            raise ValueError("inconsistent unique-paper counts")
    ci = raw.get("mh_ci95")
    if ci is not None:
        if (not isinstance(ci, list) or len(ci) != 2
                or any(isinstance(v, bool) or not isinstance(v, (int, float))
                       or not math.isfinite(v) or v < 0 for v in ci)
                or ci[0] > ci[1]):
            raise ValueError("mh_ci95 must be null or two ordered finite values")
        if _mh(raw["strata"]) is None:
            raise ValueError("bootstrap interval exists for an undefined pooled ratio")
    meta = raw.get("bootstrap")
    if meta is not None:
        if not isinstance(meta, dict):
            raise ValueError("bootstrap metadata must be an object")
        for key in ("resamples", "seed", "finite_resamples", "undefined_resamples"):
            _count(meta.get(key), f"bootstrap.{key}")
        if meta["finite_resamples"] + meta["undefined_resamples"] != meta["resamples"]:
            raise ValueError("bootstrap resample counts do not add up")
        if bool(meta["finite_resamples"]) != (ci is not None):
            raise ValueError("bootstrap interval disagrees with finite resample count")
        if (meta.get("sampling_unit") != "entity pair"
                or meta.get("row_order") != "sorted raw identifier pairs"):
            raise ValueError("unsupported bootstrap sampling unit or row order")
    return totals


def _percent(numerator, denominator):
    if not denominator:
        return "not estimable (no observations)"
    value = 100 * numerator / denominator
    return f"{value:.6f}%" if 0 < value < 0.01 else f"{value:.2f}%"


def _formatted_ratio(value):
    return "not estimable" if value is None else f"{value:.2f}x"


def render(raw: dict) -> str:
    """Render a validated aggregate snapshot without reading or changing inputs."""
    totals = _validate(raw)
    an, ac = totals["amb"]
    cn, cc = totals["clean"]
    crude = _ratio(an, ac, cn, cc)
    pooled = _mh(raw["strata"])
    flagged = raw["conflicting_pairs"]
    overlap_pairs = raw["pairs_with_self_contradiction"]
    overlap = raw["self_contradicting_assertions"]
    incidences = raw["total_assertions_in_conflicts"]
    pair_papers = incidences - overlap
    lines = [
        "# Extracted direction overlap and ambiguity association (#ATLAS-CONTRA-Q)", "",
        "Generated by `scripts/atlas_contradiction_quality.py`. Reproduce this report",
        "with `python scripts/atlas_contradiction_quality.py --render-only`; that reads",
        "the committed JSON and leaves it unchanged. These are structural diagnostics,",
        "not measurements of extraction precision, recall, or biological disagreement.", "",
        "## Population and counting units", "",
        f"A pair is eligible with at least {MIN_TOTAL} distinct (pair, PMID, direction)",
        f"incidences. It is flagged when each direction has at least {MIN_WEAK} distinct PMIDs.",
        "Repeated rows for the same pair, PMID and direction count once. The same PMID",
        "can contribute to many pairs and to both directions of one pair.", "",
        "This diagnostic uses raw identifiers and only `positive_correlate` and",
        "`negative_correlate`. The contradiction queue instead uses corrected identifiers",
        "and totals over all relation predicates. Its counts and this diagnostic's counts",
        "therefore describe different cohorts; this is not an accuracy audit of that queue.", "",
        "## Same-paper overlap within flagged pairs", "",
        "| quantity | count | share |", "|---|---|---|",
        f"| flagged pairs | {flagged:,} | |",
        f"| flagged pairs with at least one PMID in both directions | {overlap_pairs:,} | "
        f"{_percent(overlap_pairs, flagged)} |",
        f"| distinct (pair, PMID, direction) incidences in flagged pairs | {incidences:,} | |",
        f"| distinct (pair, PMID) incidences in flagged pairs | {pair_papers:,} | |",
        f"| (pair, PMID) incidences appearing in both directions | {overlap:,} | "
        f"{_percent(overlap, pair_papers)} |", "",
        "The pair-PMID denominator is the directional-incidence count minus the overlap:",
        f"{incidences:,} - {overlap:,} = {pair_papers:,}. Each overlapping pair-PMID is counted",
        "twice in the directional total and once in the pair-PMID total.", "",
    ]
    if "unique_papers_in_conflicts" in raw:
        papers = raw["unique_papers_in_conflicts"]
        both_papers = raw["unique_papers_with_overlap"]
        lines += [
            f"Across flagged pairs there are {papers:,} unique PMIDs; {both_papers:,} have",
            "both directions for at least one flagged pair. These separate paper counts",
            "deduplicate PMIDs across pairs; the table above counts pair-PMID incidences.", "",
        ]
    else:
        lines += [
            "Unique-paper counts were not retained in this historical aggregate. The",
            "pair-PMID totals cannot supply a unique-paper denominator.", "",
        ]
    lines += [
        "Overlap can reflect distinct tissues, conditions, or claims within one paper;",
        "it is not itself proof of extraction error. Low overlap does not establish",
        "extraction accuracy, and opposite outputs from different papers do not establish",
        "scientific disagreement. Pairs below the flag thresholds were not tested here.", "",
        "## Association with measured identifier ambiguity", "",
        "A pair is marked as involving measured ambiguity when either raw identifier",
        "appears as the top or runner-up identifier in the ambiguity scan's sense rows.",
        "This does not show that the particular mentions forming that pair were conflated.", "",
        "| eligible pairs | pairs | flagged | flag rate |", "|---|---|---|---|",
        f"| involving a measured sense collision | {an:,} | {ac:,} | {_percent(ac, an)} |",
        f"| no measured sense collision | {cn:,} | {cc:,} | {_percent(cc, cn)} |", "",
        f"Crude risk ratio: **{_formatted_ratio(crude)}**.", "",
        "All stored assertion-count strata are shown, including sparse groups:", "",
        "| directional incidences | measured collision: flagged/pairs | "
        "no measured collision: flagged/pairs | risk ratio |", "|---|---|---|---|",
    ]
    ratios = []
    for bucket in sorted(raw["strata"], key=int):
        b = int(bucket)
        s = raw["strata"][bucket]
        sn, sc = s["amb"]
        tn, tc = s["clean"]
        ratio = _ratio(sn, sc, tn, tc)
        if ratio is not None:
            ratios.append(ratio)
        label = f"{2**b}+" if b == MAX_BUCKET else f"{2**b}-{2**(b+1)-1}"
        lines.append(f"| {label} | {sc}/{sn} ({_percent(sc, sn)}) | "
                     f"{tc}/{tn} ({_percent(tc, tn)}) | {_formatted_ratio(ratio)} |")
    lines += ["", f"Mantel-Haenszel pooled risk ratio: **{_formatted_ratio(pooled)}**.", ""]
    if any(a > b for a, b in zip(ratios, ratios[1:])):
        lines += ["The estimable stratum ratios do not increase monotonically with assertion count.", ""]
    lines += [
        "A risk ratio requires both groups and a nonzero comparison-group flag rate.",
        "Strata missing either group do not contribute to the pooled estimate; a zero",
        "pooled denominator makes that estimate undefined. Undefined estimates are",
        "reported as not estimable rather than as zero or an arbitrary large ratio.", "",
    ]
    ci = raw.get("mh_ci95")
    meta = raw.get("bootstrap")
    if meta is None:
        if ci is not None:
            lines += [f"Stored historical 95% pair-bootstrap interval: {ci[0]:.2f}-{ci[1]:.2f}."]
        else:
            lines += ["No historical bootstrap interval is available."]
        lines += [
            "The historical interval is retained, not recomputed: the original pair order",
            "was not saved, and the legacy sampler used unordered set iteration. A fixed",
            "random seed alone cannot reproduce its endpoints. Finite and undefined",
            "resample counts were not retained.", "",
        ]
    else:
        if ci is not None:
            lines += [f"95% pair-bootstrap percentile interval: {ci[0]:.2f}-{ci[1]:.2f}."]
        else:
            lines += ["Bootstrap interval: not estimable (no finite bootstrap estimates)."]
        lines += [
            f"Bootstrap resamples: {meta['resamples']:,}; finite: {meta['finite_resamples']:,};",
            f"undefined: {meta['undefined_resamples']:,}. Seed: {meta['seed']}. Fresh runs",
            "order raw identifier pairs before sampling. Percentiles use only finite",
            "resamples; discarding undefined estimates can affect their interpretation.", "",
        ]
    lines += [
        "The bootstrap samples pairs as if independent. Pairs can share PMIDs and",
        "entities, so the interval does not account for that dependence or for extraction",
        "and ambiguity-label error. It is not a validated uncertainty interval for",
        "biological disagreement or for a causal effect of ambiguity.", "",
        "## Interpretation and remaining evidence work", "",
        "The pooled ratio describes an association after grouping directional-incidence",
        "counts into broad bins. It does not establish that ambiguity caused extra flags",
        "or rule out popularity, context, or other confounding. The estimates cannot be",
        "used to subtract an attributed fraction of conflicts.", "",
        "No measured collision means only that this ambiguity scan did not flag the",
        "identifier. The scan covers selected surface forms; unmeasured ambiguity may",
        "remain, and the direction of resulting bias is not established here.", "",
        "Sentence-level claims, source context, independent adjudication, and an explicit",
        "sampling frame are still required to estimate extraction quality. Abstract versus",
        "full-text coverage and census-stream coverage cannot be recovered from these",
        "aggregate counts. The `treat` versus `cause` comparison is outside this analysis.",
    ]
    return "\n".join(lines) + "\n"


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--render-only", action="store_true",
                        help="render committed aggregate JSON without external inputs or rewriting it")
    args = parser.parse_args(argv)
    try:
        if args.render_only:
            raw = json.loads(OUT_JSON.read_text(encoding="utf-8"))
        else:
            contested = load_contested(json.loads(SCAN.read_text(encoding="utf-8")))
            print("reading directional relations ...", flush=True)
            pos, neg = load_directional(atlas_root())
            raw = analyze(pos, neg, contested)
        # Prepare and validate both products before replacing either output.
        outputs = [(OUT_MD, render(raw))]
        if not args.render_only:
            outputs.insert(0, (OUT_JSON, json.dumps(raw, indent=2, allow_nan=False) + "\n"))
        for path, payload in outputs:
            path.write_text(payload, encoding="utf-8")
    except (OSError, EOFError, UnicodeError, ValueError) as exc:
        print(f"atlas contradiction quality: {exc}", file=sys.stderr)
        return 1
    for path, _ in outputs:
        print(f"wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
