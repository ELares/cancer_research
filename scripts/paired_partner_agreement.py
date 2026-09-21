#!/usr/bin/env python3
"""Compare frozen keyword and MeSH partner counts on explicitly paired PMIDs.

Only committed index records, article frontmatter and the descriptor map are
read. No tagger, article body, download or raw census is used. Empty MeSH lists
are missing comparator information, never negative mechanism labels. Frozen
keyword labels can themselves use MeSH metadata: these arms are not independent.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
import os
import re
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path

import yaml

REPO = Path(__file__).resolve().parent.parent
INDEX = REPO / "corpus/INDEX.jsonl"
ARTICLES = REPO / "corpus/by-pmid"
MAP = REPO / "analysis/mesh-mechanism-map.yaml"
OUT_JSON = REPO / "analysis/paired-partner-agreement.json"
OUT_MD = REPO / "analysis/paired-partner-agreement.md"
TOP_K = 3
MIN_TARGET = 30
MIN_POSITIVE_PARTNERS = 3
SCHEMA_VERSION = 1


class InputError(ValueError):
    """Invalid or incomplete inputs must never replace existing reports."""


class _UniqueLoader(yaml.SafeLoader):
    pass


def _unique_mapping(loader, node, deep=False):
    out = {}
    for key_node, value_node in node.value:
        key = loader.construct_object(key_node, deep=deep)
        try:
            hash(key)
        except TypeError as exc:
            raise InputError("unhashable YAML key") from exc
        if key in out:
            raise InputError(f"duplicate YAML key: {key}")
        out[key] = loader.construct_object(value_node, deep=deep)
    return out


_UniqueLoader.add_constructor(
    yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, _unique_mapping)


def _json_object(pairs):
    out = {}
    for key, value in pairs:
        if key in out:
            raise InputError(f"duplicate JSON key: {key}")
        out[key] = value
    return out


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def pmid_hash(pmids) -> str:
    """Hash sorted string PMIDs, newline-delimited with a final newline."""
    return sha256(("\n".join(sorted(pmids)) + "\n").encode("utf-8"))


def _canonical_labels(value, context: str) -> list[str]:
    if not isinstance(value, list) or any(
            not isinstance(x, str) or not x.strip() for x in value):
        raise InputError(f"{context}: expected a list of nonempty strings")
    return sorted({x.strip().lower() for x in value})


def _pmid(value, context: str) -> str:
    if isinstance(value, bool) or not isinstance(value, (str, int)):
        raise InputError(f"{context}: missing or invalid PMID")
    result = str(value)
    if not re.fullmatch(r"0|[1-9][0-9]*", result):
        raise InputError(f"{context}: invalid PMID {result!r}")
    return result


def _frontmatter(path: Path) -> tuple[dict, str]:
    """Read and hash only the frontmatter, requiring its closing delimiter."""
    with path.open("rb") as fh:
        first = fh.readline()
        if first.strip() != b"---":
            raise InputError(f"{path.name}: missing frontmatter opening delimiter")
        lines = [first]
        for line in fh:
            lines.append(line)
            if line.strip() == b"---":
                break
        else:
            raise InputError(f"{path.name}: missing frontmatter closing delimiter")
    raw = b"".join(lines)
    data = yaml.load(b"".join(lines[1:-1]).decode("utf-8"), Loader=_UniqueLoader)
    if not isinstance(data, dict):
        raise InputError(f"{path.name}: frontmatter must be a mapping")
    return data, sha256(raw)


def _input_revision(paths: list[Path]) -> str | None:
    """Pin the last commit touching these inputs, unaffected by report commits."""
    try:
        relative = [str(p.resolve().relative_to(REPO)) for p in paths]
        result = subprocess.run(
            ["git", "log", "-1", "--format=%H", "--", *relative],
            cwd=REPO, capture_output=True, text=True, check=True)
        return result.stdout.strip() or None
    except (OSError, ValueError, subprocess.CalledProcessError):
        return None


def load_snapshot(index_path=INDEX, article_dir=ARTICLES, map_path=MAP) -> dict:
    """Validate the complete PMID join and retain operational labels only."""
    index_path, article_dir, map_path = map(Path, (index_path, article_dir, map_path))
    try:
        return _load_snapshot(index_path, article_dir, map_path)
    except (yaml.YAMLError, UnicodeError, json.JSONDecodeError) as exc:
        raise InputError(f"invalid metadata: {exc}") from exc


def _load_snapshot(index_path: Path, article_dir: Path, map_path: Path) -> dict:
    if not index_path.is_file() or not article_dir.is_dir() or not map_path.is_file():
        raise InputError("missing index, article directory or descriptor map")
    index_bytes, map_bytes = index_path.read_bytes(), map_path.read_bytes()
    mapping = yaml.load(map_bytes.decode("utf-8"), Loader=_UniqueLoader)
    if not isinstance(mapping, dict) or not isinstance(mapping.get("mechanisms"), dict):
        raise InputError("descriptor map requires a mechanisms mapping")
    descriptors, map_metadata, empty_groups = {}, {}, []
    canonical_keys = set()
    for name, item in mapping["mechanisms"].items():
        if not isinstance(name, str) or not name.strip() or not isinstance(item, dict):
            raise InputError("descriptor map has an invalid mechanism entry")
        canonical = name.strip().lower()
        if canonical in canonical_keys:
            raise InputError(f"duplicate canonical mechanism in map: {canonical}")
        canonical_keys.add(canonical)
        terms = _canonical_labels(item.get("descriptors"), f"map {canonical}.descriptors")
        if terms:
            descriptors[canonical] = terms
            map_metadata[canonical] = {
                "note": str(item.get("note", "")).strip(),
                "proxy_confounded": bool(item.get("proxy_confounded", False)),
            }
        else:
            empty_groups.append(canonical)
    if len(descriptors) < 2:
        raise InputError("descriptor map must contain at least two nonempty groups")

    index, frozen_vocabulary = {}, set()
    for lineno, line in enumerate(index_bytes.decode("utf-8").splitlines(), 1):
        if not line.strip():
            continue
        row = json.loads(line, object_pairs_hook=_json_object)
        if not isinstance(row, dict):
            raise InputError(f"index line {lineno}: expected an object")
        pmid = _pmid(row.get("pmid"), f"index line {lineno}")
        if pmid in index:
            raise InputError(f"duplicate index PMID: {pmid}")
        labels = _canonical_labels(row.get("mechanisms"), f"index {pmid}.mechanisms")
        index[pmid] = labels
        frozen_vocabulary.update(labels)
    if not index:
        raise InputError("index contains no records")
    files = sorted(article_dir.glob("*.md"))
    filenames = {_pmid(p.stem, p.name) for p in files}
    if len(filenames) != len(files):
        raise InputError("duplicate article filename PMID")
    if filenames != set(index):
        raise InputError(
            "index/article PMID sets differ: "
            f"{len(set(index) - filenames)} missing article files; "
            f"{len(filenames - set(index))} extra article files")

    universe = sorted(descriptors)
    universe_set = set(universe)
    articles, frontmatter_hashes, seen = [], [], set()
    for path in files:
        frontmatter, digest = _frontmatter(path)
        pmid = _pmid(frontmatter.get("pmid"), path.name)
        if pmid in seen:
            raise InputError(f"duplicate frontmatter PMID: {pmid}")
        seen.add(pmid)
        if pmid != path.stem or pmid not in index:
            raise InputError(f"filename/frontmatter/index PMID mismatch: {path.name} / {pmid}")
        frozen = _canonical_labels(frontmatter.get("mechanisms"), f"{path.name}.mechanisms")
        if frozen != index[pmid]:
            raise InputError(f"index/frontmatter mechanism mismatch: {pmid}")
        terms = _canonical_labels(frontmatter.get("mesh_terms"), f"{path.name}.mesh_terms")
        observed = bool(terms)
        term_set = set(terms)
        mesh = sorted(k for k, values in descriptors.items() if term_set.intersection(values))
        articles.append({
            "pmid": pmid,
            "keyword": frozen,
            "mesh": mesh if observed else None,
        })
        frontmatter_hashes.append({"pmid": pmid, "sha256": digest})
    if not any(r["mesh"] is not None for r in articles):
        raise InputError("no articles have nonempty MeSH; paired cohort is empty")
    hashes_payload = "".join(f"{r['pmid']}\t{r['sha256']}\n" for r in frontmatter_hashes)
    return {
        "schema_version": SCHEMA_VERSION,
        "vocabulary": universe,
        "descriptors": {k: descriptors[k] for k in universe},
        "descriptor_metadata": {k: map_metadata[k] for k in universe},
        "excluded_frozen_mechanisms": sorted(frozen_vocabulary - universe_set),
        "empty_descriptor_groups": sorted(empty_groups),
        "records": articles,
        "provenance": {
            "index_path": "corpus/INDEX.jsonl",
            "index_sha256": sha256(index_bytes),
            "descriptor_map_path": "analysis/mesh-mechanism-map.yaml",
            "descriptor_map_sha256": sha256(map_bytes),
            "input_last_changed_commit": _input_revision([index_path, article_dir, map_path]),
            "map_semantics_sha256": sha256(json.dumps(
                {k: descriptors[k] for k in universe}, sort_keys=True,
                separators=(",", ":")).encode("utf-8")),
            "all_pmid_sha256": pmid_hash(r["pmid"] for r in articles),
            "observed_mesh_pmid_sha256": pmid_hash(
                r["pmid"] for r in articles if r["mesh"] is not None),
            "missing_mesh_pmid_sha256": pmid_hash(
                r["pmid"] for r in articles if r["mesh"] is None),
            "article_frontmatter_sha256": frontmatter_hashes,
            "article_frontmatter_manifest_sha256": sha256(hashes_payload.encode("utf-8")),
            "frontmatter_manifest_hash_rule": "sorted PMID + TAB + frontmatter SHA-256 + newline",
            "pmid_set_hash_rule": "lexicographically sorted string PMIDs joined by newline, with final newline",
        },
    }


def average_ranks(values) -> list[float]:
    """Return one-based ascending ranks, averaging the occupied ranks in ties."""
    order = sorted(range(len(values)), key=lambda i: values[i])
    ranks = [0.0] * len(values)
    start = 0
    while start < len(order):
        end = start + 1
        while end < len(order) and values[order[end]] == values[order[start]]:
            end += 1
        rank = (start + 1 + end) / 2
        for i in order[start:end]:
            ranks[i] = rank
        start = end
    return ranks


def spearman(x, y) -> float | None:
    """Pearson correlation of average ranks; constant vectors are undefined."""
    if len(x) != len(y):
        raise ValueError("Spearman vectors must have equal lengths")
    if len(x) < 2:
        return None
    rx, ry = average_ranks(x), average_ranks(y)
    mean = (len(x) + 1) / 2
    covariance = math.fsum((a - mean) * (b - mean) for a, b in zip(rx, ry))
    vx = math.fsum((a - mean) ** 2 for a in rx)
    vy = math.fsum((b - mean) ** 2 for b in ry)
    return covariance / math.sqrt(vx * vy) if vx and vy else None


def top_k_membership(counts: dict[str, int], k=TOP_K) -> dict[str, float] | None:
    """Fractional top-k membership, shared equally across a boundary tie.

    Only positive counts can enter. A boundary tie splits the remaining slots
    evenly; zero-count partners never fill an under-supported top-k list.
    """
    if isinstance(k, bool) or not isinstance(k, int) or k < 1:
        raise ValueError("top-k must be a positive integer")
    if any(isinstance(v, bool) or not isinstance(v, int) or v < 0 for v in counts.values()):
        raise ValueError("partner counts must be nonnegative integers")
    positive = sorted((v for v in counts.values() if v > 0), reverse=True)
    if len(positive) < k:
        return None
    boundary = positive[k - 1]
    above = sum(v > boundary for v in positive)
    tied = positive.count(boundary)
    fraction = (k - above) / tied
    return {name: (1.0 if count > boundary else fraction if count == boundary else 0.0)
            for name, count in sorted(counts.items())}


def top_k_overlap(a: dict[str, int], b: dict[str, int], k=TOP_K) -> float | None:
    """Shared fractional top-k membership divided by k (range zero to one)."""
    if set(a) != set(b):
        raise ValueError("top-k comparisons require the same candidate partners")
    wa, wb = top_k_membership(a, k), top_k_membership(b, k)
    if wa is None or wb is None:
        return None
    # Python's built-in float summation changed in 3.12. Explicit fsum keeps
    # frozen artifacts reproducible on all supported Python versions.
    return math.fsum(min(wa[name], wb[name]) for name in wa) / k


def _validate_snapshot(snapshot: dict) -> tuple[list[str], list[dict]]:
    if not isinstance(snapshot, dict) or snapshot.get("schema_version") != SCHEMA_VERSION:
        raise InputError("unsupported or missing snapshot schema_version")
    vocabulary = snapshot.get("vocabulary")
    if (_canonical_labels(vocabulary, "snapshot vocabulary") != vocabulary
            or len(vocabulary) < 2):
        raise InputError("snapshot vocabulary must be sorted, unique and canonical")
    records = snapshot.get("records")
    if not isinstance(records, list) or not records:
        raise InputError("snapshot requires nonempty article records")
    universe, seen = set(vocabulary), set()
    for row in records:
        if not isinstance(row, dict):
            raise InputError("snapshot article must be an object")
        pmid = _pmid(row.get("pmid"), "snapshot article")
        if not isinstance(row["pmid"], str):
            raise InputError("snapshot PMID must be a string")
        if pmid in seen:
            raise InputError(f"duplicate snapshot PMID: {pmid}")
        seen.add(pmid)
        if _canonical_labels(row.get("keyword"), "snapshot keyword") != row["keyword"]:
            raise InputError("snapshot keyword labels must be sorted, unique and canonical")
        if "mesh" not in row:
            raise InputError("snapshot article is missing mesh observation state")
        if row["mesh"] is not None:
            if (_canonical_labels(row["mesh"], "snapshot mesh") != row["mesh"]
                    or not set(row["mesh"]).issubset(universe)):
                raise InputError("snapshot MeSH labels must be canonical vocabulary members")
    if not any(r["mesh"] is not None for r in records):
        raise InputError("snapshot has no MeSH-observed paired records")
    # Stored identity claims must remain attached to the operational labels
    # that actually regenerate the report. Never trust stale derived totals.
    provenance = snapshot.get("provenance", {})
    if not isinstance(provenance, dict):
        raise InputError("snapshot provenance must be a mapping")
    for key, pmids in (
            ("all_pmid_sha256", [r["pmid"] for r in records]),
            ("observed_mesh_pmid_sha256", [r["pmid"] for r in records if r["mesh"] is not None]),
            ("missing_mesh_pmid_sha256", [r["pmid"] for r in records if r["mesh"] is None])):
        if key in provenance and provenance[key] != pmid_hash(pmids):
            raise InputError(f"snapshot provenance mismatch: {key}")
    return vocabulary, sorted(records, key=lambda r: r["pmid"])


def _arm_summary(labels: list[set[str]]) -> dict:
    tagged = sum(bool(s) for s in labels)
    multi = sum(len(s) >= 2 for s in labels)
    return {
        "n": len(labels),
        "tagged": tagged,
        "untagged": len(labels) - tagged,
        "multi_tagged": multi,
        "tagged_denominator": len(labels),
        "multi_tagged_denominator": tagged,
        "multi_tagged_fraction": multi / tagged if tagged else None,
    }


def _profile(mechanism, candidates, records) -> dict:
    keyword, mesh = Counter(), Counter()
    nk = nm = both = 0
    for kw, ms in records:
        if mechanism in kw:
            nk += 1
            keyword.update(kw - {mechanism})
        if mechanism in ms:
            nm += 1
            mesh.update(ms - {mechanism})
        both += mechanism in kw and mechanism in ms
    a = {p: keyword[p] for p in candidates}
    b = {p: mesh[p] for p in candidates}
    positive = {"keyword": sum(v > 0 for v in a.values()),
                "mesh": sum(v > 0 for v in b.values())}
    reasons = []
    for arm, n in (("keyword", nk), ("mesh", nm)):
        if n < MIN_TARGET:
            reasons.append(f"{arm}: target count {n} < {MIN_TARGET}")
        if positive[arm] < MIN_POSITIVE_PARTNERS:
            reasons.append(f"{arm}: positive partners {positive[arm]} < {MIN_POSITIVE_PARTNERS}")
    rho = spearman(list(a.values()), list(b.values()))
    return {
        "n": len(records),
        "target_counts": {"keyword": nk, "mesh": nm, "both": both},
        "positive_partners": positive,
        "zero_partners": {arm: len(candidates) - n for arm, n in positive.items()},
        "partner_counts": {"keyword": a, "mesh": b},
        "top_k_membership": {"keyword": top_k_membership(a), "mesh": top_k_membership(b)},
        "scored": not reasons,
        "unscored_reasons": reasons,
        "spearman": rho if not reasons else None,
        "spearman_reason": ("support floor not met" if reasons else
                            "constant partner-count vector" if rho is None else None),
        "top_k_overlap": top_k_overlap(a, b) if not reasons else None,
    }


def assemble(snapshot: dict) -> dict:
    """Recompute every count, ranking and support decision from article labels."""
    vocabulary, records = _validate_snapshot(snapshot)
    universe = set(vocabulary)
    observed = [r for r in records if r["mesh"] is not None]
    missing = [r for r in records if r["mesh"] is None]
    paired = [(set(r["keyword"]) & universe, set(r["mesh"])) for r in observed]
    missing_keywords = [set(r["keyword"]) for r in missing]
    all_keywords = set().union(*(set(r["keyword"]) for r in records))
    pair_counts = {arm: Counter() for arm in ("keyword", "mesh")}
    for kw, ms in paired:
        pair_counts["keyword"].update(itertools.combinations(sorted(kw), 2))
        pair_counts["mesh"].update(itertools.combinations(sorted(ms), 2))
    rows = []
    for focal in vocabulary:
        candidates = [name for name in vocabulary if name != focal]
        row = {"mechanism": focal, **_profile(focal, candidates, paired)}
        shared = [(kw, ms) for kw, ms in paired if focal in kw and focal in ms]
        row["shared_focal"] = _profile(focal, candidates, shared)
        row["shared_focal"]["pmid_sha256"] = pmid_hash(
            r["pmid"] for r in observed if focal in r["keyword"] and focal in r["mesh"])
        rows.append(row)
    return {
        "schema_version": SCHEMA_VERSION,
        "snapshot": snapshot,
        "method": {
            "top_k": TOP_K,
            "minimum_target_count_each_arm": MIN_TARGET,
            "minimum_positive_partners_each_arm": MIN_POSITIVE_PARTNERS,
            "candidate_partners_per_focal": len(vocabulary) - 1,
            "spearman": "average ascending ranks over all candidate partners, including zeros",
            "top_k_ties": "fractional boundary membership; sum(min(keyword weight, MeSH weight))/k",
            "support_rule": "both target counts >= minimum and both positive-partner counts >= minimum; constant ranks undefined",
            "interpretation": "descriptive reporting floors, not precision guarantees or independent biological validation",
        },
        "cohorts": {
            "all": {"n": len(records), "pmid_sha256": pmid_hash(r["pmid"] for r in records)},
            "observed_mesh": {"n": len(observed), "pmid_sha256": pmid_hash(r["pmid"] for r in observed)},
            "missing_mesh": {"n": len(missing), "pmid_sha256": pmid_hash(r["pmid"] for r in missing)},
        },
        "primary": {
            "n": len(observed),
            "arms": {"keyword": _arm_summary([kw for kw, ms in paired]),
                     "mesh": _arm_summary([ms for kw, ms in paired])},
            "pairs": [{"mechanisms": [a, b], "keyword": pair_counts["keyword"][a, b],
                       "mesh": pair_counts["mesh"][a, b]}
                      for a, b in itertools.combinations(vocabulary, 2)],
            "rows": rows,
        },
        "missing_mesh": {
            "n": len(missing),
            "keyword": {"all": _arm_summary(missing_keywords),
                        "common": _arm_summary([s & universe for s in missing_keywords])},
            "mechanism_counts": {name: sum(name in s for s in missing_keywords)
                                 for name in sorted(all_keywords)},
        },
        "excluded_frozen_mechanisms": {
            name: {"all": sum(name in r["keyword"] for r in records),
                   "observed_mesh": sum(name in r["keyword"] for r in observed),
                   "missing_mesh": sum(name in r["keyword"] for r in missing)}
            for name in sorted(all_keywords - universe)},
    }


def _number(value) -> str:
    return f"{value:.3f}" if value is not None else "not scored"


def _coverage_text(summary: dict) -> str:
    fraction = summary["multi_tagged_fraction"]
    percentage = f"{100 * fraction:.2f}%" if fraction is not None else "undefined"
    return (f"{summary['tagged']:,}/{summary['tagged_denominator']:,} articles tagged; "
            f"{summary['multi_tagged']:,}/{summary['multi_tagged_denominator']:,} "
            f"tagged articles have at least two labels ({percentage})")


def _top_text(row: dict, arm: str) -> str:
    weights = row["top_k_membership"][arm]
    if weights is None:
        return "fewer than three positive partners"
    counts = row["partner_counts"][arm]
    return "; ".join(
        f"{name} ({counts[name]}; weight {weight:.3g})"
        for name, weight in sorted(weights.items(), key=lambda p: (-counts[p[0]], p[0]))
        if weight > 0)


def render(result: dict) -> str:
    """Render only derived results; no corpus or tagger reads occur here."""
    snapshot, primary = result["snapshot"], result["primary"]
    vocabulary, missing = snapshot["vocabulary"], result["missing_mesh"]
    provenance = snapshot.get("provenance", {})
    lines = [
        "# Paired frozen-keyword versus MeSH partner agreement",
        "",
        "Generated by `scripts/paired_partner_agreement.py` from committed frozen "
        "keyword labels and MeSH article frontmatter. No articles are retagged, "
        "and no download or raw census is required.",
        "",
        f"The explicit PMID join contains **{result['cohorts']['all']['n']:,} articles**. "
        f"The primary comparison retains all **{primary['n']:,} articles with nonempty "
        f"MeSH**, including articles with zero labels in either arm. The "
        f"**{missing['n']:,} articles with empty MeSH** lack comparator information "
        "and are excluded from both primary arms; they are not MeSH-negative records.",
        "",
        f"Both arms use the same {len(vocabulary)} nonempty descriptor-map groups. "
        f"For each focal mechanism, its {len(vocabulary) - 1} candidate partners "
        "are every other common group, including zero-count partners. Each unordered "
        "mechanism pair is counted once per article. Counts describe operational "
        "co-labeling, not tested treatment combinations.",
        "",
        "Frozen production keyword labels may use MeSH metadata as well as title, "
        "abstract and other annotations. These labeling methods are **not independent**. "
        "Agreement is neither independent biological validation nor evidence of "
        "census-wide ranking stability; it applies to this retrieved corpus and "
        "these operational definitions. Descriptor breadth and confounded proxies "
        "remain limitations, including Ultrasonic Therapy for sonodynamic therapy, "
        "DNA Methylation for epigenetic therapy, and mRNA Vaccines for cancer vaccination.",
        "",
        "## Methods fixed before inspecting the comparison",
        "",
        "Mechanism keys are lowercased (including mRNA-vaccine to mrna-vaccine). "
        "MeSH terms match mapped descriptor strings exactly, ignoring case. Multiple "
        "matching descriptors and duplicate tags count only once per mechanism/article. "
        "All nonempty map groups are retained even when their observed count is zero; "
        "empty or absent descriptor groups never become zero-valued candidate partners.",
        "",
        f"The reporting floors are at least **{MIN_TARGET} focal articles in each arm** "
        f"and **{MIN_POSITIVE_PARTNERS} positive-count partners in each arm**. These "
        "are descriptive reporting floors, not precision guarantees. Sparse rows retain "
        "their counts but neither agreement score is reported. A constant partner-count "
        "vector makes Spearman undefined even when support floors pass; it does not "
        "prevent a supported top-three overlap.",
        "",
        "Spearman correlates ascending average ranks of the complete count vectors, "
        "including tied zeros. Top-three membership uses positive counts only: "
        "partners above the third count receive weight 1, all partners tied at the "
        "boundary share the remaining slots equally, and all others receive 0. "
        "Weights sum to three. Overlap is the sum of the smaller arm-specific "
        "weight for each partner divided by three. It is fractional membership "
        "overlap, not expected overlap under random tie breaking. Partner names "
        "only order the display. No p-values, confidence intervals or bootstrap "
        "inference are reported.",
        "",
        "## Cohort coverage",
        "",
        f"- Primary frozen keywords, common groups: {_coverage_text(primary['arms']['keyword'])}.",
        f"- Primary MeSH, common groups: {_coverage_text(primary['arms']['mesh'])}.",
        f"- Missing-MeSH keywords, common groups: {_coverage_text(missing['keyword']['common'])}.",
        f"- Missing-MeSH keywords, all frozen groups: {_coverage_text(missing['keyword']['all'])}.",
        "",
        "The multi-label denominator is articles tagged at least once by that arm; "
        "these arm-specific denominators are coverage diagnostics, not the denominator "
        "of the paired rank comparison. Its article cohort remains the same in both arms.",
        "",
        "Excluded frozen keyword categories: "
        + ", ".join(f"`{name}`" for name in sorted(result["excluded_frozen_mechanisms"])) + ".",
        "",
        "## Primary partner agreement",
        "",
        "K = frozen keywords; M = MeSH. Positive and zero columns count candidate "
        "partners. Both counts the articles whose focal label is assigned by both methods. "
        "All rows are shown, with no pooled stable/unstable verdict.",
        "",
        "| Focal mechanism | Target K | Target M | Both | Positive K/M | Zero K/M | Spearman | Top-3 overlap |",
        "|---|--:|--:|--:|--:|--:|--:|--:|",
    ]
    for row in primary["rows"]:
        targets, positive, zero = row["target_counts"], row["positive_partners"], row["zero_partners"]
        lines.append(
            f"| {row['mechanism']} | {targets['keyword']:,} | {targets['mesh']:,} | "
            f"{targets['both']:,} | {positive['keyword']}/{positive['mesh']} | "
            f"{zero['keyword']}/{zero['mesh']} | {_number(row['spearman'])} | "
            f"{_number(row['top_k_overlap'])} |")
    for row in primary["rows"]:
        if row["unscored_reasons"]:
            lines.append(f"\n- `{row['mechanism']}`: " + "; ".join(row["unscored_reasons"]) + ".")
        elif row["spearman_reason"]:
            lines.append(f"\n- `{row['mechanism']}` Spearman: {row['spearman_reason']}.")
    lines += [
        "",
        "## Shared-focal sensitivity",
        "",
        "For each focal mechanism separately, restrict to articles where both methods "
        "assign that focal label, then recompute both partner vectors and their support. "
        "The same article membership is used for both arms within each row. This "
        "conditions on agreement about the target and asks a narrower question about "
        "partner labeling; it does not estimate the full-cohort comparison. N is also "
        "the focal target count in each arm.",
        "",
        "| Focal mechanism | Shared N | Positive K/M | Zero K/M | Spearman | Top-3 overlap |",
        "|---|--:|--:|--:|--:|--:|",
    ]
    for original in primary["rows"]:
        row = original["shared_focal"]
        positive, zero = row["positive_partners"], row["zero_partners"]
        lines.append(
            f"| {original['mechanism']} | {row['n']:,} | "
            f"{positive['keyword']}/{positive['mesh']} | "
            f"{zero['keyword']}/{zero['mesh']} | {_number(row['spearman'])} | "
            f"{_number(row['top_k_overlap'])} |")
    for original in primary["rows"]:
        row = original["shared_focal"]
        if row["unscored_reasons"]:
            lines.append(f"\n- `{original['mechanism']}`: " + "; ".join(row["unscored_reasons"]) + ".")
        elif row["spearman_reason"]:
            lines.append(f"\n- `{original['mechanism']}` Spearman: {row['spearman_reason']}.")
    lines += [
        "",
        "## Primary top-three membership",
        "",
        "Each entry gives its raw co-label count and fractional membership weight. "
        "Membership is descriptive even where a focal article support floor is not met; "
        "the agreement table controls which comparisons receive scores. The JSON "
        "stores every candidate partner count and the complete unordered pair table.",
        "",
        "| Focal mechanism | Frozen keyword partners (count; weight) | MeSH partners (count; weight) |",
        "|---|---|---|",
    ]
    for row in primary["rows"]:
        lines.append(f"| {row['mechanism']} | {_top_text(row, 'keyword')} | {_top_text(row, 'mesh')} |")
    lines += [
        "",
        "## Provenance and regeneration",
        "",
        "The JSON snapshot stores each PMID and its operational frozen keyword "
        "and mapped MeSH labels, with `null` for unobserved MeSH. It stores no article "
        "title, abstract or body. Frontmatter hashes pin the metadata actually read; "
        "the descriptor map's byte hash and normalized semantics hash pin the map.",
        "",
        "```sh",
        "python scripts/paired_partner_agreement.py",
        "python scripts/paired_partner_agreement.py --render-only",
        "```",
        "",
        "The first command validates the complete index/file/frontmatter PMID join "
        "and frozen-label agreement before writing these two new artifacts. The "
        "second recomputes every derived count and score from the stored snapshot "
        "without reading the index, article files, descriptor map or census. Inputs "
        "are never modified. Invalid or empty input leaves existing reports unchanged.",
        "",
        "| PMID cohort | N | SHA-256 |",
        "|---|--:|---|",
    ]
    for name in ("all", "observed_mesh", "missing_mesh"):
        cohort = result["cohorts"][name]
        lines.append(f"| {name} | {cohort['n']:,} | `{cohort['pmid_sha256']}` |")
    lines += ["", "PMID hashes sort string identifiers lexicographically, join them "
              "with newlines and append a final newline.", ""]
    for key in ("index_sha256", "descriptor_map_sha256", "map_semantics_sha256",
                "article_frontmatter_manifest_sha256", "input_last_changed_commit"):
        if key in provenance:
            lines.append(f"- {key}: `{provenance[key]}`")
    lines.append("")
    return "\n".join(lines)


def _write_reports(result: dict) -> None:
    """Finish validation and serialization before replacing either output."""
    payloads = [(OUT_JSON, json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n"),
                (OUT_MD, render(result))]
    staged = []
    try:
        for path, payload in payloads:
            with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=path.parent,
                                             prefix=f".{path.name}.", delete=False) as fh:
                staged.append((Path(fh.name), path))
                fh.write(payload)
        for temporary, destination in staged:
            os.replace(temporary, destination)
    finally:
        for temporary, destination in staged:
            temporary.unlink(missing_ok=True)


def _split_stored(stored: dict) -> tuple[dict]:
    """Recover the operational-label input for generic artifact regeneration."""
    return (stored["snapshot"],)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--render-only", action="store_true",
                        help="recompute from the stored operational-label snapshot")
    args = parser.parse_args(argv)
    try:
        if args.render_only:
            existing = json.loads(OUT_JSON.read_text(encoding="utf-8"), object_pairs_hook=_json_object)
            if not isinstance(existing, dict) or "snapshot" not in existing:
                raise InputError("existing report lacks its operational-label snapshot")
            result = assemble(*_split_stored(existing))
        else:
            result = assemble(load_snapshot(INDEX, ARTICLES, MAP))
    except (OSError, ValueError, yaml.YAMLError) as exc:
        print(f"Paired agreement failed: {exc}. Existing reports were not changed.", file=sys.stderr)
        return 2
    try:
        _write_reports(result)
    except (OSError, ValueError) as exc:
        print(f"Paired agreement publication failed: {exc}. Check both report files before use.",
              file=sys.stderr)
        return 2
    print(f"Wrote {OUT_JSON} and {OUT_MD}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
