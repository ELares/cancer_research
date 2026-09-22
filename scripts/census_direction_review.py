#!/usr/bin/env python3
"""Prepare an unreviewed direction packet and a text-free provenance report.

The local packet retains the complete two-intersection union and blank review
worksheets. Only the readiness JSON and Markdown belong in the repository.
Offline verification checks that packet against the recorded code; it cannot
revalidate the source census without its shards or authenticate an adjudicator.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import shutil
import subprocess
import tempfile

from atlas_baseline import atlas_root
from census_adjudication import candidate_csv
import census_hypoxia_direction as hypoxia
from census_snapshot import read_snapshot, validate_snapshot
import census_thesis_direction as thesis


REPO = Path(__file__).resolve().parent.parent
OUT_JSON = REPO / "analysis/census-direction-review.json"
OUT_MD = REPO / "analysis/census-direction-review.md"
SCANNERS = {"hypoxia": hypoxia, "thesis": thesis}
CLASS_LABELS = {
    "hypoxia": ("protects", "sensitises", "both", "neither"),
    "thesis": ("exploit", "obstacle", "both", "neither"),
}
SOURCE_PATHS = (
    "scripts/census_direction_review.py", "scripts/census_snapshot.py",
    "scripts/census_hypoxia_direction.py", "scripts/census_thesis_direction.py",
    "scripts/census_adjudication.py", "scripts/census_input.py",
    "scripts/atlas_baseline.py",
)
PAYLOAD_PATHS = (
    "intersection-records.jsonl",
    "hypoxia/census-hypoxia-direction.json",
    "hypoxia/census-hypoxia-direction.md", "hypoxia/candidates.csv",
    "thesis/census-thesis-direction.json",
    "thesis/census-thesis-direction.md", "thesis/candidates.csv",
)
BUNDLE_PATHS = (*PAYLOAD_PATHS, "readiness.json", "readiness.md")


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _json_bytes(value) -> bytes:
    return (json.dumps(value, indent=1, ensure_ascii=True, allow_nan=False)
            + "\n").encode("utf-8")


def _source_hashes() -> dict:
    return {name: _sha256((REPO / name).read_bytes()) for name in SOURCE_PATHS}


def _git_context() -> dict:
    try:
        head = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=REPO, capture_output=True,
            text=True, check=True).stdout.strip()
        status = subprocess.run(
            ["git", "status", "--porcelain"], cwd=REPO, capture_output=True,
            text=True, check=True).stdout
        return {"git_head": head, "git_dirty": bool(status)}
    except (OSError, subprocess.CalledProcessError):
        return {"git_head": None, "git_dirty": None}


def _selected(record: dict) -> bool:
    mesh = record.get("mesh")
    if mesh is not None and (not isinstance(mesh, list)
                             or any(not isinstance(item, str) for item in mesh)):
        raise ValueError("record mesh must be a list of text descriptors or null")
    for field in ("title", "abstract"):
        value = record.get(field)
        if value is not None and not isinstance(value, str):
            raise ValueError(f"record {field} must be text or null")
    return any(scanner.in_scope(record) for scanner in SCANNERS.values())


def _analysis_files(records: list[dict]) -> tuple[dict, dict[str, bytes]]:
    summaries, files = {}, {}
    for name, scanner in SCANNERS.items():
        raw = scanner.scan_records(records)
        data = scanner.assemble(raw)
        files[f"{name}/census-{name}-direction.json"] = _json_bytes(data)
        files[f"{name}/census-{name}-direction.md"] = scanner.render(data).encode("utf-8")
        files[f"{name}/candidates.csv"] = candidate_csv(raw["cohort"]).encode("utf-8")
        missing = {label: {"title": 0, "abstract": 0, "both": 0}
                   for label in raw["counts"]}
        for record in records:
            if scanner.in_scope(record):
                label = scanner.candidate_label(record)
                title = not str(record.get("title") or "").strip()
                abstract = not str(record.get("abstract") or "").strip()
                missing[label]["title"] += int(title)
                missing[label]["abstract"] += int(abstract)
                missing[label]["both"] += int(title and abstract)
        summaries[name] = {
            "intersection_records": data["total"],
            "counts": data["counts"],
            "candidate_records": data["classified"],
            "excluded_records": data["unclassified"],
            "cohort_sha256": raw["cohort"]["cohort_sha256"],
            "missing_text": missing,
        }
    return summaries, files


def render(report: dict) -> str:
    snapshot, provenance = report["snapshot"], report["provenance"]
    lines = [
        "# Direction candidates prepared for independent review\n",
        "**Unreviewed: no biological decisions have been assigned.** The local "
        "packet contains complete blank worksheets for the singly classified "
        "hypoxia and thesis candidates. Historical reports and labels remain "
        "separate.\n",
        f"Captured {provenance['created_utc']} with Python "
        f"{provenance['python_version']}. The indexed `records/` snapshot has "
        f"{snapshot['observed_records']:,} parsed records across "
        f"{len(snapshot['shards']):,} shards; its manifest declares "
        f"{snapshot['declared_records']:,}. Every sorted shard was read "
        "(`stride=1`). This validates the acquired local snapshot, not its "
        "coverage of literature published since acquisition.\n",
        f"The local packet retains {report['union_records']:,} complete records "
        "in the union of the two descriptor intersections. This public report "
        "contains counts and fingerprints without article text.\n",
        "The indexed stream can include C04 and adjacent records; the direction "
        "scanners do not add a C04-only restriction. Snapshot cancer-basis "
        "counts: " + ", ".join(f"`{key}` {value:,}" for key, value in
                                 sorted(snapshot["cancer_basis"].items())) + ".\n",
        "## Candidate coverage\n",
        "| Analysis | Intersection | Singly classified candidates | Both/neither excluded |",
        "|---|---:|---:|---:|",
    ]
    for name in SCANNERS:
        summary = report["analyses"][name]
        lines.append(f"| {name} | {summary['intersection_records']:,} | "
                     f"{summary['candidate_records']:,} | {summary['excluded_records']:,} |")
    lines.extend([
        "", "All four lexical classes are retained in the intersection. Only "
        "singly classified candidates enter the exported worksheets; these "
        "counts do not measure biological direction or classifier accuracy. "
        "Even complete candidate adjudication cannot establish recall among "
        "excluded or unindexed articles.\n",
        "| Analysis | Lexical class | Records | Missing title | Missing abstract | Missing both |",
        "|---|---|---:|---:|---:|---:|",
    ])
    for name in SCANNERS:
        summary = report["analyses"][name]
        for label in CLASS_LABELS[name]:
            count = summary["counts"][label]
            missing = summary["missing_text"][label]
            lines.append(f"| {name} | {label} | {count:,} | {missing['title']:,} | "
                         f"{missing['abstract']:,} | {missing['both']:,} |")
    lines.extend([
        "", "Missing text includes empty and whitespace-only fields. Missing "
        "both is included in each separate missing-field count. Exported "
        "evidence must remain unchanged: record later source checks and any "
        "additional evidence in review reasons or a sidecar keyed by the "
        "record fingerprint.\n",
        "## Reproduction and review\n",
        "Generate a new local packet outside the repository using the same "
        "acquired census and code:\n",
        "```bash",
        'python scripts/census_direction_review.py --atlas-root "$FERRO_ATLAS_ROOT" \\\n  --bundle-dir /path/outside/repository/direction-review',
        "python scripts/census_direction_review.py \\\n  --verify-bundle /path/outside/repository/direction-review",
        "python scripts/census_direction_review.py --render-only",
        "```\n",
        "The packet directory must not already exist. All parsing, hashing, "
        "serialization and rendering finish before publication. The packet "
        "uses a sibling staging directory and rename; the subsequent public "
        "JSON and Markdown writes are not an atomic pair against disk failure "
        "or interruption.\n",
        "Offline verification checks the fixed packet paths and hashes, "
        "recorded source-code fingerprints, and all regenerated scanner "
        "reports, blank worksheets, cohort fingerprints and readiness counts "
        "against the retained intersection records. It does not reread the "
        "full census, prove that the retained union is complete, authenticate "
        "provenance, or assess scientific labels. Full input verification "
        "requires the original shards and a new build.\n",
        "Independent readers and source checks remain outstanding. Preserve "
        "separate reviewer decisions, identities, evidence sources and "
        "locators, and documented disagreement resolution. Only completed "
        "identity-matched candidate worksheets can be imported using "
        "`docs/CENSUS_DIRECTION_ADJUDICATION.md`; historical title-only labels "
        "cannot be carried forward.\n",
        "## Provenance\n",
        f"Source URL recorded by the snapshot: `{snapshot['source_url']}`.\n",
        f"Source inventory SHA-256: `{snapshot['inventory_sha256']}`.\n",
        f"Acquisition manifest SHA-256: `{snapshot['manifest_sha256']}`.\n",
        f"Baseline manifest SHA-256: `{snapshot['baseline_manifest_sha256']}`.\n",
        f"Git context: `{provenance['git_head']}`; working tree dirty: "
        f"`{provenance['git_dirty']}`. Source hashes below are authoritative "
        "for the code used; the Git head alone may predate uncommitted code.\n",
        "| Analysis | Cohort SHA-256 |", "|---|---|",
    ])
    for name in SCANNERS:
        summary = report["analyses"][name]
        lines.append(f"| {name} | `{summary['cohort_sha256']}` |")
    lines.extend(["", "| Source file | SHA-256 |", "|---|---|"])
    for name in SOURCE_PATHS:
        digest = provenance["source_sha256"][name]
        lines.append(f"| `{name}` | `{digest}` |")
    lines.extend(["", "| Local packet payload | SHA-256 |", "|---|---|"])
    for name in PAYLOAD_PATHS:
        digest = report["payload_sha256"][name]
        lines.append(f"| `{name}` | `{digest}` |")
    lines.append("")
    return "\n".join(lines)


def prepare_bundle(root: Path) -> tuple[dict, dict[str, bytes]]:
    """Read one snapshot and prepare every artifact without replacing outputs."""
    try:
        source_hashes = _source_hashes()
        snapshot, records = read_snapshot(Path(root), _selected)
        summaries, files = _analysis_files(records)
        files["intersection-records.jsonl"] = "".join(
            json.dumps(record, sort_keys=True, ensure_ascii=False, allow_nan=False)
            + "\n" for record in records).encode("utf-8")
        report = {
            "schema_version": 1, "status": "unreviewed", "snapshot": snapshot,
            "provenance": {
                "created_utc": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
                "python_version": platform.python_version(), **_git_context(),
                "source_sha256": source_hashes,
            },
            "union_records": len(records), "analyses": summaries,
            "payload_sha256": {name: _sha256(files[name]) for name in PAYLOAD_PATHS},
        }
        files["readiness.json"] = _json_bytes(report)
        files["readiness.md"] = render(report).encode("utf-8")
        if _source_hashes() != source_hashes:
            raise SystemExit("Review source code changed during preparation; no outputs were published.")
        return report, files
    except (OSError, ValueError, TypeError, KeyError, AttributeError) as exc:
        raise SystemExit(f"Cannot prepare direction review packet: {exc}") from exc


def verify_bundle(bundle: Path) -> dict:
    """Verify retained evidence offline, without claiming full-census coverage."""
    bundle = Path(bundle)
    try:
        entries = list(bundle.rglob("*"))
        if any(path.is_symlink() for path in entries):
            raise ValueError("packet entries must not be symbolic links")
        paths = {path.relative_to(bundle).as_posix() for path in entries if path.is_file()}
        if paths != set(BUNDLE_PATHS):
            raise ValueError("packet must contain exactly the fixed review payload and readiness files")
        files = {name: (bundle / name).read_bytes() for name in BUNDLE_PATHS}
        report = json.loads(files["readiness.json"])
        if report["schema_version"] != 1 or report["status"] != "unreviewed":
            raise ValueError("unsupported schema or packet is not marked unreviewed")
        if report["provenance"]["source_sha256"] != _source_hashes():
            raise ValueError("source code fingerprints differ from the recorded build")
        validate_snapshot(report["snapshot"])
        expected_hashes = {name: _sha256(files[name]) for name in PAYLOAD_PATHS}
        if report["payload_sha256"] != expected_hashes:
            raise ValueError("payload fingerprints do not match readiness.json")
        records = [json.loads(line) for line in files["intersection-records.jsonl"].decode("utf-8").splitlines()]
        if any(not isinstance(record, dict) or not _selected(record) for record in records):
            raise ValueError("intersection sidecar contains an invalid or out-of-scope record")
        summaries, regenerated = _analysis_files(records)
        if report["union_records"] != len(records) or report["analyses"] != summaries:
            raise ValueError("readiness counts or cohorts differ from the retained intersection records")
        for name, content in regenerated.items():
            if files[name] != content:
                raise ValueError(f"regenerated {name} differs from the packet")
        if files["readiness.md"] != render(report).encode("utf-8"):
            raise ValueError("readiness Markdown differs from its JSON report")
        return report
    except (OSError, ValueError, TypeError, KeyError, AttributeError) as exc:
        raise SystemExit(f"Cannot verify direction review packet: {exc}") from exc


def _exists(path: Path) -> bool:
    return path.exists() or path.is_symlink()


def publish_bundle(bundle: Path, files: dict[str, bytes]) -> None:
    """Publish a complete packet through a sibling staging directory rename."""
    bundle = Path(bundle)
    if _exists(bundle):
        raise SystemExit(f"Review packet directory already exists: {bundle}")
    if set(files) != set(BUNDLE_PATHS):
        raise SystemExit("Cannot publish an incomplete direction review packet.")
    bundle.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{bundle.name}-", dir=bundle.parent))
    try:
        for name, content in files.items():
            target = staging / name
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(content)
        if _exists(bundle):
            raise SystemExit(f"Review packet directory appeared during preparation: {bundle}")
        staging.rename(bundle)
    finally:
        if staging.exists():
            shutil.rmtree(staging)


def _check_bundle_path(bundle: Path, root: Path) -> None:
    target = bundle.resolve()
    for protected in (root.resolve(), OUT_JSON.resolve(), OUT_MD.resolve()):
        if target == protected or target in protected.parents or protected in target.parents:
            raise SystemExit("Review packet directory must not overlap census input or public report paths.")
    if _exists(bundle):
        raise SystemExit(f"Review packet directory already exists: {bundle}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--atlas-root", type=Path, default=atlas_root())
    parser.add_argument("--bundle-dir", type=Path)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--render-only", action="store_true")
    mode.add_argument("--verify-bundle", type=Path)
    args = parser.parse_args()
    if args.bundle_dir and (args.render_only or args.verify_bundle):
        parser.error("--bundle-dir cannot be combined with an offline mode")
    if args.verify_bundle:
        report = verify_bundle(args.verify_bundle)
        print(f"Verified unreviewed packet: {report['union_records']} retained intersection records; "
              "full census shards were not revalidated.")
        return 0
    if args.render_only:
        markdown = render(json.loads(OUT_JSON.read_text(encoding="utf-8")))
        OUT_MD.parent.mkdir(parents=True, exist_ok=True)
        OUT_MD.write_text(markdown, encoding="utf-8")
        print(f"Wrote {OUT_MD}")
        return 0
    if args.bundle_dir is None:
        parser.error("--bundle-dir is required to build a new local packet")
    _check_bundle_path(args.bundle_dir, args.atlas_root)
    report, files = prepare_bundle(args.atlas_root)
    json_text = json.dumps(report, indent=1, ensure_ascii=True, allow_nan=False) + "\n"
    md_text = files["readiness.md"].decode("utf-8")
    publish_bundle(args.bundle_dir, files)
    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_MD.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json_text, encoding="utf-8")
    OUT_MD.write_text(md_text, encoding="utf-8")
    print(f"Prepared unreviewed packet at {args.bundle_dir}; "
          f"{report['union_records']} retained intersection records.")
    print(f"Wrote {OUT_JSON} and {OUT_MD}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
