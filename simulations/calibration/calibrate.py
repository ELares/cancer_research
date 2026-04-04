#!/usr/bin/env python3
"""
Evaluate simulation parameters against published calibration targets.

Reads existing simulation output files (no recompilation needed) and
compares them to the target values defined in targets.yaml.

Usage:
    python calibrate.py --evaluate
"""

import argparse
import csv
import json
import sys
from pathlib import Path

import yaml

SCRIPT_DIR = Path(__file__).resolve().parent
SIM_ROOT = SCRIPT_DIR.parent
TARGETS_FILE = SCRIPT_DIR / "targets.yaml"
REPORT_FILE = SCRIPT_DIR / "calibration_report.md"


def load_targets() -> list[dict]:
    with open(TARGETS_FILE, encoding="utf-8") as f:
        data = yaml.safe_load(f)
    return data["targets"]


def _resolve_output_path(target: dict) -> Path:
    """Resolve the output file path from a target definition.

    sim-original outputs live in simulations/ root; all others live in
    simulations/output/.  The target's output_file is always relative to
    the appropriate base.
    """
    output_file = target["output_file"]
    candidate = SIM_ROOT / "output" / output_file
    if candidate.exists():
        return candidate
    return SIM_ROOT / output_file


def extract_sim_original(target: dict) -> float | None:
    """Extract observable from sim-original JSON output."""
    path = _resolve_output_path(target)
    if not path.exists():
        return None
    with open(path) as f:
        results = json.load(f)

    ext = target["extraction"]
    for r in results:
        phenotype = r.get("phenotype", "")
        if ext.get("phenotype_contains") and ext["phenotype_contains"] not in phenotype:
            continue
        if ext.get("phenotype_excludes") and ext["phenotype_excludes"] in phenotype:
            continue
        if ext.get("treatment") and r.get("treatment") != ext["treatment"]:
            continue
        return r.get(ext["field"])
    return None


def extract_spatial_csv(target: dict) -> float | None:
    """Extract observable from spatial depth_kill_curves.csv."""
    path = _resolve_output_path(target)
    if not path.exists():
        return None

    ext = target["extraction"]
    treatment = ext.get("treatment", "")
    depth_range = ext.get("depth_um_range", [0, 999999])
    field = ext.get("field", "death_rate")
    aggregation = ext.get("aggregation", "mean")

    values = []
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            if row.get("treatment") != treatment:
                continue
            depth = float(row.get("depth_um", 0))
            if depth < depth_range[0] or depth > depth_range[1]:
                continue
            n_cells = int(row.get("n_cells", 0))
            if n_cells == 0:
                continue
            values.append(float(row[field]))

    if not values:
        return None
    if aggregation == "min":
        return min(values)
    if aggregation == "max":
        return max(values)
    return sum(values) / len(values)


def extract_window_csv(target: dict) -> float | None:
    """Extract observable from vulnerability window CSV."""
    path = _resolve_output_path(target)
    if not path.exists():
        return None

    ext = target["extraction"]
    treatment = ext.get("treatment", "")
    timepoint = ext.get("timepoint_hours", 0)
    field = ext.get("field", "death_rate")

    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            if row.get("treatment") != treatment:
                continue
            if abs(float(row.get("timepoint_hours", 0)) - timepoint) > 0.5:
                continue
            return float(row[field])
    return None


def extract_invivo_json(target: dict) -> float | None:
    """Extract derived protection factor from invivo comparison JSON."""
    path = _resolve_output_path(target)
    if not path.exists():
        return None

    with open(path) as f:
        results = json.load(f)

    ext = target["extraction"]
    contexts = ext.get("contexts", [])
    treatment = ext.get("treatment", "")
    values_by_context = {}

    for r in results:
        phenotype = r.get("phenotype", "")
        if ext.get("phenotype_contains") and ext["phenotype_contains"] not in phenotype:
            continue
        if ext.get("phenotype_excludes") and ext["phenotype_excludes"] in phenotype:
            continue
        if r.get("treatment") != treatment:
            continue
        ctx = r.get("context", "")
        if ctx in contexts:
            values_by_context[ctx] = r.get("death_rate", 0)

    if "2d" in values_by_context and "invivo" in values_by_context:
        invivo_val = values_by_context["invivo"]
        if invivo_val > 0:
            return values_by_context["2d"] / invivo_val
    return None


EXTRACTORS = {
    "sim-original": extract_sim_original,
    "sim-spatial": extract_spatial_csv,
    "sim-window": extract_window_csv,
    "sim-invivo": extract_invivo_json,
}


def evaluate_target(target: dict) -> dict:
    """Evaluate a single calibration target against simulation output."""
    binary = target.get("binary", "")

    # Route to the right extractor
    if binary in EXTRACTORS:
        observed = EXTRACTORS[binary](target)
    else:
        observed = None

    if observed is None:
        return {
            "id": target["id"],
            "status": "SKIP",
            "reason": f"Output file not found or observable not extractable",
            "target": target.get("target_value"),
            "observed": None,
            "residual": None,
        }

    target_value = target["target_value"]
    tolerance = target["tolerance"]
    comparator = target.get("target_comparator")

    if comparator == "<":
        passed = observed < target_value
        residual = observed - target_value
    elif comparator == ">":
        passed = observed > target_value
        residual = target_value - observed
    else:
        residual = observed - target_value
        passed = abs(residual) <= tolerance

    return {
        "id": target["id"],
        "status": "PASS" if passed else "FAIL",
        "target": target_value,
        "observed": round(observed, 6),
        "residual": round(residual, 6),
        "tolerance": tolerance,
        "comparator": comparator or "==",
        "confidence": target.get("confidence", ""),
        "description": target.get("description", ""),
    }


def generate_report(results: list[dict]) -> str:
    """Generate markdown calibration report."""
    lines = ["# Calibration Report\n"]
    lines.append("Comparison of current simulation outputs against published calibration targets.\n")

    passed = sum(1 for r in results if r["status"] == "PASS")
    failed = sum(1 for r in results if r["status"] == "FAIL")
    skipped = sum(1 for r in results if r["status"] == "SKIP")
    lines.append(f"**Results: {passed} PASS, {failed} FAIL, {skipped} SKIP out of {len(results)} targets.**\n")

    lines.append("| Target | Status | Observed | Target | Residual | Tolerance | Confidence |")
    lines.append("|--------|--------|----------|--------|----------|-----------|------------|")
    for r in results:
        obs = f"{r['observed']:.4f}" if r["observed"] is not None else "N/A"
        res = f"{r['residual']:+.4f}" if r["residual"] is not None else "N/A"
        tgt = r.get("comparator", "==")
        if tgt in ("<", ">"):
            tgt_str = f"{tgt} {r['target']}"
        else:
            tgt_str = f"{r['target']}"
        lines.append(
            f"| **{r['id']}** | {r['status']} | {obs} | {tgt_str} | {res} | {r['tolerance']} | {r['confidence']} |"
        )

    lines.append("\n## Details\n")
    for r in results:
        emoji = {"PASS": "+", "FAIL": "!", "SKIP": "?"}[r["status"]]
        lines.append(f"- [{emoji}] **{r['id']}**: {r['description']}")
        if r["status"] == "SKIP":
            lines.append(f"  - Skipped: {r.get('reason', 'unknown')}")
        elif r["status"] == "FAIL":
            lines.append(f"  - Observed {r['observed']} vs target {r['target']} (residual {r['residual']:+.4f}, tolerance {r['tolerance']})")

    lines.append("\n## Interpretation\n")
    lines.append("- PASS means the current simulation output is within tolerance of the published target.")
    lines.append("- FAIL means the output is outside tolerance and the parameter may need adjustment.")
    lines.append("- SKIP means the required simulation output file was not found locally.")
    lines.append("- See `parameter_provenance.md` for the source and confidence of each parameter.")
    lines.append("- See `targets.yaml` for the full target definitions including source PMIDs.")
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser(description="Evaluate simulation calibration targets")
    parser.add_argument("--evaluate", action="store_true", help="Run evaluation against current outputs")
    args = parser.parse_args()

    if not args.evaluate:
        parser.print_help()
        print("\nUse --evaluate to compare current simulation outputs against calibration targets.")
        return

    if not TARGETS_FILE.exists():
        print(f"ERROR: Targets file not found: {TARGETS_FILE}")
        sys.exit(1)

    targets = load_targets()
    print(f"Loaded {len(targets)} calibration targets from {TARGETS_FILE.name}\n")

    results = []
    for target in targets:
        result = evaluate_target(target)
        status_mark = {"PASS": "PASS", "FAIL": "FAIL", "SKIP": "SKIP"}[result["status"]]
        obs_str = f"{result['observed']:.4f}" if result["observed"] is not None else "N/A"
        print(f"  [{status_mark}] {result['id']}: observed={obs_str}, target={result['target']}")
        results.append(result)

    report = generate_report(results)
    REPORT_FILE.write_text(report, encoding="utf-8")
    print(f"\nReport written to {REPORT_FILE}")

    passed = sum(1 for r in results if r["status"] == "PASS")
    total = len(results)
    print(f"Overall: {passed}/{total} targets passed")


if __name__ == "__main__":
    main()
