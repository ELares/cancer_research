#!/usr/bin/env python3
"""Generate the 2D-vs-3D TME comparison table from sim-tme + sim-tme-3d outputs.

Closes #195 AC: "Comparison table: 2D vs 3D for every TME feature."

**Scale-mismatch caveat (load-bearing)**: sim-tme uses a 500x500 grid
(tumor radius ~4500 um, ~159k tumor cells). sim-tme-3d uses a 60^3 grid
(tumor radius ~540 um, ~12k tumor cells). The 8x linear scale difference
means absolute kill counts are NOT directly comparable. This script
reports RATIOS (e.g., RSL3 hypoxic kill / RSL3 normoxic kill) which are
dimensionally meaningful at different scales.

Inputs:
- output/tme/tme_summary.json       (sim-tme, 2D)
- output/tme-3d/summary.json        (sim-tme-3d, 3D)

Outputs:
- output/tme-3d/comparison_2d_vs_3d.csv  — per-condition side-by-side
- output/tme-3d/key_questions.txt        — the 4 issue-#195 answers
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
SIM_TME_JSON = REPO_ROOT / "output" / "tme" / "tme_summary.json"
SIM_TME_3D_JSON = REPO_ROOT / "output" / "tme-3d" / "summary.json"
OUT_DIR = REPO_ROOT / "output" / "tme-3d"

# Conditions we compare. Each maps a logical key to the lookup criteria
# in both the 2D and 3D JSONs.
COMPARISON_KEYS = [
    # (key_name, treatment, o2_condition, immune_mode)
    ("baseline_control",    "Control", "uniform",  "off"),
    ("baseline_rsl3",       "RSL3",    "uniform",  "off"),
    ("baseline_sdt",        "SDT",     "uniform",  "off"),
    ("o2_120_control",      "Control", "gradient", "off"),  # o2_lambda_um=120
    ("o2_120_rsl3",         "RSL3",    "gradient", "off"),
    ("o2_120_sdt",          "SDT",     "gradient", "off"),
]


def load_2d(path: Path) -> list[dict]:
    """sim-tme writes a bare array of ConditionResult."""
    with open(path) as f:
        return json.load(f)


def load_3d(path: Path) -> list[dict]:
    """sim-tme-3d writes {grid_dim, ..., conditions: [...]}."""
    with open(path) as f:
        data = json.load(f)
    return data["conditions"]


def find_condition(conditions: list[dict], treatment: str, o2_condition: str,
                   immune_mode: str, o2_lambda: float | None = None) -> dict | None:
    """Find a condition by (treatment, o2_condition, immune_mode, optional λ).

    Returns None if not found (a condition may be skipped in one binary).
    """
    for c in conditions:
        if c.get("treatment") != treatment:
            continue
        # 2D uses "gradient_80um" etc; 3D uses "gradient" + o2_lambda_um field.
        cond_o2 = c.get("o2_condition", "")
        if o2_condition == "uniform" and cond_o2 != "uniform":
            continue
        if o2_condition == "gradient":
            # Accept either "gradient" (3D) or "gradient_*" (2D).
            if not cond_o2.startswith("gradient"):
                continue
            if o2_lambda is not None and c.get("o2_lambda_um") != o2_lambda:
                continue
        if c.get("immune_mode") != immune_mode:
            continue
        return c
    return None


def write_comparison_csv(out_path: Path, rows_2d: list[dict],
                         rows_3d: list[dict]) -> None:
    """Write a side-by-side CSV of key conditions.

    Columns: condition_key, metric, 2d_value, 3d_value, ratio_3d_over_2d.
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["condition_key,metric,2d_value,3d_value,ratio_3d_over_2d"]
    metrics = [
        "total_tumor",
        "total_dead",
        "overall_kill_rate",
        "normoxic_kill_rate",
        "transition_kill_rate",
        "hypoxic_kill_rate",
    ]
    for key_name, tx, o2c, im in COMPARISON_KEYS:
        # λ=120 specifically for the gradient cases.
        lam = 120.0 if o2c == "gradient" else None
        r2 = find_condition(rows_2d, tx, o2c, im, lam)
        r3 = find_condition(rows_3d, tx, o2c, im, lam)
        for metric in metrics:
            v2 = r2.get(metric) if r2 else None
            v3 = r3.get(metric) if r3 else None
            ratio = ""
            if (isinstance(v2, (int, float)) and isinstance(v3, (int, float))
                    and v2 != 0):
                ratio = f"{v3 / v2:.4f}"
            v2s = f"{v2}" if v2 is not None else ""
            v3s = f"{v3}" if v3 is not None else ""
            lines.append(f"{key_name},{metric},{v2s},{v3s},{ratio}")
    out_path.write_text("\n".join(lines) + "\n")
    print(f"Wrote {out_path}")


def answer_key_questions(rows_2d: list[dict], rows_3d: list[dict],
                          out_path: Path) -> None:
    """Compute the 4 manuscript-keystone answers from the issue.

    Each answer compares a RATIO (so it's scale-invariant) between 2D and 3D.
    """
    lines = []
    lines.append("=" * 72)
    lines.append("sim-tme-3d / sim-tme: 4 key questions from issue #195")
    lines.append("=" * 72)
    lines.append("")
    lines.append("Caveat: 2D uses 500x500 grid (radius ~4500um, ~159k tumor cells);")
    lines.append("3D uses 60^3 grid (radius ~540um, ~12k tumor cells).")
    lines.append("Compare RATIOS not absolute counts.")
    lines.append("")

    # --- Q1: Hypoxia RSL3 collapse ---
    # 2D baseline-RSL3 vs RSL3-at-lambda-120 (the canonical hypoxia case)
    def kill_rate(rows, tx, o2c, lam):
        c = find_condition(rows, tx, o2c, "off", lam)
        return c["overall_kill_rate"] if c else None

    base_2d = kill_rate(rows_2d, "RSL3", "uniform", None)
    hyp_2d = kill_rate(rows_2d, "RSL3", "gradient", 120.0)
    base_3d = kill_rate(rows_3d, "RSL3", "uniform", None)
    hyp_3d = kill_rate(rows_3d, "RSL3", "gradient", 120.0)
    lines.append("Q1: Hypoxia RSL3 collapse — RSL3 kill: baseline vs O2-gradient")
    if base_2d and hyp_2d and base_3d and hyp_3d:
        ratio_2d = hyp_2d / base_2d
        ratio_3d = hyp_3d / base_3d
        lines.append(
            f"  2D: baseline={base_2d:.4f} → λ=120={hyp_2d:.4f} (ratio={ratio_2d:.3f}×)")
        lines.append(
            f"  3D: baseline={base_3d:.4f} → λ=120={hyp_3d:.4f} (ratio={ratio_3d:.3f}×)")
        if ratio_3d < ratio_2d:
            lines.append("  → 3D collapse is MORE pronounced than 2D (smaller ratio).")
        elif ratio_3d > ratio_2d:
            lines.append("  → 3D collapse is LESS pronounced than 2D (larger ratio).")
        else:
            lines.append("  → ratios are essentially equal.")
    else:
        lines.append("  (incomplete data — check both summaries)")
    lines.append("")

    # --- Q2: Immune 104:1 ratio (SDT immune kills vs RSL3 immune kills) ---
    def immune_kills(rows, tx):
        c = find_condition(rows, tx, "gradient", "on", 120.0)
        return c.get("immune_kills") if c else None

    sdt_imm_2d = immune_kills(rows_2d, "SDT")
    rsl3_imm_2d = immune_kills(rows_2d, "RSL3")
    sdt_imm_3d = immune_kills(rows_3d, "SDT")
    rsl3_imm_3d = immune_kills(rows_3d, "RSL3")
    lines.append("Q2: Immune SDT-vs-RSL3 ratio (the 104:1 manuscript claim)")
    if sdt_imm_2d and rsl3_imm_2d and rsl3_imm_2d > 0:
        lines.append(f"  2D: SDT={sdt_imm_2d}, RSL3={rsl3_imm_2d} (ratio={sdt_imm_2d / rsl3_imm_2d:.1f}×)")
    else:
        lines.append("  2D: (immune_on data missing)")
    if sdt_imm_3d and rsl3_imm_3d and rsl3_imm_3d > 0:
        lines.append(f"  3D: SDT={sdt_imm_3d}, RSL3={rsl3_imm_3d} (ratio={sdt_imm_3d / rsl3_imm_3d:.1f}×)")
    elif sdt_imm_3d is not None and rsl3_imm_3d is not None:
        lines.append(f"  3D: SDT={sdt_imm_3d}, RSL3={rsl3_imm_3d} (cannot compute ratio — RSL3=0)")
    else:
        lines.append("  3D: (immune_on data missing)")
    lines.append("")

    # --- Q3: Stromal shielding impact (stromal_adjacent kill rate vs overall) ---
    def stromal_data(rows, tx):
        # Find the stromal-on condition.
        for c in rows:
            if (c.get("treatment") == tx
                    and c.get("stromal_mode") == "stromal_on"
                    and c.get("o2_lambda_um") == 120.0):
                return c
        return None

    s2d = stromal_data(rows_2d, "RSL3")
    s3d = stromal_data(rows_3d, "RSL3")
    lines.append("Q3: Stromal shielding impact on RSL3")
    if s2d and s3d:
        # The "impact" is: overall - stromal_adjacent (the shielded cells survive more)
        overall_2d = s2d["overall_kill_rate"]
        adj_2d = s2d.get("stromal_adjacent_kill_rate", 0)
        overall_3d = s3d["overall_kill_rate"]
        adj_3d = s3d.get("stromal_adjacent_kill_rate", 0)
        impact_2d = (overall_2d - adj_2d) if adj_2d is not None else 0
        impact_3d = (overall_3d - adj_3d) if adj_3d is not None else 0
        lines.append(
            f"  2D: overall_kill={overall_2d:.4f}, stromal_adjacent_kill={adj_2d:.4f}, gap={impact_2d:.4f}")
        lines.append(
            f"  3D: overall_kill={overall_3d:.4f}, stromal_adjacent_kill={adj_3d:.4f}, gap={impact_3d:.4f}")
        lines.append(f"  → 3D boundary-fraction is larger (~3/R vs 2/R), so the SHIELDING fraction differs.")
    else:
        lines.append("  (incomplete data)")
    lines.append("")

    # --- Q4: pH ion trapping RSL3 reduction ---
    def ph_data(rows, tx):
        for c in rows:
            if (c.get("treatment") == tx
                    and c.get("ph_mode") == "ph_on"
                    and c.get("o2_lambda_um") == 120.0):
                return c
        return None

    rsl3_no_ph_2d = find_condition(rows_2d, "RSL3", "gradient", "off", 120.0)
    rsl3_ph_2d = ph_data(rows_2d, "RSL3")
    rsl3_no_ph_3d = find_condition(rows_3d, "RSL3", "gradient", "off", 120.0)
    rsl3_ph_3d = ph_data(rows_3d, "RSL3")
    lines.append("Q4: pH ion-trapping reduction for RSL3 (kill with vs without pH)")
    if rsl3_no_ph_2d and rsl3_ph_2d and rsl3_no_ph_3d and rsl3_ph_3d:
        no_2d = rsl3_no_ph_2d["overall_kill_rate"]
        ph_2d = rsl3_ph_2d["overall_kill_rate"]
        no_3d = rsl3_no_ph_3d["overall_kill_rate"]
        ph_3d = rsl3_ph_3d["overall_kill_rate"]
        red_2d = (1 - ph_2d / no_2d) if no_2d > 0 else 0
        red_3d = (1 - ph_3d / no_3d) if no_3d > 0 else 0
        lines.append(
            f"  2D: no-pH={no_2d:.4f}, pH={ph_2d:.4f} (reduction={red_2d * 100:.1f}%)")
        lines.append(
            f"  3D: no-pH={no_3d:.4f}, pH={ph_3d:.4f} (reduction={red_3d * 100:.1f}%)")
    else:
        lines.append("  (incomplete data)")
    lines.append("")
    lines.append("=" * 72)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(lines) + "\n")
    print(f"Wrote {out_path}")
    print()
    print("\n".join(lines))


def main() -> int:
    if not SIM_TME_JSON.exists():
        print(f"ERROR: missing {SIM_TME_JSON} (run sim-tme first)", file=sys.stderr)
        return 1
    if not SIM_TME_3D_JSON.exists():
        print(f"ERROR: missing {SIM_TME_3D_JSON} (run sim-tme-3d first)", file=sys.stderr)
        return 1

    rows_2d = load_2d(SIM_TME_JSON)
    rows_3d = load_3d(SIM_TME_3D_JSON)
    print(f"Loaded {len(rows_2d)} 2D conditions, {len(rows_3d)} 3D conditions")

    write_comparison_csv(OUT_DIR / "comparison_2d_vs_3d.csv", rows_2d, rows_3d)
    answer_key_questions(rows_2d, rows_3d, OUT_DIR / "key_questions.txt")
    return 0


if __name__ == "__main__":
    sys.exit(main())
