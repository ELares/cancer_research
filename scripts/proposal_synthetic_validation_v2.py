#!/usr/bin/env python3
"""Fixed second synthetic study; preserves the first failed study unchanged."""

import argparse
import importlib
import json
from pathlib import Path

import synthetic_proposal_study as study

ROOT = Path(__file__).resolve().parents[1]
OUT_JSON = ROOT / "analysis" / "calibration" / "proposal-synthetic-validation-v2.json"
OUT_MD = ROOT / "analysis" / "calibration" / "proposal-synthetic-validation-v2.md"
SPEC = study.StudySpec(
    revision_id="local-move-v2", strategy_module="resample_move_local",
    positive_seeds=(2026092401, 2026092402, 2026092403),
    source_paths=("scripts/proposal_synthetic_validation_v2.py", "scripts/synthetic_proposal_study.py",
                  "scripts/resample_move_local.py", "scripts/proposal_synthetic_validation.py",
                  "scripts/resample_move.py", "scripts/bounded_proposal.py", "scripts/importance_sampling.py"),
    prior_study_path="analysis/calibration/proposal-synthetic-validation.json",
    prior_study_sha256="1b045df4b719b1131c9c15256e0a803959c4b57a7df787e17cbf26e311ae3e82",
)


def strategy():
    # Only this checked-in specification chooses code; archived JSON cannot.
    return importlib.import_module(SPEC.strategy_module)


def assemble(stored):
    return study.assemble(stored, SPEC, strategy())


def render(result):
    return study.render(result)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--render-only", action="store_true")
    args = parser.parse_args()
    if args.render_only:
        result = assemble(json.loads(OUT_JSON.read_text()))
    else:
        result = study.run_study(SPEC, strategy())
    json_text = json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n"
    md_text = render(result)
    OUT_JSON.write_text(json_text)
    OUT_MD.write_text(md_text)
    print(f"wrote {OUT_MD}; all checks passed: {result['passed']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
