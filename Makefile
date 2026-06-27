# One-command reproduction + common developer tasks (#542).
#
#   make reproduce   # the headline target: full test suite + a regenerated result
#   make test        # Python + Rust test suites
#   make bindings    # build the PyO3 Python bindings (maturin)
#   make manifest    # regenerate + verify MANIFEST.sha256 (needs git)
#
# `make reproduce` is what the Dockerfile runs. The Rust suite includes the
# calibration-regression gate (tests/test_calibration_regression.py wired through
# cargo? no -> it is a pytest), so the full reproduction runs BOTH suites.
.PHONY: reproduce test test-python test-rust bindings manifest headline help

help:
	@echo "Targets:"
	@echo "  reproduce  - run both test suites + regenerate a data-anchored headline"
	@echo "  test       - Python (pytest) + Rust (cargo test --workspace)"
	@echo "  bindings   - build the PyO3 Python bindings via maturin"
	@echo "  manifest   - regenerate + verify MANIFEST.sha256 (requires git)"

reproduce: test headline
	@echo ""
	@echo "==> Reproduction complete: both test suites pass and a data-anchored"
	@echo "    headline (ferroptotic trigger-wave speed, #482) was regenerated."

test: test-python test-rust

test-python:
	python3 -m pytest tests/ -q

test-rust:
	cd simulations && cargo test --workspace

bindings:
	cd simulations && maturin develop -m ferroptosis-python/Cargo.toml --release

# A tangible, dependency-light reproduction: re-derive a result that is anchored
# to published measurement (pure stdlib, no compiled extension needed).
headline:
	@echo "==> Regenerating a data-anchored headline (ferroptotic trigger-wave, #482):"
	python3 scripts/validate_trigger_wave.py

# Needs git (generate_release_manifest.py uses `git ls-files`).
manifest:
	python3 scripts/generate_release_manifest.py
	@git diff --quiet -- MANIFEST.sha256 \
		&& echo "MANIFEST.sha256 is current." \
		|| (echo "ERROR: MANIFEST.sha256 drifted — commit the regenerated file." && exit 1)
