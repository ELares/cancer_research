#!/usr/bin/env python3
"""Make matplotlib PDF output reproducible.

A PDF written by matplotlib embeds `/CreationDate`, so regenerating a figure
whose data has not changed still produces different bytes. That has one
consequence worth fixing: **figure freshness cannot be checked**. Every other
committed artifact in this repo is now gated by regenerating it and comparing
(`tests/test_artifact_freshness.py`), and figures were the one class where the
comparison always failed, so a genuinely stale figure looked exactly like a
fresh one.

Measured before fixing, on `fig28_census_capture.pdf`: regenerating produced a
file of identical length differing at exactly one offset, inside
`/CreationDate (D:2026...)`. PNGs already reproduce -- matplotlib writes no
timestamp into them -- so this is a PDF-only problem.

WHY A WRAPPER RATHER THAN 59 CALL-SITE EDITS. `metadata={"CreationDate": None}`
has to reach every `savefig` that writes a PDF, and the two figure scripts have
59 of them. Editing each one leaves the next one somebody adds silently
nondeterministic again, which is the failure this exists to end. Wrapping the
method once covers the call sites that exist and the ones that do not yet.

It is deliberately NARROW: it only fills in a `CreationDate` that the caller
did not set, only for `.pdf` targets, and it leaves every other argument alone.
A caller that wants its own metadata still gets it.
"""
from pathlib import Path

from matplotlib.figure import Figure

_PATCHED = False


def _is_pdf(fname) -> bool:
    if isinstance(fname, (str, Path)):
        return str(fname).lower().endswith(".pdf")
    return False


def make_figures_deterministic() -> None:
    """Idempotently wrap `Figure.savefig` so PDFs carry no creation date."""
    global _PATCHED
    if _PATCHED:
        return
    original = Figure.savefig

    def savefig(self, fname, *args, **kwargs):
        if _is_pdf(fname):
            meta = dict(kwargs.get("metadata") or {})
            meta.setdefault("CreationDate", None)
            kwargs["metadata"] = meta
        return original(self, fname, *args, **kwargs)

    savefig.__doc__ = original.__doc__
    savefig._deterministic_wrapper = True          # so a guard can find it
    Figure.savefig = savefig
    _PATCHED = True
