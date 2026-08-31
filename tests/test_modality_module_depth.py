"""Guards for `analysis/modality-module-depth.{md,json}`.

This page exists because a sentence in Chapter 6 was a measurement written as
an assertion, was true when written, and stopped being true silently. So the
guards are about the two ways it could go wrong again: the count drifting from
the crate, and the page being read as progress.
"""
import importlib.util
import json
import re
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
MD = REPO / "analysis/modality-module-depth.md"
JSON_ = REPO / "analysis/modality-module-depth.json"
CORE = REPO / "simulations/ferroptosis-core/src"
MANUSCRIPT = REPO / "article/drafts/v1.md"


def _load():
    spec = importlib.util.spec_from_file_location(
        "modality_module_depth", REPO / "scripts/modality_module_depth.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


MD_ = _load()


@pytest.fixture(scope="module")
def d():
    return json.loads(JSON_.read_text())


def test_the_counts_are_recomputed_not_stored(d):
    """A hand-edited number is the whole failure mode this page addresses."""
    live = MD_.assemble(MD_.scan())
    for k in ("dedicated_modules", "dedicated_pub_fns", "dedicated_code_lines",
              "ferroptosis_engine_modules", "engine_pub_fns",
              "engine_code_lines"):
        assert live[k] == d[k], f"{k} drifted: live {live[k]} vs stored {d[k]}"


def test_every_named_module_exists_and_every_module_is_classified(d):
    """A module missing from all three buckets would be uncounted, and one
    named that does not exist would be counted twice over."""
    for m in d["dedicated"] + d["shared"]:
        assert (CORE / f"{m['module']}.rs").exists(), m["module"]
        assert m["serves"].strip(), m["module"]
    on_disk = {p.stem for p in CORE.glob("*.rs") if p.stem != "lib"}
    named = {m["module"] for m in d["dedicated"] + d["shared"]}
    assert named <= on_disk, f"named but absent: {sorted(named - on_disk)}"
    assert d["total_modules"] == len(on_disk)
    assert (d["ferroptosis_engine_modules"] + len(named)) == len(on_disk), (
        "the buckets do not partition the crate, so some module is counted "
        "twice or not at all")


def test_the_ferroptosis_engine_is_still_the_larger_body_of_work(d):
    """The page's whole point is honesty about the gap.

    If the arms ever DO overtake the engine this fails, and the chapter's
    paragraph has to be rewritten rather than the number quietly updated --
    which is the failure that produced this page.
    """
    assert d["engine_code_lines"] > d["dedicated_code_lines"], (
        "the modality arms now exceed the ferroptosis engine in code; "
        "Chapter 6's paragraph claims the opposite and must be re-derived")
    assert d["engine_pub_fns"] > d["dedicated_pub_fns"]
    assert "still the larger body of work" in MD.read_text()


def test_the_page_refuses_to_be_read_as_progress(d):
    """Lines of code are the weakest possible evidence of depth. The page says
    so, and says what it does not measure."""
    md = MD.read_text()
    for frag in ("weakest possible evidence of depth",
                 "A count going up is not progress on its own",
                 "A module can be large and wrong",
                 "none of this appears in a figure or a claim"):
        assert frag in md, f"the page no longer says: {frag}"


def test_the_chapter_quotes_the_measurement_and_not_the_old_sentence(d):
    """The sentence this page replaced must be GONE, not merely contradicted
    somewhere else -- a retraction that leaves the original standing is how
    this repository has shipped two contradictory claims before."""
    text = MANUSCRIPT.read_text()
    assert "one function and a configuration struct each" in text, (
        "the chapter no longer records what the old sentence said, so a "
        "reader cannot tell the claim was retracted")
    # ...but only as a QUOTED retraction, never as a live claim.
    live = re.search(
        r"\*\*Most of these arms are one function and a configuration struct\.\*\*",
        text)
    assert live is None, (
        "the retracted sentence is still standing as a live claim")
    for k in ("dedicated_modules", "dedicated_pub_fns", "ferroptosis_engine_modules"):
        assert str(d[k]) in text, (
            f"the chapter does not quote the measured {k}, so it can drift "
            "from the artifact again")
