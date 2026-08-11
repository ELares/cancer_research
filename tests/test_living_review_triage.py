"""The living-review triage must not report a tagging artifact as a trend.

The living-review workflow uploads a monthly delta and never commits, leaving a
human to decide whether anything in it matters. Two windows accumulated unread,
so the triage exists to make that decision cheap.

The trap it was built around is real and fires on the very first delta: comparing
the new window's mechanism shares against the frozen corpus reports `cuproptosis`
at 0.0% frozen against 4.0% in the window, which reads as explosive emergence. It
is not. The frozen index was tagged with an older mechanism vocabulary, and 38 of
its records already name cuproptosis on a different tag axis.

These guards pin the separation, not the numbers.
"""

import json
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "scripts"))

import living_review_triage as lrt  # noqa: E402


def _frozen():
    return lrt.load_jsonl(REPO_ROOT / "corpus" / "INDEX.jsonl")


def test_the_frozen_vocabulary_is_smaller_than_the_current_config():
    """The premise: the frozen tags are a stale snapshot of the keyword list.

    If these ever agree, the confound is gone and the triage's separation is
    unnecessary -- which would be good news, and worth noticing rather than
    carrying dead machinery.
    """
    import config

    frozen_vocab = lrt.frozen_vocabulary(_frozen())
    current = set(getattr(config, "MECHANISM_KEYWORDS", {}) or {})
    if not current:
        return
    missing = current - frozen_vocab
    assert missing, (
        "the frozen index now covers every mechanism in the config, so the "
        "not-comparable split this triage performs is no longer needed")


def test_a_mechanism_outside_the_frozen_vocabulary_is_not_called_movement():
    """The whole point: a structurally-zero baseline is not a denominator."""
    frozen_vocab = lrt.frozen_vocabulary(_frozen())
    assert "cuproptosis" not in frozen_vocab, (
        "cuproptosis is now in the frozen vocabulary; re-check whether the "
        "confound this guards against still exists")


def test_the_untagged_search_looks_where_the_data_actually_is():
    """The frozen index has no abstract and no MeSH terms.

    An earlier version searched `abstract`, which does not exist there, and so
    reported 5 mentions instead of 38 -- a missing field reading as a small
    count. The search must cover the tag axes that carry the evidence.
    """
    froz = _frozen()
    assert "abstract" not in froz[0], (
        "the frozen index gained an abstract field; the search fields should be "
        "revisited rather than left as they are")
    n, where = lrt.untagged_mentions(froz, "cuproptosis")
    assert n >= 20, (
        f"only {n} frozen records found naming cuproptosis; the search is "
        "probably looking at a field that does not exist")
    assert "pathway_targets" in where, (
        "the evidence for the vocabulary gap lives in pathway_targets and the "
        "search is no longer finding it")


def test_the_triage_runs_and_separates_the_two_classes():
    """End to end on a synthetic delta, so it does not need the artifact."""
    import tempfile

    rows = [
        {"pmid": "1", "mechanisms": ["cuproptosis"], "pub_types": ["Journal Article"]},
        {"pmid": "2", "mechanisms": ["immunotherapy"],
         "pub_types": ["Randomized Controlled Trial"]},
    ]
    with tempfile.TemporaryDirectory() as d:
        delta = Path(d) / "index.jsonl"
        delta.write_text("\n".join(json.dumps(r) for r in rows) + "\n")
        out = Path(d) / "report.md"
        rc = subprocess.run(
            [sys.executable, str(REPO_ROOT / "scripts" / "living_review_triage.py"),
             "--delta", str(delta), "--out", str(out)],
            capture_output=True, text=True)
        assert rc.returncode == 0, rc.stderr
        txt = out.read_text()
    assert "NOT valid" in txt, "the not-comparable section is missing"
    assert "cuproptosis" in txt
    # immunotherapy is in the frozen vocabulary, so it must NOT be in the
    # not-comparable table.
    head = txt[:txt.index("## Comparable mechanisms")]
    assert "immunotherapy" not in head, (
        "a mechanism the frozen corpus was tagged with is being reported as "
        "not comparable")


def test_the_report_states_what_it_cannot_settle():
    """A five-week window against a decade is not a base rate."""
    import tempfile

    with tempfile.TemporaryDirectory() as d:
        delta = Path(d) / "index.jsonl"
        delta.write_text(json.dumps(
            {"pmid": "1", "mechanisms": ["immunotherapy"], "pub_types": []}) + "\n")
        out = Path(d) / "r.md"
        subprocess.run(
            [sys.executable, str(REPO_ROOT / "scripts" / "living_review_triage.py"),
             "--delta", str(delta), "--out", str(out)],
            capture_output=True, text=True, check=True)
        txt = " ".join(out.read_text().split())
    assert "not a like-for-like" in txt
    assert "indexing lag" in txt
    assert "says where to look, not what is true" in txt
