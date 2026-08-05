"""Guards for the off-by-default authority-name filter (#628).

It changes which entities the co-mention layer can see, so the properties that
matter are: it does nothing unless asked, it touches only MeSH, it uses the
normalised comparison the measurements were made with, and it fails OPEN rather
than silently emptying the map when the label table is absent.
"""

import os
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "scripts"))

import atlas_comention as ac  # noqa: E402


def test_the_normalised_comparison_matches_how_mesh_and_papers_differ():
    """MeSH says `Breast Neoplasms` and stores `Cancer of Breast`; papers write
    `breast cancer`. Without all three normalisations the rule deletes three
    quarters of the cancer vocabulary."""
    b = ac._authority_bag
    assert b("Breast Neoplasms") == b("breast cancer")   # synonym + plural
    assert b("Cancer of Breast") == b("breast cancer")   # stopword + inversion
    assert b("Lung Neoplasms") == b("lung cancer")
    # It must not collapse genuinely different concepts.
    assert b("lung cancer") != b("breast cancer")
    assert b("Malignant Neoplasm of Breast") != b("breast cancer")


def test_the_filter_is_on_by_default_with_a_working_escape_hatch(monkeypatch):
    """Promoted on measured evidence (#628); `=0` still turns it off.

    An escape hatch that does not work is worse than none, because the first
    person to need it will assume the layer is unfiltered when it is not.
    """
    monkeypatch.setattr(ac, "_authority_labels",
                        lambda: {"MESH:D016609": ["Neoplasms, Second Primary"]})
    idx = {"alias": {"gpx4": "2879", "treatment": "MESH:D016609"},
           "alias_support": {}, "ident_mentions": {}}

    monkeypatch.delenv("FERRO_COMENTION_AUTHORITY", raising=False)
    out, stats = ac.build_alias_map(idx)
    assert "treatment" not in out, "the filter is not on by default"
    assert stats["dropped_not_a_name"] == 1

    monkeypatch.setenv("FERRO_COMENTION_AUTHORITY", "0")
    out, stats = ac.build_alias_map(idx)
    assert set(out) == {"gpx4", "treatment"}, "the escape hatch does not disable it"
    assert stats["dropped_not_a_name"] == 0


def test_it_fails_open_when_the_label_table_is_missing(monkeypatch, capsys):
    """Silently emptying the MeSH half of the map would be far worse than
    doing nothing, so an absent table must warn and pass everything through."""
    monkeypatch.setattr(ac, "_authority_labels", lambda: {})
    monkeypatch.delenv("FERRO_COMENTION_AUTHORITY", raising=False)
    idx = {"alias": {"treatment": "MESH:D016609"}, "alias_support": {},
           "ident_mentions": {}}
    out, stats = ac.build_alias_map(idx)
    assert set(out) == {"treatment"}, "the filter dropped forms with no table"
    assert stats["dropped_not_a_name"] == 0
    assert "proceeding UNFILTERED" in capsys.readouterr().err


def test_it_touches_only_mesh_identifiers(monkeypatch):
    """Genes are already 75% precise; the rule removes no false positives there
    and still costs true ones, so applying it to them is a measured mistake."""
    monkeypatch.setattr(ac, "_authority_labels",
                        lambda: {"MESH:D001943": ["Breast Neoplasms"],
                                 "MESH:D016609": ["Neoplasms, Second Primary"],
                                 "2879": ["GPX4", "glutathione peroxidase 4"]})
    monkeypatch.delenv("FERRO_COMENTION_AUTHORITY", raising=False)
    idx = {"alias": {"breast cancer": "MESH:D001943",   # a name -> kept
                     "treatment": "MESH:D016609",       # not a name -> dropped
                     "phgpx": "2879"},                  # gene, not a listed name
           "alias_support": {}, "ident_mentions": {}}
    out, stats = ac.build_alias_map(idx)
    assert "breast cancer" in out, "a genuine MeSH name was rejected"
    assert "treatment" not in out, "a generic word survived the filter"
    assert "phgpx" in out, "a gene form was filtered; the rule is MeSH-only"
    assert stats["dropped_not_a_name"] == 1


def test_an_unlabelled_mesh_identifier_is_kept_not_dropped():
    """Absence of a label is not evidence the form is wrong. 244 alias forms
    resolve to identifiers the table does not cover."""
    import types

    monkey = types.SimpleNamespace()
    old = ac._authority_labels
    ac._authority_labels = lambda: {"MESH:D001943": ["Breast Neoplasms"]}
    os.environ.pop("FERRO_COMENTION_AUTHORITY", None)
    try:
        idx = {"alias": {"whatever": "MESH:D999999"}, "alias_support": {},
               "ident_mentions": {}}
        out, stats = ac.build_alias_map(idx)
        assert "whatever" in out
        assert stats["dropped_not_a_name"] == 0
    finally:
        ac._authority_labels = old
        os.environ.pop("FERRO_COMENTION_AUTHORITY", None)


def test_the_curated_redirects_survive_the_filter():
    """The check runs AFTER the DOMAIN_SENSE redirect, and all five redirects
    target genes, so a hand-curated correction cannot be silently undone."""
    from atlas_ambiguity import DOMAIN_SENSE

    from build_label_source import load_table
    labels = load_table()
    if not labels:
        pytest.skip("authority table not built")
    for form, v in DOMAIN_SENSE.items():
        ident = v["id"] if isinstance(v, dict) else v[0]
        assert not ident.startswith("MESH:"), (
            f"{form} redirects to a MeSH identifier; the filter could undo it")
        names = labels.get(ident, [])
        assert any(ac._authority_bag(n) == ac._authority_bag(form) for n in names), (
            f"{form} is not an authority name of its redirect target {ident}")
