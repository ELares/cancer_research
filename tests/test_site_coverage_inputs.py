"""Exercise real site scans and all three outputs with tiny isolated inputs."""
import gzip
import importlib.util
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest


REPO = Path(__file__).resolve().parent.parent
SENTINELS = {
    "OUT_JSON": b'{"existing": "published site counts"}\n',
    "OUT_MD": b"Existing site interpretation.\n",
    "OUT_MAP": b"Existing site descriptor map.\n",
}
LUNG = {"pmid": "1", "mesh": ["Lung Neoplasms"], "cancer_basis": "C04"}
GENERIC = {"pmid": "2", "mesh": ["Neoplasms"], "cancer_basis": "C04"}


def _write_shard(path, records):
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8") as stream:
        for record in records:
            stream.write(json.dumps(record) + "\n")


def _snapshot(module):
    return {attr: getattr(module, attr).read_bytes() for attr in SENTINELS}


@pytest.fixture
def scanner(tmp_path, monkeypatch):
    root = tmp_path / "external-atlas"
    monkeypatch.setenv("FERRO_ATLAS_ROOT", str(root))
    monkeypatch.syspath_prepend(str(REPO / "scripts"))
    spec = importlib.util.spec_from_file_location(
        "site_coverage_inputs", REPO / "scripts" / "atlas_site_coverage.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    # Keep regressions against the old hardcoded reader isolated too. The
    # public-root test below deliberately makes this legacy location differ.
    monkeypatch.setattr(module, "ATLAS", root)
    for attr, content in SENTINELS.items():
        path = tmp_path / attr
        path.write_bytes(content)
        monkeypatch.setattr(module, attr, path)
    monkeypatch.setattr(sys, "argv", ["atlas_site_coverage.py"])
    return SimpleNamespace(module=module, root=root, records=root / "records")


@pytest.mark.parametrize("state", ["missing", "empty-directory", "empty-shards"])
def test_unavailable_census_preserves_all_three_outputs(scanner, state):
    if state != "missing":
        scanner.records.mkdir(parents=True)
    if state == "empty-shards":
        _write_shard(scanner.records / "a.jsonl.gz", [])
        _write_shard(scanner.records / "b.jsonl.gz", [])

    with pytest.raises(SystemExit) as exc:
        scanner.module.main()

    assert _snapshot(scanner.module) == SENTINELS
    assert str(scanner.records) in str(exc.value)
    assert "FERRO_ATLAS_ROOT" in str(exc.value)
    assert "--render-only" in str(exc.value)


@pytest.mark.parametrize("damage", ["json", "gzip", "truncated-gzip"])
def test_late_corruption_preserves_all_three_outputs(scanner, damage):
    _write_shard(scanner.records / "a.jsonl.gz", [LUNG])
    tail = scanner.records / "z.jsonl.gz"
    if damage == "json":
        with gzip.open(tail, "wt", encoding="utf-8") as stream:
            stream.write('{"pmid":\n')
    elif damage == "gzip":
        tail.write_bytes(b"not a gzip file")
    else:
        _write_shard(tail, [GENERIC])
        tail.write_bytes(tail.read_bytes()[:-6])

    with pytest.raises((ValueError, OSError, EOFError)):
        scanner.module.main()

    assert _snapshot(scanner.module) == SENTINELS


@pytest.mark.parametrize("record", [
    [], {"mesh": "Lung Neoplasms"}, {"mesh": [42]},
    {"mesh": {"Lung Neoplasms": True}}, {"mesh": [], "cancer_basis": []},
])
def test_malformed_record_shape_preserves_all_three_outputs(scanner, record):
    _write_shard(scanner.records / "a.jsonl.gz", [LUNG, record])

    with pytest.raises(ValueError, match="Invalid census record"):
        scanner.module.main()

    assert _snapshot(scanner.module) == SENTINELS


@pytest.mark.parametrize("renderer", ["render", "render_map"])
def test_render_failure_preserves_all_three_outputs(scanner, monkeypatch, renderer):
    _write_shard(scanner.records / "a.jsonl.gz", [LUNG])

    def fail_render(_data):
        raise RuntimeError("intentional renderer failure")

    monkeypatch.setattr(scanner.module, renderer, fail_render)
    with pytest.raises(RuntimeError, match="intentional renderer failure"):
        scanner.module.main()

    assert _snapshot(scanner.module) == SENTINELS


def test_invalid_optional_manifest_preserves_all_three_outputs(scanner):
    _write_shard(scanner.records / "a.jsonl.gz", [LUNG])
    (scanner.root / "unindexed-manifest.json").write_text("{")

    with pytest.raises(ValueError):
        scanner.module.main()

    assert _snapshot(scanner.module) == SENTINELS


def test_public_root_controls_records_and_optional_manifests(scanner, monkeypatch):
    module = scanner.module
    monkeypatch.setattr(module, "ATLAS", scanner.root / "unused-legacy-root")
    _write_shard(scanner.records / "a.jsonl.gz", [LUNG, GENERIC])
    _write_shard(scanner.root / "records_unindexed" / "ignored.jsonl.gz", [LUNG])
    _write_shard(scanner.root / "records_updates" / "ignored.jsonl.gz", [LUNG])
    (scanner.root / "unindexed-manifest.json").write_text(
        json.dumps({"files": {"fixture": {"cancer_text": 3}}}))
    (scanner.root / "manifest-c04only.json").write_text(
        json.dumps({"files": {"fixture": {"cancer": 2}}}))

    module.main()

    result = json.loads(module.OUT_JSON.read_text())
    assert result["census"] == 2
    assert result["assigned"] == 1
    assert result["excluded_streams"] == {"text_matched_no_mesh": 3, "c04_core": 2}
    assert dict(result["sites"])["lung"] == 1
    assert len(result["sites"]) == result["n_sites"] == len(module.SITES)
    assert module.OUT_MD.read_text() == module.render(result)
    assert module.OUT_MAP.read_text() == module.render_map(module.deep_map())
    assert "**20.0%**" in module.OUT_MD.read_text()


def test_zero_assignments_are_valid_and_findings_are_conditional(scanner):
    module = scanner.module
    _write_shard(scanner.records / "a.jsonl.gz", [GENERIC])

    module.main()

    result = json.loads(module.OUT_JSON.read_text())
    md = module.OUT_MD.read_text()
    assert result["census"] == 1 and result["assigned"] == 0
    assert result["adjacent_basis"] == 0
    assert result["unassigned"]["generic_neoplasms"] == 1
    assert set(dict(result["sites"]).values()) == {0}
    assert "| lung | 0 | 0 | n/a (zero shallow) | n/a -> n/a |" in md
    assert "No article in this input matches the shallow site list" in md
    assert "n/a (no assigned articles)" in md
    assert "no finite largest/smallest ratio" in md
    assert "is the largest single bucket" not in md
    assert "sums to more than the assigned total" not in md
    assert "The gap between the two lists is not uniform" not in md
    assert "not the census's limit" not in md


def test_deep_only_site_is_visible_without_an_invented_ratio_or_shallow_rank(scanner):
    module = scanner.module
    dm = module.deep_map()
    deep_only = sorted(dm["lung"]["deep"] - module.SITES["lung"])[0]
    _write_shard(scanner.records / "a.jsonl.gz", [
        {"mesh": [deep_only], "cancer_basis": "C04"},
    ])

    module.main()

    result = json.loads(module.OUT_JSON.read_text())
    assert result["assigned"] == 0
    assert result["variants"]["deep"]["assigned"] == 1
    assert result["unassigned"]["same_sites_deeper"] == 1
    assert "| lung | 0 | 1 | n/a (zero shallow) | n/a -> 1 |" in module.OUT_MD.read_text()


def test_missing_optional_metadata_uses_actual_c04_membership(scanner):
    module = scanner.module
    _write_shard(scanner.records / "a.jsonl.gz", [
        {"mesh": ["Neoplasms"]}, {"mesh": ["Lung Neoplasms"]}, {"mesh": None}, {},
    ])

    result = module.scan()

    assert result["census"] == 4 and result["assigned"] == 1
    assert result["adjacent_basis"] == 2
    assert result["unassigned"]["total"] == 3
    assert result["unassigned"]["generic_neoplasms"] == 1
    assert result["unassigned"]["no_c04_descriptor"] == 2
    assert "the remaining 0 (0.0%)" in module.render(result)
    assert _snapshot(module) == SENTINELS


def test_deep_only_site_can_change_ranks_when_existing_ratios_are_uniform(scanner):
    module = scanner.module
    deep_only = sorted(module.deep_map()["breast"]["deep"] - module.SITES["breast"])[0]
    _write_shard(scanner.records / "a.jsonl.gz", [
        LUNG, {"mesh": [deep_only]}, {"mesh": [deep_only]},
    ])

    result = module.scan()
    md = module.render(result)

    assert result["assigned"] == 1 and result["variants"]["deep"]["assigned"] == 3
    assert "| lung | 1 | 1 | 1.00x | **1 -> 2** |" in md
    assert "| breast | 0 | 2 | n/a (zero shallow) | n/a -> 1 |" in md
    assert "1 of 1 rankable sites change rank, `lung` 1 -> 2" in md


def test_sparse_ties_are_reproducible_and_share_ranks(scanner):
    module = scanner.module
    rows = [
        {"mesh": ["Breast Neoplasms", "Lung Neoplasms"], "cancer_basis": "C04"},
        {"mesh": ["Neoplasms", "Z fixture", "A fixture"], "cancer_basis": "C04"},
    ]
    shard = scanner.records / "a.jsonl.gz"
    _write_shard(shard, rows)
    first = module.scan()
    _write_shard(shard, [{**r, "mesh": list(reversed(r["mesh"]))} for r in reversed(rows)])
    second = module.scan()

    assert first == second
    md = module.render(first)
    assert md == module.render(json.loads(json.dumps(first, sort_keys=True)))
    for site in ("breast", "lung"):
        assert f"| {site} | 1 | 1 | 1.00x | 1 -> 1 |" in md
    assert "Tied positive counts share a rank" in md
    assert "ratio is uniformly 1.00x" in md
    assert "sums to more than the assigned total" in md
    assert first["unassigned"]["top_descriptors"] == [
        ["a fixture", 1], ["neoplasms", 1], ["z fixture", 1],
    ]
    assert _snapshot(module) == SENTINELS


@pytest.mark.parametrize("failure", [False, True], ids=["success", "render-failure"])
def test_render_only_needs_no_raw_input_or_taxonomy(scanner, monkeypatch, failure):
    module = scanner.module
    original = (REPO / "analysis" / "atlas-site-coverage.json").read_bytes()
    module.OUT_JSON.write_bytes(original)
    before = _snapshot(module)
    monkeypatch.setattr(sys, "argv", ["atlas_site_coverage.py", "--render-only"])
    monkeypatch.setattr(module, "C04_TSV", scanner.root / "missing-descriptors.tsv")
    monkeypatch.setattr(module, "C04_TREE", scanner.root / "missing-tree.tsv")

    def fail_scan(*_args):
        raise AssertionError("render-only reopened the census or taxonomy")

    monkeypatch.setattr(module, "scan", fail_scan)
    monkeypatch.setattr(module, "deep_map", fail_scan)
    assert not scanner.records.exists()
    if failure:
        def fail_render(_data):
            raise RuntimeError("intentional renderer failure")
        monkeypatch.setattr(module, "render", fail_render)
        with pytest.raises(RuntimeError, match="intentional renderer failure"):
            module.main()
        assert _snapshot(module) == before
    else:
        module.main()
        assert module.OUT_MD.read_text() == module.render(json.loads(original))
        assert module.OUT_MD.read_bytes() != before["OUT_MD"]
        assert module.OUT_JSON.read_bytes() == original
        assert module.OUT_MAP.read_bytes() == before["OUT_MAP"]


@pytest.mark.parametrize('manifest', [
    {}, {'files': []}, {'files': {'fixture': {}}},
    {'files': {'fixture': {'cancer_text': -1}}},
    {'files': {'fixture': {'cancer_text': True}}},
    {'files': {'fixture': {'cancer_text': 0.5}}},
])
def test_invalid_manifest_population_preserves_all_outputs(scanner, manifest):
    _write_shard(scanner.records / 'a.jsonl.gz', [LUNG, GENERIC])
    (scanner.root / 'unindexed-manifest.json').write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match='manifest'):
        scanner.module.main()
    assert _snapshot(scanner.module) == SENTINELS


def test_core_manifest_must_match_the_scanned_population(scanner):
    _write_shard(scanner.records / 'a.jsonl.gz', [LUNG, LUNG, {'mesh': []}])
    (scanner.root / 'manifest-c04only.json').write_text(
        json.dumps({'files': {'fixture': {'cancer': 1}}}))
    with pytest.raises(ValueError, match='C04'):
        scanner.module.main()
    assert _snapshot(scanner.module) == SENTINELS


def test_valid_zero_manifest_populations_are_not_missing_inputs(scanner):
    _write_shard(scanner.records / 'a.jsonl.gz', [{'mesh': []}])
    (scanner.root / 'unindexed-manifest.json').write_text(
        json.dumps({'files': {'fixture': {'cancer_text': 0}}}))
    (scanner.root / 'manifest-c04only.json').write_text(json.dumps({'files': {}}))
    scanner.module.main()
    result = json.loads(scanner.module.OUT_JSON.read_text())
    assert result['excluded_streams'] == {'text_matched_no_mesh': 0, 'c04_core': 0}


@pytest.mark.parametrize('excluded', [
    {'text_matched_no_mesh': -1}, {'text_matched_no_mesh': True},
    {'c04_core': 1},
])
def test_offline_impossible_population_preserves_all_outputs(scanner, monkeypatch, excluded):
    module = scanner.module
    raw = json.loads((REPO / 'analysis/atlas-site-coverage.json').read_text())
    raw['excluded_streams'] = excluded
    module.OUT_JSON.write_text(json.dumps(raw))
    before = _snapshot(module)
    monkeypatch.setattr(sys, 'argv', ['atlas_site_coverage.py', '--render-only'])
    with pytest.raises(ValueError):
        module.main()
    assert _snapshot(module) == before


def test_generic_descriptor_does_not_imply_site_free_literature(scanner):
    _write_shard(scanner.records / 'a.jsonl.gz', [
        {'mesh': ['Neoplasms', 'Bone Neoplasms']}, {'mesh': None},
    ])
    scanner.module.main()
    md = scanner.module.OUT_MD.read_text()
    assert 'no record here can lack MeSH' not in md
    assert 'Custom inputs may contain records without MeSH' in md
    assert 'They may also name sites outside this list' in md
    assert 'is the reading the original sentence described' not in md
