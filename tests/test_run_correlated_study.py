"""Operational preflight uses toy records and never runs a study pilot."""

import copy
import gzip
import json
from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import run_correlated_study as runner  # noqa: E402


@pytest.fixture
def environment(monkeypatch, tmp_path):
    study = runner.study
    archive_dir = tmp_path / 'archives'
    archive_dir.mkdir()
    monkeypatch.setattr(study, 'ARCHIVES', archive_dir)
    monkeypatch.setattr(study, 'LEARNED_SEEDS', {'toy': (41101, 41102, 41103)})
    monkeypatch.setattr(study, 'OUT_JSON', tmp_path / 'report.json')
    monkeypatch.setattr(study, 'OUT_MD', tmp_path / 'report.md')
    study.OUT_JSON.write_text('original JSON')
    study.OUT_MD.write_text('original Markdown')
    runtime = {'python': '3.14.6', 'numpy': '2.4.0', 'scipy': '1.17.0', 'bit_generator': 'PCG64'}
    commit = 'a' * 40
    monkeypatch.setattr(study, 'committed_sources', lambda: commit)
    monkeypatch.setattr(runner, 'current_runtime', lambda: copy.deepcopy(runtime))
    assessed, delegated = [], []

    def assess(raw):
        assessed.append(raw['seed'])
        study._equal(raw['specification'], {'revision': 'toy'}, 'archive specification')
        if not raw['valid']:
            raise ValueError('corrupt toy archive')

    def delegate():
        delegated.append(list(sys.argv[1:]))
        return 17

    def forbidden(*args, **kwargs):
        pytest.fail('preflight tests must never generate a planned run or train a pilot')

    monkeypatch.setattr(study, 'assess_archive', assess)
    monkeypatch.setattr(study, 'main', delegate)
    monkeypatch.setattr(study, 'run_learned', forbidden)
    monkeypatch.setattr(study.strategy, 'train_pilot', forbidden)

    def write(seed, **changes):
        raw = {'fixture': 'toy', 'seed': seed, 'specification': {'revision': 'toy'},
               'runtime': copy.deepcopy(runtime), 'implementation_commit': commit, 'valid': True}
        raw.update(changes)
        path = study.archive_path(archive_dir, 'toy', seed)
        path.write_bytes(gzip.compress(json.dumps(raw).encode(), mtime=0))
        return path

    def invoke(*args):
        monkeypatch.setattr(sys, 'argv', ['run_correlated_study.py', *args])
        return runner.main()

    def snapshot():
        return {path.relative_to(tmp_path): path.read_bytes()
                for path in tmp_path.rglob('*') if path.is_file()}

    return write, invoke, assessed, delegated, runtime, snapshot


def test_missing_first_run_does_not_start_before_corrupt_later_archive(environment):
    write, invoke, assessed, delegated, _, snapshot = environment
    write(41102, valid=False)
    before = snapshot()
    with pytest.raises(ValueError, match='corrupt toy archive'):
        invoke('--resume')
    assert assessed == [41102]
    assert delegated == []
    assert snapshot() == before


@pytest.mark.parametrize('fault', ['identity', 'specification', 'runtime', 'implementation_commit'])
def test_all_existing_records_are_validated_before_delegation(environment, fault):
    write, invoke, _, delegated, runtime, snapshot = environment
    write(41102)
    changes = {
        'identity': {'fixture': 'another'},
        'specification': {'specification': {'revision': 'other'}},
        'runtime': {'runtime': {**runtime, 'numpy': 'other'}},
        'implementation_commit': {'implementation_commit': 'b' * 40},
    }
    write(41103, **changes[fault])
    before = snapshot()
    with pytest.raises(ValueError, match='identity|specification|runtimes|implementation commits'):
        invoke('--resume')
    assert delegated == []
    assert snapshot() == before


@pytest.mark.parametrize('field', ['python', 'numpy', 'scipy', 'bit_generator'])
def test_partial_resume_requires_matching_current_runtime(environment, monkeypatch, field):
    write, invoke, assessed, delegated, runtime, snapshot = environment
    write(41101)
    monkeypatch.setattr(runner, 'current_runtime', lambda: {**runtime, field: 'other'})
    before = snapshot()
    with pytest.raises(ValueError, match='Current runtime differs'):
        invoke('--resume')
    assert assessed == [41101]
    assert delegated == []
    assert snapshot() == before


def test_partial_resume_requires_matching_current_implementation_commit(environment, monkeypatch):
    write, invoke, _, delegated, _, snapshot = environment
    write(41101)
    monkeypatch.setattr(runner.study, 'committed_sources', lambda: 'b' * 40)
    before = snapshot()
    with pytest.raises(ValueError, match='Current implementation commit differs'):
        invoke('--resume')
    assert delegated == []
    assert snapshot() == before


def test_complete_archives_allow_a_different_current_runtime_and_commit(environment, monkeypatch):
    write, invoke, assessed, delegated, _, _ = environment
    for seed in (41103, 41101, 41102):
        write(seed)
    monkeypatch.setattr(runner, 'current_runtime', lambda: pytest.fail('no new-run runtime required'))
    monkeypatch.setattr(runner.study, 'committed_sources', lambda: 'b' * 40)
    assert invoke('--resume') == 17
    assert assessed == [41101, 41102, 41103]
    assert delegated == [['--resume']]


def test_default_rejects_existing_later_seed_before_any_new_run(environment):
    write, invoke, assessed, delegated, _, snapshot = environment
    write(41103)
    before = snapshot()
    with pytest.raises(FileExistsError, match='use --resume'):
        invoke()
    assert assessed == delegated == []
    assert snapshot() == before


@pytest.mark.parametrize('args', [(), ('--resume',)])
def test_unexpected_archives_are_rejected_before_any_new_run(environment, args):
    _, invoke, assessed, delegated, _, snapshot = environment
    (runner.study.ARCHIVES / 'unexpected.json.gz').write_bytes(b'unexpected')
    before = snapshot()
    with pytest.raises(ValueError, match='Unexpected study archives'):
        invoke(*args)
    assert assessed == delegated == []
    assert snapshot() == before


def test_valid_partial_resume_validates_all_existing_records_and_delegates(environment):
    write, invoke, assessed, delegated, _, snapshot = environment
    write(41103)
    write(41102)
    before = snapshot()
    assert invoke('--resume') == 17
    assert assessed == [41102, 41103]
    assert delegated == [['--resume']]
    assert snapshot() == before


@pytest.mark.parametrize('args', [(), ('--resume',)])
def test_empty_archive_set_can_start_with_no_existing_provenance(environment, args):
    _, invoke, assessed, delegated, _, _ = environment
    assert invoke(*args) == 17
    assert assessed == []
    assert delegated == [list(args)]


def test_render_only_delegates_without_production_preflight(environment, monkeypatch):
    _, invoke, assessed, delegated, _, _ = environment
    monkeypatch.setattr(runner, 'preflight', lambda **kwargs: pytest.fail('render-only cannot start runs'))
    assert invoke('--render-only') == 17
    assert assessed == []
    assert delegated == [['--render-only']]


@pytest.mark.parametrize('args', [(), ('--resume',)])
def test_source_freeze_failure_prevents_delegation(environment, monkeypatch, args):
    _, invoke, _, delegated, _, snapshot = environment

    def uncommitted():
        raise ValueError('commit the protocol and numerical sources before evaluation')

    monkeypatch.setattr(runner.study, 'committed_sources', uncommitted)
    before = snapshot()
    with pytest.raises(ValueError, match='commit the protocol'):
        invoke(*args)
    assert delegated == []
    assert snapshot() == before


def test_modes_remain_mutually_exclusive(environment):
    _, invoke, assessed, delegated, _, _ = environment
    with pytest.raises(SystemExit) as exc:
        invoke('--resume', '--render-only')
    assert exc.value.code == 2
    assert assessed == delegated == []
