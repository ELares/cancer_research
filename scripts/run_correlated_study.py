#!/usr/bin/env python3
"""Preflight archive preservation before invoking the frozen covariance study.

The original driver remains part of the experiment's hashed implementation.
This operational entry point validates all existing records before allowing
that driver to fill a missing archive slot. It does not change any numerical
calculation, experimental input, or replay rule.
"""

import argparse
import gzip
import json

import proposal_correlated_study as study


def current_runtime():
    return {
        'python': study.platform.python_version(),
        'numpy': study.np.__version__,
        'scipy': study.scipy.__version__,
        'bit_generator': type(study.np.random.default_rng(0).bit_generator).__name__,
    }


def preflight(*, resume=False):
    """Reject incompatible or invalid existing records before any new run."""
    planned = {
        study.archive_path(study.ARCHIVES, fixture, seed): (fixture, seed)
        for fixture, seeds in study.LEARNED_SEEDS.items() for seed in seeds
    }
    found = set(study.ARCHIVES.glob('*.json.gz'))
    extras = found - planned.keys()
    if extras:
        names = ', '.join(sorted(path.name for path in extras))
        raise ValueError(f'Unexpected study archives: {names}; inspect the archive set before running')
    existing = [path for path in planned if path.exists()]
    if existing and not resume:
        raise FileExistsError(
            f'{existing[0].name} already exists; use --resume to validate and retain existing runs')

    commit = study.committed_sources()
    runtime, archived_commit = None, None
    for path in existing:
        raw = json.loads(gzip.decompress(path.read_bytes()))
        study._equal([raw['fixture'], raw['seed']], list(planned[path]), 'existing archive identity')
        study.assess_archive(raw)
        if runtime is None:
            runtime, archived_commit = raw['runtime'], raw['implementation_commit']
        study._equal(raw['runtime'], runtime, 'existing archive runtimes')
        study._equal(raw['implementation_commit'], archived_commit, 'existing implementation commits')

    if existing and len(existing) < len(planned):
        if current_runtime() != runtime:
            raise ValueError(
                'Current runtime differs from existing study archives; resume with their recorded '
                f'Python/NumPy/SciPy/bit-generator environment: {runtime}. '
                'No new runs were started; preserve the existing records.')
        if commit != archived_commit:
            raise ValueError(
                'Current implementation commit differs from existing study archives; '
                f'resume from their recorded source commit {archived_commit}. '
                'No new runs were started; preserve the existing records.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--render-only', action='store_true')
    mode.add_argument('--resume', action='store_true')
    args = parser.parse_args()
    if not args.render_only:
        preflight(resume=args.resume)
    return study.main()


if __name__ == '__main__':
    raise SystemExit(main())
