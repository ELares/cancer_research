"""Tag-level precision/recall cannot sign bias in an article-level co-tag rate.

The former guards required an unjustified direction from precision < recall.
These counterexamples preserve even per-label error totals while giving
opposite biases. The manuscript must retain this distinction at both sites.
"""
import json
import re
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
MANUSCRIPT = REPO / 'article/drafts/v1.md'


def _rate(rows):
    return sum(len(row) >= 2 for row in rows) / sum(bool(row) for row in rows)


def _errors(truth, predicted, label):
    return (
        sum(label in t and label in p for t, p in zip(truth, predicted)),
        sum(label not in t and label in p for t, p in zip(truth, predicted)),
        sum(label in t and label not in p for t, p in zip(truth, predicted)),
    )


def _section(number):
    text = MANUSCRIPT.read_text()
    return re.search(rf'^### {re.escape(number)} .*?\n(.*?)(?=^### |\Z)',
                     text, re.M | re.S)[1]


def test_identical_per_label_precision_and_recall_allow_opposite_rate_bias():
    truth = list(map(set, ['AB', 'AB', 'AB', 'A', 'B', '', '', '', '']))
    inflated = list(map(set, ['AB', 'AB', 'AB', '', '', 'AB', 'AB', '', '']))
    deflated = list(map(set, ['B', 'A', 'AB', 'A', 'B', 'A', 'A', 'B', 'B']))
    for label in 'AB':
        assert _errors(truth, inflated, label) == (3, 2, 1)
        assert _errors(truth, deflated, label) == (3, 2, 1)
    precision, recall = 3 / (3 + 2), 3 / (3 + 1)
    assert precision < recall
    assert _rate(deflated) == pytest.approx(1 / 9)
    assert _rate(truth) == pytest.approx(3 / 5)
    assert _rate(inflated) == 1


def test_reported_precision_still_names_its_source_measurement():
    source = (REPO / 'analysis/mechanism-precision-report.md').read_text()
    precision = re.search(r'\*\*Overall strict precision: \d+/\d+ = ([\d.]+)%\*\*', source)[1]
    assert f'{precision}% strict' in _section('3.6')
    recall = json.loads((REPO / 'analysis/mechanism-recall.json').read_text())
    assert 0 < recall['aggregates']['volume_weighted_recall_leakage_free'] <= 1


def test_limitation_requires_joint_error_information_not_a_signed_bias():
    text = ' '.join(_section('3.6').split())
    assert 'Aggregate precision and recall do not determine the direction of bias' in text
    assert 'Article-level joint error information is needed' in text
    assert 'tagged-article denominator' in text
    assert 'more likely inflated than deflated' not in text
    assert 'precision error is amplified and the recall error is not' not in text


def test_rate_comparison_does_not_claim_a_causal_decomposition():
    text = ' '.join(_section('3.13').split())
    assert 'do not identify a causal decomposition' in text
    assert 'not independent' in text
    assert 'generous reading' not in text
    assert 'entirely explained' not in text
    assert 'full-text keyword matcher' not in text
