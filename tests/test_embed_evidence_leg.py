"""Tests for the #411 semantic-retrieval evidence-detection leg.

The embedding RECOMPUTE (scripts/embed_gold_set.py) needs sentence-transformers
and is not run in CI; these tests cover the pure k-NN / combine logic and the
committed artifacts (the embedding vectors + the eval result), so CI stays offline.
"""

import json
import sys
from pathlib import Path

import numpy as np
import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPTS_DIR = REPO_ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS_DIR))

import embed_evidence_leg as leg  # noqa: E402

NPZ = REPO_ROOT / "analysis" / "evidence-gold-embeddings.npz"
EVAL_JSON = REPO_ROOT / "analysis" / "evidence-gold-embedding-eval.json"


# --------------------------------------------------------------------------
# Pure logic
# --------------------------------------------------------------------------


def test_knn_evidence_pred_separates_clusters():
    # Two well-separated clusters: positives near [1,0], negatives near [0,1].
    pos = np.array([[1.0, 0.0], [0.99, 0.01], [0.98, 0.02], [0.97, 0.03]])
    neg = np.array([[0.0, 1.0], [0.01, 0.99], [0.02, 0.98], [0.03, 0.97]])
    emb = np.vstack([pos, neg])
    emb = emb / np.linalg.norm(emb, axis=1, keepdims=True)
    gold = np.array([1, 1, 1, 1, 0, 0, 0, 0])
    pred = leg.knn_evidence_pred(emb, gold, k=3)
    assert list(pred) == list(gold)  # LOO k-NN recovers the clusters


def test_knn_excludes_self():
    # A lone positive surrounded by negatives must NOT predict itself positive
    # (self is excluded), so its neighbors (all negative) drive it to 0.
    emb = np.array([[1.0, 0.0], [0.0, 1.0], [0.0, 1.0], [0.0, 1.0]])
    emb = emb / np.linalg.norm(emb, axis=1, keepdims=True)
    gold = np.array([1, 0, 0, 0])
    pred = leg.knn_evidence_pred(emb, gold, k=3)
    assert pred[0] == 0


def test_combine_logic():
    assert leg.combine(["preclinical-invivo", "", ""], [0, 1, 0]) == [
        "preclinical-invivo", "embedding-detected", ""
    ]
    # an existing leg prediction is never overwritten by the embedding
    assert leg.combine(["clinical-other"], [1]) == ["clinical-other"]


def test_load_gold_v2_is_balanced():
    rows = leg.load_gold_v2()
    assert len(rows) > 200
    pos = sum(1 for _p, g in rows if g not in ("", "none-applicable"))
    frac = pos / len(rows)
    assert 0.2 < frac < 0.45  # balanced (NOT the 87%-positive v1 set)


# --------------------------------------------------------------------------
# Committed artifacts
# --------------------------------------------------------------------------


def test_committed_embeddings_shape_and_normalized():
    data = np.load(NPZ, allow_pickle=True)
    emb = data["embeddings"]
    assert emb.shape[0] == len(data["pmids"]) > 200
    assert emb.shape[1] == 768  # SPECTER dimension
    norms = np.linalg.norm(emb.astype("float64"), axis=1)
    assert np.allclose(norms, 1.0, atol=1e-3)  # row-normalized for cosine


def test_committed_eval_recall_ladder():
    r = json.loads(EVAL_JSON.read_text())
    b, mh, me = r["baseline"], r["mesh"], r["mesh_plus_embedding"]
    # recall ladder: keyword <= +MeSH <= +MeSH+embedding
    assert b["recall"] <= mh["recall"] <= me["recall"]
    # the documented v1 recalls are reproduced (shared positives)
    assert b["recall"] == pytest.approx(0.552, abs=0.01)
    assert mh["recall"] == pytest.approx(0.678, abs=0.01)


def test_embedding_recovers_residue_without_precision_drop():
    """The load-bearing #411 finding: the semantic leg recovers a real chunk of the
    MeSH-residual false-negatives, and (on the balanced v2 set) does so WITHOUT a
    precision drop, because its own precision exceeds the keyword baseline's."""
    r = json.loads(EVAL_JSON.read_text())
    assert r["residue_total"] > 0
    assert r["residue_recovered"] >= 10  # recovers a substantial chunk of the residue
    mh, me = r["mesh"], r["mesh_plus_embedding"]
    assert me["recall"] > mh["recall"] + 0.1  # meaningful recall lift
    assert me["precision"] >= mh["precision"] - 0.02  # precision does not drop materially
