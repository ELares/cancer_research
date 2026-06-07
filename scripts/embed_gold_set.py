#!/usr/bin/env python3
"""Precompute SPECTER embeddings for the evidence gold set (#411 / #346 leg).

#346's MeSH-descriptor fallback lifted gold-set evidence-detection recall from
55.2% to 67.8%, but hit a hard floor: of the 39 baseline false-negatives, 9 carry
empty `mesh_terms` and more carry non-discriminative MeSH, so MeSH cannot reach
that residue. #411 adds a semantic-retrieval leg to recover part of it.

This is the RECOMPUTE step (option B dependency story): it needs a biomedical
document-embedding model (SPECTER) and is run LOCALLY, not in CI. It embeds each
gold-set record's title+abstract and writes the embedding VECTORS to
`analysis/evidence-gold-embeddings.npz` (committed, ~1 MB). The scoring/evaluation
(`scripts/embed_evidence_leg.py`) and CI then read those committed vectors with
numpy only, so the offline/pinned-deps CI contract is preserved (sentence-
transformers / torch are never added to requirements-lock.txt).

Run (in a venv with sentence-transformers installed):
    pip install sentence-transformers
    python scripts/embed_gold_set.py
Model: allenai-specter (768-d, designed to embed scientific title+abstract for
document similarity).
"""

import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPTS_DIR = REPO_ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS_DIR))

import retag_gold_set as rg  # noqa: E402  (corpus paths, MeSH harness)
import tag_articles  # noqa: E402
from article_io import load_article  # noqa: E402
from embed_evidence_leg import GOLD_V2, load_gold_v2  # noqa: E402  (shared v2 loader)


OUT_NPZ = REPO_ROOT / "analysis" / "evidence-gold-embeddings.npz"
MODEL_NAME = "allenai-specter"


def gold_texts():
    """(pmids, texts) for the v2 gold records present in the corpus. v2 (32%
    positive, 183 negatives) is balanced enough to measure a detection leg's
    precision, unlike v1 (87% positive). Text is title+abstract (no MeSH), so the
    embedding never sees the controlled vocabulary the keyword/MeSH legs use."""
    pmids, texts = [], []
    for pmid, _gold in load_gold_v2():
        path = rg.PMID_DIR / f"{pmid}.md"
        if not path.exists():
            continue
        fm, body = load_article(path)
        # include_metadata=False => title + abstract only (the #412 leakage-free path)
        texts.append(tag_articles.get_searchable_text(fm, body, include_metadata=False))
        pmids.append(pmid)
    return pmids, texts


def main() -> int:
    try:
        from sentence_transformers import SentenceTransformer
    except ImportError:
        raise SystemExit(
            "sentence-transformers is required for the RECOMPUTE step only "
            "(not for CI / the offline eval). `pip install sentence-transformers` in a venv."
        )
    pmids, texts = gold_texts()
    print(f"embedding {len(pmids)} gold records with {MODEL_NAME} ...")
    model = SentenceTransformer(MODEL_NAME)
    emb = model.encode(texts, convert_to_numpy=True, normalize_embeddings=True, batch_size=32)
    np.savez_compressed(
        OUT_NPZ,
        pmids=np.array(pmids),
        embeddings=emb.astype("float32"),
        model=np.array(MODEL_NAME),
    )
    print(f"wrote {OUT_NPZ.relative_to(REPO_ROOT)}  shape={emb.shape}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
