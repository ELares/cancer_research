#!/usr/bin/env python3
"""Measure the v2 evidence tagger against the v1 baseline (#TAGGER-V2).

Runs the REAL production path (`get_searchable_text` -> `match_evidence_level`)
twice per evaluation set, once with `FERRO_EVIDENCE_V2` off and once on, and
writes `analysis/evidence-v2-eval.md`.

Three evaluation sets, because the honest answer depends on which labels you
trust:

  DEV        100 records carrying the v1 human labels. The v2 keyword lists were
             developed against these, so this number is optimistic by
             construction and is reported only for completeness.
  HELDOUT    170 records with NO human label, never inspected during
             development. Labels come from the independent full-text relabeling
             pass in analysis/evidence-gold-set-v3-fulltext.csv. This is the
             generalization test.
  CONSENSUS  Records where the v1 human label and the independent relabel AGREE
             (77 of 100). The most defensible labels available, and the
             conservative headline.

Offline and deterministic. Re-exec's itself in a subprocess per configuration
because the feature flags are read at import time.

Usage:
    python scripts/evaluate_evidence_v2.py
"""

import collections
import csv
import json
import os
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from config import PROJECT_ROOT

HUMAN_LABELS = PROJECT_ROOT / "analysis" / "evidence-gold-labels-v1.csv"
FULLTEXT_LABELS = PROJECT_ROOT / "analysis" / "evidence-gold-set-v3-fulltext.csv"
PMID_DIR = PROJECT_ROOT / "corpus" / "by-pmid"
OUTPUT = PROJECT_ROOT / "analysis" / "evidence-v2-eval.md"

CLINICAL = {"phase3-clinical", "phase2-clinical", "phase1-clinical", "clinical-other"}
LEVELS = ["phase3-clinical", "phase2-clinical", "phase1-clinical", "clinical-other",
          "preclinical-invivo", "preclinical-invitro", "theoretical", "none-applicable"]

_WORKER = "--_predict-worker"


def _predict_worker(pmid_file: str) -> None:
    """Child process: predict every PMID under the ambient flag settings."""
    import tag_articles
    from article_io import load_article

    out = {}
    for pmid in json.load(open(pmid_file)):
        path = PMID_DIR / f"{pmid}.md"
        if not path.exists():
            continue
        fm, body = load_article(path)
        out[pmid] = tag_articles.match_evidence_level(fm, tag_articles.get_evidence_text(fm, body))
    sys.stdout.write("@@JSON@@" + json.dumps(out))


def predict(pmids, env_extra: dict, tmp: Path) -> dict:
    tmp.write_text(json.dumps(list(pmids)), encoding="utf-8")
    env = dict(os.environ)
    env.update(env_extra)
    proc = subprocess.run([sys.executable, str(Path(__file__).resolve()), _WORKER, str(tmp)],
                          capture_output=True, text=True, env=env)
    if proc.returncode != 0:
        raise SystemExit(f"prediction worker failed:\n{proc.stderr[-2000:]}")
    marker = proc.stdout.find("@@JSON@@")
    return json.loads(proc.stdout[marker + len("@@JSON@@"):])


def norm(pred: str) -> str:
    """Repo convention: an empty prediction is the explicit none-applicable class."""
    return pred or "none-applicable"


def wilson(k: int, n: int, z: float = 1.96):
    if n == 0:
        return 0.0, 0.0
    p = k / n
    denom = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / denom
    half = z * ((p * (1 - p) / n + z * z / (4 * n * n)) ** 0.5) / denom
    return max(0.0, centre - half), min(1.0, centre + half)


def metrics(pairs):
    n = len(pairs)
    exact = sum(1 for g, p in pairs if norm(p) == g)
    tp = sum(1 for g, p in pairs if g != "none-applicable" and p)
    fp = sum(1 for g, p in pairs if g == "none-applicable" and p)
    fn = sum(1 for g, p in pairs if g != "none-applicable" and not p)
    prec = tp / (tp + fp) if tp + fp else 0.0
    rec = tp / (tp + fn) if tp + fn else 0.0
    f1 = 2 * prec * rec / (prec + rec) if prec + rec else 0.0

    def coarse(x):
        if x == "none-applicable":
            return "none"
        return "clinical" if x in CLINICAL else "preclinical"

    co = sum(1 for g, p in pairs if coarse(g) == coarse(norm(p)))
    return dict(n=n, exact_k=exact, exact=exact / n, coarse=co / n,
                prec=prec, rec=rec, f1=f1)


def load_label_sets():
    human = {r["pmid"]: r["gold_evidence_level"] for r in csv.DictReader(open(HUMAN_LABELS))}
    with open(FULLTEXT_LABELS, encoding="utf-8") as fh:
        rows = list(csv.DictReader(l for l in fh if not l.startswith("#")))
    fulltext = {r["pmid"]: r["tier_fulltext"] for r in rows}
    dev = dict(human)
    heldout = {p: t for p, t in fulltext.items() if p not in human}
    consensus = {p: human[p] for p in human if fulltext.get(p) == human[p]}
    return dev, heldout, consensus, human, fulltext


def main() -> None:
    dev, heldout, consensus, human, fulltext = load_label_sets()
    tmp = PROJECT_ROOT / "analysis" / ".evidence_v2_pmids.json"
    lines = ["# Evidence tagger v2 evaluation (#TAGGER-V2)", "",
             "Generated by `scripts/evaluate_evidence_v2.py`. Offline and deterministic.", "",
             "`v2` = `FERRO_EVIDENCE_V2=1`: section-scoped full text (Methods/Results only),",
             "expanded evidence vocabularies, in-vivo reagent guard, opinion-pub-type veto,",
             "title-shape review guard, the theoretical-dominance rule, and a dedicated "
             "evidence prose channel. Off by default, so",
             "the frozen corpus and every manuscript number are unchanged.", ""]

    n_both = sum(1 for p in human if p in fulltext)
    n_agree = sum(1 for p in human if fulltext.get(p) == human[p])
    lines += [f"Annotator agreement between the v1 human labels and the independent full-text",
              f"relabel: **{n_agree}/{n_both} = {n_agree / n_both:.1%}**. Neither column is ground truth;",
              "treat this as the practical ceiling for the task.", ""]

    try:
        for name, gold, note in [
            ("DEV (v1 human labels; v2 was TUNED on these)", dev,
             "Optimistic by construction. Reported for completeness only."),
            ("HELDOUT (independent relabel, never inspected during development)", heldout,
             "The generalization test."),
            ("CONSENSUS (human and independent relabel agree)", consensus,
             "The conservative headline."),
        ]:
            res, preds = {}, {}
            for cfg, env in [("v1 baseline", {}), ("v2", {"FERRO_EVIDENCE_V2": "1"})]:
                preds[cfg] = predict(gold.keys(), env, tmp)
                res[cfg] = metrics([(gold[p], preds[cfg][p]) for p in preds[cfg]])

            b, v = res["v1 baseline"], res["v2"]
            lines += [f"## {name}", "", f"n = {b['n']}. {note}", "",
                      "| config | exact-label | 95% CI | coarse | precision | recall | binary F1 |",
                      "|---|---|---|---|---|---|---|"]
            for cfg in ("v1 baseline", "v2"):
                m = res[cfg]
                lo, hi = wilson(m["exact_k"], m["n"])
                lines.append(f"| {cfg} | {m['exact']:.1%} | [{lo:.1%}, {hi:.1%}] | {m['coarse']:.1%} | "
                             f"{m['prec']:.1%} | {m['rec']:.1%} | {m['f1']:.3f} |")
            eb, ev = 1 - b["exact"], 1 - v["exact"]
            lines += ["", f"- exact-label error **{eb:.1%} -> {ev:.1%}** "
                          f"(**{eb / ev:.2f}x** reduction)" if ev else "",
                      f"- binary F1 {b['f1']:.3f} -> {v['f1']:.3f}", ""]

            lines += ["| tier | n | recall v1 | recall v2 |", "|---|---|---|---|"]
            bp, vp = preds["v1 baseline"], preds["v2"]
            for lv in LEVELS:
                idx = [p for p in vp if gold[p] == lv]
                if not idx:
                    continue
                rb = sum(1 for p in idx if norm(bp[p]) == lv) / len(idx)
                rv = sum(1 for p in idx if norm(vp[p]) == lv) / len(idx)
                lines.append(f"| {lv} | {len(idx)} | {rb:.1%} | {rv:.1%} |")
            conf = collections.Counter((gold[p], norm(vp[p])) for p in vp if norm(vp[p]) != gold[p])
            if conf:
                lines += ["", "Remaining v2 errors (gold -> predicted):", ""]
                lines += [f"- {c} x `{g}` -> `{p}`" for (g, p), c in conf.most_common(8)]
            lines.append("")
    finally:
        tmp.unlink(missing_ok=True)

    OUTPUT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"wrote {OUTPUT}")
    print("\n".join(l for l in lines if l.startswith(("|", "- exact", "## ", "Annotator"))))


if __name__ == "__main__":
    if len(sys.argv) > 2 and sys.argv[1] == _WORKER:
        _predict_worker(sys.argv[2])
    else:
        main()
