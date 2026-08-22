#!/usr/bin/env python3
"""How much precedent each reagent in the P1 protocol actually has.

`analysis/p1-wetlab-protocol.md` proposes the experiment this project calls its
single biggest credibility step: a GPX4-inhibitor by FSP1-inhibitor dose matrix
in an FSP1-low persister line. It names RSL3 (or ML162) against iFSP1 (or
brequinar), and presents the two arms symmetrically.

THE CENSUS CAN SAY WHETHER THEY ARE SYMMETRIC IN PRECEDENT, and they are not.
This matters for the thing the protocol is for: it exists to be handed to a
collaborator, and a collaborator sourcing reagents, choosing doses and
comparing results needs to know which half of the experiment has a literature
behind it and which half is close to first-in-field.

WHAT A MENTION COUNT MEASURES. How often the ferroptosis literature names a
compound in a title or abstract -- so it bounds how much comparative work
exists, which is the collaborator's question. It says nothing about whether a
compound is GOOD. A new and excellent inhibitor is rare in the literature for
the same reason a poor one is, and the fix for that is time rather than a
different count. Nor is a thin literature a reason not to run an experiment;
it is a reason to say so in the protocol.

Mentions are counted in title and abstract only, so a compound used in Methods
and never named in the abstract is undercounted. That biases every row the same
way and does not disturb the comparison between rows, which is what carries
here.
"""
import argparse
import gzip
import json
import re
from collections import Counter
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
RECORDS = REPO / "corpus/atlas/records"
PROTOCOL = REPO / "analysis/p1-wetlab-protocol.md"
OUT_MD = REPO / "analysis/census-protocol-precedent.md"
OUT_JSON = REPO / "analysis/census-protocol-precedent.json"

FERRO = "ferroptosis"
# Grouped by the role the protocol assigns, so the comparison is between the
# two ARMS of the experiment rather than between arbitrary compounds.
REAGENTS = {
    "GPX4 inhibitor (arm 1)": {
        "RSL3": r"\brsl3\b", "ML162": r"\bml162\b", "ML210": r"\bml210\b",
        "FIN56": r"\bfin56\b", "FINO2": r"\bfino2\b",
        "withaferin A": r"\bwithaferin a\b", "altretamine": r"\baltretamine\b",
    },
    "FSP1 / DHODH inhibitor (arm 2)": {
        "iFSP1": r"\bifsp1\b", "FSEN1": r"\bfsen1\b", "viFSP1": r"\bvifsp1\b",
        "icFSP1": r"\bicfsp1\b", "brequinar": r"\bbrequinar\b",
    },
    "system Xc- inhibitor (not in P1)": {
        "erastin": r"\berastin\b", "IKE": r"\b(ike|imidazole ketone erastin)\b",
        "sulfasalazine": r"\bsulfasalazine\b", "sorafenib": r"\bsorafenib\b",
    },
}
# The protocol names a persister system rather than a line. These are the
# lines the field actually uses, so a collaborator can see whether the
# system they have is one others have published in.
CELL_LINES = ["hepg2", "a549", "4t1", "hela", "mcf-7", "hct116", "huh7",
              "u87", "pc-3", "ht-1080", "panc-1", "b16", "pc9", "h1975",
              "hcc827", "sk-hep-1"]
# Below this a compound has essentially no comparative literature.
THIN = 50


def scan(stride: int = 1) -> dict:
    pats = {g: {k: re.compile(v) for k, v in d.items()}
            for g, d in REAGENTS.items()}
    line_pats = {k: re.compile(rf"\b{re.escape(k)}\b") for k in CELL_LINES}
    counts = {g: Counter() for g in REAGENTS}
    lines = Counter()
    n = 0
    for f in sorted(RECORDS.glob("*.jsonl.gz"))[::stride]:
        with gzip.open(f, "rt", encoding="utf-8") as fh:
            for line in fh:
                r = json.loads(line)
                ms = {m.lower() for m in (r.get("mesh") or [])}
                if FERRO not in ms:
                    continue
                n += 1
                blob = f"{r.get('title') or ''} {r.get('abstract') or ''}".lower()
                for g, ps in pats.items():
                    for k, p in ps.items():
                        if p.search(blob):
                            counts[g][k] += 1
                for k, p in line_pats.items():
                    if p.search(blob):
                        lines[k] += 1
    return {"ferroptosis_articles": n, "thin_threshold": THIN,
            "reagents": {g: dict(c) for g, c in counts.items()},
            "cell_lines": dict(lines)}


def _protocol_named() -> set:
    """Compounds the protocol actually names, read from it rather than typed."""
    txt = PROTOCOL.read_text(encoding="utf-8").lower()
    named = set()
    for group in REAGENTS.values():
        for k in group:
            if re.search(rf"\b{re.escape(k.lower())}\b", txt):
                named.add(k)
    return named


def assemble(d: dict) -> dict:
    named = _protocol_named()
    n = d["ferroptosis_articles"]
    groups = []
    for g, c in d["reagents"].items():
        rows = [{"reagent": k, "articles": v,
                 "share": round(100 * v / n, 2) if n else None,
                 "in_protocol": k in named,
                 "thin": v < d["thin_threshold"]}
                for k, v in sorted(c.items(), key=lambda kv: -kv[1])]
        groups.append({"role": g, "rows": rows,
                       "protocol_total": sum(r["articles"] for r in rows
                                             if r["in_protocol"]),
                       "group_total": sum(r["articles"] for r in rows)})
    out = dict(d)
    out["groups"] = groups
    out["protocol_named"] = sorted(named)
    # THE COMPARISON THAT MATTERS: the two arms of one experiment.
    arm1 = next(g for g in groups if g["role"].endswith("(arm 1)"))
    arm2 = next(g for g in groups if g["role"].endswith("(arm 2)"))
    out["arm1_precedent"] = arm1["protocol_total"]
    out["arm2_precedent"] = arm2["protocol_total"]
    out["arm_asymmetry"] = (round(arm1["protocol_total"] / arm2["protocol_total"], 1)
                            if arm2["protocol_total"] else None)
    out["thin_arm"] = (arm2["role"] if out["arm_asymmetry"]
                       and out["arm_asymmetry"] > 1 else arm1["role"])
    out["cell_line_rows"] = [
        {"line": k, "articles": v, "share": round(100 * v / n, 2)}
        for k, v in sorted(d["cell_lines"].items(), key=lambda kv: -kv[1])]
    return out


def render(d: dict) -> str:
    L = ["# What precedent the P1 protocol's reagents have\n"]
    L.append(
        f"Generated by `scripts/census_protocol_precedent.py` over the "
        f"{d['ferroptosis_articles']:,} census articles carrying the "
        f"`Ferroptosis` descriptor. `analysis/p1-wetlab-protocol.md` proposes "
        f"the experiment this project calls its single biggest credibility "
        f"step, and names its two arms symmetrically. They are not symmetric "
        f"in precedent.\n"
    )
    L.append(
        f"**Arm 1 ({d['arm1_precedent']:,} articles) against arm 2 "
        f"({d['arm2_precedent']:,}) — a factor of {d['arm_asymmetry']}.** The "
        f"thin side is the *{d['thin_arm']}*, and the protocol presents it as "
        f"an equal partner.\n"
    )
    for g in d["groups"]:
        L.append(f"### {g['role']}\n")
        L.append("| reagent | articles | share of ferroptosis literature | in P1 |")
        L.append("|---|--:|--:|---|")
        for r in g["rows"]:
            mark = "**yes**" if r["in_protocol"] else "—"
            thin = " *" if r["thin"] and r["in_protocol"] else ""
            L.append(f"| {r['reagent']}{thin} | {r['articles']:,} | "
                     f"{r['share']}% | {mark} |")
        L.append("")
    L.append(
        f"\\* named in the protocol and mentioned in fewer than "
        f"{d['thin_threshold']} articles.\n"
    )
    L.append("## What this means for the protocol\n")
    L.append(
        "It does not mean the experiment is wrong. FSP1 inhibitors are recent "
        "-- iFSP1 comes from the 2019 papers this project builds on -- and a "
        "compound is rare in the literature when it is new for the same reason "
        "it is rare when it is poor. A count cannot tell those apart, and the "
        "fix is time rather than a different measurement.\n"
    )
    L.append(
        "What it means is that a collaborator handed this protocol will find "
        "abundant precedent for choosing an RSL3 dose and almost none for "
        "choosing an iFSP1 one, and the protocol should say so where it names "
        "the reagents rather than leaving it to be discovered at the bench. "
        "The same applies to the surrogates it offers: ML162 is presented as a "
        "more stable stand-in for RSL3, and it appears in a fiftieth as many "
        "articles, so a result obtained with it has correspondingly less to be "
        "compared against.\n"
    )
    L.append("## Cell lines the field uses\n")
    L.append(
        "The protocol specifies a persister SYSTEM rather than a line, which "
        "is the right level of specification for a state-dependent hypothesis. "
        "These are the lines the ferroptosis literature actually names, so a "
        "collaborator can see whether the system available to them is one "
        "others have published in.\n"
    )
    L.append("| cell line | articles | share |")
    L.append("|---|--:|--:|")
    for r in d["cell_line_rows"][:10]:
        L.append(f"| {r['line']} | {r['articles']:,} | {r['share']}% |")
    L.append("")
    L.append("## What a mention count is not\n")
    L.append(
        "Evidence about quality, and not a reason to avoid a reagent. It "
        "bounds how much COMPARATIVE work exists, which is the collaborator's "
        "question when setting a dose or interpreting a result. Mentions are "
        "counted in title and abstract only, so a compound used in Methods and "
        "never named in the abstract is undercounted -- which biases every row "
        "the same way and leaves the comparison between rows, the thing that "
        "carries here, intact.\n"
    )
    return "\n".join(L)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--render-only", action="store_true")
    a = ap.parse_args()
    d = assemble(json.loads(OUT_JSON.read_text()) if a.render_only
                 else scan(a.stride))
    OUT_JSON.write_text(json.dumps(d, indent=1) + "\n")
    OUT_MD.write_text(render(d))
    print(f"wrote {OUT_MD}")
    print(f"  arm1 {d['arm1_precedent']:,} vs arm2 {d['arm2_precedent']:,} "
          f"= {d['arm_asymmetry']}x, thin side: {d['thin_arm']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
