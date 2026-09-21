#!/usr/bin/env python3
"""Which synergy metric does the ferroptosis literature name?

`analysis/p1-wetlab-protocol.md` states its falsification threshold in
Chou-Talalay combination-index terms: refute if CI > 0.8. A pre-stated
threshold is only as useful as the comparison it enables -- a collaborator
who measures CI = 0.74 wants to know what the published CIs look like. The
protocol also names Bliss independence, which uses a different reference model.

SO THE QUESTION IS NOT WHETHER THE PROTOCOL IS RIGOROUS. It is whether its
number will be comparable to anything, and which metric the field names when
it names one.

THE SUBJECT ARM IS TOO THIN TO ANSWER THE SECOND HALF, and that is reported
rather than papered over: 8 combination papers name Chou-Talalay and 0 name
anything else, which is not a ranking. So the same rules run over a strided
sample of the WHOLE census as a control, where the counts are large enough to
separate the metrics -- and the ranking is read from there, with the thin
subject arm shown beside it as agreement rather than as evidence.

WHAT THIS INSTRUMENT CAN AND CANNOT DO, and the limit is severe enough to set
the framing. A synergy metric may be computed in a figure and quoted in
Results without appearing in the abstract. What this search measures is
which metric names appear in titles and abstracts, not computation rates.

The relative counts order observed mentions, not underlying computation rates.
Different metric names and bare claim words may have different probabilities
of appearing in abstracts. Those reporting rates are not measured here, so the
claim-to-metric ratio does not bound the corresponding gap in full papers.
"""
import argparse
import json
import re
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from census_input import atlas_root, iter_census_records  # noqa: E402

REPO = Path(__file__).resolve().parent.parent
RECORDS = atlas_root() / "records"
PROTOCOL = REPO / "analysis/p1-wetlab-protocol.md"
OUT_MD = REPO / "analysis/census-synergy-metrics.md"
OUT_JSON = REPO / "analysis/census-synergy-metrics.json"

FERRO = "ferroptosis"

# Named metrics: each is a defined quantity on a stated scale, so two papers
# reporting one can be compared. Grouped by metric, not by spelling.
METRICS = {
    "Chou-Talalay combination index": (
        r"\bchou[- ]talalay\b|\bcombination index\b|\bcompusyn\b|\bcalcusyn\b"),
    "Bliss independence": r"\bbliss\b",
    "Loewe additivity": r"\bloewe\b",
    "HSA (highest single agent)": r"\bhighest single agent\b|\bhsa (model|score|synergy)\b",
    "ZIP (zero interaction potency)": r"\bzero interaction potency\b|\bzip (model|score|synergy)\b",
    "isobologram": r"\bisobologram?\b|\bisobolographic\b",
    "dose reduction index": r"\bdose[- ]reduction index\b",
}
# Claim words that assert an interaction without naming a scale.
CLAIM = r"\bsynergis(?:m|tic|tically)\b|\bsynergy\b|\bsupra-?additive\b|\bpotentiat"
# The combination subset: a paper not about a combination has no occasion to
# report a synergy metric, so a share over all ferroptosis articles would be
# diluted by a denominator that was never eligible.
COMBO = (r"\bcombination therapy\b|\bcombined with\b|\bin combination\b|"
         r"\bco-?treatment\b|\bcombination treatment\b|\bdual (?:inhibit|target)")
# The control arm samples every 40th shard. Shards are CHRONOLOGICAL, so a
# prefix would sample one era; a stride spreads the sample across the range.
CONTROL_STRIDE = 40
# Below this the subject arm cannot separate two metrics and says so.
RANKABLE_MIN = 30


def scan_control() -> dict:
    """The same rules over a strided census sample, where counts are large
    enough to rank the metrics the subject arm cannot separate."""
    mp = {k: re.compile(v) for k, v in METRICS.items()}
    c = Counter()
    n = 0
    for r in iter_census_records(RECORDS, CONTROL_STRIDE):
        n += 1
        blob = f"{r.get('title') or ''} {r.get('abstract') or ''}".lower()
        for k, p in mp.items():
            if p.search(blob):
                c[k] += 1
    return {"control_records": n, "control_stride": CONTROL_STRIDE,
            "control_metric": dict(c)}


def scan(stride: int = 1) -> dict:
    mp = {k: re.compile(v) for k, v in METRICS.items()}
    claim = re.compile(CLAIM)
    combo = re.compile(COMBO)
    metric_all, metric_combo = Counter(), Counter()
    n = n_combo = n_claim = n_claim_combo = 0
    any_metric = any_metric_combo = 0
    claim_only_combo = 0
    for r in iter_census_records(RECORDS, stride):
        if FERRO not in {m.lower() for m in (r.get("mesh") or [])}:
            continue
        n += 1
        blob = f"{r.get('title') or ''} {r.get('abstract') or ''}".lower()
        is_combo = bool(combo.search(blob))
        has_claim = bool(claim.search(blob))
        hit = False
        if is_combo:
            n_combo += 1
        if has_claim:
            n_claim += 1
            if is_combo:
                n_claim_combo += 1
        for k, p in mp.items():
            if p.search(blob):
                hit = True
                metric_all[k] += 1
                if is_combo:
                    metric_combo[k] += 1
        if hit:
            any_metric += 1
            if is_combo:
                any_metric_combo += 1
        elif is_combo and has_claim:
            claim_only_combo += 1
    return {"ferroptosis_articles": n, "combination_articles": n_combo,
            "claim_articles": n_claim, "claim_combination": n_claim_combo,
            "any_metric": any_metric, "any_metric_combination": any_metric_combo,
            "claim_only_combination": claim_only_combo,
            "metric_all": dict(metric_all), "metric_combination": dict(metric_combo),
            **scan_control()}


def _protocol_metrics() -> list:
    """Which metrics the protocol names, read from it rather than typed."""
    txt = PROTOCOL.read_text(encoding="utf-8").lower()
    return [k for k, v in METRICS.items() if re.search(v, txt)]


def assemble(d: dict) -> dict:
    out = dict(d)
    named = _protocol_metrics()
    nc = d["combination_articles"]
    rows = [{"metric": k,
             "articles": d["metric_all"].get(k, 0),
             "combination_articles": d["metric_combination"].get(k, 0),
             "share_of_combination": (round(100 * d["metric_combination"].get(k, 0) / nc, 2)
                                      if nc else None),
             "in_protocol": k in named}
            for k in METRICS]
    for r in rows:
        r["control_articles"] = d["control_metric"].get(r["metric"], 0)
    # Ranked by the CONTROL, which is the arm with the counts to rank with.
    rows.sort(key=lambda r: (-r["control_articles"], -r["combination_articles"]))
    out["rows"] = rows
    out["protocol_named"] = named
    # THE SUBJECT ARM'S POWER, derived rather than assumed: with this few
    # named-metric articles it cannot separate two metrics, so the ranking is
    # read off the control and the subject arm is shown as agreement.
    out["subject_rankable"] = d["any_metric"] >= RANKABLE_MIN
    out["rankable_min"] = RANKABLE_MIN
    out["leading_metric"] = (rows[0]["metric"] if rows and
                             rows[0]["control_articles"] > rows[1]["control_articles"]
                             else None)
    out["leading_control_articles"] = rows[0]["control_articles"] if rows else 0
    out["runner_up_control_articles"] = rows[1]["control_articles"] if len(rows) > 1 else 0
    out["protocol_names_the_leader"] = bool(out["leading_metric"] and rows[0]["in_protocol"])
    # Agreement between the arms is a weaker claim than a ranking and is the
    # one the subject arm can carry: does its top metric match the control's?
    subj = sorted(rows, key=lambda r: -r["combination_articles"])
    out["subject_top_metric"] = (subj[0]["metric"] if
                                 subj[0]["combination_articles"] > subj[1]["combination_articles"]
                                 else None)
    out["arms_agree"] = bool(out["leading_metric"] and
                             out["subject_top_metric"] == out["leading_metric"])
    # THE RATIO THE PROTOCOL CARES ABOUT: among combination papers that assert
    # an interaction, how many put a named scale beside the claim?
    out["claim_with_metric_share"] = (
        round(100 * (d["claim_combination"] - d["claim_only_combination"])
              / d["claim_combination"], 2)
        if d["claim_combination"] else None)
    out["combination_share"] = (round(100 * nc / d["ferroptosis_articles"], 2)
                                if d["ferroptosis_articles"] else None)
    return out


def render(d: dict) -> str:
    def percent(value):
        return "unavailable" if value is None else f"{value}%"

    L = ["# Which synergy metric the ferroptosis literature names\n"]
    L.append(
        f"Generated by `scripts/census_synergy_metrics.py` over the "
        f"{d['ferroptosis_articles']:,} census articles carrying the "
        f"`Ferroptosis` descriptor, of which {d['combination_articles']:,} "
        f"({percent(d['combination_share'])}) describe a combination. "
        f"`analysis/p1-wetlab-protocol.md` states its falsification threshold "
        f"in one metric's terms, and a threshold is only as useful as the "
        f"comparison it enables.\n")
    lead = d["leading_metric"]
    a, b = d["leading_control_articles"], d["runner_up_control_articles"]
    if not d["subject_rankable"]:
        L.append(
            f"**The subject arm cannot rank these metrics.** Only {d['any_metric']} "
            f"of the {d['ferroptosis_articles']:,} ferroptosis articles name a metric, "
            f"below the {d['rankable_min']} this analysis requires before reading "
            "an order from that arm.\n")
    else:
        L.append(
            f"The subject arm names a metric in {d['any_metric']:,} articles and "
            f"meets the {d['rankable_min']}-article threshold. This count threshold "
            "does not establish an ordering with measured uncertainty.\n")
    L.append(
        f"The control samples {d['control_records']:,} census records at 1-in-"
        f"{d['control_stride']} shards.")
    if lead is None:
        L.append("The control has no unique leading metric: its highest "
                 "counts are tied or all zero. No agreement between "
                 "leading metrics can be established.\n")
    else:
        L.append(
            f"The largest observed control count is {a:,} for {lead}, "
            f"against {b:,} for the next. This is a descriptive ordering; "
            "the chronological shard sample and overlapping metric matches "
            "do not supply an independent-count uncertainty model. "
            + ("The subject combination counts have no unique leader.\n"
               if d['subject_top_metric'] is None else
               f"The subject combination leader {'AGREES' if d['arms_agree'] else 'DISAGREES'} "
               "with it. This is descriptive agreement, the most it can support.\n"))
    if lead is None:
        L.append("Without a unique control leader, the protocol's named "
                 "metrics cannot be compared with a leading scale.\n")
    elif d["protocol_names_the_leader"]:
        L.append("The protocol names the metric with the largest observed "
                 "control count.\n")
    else:
        L.append("The protocol does not name the metric with the largest "
                 "observed control count. Its threshold therefore uses a "
                 "different scale from that observed leader.\n")
    L.append("| metric | census control | ferroptosis combination papers | share | in P1 |")
    L.append("|---|--:|--:|--:|---|")
    for r in d["rows"]:
        mark = "**yes**" if r["in_protocol"] else "\u2014"
        L.append(f"| {r['metric']} | {r['control_articles']:,} | "
                 f"{r['combination_articles']:,} | "
                 f"{percent(r['share_of_combination'])} | {mark} |")
    L.append("")
    L.append(
        f"The control column is a count over {d['control_records']:,} records "
        f"of all kinds, so it is not a share of anything comparable to the "
        f"ferroptosis column beside it -- it is there to ORDER the metrics, "
        f"not to size them.\n")
    L.append("## The claim-to-metric gap\n")
    with_metric = d["claim_combination"] - d["claim_only_combination"]
    if d["claim_combination"]:
        L.append(
            f"{d['claim_combination']:,} combination papers assert an interaction "
            f"in words -- *synergistic*, *synergy*, *potentiates* -- and "
            f"{with_metric:,} of them "
            f"({percent(d['claim_with_metric_share'])}) name a metric anywhere a reader "
            f"can see it. The remaining {d['claim_only_combination']:,} state the "
            f"conclusion without the scale it was measured on.\n")
    else:
        L.append("No combination papers in the selected subject records "
                 "assert an interaction using the claim vocabulary. The "
                 "claim-with-metric share is unavailable.\n")
    if d["claim_combination"]:
        L.append(
            "This share measures named metrics alongside claim words in titles "
            "and abstracts. It does not measure whether a metric was computed "
            "elsewhere in the paper. Without measured reporting rates for both "
            "vocabularies, these counts establish neither the true computation "
            "rate nor a bound or bias direction for the full-paper gap.\n")
    L.append("## What this means for the protocol\n")
    bliss = next((r for r in d["rows"] if r["metric"].startswith("Bliss")), None)
    if (bliss and bliss["in_protocol"] and d["protocol_names_the_leader"] and
            lead != bliss["metric"]):
        L.append(
            f"The protocol names both Bliss and {lead}. Bliss is named in {bliss['control_articles']:,} "
            f"control records and {lead} in {d['leading_control_articles']:,}. "
            "These observed mention counts describe the available abstract text, "
            "not which metric is preferable or how often each was computed.\n")
    L.append(
        "Report the raw dose-response matrix alongside whichever index is "
        "computed. An index is a number on a scale, and the scales are not "
        "interconvertible without the underlying data -- a Chou-Talalay CI "
        "and a Bliss excess answer different questions about the same "
        "experiment. A collaborator who publishes the matrix lets a reader "
        "holding the other metric recompute; one who publishes only the index "
        "does not.\n")
    L.append("## What this does not measure\n")
    L.append(
        "Whether a metric was computed. Title and abstract searches may miss "
        "metrics named elsewhere or match a name without its calculation. "
        "Metric-specific reporting rates are unknown, so observed mention "
        "ordering need not match the ordering of computation rates. "
        "This analysis also says nothing about whether a reported synergy is real: "
        "naming a scale is not the same as using it well.\n")
    return "\n".join(L)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--render-only", action="store_true")
    a = ap.parse_args()
    d = assemble(json.loads(OUT_JSON.read_text()) if a.render_only
                 else scan(a.stride))
    json_text = json.dumps(d, indent=1) + "\n"
    md_text = render(d)
    OUT_JSON.write_text(json_text)
    OUT_MD.write_text(md_text)
    print(f"wrote {OUT_MD}")
    print(f"  leader: {d['leading_metric']} (protocol names it: "
          f"{d['protocol_names_the_leader']}); claim-with-metric "
          f"{d['claim_with_metric_share']}%")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
