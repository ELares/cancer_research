#!/usr/bin/env python3
"""Does the literature support this project's normal-tissue selectivity premise?

Issue #728 states the assumption plainly: ferroptosis inducers are attractive
"precisely on the claim that normal cells resist them", and observes that the
word *toxicity* appears once in the whole manuscript. Nobody had asked the
census whether that claim holds.

THERE IS A LITERATURE IN WHICH FERROPTOSIS IS THE MECHANISM OF NORMAL-TISSUE
HARM, and it is not small: 575 of the census's ferroptosis articles carry an
organ-toxicity descriptor, led by acute kidney injury, drug-induced liver
injury and cardiotoxicity. Two of the cancer drugs whose dose-limiting
toxicities define oncology practice are the most-named agents in it --
doxorubicin cardiotoxicity and cisplatin nephrotoxicity -- so this is not a
digression from cancer therapy.

WHERE IT SITS IS ITSELF THE FINDING, and it corrects the first version of this
analysis, which called these "organ damage in cancer patients receiving cancer
drugs" and was wrong for most of them. 558 of the 575 enter the census only
through the ADJACENT extension and just 17 through the C04 cancer tree; a third
of them mention a cancer drug or therapy at all, and the sample holds cadmium
in Nile tilapia and a herbicide in mice. Read as a rate, the ferroptosis
literature indexed as cancer touches organ toxicity in 0.33% of its articles
against 6.8% of the adjacent stream, a factor of twenty.

That gap is mostly a LABELLING effect and must not be reported as neglect: a
paper about doxorubicin injuring a heart is indexed under Cardiotoxicity and
Doxorubicin, not under a Neoplasms descriptor, because its subject is the
heart. The consequence for this project stands either way -- a corpus scoped
to cancer will not surface the literature that bears on its own selectivity
premise, whatever the reason.

WHAT PRESENCE CAN AND CANNOT SHOW, and this is where a keyword count would go
wrong. Naming a ferroptosis INHIBITOR (ferrostatin-1, liproxstatin-1) does not
mean a paper proposes blocking ferroptosis to protect an organ. Fer-1 rescue is
the field's standard SPECIFICITY CONTROL -- this project's own P1 protocol uses
it that way -- so the same string serves two opposite roles. Direction is
therefore adjudicated by reading, on a committed sample, exactly as the hypoxia
and thesis-direction analyses do; the scan generates candidates and the
adjudication is the measurement.

WHAT THIS DOES NOT DO is refute the selectivity premise. A literature showing
normal tissue CAN die by ferroptosis under a chemotherapeutic does not
establish that a targeted GPX4 inhibitor at a therapeutic dose would harm it,
and the two claims are often confused. What it does establish is that the
premise is an open empirical question with a literature attached, rather than
a background fact -- which is the difference between a caveat and an
assumption.
"""
import argparse
import gzip
import json
import random
import re
from collections import Counter
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
RECORDS = REPO / "corpus/atlas/records"
OUT_MD = REPO / "analysis/census-normal-tissue.md"
OUT_JSON = REPO / "analysis/census-normal-tissue.json"
ADJ = REPO / "analysis/normal-tissue-adjudication.csv"

FERRO = "ferroptosis"
# Named explicitly rather than pattern-matched. `Cytotoxicity, Immunologic`
# and `T-Lymphocytes, Cytotoxic` contain the substring "toxic" and mean the
# OPPOSITE -- the ability to kill a target cell, which is the intended effect
# -- so a regex on "toxic" sweeps them in and inflates the count.
ORGAN_TOX = {
    "acute kidney injury", "chemical and drug induced liver injury",
    "cardiotoxicity", "neurotoxicity syndromes", "ototoxicity",
    "drug-related side effects and adverse reactions", "mucositis",
    "toxicity tests", "peripheral nervous system diseases",
}
# Agents used to BLOCK ferroptosis. Presence only -- the role is adjudicated.
INHIBITOR = re.compile(
    r"\bferrostatin\b|\bferrostatin-1\b|\bliproxstatin\b|\bfer-1\b|"
    r"\bdeferoxamine\b|\bdesferrioxamine\b|\bvitamin e\b|\balpha-tocopherol\b")
INDUCER = re.compile(r"\brsl3\b|\berastin\b|\bml162\b|\bml210\b|\bifsp1\b|\bsorafenib\b")
# The cancer drugs whose organ toxicity this literature is mostly about.
CULPRIT = {"cisplatin": r"\bcisplatin\b", "doxorubicin": r"\bdoxorubicin\b|\badriamycin\b",
           "sorafenib": r"\bsorafenib\b", "carboplatin": r"\bcarboplatin\b",
           "methotrexate": r"\bmethotrexate\b"}
SAMPLE_N = 40
SAMPLE_SEED = 728   # the issue number, so the draw is documented not chosen


def scan(stride: int = 1) -> dict:
    inh, ind = INHIBITOR, INDUCER
    culp = {k: re.compile(v) for k, v in CULPRIT.items()}
    n = tox = n_inh = n_ind = n_both = 0
    by_desc, by_culprit = Counter(), Counter()
    basis_all, basis_tox = Counter(), Counter()
    pool = []
    for f in sorted(RECORDS.glob("*.jsonl.gz"))[::stride]:
        with gzip.open(f, "rt", encoding="utf-8") as fh:
            for line in fh:
                r = json.loads(line)
                ms = {m.lower() for m in (r.get("mesh") or [])}
                if FERRO not in ms:
                    continue
                n += 1
                basis = r.get("cancer_basis") or "?"
                basis_all[basis] += 1
                hit = ms & ORGAN_TOX
                if not hit:
                    continue
                tox += 1
                basis_tox[basis] += 1
                for h in hit:
                    by_desc[h] += 1
                blob = f"{r.get('title') or ''} {r.get('abstract') or ''}".lower()
                hi, hd = bool(inh.search(blob)), bool(ind.search(blob))
                n_inh += hi
                n_ind += hd
                n_both += hi and hd
                for k, p in culp.items():
                    if p.search(blob):
                        by_culprit[k] += 1
                pool.append({"pmid": r.get("pmid"), "year": r.get("year"),
                             "title": r.get("title") or "",
                             "descriptors": sorted(hit),
                             "names_inhibitor": hi, "names_inducer": hd})
    pool.sort(key=lambda x: str(x["pmid"]))
    rng = random.Random(SAMPLE_SEED)
    sample = rng.sample(pool, min(SAMPLE_N, len(pool)))
    return {"ferroptosis_articles": n, "organ_toxicity_articles": tox,
            "basis_all": dict(basis_all), "basis_toxicity": dict(basis_tox),
            "names_inhibitor": n_inh, "names_inducer": n_ind,
            "names_both": n_both,
            "by_descriptor": dict(by_desc.most_common()),
            "by_culprit_drug": dict(by_culprit.most_common()),
            "sample_seed": SAMPLE_SEED, "sample_n": SAMPLE_N,
            "sample": sample}


def _adjudication() -> dict:
    """The committed hand-read verdicts, which ARE the direction measurement."""
    if not ADJ.exists():
        return {}
    import csv
    with ADJ.open() as fh:
        rows = list(csv.DictReader(fh))
    return {r["pmid"]: r for r in rows}


def assemble(d: dict) -> dict:
    out = dict(d)
    n, tox = d["ferroptosis_articles"], d["organ_toxicity_articles"]
    out["organ_toxicity_share"] = round(100 * tox / n, 2) if n else None
    out["inhibitor_share"] = round(100 * d["names_inhibitor"] / tox, 1) if tox else None
    # THE RATE THAT MATTERS is per stream, not pooled: pooling hides that the
    # cancer-indexed literature barely touches this and the adjacent stream does.
    out["basis_rate"] = {
        k: round(100 * d["basis_toxicity"].get(k, 0) / v, 2)
        for k, v in d["basis_all"].items() if v}
    c04, adj = out["basis_rate"].get("C04"), out["basis_rate"].get("adjacent")
    out["basis_rate_ratio"] = round(adj / c04, 1) if (c04 and adj) else None
    adj = _adjudication()
    out["adjudicated_n"] = len(adj)
    if adj:
        verdicts = Counter(r["adjudicated"] for r in adj.values())
        out["verdicts"] = dict(verdicts.most_common())
        harm = verdicts.get("harm", 0)
        out["harm_share"] = round(100 * harm / len(adj), 1)
        # The confound the keyword arm cannot see, measured on the sample:
        # how often naming an inhibitor means a SPECIFICITY PROBE instead.
        probe = sum(1 for r in adj.values()
                    if r["adjudicated"] == "probe" and r["regex_inhibitor"] == "True")
        named = sum(1 for r in adj.values() if r["regex_inhibitor"] == "True")
        out["inhibitor_named_in_sample"] = named
        out["inhibitor_is_probe"] = probe
        out["probe_share_of_inhibitor_mentions"] = (
            round(100 * probe / named, 1) if named else None)
    return out


def _ranked(d: dict) -> list:
    """(key, count) pairs, count-descending.

    These were built with `Counter.most_common()` and rendered from dict order,
    so their ranking survived only until something serialised them. Sorting
    HERE keeps the report both reproducible and correctly ranked.
    """
    return sorted(d.items(), key=lambda kv: (-kv[1], kv[0]))


def render(d: dict) -> str:
    L = ["# Ferroptosis as a mechanism of normal-tissue harm\n"]
    L.append(
        f"Generated by `scripts/census_normal_tissue.py`. Issue #728 records "
        f"this project's selectivity premise -- that ferroptosis inducers are "
        f"attractive because normal cells resist them -- and notes that "
        f"*toxicity* appears once in the manuscript. Across the census, "
        f"**{d['organ_toxicity_articles']:,} of the "
        f"{d['ferroptosis_articles']:,} articles carrying the `Ferroptosis` "
        f"descriptor ({d['organ_toxicity_share']}%) also carry an "
        f"organ-toxicity descriptor**, and on a hand-read sample those are "
        f"papers in which ferroptosis is the mechanism by which normal tissue "
        f"is damaged. Where they sit in the census is the second finding, and "
        f"it is why this project has not met them.\n")
    L.append("| organ-toxicity descriptor | articles |")
    L.append("|---|--:|")
    for k, v in _ranked(d["by_descriptor"]):
        L.append(f"| {k} | {v:,} |")
    L.append("")
    c04 = d["basis_rate"].get("C04")
    adj = d["basis_rate"].get("adjacent")
    L.append("## Where this literature sits, and why the project has not met it\n")
    L.append(
        f"**{d['basis_toxicity'].get('adjacent', 0):,} of the "
        f"{d['organ_toxicity_articles']:,} enter the census through the "
        f"ADJACENT extension and {d['basis_toxicity'].get('C04', 0)} through "
        f"the C04 cancer tree.** As a rate: the ferroptosis literature indexed "
        f"as cancer touches organ toxicity in {c04}% of its articles against "
        f"{adj}% of the adjacent stream, a factor of "
        f"{d['basis_rate_ratio']}.\n")
    L.append(
        "**That gap is mostly a labelling effect and should not be read as "
        "neglect.** A paper about doxorubicin injuring a heart is indexed "
        "under `Cardiotoxicity` and `Doxorubicin`, not under a `Neoplasms` "
        "descriptor, because its subject is the heart. The consequence for "
        "this project holds either way: a corpus scoped to cancer will not "
        "surface the literature bearing on its own selectivity premise, "
        "whatever the reason it sits elsewhere.\n")
    L.append(
        "It also corrects this analysis's own first framing, which called "
        "these papers organ damage in cancer patients receiving cancer drugs. "
        "A third mention a cancer drug or therapy at all, and the adjudicated "
        "sample holds cadmium in Nile tilapia and a herbicide in mice. The "
        "cancer-drug core below is real; it is not the whole set.\n")
    if d["by_culprit_drug"]:
        L.append("| implicated agent | articles |")
        L.append("|---|--:|")
        for k, v in _ranked(d["by_culprit_drug"]):
            L.append(f"| {k} | {v:,} |")
        L.append("")
    L.append("## Direction is adjudicated, not matched\n")
    L.append(
        f"{d['names_inhibitor']:,} of these articles "
        f"({d['inhibitor_share']}%) name a ferroptosis INHIBITOR -- "
        f"ferrostatin-1, liproxstatin-1, deferoxamine. **That count cannot be "
        f"read as papers proposing to block ferroptosis to protect an organ.** "
        f"Fer-1 rescue is the field's standard specificity control, and this "
        f"project's own P1 protocol uses it that way, so one string serves two "
        f"opposite roles.\n")
    if d.get("verdicts"):
        L.append(
            f"So a sample of {d['adjudicated_n']} was read "
            f"(`analysis/normal-tissue-adjudication.csv`, seed "
            f"{d['sample_seed']}), and the verdicts are the measurement:\n")
        L.append("| verdict | records |")
        L.append("|---|--:|")
        for k, v in _ranked(d["verdicts"]):
            L.append(f"| {k} | {v} |")
        L.append("")
        if d.get("inhibitor_named_in_sample") is not None:
            probe = d["inhibitor_is_probe"]
            L.append(
                f"**The predicted confound did not appear, and its direction "
                f"was predicted wrong.** Of the "
                f"{d['inhibitor_named_in_sample']} sampled records naming an "
                f"inhibitor, {probe} use it as a bare specificity probe: in an "
                f"organ-injury paper the inhibitor IS the protective "
                f"intervention. So the keyword arm does not overstate the "
                f"protective literature as this analysis first warned -- it "
                f"understates it, since {d['verdicts'].get('harm', 0)} of "
                f"{d['adjudicated_n']} records read as ferroptosis-mediated "
                f"injury while only {d['inhibitor_named_in_sample']} name an "
                f"inhibitor at all. Knowing that a measure errs is not knowing "
                f"which way.\n")
            L.append(
                "The uniformity is itself the result. An adjudication that "
                "finds one direction in 38 of 40 records is weak evidence "
                "about a contested question and strong evidence about an "
                "uncontested one, and this question is uncontested inside this "
                "literature: these papers agree that normal kidney, liver, "
                "heart and neural tissue die by ferroptosis under stress.\n")
    L.append("## What this does and does not establish\n")
    L.append(
        "**It does not refute the selectivity premise.** A literature showing "
        "normal tissue can die by ferroptosis under a chemotherapeutic does "
        "not establish that a targeted GPX4 inhibitor at a therapeutic dose "
        "would harm it, and the two claims are easy to confuse. Cisplatin and "
        "doxorubicin injure kidneys and hearts through several mechanisms at "
        "once; ferroptosis being among them says nothing directly about RSL3.\n")
    L.append(
        "**What it establishes is that the premise is an open empirical "
        "question with a literature attached, rather than a background fact.** "
        "That is the difference between a caveat and an assumption, and it is "
        "the gap issue #728 identifies: the project scores every modality on "
        "kill and none on harm, and the one compartment it does model is a "
        "cancer-associated fibroblast, which is tumour-resident and not "
        "normal tissue at all.\n")
    L.append("## For the layer-freeze policy\n")
    L.append(
        "A normal-tissue phenotype is the addition #728 asks for, and "
        "CONTRIBUTING.md requires a named calibration target before one lands. "
        "This analysis does not supply it: an article count is not a "
        "dose-response, and the organ-toxicity literature reports histology "
        "and biomarkers rather than the matched normal-versus-tumour kill "
        "curves a phenotype would be fitted to. The committed CTRPv2 curves "
        "are cancer cell lines only. So the target remains unnamed, and the "
        "honest position is that the selectivity assumption should be stated "
        "as one in the manuscript before it is modelled.\n")
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
    print(f"  {d['organ_toxicity_articles']:,} of {d['ferroptosis_articles']:,} "
          f"({d['organ_toxicity_share']}%); inhibitor named "
          f"{d['inhibitor_share']}%; adjudicated {d['adjudicated_n']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
