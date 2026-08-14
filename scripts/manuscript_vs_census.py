#!/usr/bin/env python3
"""Testing the manuscript's own quantitative claims against the census.

WHY THIS EXISTS
---------------
The manuscript's numbers were computed over a keyword-retrieved corpus, or over
live PubMed queries run once. The census is the whole indexed cancer
literature with expert MeSH labels. Several claims that the frozen corpus could
not test are now testable, and a project whose stated principle is "let the
evidence lead" should test its own published claims first.

WHAT THIS IS NOT
----------------
Not a replication and not a refutation exercise. Two of the claims tested here
SURVIVE and one of them is understated by the manuscript, which is reported as
prominently as any failure would be. A document that only ever finds its own
work wanting is as unreliable as one that only ever confirms it.

THE SCOPE DIFFERENCE THAT DECIDES HOW EVERY ROW READS
------------------------------------------------------
The census is cancer-restricted: MeSH tree C04 plus a set of adjacent
descriptors.
The manuscript's modality counts came from UNRESTRICTED PubMed queries, so they
include ferroptosis work outside oncology. Absolute counts therefore differ
several-fold BY CONSTRUCTION and are not evidence of anything. Only ratios
within a table, and directions, are comparable. Every row below states which it
is testing.

A SECOND SCOPE LIMIT, ON THE OTHER SIDE. A claim is only testable here if its
concept has a MeSH descriptor. TTFields has none, so the manuscript's zero for
it can be neither confirmed nor contradicted, and rows whose census count is
too small to support a ratio are reported as unmeasurable rather than as a
finding.

Usage:
    python scripts/manuscript_vs_census.py
"""

import collections
import gzip
import json
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
RECORDS = PROJECT_ROOT / "corpus" / "atlas" / "records"
MANUSCRIPT = PROJECT_ROOT / "article" / "drafts" / "v1.md"
OUT_MD = PROJECT_ROOT / "analysis" / "manuscript-vs-census.md"
OUT_JSON = PROJECT_ROOT / "analysis" / "manuscript-vs-census.json"

FERROPTOSIS = "Ferroptosis"

# The modality table in section 8.2, with the descriptor each row maps to and
# the figure the manuscript states. `None` means the concept has no MeSH
# descriptor and the row cannot be tested at all.
MODALITIES = [
    ("PDT", "Photochemotherapy", 355, 67),
    ("SDT", "Ultrasonic Therapy", 121, 25),
    ("IRE", "Electroporation", 15, None),
    ("HIFU", "High-Intensity Focused Ultrasound Ablation", 3, None),
    ("TTFields", None, 0, 0),
]

# A ratio below this many articles on either side is not a ratio.
MIN_FOR_A_RATIO = 10

# DESCRIPTOR SENSITIVITY. MeSH gives PDT both a procedure and an agent
# descriptor while SDT has only a procedure one plus generic physics terms, so
# genuine sonodynamic papers land outside `Ultrasonic Therapy`. A single
# descriptor pair therefore ships a point estimate whose margin depends on a
# choice the analysis made. Each variant is computed and the range reported.
DESCRIPTOR_VARIANTS = {
    "as the manuscript frames them": (
        ["Photochemotherapy"], ["Ultrasonic Therapy"]),
    "SDT widened": (
        ["Photochemotherapy"],
        ["Ultrasonic Therapy", "Ultrasonic Waves", "Sonication", "Ultrasonics",
         "Microbubbles"]),
    "PDT widened": (
        ["Photochemotherapy", "Photosensitizing Agents",
         "Low-Level Light Therapy", "Phototherapy"],
        ["Ultrasonic Therapy"]),
    "both widened": (
        ["Photochemotherapy", "Photosensitizing Agents",
         "Low-Level Light Therapy", "Phototherapy"],
        ["Ultrasonic Therapy", "Ultrasonic Waves", "Sonication", "Ultrasonics",
         "Microbubbles"]),
}

# THE INVALIDATOR THIS ANALYSIS NAMED AND DID NOT MEASURE. Both descriptors are
# broader than the modality they stand for, so the ratio is only as good as
# their RELATIVE over-estimation. Measured by asking whether each record's own
# title and abstract discuss the modality and a tumour.
PDT_WORDS = ("photodynamic", "photosensitiz", "photosensitis", "photo-activat")
SDT_WORDS = ("sonodynamic", "sono-dynamic", "sonosensitiz", "sonosensitis",
             "ultrasound-activat", "ultrasound activat")
TUMOUR_WORDS = ("tumor", "tumour", "cancer", "carcinoma", "neoplas", "melanoma",
                "glioma", "sarcoma", "leukemi", "lymphoma", "metasta")

# The growth claim in section 3.7, and the two confounds the census can rule
# out or in. The corpus was retrieved ONCE with fixed queries, so the only
# year-dependent retrieval effect is how much of each year is open access.
GROWTH_START, GROWTH_END = 2015, 2025
CORPUS_GROWTH_START, CORPUS_GROWTH_END = 38, 1167


def understates(census_ratio, manuscript_ratio, measurable) -> bool:
    """Does the census show a LARGER ratio than the manuscript argued from?

    A function, not an expression inlined into the result dict, so a test can
    hand it inputs where the answer must be False. A guard that recomputes the
    verdict from the same fields the generator wrote cannot tell a live
    comparison from a hardcoded True while the claim happens to hold.
    """
    if not (measurable and census_ratio and manuscript_ratio):
        return False
    return census_ratio > manuscript_ratio


def outgrew(corpus, field) -> bool:
    """Did the corpus grow faster than the field it was drawn from?"""
    return bool(corpus and field and corpus > field)


def ratio_is_measurable(a, b, floor=None) -> bool:
    """Are both sides large enough that their ratio means anything?"""
    return min(a, b) >= (MIN_FOR_A_RATIO if floor is None else floor)


def scan():
    """One pass over the census records."""
    want = {d for _, d, _, _ in MODALITIES if d}
    want |= {d for a, b in DESCRIPTOR_VARIANTS.values() for d in a + b}
    on_modality = collections.Counter()
    variant_hits = collections.defaultdict(set)
    mod_total = collections.Counter()
    mod_ferro = collections.Counter()
    mod_ferro_icd = collections.Counter()
    year_total = collections.Counter()
    year_oa = collections.Counter()
    n = ferro = 0
    for f in sorted(RECORDS.glob("*.jsonl.gz")):
        with gzip.open(f, "rt") as fh:
            for line in fh:
                try:
                    r = json.loads(line)
                except json.JSONDecodeError:
                    continue
                mesh = set(r.get("mesh") or ())
                y = r.get("year")
                if y:
                    try:
                        y = int(y)
                    except (TypeError, ValueError):
                        y = None
                if y and GROWTH_START <= y <= GROWTH_END:
                    year_total[y] += 1
                    if r.get("pmcid"):
                        year_oa[y] += 1
                if not mesh:
                    continue
                n += 1
                is_f = FERROPTOSIS in mesh
                ferro += is_f
                is_icd = "Immunogenic Cell Death" in mesh
                hit = want & mesh
                for d in hit:
                    mod_total[d] += 1
                    if is_f:
                        mod_ferro[d] += 1
                        mod_ferro_icd[d] += is_icd
                if is_f and hit:
                    pm = r.get("pmid")
                    for d in hit:
                        variant_hits[d].add(pm)
                    # the relative-over-estimation measurement
                    txt = ((r.get("title") or "") + " " +
                           (r.get("abstract") or "")).lower()
                    tumour = any(w in txt for w in TUMOUR_WORDS)
                    if "Photochemotherapy" in hit and tumour and \
                            any(w in txt for w in PDT_WORDS):
                        on_modality["PDT"] += 1
                    if "Ultrasonic Therapy" in hit and tumour and \
                            any(w in txt for w in SDT_WORDS):
                        on_modality["SDT"] += 1
    return {"records_with_mesh": n, "ferroptosis_records": ferro,
            "on_modality": on_modality, "variant_hits": variant_hits,
            "mod_total": mod_total, "mod_ferro": mod_ferro,
            "mod_ferro_icd": mod_ferro_icd,
            "year_total": year_total, "year_oa": year_oa}


def main() -> int:
    if not RECORDS.exists():
        print(f"census records absent at {RECORDS}; nothing to test")
        return 1
    s = scan()

    rows = []
    for name, desc, ms_ferro, ms_icd in MODALITIES:
        rows.append({
            "modality": name, "descriptor": desc,
            "manuscript_ferroptosis": ms_ferro,
            "manuscript_ferroptosis_icd": ms_icd,
            "census_total": s["mod_total"].get(desc, 0) if desc else None,
            "census_ferroptosis": s["mod_ferro"].get(desc, 0) if desc else None,
            "census_ferroptosis_icd": s["mod_ferro_icd"].get(desc, 0) if desc else None,
            "testable": bool(desc),
        })
    by = {r["modality"]: r for r in rows}
    pdt, sdt = by["PDT"], by["SDT"]

    def ratio(a, b):
        return round(a / b, 2) if b else None

    pdt_ferro_raw = pdt["census_ferroptosis"]
    sdt_ferro_raw = sdt["census_ferroptosis"]
    ms_ratio = ratio(pdt["manuscript_ferroptosis"], sdt["manuscript_ferroptosis"])
    cs_ratio = ratio(pdt["census_ferroptosis"], sdt["census_ferroptosis"])
    ratio_measurable = ratio_is_measurable(pdt["census_ferroptosis"],
                                           sdt["census_ferroptosis"])
    icd_measurable = ratio_is_measurable(pdt["census_ferroptosis_icd"],
                                         sdt["census_ferroptosis_icd"])

    # the variant range, and whether the direction survives the SDT-most-
    # favourable descriptor set rather than only the one this analysis chose
    variants = []
    for label, (pd_, sd_) in DESCRIPTOR_VARIANTS.items():
        a = len(set().union(*(s["variant_hits"][d] for d in pd_)))
        b = len(set().union(*(s["variant_hits"][d] for d in sd_)))
        variants.append({"variant": label, "pdt": a, "sdt": b,
                         "ratio": round(a / b, 2) if b else None,
                         "descriptors_pdt": pd_, "descriptors_sdt": sd_})
    v_ratios = [v["ratio"] for v in variants if v["ratio"]]
    om = s["on_modality"]
    pdt_on = round(100.0 * om["PDT"] / max(pdt_ferro_raw, 1), 1)
    sdt_on = round(100.0 * om["SDT"] / max(sdt_ferro_raw, 1), 1)

    yt, yo = s["year_total"], s["year_oa"]
    field = ratio(yt.get(GROWTH_END, 0), yt.get(GROWTH_START, 0))
    oa = ratio(yo.get(GROWTH_END, 0), yo.get(GROWTH_START, 0))
    corpus = ratio(CORPUS_GROWTH_END, CORPUS_GROWTH_START)

    res = {
        "records_with_mesh": s["records_with_mesh"],
        "ferroptosis_records": s["ferroptosis_records"],
        "modality_table": {
            "rows": rows,
            "manuscript_pdt_sdt_ratio": ms_ratio,
            "census_pdt_sdt_ratio": cs_ratio,
            "ratio_is_measurable": ratio_measurable,
            "census_exceeds_manuscript": understates(
                cs_ratio, ms_ratio, ratio_measurable),
            "icd_column_is_measurable": icd_measurable,
            "untestable_rows": [r["modality"] for r in rows if not r["testable"]],
            "rows_below_ratio_floor": [
                r["modality"] for r in rows
                if r["testable"] and (r["census_ferroptosis"] or 0) < MIN_FOR_A_RATIO],
            "min_for_a_ratio": MIN_FOR_A_RATIO,
            # The descriptor choice is a choice; here is what it is worth.
            "descriptor_variants": variants,
            "ratio_range": [min(v_ratios), max(v_ratios)] if v_ratios else None,
            "direction_holds_under_every_variant": bool(
                v_ratios and ms_ratio and min(v_ratios) > 1.0),
            "understatement_holds_under_every_variant": bool(
                v_ratios and ms_ratio and min(v_ratios) > ms_ratio),
            # The invalidator this analysis names: is the over-estimation
            # SYMMETRIC? If PDT's descriptor is proportionally broader than
            # SDT's, the ratio is inflated and the verdict is unearned.
            "on_modality_and_tumour": {
                "pdt_pct": pdt_on, "sdt_pct": sdt_on,
                "pdt_n": om["PDT"], "sdt_n": om["SDT"],
                "gap_points": round(abs(pdt_on - sdt_on), 1),
                "symmetric_within_5_points": abs(pdt_on - sdt_on) <= 5.0,
                "filtered_ratio": round(om["PDT"] / om["SDT"], 2) if om["SDT"] else None,
            },
            "direction_holds": bool(
                ratio_measurable and cs_ratio and cs_ratio > 1.0),
        },
        "growth": {
            "start_year": GROWTH_START, "end_year": GROWTH_END,
            "census_start": yt.get(GROWTH_START, 0),
            "census_end": yt.get(GROWTH_END, 0),
            "census_growth": field,
            "open_access_start": yo.get(GROWTH_START, 0),
            "open_access_end": yo.get(GROWTH_END, 0),
            "open_access_growth": oa,
            "corpus_start": CORPUS_GROWTH_START,
            "corpus_end": CORPUS_GROWTH_END,
            "corpus_growth": corpus,
            "corpus_exceeds_field": outgrew(corpus, field),
            "unexplained_by_availability": round(corpus / oa, 2) if oa else None,
            "per_year": [{"year": y, "census": yt.get(y, 0), "open_access": yo.get(y, 0),
                          "oa_share": round(100.0 * yo.get(y, 0) / yt[y], 1)}
                         for y in sorted(yt)],
        },
    }
    OUT_JSON.write_text(json.dumps(res, indent=2, sort_keys=True) + "\n")
    OUT_MD.write_text(render(res), encoding="utf-8")
    print(f"wrote {OUT_MD}\nwrote {OUT_JSON}")
    print(f"  section 8.2 PDT:SDT  manuscript {ms_ratio} -> census {cs_ratio}")
    print(f"  section 3.7 growth   corpus {corpus} vs field {field} "
          f"vs open access {oa}")
    return 0


def _headline(r: dict) -> str:
    """One sentence stating the outcome, computed from the verdicts."""
    mt, g = r["modality_table"], r["growth"]
    held = [n for n, v in (("section 8.2", mt["census_exceeds_manuscript"]
                            or mt["direction_holds"]),
                           ("section 3.7", g["corpus_exceeds_field"])) if v]
    failed = [n for n in ("section 8.2", "section 3.7") if n not in held]
    parts = []
    if held:
        parts.append(f"{'Both' if len(held) == 2 else held[0].capitalize()} "
                     + ("claims tested here survive" if len(held) == 2
                        else "survives"))
    if failed:
        parts.append(("but " if held else "")
                     + f"{' and '.join(failed)} does not")
    tail = (" and the manuscript understates one of them."
            if mt["census_exceeds_manuscript"] else ".")
    return ("Nothing tested here survives." if not held
            else ", ".join(parts) + tail)


def render(r: dict) -> str:
    m, g = r["modality_table"], r["growth"]
    pdt = next(x for x in m["rows"] if x["modality"] == "PDT")
    sdt = next(x for x in m["rows"] if x["modality"] == "SDT")
    L = [
        "# The manuscript's own claims, tested against the census", "",
        "Generated by `scripts/manuscript_vs_census.py`.", "",
        "The manuscript's numbers were computed over a keyword-retrieved corpus,",
        "or over live PubMed queries run once. The census is the whole indexed",
        "cancer literature with expert MeSH labels, so several claims the frozen",
        "corpus could not test now are testable.", "",
        # DERIVED. This sentence was a static literal ungated by either verdict,
        # so forcing the comparison to fail produced a document whose section
        # body said the ratio was smaller while its opening line still said the
        # claim survived -- and the guard written to keep the outcome prominent
        # asserted the fixed words, pinning the endorsement in place.
        f"**{_headline(r)}** That is reported as prominently as the opposite",
        "would be: a document that only ever finds its own work wanting is as",
        "unreliable as one that only ever confirms it.", "",
        "## The scope difference that decides how every row reads", "",
        "The census is cancer-restricted (MeSH tree C04 plus a set of adjacent",
        "descriptors). The manuscript's modality counts came from UNRESTRICTED",
        "PubMed queries and include ferroptosis work outside oncology. Absolute",
        "counts therefore differ several-fold **by construction** and are",
        "evidence of nothing. Only ratios within a table, and directions, are",
        "comparable.", "",
        "## Section 8.2: the modality table", "",
        "This is the thesis chapter's only quantitative table and it carries the",
        "manuscript's central self-correction, that SDT is not the unique",
        "ferroptosis-ICD modality and PDT is the established one.", "",
        "| modality | descriptor | census, all | census, x Ferroptosis | manuscript |",
        "|---|---|--:|--:|--:|"]
    for x in m["rows"]:
        if not x["testable"]:
            L.append(f"| {x['modality']} | _no MeSH descriptor_ | -- | -- | "
                     f"{x['manuscript_ferroptosis']} |")
            continue
        L.append(f"| {x['modality']} | {x['descriptor']} | "
                 f"{x['census_total']:,} | {x['census_ferroptosis']:,} | "
                 f"{x['manuscript_ferroptosis']} |")
    L += ["",
          f"The manuscript states the ratio as **{m['manuscript_pdt_sdt_ratio']}:1**",
          "and reads it as PDT dominating SDT by approximately three to one.", ""]
    if m["ratio_is_measurable"]:
        L += [f"On the census it is **{m['census_pdt_sdt_ratio']}:1** "
              f"({pdt['census_ferroptosis']} against {sdt['census_ferroptosis']}).",
              ""]
        if m["census_exceeds_manuscript"]:
            L += ["**The ferroptosis-count leg of the claim survives, on a larger",
                  "ratio than the manuscript argued from.** That is the whole of",
                  "what a count ratio can establish, and the manuscript's claim",
                  "is broader than it.", "",
                  "NOT confirmed here, and each for a different reason: the",
                  "manuscript states the ratio holds on every ferroptosis AND ICD",
                  "metric, and the ICD leg is unmeasurable at census scale (see",
                  "below); the depth-penetration argument is tissue optics, which",
                  "no publication count bears on; and 'more clinical experience,",
                  "approved photosensitizers' is not a literature-count claim",
                  "either. A count ratio speaks to attention, not to any of",
                  "those.", ""]
        else:
            L += ["**The census ratio is smaller than the manuscript's**, so the",
                  "strength of the claim needs revisiting even though its",
                  "direction holds.", ""]
    else:
        L += ["The census cannot decide the ratio: one side falls below the",
              f"{m['min_for_a_ratio']}-article floor this analysis uses.", ""]
    om, vs = m["on_modality_and_tumour"], m["descriptor_variants"]
    L += ["### The descriptor choice, and what it is worth", "",
          "MeSH gives PDT both a procedure and an agent descriptor while SDT has",
          "only a procedure one plus generic physics terms, so genuine",
          "sonodynamic papers land outside `Ultrasonic Therapy`. A single pair",
          "therefore ships a point estimate whose margin depends on a choice",
          "this analysis made. Every defensible variant:", "",
          "| descriptor set | PDT | SDT | ratio |", "|---|--:|--:|--:|"]
    for v in vs:
        L.append(f"| {v['variant']} | {v['pdt']} | {v['sdt']} | {v['ratio']} |")
    L += ["",
          f"The ratio ranges {m['ratio_range'][0]} to {m['ratio_range'][1]}."
          + (" **Every variant exceeds the manuscript's**, so neither the"
             " direction nor the understatement rests on the descriptor pair"
             " chosen here."
             if m["understatement_holds_under_every_variant"] else
             " **Not every variant exceeds the manuscript's**, so the"
             " conclusion depends on which descriptors are used and the"
             " point estimate should not be read alone."), "",
          "### Is the over-estimation symmetric?", "",
          "Both descriptors are broader than the modality they name, so the",
          "ratio is only as good as their RELATIVE over-estimation. An earlier",
          "version stated that and never measured it, which left a named",
          "invalidator sitting beside an unconditional verdict. Measured by",
          "asking whether each record's own title and abstract discuss the",
          "modality and a tumour:", "",
          "| | on modality and tumour | of |", "|---|--:|--:|",
          f"| PDT | {om['pdt_n']} ({om['pdt_pct']}%) | "
          f"{next(x['census_ferroptosis'] for x in m['rows'] if x['modality'] == 'PDT')} |",
          f"| SDT | {om['sdt_n']} ({om['sdt_pct']}%) | "
          f"{next(x['census_ferroptosis'] for x in m['rows'] if x['modality'] == 'SDT')} |",
          "",
          (f"The gap is {om['gap_points']} points, so the over-estimation is "
           f"symmetric and the ratio survives it: filtering to on-modality "
           f"records gives {om['filtered_ratio']} against the raw "
           f"{m['census_pdt_sdt_ratio']}."
           if om["symmetric_within_5_points"] else
           f"The gap is {om['gap_points']} points, so one descriptor is "
           "materially broader than the other and the raw ratio is inflated. "
           f"The filtered ratio is {om['filtered_ratio']}."), "",
          "### What this table cannot test", "",
          f"* **{', '.join(m['untestable_rows'])}** has no MeSH descriptor, so the",
          "  manuscript's figure for it can be neither confirmed nor",
          "  contradicted here.",
          f"* **{', '.join(m['rows_below_ratio_floor'])}** fall below the",
          f"  {m['min_for_a_ratio']}-article floor at census scale, so their rows",
          "  carry no census signal and are reported as unmeasurable rather than",
          "  as a finding.",
          "* **The ICD column is not measurable.** Both sides are far below the",
          "  floor, so no ratio is computed from it. The manuscript's ICD figures",
          "  came from a broader keyword query, not a descriptor intersection.",
          "* `Ultrasonic Therapy` is broader than sonodynamic therapy and",
          "  `Photochemotherapy` is broader than tumour PDT, which the manuscript",
          "  already states. Both legs are over-estimates and the ratio is only",
          "  as good as their relative over-estimation.", "",
          "## Section 3.7: the growth claim", "",
          "The manuscript reports the corpus growing from",
          f"{g['corpus_start']} full-text articles in {g['start_year']} to",
          f"{g['corpus_end']} in {g['end_year']}, and attributes it to the",
          "maturation of checkpoint immunotherapy, cell therapy and RNA-platform",
          "attention. The census can rule two confounds in or out.", "",
          "| | " f"{g['start_year']}" " | " f"{g['end_year']}" " | growth |",
          "|---|--:|--:|--:|",
          f"| cancer literature (census) | {g['census_start']:,} | "
          f"{g['census_end']:,} | x{g['census_growth']} |",
          f"| ...of it open access | {g['open_access_start']:,} | "
          f"{g['open_access_end']:,} | x{g['open_access_growth']} |",
          f"| **the manuscript's corpus** | {g['corpus_start']} | "
          f"{g['corpus_end']:,} | **x{g['corpus_growth']}** |", ""]
    if g["corpus_exceeds_field"]:
        L += ["**The claim survives, and the census strengthens it by removing",
              "the two obvious confounds.** The corpus was retrieved once with",
              "fixed queries, so the only year-dependent retrieval effect is how",
              "much of each year is open access. Neither the field's growth nor",
              "the rise in availability comes close to accounting for the",
              f"corpus's, which is **x{g['unexplained_by_availability']}** larger",
              "than availability growth alone.", "",
              "So the mechanisms the corpus tracks did grow far faster than",
              "cancer literature as a whole, which is what the manuscript",
              "attributes the rise to.", ""]
    else:
        L += ["**The corpus did not outgrow the field**, so the attribution in",
              "section 3.7 needs revisiting.", ""]
    L += ["| year | census | open access | share |", "|--:|--:|--:|--:|"]
    for y in g["per_year"]:
        L.append(f"| {y['year']} | {y['census']:,} | {y['open_access']:,} | "
                 f"{y['oa_share']}% |")
    L += ["", "## What this analysis cannot say", "",
          "* **A MeSH descriptor is not the concept.** Every row inherits the",
          "  breadth of the descriptor it maps to, and two of them are broader",
          "  than the modality they stand for.",
          "* **MeSH indexing lags**, so the most recent years are undercounted",
          "  and every growth figure here is a lower bound.",
          "* **This tests two claims.** The manuscript makes many more, and",
          "  nothing here speaks to those.",
          "* **Surviving a census test is not validation.** It removes two",
          "  specific confounds from one claim and shows a ratio holds on",
          "  independent labels. Neither makes the underlying biology true.", ""]
    return "\n".join(L) + "\n"


if __name__ == "__main__":
    raise SystemExit(main())
