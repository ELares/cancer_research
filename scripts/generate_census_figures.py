#!/usr/bin/env python3
"""Figures for the census layer (#CENSUS-FIG).

WHY THIS EXISTS
---------------
The census produced the finding that bounds the widest set of manuscript claims
-- the frozen corpus is a 213-fold non-uniform sample of the cancer literature --
and it had no figure. Every relative prevalence the manuscript computes over the
4,830-article corpus inherits that unevenness, so a reader deciding how much to
trust a cross-mechanism comparison needs to see its shape, not a sentence about
it.

It is a separate script from `generate_figures.py` because the inputs are
different in kind: this reads only committed atlas JSON, so it runs offline in
seconds, while `generate_figures.py` loads the full corpus and several gitignored
simulation outputs.

FIGURES
-------
fig28  Two panels. (A) per-mechanism capture -- the frozen corpus's share of the
       census for that mechanism -- on a log axis, because a linear one collapses
       everything below 10% into the axis. (B) what that unevenness does to the
       manuscript's central corpus claim.

Five more replace corpus-derived figures whose CAPTIONS were rewritten for the
census. That gap is worth naming rather than quietly closing: a caption
describing one measurement above an image plotting another is worse than
leaving both stale, because a reader checks the caption and trusts the picture.

fig2c   census publication volume by year, with the retrieved corpus overlaid on
        a second axis -- the 31-fold rise against the census's 1.10-fold IS the
        retrieval effect, so the two curves are the argument.
fig9c   study-design composition from NLM publication types and check tags, with
        the UNDETERMINED share drawn rather than dropped, since it is the
        largest class and omitting it would imply the census classifies
        everything.
fig14c  mechanism class by anatomical site: enrichment against each site's own
        share, physical against pharmacological, so a site's general prominence
        cancels and only a disagreement between classes carries.
fig15c  the ten most frequent mechanism pairs. Counts, never a RATE -- Section
        3.13 shows the rate is a property of the labelling instrument.
fig16c  clinical-trial share against volume on a log axis, the replacement for a
        weighted composite whose ranking moved seven places under a defensible
        reweighting.

Usage:
    python scripts/generate_census_figures.py
"""

import json
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))

from atlas_landscape import PHARMACOLOGICAL, PHYSICAL  # noqa: E402
from config import PROJECT_ROOT  # noqa: E402

FIG_DIR = PROJECT_ROOT / "article" / "figures"
LANDSCAPE = PROJECT_ROOT / "analysis" / "atlas-landscape.json"
ANALYSIS = PROJECT_ROOT / "analysis"
GROWTH = ANALYSIS / "census-mechanism-growth.json"
DESIGN = ANALYSIS / "census-evidence-design.json"
SITES = ANALYSIS / "census-mechanism-sites.json"
PROFILE = ANALYSIS / "census-mechanism-profile.json"
RATIO = ANALYSIS / "atlas-modality-ratio.json"
MATRIX = ANALYSIS / "census-mechanism-cancer-matrix.json"
# The manuscript's own keyword-method figure, for the fig1c reference line. It
# is the thing the census ratios STRADDLE, so it is the panel's whole point.
MANUSCRIPT_RATIO = 9.1

# The retrieved corpus's own volume, for the fig2c overlay. Read from the
# frozen index rather than typed, so the contrast cannot be drawn against a
# number nobody can check.
FROZEN_INDEX = PROJECT_ROOT / "corpus" / "INDEX.jsonl"

plt.rcParams.update({
    "font.size": 11, "font.family": "serif", "axes.titlesize": 13,
    "axes.labelsize": 12, "figure.dpi": 300, "savefig.dpi": 300,
    "savefig.bbox": "tight",
})


def fig28_census_capture():
    """The frozen corpus's capture of the census, per mechanism."""
    d = json.loads(LANDSCAPE.read_text())
    # A mechanism with no census articles has no defined capture. Dropping it is
    # the honest move: 0/0 plotted as 0% would read as "never captured". A
    # mechanism present in the census but absent from the corpus is dropped for
    # the same reason -- capture 0 has no position on a log axis, and it would
    # otherwise divide by zero in the spread below.
    rows, dropped = [], []
    for r in d["rows"]:
        if r["mesh_census"] > 0 and r["mesh_frozen"] > 0:
            rows.append({**r, "capture": r["mesh_frozen"] / r["mesh_census"]})
        else:
            dropped.append(r["mechanism"])
    rows.sort(key=lambda r: r["capture"])

    fig, (ax, bx) = plt.subplots(
        1, 2, figsize=(14, 6.6), gridspec_kw={"width_ratios": [2, 1.15]})

    names = [r["mechanism"] for r in rows]
    caps = [100 * r["capture"] for r in rows]
    colors = ["#b5651d" if r["mechanism"] in PHYSICAL else "#4a6f8a" for r in rows]

    # A dot-and-stem plot, NOT bars. A bar encodes magnitude by length measured
    # from zero, and a log axis has no zero -- matplotlib picks the left limit,
    # so bar lengths would encode that arbitrary choice rather than the data. At
    # the default limit the longest:shortest bar reads about 21x against a true
    # 213x, flattening the finding by an order of magnitude. Position on a log
    # axis is the honest encoding, and the stems are a reading aid only.
    y = range(len(rows))
    ax.hlines(y, min(caps) * 0.55, caps, color=colors, alpha=0.35, lw=1.6)
    ax.scatter(caps, y, color=colors, s=58, zorder=3)
    ax.set_yticks(list(y))
    ax.set_yticklabels(names)
    ax.set_xscale("log")
    ax.set_xlim(min(caps) * 0.5, max(caps) * 2.6)
    ax.set_xlabel("share of the census this corpus captures (%, log scale)")
    ax.set_title("A. The frozen corpus samples the literature unevenly")
    ax.grid(axis="x", alpha=0.18, which="both")
    ax.set_axisbelow(True)

    for i, c in enumerate(caps):
        ax.text(c * 1.13, i, f"{c:.2f}%", va="center", fontsize=8.5, color="#333")

    spread = caps[-1] / caps[0]
    ax.set_ylim(-2.4, len(rows) - 0.4)
    ax.annotate(
        f"{spread:.0f}-fold spread:  {names[0]} {caps[0]:.2f}%   to   "
        f"{names[-1]} {caps[-1]:.1f}%",
        xy=(0.02, -1.55), xycoords=("axes fraction", "data"), ha="left",
        va="center", fontsize=9.5,
        bbox=dict(boxstyle="round,pad=0.45", fc="#f4efe6", ec="#b5651d", lw=1.1))
    ax.legend(handles=[
        plt.Rectangle((0, 0), 1, 1, fc="#b5651d"),
        plt.Rectangle((0, 0), 1, 1, fc="#4a6f8a")],
        labels=["physical modality", "other mechanism"],
        loc="lower right", fontsize=9, framealpha=0.95)

    # Panel B: all THREE arms of the committed design. Showing only the first and
    # last attributes the whole move to corpus selection, and it does not
    # reconcile -- 9.1 x 3.3 is 30, not 17.6. The middle arm is where the missing
    # factor lives, and it runs the OTHER way: holding the corpus fixed and
    # switching keyword labels for MeSH CUTS the ratio, so by the more
    # independent labelling method the manuscript overstates. Recomputed here
    # from the committed rows rather than quoted, so it cannot drift from panel A.
    R = {r["mechanism"]: r for r in rows}
    tot = lambda ks, c: sum(R[k][c] for k in ks if k in R and R[k].get(c))  # noqa: E731
    arms = [
        ("this corpus,\nkeyword labels\n(what the\nmanuscript reports)",
         tot(PHARMACOLOGICAL, "keyword_frozen") / tot(PHYSICAL, "keyword_frozen"), "#4a6f8a"),
        ("same articles,\nMeSH labels\n(method effect\nalone)",
         tot(PHARMACOLOGICAL, "mesh_frozen") / tot(PHYSICAL, "mesh_frozen"), "#9aaebd"),
        ("the census,\nMeSH labels",
         tot(PHARMACOLOGICAL, "mesh_census") / tot(PHYSICAL, "mesh_census"), "#2f4f4f"),
    ]
    cap_phys = 100 * tot(PHYSICAL, "mesh_frozen") / tot(PHYSICAL, "mesh_census")
    cap_pharm = 100 * tot(PHARMACOLOGICAL, "mesh_frozen") / tot(PHARMACOLOGICAL, "mesh_census")

    vals = [a[1] for a in arms]
    bx.bar(range(3), vals, color=[a[2] for a in arms], width=0.6)
    bx.set_xticks(range(3))
    bx.set_xticklabels([a[0] for a in arms], fontsize=8.5)
    bx.set_ylabel("pharmacological : physical articles")
    bx.set_title("B. Method and corpus pull in opposite directions")
    for x, v in enumerate(vals):
        bx.text(x, v + max(vals) * 0.03, f"{v:.1f} : 1", ha="center",
                fontweight="bold", fontsize=11.5)
    bx.set_ylim(0, max(vals) * 1.62)
    bx.annotate(
        f"Switching keyword labels for MeSH on the SAME\n"
        f"articles CUTS the ratio ({vals[1]/vals[0]:.2f}x). Widening to the\n"
        f"census then RAISES it {vals[2]/vals[1]:.1f}x, because the corpus\n"
        f"over-samples physical modalities {cap_phys/cap_pharm:.1f}x\n"
        f"({cap_phys:.2f}% vs {cap_pharm:.2f}% capture).\n"
        f"Net {vals[2]/vals[0]:.1f}x: the claim survives, but the\n"
        f"corpus effect is larger than the net move shows.",
        xy=(0.5, 0.875), xycoords="axes fraction", ha="center", va="center",
        fontsize=8.2,
        bbox=dict(boxstyle="round,pad=0.4", fc="#eef2f5", ec="#4a6f8a", lw=1))

    fig.suptitle(
        "Sampling bias in the frozen corpus, measured against the full cancer census "
        f"({d['census_records']:,} articles)", fontsize=13.5, y=0.985)
    fig.tight_layout(rect=[0, 0, 1, 0.955])
    for ext in ("pdf", "png"):
        fig.savefig(FIG_DIR / f"fig28_census_capture.{ext}")
    plt.close(fig)
    print(f"fig28_census_capture: {len(rows)} mechanisms, {spread:.0f}x spread, "
          f"{vals[0]:.1f}:1 -> {vals[1]:.1f}:1 -> {vals[2]:.1f}:1"
          + (f"; dropped: {', '.join(dropped)}" if dropped else ""))
    return {"spread": spread, "ratios": vals, "dropped": dropped,
            "capture": {r["mechanism"]: r["capture"] for r in rows},
            # Reported so a guard can hold the encoding: a linear axis collapses
            # everything under 10% into the axis and destroys the finding.
            "xscale": ax.get_xscale()}


def _load(path):
    return json.loads(path.read_text(encoding="utf-8"))


def _corpus_year_counts():
    """Year histogram of the retrieved corpus, or None if it is absent.

    Read rather than typed. Fail-soft: without it fig2c draws the census alone
    and says so, because a contrast panel missing one of its two curves must
    not silently render as a single-series chart that looks complete.
    """
    if not FROZEN_INDEX.exists():
        return None
    counts = {}
    for line in FROZEN_INDEX.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        y = json.loads(line).get("year")
        if isinstance(y, int):
            counts[y] = counts.get(y, 0) + 1
    return counts


def fig2c_census_volume():
    """Census volume by year against the corpus's own, on twin axes.

    TWIN AXES ARE THE HONEST CHOICE HERE and also the risky one: two series
    four orders of magnitude apart cannot share a scale, but twin axes let a
    small series be drawn as tall as a large one. So the point of the panel is
    the SHAPE difference, and each axis is labelled with its own total to stop
    the heights being read against each other.
    """
    g = _load(GROWTH)
    field = {int(k): v for k, v in g["field_by_year"].items()}
    years = [y for y in sorted(field) if 1990 <= y <= g["end_year"]]
    corpus = _corpus_year_counts()

    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(years, [field[y] for y in years], color="#1565C0", linewidth=2.2,
            label=f"census (x{g['field_growth']} over "
                  f"{g['start_year']}-{g['end_year']})")
    ax.set_xlabel("year")
    ax.set_ylabel("census articles indexed", color="#1565C0")
    ax.tick_params(axis="y", labelcolor="#1565C0")
    ax.grid(alpha=0.25, linewidth=0.6)

    if corpus:
        cy = [y for y in years if y in corpus]
        ax2 = ax.twinx()
        ax2.plot(cy, [corpus[y] for y in cy], color="#E65100", linewidth=2.2,
                 linestyle="--", label="retrieved corpus (right axis)")
        ax2.set_ylabel("retrieved-corpus articles", color="#E65100")
        ax2.tick_params(axis="y", labelcolor="#E65100")
        lines = ax.get_lines() + ax2.get_lines()
        ax.legend(lines, [l.get_label() for l in lines], loc="upper left",
                  frameon=False, fontsize=9)
        note = ("Separate axes: the two series differ by orders of magnitude, so "
                "their HEIGHTS are not comparable and their SHAPES are.")
    else:
        ax.legend(loc="upper left", frameon=False, fontsize=9)
        note = ("The retrieved-corpus overlay is absent because "
                "corpus/INDEX.jsonl is not present.")

    ax.set_title("Publication volume: the census, and the retrieval that "
                 "looked like growth")
    fig.text(0.5, -0.03, note, ha="center", fontsize=8.5, style="italic",
             color="#455A64")
    fig.savefig(FIG_DIR / "fig2c_census_volume.pdf")
    fig.savefig(FIG_DIR / "fig2c_census_volume.png")
    plt.close(fig)
    print("  fig2c_census_volume")


def fig9c_design_composition():
    """Study-design classes, with the undetermined share DRAWN.

    It is the largest class. A chart showing only the classified records would
    imply the census assigns a design to everything, which is the single most
    misleading thing this panel could do.
    """
    d = _load(DESIGN)
    classes = d["classes"]
    order = ["trial", "clinical-other", "animal-model", "cell-culture",
             "animal-other", "non-primary", "undetermined"]
    labels = {"trial": "clinical trial", "clinical-other": "patient study,\nno trial type",
              "animal-model": "animal model", "cell-culture": "cell culture",
              "animal-other": "animal, no human\ncheck tag",
              "non-primary": "review, editorial,\ncomment", "undetermined": "undetermined"}
    colors = {"trial": "#2E7D32", "clinical-other": "#66BB6A",
              "animal-model": "#FB8C00", "cell-culture": "#FFB74D",
              "animal-other": "#FFE0B2", "non-primary": "#90A4AE",
              "undetermined": "#CFD8DC"}
    vals = [classes[k] for k in order]
    total = d["census"]

    fig, ax = plt.subplots(figsize=(9, 5))
    bars = ax.barh([labels[k] for k in order], vals,
                   color=[colors[k] for k in order], edgecolor="#37474F",
                   linewidth=0.7)
    ax.invert_yaxis()
    ax.set_xlabel(f"census articles (of {total:,})")
    ax.set_title("Study design, from NLM publication types and MeSH check tags")
    for b, v in zip(bars, vals):
        ax.text(b.get_width() + total * 0.008, b.get_y() + b.get_height() / 2,
                f"{v:,}  ({100 * v / total:.1f}%)", va="center", fontsize=9)
    ax.set_xlim(0, max(vals) * 1.28)
    ax.grid(axis="x", alpha=0.25, linewidth=0.6)
    trial_share_all = 100 * classes["trial"] / total
    trial_share_cls = 100 * classes["trial"] / d["classifiable"]
    fig.text(0.5, -0.04,
             f"Trial share has two denominators and both are shown: "
             f"{trial_share_all:.1f}% of the census, {trial_share_cls:.1f}% of "
             f"the {d['classifiable']:,} records carrying any design label. "
             f"The undetermined class is drawn, not dropped.",
             ha="center", fontsize=8.5, style="italic", color="#455A64")
    fig.savefig(FIG_DIR / "fig9c_design_composition.pdf")
    fig.savefig(FIG_DIR / "fig9c_design_composition.png")
    plt.close(fig)
    print("  fig9c_design_composition")


def fig14c_class_by_site():
    """Enrichment per site, physical against pharmacological.

    Plotted as enrichment rather than as counts because a count chart
    reproduces the ordering of the SITES, not of the modality. The 1.0 line is
    drawn because it is the only value that means anything on its own.
    """
    d = _load(SITES)
    rows = sorted(d["rows"], key=lambda r: -r["physical_enrichment"])
    sites = [r["site"] for r in rows]
    y = range(len(rows))
    fig, ax = plt.subplots(figsize=(9, 7))
    ax.barh([i - 0.2 for i in y], [r["physical_enrichment"] for r in rows],
            height=0.4, color="#E65100", label="physical class")
    ax.barh([i + 0.2 for i in y], [r["pharmacological_enrichment"] for r in rows],
            height=0.4, color="#1565C0", label="pharmacological class")
    ax.axvline(1.0, color="#37474F", linewidth=1.3)
    ax.set_yticks(list(y))
    ax.set_yticklabels(sites)
    ax.invert_yaxis()
    ax.set_xlabel("enrichment against the site's own share of site-assigned records")
    ax.set_title("Where each class of modality sits, by anatomical site")
    ax.legend(frameon=False, fontsize=9, loc="lower right")
    ax.grid(axis="x", alpha=0.25, linewidth=0.6)
    opposed = set(d["opposed_sites"])
    for i, r in enumerate(rows):
        if r["site"] in opposed:
            ax.get_yticklabels()[i].set_fontweight("bold")
    # WRAPPED. `bbox_inches="tight"` expands the canvas to fit any single-line
    # text, so an unwrapped footnote silently stretched this panel to a 2.2:1
    # aspect and shrank every bar.
    fig.text(0.5, -0.03,
             f"1.0 is the site's own weight. Bold labels mark the "
             f"{len(opposed)} sites where the two classes move in OPPOSITE\n"
             f"directions -- the reading that does not depend on how much a "
             f"site is written about.\nThe physical class holds "
             f"{len(d['physical_members'])} mechanisms and omits radiotherapy, "
             f"its largest real member.",
             ha="center", va="top", fontsize=8.5, style="italic",
             color="#455A64")
    fig.savefig(FIG_DIR / "fig14c_class_by_site.pdf")
    fig.savefig(FIG_DIR / "fig14c_class_by_site.png")
    plt.close(fig)
    print("  fig14c_class_by_site")


def fig15c_mechanism_pairs():
    """The ten most frequent mechanism pairs. Counts, never a rate."""
    d = _load(PROFILE)
    seen, pairs = set(), []
    for r in d["rows"]:
        for p in r["top_partners"]:
            key = tuple(sorted((r["mechanism"], p["mechanism"])))
            if key not in seen:
                seen.add(key)
                pairs.append((key, p["n"]))
    pairs.sort(key=lambda kv: -kv[1])
    top = pairs[:10]
    fig, ax = plt.subplots(figsize=(9, 5.5))
    labels = [f"{a} + {b}" for (a, b), _ in top]
    vals = [n for _, n in top]
    bars = ax.barh(labels, vals, color="#5E35B1", edgecolor="#311B92",
                   linewidth=0.7)
    ax.invert_yaxis()
    ax.set_xlabel("census articles carrying both mechanism descriptors")
    ax.set_title("Where mechanisms co-occur")
    for b, v in zip(bars, vals):
        ax.text(b.get_width() + max(vals) * 0.01,
                b.get_y() + b.get_height() / 2, f"{v:,}", va="center", fontsize=9)
    ax.set_xlim(0, max(vals) * 1.14)
    ax.grid(axis="x", alpha=0.25, linewidth=0.6)
    fig.text(0.5, -0.04,
             "Counts, not a rate. Co-tagging records that two vocabularies "
             "appear on one article -- not that two mechanisms were tested in "
             "combination -- and the co-occurrence RATE is a property of the "
             "labelling instrument rather than of the field.",
             ha="center", fontsize=8.5, style="italic", color="#455A64")
    fig.savefig(FIG_DIR / "fig15c_mechanism_pairs.pdf")
    fig.savefig(FIG_DIR / "fig15c_mechanism_pairs.png")
    plt.close(fig)
    print("  fig15c_mechanism_pairs")


def fig16c_trial_share():
    """Trial share against volume. The two axes are the point.

    Volume on a log axis because the mechanisms span nearly two orders of
    magnitude, and because a linear axis would put every small mechanism on the
    spine -- including the two whose trial share is highest.
    """
    d = _load(PROFILE)
    rows = [r for r in d["rows"] if r["census"] >= 200]
    fig, ax = plt.subplots(figsize=(9, 6))
    xs = [r["census"] for r in rows]
    ys = [r["trial_share"] for r in rows]
    ax.scatter(xs, ys, s=70, color="#00695C", edgecolor="#003D33", zorder=3)
    # Offset each label away from its nearest already-placed neighbour rather
    # than using one fixed offset: the two highest-share mechanisms sit within
    # a label's height of each other and collided.
    placed = []
    for r in sorted(rows, key=lambda r: -r["trial_share"]):
        x, y = r["census"], r["trial_share"]
        dx, dy = 7, 4
        for px, py, ptext in placed:
            # The threshold has to account for the LABEL's width, not the
            # point's position: at 0.22 decades `bispecific-antibody` still
            # ran into `antibody-drug-conjugate` 0.24 decades away, because a
            # long label reaches far past its own marker. Scale the exclusion
            # with the neighbouring label's length.
            reach = 0.14 + 0.016 * len(ptext)
            if abs(math.log10(x) - math.log10(px)) < reach and abs(y - py) < 0.6:
                dx, dy = -9, -13
                break
        ax.annotate(r["mechanism"], (x, y), textcoords="offset points",
                    xytext=(dx, dy), fontsize=8.5, color="#263238",
                    ha="left" if dx > 0 else "right")
        placed.append((x, y, r["mechanism"]))
    ax.set_xscale("log")
    ax.set_xlabel("census articles (log)")
    ax.set_ylabel("share carrying an NLM clinical-trial publication type (%)")
    ax.set_title("Maturity does not follow volume")
    ax.grid(alpha=0.25, linewidth=0.6)
    hifu = next((r for r in rows if r["mechanism"] == "hifu"), None)
    cart = next((r for r in rows if r["mechanism"] == "car-t"), None)
    if hifu and cart and hifu["trial_share"] > cart["trial_share"]:
        note = (f"HIFU sits at {hifu['trial_share']}% against CAR-T's "
                f"{cart['trial_share']}% on {cart['census'] / hifu['census']:.0f}x "
                f"the volume, so `physical modality` is not a maturity class.")
    else:
        note = ("Volume and trial share are plotted on separate axes because "
                "they do not track each other.")
    fig.text(0.5, -0.02, note, ha="center", fontsize=8.5, style="italic",
             color="#455A64")
    fig.savefig(FIG_DIR / "fig16c_trial_share.pdf")
    fig.savefig(FIG_DIR / "fig16c_trial_share.png")
    plt.close(fig)
    print("  fig16c_trial_share")


def fig1c_ratio_straddle():
    """The pharmacological:physical ratio under four class definitions.

    THE FIGURE EXISTS BECAUSE THE NUMBER IS NOT ONE NUMBER. The manuscript's
    central corpus claim was reported as a single ratio, recomputed on the
    census as a single larger one, and read as "the census understates the
    manuscript's case". Restricting BOTH classes symmetrically -- to
    descriptors naming a therapy rather than a process or a material -- gives
    figures BELOW the manuscript's, so under two readings of three the census
    overstates it instead.

    Plotting one bar per definition against the manuscript's line is the only
    honest presentation: a single bar would be picking one of four, and which
    one gets picked decides the direction of the conclusion.
    """
    d = _load(RATIO)
    lc = d["landscape_composition"]
    bars = [
        ("both curated classes\n(the headline)", lc["ratio"]),
        ("therapy-naming descriptors,\nrestriction's own criterion", lc["criterion_restored_ratio"]),
        ("therapy-naming, as the\nmaturity table defines it", lc["landscape_own_ratio"]),
        ("therapy-naming, intersected with\nthe curated pharmacological list", lc["precise_ratio"]),
    ]
    fig, ax = plt.subplots(figsize=(9.5, 5.5))
    labels = [b[0] for b in bars]
    vals = [b[1] for b in bars]
    colours = ["#1565C0" if v > MANUSCRIPT_RATIO else "#E65100" for v in vals]
    rects = ax.barh(labels, vals, color=colours, edgecolor="#263238", linewidth=0.7)
    ax.invert_yaxis()
    ax.axvline(MANUSCRIPT_RATIO, color="#B71C1C", linewidth=1.8, linestyle="--")
    # INSIDE the axes. Placed above the top bar in data coordinates it landed
    # on top of the title, because the bar axis is inverted and -0.72 is off
    # the top of the plot rather than below it.
    ax.text(MANUSCRIPT_RATIO + 0.2, len(bars) - 0.55,
            f"earlier keyword method: {MANUSCRIPT_RATIO}:1",
            color="#B71C1C", fontsize=9, va="center", ha="left")
    for r, v in zip(rects, vals):
        ax.text(r.get_width() + 0.25, r.get_y() + r.get_height() / 2,
                f"{v:.2f}:1", va="center", fontsize=9.5, fontweight="bold")
    ax.set_xlim(0, max(vals) * 1.16)
    ax.set_xlabel("pharmacological : physical, by census article count")
    ax.set_title("The volume ratio is a choice about class membership")
    ax.grid(axis="x", alpha=0.25, linewidth=0.6)
    above = sum(1 for v in vals if v > MANUSCRIPT_RATIO)
    fig.text(0.5, -0.06,
             f"Blue exceeds the earlier figure, orange falls below it: "
             f"{above} of {len(vals)} definitions above, {len(vals) - above} "
             f"below.\nThe DIRECTION of the imbalance survives every reading; "
             f"the claim that the census understates the earlier case does not.",
             ha="center", va="top", fontsize=8.5, style="italic", color="#455A64")
    fig.savefig(FIG_DIR / "fig1c_ratio_straddle.pdf")
    fig.savefig(FIG_DIR / "fig1c_ratio_straddle.png")
    plt.close(fig)
    print("  fig1c_ratio_straddle")


def fig5c_mechanism_site_matrix():
    """Mechanism by anatomical site, coloured by observed over expected.

    COLOURED BY RATIO, NOT BY COUNT, and the difference is the figure's whole
    argument. A count heatmap reproduces the product of the two marginals: the
    biggest mechanism crossed with the biggest site is always the brightest
    cell, so the picture shows which rows and columns are large and says
    nothing about where a literature is unusually thick or thin.

    Cells whose expectation falls below the interpretability floor are drawn in
    grey rather than coloured. A ratio computed on a handful of articles is
    noise, and colouring it would put the loudest colours on the least
    reliable cells -- the failure mode of every heatmap built from sparse
    counts.
    """
    import numpy as np
    from matplotlib.colors import TwoSlopeNorm

    d = _load(MATRIX)
    mechs = [m for m, _ in sorted(d["mechanism_totals"].items(),
                                  key=lambda kv: -kv[1])]
    sites = [s for s, _ in sorted(d["site_totals"].items(),
                                  key=lambda kv: -kv[1])]
    cell = {(r["mechanism"], r["site"]): r for r in d["rows"]}

    M = np.full((len(mechs), len(sites)), np.nan)
    for i, m in enumerate(mechs):
        for j, s in enumerate(sites):
            r = cell.get((m, s))
            if r and r["interpretable"] and r["ratio"]:
                M[i, j] = r["ratio"]

    fig, ax = plt.subplots(figsize=(11, 8))
    # Diverging around 1.0 -- the only value that means anything on its own --
    # on a log scale, because a ratio of 4 and a ratio of 1/4 are the same
    # size of departure and a linear scale would render them as 3 and 0.75.
    finite = M[np.isfinite(M)]
    hi = float(np.nanmax(finite)) if finite.size else 2.0
    norm = TwoSlopeNorm(vmin=np.log10(0.1), vcenter=0.0, vmax=np.log10(hi))
    im = ax.imshow(np.log10(M), cmap="RdBu_r", norm=norm, aspect="auto")
    ax.set_facecolor("#E0E0E0")

    ax.set_xticks(range(len(sites)))
    ax.set_xticklabels(sites, rotation=45, ha="right", fontsize=9)
    ax.set_yticks(range(len(mechs)))
    ax.set_yticklabels(mechs, fontsize=9)
    ax.set_title("Where each mechanism's literature is thicker or thinner "
                 "than its own marginals predict")
    for i in range(len(mechs)):
        for j in range(len(sites)):
            if np.isfinite(M[i, j]) and (M[i, j] >= 3 or M[i, j] <= 0.2):
                # White on the saturated ends: dark text on a dark cell is the
                # label that only the person who wrote it can read, and the
                # extreme cells are exactly the ones being labelled.
                dark = M[i, j] >= 0.75 * hi or M[i, j] <= 0.15
                ax.text(j, i, f"{M[i, j]:.1f}", ha="center", va="center",
                        fontsize=7.5,
                        color="#FFFFFF" if dark else "#212121")
    cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cb.set_label("observed / expected (log scale)", fontsize=9)
    cb.set_ticks([np.log10(v) for v in (0.1, 0.25, 0.5, 1, 2, 4, 8)])
    cb.set_ticklabels(["0.1", "0.25", "0.5", "1", "2", "4", "8"])
    fig.text(0.5, -0.02,
             f"Grey = expectation below {d['min_expected']:.0f} articles, where "
             f"a ratio describes a handful. Labels shown only at 3x or 0.2x.\n"
             f"Colour is NOT article count: a count heatmap reproduces the "
             f"product of the marginals and shows only which rows and columns "
             f"are large.",
             ha="center", va="top", fontsize=8.2, style="italic",
             color="#455A64")
    fig.savefig(FIG_DIR / "fig5c_mechanism_site_matrix.pdf", bbox_inches="tight")
    fig.savefig(FIG_DIR / "fig5c_mechanism_site_matrix.png", bbox_inches="tight")
    plt.close(fig)
    print("  fig5c_mechanism_site_matrix")


CENSUS_FIGURES = [
    (RATIO, fig1c_ratio_straddle),
    (MATRIX, fig5c_mechanism_site_matrix),
    (LANDSCAPE, fig28_census_capture),
    (GROWTH, fig2c_census_volume),
    (DESIGN, fig9c_design_composition),
    (SITES, fig14c_class_by_site),
    (PROFILE, fig15c_mechanism_pairs),
    (PROFILE, fig16c_trial_share),
]


def main() -> int:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    missing = sorted({p.name for p, _ in CENSUS_FIGURES if not p.exists()})
    if missing:
        print(f"missing committed artifacts: {', '.join(missing)}; run the "
              f"matching scripts/census_*.py or scripts/atlas_landscape.py",
              file=sys.stderr)
        return 1
    for _, fn in CENSUS_FIGURES:
        fn()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
