#!/usr/bin/env python3
"""Compare descriptor agreement with symmetric modality text rules.

Recall and precision use title/abstract matches as the reference, not clinical
or manually adjudicated ground truth. Reports derive their conclusions from
raw counts, including equal, reversed, zero and unavailable comparisons.
Historical log-ratio intervals retain their original approximation: descriptor
and text totals share a within-arm covariance term, but cross-arm overlap is
not recorded and is not corrected. These are model intervals, not estimates
of error from sampling an otherwise completely enumerated census.

The indexed census root follows FERRO_ATLAS_ROOT. Missing or unreadable input
must fail before either report is replaced. A readable stream with no matching
subject or modality is a valid result. --render-only requires only saved counts
and never rewrites them.

Usage:
    python scripts/atlas_descriptor_recall.py
    python scripts/atlas_descriptor_recall.py --render-only
"""

import argparse
from copy import deepcopy
import json
import math
import re
import statistics
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from census_input import atlas_root, iter_census_records  # noqa: E402

PROJECT_ROOT = Path(__file__).resolve().parent.parent
ATLAS = atlas_root()
OUT_MD = PROJECT_ROOT / "analysis" / "atlas-descriptor-recall.md"
OUT_JSON = PROJECT_ROOT / "analysis" / "atlas-descriptor-recall.json"

SUBJECT = "ferroptosis"

# One shape of text rule for every arm: the full word, its hyphenation, and the
# conventional acronym. Asymmetry between arms here would reproduce the defect
# this analysis exists to correct, so the shape is held fixed and only the stem
# changes.
ARMS = {
    "PDT": {
        "label": "photodynamic therapy",
        "descriptors": {"photochemotherapy"},
        "text": r"photodynamic|photo-dynamic|\bPDT\b",
    },
    "SDT": {
        "label": "sonodynamic therapy",
        "descriptors": {"ultrasonic therapy"},
        "text": r"sonodynamic|sono-dynamic|\bSDT\b",
    },
}
# The pair the manuscript's Section 8.2 states, and its stated value.
RATIO_PAIR = ("PDT", "SDT")
MANUSCRIPT_RATIO = 2.93
# The interval level, interpolated wherever it is printed. A "95%" typed into
# prose is still a number a sentence can outlive, and this file's own subject
# is exactly that failure.
CONF = 0.95
# DERIVED from CONF, never typed beside it. Hardcoding Z let a
# mutation print a 68% interval under a "95% CI" label, because the
# label was interpolated and the arithmetic was not.
Z = statistics.NormalDist().inv_cdf(1 - (1 - CONF) / 2)


def scan() -> dict:
    pats = {k: re.compile(v["text"], re.I) for k, v in ARMS.items()}
    # the descriptor name travels WITH the counts, so a consumer quoting this
    # artifact names the right descriptor instead of falling back to the arm
    # key and printing "`PDT` recalls 80.2% of PDT papers"
    stat = {k: {"text": 0, "descriptor": 0, "both": 0,
                "descriptors": sorted(v["descriptors"]),
                "label": v["label"],
                "descriptor_not_text": []} for k, v in ARMS.items()}
    n = 0
    for r in iter_census_records(ATLAS / "records"):
        # Validate before filtering: malformed non-subject records cannot make
        # a corrupt stream look like a valid subject with no matches.
        if not isinstance(r, dict):
            raise ValueError("Census record must be a JSON object")
        mesh_values = r.get("mesh")
        if mesh_values is not None and (
            not isinstance(mesh_values, list)
            or any(not isinstance(value, str) for value in mesh_values)
        ):
            raise ValueError("Census mesh must be a list of strings or null")
        for field in ("title", "abstract"):
            if r.get(field) is not None and not isinstance(r[field], str):
                raise ValueError(f"Census {field} must be a string or null")
        mesh = {m.lower() for m in (mesh_values or [])}
        if SUBJECT not in mesh:
            continue
        n += 1
        blob = (r.get("title") or "") + " " + (r.get("abstract") or "")
        for k, arm in ARMS.items():
            t = bool(pats[k].search(blob))
            d = bool(mesh & arm["descriptors"])
            stat[k]["text"] += t
            stat[k]["descriptor"] += d
            if t and d:
                stat[k]["both"] += 1
            elif d and not t and len(stat[k]["descriptor_not_text"]) < 12:
                stat[k]["descriptor_not_text"].append(
                    {"pmid": r.get("pmid"),
                     "title": (r.get("title") or "")[:110]})

    return assemble({
        "subject": SUBJECT,
        "subject_articles": n,
        "arms": stat,
        "ratio_pair": list(RATIO_PAIR),
        "manuscript_ratio": MANUSCRIPT_RATIO,
    })


def assemble(d: dict) -> dict:
    """Recompute derived results from raw counts without mutating the input."""
    d = deepcopy(d)
    stat = d["arms"]
    n = d["subject_articles"]
    if not isinstance(stat, dict) or set(stat) != set(ARMS):
        raise ValueError("arms must contain the configured modality pair")
    pair = d["ratio_pair"]
    if (not isinstance(pair, (list, tuple)) or len(pair) != 2
            or any(not isinstance(name, str) or name not in stat for name in pair)
            or pair[0] == pair[1]):
        raise ValueError("ratio_pair must name distinct known arms")
    manuscript_ratio = d["manuscript_ratio"]
    if (type(manuscript_ratio) not in (int, float)
            or not math.isfinite(manuscript_ratio) or manuscript_ratio <= 0):
        raise ValueError("manuscript_ratio must be finite and positive")
    if type(n) is not int or n < 0:
        raise ValueError("subject_articles must be a nonnegative integer")
    for name, counts in stat.items():
        for field in ("text", "descriptor", "both"):
            value = counts[field]
            if type(value) is not int or not 0 <= value <= n:
                raise ValueError(f"{name} {field} must be an integer within subject_articles")
        if counts["both"] > min(counts["text"], counts["descriptor"]):
            raise ValueError(f"{name} overlap exceeds its marginal counts")
        if counts["text"] + counts["descriptor"] - counts["both"] > n:
            raise ValueError(f"{name} union exceeds subject_articles")

    for k, s in stat.items():
        s["recall"] = s["both"] / s["text"] if s["text"] else None
        s["precision"] = s["both"] / s["descriptor"] if s["descriptor"] else None

    a, b = d["ratio_pair"]
    by_desc = (stat[a]["descriptor"] / stat[b]["descriptor"]
               if stat[b]["descriptor"] else None)
    by_text = stat[a]["text"] / stat[b]["text"] if stat[b]["text"] else None
    recalls = [s["recall"] for s in stat.values() if s["recall"] is not None]

    # Retain historical deterministic log-ratio and Wilson approximations.
    def _logratio_ci(x, y):
        """Interval for x/y treating both as Poisson counts."""
        if not x or not y:
            return None
        se = math.sqrt(1 / x + 1 / y)
        pt = math.log(x / y)
        return [math.exp(pt - Z * se), math.exp(pt + Z * se)]

    def _wilson(k, n):
        if not n:
            return None
        p, z = k / n, Z
        d = 1 + z * z / n
        c = p + z * z / (2 * n)
        m = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
        return [(c - m) / d, (c + m) / d]

    for k, s in stat.items():
        s["recall_ci"] = _wilson(s["both"], s["text"])
        s["precision_ci"] = _wilson(s["both"], s["descriptor"])

    text_ci = _logratio_ci(stat[a]["text"], stat[b]["text"])
    desc_ci = _logratio_ci(stat[a]["descriptor"], stat[b]["descriptor"])
    # Relative movement between descriptor and text counting routes.
    inflation = (by_desc / by_text
                 if by_desc is not None and by_text is not None and by_text > 0
                 else None)
    # Historical approximation: paired descriptor/text counts WITHIN each
    # arm. Cross-arm overlaps are not retained, so their covariance is absent.
    # Preserve the estimator for historical comparability; report its scope.
    var = 0.0
    for k in (a, b):
        D, T, B = stat[k]["descriptor"], stat[k]["text"], stat[k]["both"]
        if not (D and T):
            var = None
            break
        var += 1 / D + 1 / T - 2 * B / (D * T)
    infl_se = math.sqrt(var) if var and var > 0 else None
    infl_ci = ([math.exp(math.log(inflation) - Z * infl_se),
                math.exp(math.log(inflation) + Z * infl_se)]
               if (inflation and infl_se) else None)
    # Descriptive max/min recall ratio, not a test of recall parity.
    ra = (max(recalls) / min(recalls)) if len(recalls) > 1 and min(recalls) else None
    d.update({
        "arms": stat,
        "ratio_by_descriptor": by_desc,
        "ratio_by_descriptor_ci": desc_ci,
        "ratio_by_text": by_text,
        "ratio_by_text_ci": text_ci,
        "recall_asymmetry": ra,
        "descriptor_inflation": inflation,
        "descriptor_inflation_ci": infl_ci,
        "symmetric_ratio_covers_manuscript": bool(
            text_ci and text_ci[0] <= d["manuscript_ratio"] <= text_ci[1]),
        "descriptor_route_significantly_inflates": bool(
            infl_ci and infl_ci[0] > 1.0),
    })
    return d


def _percent(value):
    return "n/a" if value is None else f"{100 * value:.1f}%"


def _ratio(value):
    return "n/a" if value is None else f"{value:.2f}:1"


def render(d: dict) -> str:
    d = assemble(d)
    a, b = d["ratio_pair"]
    A, B = d["arms"][a], d["arms"][b]

    def ci(value):
        return f"[{value[0]:.2f}, {value[1]:.2f}]" if value is not None else "n/a"

    L = ["# Modality descriptor agreement with title and abstract text", "",
         "*Generated by `scripts/atlas_descriptor_recall.py` from raw counts. "
         "Render-only regeneration does not rescan the census.*", "",
         f"Over the {d['subject_articles']:,} census articles indexed "
         f"`{d['subject'].title()}`, measured with one text rule applied "
         "identically to both arms:", "",
         "| arm | descriptor | text says so | descriptor | both | recall | precision |",
         "|---|---|--:|--:|--:|--:|--:|"]
    for k in (a, b):
        s = d["arms"][k]
        descriptors = ", ".join(f"`{name.title()}`" for name in s["descriptors"])
        L.append(f"| {k} | {descriptors} | {s['text']:,} | "
                 f"{s['descriptor']:,} | {s['both']:,} | "
                 f"**{_percent(s['recall'])}** | {_percent(s['precision'])} |")
    L += ["", "Here **recall** is the fraction of text matches carrying the "
          "descriptor; **precision** is the fraction of descriptor matches "
          "also matching the text rule. Both measure **text agreement**, "
          "not accuracy against independently adjudicated ground truth. "
          "`n/a` marks an undefined denominator or an unavailable interval; "
          "a measured zero remains zero.", ""]
    if not d["subject_articles"]:
        L += ["The census was readable, but no article matched the subject "
              "descriptor. Modality rates and ratios are unavailable for "
              "this empty subject subset.", ""]

    recalls = [s["recall"] for s in (A, B)]
    if any(value is None for value in recalls):
        L += ["The available text matches do not define recall for both arms; "
              "a recall difference cannot be assessed.", ""]
    elif recalls[0] == recalls[1]:
        L += [f"The observed recalls are equal at {_percent(recalls[0])}; "
              "these counts show no recall asymmetry.", ""]
    else:
        higher, lower = (a, b) if recalls[0] > recalls[1] else (b, a)
        L += [f"Observed text-agreement recall is higher for {higher} than "
              f"for {lower}. " + (
                  f"The recalls differ by a factor of **{d['recall_asymmetry']:.2f}**."
                  if d["recall_asymmetry"] is not None else
                  "The lower recall is zero, so their multiplicative gap is undefined."
              ) + " This is a descriptive comparison, not a test of recall parity.", ""]

    L += ["## What that does to the ratio the manuscript states", "",
          f"The ratio below is {a}:{b}.", "",
          f"| how the two are counted | ratio | {100*CONF:g}% approximate CI |",
          "|---|--:|--:|",
          "| by descriptor (the route used in `manuscript-vs-census.md`) | "
          f"**{_ratio(d['ratio_by_descriptor'])}** | {ci(d['ratio_by_descriptor_ci'])} |",
          f"| by one text rule applied to both | **{_ratio(d['ratio_by_text'])}** "
          f"| {ci(d['ratio_by_text_ci'])} |",
          f"| as the manuscript states it | {_ratio(d['manuscript_ratio'])} | |", ""]
    tc, dc = d["ratio_by_text_ci"], d["ratio_by_descriptor_ci"]
    if tc is not None and dc is not None:
        overlap = max(tc[0], dc[0]) <= min(tc[1], dc[1])
        L += ["The marginal intervals " + ("overlap" if overlap else "do not overlap") +
              ". Their overlap alone is not a test of the difference between "
              "routes, since both routes count the same articles.", ""]

    L += ["### What the data supports", ""]
    if tc is None:
        L += ["The text-route interval is unavailable, so this approximation "
              "cannot assess whether the manuscript understates or overstates "
              "the text-matched ratio.", ""]
    elif d["symmetric_ratio_covers_manuscript"]:
        L += [f"The manuscript's {_ratio(d['manuscript_ratio'])} sits **inside** "
              f"the symmetric interval {ci(tc)}, so this model cannot distinguish "
              "the census text ratio from the manuscript's ratio. It does not "
              "establish agreement. An UNDERSTATEMENT verdict is not supported "
              "by this comparison.", ""]
    elif tc[0] > d["manuscript_ratio"]:
        L += [f"The symmetric interval {ci(tc)} lies above the manuscript's "
              f"{_ratio(d['manuscript_ratio'])}. Under this model the manuscript "
              "understates the text-matched ratio in this census subset.", ""]
    else:
        L += [f"The symmetric interval {ci(tc)} lies below the manuscript's "
              f"{_ratio(d['manuscript_ratio'])}. Under this model the manuscript "
              "overstates the text-matched ratio in this census subset.", ""]

    inflation, ic = d["descriptor_inflation"], d["descriptor_inflation_ci"]
    if inflation is None:
        L += ["The descriptor/text ratio multiplier is undefined because the "
              "required ratios cannot both support division. No directional "
              "claim about descriptor-route inflation follows.", ""]
    elif inflation == 1:
        L += ["The descriptor and text routes give the same observed ratio "
              f"(multiplier **{inflation:.2f}x**); there is no observed "
              "descriptor-route inflation or reduction.", ""]
    else:
        direction = "higher" if inflation > 1 else "lower"
        L += [f"The descriptor route gives a {direction} point ratio than the "
              f"text route, with multiplier **{inflation:.2f}x**.", ""]
    if ic is None:
        L += ["The multiplier interval is unavailable when a count is zero "
              "or the estimated log variance is zero. This does not establish "
              "population equality or a statistically resolved difference.", ""]
    elif d["descriptor_route_significantly_inflates"]:
        L += [f"The multiplier interval {ci(ic)} excludes parity above it: "
              "the descriptor route inflates this ratio relative to the text "
              "route under the stated approximation.", ""]
    elif ic[1] < 1:
        L += [f"The multiplier interval {ci(ic)} excludes parity below it: "
              "the descriptor route reduces this ratio relative to the text "
              "route under the stated approximation.", ""]
    else:
        L += [f"The multiplier interval {ci(ic)} includes parity; a directional "
              "difference between routes is not resolved by this approximation.", ""]

    L += ["## Descriptor-only and text-only matches", ""]
    for k, s in ((a, A), (b, B)):
        descriptor_only = s["descriptor"] - s["both"]
        text_only = s["text"] - s["both"]
        relation = ("fewer" if text_only > descriptor_only else
                    "more" if text_only < descriptor_only else "as many")
        comparison = "as" if relation == "as many" else "than"
        L += [f"For {k}, {descriptor_only:,} records match only the descriptor "
              f"and {text_only:,} match only the text rule. The descriptor "
              f"therefore counts {relation} records {comparison} the text rule.", ""]
    L += ["These disagreements do not identify false positives or false "
          "negatives without adjudication. A descriptor can correctly label "
          "a paper whose title and abstract use other words; a text match can "
          "also be off-topic. Neither route supplies upper or lower bounds "
          "on the true modality count or ratio.", ""]
    examples = B.get("descriptor_not_text", [])
    if examples:
        L += [f"Examples of {b} descriptor matches without a text match "
              "(the stored list may be truncated):", ""]
        L += [f"* {r['pmid']} — {r['title']}" for r in examples]
        L += [""]

    L += ["## Interpretation and uncertainty", "",
          "A precision comparison alone cannot establish equal recall. "
          "Changing descriptor sets also does not by itself validate their "
          "agreement with a text rule. Symmetry of the regex structure does "
          "not establish equal sensitivity or specificity for the modalities.", "",
          f"The {100*CONF:g}% intervals describe a statistical model, not "
          "sampling error from observing only part of this enumerated census. "
          "Marginal log-ratio intervals treat arm counts as independent "
          "Poisson variables. The descriptor/text multiplier includes "
          "descriptor-text covariance within each arm, but does not account "
          "for cross-arm overlap: those joint counts are absent from the "
          "saved artifact. It is not a complete article-level paired "
          "uncertainty analysis. Wilson intervals in the JSON likewise "
          "assume a binomial model for text agreement. Sparse or boundary "
          "counts can make these approximations unreliable; unavailable "
          "intervals are not evidence of equivalence. None of these intervals "
          "captures taxonomy, retrieval, text-rule or indexing error.", "",
          "## What this does not claim", "",
          "* Not that the manuscript is wrong or confirmed: the comparison "
          "concerns text and descriptor counts under the stated model, not "
          "the underlying therapeutic literature's true ratio.",
          "* Not that the text rule is ground truth. It has errors in both "
          "directions and is applied identically to both arms.",
          "* Not that either descriptor is intrinsically inaccurate: "
          "descriptor/text disagreements require independent review.", ""]
    return "\n".join(L) + "\n"


def _write_reports(outputs):
    """Stage every output before replacing any destination.

    Replacement is atomic per file, not a transaction across the report pair.
    """
    staged = []
    try:
        for destination, content in outputs:
            with tempfile.NamedTemporaryFile(
                mode="w", encoding="utf-8", dir=destination.parent,
                prefix=f".{destination.name}.", delete=False,
            ) as stream:
                temporary = Path(stream.name)
                staged.append((temporary, destination))
                stream.write(content)
            # Keep an existing report's permissions when replacing it.
            if destination.exists():
                temporary.chmod(destination.stat().st_mode & 0o777)
        for temporary, destination in staged:
            temporary.replace(destination)
    finally:
        for temporary, _ in staged:
            temporary.unlink(missing_ok=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--render-only", action="store_true")
    args = ap.parse_args()
    if args.render_only:
        d = assemble(json.loads(OUT_JSON.read_text(encoding="utf-8")))
    else:
        d = scan()
    markdown = render(d)
    outputs = []
    if not args.render_only:
        outputs.append((OUT_JSON, json.dumps(
            d, indent=1, sort_keys=True, allow_nan=False) + "\n"))
    outputs.append((OUT_MD, markdown))
    _write_reports(outputs)
    for destination, _ in outputs:
        print(f"wrote {destination}")
    for k, s in d["arms"].items():
        print(f"  {k:4s} recall {_percent(s['recall'])}  "
              f"precision {_percent(s['precision'])}")
    print(f"  ratio by descriptor {_ratio(d['ratio_by_descriptor'])}  "
          f"by text {_ratio(d['ratio_by_text'])}  "
          f"manuscript {_ratio(d['manuscript_ratio'])}")


if __name__ == "__main__":
    main()
