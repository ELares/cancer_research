#!/usr/bin/env python3
"""How much did the ABC posteriors actually learn? Measured against a null.

THE QUESTION
------------
`abc_joint_posterior.py` formerly flagged a parameter "unconstrained" when its
95% posterior interval occupied at least 0.6 of the prior width. That legacy
threshold was a bare constant. It took no
account of how many draws were accepted, and that is the thing that decides what
an UNINFORMATIVE posterior looks like.

With 30 accepted draws, the 2.5-97.5 percentile range of samples drawn from the
prior and nothing else covers about 0.90 of the prior width -- not 1.0, because
30 points rarely reach the corners. The legacy cutoff could therefore label a
contracted interval unconstrained. The joint generator now uses a null threshold
matched to the accepted sample size; this report retains 0.6 for comparison.

THE TEST
--------
For each posterior, draw `n_accepted` uniform samples, take the same 2.5/97.5
quantiles, express the span as a fraction of the prior width, and repeat. That is
the null: the distribution of apparent constraint when the data says nothing. A
parameter counts as informed only if its observed width falls below what noise
produces at that sample size.

This is cheap and needs no compiled extension: it reads the committed artifacts
and does not re-run the ABC.

WHAT IT IS NOT
--------------
It measures marginal width contraction, not whether the model or the fit is
right. A parameter can be sharply constrained by a badly specified model; a
width test also misses shifts in location or changes in dependence.
Nor does it rescue a posterior whose accepted draws are all worse than a vector
already sitting in the repository -- that is a separate defect and is recorded in
`analysis/headline-at-fitted-cascade.md`.

Usage:
    python scripts/abc_posterior_information.py
"""

import json
import random
import statistics
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
CAL = PROJECT_ROOT / "analysis" / "calibration"
OUT = CAL / "abc-information-content.md"
OUT_JSON = CAL / "abc-information-content.json"

# Historical comparator only; the joint generator now uses a sample-size null.
LEGACY_THRESHOLD = 0.6
REPLICATES = 20000
SEED = 20260812


def null_widths(n_accepted: int, replicates: int = REPLICATES) -> list:
    """Apparent posterior width when the data is uninformative.

    Uniform draws on [0,1] stand in for any prior box: the width fraction is
    scale-free, so the null does not depend on the parameter's actual range.
    """
    rng = random.Random(SEED)
    out = []
    for _ in range(replicates):
        u = sorted(rng.random() for _ in range(n_accepted))
        out.append(_pct(u, 97.5) - _pct(u, 2.5))
    return out


def _pct(sorted_vals: list, p: float) -> float:
    """Linear-interpolation percentile, matching numpy's default."""
    if not sorted_vals:
        return 0.0
    k = (len(sorted_vals) - 1) * p / 100.0
    lo, hi = int(k), min(int(k) + 1, len(sorted_vals) - 1)
    return sorted_vals[lo] + (sorted_vals[hi] - sorted_vals[lo]) * (k - lo)


def assess(path: Path) -> dict:
    d = json.loads(path.read_text())
    post = d.get("posterior") or {}
    n = d.get("n_accepted")
    minimum = d.get("min_posterior")
    try:
        artifact = str(path.relative_to(PROJECT_ROOT))
    except ValueError:
        artifact = str(path)
    base = {"artifact": artifact, "n_draws": d.get("n_draws"), "n_accepted": n}
    if minimum is not None:
        base["min_posterior"] = minimum
    # A sample-size null cannot turn too few accepted draws into an assessable
    # posterior. Check the declared minimum even if its status flag is stale.
    if d.get("underpowered") or (minimum is not None and n is not None and n < minimum):
        return {**base, "underpowered": True,
                "unassessable": (
                    f"underpowered: {n} accepted draws; minimum required {minimum}"
                    if minimum is not None else
                    f"underpowered: {n} accepted draws; source marks the run underpowered")}
    if not post or not n:
        return {**base, "unassessable": "no posterior or no accepted count"}
    null = sorted(null_widths(n))
    # The #500 artifact stores the width fraction; the #332 one stores only the
    # quantiles plus its priors, so derive it there rather than skipping the
    # artifact silently -- an unparsed input reported as "0 of 0 informed" reads
    # exactly like a clean result.
    priors = d.get("priors") or {}
    rows = {}
    for name, v in post.items():
        w = v.get("posterior_width_frac_of_prior")
        if w is None:
            pr = priors.get(name)
            if isinstance(pr, dict):
                lo_, hi_ = pr.get("low"), pr.get("high")
            elif isinstance(pr, (list, tuple)) and len(pr) == 2:
                lo_, hi_ = pr
            else:
                lo_ = hi_ = None
            if lo_ is None or hi_ is None or hi_ <= lo_:
                continue
            w = round((v["q97_5"] - v["q2_5"]) / (hi_ - lo_), 3)
        # Where this width sits in the null: a high percentile means the width
        # is unremarkable for noise, i.e. the data did not move it.
        below = sum(1 for x in null if x < w)
        pctile = 100.0 * below / len(null)
        rows[name] = {
            "width_frac_of_prior": w,
            "null_percentile": round(pctile, 1),
            "informed": pctile <= 5.0,
            "legacy_flag_unconstrained": w >= LEGACY_THRESHOLD,
        }
    return {
        **base,
        "null_median_width": round(statistics.median(null), 3),
        "null_p5_width": round(_pct(null, 5), 3),
        "parameters": rows,
    }


def render(reports: list) -> str:
    L = ["# How much did the ABC posteriors actually learn?", "",
         "Generated by `scripts/abc_posterior_information.py`. Reads the committed",
         "artifacts; does not re-run the ABC.", "",
         "Runs marked underpowered, or below their declared minimum accepted",
         "count, are not assessed. Their diagnostic intervals cannot support",
         "parameter-information or contraction claims.", "",
         "## The method", "",
         "This audit reads or derives each parameter's 95% interval width as a",
         "fraction of its prior width. The former joint-generator rule flagged it",
         f"*unconstrained* at **{LEGACY_THRESHOLD}** or above. That legacy comparator",
         "is retained in the table; it is no longer the joint generator's rule.",
         "A fixed width cutoff takes no account of how many draws",
         "were accepted — which is the thing that decides what an uninformative",
         "posterior looks like.", "",
         "So the null is measured instead: draw `n_accepted` samples from the prior",
         "and nothing else, take the same 2.5/97.5 quantiles, and express the span",
         "as a fraction of the prior width. A parameter counts as **informed** only",
         "if its width is at or below the 5th percentile of that null. The current",
         "joint generator uses the same criterion with its own simulated null;",
         "small boundary differences can arise from Monte Carlo and rounding.", ""]

    for r in reports:
        if not r.get("parameters"):
            L += [f"## `{r['artifact']}`", "",
                  f"Not assessed: {r.get('unassessable', 'no parameters parsed')}.", ""]
            continue
        L += [f"## `{r['artifact']}`", "",
              f"{r['n_draws']} draws, **{r['n_accepted']} accepted**. With that many",
              f"draws an uninformative posterior still shows a median width of",
              f"**{r['null_median_width']}** of the prior, and only "
              f"{r['null_p5_width']} at its 5th percentile. The table applies this",
              "sample-size-dependent width test.", "",
              "| parameter | width / prior | null percentile | informed? | legacy 0.6 would flag |",
              "|---|--:|--:|---|---|"]
        for name, v in r["parameters"].items():
            L.append(f"| `{name}` | {v['width_frac_of_prior']:.3f} | "
                     f"{v['null_percentile']:.1f}% | "
                     f"{'**yes**' if v['informed'] else 'no'} | "
                     f"{'unconstrained' if v['legacy_flag_unconstrained'] else '—'} |")
        informed = [n for n, v in r["parameters"].items() if v["informed"]]
        prior_only = [n for n, v in r["parameters"].items() if not v["informed"]]
        mislabelled = [n for n, v in r["parameters"].items()
                       if v["informed"] and v["legacy_flag_unconstrained"]]
        L += ["",
              f"**{len(informed)} of {len(r['parameters'])} parameters are informed "
              f"under this width criterion**: {', '.join(f'`{n}`' for n in informed) or 'none'}.", ""]
        if prior_only:
            L += [f"**{len(prior_only)} show no detected marginal width contraction**: "
                  + ", ".join(f"`{n}`" for n in prior_only)
                  + ". This does not establish that their full distributions equal",
                  "the prior or that they have no effect on model predictions.", ""]
        if mislabelled:
            L += [f"**The legacy {LEGACY_THRESHOLD} threshold would mislabel "
                  f"{len(mislabelled)} of them.** "
                  + ", ".join(f"`{n}`" for n in mislabelled)
                  + " would be flagged *unconstrained* despite passing this null-based",
                  "contraction check. This is a comparison with the retired flag,",
                  "not a claim that the current joint generator still applies it.", ""]

    L += ["## What this does not say", "",
          "* It measures marginal width contraction, not whether the model",
          "  or the fit is right. A parameter can be sharply constrained by a badly",
          "  specified model. Failing this width test does not rule out shifts in",
          "  location or changes in parameter dependence.",
          "* The separate `analysis/headline-at-fitted-cascade.md` substitution test",
          "  found recorded historical vectors inadmissible in the spatial model.",
          "  It did not test the corrected posterior; this width analysis also makes",
          "  no admissibility claim about that new fit.",
          "* The null assumes the accepted draws would otherwise be uniform on the",
          "  prior box, which is what these runs sample. It would need rederiving",
          "  for a non-uniform prior.", ""]

    return "\n".join(L) + "\n"


def main() -> int:
    reports = [r for r in (assess(CAL / "joint-posterior.json"),
                           assess(CAL / "abc-posterior.json")) if r]
    if not reports:
        print("no ABC artifacts found")
        return 1
    OUT.write_text(render(reports), encoding="utf-8")
    OUT_JSON.write_text(json.dumps(reports, indent=2, sort_keys=True) + "\n")
    print(f"wrote {OUT}\nwrote {OUT_JSON}")
    for r in reports:
        if not r.get("parameters"):
            print(f"  {r['artifact']}: NOT ASSESSED "
                  f"({r.get('unassessable', 'no parameters parsed')})")
            continue
        inf = sum(1 for v in r["parameters"].values() if v["informed"])
        mis = sum(1 for v in r["parameters"].values()
                  if v["informed"] and v["legacy_flag_unconstrained"])
        print(f"  {r['artifact']}: {inf}/{len(r['parameters'])} informed, "
              f"{mis} would be mislabelled by the legacy {LEGACY_THRESHOLD} threshold")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
