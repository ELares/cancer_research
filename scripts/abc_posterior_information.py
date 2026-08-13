#!/usr/bin/env python3
"""How much did the ABC posteriors actually learn? Measured against a null.

THE QUESTION
------------
`abc_joint_posterior.py` reports, per parameter, the fraction of the prior width
its 95% posterior interval still occupies, and flags a parameter "unconstrained"
when that fraction is at least 0.6. The threshold is a bare constant. It takes no
account of how many draws were accepted, and that is the thing that decides what
an UNINFORMATIVE posterior looks like.

With 30 accepted draws, the 2.5-97.5 percentile range of samples drawn from the
prior and nothing else covers about 0.90 of the prior width -- not 1.0, because
30 points rarely reach the corners. So 0.6 is far below anything noise produces,
and every parameter the run flagged was flagged for being WELL constrained.

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
It measures whether the data moved each parameter, not whether the model or the
fit is right. A parameter can be sharply constrained by a badly specified model.
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

# The threshold the generator uses today, kept here so the report can say what it
# is comparing against rather than describing it in prose.
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
    if not post or not n:
        return {"artifact": str(path.relative_to(PROJECT_ROOT)),
                "unassessable": "no posterior or no accepted count"}
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
        "artifact": str(path.relative_to(PROJECT_ROOT)),
        "n_draws": d.get("n_draws"), "n_accepted": n,
        "null_median_width": round(statistics.median(null), 3),
        "null_p5_width": round(_pct(null, 5), 3),
        "parameters": rows,
    }


def main() -> int:
    reports = [r for r in (assess(CAL / "joint-posterior.json"),
                           assess(CAL / "abc-posterior.json")) if r]
    if not reports:
        print("no ABC artifacts found")
        return 1

    L = ["# How much did the ABC posteriors actually learn?", "",
         "Generated by `scripts/abc_posterior_information.py`. Reads the committed",
         "artifacts; does not re-run the ABC.", "",
         "## The method", "",
         "Each ABC run reports, per parameter, the fraction of the prior width its",
         "95% posterior interval still occupies, and flags a parameter",
         f"*unconstrained* when that fraction is at least **{LEGACY_THRESHOLD}**.",
         "That threshold is a bare constant and takes no account of how many draws",
         "were accepted — which is the thing that decides what an uninformative",
         "posterior looks like.", "",
         "So the null is measured instead: draw `n_accepted` samples from the prior",
         "and nothing else, take the same 2.5/97.5 quantiles, and express the span",
         "as a fraction of the prior width. A parameter counts as **informed** only",
         "if its observed width falls below the 5th percentile of that null.", ""]

    for r in reports:
        if not r.get("parameters"):
            L += [f"## `{r['artifact']}`", "",
                  f"Not assessed: {r.get('unassessable', 'no parameters parsed')}.", ""]
            continue
        L += [f"## `{r['artifact']}`", "",
              f"{r['n_draws']} draws, **{r['n_accepted']} accepted**. With that many",
              f"draws an uninformative posterior still shows a median width of",
              f"**{r['null_median_width']}** of the prior, and only "
              f"{r['null_p5_width']} at its 5th percentile — so anything at or below",
              f"{r['null_p5_width']} is doing real work.", "",
              "| parameter | width / prior | null percentile | informed? | old flag |",
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
              f"by the data**: {', '.join(f'`{n}`' for n in informed) or 'none'}.", ""]
        if prior_only:
            L += [f"**{len(prior_only)} are indistinguishable from the prior**: "
                  + ", ".join(f"`{n}`" for n in prior_only)
                  + ". Their reported credible intervals are the prior's, and should",
                  "not be read as inferred.", ""]
        if mislabelled:
            L += [f"**The {LEGACY_THRESHOLD} threshold mislabels "
                  f"{len(mislabelled)} of them.** "
                  + ", ".join(f"`{n}`" for n in mislabelled)
                  + " are flagged *unconstrained* by the generator while sitting at",
                  "the very bottom of the uninformative null — they are among the",
                  "best-determined parameters in the run. The flag understates the",
                  "analysis's own result, and it does so for exactly the cascade",
                  "parameters the manuscript quotes credible intervals for.", ""]

    L += ["## What this does not say", "",
          "* It measures whether the data moved a parameter, not whether the model",
          "  or the fit is right. A parameter can be sharply constrained by a badly",
          "  specified model.",
          "* It does not rescue a posterior whose accepted draws are all worse than",
          "  a vector already in the repository, nor one whose medians produce an",
          "  inadmissible model — see `analysis/headline-at-fitted-cascade.md`.",
          "* The null assumes the accepted draws would otherwise be uniform on the",
          "  prior box, which is what these runs sample. It would need rederiving",
          "  for a non-uniform prior.", ""]

    OUT.write_text("\n".join(L) + "\n", encoding="utf-8")
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
              f"{mis} mislabelled by the {LEGACY_THRESHOLD} threshold")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
