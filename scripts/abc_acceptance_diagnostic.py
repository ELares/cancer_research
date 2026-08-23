#!/usr/bin/env python3
"""Why the joint ABC accepts parameter vectors worse than one already committed.

THE OBSERVATION
---------------
`analysis/calibration/joint-posterior.json` reports a posterior over 30 accepted
draws with an acceptance threshold of 0.35 in the joint distance. But the vector
this repository had already committed -- the #330 cascade plus the #502
shared-switch erastin parameters -- scores 0.2202 on that same distance, and the
posterior MEDIAN scores 0.2413. The inference returned something worse than its
own starting material.

WHAT IT IS NOT (both tested here, not assumed)
----------------------------------------------
1. A TRUNCATED PRIOR. The committed vector sits exactly on two prior bounds --
   `k_erastin` at its low bound of 3.0 and `hill` at its high bound of 6.0 --
   which looks like the box clipping the optimum. It is not: pushing `k_erastin`
   below 3.0 makes the fit monotonically worse, and `hill` is INERT, changing the
   distance by nothing at all between 6 and 10. The box is not in the way.

2. AN UNDER-SAMPLED POSTERIOR. More draws would help, but the defect is not that
   30 is a small number. It is what the acceptance rule does with any number.

WHAT IT IS
----------
The acceptance is a fixed FRACTION -- `n_accept = n_draws * 0.02` -- so the run
always accepts its best 2% no matter how bad they are. The reported epsilon is
therefore an OUTPUT, whatever the 30th-best draw happened to score, and never a
criterion anything had to meet. A rejection ABC built this way has no floor: hand
it draws that are uniformly terrible and it will still return 2% of them and
label the result a posterior.

Here that shows up as a measurable gap. Uniform draws essentially never reach the
good region in 7 dimensions -- 0 of 300 beat the committed vector in the run
below -- so the 2% cut lands in a shell well above what is achievable, and the
posterior median inherits the shell rather than the optimum.

Both things are true at once, and the second does not cancel the first: the
accepted set IS narrower than the prior for five of seven parameters
(`analysis/calibration/abc-information-content.md`), so the data does move them.
It is centred in the wrong place.

Usage:
    python scripts/abc_acceptance_diagnostic.py            # 300 draws, ~10 s
    python scripts/abc_acceptance_diagnostic.py --draws 1500
"""

import argparse
import importlib.util
import json
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
ART = PROJECT_ROOT / "analysis" / "calibration" / "joint-posterior.json"
OUT_MD = PROJECT_ROOT / "analysis" / "calibration" / "abc-acceptance-diagnostic.md"
OUT_JSON = PROJECT_ROOT / "analysis" / "calibration" / "abc-acceptance-diagnostic.json"

# The vector already in the repository when the ABC ran: the #330 CTRPv2 cascade
# fit, plus #502's shared-switch erastin parameters, with the remaining two at
# their Params::default() values (the ABC's own _shared() reads all four).
# The run that EXHIBITED the defect, kept as a record. These cannot be recomputed
# -- the artifact they describe has been replaced by the fixed run -- so they are
# stated as history and the guards check them as constants, not as measurements.
HISTORICAL = {
    "n_draws": 1500, "n_accepted": 30,
    "epsilon": 0.35, "posterior_median_distance": 0.2413,
    "rule": "fixed 2% quantile",
}

COMMITTED = {"lp_propagation": 0.7, "lp_rate": 0.4, "gpx4_rate": 0.30,
             "gsh_scav_efficiency": 0.5, "k_um": 0.25,
             "k_erastin": 3.0, "hill": 6.0}


def _abc():
    spec = importlib.util.spec_from_file_location(
        "abc_joint_posterior", PROJECT_ROOT / "scripts" / "abc_joint_posterior.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def _roundtrip(d: dict) -> dict:
    """Render from what the artifact WILL contain, not from the live dict.

    The JSON is written with `sort_keys=True`, so a dict rendered in insertion
    order produces a document that can never be reproduced from its own
    artifact -- the row ordering differs.

    ROUND-TRIPPING IS NOT ENOUGH ON ITS OWN, and assuming it was regressed a
    published finding here: any ordering that CARRIED MEANING has to be
    re-established inside the renderer, because sorting the input replaces a
    rank order with an alphabetical one. Every table below that had a
    meaningful order now sorts explicitly.
    """
    return json.loads(json.dumps(d, sort_keys=True))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--draws", type=int, default=300)
    args = ap.parse_args()

    try:
        import numpy as np
    except ImportError:
        print("numpy required", file=sys.stderr)
        return 1
    sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
    abc = _abc()
    ck = abc.ck

    art = json.loads(ART.read_text())
    c = art["curves"]
    rd, ed = c["rsl3_doses_um"], c["erastin_doses_um"]
    er, ee = c["empirical_rsl3_ml162"], c["empirical_erastin"]

    def dist(p):
        return (ck.rmse(abc.model_rsl3(rd, p), er)
                + ck.rmse(abc.model_erastin(ed, p), ee))

    committed = dist(COMMITTED)
    median_vec = {k: v["median"] for k, v in art["posterior"].items()}
    median_d = dist(median_vec)
    eps = art.get("epsilon_joint_distance")

    # 1. Is the prior truncating? Step outside the two bounds the committed
    #    vector sits on and see whether the fit improves.
    outside = {"k_erastin": {}, "hill": {}}
    for v in (3.0, 2.0, 1.0):
        outside["k_erastin"][v] = dist(dict(COMMITTED, k_erastin=v))
    for v in (6.0, 8.0, 10.0):
        outside["hill"][v] = dist(dict(COMMITTED, hill=v))

    # 2. How often does uniform sampling reach the good region at all?
    names = [p[0] for p in abc.PRIORS]
    lows = np.array([p[1] for p in abc.PRIORS])
    highs = np.array([p[2] for p in abc.PRIORS])
    rng = np.random.default_rng(7)
    ds = np.array([dist(dict(zip(names, rng.uniform(lows, highs))))
                   for _ in range(args.draws)])

    res = {
        "committed_vector": COMMITTED,
        "committed_distance": round(float(committed), 4),
        "posterior_median_distance": round(float(median_d), 4),
        "reported_epsilon": eps,
        "epsilon_excess_over_committed": round(float(eps / committed - 1.0), 4) if eps else None,
        "prior_truncation_test": {k: {str(a): round(float(b), 4) for a, b in v.items()}
                                  for k, v in outside.items()},
        "n_accepted_now": art.get("n_accepted"),
        "sampling": {
            "draws": args.draws,
            "draws_of_record": art.get("n_draws"),
            "best": round(float(ds.min()), 4),
            "quantile_2pct": round(float(np.quantile(ds, 0.02)), 4),
            "n_beating_committed": int((ds < committed).sum()),
            "frac_inside_epsilon": round(float((ds < eps).mean()), 4) if eps else None,
        },
    }
    OUT_JSON.write_text(json.dumps(res, indent=2, sort_keys=True) + "\n")
    OUT_MD.write_text(render(_roundtrip(res)), encoding="utf-8")
    print(f"wrote {OUT_MD}\nwrote {OUT_JSON}")
    print(f"  committed {res['committed_distance']}  median "
          f"{res['posterior_median_distance']}  eps {eps}")
    print(f"  {res['sampling']['n_beating_committed']}/{args.draws} draws beat the "
          "committed vector")
    return 0


# The bounds the committed vector sits on. Each sweep steps OUTWARD from its
# bound (`k_erastin` 3->2->1, `hill` 6->8->10), which is the order the prose
# below describes, so the table has to lead with the bound.
BOUNDS = {"k_erastin": 3.0, "hill": 6.0}


def _outward(d: dict, param: str) -> list:
    """(value, distance) pairs ordered OUTWARD FROM THE BOUND.

    Three orderings are possible here and two are wrong. The swept values are
    dict KEYS, so JSON makes them strings and a round-tripped artifact renders
    them lexicographically -- 10.0, 6.0, 8.0 -- directly above prose about what
    happens between 6 and 10. Sorting numerically ASCENDING fixes `hill` and
    silently reverses `k_erastin`, putting `(prior low bound)` on the last row
    under a sentence that says "pushing below its bound", which is what the
    first attempt at this shipped.

    Distance from the bound is the sequence the sweep actually walks.
    """
    b = BOUNDS[param]
    return sorted(d.items(), key=lambda kv: abs(float(kv[0]) - b))


def render(r: dict) -> str:
    s = r["sampling"]
    ke = r["prior_truncation_test"]["k_erastin"]
    hl = r["prior_truncation_test"]["hill"]
    hill_inert = len(set(hl.values())) == 1
    ke_worse = all(v >= ke["3.0"] for v in ke.values())
    fixed = r["posterior_median_distance"] < r["committed_distance"]
    h = HISTORICAL

    L = ["# The ABC acceptance rule: the defect, and its fix", "",
         "Generated by `scripts/abc_acceptance_diagnostic.py`.", ""]

    if fixed:
        L += ["## Status: RESOLVED", "",
              "The acceptance rule has been changed from a fixed quantile to a",
              "tolerance anchored to a reachable distance, and the posterior it",
              "produces now beats the reference it used to lose to. This document",
              "records what was wrong, how it was diagnosed, and what the fix",
              "bought — the defect itself is no longer present.", "",
              "| | before | after |", "|---|--:|--:|",
              f"| acceptance rule | {h['rule']} | tolerance |",
              f"| draws | {h['n_draws']:,} | {r['sampling']['draws_of_record']:,} |",
              f"| accepted | {h['n_accepted']} | {r['n_accepted_now']} |",
              f"| epsilon | {h['epsilon']} (an output) | "
              f"{r['reported_epsilon']} (a criterion) |",
              f"| posterior median distance | {h['posterior_median_distance']} | "
              f"**{r['posterior_median_distance']}** |",
              f"| reference (committed vector) | {r['committed_distance']} | "
              f"{r['committed_distance']} |", "",
              f"The median went from **losing** to the reference by "
              f"{100*(h['posterior_median_distance']/r['committed_distance']-1):.0f}% to "
              f"**beating** it by "
              f"{100*(1-r['posterior_median_distance']/r['committed_distance']):.0f}%.", ""]
    else:
        L += ["## Status: PRESENT", "",
              f"The posterior median ({r['posterior_median_distance']}) fits worse",
              f"than the committed reference ({r['committed_distance']}), with an",
              f"acceptance threshold of {r['reported_epsilon']}.", ""]

    L += ["## What the defect was", "",
          f"Acceptance was a fixed FRACTION — the run kept its best 2% however bad",
          "they were — so the reported epsilon was an **output**, whatever the last",
          "accepted draw happened to score, and never a criterion anything had to",
          "meet. A rejection ABC built that way has no floor: hand it uniformly",
          "terrible draws and it returns 2% of them and calls the result a",
          "posterior.", "",
          "## It was not a truncated prior", "",
          "The reference vector sits exactly on two prior bounds — `k_erastin` at",
          "its low 3.0, `hill` at its high 6.0 — which looks like the box clipping",
          "the optimum. Stepping outside says otherwise.", "",
          "| parameter | value | joint distance |", "|---|--:|--:|"]
    for v, d in _outward(ke, "k_erastin"):
        L.append(f"| `k_erastin` | {v} | {d}{'  (prior low bound)' if v == '3.0' else ''} |")
    for v, d in _outward(hl, "hill"):
        L.append(f"| `hill` | {v} | {d}{'  (prior high bound)' if v == '6.0' else ''} |")
    L += ["",
          ("Pushing `k_erastin` below its bound makes the fit monotonically worse."
           if ke_worse else "`k_erastin` improves outside its bound — the prior IS truncating."),
          ("And `hill` is **inert**: the distance does not move at all between 6 and "
           "10. It is not weakly identified, it has no effect — which is why it is "
           "the one parameter the information-content analysis still finds "
           "indistinguishable from its prior."
           if hill_inert else "`hill` does change the fit outside its bound."), "",
          "## How rare the good region is", "",
          f"Over {s['draws']:,} uniform draws:", "",
          f"* best: **{s['best']}**;",
          f"* draws beating the reference {r['committed_distance']}: "
          f"**{s['n_beating_committed']} of {s['draws']:,}** "
          f"(about {s['n_beating_committed']/max(s['draws'],1):.1e} per draw);", "",
          "That rate is why the old rule failed and why the fix needed more draws",
          "rather than fewer: the region exists and is reachable, but a 1,500-draw",
          "run essentially never lands in it, so a quantile cut over those draws",
          "lands in a shell well above what is achievable.", "",
          "## What the fix does not claim", "",
          "* The tolerance is anchored to a hand-tuned reference vector. That sets",
          "  the bar; it never enters the posterior. The alternative is a bar set by",
          "  whatever the sampler happened to draw, which is what produced the",
          "  defect.",
          "* A better-fitting posterior is not a validated one. It fits the two",
          "  in-vitro dose-response panels better; it says nothing about the in-vivo",
          "  regime, where substituting these values remains inadmissible",
          "  (`analysis/headline-at-fitted-cascade.md`).", ""]
    return "\n".join(L) + "\n"


if __name__ == "__main__":
    raise SystemExit(main())
