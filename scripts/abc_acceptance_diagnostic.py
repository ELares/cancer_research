#!/usr/bin/env python3
"""Check the joint ABC's reference, acceptance rule and local boundary probes.

The original fixed-fraction run accepted its best 2% irrespective of fit. Its
coordinate-wise posterior median scored 0.2413 against a known reference at
0.2202 on the same original targets. That history is retained separately from
the current-target comparison: the dose-support correction removes the
unsupported 100 µM erastin point and changes the retained cohorts, so comparing
old and new distances would confound the acceptance rule with a changed target.

Current probes score the reference and median vector against the same stored
targets, vary two coordinates outside their prior bounds, and estimate how often
fresh prior draws improve on the reference. These are fit diagnostics, not proof
of a global optimum, mechanistic identifiability or biological validity.

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
    "reference_distance": 0.2202,
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
    if args.draws < 1:
        ap.error("--draws must be positive")

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
        "acceptance_rule": art.get("acceptance_rule"),
        "target_source": art.get("target_source"),
        "target_dose_grids_um": {"ML162": rd, "ERASTIN": ed},
        "committed_vector": COMMITTED,
        "committed_distance": round(float(committed), 4),
        "posterior_median_distance": round(float(median_d), 4),
        "reported_epsilon": eps,
        "epsilon_excess_over_committed": round(float(eps / committed - 1.0), 4) if eps else None,
        "prior_truncation_test": {k: {str(a): round(float(b), 4) for a, b in v.items()}
                                  for k, v in outside.items()},
        "n_accepted_now": art.get("n_accepted"),
        "min_posterior": art.get("min_posterior"),
        "underpowered": art.get("underpowered"),
        "sampling": {
            "draws": args.draws,
            "draws_of_record": art.get("n_draws"),
            "seed": 7,
            "best": round(float(ds.min()), 4),
            "quantile_2pct": round(float(np.quantile(ds, 0.02)), 4),
            "n_beating_committed": int((ds < committed).sum()),
            "n_inside_epsilon": int((ds <= eps).sum()) if eps is not None else None,
            "frac_inside_epsilon": round(float((ds <= eps).mean()), 4) if eps is not None else None,
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
    hill_unchanged = len(set(hl.values())) == 1
    ke_no_improvement = all(v >= ke["3.0"] for v in ke.values())
    median_better = r["posterior_median_distance"] < r["committed_distance"]
    minimum = r.get("min_posterior")
    count = r.get("n_accepted_now")
    # A recorded shortfall or a count below the stated minimum is sufficient;
    # absent metadata must not silently turn into a claim of adequate sampling.
    underpowered = (r.get("underpowered") is True
                    or (minimum is not None and count is not None and count < minimum))
    adequacy_recorded = (minimum is not None and count is not None
                         and r.get("underpowered") is not None)
    h = HISTORICAL

    L = ["# The ABC acceptance rule: the defect, and its fix", "",
         "Generated by `scripts/abc_acceptance_diagnostic.py`.", ""]

    if r.get("acceptance_rule") == "tolerance":
        status = ("ACCEPTANCE RULE FIXED; POSTERIOR UNDERPOWERED" if underpowered
                  else "ACCEPTANCE RULE FIXED" if adequacy_recorded
                  else "ACCEPTANCE RULE FIXED; SAMPLING ADEQUACY UNKNOWN")
        L += [f"## Status: {status}", "",
              "The joint run uses a reference-based tolerance instead of filling a",
              "fixed acceptance quota. This resolves the acceptance-rule defect;",
              "it does not establish that the accepted sample supports reliable",
              "posterior summaries.", ""]
    else:
        L += ["## Status: TOLERANCE RULE NOT RECORDED", "",
              "This diagnostic does not record a tolerance-based acceptance rule.",
              "An older artifact may lack that metadata; regenerate before using",
              "this document to establish the current rule.", ""]

    L += ["## Historical run: original targets", "",
          f"The {h['rule']} run drew {h['n_draws']:,} vectors and accepted",
          f"{h['n_accepted']}, with epsilon {h['epsilon']} determined by that quota.",
          f"Its coordinate-wise posterior median scored {h['posterior_median_distance']}",
          f"against a reference distance of {h['reference_distance']} on the SAME",
          "original targets. Those targets used unfiltered cohorts and included",
          "the unsupported 100 µM erastin point.", "",
          "Historical and corrected-target distances are not directly comparable.",
          "Removing an extrapolated target and changing the retained cohort changes",
          "the objective; a lower distance does not by itself demonstrate better",
          "biological validity. Historical reference values are preserved above,",
          "not recomputed using the current targets.", "",
          "## Current-target check", ""]
    grids = r.get("target_dose_grids_um")
    if grids:
        for name, doses in grids.items():
            L.append(f"* {name} dose grid (µM): {doses}.")
    else:
        L.append("Dose-grid metadata is absent in this older diagnostic artifact.")
    if r.get("target_source"):
        L.append(f"* Input CSV SHA256: `{r['target_source']['sha256']}`.")
    L += ["",
          f"The run drew {s['draws_of_record']:,} vectors and accepted",
          f"{r['n_accepted_now']}; its recorded epsilon is {r['reported_epsilon']}.", ""]
    if underpowered:
        L += [f"**Posterior underpowered:** {count} accepted draws; recorded minimum",
              f"{minimum if minimum is not None else 'not available'}.",
              "The coordinate medians, interval endpoints and parameter-draw bands",
              "are exploratory summaries of this small accepted sample, not reliable",
              "posterior estimates. A favorable median-vector distance cannot resolve",
              "this sampling shortfall.", ""]
    elif adequacy_recorded:
        L += [f"The accepted count meets the recorded minimum of {minimum}.",
              "Meeting that minimum alone does not establish Monte Carlo stability",
              "or posterior calibration.", ""]
    else:
        L += ["Sampling-adequacy metadata is incomplete. This diagnostic cannot",
              "establish whether the recorded minimum accepted count was met.", ""]
    L += [
          f"On these targets the reference scores **{r['committed_distance']}**",
          f"and the coordinate-wise accepted-draw median scores **{r['posterior_median_distance']}**.",
          ("That single constructed vector fits better than the reference on the current targets."
           if median_better else
           "That single constructed vector does not fit better than the reference on the current targets."),
          "A coordinate-wise median need not be an accepted vector, and its distance",
          "is not the median distance of accepted draws. This check alone neither",
          "establishes nor refutes compliance with the acceptance tolerance, and",
          "does not establish that the posterior is adequately sampled.", ""]

    L += ["## What the defect was", "",
          f"Acceptance was a fixed FRACTION — the run kept its best 2% however bad",
          "they were — so the reported epsilon was an **output**, whatever the last",
          "accepted draw happened to score, and never a criterion anything had to",
          "meet. A rejection ABC built that way has no floor: hand it uniformly",
          "terrible draws and it returns 2% of them and calls the result a",
          "posterior.", "",
          "## Boundary probes on the current targets", "",
          "The reference vector sits exactly on two prior bounds — `k_erastin` at",
          "its low 3.0, `hill` at its high 6.0. These probes change one coordinate",
          "at that reference; they do not search for a global optimum outside the prior.", "",
          "| parameter | value | joint distance |", "|---|--:|--:|"]
    for v, d in _outward(ke, "k_erastin"):
        L.append(f"| `k_erastin` | {v} | {d}{'  (prior low bound)' if v == '3.0' else ''} |")
    for v, d in _outward(hl, "hill"):
        L.append(f"| `hill` | {v} | {d}{'  (prior high bound)' if v == '6.0' else ''} |")
    L += ["",
          ("The tested lower `k_erastin` values do not improve the reported distance."
           if ke_no_improvement else
           "A tested `k_erastin` value outside the prior improves the reported distance."),
          ("The reported `hill` distances are unchanged at the probed values 6, 8 and 10. "
           "This does not establish that the parameter is inert throughout its prior."
           if hill_unchanged else "The reported `hill` distances differ across the tested values."), "",
          "## Fresh prior-search diagnostic", "",
          f"Over {s['draws']:,} uniform draws:", "",
          f"* best: **{s['best']}**;",
          f"* draws beating the reference {r['committed_distance']}: "
          f"**{s['n_beating_committed']} of {s['draws']:,}** "
          f"(about {s['n_beating_committed']/max(s['draws'],1):.1e} per draw);"]
    if s.get("n_inside_epsilon") is not None:
        L += [f"* draws at or below the recorded tolerance: **{s['n_inside_epsilon']} "
              f"of {s['draws']:,}**."]
    L += ["",
          "This estimates how frequently uniform prior draws beat the reference on",
          "the current targets. It does not reconstruct the sampling rate on the",
          "different historical target set. Zero observed hits would not prove",
          "that no better-fitting region exists. These additional draws diagnose",
          "search difficulty; they are not added to the joint posterior and do not",
          "cure a shortfall in its accepted sample.", "",
          "## What the fix does not claim", "",
          "* The tolerance is anchored to a hand-tuned reference vector. That sets",
          "  the bar; it never enters the posterior. The alternative is a bar set by",
          "  whatever the sampler happened to draw, which is what produced the",
          "  defect.",
          "* A better-fitting posterior is not a validated one. This is a comparison",
          "  against two in-vitro fitted-curve summaries. The separate",
          "  `analysis/headline-at-fitted-cascade.md` substitution test found recorded",
          "  historical vectors inadmissible in the spatial model. It did not test",
          "  the corrected posterior and cannot establish its admissibility.", ""]
    return "\n".join(L) + "\n"


if __name__ == "__main__":
    raise SystemExit(main())
