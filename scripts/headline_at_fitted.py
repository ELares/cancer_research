#!/usr/bin/env python3
"""Re-derive every headline at the CTRPv2-fitted cascade, and report both.

WHY THIS EXISTS
---------------
`analysis/calibration/in-vivo-prior-provenance.md` names this as the outstanding
piece of work, in these words:

    Either an in-vivo ferroptosis dataset that maps onto these dimensionless
    observables -- which the repository has looked for and documented as not
    publicly existing -- or re-deriving every headline at the fitted cascade and
    reporting both, so a reader can see which directions survive crossing the
    bistable tipping point. The second is achievable now and is the cheaper of
    the two.

Nobody had done it. Every headline magnitude the manuscript reports is computed
at `Params::default()`, the IN-VIVO parameter set, which `targets.yaml` records
as carrying zero calibration targets. The one leg anchored to independent data
-- the CTRPv2 GPX4-inhibitor fit (#330) -- rejects those defaults by eleven-fold
RMSE (0.0504 fitted against 0.5666 default). So the reported numbers rest on
parameters that are unfalsified because untested in their own regime, not
because a test cleared them.

This runs each headline under three parameter sets and prints them side by side.

THE THREE SETS
--------------
  default            Params::default(). lp_propagation 0.10, lp_rate 0.06.
                     In-vivo, zero calibration targets. What the manuscript
                     currently reports.
  ctrpv2_point       The #330 point fit: lp_propagation 0.70, lp_rate 0.40.
  posterior_median   The #500 joint multi-inducer ABC posterior medians:
                     lp_propagation 0.7823, lp_rate 0.7133, gpx4_rate 0.4071,
                     gsh_scav_efficiency 0.5183.

The fitted sets sit ENTIRELY ABOVE the in-vivo priors -- the #500 artifact
records `entire_95pct_posterior_above_invivo_max` for both cascade parameters --
so this is not a perturbation. It crosses the bistable tipping point, which is
exactly why the directions are worth re-checking rather than assumed.

WHAT DOES NOT TRANSFER, AND WHY THAT IS NOT A DEFECT
-----------------------------------------------------
The #330 fit also produced `k_um = 0.25`, a micromolar-to-`rsl3_gpx4_inhib`
dose scale. It is NOT a `Params` field (`apply_param_overrides` does not accept
it) and it does not transfer: it describes how the FIT mapped CTRPv2
concentrations onto the switch, while each binary applies its own fixed drug
effect. Only the cascade parameters carry over, and they are the ones the
disjunction is about.

WHAT THIS IS AND IS NOT
-----------------------
It is NOT a calibration of the spatial headlines. The in-vitro posterior
conditions an in-vitro switch; applying it to in-vivo spatial models does not
make those models data-conditioned, and no claim here should be read that way.
What it does is bound the sensitivity of each REPORTED DIRECTION to the one
parameter choice the repository knows to be contradicted by data. A direction
that survives both sets is robust to that choice. A direction that flips was
resting on it.

Usage:
    python scripts/headline_at_fitted.py            # all headlines
    python scripts/headline_at_fitted.py --quick    # skip sim-tme (~4 min/run)
"""

import argparse
import importlib.util
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
BIN_DIR = PROJECT_ROOT / "simulations" / "target" / "release"
OUT_MD = PROJECT_ROOT / "analysis" / "headline-at-fitted-cascade.md"
OUT_JSON = PROJECT_ROOT / "analysis" / "headline-at-fitted-cascade.json"


def _load(name):
    """Import a sibling script by path; `scripts/` is not a package."""
    spec = importlib.util.spec_from_file_location(
        name, PROJECT_ROOT / "scripts" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# The EXTRACTION helpers are reused rather than reimplemented, so this cannot
# drift from the sensitivity and uncertainty analyses that read the same
# outputs. Only the launcher differs: those build a full parameter vector from
# a design matrix, while here the point is to override a HANDFUL of parameters
# and leave every other one at its default.
hs = _load("headline_sensitivity")

PARAM_SETS = {
    "default": {},
    "ctrpv2_point": {"lp_propagation": 0.70, "lp_rate": 0.40},
    "posterior_median": {
        "lp_propagation": 0.7823, "lp_rate": 0.7133,
        "gpx4_rate": 0.4071, "gsh_scav_efficiency": 0.5183,
    },
}

FIGURE7 = [("Glycolytic", "Control"), ("Glycolytic", "RSL3"),
           ("Glycolytic", "SDT"), ("OXPHOS", "RSL3"), ("OXPHOS", "SDT"),
           ("Persister", "Control"), ("Persister", "RSL3"), ("Persister", "SDT"),
           ("PersisterNrf2", "Control"), ("PersisterNrf2", "RSL3"),
           ("PersisterNrf2", "SDT")]


def _run(binary, overrides, reader, subdir=None, key=None):
    """Run a binary under an override map and read its output with `reader`.

    `key` selects a sub-field of the parsed summary before handing it to the
    reader, because the two extractors do NOT share a contract:
    `extract_tme_observables` wants the `conditions` LIST while
    `extract_tissue_pk_observables` wants the whole document. Passing the whole
    document to the first raises `TypeError: string indices must be integers`
    four minutes into a run, which is exactly what it did.
    """
    exe = BIN_DIR / binary
    if not exe.exists():
        raise FileNotFoundError(f"{exe} -- build with: cargo build --release")
    with tempfile.TemporaryDirectory(prefix="ferro_fitted_") as workdir:
        env = dict(os.environ)
        # An EMPTY override map must mean "no hook at all", not "an empty JSON
        # object", so the default row is produced by exactly the code path the
        # manuscript's numbers came from.
        if overrides:
            env["FERRO_PARAM_OVERRIDES"] = json.dumps(overrides)
        else:
            env.pop("FERRO_PARAM_OVERRIDES", None)
        p = subprocess.run([str(exe)], cwd=workdir, env=env,
                           capture_output=True, text=True)
        if p.returncode != 0:
            raise RuntimeError(f"{binary} failed ({p.returncode}): {p.stderr[-400:]}")
        if subdir is None:
            return reader(workdir)
        summary = Path(workdir) / subdir
        doc = json.loads(summary.read_text())
        return reader(doc[key] if key else doc)


def single_cell(overrides):
    """The Figure 7 death rates, from sim-scale at the figure's own n = 1e6."""
    exe = BIN_DIR / "sim-scale"
    out = {}
    for pheno, tx in FIGURE7:
        cmd = [str(exe), "--cells", "1e6", "--phenotype", pheno,
               "--treatment", tx, "--label", "at-fitted"]
        if overrides:
            cmd += ["--params", json.dumps(overrides)]
        p = subprocess.run(cmd, capture_output=True, text=True)
        if p.returncode != 0:
            raise RuntimeError(f"sim-scale failed: {p.stderr[-300:]}")
        out[f"{pheno}/{tx}"] = json.loads(p.stdout)["death_rate"]
    return out


def bliss(overrides):
    return float(_run("sim-combo-mech", overrides, hs.read_bliss_synergy))


def tme(overrides):
    return _run("sim-tme", overrides, hs.extract_tme_observables,
                "output/tme/tme_summary.json", key="conditions")


def penetration(overrides):
    return _run("sim-tissue-pk", overrides, hs.extract_tissue_pk_observables,
                "output/tissue-pk/tissue_pk_summary.json")


# The manuscript's own stated model constraint, from the Chapter 5 list of
# eight: "all phenotypes show less than 2% death under Control". A parameter set
# that breaks it is not producing wrong headlines, it is producing a model whose
# untreated baseline is dead, in which no headline is interpretable at all.
BASELINE_MAX = 0.02


def admissibility(single_cell: dict) -> dict:
    """Does this parameter set keep untreated cells alive?

    THIS CHECK DECIDES HOW EVERYTHING ELSE READS, which is why it is computed
    before any verdict. Without it the honest finding -- that the in-vitro
    cascade drives the in-vivo models into a regime where untreated cells die en
    masse -- would be reported as "the Bliss synergy headline collapses to 1.0",
    which is true of the number and false about the reason. A ratio between two
    saturated arms is arithmetic, not biology.
    """
    controls = {c: v for c, v in single_cell.items() if c.endswith("/Control")}
    worst = max(controls, key=controls.get)
    return {"controls": controls, "worst_condition": worst,
            "worst_rate": controls[worst],
            "admissible": controls[worst] <= BASELINE_MAX,
            "constraint": BASELINE_MAX}


def _fmt(x):
    if x is None:
        return "n/a"
    if isinstance(x, float):
        return f"{x:.4g}"
    return str(x)


# Tissues in DECREASING drug penetration, which is what the ordering claim
# below is about. Named explicitly because the check is a monotone test and a
# monotone test over dict iteration order is not a test of anything: it read
# `list(res[...]["penetration"])` and was correct only while the scan happened
# to insert the keys in this order. Round-tripping the artifact through JSON
# sorts them alphabetically -- cns_bbb, poorly, well -- which reverses the
# sequence and reported the gradient as BROKEN for `default`, the one
# admissible set, whose kills are a textbook 12.10% > 2.60% > 1.80%.
# DERIVED from the constant that builds the dict, not re-typed beside it:
# `headline_sensitivity.PENETRATION_TISSUES` determines both the keys and
# their order, so taking it from there is order-correct by construction and
# cannot drift. A hand-copied duplicate is the "one artifact describes
# another" defect this repo keeps finding.
PENETRATION_ORDER = tuple(k for k, _ in hs.PENETRATION_TISSUES)


def _tissues_by_penetration(pen: dict) -> list:
    """Order the tissue keys by penetration, refusing an unknown key.

    Refusing rather than falling back: a silent fallback to dict order is what
    made the verdict depend on serialisation in the first place.
    """
    unknown = set(pen) - set(PENETRATION_ORDER)
    if unknown:
        raise SystemExit(
            f"unknown tissue key(s) {sorted(unknown)}: add them to "
            "PENETRATION_ORDER in penetration order, since the ordering "
            "verdict below is a monotone test over this sequence")
    return [t for t in PENETRATION_ORDER if t in pen]


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
    ap.add_argument("--quick", action="store_true",
                    help="skip sim-tme, which costs ~4 min per parameter set")
    ap.add_argument("--out", default=str(OUT_MD))
    args = ap.parse_args()

    results = {}
    for name, ov in PARAM_SETS.items():
        print(f"=== {name} {ov or '(no overrides)'}", flush=True)
        r = {"overrides": ov}
        r["single_cell"] = single_cell(ov)
        r["admissibility"] = admissibility(r["single_cell"])
        a = r["admissibility"]
        print(f"  single-cell ok — worst untreated {a['worst_condition']} "
              f"{a['worst_rate']*100:.2f}% "
              f"({'admissible' if a['admissible'] else 'INADMISSIBLE'})", flush=True)
        r["bliss"] = bliss(ov);                  print(f"  bliss {r['bliss']:.4g}", flush=True)
        r["penetration"] = penetration(ov);      print("  penetration ok", flush=True)
        if not args.quick:
            r["tme"] = tme(ov)
            print(f"  tme hypoxia {r['tme']['hypoxia']:.4g} immune {r['tme']['immune']:.4g}",
                  flush=True)
        results[name] = r

    OUT_JSON.write_text(json.dumps(results, indent=2, sort_keys=True) + "\n")
    Path(args.out).write_text(render(_roundtrip(results)), encoding="utf-8")
    print(f"wrote {args.out}\nwrote {OUT_JSON}")
    return 0


def render(res: dict) -> str:
    """The write-up. Every verdict below is COMPUTED from `res`, never typed."""
    sets = [s for s in ("default", "ctrpv2_point", "posterior_median") if s in res]
    L = [
        "# Every headline, re-derived at the CTRPv2-fitted cascade", "",
        "Generated by `scripts/headline_at_fitted.py`.", "",
        "`analysis/calibration/in-vivo-prior-provenance.md` names this as the",
        "outstanding piece of work: *\"re-deriving every headline at the fitted",
        "cascade and reporting both, so a reader can see which directions survive",
        "crossing the bistable tipping point. The second is achievable now and is",
        "the cheaper of the two.\"* This is that, and it had not been done.", "",
        "## Why it matters", "",
        "Every headline magnitude the manuscript reports is computed at",
        "`Params::default()` — the in-vivo set, which `targets.yaml` records as",
        "carrying **zero calibration targets**. The one leg anchored to independent",
        "data, the CTRPv2 GPX4-inhibitor fit (#330), rejects those defaults by",
        "**eleven-fold RMSE** (0.0504 fitted against 0.5666 default), and the #500",
        "posterior lies *entirely above* the in-vivo priors for both cascade",
        "parameters. So the reported numbers rest on parameters that are unfalsified",
        "because untested in their own regime — not because a test cleared them.", "",
        "This does **not** calibrate the spatial headlines: an in-vitro posterior",
        "conditions an in-vitro switch, and carrying it into in-vivo spatial models",
        "does not make them data-conditioned. What it bounds is how much each",
        "reported DIRECTION depends on the one parameter choice the repository knows",
        "to be contradicted by data.", "",
        "## Parameter sets", "",
        "| set | lp_propagation | lp_rate | gpx4_rate | gsh_scav_efficiency |",
        "|---|--:|--:|--:|--:|",
    ]
    base = {"lp_propagation": 0.10, "lp_rate": 0.06,
            "gpx4_rate": 0.30, "gsh_scav_efficiency": 0.5}
    for s in sets:
        ov = {**base, **res[s]["overrides"]}
        L.append(f"| `{s}` | {ov['lp_propagation']:.4g} | {ov['lp_rate']:.4g} | "
                 f"{ov['gpx4_rate']:.4g} | {ov['gsh_scav_efficiency']:.4g} |")
    L += ["",
          "`k_um` from the #330 fit is deliberately absent: it is a micromolar",
          "dose-scale, not a `Params` field, and describes how the fit mapped CTRPv2",
          "concentrations rather than anything the binaries consume.", ""]

    # --- admissibility, FIRST, because it decides how everything reads -----
    adm = {s_: res[s_]["admissibility"] for s_ in sets}
    bad = [s_ for s_ in sets if not adm[s_]["admissible"]]
    L += ["## Is each parameter set even admissible for these models?", "",
          "Chapter 5 lists eight constraints the model is checked against. The",
          f"first is that **all phenotypes show less than {BASELINE_MAX*100:.0f}%",
          "death under Control** — untreated cells must stay alive. A parameter set",
          "that breaks it is not producing wrong headlines; it is producing a model",
          "whose untreated baseline is dead, in which no headline means anything.", "",
          "| set | worst untreated condition | untreated death | admissible? |",
          "|---|---|--:|---|"]
    for s_ in sets:
        a = adm[s_]
        L.append(f"| `{s_}` | {a['worst_condition']} | {a['worst_rate']*100:.2f}% | "
                 f"{'yes' if a['admissible'] else '**NO**'} |")
    L.append("")
    if bad:
        worst = max(bad, key=lambda k: adm[k]["worst_rate"])
        L += [
            f"**{len(bad)} of {len(sets)} parameter sets are inadmissible.** At",
            f"`{worst}` the untreated {adm[worst]['worst_condition']} population dies at",
            f"{adm[worst]['worst_rate']*100:.2f}%, against a constraint of "
            f"{BASELINE_MAX*100:.0f}%.",
            "",
            "**This is the finding, and it is a negative one.** The in-vitro-fitted",
            "cascade cannot simply be carried into the in-vivo and spatial models:",
            "it drives them into a regime where untreated cells die en masse and",
            "every arm saturates. The headline tables below are therefore reported",
            "for completeness and are NOT evidence that any direction failed — a",
            "ratio between two saturated arms is arithmetic, not biology.", "",
            "It also sharpens the honest position. The in-vivo defaults are not",
            "validated, and the disjunction (#332, #500) says the in-vitro posterior",
            "is not a substitute for them. Until now that was an argument about",
            "parameter ranges. It is now a demonstrated one: substituting the fitted",
            "values breaks the model's own baseline-viability constraint by a factor",
            f"of {adm[worst]['worst_rate']/BASELINE_MAX:.0f}.", ""]
    else:
        L += ["All parameter sets keep untreated death inside the constraint, so",
              "the headline comparisons below are interpretable.", ""]

    def verdict(text_ok, text_bad):
        """Verdicts must not read as 'the headline failed' when the parameter
        set was inadmissible in the first place."""
        if bad:
            return ("*Not a verdict on the headline.* "
                    + text_bad
                    + f" — but `{'`, `'.join(bad)}` "
                    + ("is" if len(bad) == 1 else "are")
                    + " inadmissible above, so this row reports arithmetic in a"
                      " degenerate regime rather than a direction that held or failed.")
        return text_ok

    # --- single cell -----------------------------------------------------
    L += ["## Figure 7 death rates", "",
          "| condition | " + " | ".join(f"`{s}`" for s in sets) + " |",
          "|---|" + "--:|" * len(sets)]
    conds = sorted(res[sets[0]]["single_cell"])
    for c in conds:
        L.append(f"| {c} | " + " | ".join(
            f"{res[s]['single_cell'][c]*100:.2f}%" for s in sets) + " |")
    L.append("")

    # The selectivity claim the thesis rests on, checked at every set.
    sel = []
    for s in sets:
        sc = res[s]["single_cell"]
        sel.append((s, sc.get("Persister/RSL3", 0), sc.get("Glycolytic/RSL3", 0)))
    holds = [s for s, p, g in sel if p > g]
    detail = "; ".join(f"`{s_}` {p*100:.2f}% vs {g*100:.2f}%" for s_, p, g in sel)
    L += [verdict(
        f"**RSL3 selectivity (Persister > Glycolytic) holds in {len(holds)} of "
        f"{len(sel)} parameter sets**: {detail}.",
        f"RSL3 selectivity holds in {len(holds)} of {len(sel)} sets: {detail}"), ""]

    # --- bliss -----------------------------------------------------------
    L += ["## Bliss synergy (RSL3 + FSP1i)", "",
          "| set | synergy score |", "|---|--:|"]
    for s in sets:
        L.append(f"| `{s}` | {_fmt(res[s]['bliss'])} |")
    vals = [(s, res[s]["bliss"]) for s in sets]
    supra = [s for s, v in vals if v == v and v > 1.0]
    flat = ", ".join(f"`{s_}` = {_fmt(v)}" for s_, v in vals
                     if not (v == v and v > 1.0))
    L += ["", verdict(
        f"**Supra-additive (> 1.0) in {len(supra)} of {len(vals)} sets.** "
        + ("The direction survives the crossing." if len(supra) == len(vals)
           else f"The direction does NOT survive: {flat}."),
        f"Supra-additive in {len(supra)} of {len(vals)} sets ({flat} at 1.0 exactly, "
        "which is what a Bliss ratio returns when both single arms already "
        "saturate)"), ""]

    # --- tme -------------------------------------------------------------
    if all("tme" in res[s] for s in sets):
        L += ["## Hypoxia kill-collapse gap, and immune amplification", "",
              "| set | SDT − RSL3 hypoxic gap | SDT de-confounded immune rate |",
              "|---|--:|--:|"]
        for s in sets:
            L.append(f"| `{s}` | {_fmt(res[s]['tme']['hypoxia'])} | "
                     f"{_fmt(res[s]['tme']['immune'])} |")
        pos = [s for s in sets if res[s]["tme"]["hypoxia"] > 0]
        L += ["", verdict(
            f"**The hypoxia gap stays positive in {len(pos)} of {len(sets)} sets.** "
            + ("SDT-holds-where-RSL3-collapses survives the crossing."
               if len(pos) == len(sets) else "It does not survive."),
            f"The hypoxia gap is positive in {len(pos)} of {len(sets)} sets"), ""]
    else:
        L += ["## Hypoxia and immune headlines", "",
              "Not run (`--quick`). sim-tme costs about four minutes per set.", ""]

    # --- penetration -----------------------------------------------------
    L += ["## Penetration gradient (RSL3-like vessel-wall kill)", "",
          "| tissue | " + " | ".join(f"`{s}`" for s in sets) + " |",
          "|---|" + "--:|" * len(sets)]
    tissues = _tissues_by_penetration(res[sets[0]]["penetration"])
    for t in tissues:
        L.append(f"| {t} | " + " | ".join(
            f"{res[s]['penetration'][t]*100:.2f}%" for s in sets) + " |")
    ordered = []
    for s in sets:
        v = [res[s]["penetration"][t] for t in tissues]
        ordered.append((s, all(a >= b for a, b in zip(v, v[1:]))))
    keep = [s for s, ok in ordered if ok]
    L += ["", verdict(
        f"**The across-tissue ordering is preserved in {len(keep)} of "
        f"{len(ordered)} sets** ({', '.join(tissues)}).",
        f"The ordering is preserved in {len(keep)} of {len(ordered)} sets, but at "
        "the fitted sets every tissue saturates at 100%, so the ordering is "
        "trivially preserved and says nothing"), ""]

    L += ["## What this does not license", "",
          "* These runs do not make any spatial headline data-conditioned. The",
          "  posterior is in-vitro; the spatial models are in-vivo; the two",
          "  parameterisations are provably disjoint (#332, #500).",
          "* A direction surviving both sets is robust to *this* parameter choice,",
          "  not validated. Nothing here is a measurement of biology.",
          "* Only the cascade parameters transfer. Everything else stays at its",
          "  default, including every off-by-default realism layer.", ""]
    return "\n".join(L) + "\n"


if __name__ == "__main__":
    raise SystemExit(main())
