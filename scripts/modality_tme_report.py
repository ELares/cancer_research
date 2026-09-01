"""How each treatment arm fares under the tumour microenvironment.

WHY THIS IS THE DEEP QUESTION
-----------------------------
`analysis/modality-panel.md` answers "what does each arm do to a naive
tumour". That is the shallow question and the one a reader over-interprets.
This project's actual contribution has never been kill fractions: it is the
three RESISTANCE AXES the ferroptosis chapters establish -- hypoxia, stromal
shielding and acidic pH -- and the finding that pharmacologic and physical
modalities respond to them differently.

Every arm the coverage campaign added was untested against those axes, which
is precisely the gap between "the engine can express it" and "the engine has
something to say about it". This page closes it, and reports the two axes
where it does not.

THE AXES ARE APPLIED THROUGH THE SAME HELPERS
---------------------------------------------
Hypoxia scales exogenous-ROS yield and radiation's delivered dose through the
identical Alper-Howard-Flanders hyperbola; stroma raises the antioxidant
setpoint; acidic pH traps weak bases outside the cell. No arm gets a
mechanism-specific fudge, because a modality comparison only means anything
if the environment is identical.

Offline: reads the binary's committed sweep JSON.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
SWEEP = REPO / "simulations" / "output" / "modality-panel" / "modality_tme_sweep.json"
OUT_MD = REPO / "analysis" / "modality-tme.md"
OUT_JSON = REPO / "analysis" / "modality-tme.json"

# An axis whose worst-case relative effect is below this on EVERY arm is
# reported INERT rather than quietly listed with tiny numbers. A table of
# 0.08% differences invites a reader to believe an axis was tested when the
# configuration could not see it.
INERT_THRESHOLD = 0.02


def scan() -> dict:
    if not SWEEP.exists():
        raise SystemExit(
            f"{SWEEP} not found -- run "
            "`cargo run --release -p sim-modality-panel -- --tme-sweep`")
    return json.loads(SWEEP.read_text())


AXIS_KEYS = ("hypoxic", "stroma", "acidic", "deep", "heterogeneous")
STRATA = ("phenotype",)


def _effect(conds, arm, axis, phenotype=None):
    """Largest relative change in `arm` attributable to `axis` alone.

    Paired: every ON condition is compared with the OFF condition identical in
    every OTHER axis and in the phenotype, so an effect cannot be manufactured
    by something that happens to co-vary. The phenotype is a stratum rather
    than an axis -- comparing a persister against a glycolytic cell is not a
    treatment effect.
    """
    worst = 0.0
    pool = [c for c in conds
            if phenotype is None or c.get("phenotype") == phenotype]
    for on in [c for c in pool if c[axis]]:
        off = next((c for c in pool
                    if not c[axis]
                    and all(c[k] == on[k] for k in AXIS_KEYS if k != axis)
                    and all(c.get(k) == on.get(k) for k in STRATA)), None)
        if off is None:
            continue
        base = off["arms"].get(arm)
        got = on["arms"].get(arm)
        if base is None or got is None or base <= 0:
            continue
        # SIGNED. An earlier version took the absolute value and called every
        # result a "loss", which produced the impossible report that an arm
        # lost 121% of its kill. It had GAINED: heterogeneity doubles the
        # pharmacologic arm's kill from an acid-suppressed baseline, because
        # variance supplies a low-glutathione tail that dies even when the
        # mean cell resists. An axis that can help is not the same kind of
        # thing as one that only hurts, and collapsing the sign hid a result.
        rel = (got - base) / base
        if abs(rel) > abs(worst):
            worst = rel
    return worst


def assemble(raw: dict) -> dict:
    conds = raw["conditions"]
    arms = sorted({a for c in conds for a in c["arms"]})
    phenos = sorted({c.get("phenotype") for c in conds if c.get("phenotype")})
    axes = [("hypoxic", "hypoxia"), ("stroma", "stromal shielding"),
            ("acidic", "acidic pH"), ("deep", "depth"),
            ("heterogeneous", "clonal heterogeneity")]
    # POOLED across phenotypes for the headline, and PER PHENOTYPE beside it,
    # because which axis bites turns out to depend on the phenotype -- and a
    # single-phenotype sweep reported two axes inert for that reason alone.
    effects = {label: {a: _effect(conds, a, key) for a in arms}
               for key, label in axes}
    # WHICH CELLS ARE UNDEFINED RATHER THAN ZERO. `_effect` skips a pair whose
    # unstressed kill is 0, because the relative change has no denominator --
    # and `worst` then stays 0.0, which is indistinguishable from an axis that
    # genuinely cannot move the arm. fig32's footer asserted the second reading
    # for a whole row that was really the first.
    undefined = {}
    for ph in phenos:
        pool = [c for c in conds if c.get("phenotype") == ph]
        for _, label in axes:
            for a in arms:
                if all(c["arms"].get(a, 0.0) <= 0 for c in pool):
                    undefined.setdefault(ph, {}).setdefault(label, []).append(a)
    per_pheno = {
        ph: {label: {a: _effect(conds, a, key, ph) for a in arms}
             for key, label in axes}
        for ph in phenos
    }
    inert = [label for _, label in axes
             if max((abs(v) for v in effects[label].values()), default=0.0)
             < INERT_THRESHOLD]
    live = [label for _, label in axes if label not in inert]
    order = sorted(
        arms, key=lambda a: max((abs(effects[l][a]) for l in live), default=0.0))
    # Which axis MOVES each phenotype most -- the finding the stratification
    # exists to produce, and it needs three things the first version dropped.
    #
    # (1) THE SIGN. `max(abs(v))` reported the persister winner as a
    #     "pressure" when it is a +124% GAIN -- clonal heterogeneity RAISES the
    #     pharmacologic arm's kill, which the same section states four
    #     paragraphs later. That is the sign collapse this whole analysis
    #     retracts, surviving in the field the headline reads.
    # (2) THE TIE. On the glycolytic side `acidic pH` and `depth` are both
    #     exactly 1.0, and `max` returned whichever the hard-coded `axes` list
    #     mentions first. Swapping two tuples changed the manuscript's stated
    #     finding on identical data, so ties are reported as ties.
    # (3) THE ARM. A magnitude with no arm attached cannot be checked against
    #     the figure that draws it.
    def _dominant(ph):
        scored = []
        for _, l in axes:
            vals = per_pheno[ph][l]
            arm = max(vals, key=lambda a: abs(vals[a])) if vals else None
            scored.append((l, vals.get(arm, 0.0), arm))
        top = max(abs(v) for _, v, _ in scored)
        tied = [t for t in scored if abs(abs(t[1]) - top) < 1e-9]
        label, value, arm = tied[0]
        return {"axis": label, "value": value, "arm": arm,
                "direction": "gain" if value > 0 else "loss",
                "tied_with": sorted(t[0] for t in tied[1:])}
    dominant = {ph: _dominant(ph) for ph in phenos}
    # The amplification axis, reported separately because it is not a
    # resistance pressure -- it is what an arm's death mode EARNS.
    amp = []
    for c in conds:
        q = c.get("damp_per_death") or {}
        if not q:
            continue
        amp.append({
            "phenotype": c.get("phenotype"),
            "stressed": any(c[k] for k in AXIS_KEYS),
            "heterogeneous": c.get("heterogeneous"),
            "damp_per_death": q,
            "kill": {k: c["arms"][k] for k in q if k in c["arms"]},
        })
    base = [a for a in amp
            if not a["stressed"] and not a["heterogeneous"]]
    quality_ratio = None
    for a in base:
        q = a["damp_per_death"]
        if q.get("RSL3", 0) > 0 and q.get("SDT", 0) > 0:
            quality_ratio = {"phenotype": a["phenotype"],
                             "sdt": q["SDT"], "rsl3": q["RSL3"],
                             "ratio": q["SDT"] / q["RSL3"],
                             "kill_ratio": (a["kill"]["SDT"] / a["kill"]["RSL3"]
                                            if a["kill"].get("RSL3", 0) > 0 else None)}
            break
    return dict(raw, amplification=amp, quality_ratio=quality_ratio,
                arms=arms, phenotypes=phenos, effects=effects,
                effects_by_phenotype=per_pheno, dominant_axis=dominant,
                undefined_cells=undefined,
                inert_axes=inert, live_axes=live, robustness_order=order,
                inert_threshold=INERT_THRESHOLD)


def render(d: dict) -> str:
    conds = d["conditions"]
    arms = d["arms"]
    L = ["# How each arm fares under the tumour microenvironment", "",
         "*Generated by `scripts/modality_tme_report.py --render-only` from "
         "`sim-modality-panel --tme-sweep`. Offline; reads the binary's "
         "committed sweep.*", "",
         "The head-to-head panel asks what each arm does to a NAIVE tumour. "
         "This asks the question this project actually exists to ask: what "
         "the tumour microenvironment does to it. Every axis is applied "
         "through the same helpers the ferroptosis chapters use, so no arm "
         "gets a mechanism-specific adjustment.", "",
         "| condition | " + " | ".join(arms) + " |",
         "|---|" + "--:|" * len(arms)]
    for c in conds:
        lbl = ("hypoxic" if c["hypoxic"] else "normoxic")
        if c["stroma"]:
            lbl += " + stroma"
        if c["acidic"]:
            lbl += " + acid"
        L.append(f"| {lbl} | "
                 + " | ".join(f"{c['arms'][a] * 100:.2f}%" for a in arms) + " |")

    live, inert = d["live_axes"], d["inert_axes"]
    L += ["", "## Which axis bites depends on the phenotype, and that is the "
          "finding", ""]
    dom = d.get("dominant_axis", {})
    if len(dom) > 1:
        def _one(ph, r):
            tie = (" — tied exactly with " + " and ".join(r["tied_with"])
                   + ", so which one is named is arbitrary" if r["tied_with"]
                   else "")
            return (f"**{ph}** — {r['axis']}, a {r['direction']} of "
                    f"{abs(r['value']) * 100:.0f}% for `{r['arm']}`{tie}")
        pairs = "; ".join(_one(ph, r) for ph, r in sorted(dom.items()))
        L += [f"The axis that moves each cell state most is not the same one: "
              f"{pairs}.", "",
              "**Read the direction, not just the size.** An earlier version "
              "of this section called both of these a *pressure* and printed "
              "the magnitude alone, because the field behind it was "
              "`max(abs(v))`. One of the two is a GAIN, which makes "
              "\"dominant pressure\" exactly wrong for it — the same sign "
              "collapse this page retracts elsewhere, surviving in the field "
              "the headline reads. Where two axes tie exactly, the winner was "
              "whichever a hard-coded list mentioned first, so swapping two "
              "entries changed the stated finding on identical data; ties are "
              "reported as ties now.", "",
              "**That is why the first version of this sweep reported two "
              "axes INERT.** It ran one phenotype. At the glycolytic state "
              "the delivered arm kills essentially nothing, so ion trapping "
              "had nothing to scale and the antioxidant buffer was swamped by "
              "an order-of-magnitude ROS insult — and both axes looked dead "
              "when what was dead was the configuration's ability to see "
              "them. Running the persister state as well, the state this "
              "project's thesis is actually about, makes both bite.", "",
              "An axis reported inert is a statement about the run, not about "
              "the biology, and the distinction is easy to lose. This page "
              "now stratifies rather than pooling, so an axis cannot be "
              "declared irrelevant because it was measured in the one state "
              "that cannot feel it.", ""]

    if live:
        L += ["## The ordering under the axes that bite", "",
              "**Largest relative MOVE first, and read the verb.** This header "
              "said \"largest relative loss first\" over a list whose top two "
              "entries are GAINS, because the sort is on `abs` -- the same "
              "sign collapse this page retracts elsewhere, surviving in a "
              "heading. An arm that GAINS is not an arm that is most exposed.",
              ""]
        rows = sorted(arms,
                      key=lambda a: -max(abs(d["effects"][l][a]) for l in live))
        for a in rows:
            worst_l = max(live, key=lambda l: abs(d["effects"][l][a]))
            worst = d["effects"][worst_l][a]
            verb = ("loses" if worst < 0
                    else "GAINS" if worst > 0 else "is unmoved by")
            L.append(f"* `{a}` — {verb} {abs(worst) * 100:.0f}% of its kill, "
                     f"largest to {worst_l}")
        L += ["", "That ordering was not tuned for. It follows from the "
              "mechanisms: a delivered drug loses most, because everything "
              "between the vessel and the target can stop it; an arm whose "
              "lethality depends on oxygen loses next; an arm whose dose is "
              "merely modified by oxygen loses less, which is what an "
              "enhancement ratio says a dose-modifying factor should do; and "
              "a threshold arm loses nothing, because a destroyed cell does "
              "not care about any of it.", ""]

    if inert:
        L += [f"**{', '.join(inert)} remain inert even stratified**, moving "
              f"every arm by less than {d['inert_threshold'] * 100:.0f}% in "
              "every cell state tested. That is a stronger statement than the "
              "pooled version could make, and still not a claim that the arms "
              "are resistant to them — only that this configuration cannot "
              "apply that pressure.", ""]

    gains = {a: max((d["effects"][l][a] for l in live), key=abs, default=0.0)
             for a in arms}
    helped = {a: v for a, v in gains.items() if v > 0.10}
    if helped:
        L += ["## One axis can HELP, and it helps the arm that is failing", ""]
        for a, v in sorted(helped.items(), key=lambda kv: -kv[1]):
            ax = max(live, key=lambda l: d["effects"][l][a])
            L += [f"`{a}` GAINS {v * 100:.0f}% of its kill to {ax}.", ""]
        L += ["Clonal heterogeneity is the only axis here that can raise an "
              "arm's kill rather than lower it, and it raises the "
              "pharmacologic arm from an acid-suppressed baseline. The "
              "mechanism is not subtle: widening the antioxidant setpoint "
              "while holding its MEAN fixed supplies a low-glutathione tail "
              "that dies even when the average cell resists. Variance rescues "
              "a marginal drug.", "",
              "**An earlier version of this page could not have reported "
              "that.** It took the absolute value of every change and called "
              "the result a loss, which produced the impossible line that an "
              "arm had lost 121% of its kill. Collapsing the sign did not "
              "just mislabel a number; it hid a result, because an axis that "
              "can help is not the same kind of thing as one that only "
              "hurts.", ""]

    qr = d.get("quality_ratio")
    if qr:
        L += ["## Immune amplification is a COUNT effect, not a quality effect", "",
              "The manuscript reports that sonodynamic therapy generates far "
              "more immune kills than the pharmacologic inducer, and it does. "
              "This measures WHY, and the answer refines the claim rather "
              "than confirming it.", "",
              f"In the unstressed {qr['phenotype']} state, DAMP release PER "
              f"DEATH is {qr['sdt']:.2f} for SDT against {qr['rsl3']:.2f} for "
              f"RSL3 — a ratio of **{qr['ratio']:.2f}×**. The two death modes "
              "are not very different in quality. What differs is how many "
              "cells they kill"
              + (f", and there the ratio is {qr['kill_ratio']:.2f}×"
                 if qr.get("kill_ratio") else "")
              + ". The amplification advantage is overwhelmingly a count "
              "effect.", "",
              "**And the first version of this measurement got it wrong in an "
              "instructive way.** It read lipid peroxidation at the moment of "
              "death, which returns approximately the death threshold FOR "
              "EVERY ARM by construction — death IS the threshold crossing. "
              "Measured that way both arms reported ~10.2 and the quality "
              "difference vanished entirely. It had not vanished; it was "
              "being measured before it happens. The spatial binaries read "
              "`lp_at_grace_end` rather than `lp` for exactly this reason, "
              "and the field name says so: an arm delivering exogenous ROS "
              "keeps climbing through the post-death grace period while one "
              "that merely disabled a repair enzyme does not.", "",
              "A quantity that is equal for every arm by construction is not "
              "a measurement of anything, and the tell was that it came out "
              "suspiciously close to a threshold the model defines.", ""]

    L += ["## What this does not say", "",
          (f"**{len(inert)} of the {len(axes)} axes were not tested, they were "
           "not visible.** That is a weaker statement than 'these arms are "
           "robust to them' and it is the true one. An arm can only be shown "
           "resistant to a pressure the configuration can apply."
           if inert else
           "**Every axis in this sweep moves at least one arm.** That was not "
           "always so: an earlier configuration left two of them inert and "
           "this section said so unconditionally, which meant it kept saying "
           "so after the sweep grew and nothing was inert any more. An axis "
           "reported inert is a statement about the run, and so is an axis "
           "reported live."), "",
          "**The immune arm's hypoxia response is a PREDICTION, not a "
          "measurement.** Oxygen does not enter its kill term; what changes "
          "is the suppressor field it meets, and that coupling is asserted by "
          "the model rather than fitted. A redirected T cell does not need "
          "oxygen to lyse a target, which is why the arm loses less than the "
          "ROS arms and more than ablation.", "",
          "**Ablation's flat row is the model being consistent, not a "
          "result.** A threshold arm is unaffected by every axis here by "
          "construction. Reporting it as robustness would be circular; it is "
          "in the table because leaving it out would make the other rows look "
          "like a full account of the space.", "",
          "**Every magnitude is uncalibrated** except radiation's DNA "
          "channel. `simulations/calibration/CALIBRATION_STATUS.md` carries "
          "the accounting; the ORDERING is the result, and the numbers are "
          "not.", ""]
    return "\n".join(L) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--render-only", action="store_true")
    a = ap.parse_args()
    d = assemble(json.loads(OUT_JSON.read_text()) if a.render_only else scan())
    OUT_JSON.write_text(json.dumps(d, indent=1) + "\n")
    OUT_MD.write_text(render(d))
    print(f"wrote {OUT_MD}")
    print(f"wrote {OUT_JSON}")
    print(f"  live axes: {d['live_axes']}  inert: {d['inert_axes']}")
    for a_ in d["robustness_order"]:
        worst = max((d["effects"][l][a_] for l in d["live_axes"]),
                    key=abs, default=0.0)
        verb = "loss" if worst < 0 else "GAIN"
        print(f"    {a_:22s} largest relative {verb} {abs(worst) * 100:5.1f}%")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
