"""The front door's description of the engine, bound to the engine.

README.md said the worked implementations were "ferroptosis/RSL3 biochemistry
and PDT/SDT depth physics -- there are no others" long after eight more arms
had landed, and named radiotherapy as its worked example of a modality with
"no query, no mechanism tag, no engine term" after `Treatment::Radiation`
shipped. Both were true when written. Neither was guarded, so neither moved.

WHAT THIS FILE LEARNED FROM ITS OWN FIRST VERSION. That version asserted
eleven substrings were PRESENT and called itself a check that the front door
describes the engine. A mutation sweep survived 11 of 11, and the sharpest
survivor is the one that matters: swapping which side of "against" the two
depth figures sit on -- so the README would read that the FERROPTOSIS ENGINE
is 5 modules and 366 lines against the modality arms' 27 and 3,988, the exact
inverse of the finding -- left every assertion green, because both strings
were still in the file. A number's presence is not the claim; a number is a
claim only once it is ATTACHED to a subject and a direction.

So every check here parses the sentence and asserts WHAT the figure is said
about: which side of a contrast a count sits on, which pool a share is a share
of, which arm kills which fraction. Where the sentence makes a claim about the
CODE (that two arms share one payload constant), the code is read too.

NOT a check that the prose is good. It cannot tell whether the description is
fair, only whether its figures are the live ones and are attached to the
subjects they are asserted about.
"""
import json
import re
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
README = REPO / "README.md"
ANALYSIS = REPO / "analysis"
PANEL_SRC = REPO / "simulations" / "sim-modality-panel" / "src" / "main.rs"


WORDS = {1: "one", 2: "two", 3: "three", 4: "four", 5: "five", 6: "six",
         7: "seven", 8: "eight", 9: "nine", 10: "ten", 11: "eleven",
         12: "twelve", 13: "thirteen"}

# Every live arm, and the words the README is required to name it by. The map
# is here rather than derived because an enum variant is not English; what it
# must not be is a free choice, so a missing arm fails rather than being
# skipped.
ARM_PROSE = {
    "RSL3": "ferroptosis induction",
    "PDT": "photodynamic",
    "SDT": "sonodynamic",
    "Radiation": "ionizing radiation",
    "Chemotherapy": "cytotoxic chemotherapy",
    "Immunotherapy": "checkpoint",
    "AdoptiveCell": "adoptive cell therapy",
    "OncolyticVirus": "oncolytic virus",
    "Ablation": "ablation",
    "AntibodyDrugConjugate": "antibody-drug conjugate",
}


def _j(name):
    return json.loads((ANALYSIS / name).read_text())


def _num(n: int) -> str:
    """A count as a digit OR spelled out, as one alternation.

    Prose reads better with words and artifacts carry digits; forcing one on
    the other makes the guard an editor rather than a check.
    """
    return "(?:%s|%s)" % (re.escape(f"{n:,}"), re.escape(WORDS.get(n, str(n))))


def _readme() -> str:
    # Wrapping is a layout choice and every claim here spans lines; collapsing
    # it means a reflow cannot break a guard, which is what makes the guards
    # about content.
    return re.sub(r"\s+", " ", README.read_text())


def test_the_readme_names_every_live_arm_it_counts():
    """The count and the list have to agree with the panel AND each other."""
    txt = _readme()
    binaries = sorted(p.name for p in (REPO / "simulations").glob("sim-*")
                      if p.is_dir())
    assert f"**{len(binaries)} Rust simulation binaries**" in txt, (
        f"the README does not say {len(binaries)} binaries; the tree holds "
        f"{binaries}")
    arms = [a["arm"] for a in _j("modality-panel.json")["arms"]]
    live = [a for a in arms if a != "Control"]
    assert re.search(
        r"\*\*%s treatment arms plus an untreated control\*\*" % _num(len(live)),
        txt), (f"the README does not state the live arm count ({len(live)} plus "
               f"a control); `sim-modality-panel` runs {sorted(arms)}")
    # The count is only as good as the list beneath it: a dropped arm leaves a
    # number that still matches while the sentence under-describes the engine.
    for arm in live:
        assert arm in ARM_PROSE, (
            f"{arm} is a live panel arm with no README wording registered; add "
            "it to ARM_PROSE and to the README sentence")
        assert ARM_PROSE[arm] in txt, (
            f"the README counts {len(live)} arms but never names {arm} "
            f"(expected the words {ARM_PROSE[arm]!r})")


def test_the_depth_gap_is_attached_to_the_right_side_of_the_contrast():
    """The honest half, and the one a swap silently inverts.

    Asserting both figures are PRESENT passes on a sentence claiming the
    ferroptosis engine is the small side. So the sentence is split at its own
    contrast and each figure is required on its own side of it.
    """
    txt = _readme()
    d = _j("modality-module-depth.json")
    m = re.search(
        r"the arms a modality owns outright are (?P<arms>.{0,120}?) "
        r"against the ferroptosis engine's (?P<engine>.{0,120}?), so ", txt)
    assert m, ("the README no longer contrasts the modality arms against the "
               "ferroptosis engine in the shape this guard can read; if the "
               "sentence was rewritten, re-derive the check rather than "
               "loosening it")
    want_arms = (f"{d['dedicated_modules']} modules and "
                 f"{d['dedicated_code_lines']} lines")
    want_engine = (f"{d['ferroptosis_engine_modules']} modules and "
                   f"{d['engine_code_lines']:,}")
    assert want_arms in m.group("arms"), (
        f"the modality-arm side of the contrast does not read {want_arms!r} "
        f"(it reads {m.group('arms')!r})")
    assert want_engine in m.group("engine"), (
        f"the ferroptosis-engine side of the contrast does not read "
        f"{want_engine!r} (it reads {m.group('engine')!r})")
    lo, hi = round(d["line_ratio_narrow"]), round(d["line_ratio_wide"])
    assert re.search(r"roughly %s to %s times smaller" % (_num(lo), _num(hi)),
                     txt), (
        f"the README says something other than the live {lo}-to-{hi}x gap "
        f"(exact ratios {d['line_ratio_narrow']} and {d['line_ratio_wide']})")
    assert d["engine_code_lines"] > d["dedicated_code_lines"], (
        "the modality arms now exceed the ferroptosis engine; the README's "
        "'smaller' framing has to be re-derived, not adjusted")


def test_the_delivery_finding_names_which_arm_kills_which_fraction():
    txt = _readme()
    panel = _j("modality-panel.json")
    arms = {a["arm"]: a["kill_fraction"] for a in panel["arms"]}
    adc, sdt = arms["AntibodyDrugConjugate"], arms["SDT"]
    assert (f"an antibody-drug conjugate kills **{adc * 100:.1f}%** where "
            f"sonodynamic therapy kills **{sdt * 100:.1f}%**") in txt, (
        f"the README's delivery finding does not attach the live pair to its "
        f"arms (ADC {adc * 100:.1f}%, SDT {sdt * 100:.1f}%)")
    assert adc < sdt, "the delivery argument depends on the ADC being smaller"
    # The clause that makes the contrast mean anything is a claim about the
    # SOURCE -- if the two arms ran different payloads the ratio would price a
    # chemistry difference, not transport. So the source is read.
    src = PANEL_SRC.read_text()
    for name in ("sdt_mean", "adc_mean"):
        line = re.search(r"let %s = (.+);" % name, src)
        assert line and "params.sdt_ros" in line.group(1), (
            f"{name} in sim-modality-panel no longer derives from "
            "params.sdt_ros, so the README's 'same exogenous-ROS constant' "
            "clause is false; fix the sentence, not this guard")
    assert "both driven by the same exogenous-ROS constant" in txt, (
        "the README no longer states that the two arms share their payload "
        "constant, which is what makes the ratio a transport measurement")
    free = arms["RSL3"]
    assert (f"kills **{free * 100:.2f}%** in the same table, LESS than the "
            "antibody-delivered arm") in txt, (
        f"the README does not attach the live free-RSL3 kill ({free * 100:.2f}%) "
        "to the arm it is smaller than")
    assert free < adc, (
        "the free-drug arm no longer kills less than the delivered one; that "
        "sentence is now false")


def test_the_car_t_collapse_is_attached_to_its_two_diseases():
    txt = _readme()
    ab = _j("modality-panel.json")["adoptive_barriers"]
    collapse = ab["leukaemia_kill_fraction"] / ab["solid_tumour_kill_fraction"]
    assert (f"The same CAR-T construct collapses {collapse:,.0f}-fold between "
            "two diseases") in txt, (
        f"the README does not attach the live {collapse:,.0f}-fold collapse to "
        "the construct it describes")
    assert (ab["leukaemia_kill_fraction"] > ab["solid_tumour_kill_fraction"]), (
        "the solid tumour no longer kills less than the leukaemia; 'collapses' "
        "is the wrong verb and the finding needs re-deriving")


def test_p10_names_the_pool_its_share_is_a_share_of():
    """P10 contradicts the project's own prior, so which pool shrinks is the
    whole content of it -- naming the other pool would report the belief the
    model refuted."""
    txt = _readme()
    preds = _j("modality-predictions.json")["P10"]["reach_of_negative_pool"]
    hi_k = max(preds, key=float)
    lo_k = min(preds, key=float)
    hi, lo = float(preds[hi_k]), float(preds[lo_k])
    assert hi > lo, (
        "reach no longer falls as antigen is lost; P10's sign has moved and "
        "the README sentence states the opposite of the artifact")
    # One decimal on both ends, matching PREREGISTRATION.md's own wording for
    # P10 -- the two files quoted the same number to different precisions.
    assert (f"the share of the antigen-negative pool it reaches falls from "
            f"{hi:.1%} to {lo:.1%} as antigen is lost") in txt, (
        f"the README does not attach P10's live decline ({hi:.1%} to {lo:.1%}) "
        "to the antigen-negative pool it is a share of")
    assert "starved by the escape it answers" in txt, (
        "the README no longer states P10's direction against the prior belief")


def test_the_scope_rows_name_the_modules_and_mechanisms_they_count():
    txt = _readme()
    sa = _j("scope-audit.json")
    silent = sorted(s.replace(".rs", "") for s in sa["engine_modules"]["silent"])
    m = re.search(r"Measured, %s modules \((?P<list>[^)]*)\) mention neither"
                  % _num(len(silent)), txt)
    assert m, (f"the README does not state the live silent-module count "
               f"({len(silent)}) with its list")
    named = sorted(re.findall(r"`([a-z_0-9]+)`", m.group("list")))
    assert named == silent, (
        f"the README's silent-module list is {named}, the audit's is {silent}")
    cov = _j("modality-coverage.json")
    absent = [m_ for m_ in cov["mechanisms"] if not m_.get("code_modules")]
    assert (f"**{len(absent)} now have no engine representation at all**"
            in txt), (
        f"the README does not state the live absent count ({len(absent)} of "
        f"{len(cov['mechanisms'])})")
    tiers = {}
    for row in cov["rows"]:
        tiers[row["engine_tier"]] = tiers.get(row["engine_tier"], 0) + 1
    assert (f"{tiers.get('treatment', 0)} of the {len(cov['rows'])} can "
            "be APPLIED as a treatment") in txt, (
        f"the README does not attach the live applicability split to the "
        f"treatment tier ({tiers})")
    assert (f"the other {tiers.get('modifier', 0)} are MODIFIERS") in txt, (
        f"the README does not state the live modifier count ({tiers})")
    cal = _j("modality-calibration.json")
    vc = cal["verdict_counts"]
    targeted = len(cal["arms"]) - vc.get("NO TARGET", 0)
    assert (f"Of the {targeted} arms that have a published calibration target "
            f"at all, {vc['ADMISSIBLE']} reproduce it") in txt, (
        f"the README does not attach the live verdict split to the arms that "
        f"have a target ({vc}, {len(cal['arms'])} arms)")


def test_the_readme_does_not_re_assert_the_retracted_claims():
    """Refused as PATTERNS, not as the exact sentences that were removed.

    The first version listed the four literal strings, so any rewording of the
    same claim would pass -- which is how the claims got here in the first
    place. Each pattern refuses the CLAIM.
    """
    txt = _readme()
    banned = [
        (r"(?:there are no others?|no other worked implementation)\b",
         "that ferroptosis and PDT/SDT are the only worked implementations"),
        (r"radiotherapy[^.]{0,160}?no (?:query|mechanism tag|engine term)",
         "that radiotherapy has no engine representation"),
        (r"\*\*(?!13 )\d+ Rust simulation binaries\*\*",
         "a binary count other than the live one"),
        (r"for \*\*ferroptosis and physical-ROS therapies specifically\*\*",
         "that the engine is specific to ferroptosis and physical ROS"),
        (r"different (?:payloads?|chemistr\w+)[^.]{0,80}(?:ADC|antibody-drug)",
         "that the ADC and SDT arms run different payloads"),
    ]
    for pattern, claim in banned:
        assert not re.search(pattern, txt, re.IGNORECASE), (
            f"the README asserts the retracted claim {claim} again "
            f"(matched {pattern!r})")
