# Distilled Hypotheses: What Survives Ruthless Scrutiny

## What we eliminate

**"Ferroptosis is a cancer vulnerability"** — Already a 2,045-citation review (PMID:35338310). Not novel.

**"SDT triggers ferroptosis via ROS"** — Already in dozens of papers. Not novel.

**"ICD can follow ferroptosis"** — Already established. Not novel.

**"The resistance tradeoff principle"** — As a general concept ("resistance creates new vulnerabilities"), this is known. Synthetic lethality is essentially this idea applied to genetics. The general principle is a useful framing but not a novel discovery.

**"ΔΨm as biomarker"** — Speculative, clinically unmeasurable, and the 392 articles are about basic science, not predictive biomarkers. Interesting but weak. Eliminate.

## What survives: Two hypotheses

---

### HYPOTHESIS A: SDT is categorically different from other physical modalities at the molecular level, and this difference has been invisible because nobody has compared them on ferroptosis engagement.

**Core idea in one sentence**: Among all physical cancer therapies, SDT is the only one that kills cells through a ferroptosis-dependent biochemical cascade rather than direct physical destruction, making it uniquely positioned to trigger immunogenic cell death — a distinction that has not been made in any published comparison.

**Facts it rests on** (all verified from corpus):
- SDT: 39 ferroptosis articles, 21 ICD articles, 72 GSH/GPX4 articles, 7 dual-pathway articles
- TTFields: 0 ferroptosis, 10 ICD, 0 GSH/GPX4
- HIFU: 1 ferroptosis, 4 ICD, 1 GSH/GPX4
- IRE: 1 ferroptosis, 11 ICD, 2 GSH/GPX4
- Zero papers in the corpus directly compare physical modalities on ferroptosis engagement
- 3 SDT+OXPHOS papers exist showing SDT targets mitochondrial metabolism (PMID:29555321, 35280333, 38849886)

**The hidden conceptual leap**: SDT is currently classified alongside TTFields, HIFU, and IRE as "a physical therapy." This classification is biologically wrong. SDT uses physical energy to initiate a biochemical program (ROS → GSH → GPX4 → ferroptosis); the others use physical energy for direct cell killing. This misclassification has caused SDT to be evaluated as an ablative technology (competing with HIFU, radiation) instead of as a biochemical modulator that happens to be triggered by ultrasound.

**Why it could matter if true**: SDT would have a unique clinical niche: not ablation (where it can't compete), but sub-ablative immune priming through ferroptotic ICD. No other physical modality can do this because no other physical modality engages the ferroptosis pathway.

**Why it is not obvious from the literature**: Each individual fact is published. But the SDT-ferroptosis papers are in materials chemistry / nanotechnology journals. The ICD papers are in immunology journals. The physical therapy comparison papers don't examine ferroptosis at all. The synthesis requires looking across these literatures simultaneously — which is what a corpus-wide analysis can do but individual researchers typically don't.

**Breakthrough type**: Translational framing. The biology isn't new. The reframing — "SDT is not a physical therapy, it's a physically-triggered biochemical therapy" — is new and would change what trials are designed and how SDT is evaluated.

**Biggest flaw**: The 39 SDT-ferroptosis articles are almost entirely nanoparticle engineering papers where ferroptosis induction was *designed into the nanosonosensitizer*, not an inherent property of ultrasound+sonosensitizer alone. The ferroptosis engagement may be an artifact of intentional nanoparticle design rather than an inherent SDT property. If plain SDT (without GSH-depleting nanoparticles) doesn't trigger ferroptosis, the hypothesis narrows from "SDT is unique" to "engineered nanoSDT platforms are unique" — still interesting but less fundamental.

---

### HYPOTHESIS B: OXPHOS-switched resistant cells should be selectively vulnerable to SDT-induced ferroptosis, and this specific sequence has never been tested despite both halves being well-documented.

**Core idea in one sentence**: Therapy-resistant cancer cells that survive by switching to oxidative phosphorylation acquire the exact molecular preconditions (high iron, high ROS, lipid membranes) that make them maximally susceptible to SDT-triggered ferroptosis — meaning the resistance mechanism itself becomes the vulnerability.

**Facts it rests on**:
- 61 articles document OXPHOS as a resistance mechanism
- OXPHOS requires iron-sulfur clusters (high iron demand), generates mitochondrial ROS, expands lipid membranes
- These are precisely the preconditions for ferroptosis sensitivity
- Only 4 papers in 10,413 connect OXPHOS to ferroptosis directly
- Only 2 of those 4 mention resistance
- PMID:29555321 shows blocking glycolysis sensitizes breast cancer to SDT — directly relevant
- PMID:33408185 shows mitochondrial metabolic reprogramming controls ICD induction
- Zero papers test "OXPHOS-resistant cells → SDT → ferroptosis → ICD" as a designed sequence

**The hidden conceptual leap**: Current oncology treats resistance as a problem to overcome. This hypothesis treats resistance as a *diagnostic signal* that tells you which vulnerability to exploit next. The OXPHOS switch isn't a dead end — it's an arrow pointing at ferroptosis.

**Why it could matter if true**: It would create a new therapeutic paradigm: anticipatory sequential therapy. First-line treatment → detect resistance phenotype → deploy the therapy matched to that phenotype's tradeoff. For the specific OXPHOS→ferroptosis case, it would give SDT its first rational clinical indication: post-resistance second-line immune priming.

**Why it is not obvious from the literature**: The OXPHOS-resistance and ferroptosis literatures exist in different communities (metabolism vs. cell death). 4 bridge papers out of 10,413 = 0.04% bridging rate. Neither community routinely cites the other.

**Breakthrough type**: Therapeutic + systems-level. The specific sequence is therapeutic. The deeper idea — "resistance phenotype predicts next vulnerability" — is systems-level and would generalize beyond this instance.

**Biggest flaw**: OXPHOS-resistant cells may co-upregulate NRF2/antioxidant defenses alongside OXPHOS, which would block ferroptosis even as the metabolic preconditions accumulate. If cells can have OXPHOS AND strong antioxidant defense simultaneously, the tradeoff doesn't exist and the hypothesis fails. This is testable in the isogenic cell pair experiment.

---

## Verdict

**Hypothesis A** (SDT reclassification) is an insight of moderate novelty — it synthesizes known facts into a new comparison, but relies on nanoparticle engineering literature that may not generalize to clinical SDT. Worth publishing as a perspective.

**Hypothesis B** (OXPHOS→ferroptosis→SDT sequence) is the higher-upside idea — it connects two large, separate literatures through only 4 bridge papers and proposes a testable therapeutic sequence that nobody has designed. If the isogenic cell pair experiment confirms differential ferroptosis sensitivity, it immediately justifies a clinical trial. This is the idea that could actually change what scientists test.

**The combination of A+B** is stronger than either alone: SDT is unique among physical modalities because it engages ferroptosis (A), and this uniqueness has maximal clinical value in the post-resistance OXPHOS context (B). The paper should present B as the primary hypothesis, with A as the mechanistic explanation for why SDT is the right tool.
