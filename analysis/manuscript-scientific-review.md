# Manuscript Scientific Review (deep cross-reference + external literature)

A from-scratch scientific review of the manuscript (`article/drafts/v1.md`),
cross-referencing our own corpus/analysis/simulation materials AND the live
external literature, to check that the math, physics, and biology are correct
and that our claims survive contact with what is published elsewhere. Conducted
with five independent lenses: (1) external literature that *supports* our
claims, (2) external literature that *contradicts* them, (3) math/physics
verification against the simulation code, (4) fresh-eyes internal coherence +
staleness, (5) where a figure would strengthen the argument.

Every PMID/DOI cited below was retrieved and read from PubMed / the journal page
(not inferred from a search snippet). Nothing here is a fabricated citation.

## Headline

The manuscript is unusually self-critical and survived the review well. **No
fabricated citations; the math and physics are sound; the central supporting
citation (Higuchi 2026) is real and accurately used.** The one genuinely
important scientific issue is that our **most exposed claim is the central one**:
"physical ROS *bypasses* hypoxia." That framing runs against the mainstream
sonodynamic-therapy (SDT) literature, and the manuscript previously conceded it
only in a buried caveat that the chapter-level framing then overrode. We have
now rebalanced it (changes listed below). Everything else is supporting
citations, precision fixes, and a figure roadmap.

## External literature — what SUPPORTS us (verified)

| Claim | Strongest verified support | PMID/DOI | Strength |
|-------|----------------------------|----------|----------|
| Persister cells acquire FSP1/HDAC/OXPHOS-ROS-dependent ferroptosis vulnerability (our §6.2/§6.3/§7.1 anchor) | Higuchi M et al., *Sci Adv* 2026 — citation verified exact; also shows mito-ROS is *required* for persister ferroptosis (mitoTEMPO rescues), independently supporting our §7.1 no-O2→no-trigger mechanism | PMID 41481741 | High |
| Hypoxia confers ferroptosis resistance by suppressing mitochondrial ROS (PDAC, in vivo) — supports the §7.1 RSL3-collapse mechanism | Hubbi ME … Dang CV, *Mol Cell* 2026 | PMID 41932308 | High — **added to §7.1** |
| Ferroptosis as a cancer target; RSL3 inactivates GPX4 → lipid-ROS | Zhou Q et al. *STTT* 2024 (PMID 38453898); Sui X et al. *Front Pharmacol* 2018 (PMID 30524291) | — | High |
| FSP1 is a parallel GPX4-independent axis; GPX4+FSP1 co-inhibition synergizes | Bersuker (PMID 31634900) + Doll (PMID 31634899), *Nature* 2019; Tamura et al. *Int J Mol Sci* 2024 PMID 39273151 (a prognostic/target study that also reports in-vitro RSL3+iFSP1 synergistic, ferroptosis-inhibitor-blockable cell death in gastric lines) | — | High |
| SDT/PDT induce ferroptosis + ICD markers | Theranostics 2021 (PMID 33408790, SDT-ferroptosis in hypoxic tumors); PDT-ICD PMID 31602649 | — | High (induction) |

## External literature — what CONTRADICTS / COMPLICATES us (verified)

Ranked by how damaging, most first. These drove the manuscript revisions.

1. **"Physical ROS bypasses hypoxia" contradicts the SDT field consensus.** SDT's
   ROS generation is widely treated as O2-dependent; an entire SDT subfield
   exists to *deliver oxygen* to overcome hypoxic SDT failure; the SDT mechanism
   itself is "not fully elucidated" (Dong HQ et al., *World J Clin Cases* 2023,
   **PMID 37621595**). Decisively, the **leading clinical SDT agent (SONALA-001,
   5-ALA→protoporphyrin IX) is a Type II, oxygen-dependent sensitizer** (Sanai N
   et al., *Sci Transl Med* 2025, **PMID 41296829**, DOI 10.1126/scitranslmed.ads5813,
   NCT04559685) — not the Type I O2-independent mechanism we model. → **Reframed §7.1 caveat (now the
   chapter's headline caveat), §7.5, §8.4, §10.1.**
2. **DAMP *balance*, not quantity, governs immunogenicity** (Hayashi K et al.,
   *Nat Commun* 2020, **PMID 33288764**): gemcitabine releases abundant
   stimulatory DAMPs yet fails ICD because it co-releases inhibitory PGE2. Our
   §7.2 kill-density→immunity logic assumes quantity. → **Added to §7.2 caveat.**
3. **Ferroptosis can be pro-tumorigenic** — and the strongest in-vivo evidence is
   in pancreatic cancer, a named target: GPX4 loss / high iron accelerates
   KRAS-driven PDAC via 8-OHG/STING/macrophages (Dai E et al., *Nat Commun*
   2020, **PMID 33311482**). The Wiernicki 2022 result we already cite is also
   stronger than we represented (DC vaccination with ferroptotic cells *failed*
   to protect, vs necroptotic cells). → **Expanded §7.2 caveat.**
4. **"RSL3 fails in hypoxia" is bidirectional / tumor-type-specific.** In
   clear-cell renal carcinoma, HIF-2α–HILPDA enriches PUFAs and *sensitizes*
   cells to GPX4 inhibitors (Zou Y et al., *Nat Commun* 2019, **PMID 30962421**)
   — the opposite of our PDAC-aligned mechanism. → **Added to §7.1.**

## Manuscript revisions made in this pass

All in `article/drafts/v1.md` (and `v1.tex` regenerated):
- **§7.1** — added in-vivo support for the RSL3-collapse mechanism (Hubbi/Dang
  2026) and the bidirectional caveat (Zou 2019); **rewrote the O2-independence
  caveat** to be the chapter's headline caveat, naming the Type II clinical SDT
  agent and demoting "bypasses hypoxia" to a hypothesis. Corrected the one-sided
  "hypoxic protection is likely STRONGER" line to "cuts both ways."
- **§7.2** — added the inhibitory-DAMP-balance challenge (Hayashi 2020) and the
  pro-tumorigenic-ferroptosis evidence (Dai 2020, PDAC), and strengthened the
  Wiernicki framing.
- **§7.5 / §8.4 / §10.1** — reconciled the "bypass all three barriers" /
  "directional findings" / summary language so the hypoxia leg is consistently
  flagged as the contested one (stromal + pH bypass is the stronger claim).
- **`parameter_provenance.md` + `params.rs`** — fixed the `iron_diffusion_coeff`
  mis-citation (was wrongly attributed to Jacques 2013, the optics ref; it is a
  tortuosity-reduced estimate). Moved from "Grounded" to "Assumed."

## Math / physics verification — result: SOUND

Checked against the simulation source and standard references. All equation
forms correct, units dimensionally consistent, headline numbers reproduce from
stated parameters, no sign/mass-balance errors:
- Optical (PDT): `μ_eff = √(3·μ_a·(μ_a+μ_s'))` diffusion-approx, μ_eff=0.31/mm →
  δ≈3.2 mm — matches measured 630 nm values. Correct.
- Acoustic (SDT): `I=I₀·10^(−α·f·z/10)`, α=0.7 dB/cm/MHz — correct form; value is
  high end of soft-tissue range (consensus ~0.5), i.e. *conservative* for the
  penetration thesis.
- Krogh / O2 / pH exponential fields — correct steady-state reaction-diffusion
  approximations, honestly caveated as approximations of the exact (Bessel /
  Riley) solutions.
- Bliss 1.99×: `E_A+E_B−E_A·E_B` expected, observed/expected ratio — correct
  arithmetic (0.40+0.037−0.40·0.037=0.4222; 0.841/0.4222=1.99).
- Ferroptosis ODEs, Michaelis-Menten DAMP→immune, photosensitizer PK — all
  sensible phenomenological forms with non-negativity floors; rates estimated
  (already disclosed).
- **Only issue:** the `iron_diffusion_coeff` citation error, now fixed.

## Figure roadmap (where a graph would help)

Calibration status flagged so we never plot an uncalibrated number as
quantitative. **Tier 1 (calibrated, safe to plot quantitatively):**
1. **Hypoxia kill-collapse vs O2 penetration depth** — RSL3 vs SDT line plot
   over λ = 80/100/120/150 µm. The lynchpin mechanism, currently one sentence.
   Data: `sim-tme`. (After this pass, the figure caption must carry the SDT-O2
   caveat.)
2. **Bliss synergy** (RSL3 / FSP1i / combination observed-vs-Bliss-expected bar
   chart). Data: `sim-combo-mech`. Highest value-to-effort.
3. **Depth-kill curves** PDT vs SDT (Beer-Lambert vs acoustic). Data: `sim-tme`.

**Tier 2 (medium confidence, plot with caveats):** vulnerability-window kinetics
(`sim-window`); immune-coupling spatial DAMP (note 104:1→4:1 in 3D — present the
2D calibrated value or caption the 3D discrepancy).

**Tier 3 (defer until calibrated):** pH ion-trapping (RSL3 pKa is the model's
most uncertain parameter — defer until experimental validation). 3D-suite
overlays: illustrative/mechanistic only, never quantitative (per §8.4).

**Single highest-value addition:** an integrated 2×2 "resistance-mechanism
asymmetry" figure (hypoxia / stromal / pH / immune, RSL3 vs SDT) — but it must
NOT entrench the contested hypoxia leg; build it *after* the reframe and caption
each panel with its confidence tier.

Figure generation was deliberately **not** done in this pass: the highest-value
figure embeds the very claim we just rebalanced, so settling the science first
is correct. Generating the Tier-1 figures (matplotlib scripts → sim runs →
FIGURES.yaml + generate_latex.py wiring → traceability test) is the clear next
work item.

## Open items for author judgment

- **Generate the Tier-1 figures?** (hypoxia collapse, Bliss synergy, depth-kill).
  Recommended; needs new figure scripts + FIGURES.yaml entries.
- **Persister phenotype homogeneity (minor).** Higuchi 2026 shows OXPHOS-stratified
  heterogeneity *within* the persister pool; our 16-condition matrix treats
  "persister" as one fixed FSP1-low phenotype. Worth a sentence in §5.2 if we
  want to pre-empt it.
- **Pancreatic-target framing.** Given Dai 2020 (ferroptosis accelerates PDAC in
  vivo), reconsider whether pancreatic cancer should remain a *lead* named
  context or be presented with the pro-tumor caveat attached at first mention
  (Chapter 9), not only in §7.2.
