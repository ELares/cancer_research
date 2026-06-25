# Research Roadmap (2026): What's Next

This is a handwritten interpretation note (not a generated artifact). It records a deep,
web grounded scan of the current literature and the computational oncology methods
landscape, run in June 2026, to answer a simple question: now that the corpus
consolidation and the simulation suite exist, what would make this a genuine, citable
contribution to the scientific community, and what are we missing, getting wrong, or
under representing?

The scan covered ten axes (ferroptosis biology, sonodynamic/photodynamic therapy,
immunotherapy and immunogenic cell death, computational oncology methods, drug PK and
transport, persisters and resistance evolution, under covered modalities and taxonomy,
calibration and validation, corpus and evidence methodology, and scientific contribution
strategy). It produced the issue backlog indexed at the end of this document (issues #330
through #354).

## A note on fact grounding

Every quantitative or citation claim that reaches an issue or this document was checked
against a primary source. During the scan, several machine proposed PubMed IDs resolved
to unrelated papers when verified against NCBI (a known failure mode), so the rule here
is: cite only PMIDs confirmed via NCBI esummary, or name public databases by their stable
URL. The verified references used across the backlog are listed at the bottom.

## The central theme: from uncalibrated scaffolding to validated science

`simulations/calibration/CALIBRATION_STATUS.md` is honest that nearly every layer of the
simulation suite is "uncalibrated (illustrative)": the mechanism and direction are the
claim, the magnitude is not, and the whole suite is deliberately excluded from the
manuscript's quantitative results. That honesty is a strength, but it also defines the
single highest value direction. The gap is not code maturity, it is validation maturity.
The largest cluster of high priority issues below (the validation and calibration epic)
is about turning the scaffolding into something that predicts held out data, reports
uncertainty intervals rather than point estimates, and can be cited as a method.

## Findings by axis

### Ferroptosis biology

State of the art: ferroptosis is now understood as a multi axis system. Our biochem model
captures the System xc/GSH/GPX4 axis and the FSP1/CoQ10 backup, but omits the other two
established GPX4 independent defense axes, DHODH and GCH1/BH4, and omits the lipid
remodeling escape routes (kinetic SCD1/MUFA, ether lipids via FAR1, MBOAT1/2). Our labile
iron pool is static, whereas in vivo it is dynamically set by transferrin import, ferritin
storage, and NCOA4 ferritinophagy. These omissions likely make the model overstate RSL3
monotherapy kill and understate context dependent resistance.

Contradiction worth flagging: we treat hypoxia as uniformly protective against
pharmacologic ferroptosis. The literature is more nuanced, because hypoxia also upregulates
iron scavenging (transferrin receptor up), which can increase labile iron and ferroptosis
vulnerability in some contexts. The direction of the hypoxia effect is therefore context
dependent, not monotone. Issues: #338 (DHODH/GCH1), #339 (lipid remodeling), #340 (dynamic
iron/ferritinophagy).

### Sonodynamic and photodynamic therapy

This is the most contested thread in our work, and the current literature sharpens the
caveat we already added to the manuscript. The lead clinical sonodynamic agent,
SONALA-001 (intravenous 5-ALA converted to protoporphyrin IX, activated by MR guided
focused ultrasound), is a Type II, oxygen dependent sensitizer, mechanistically like
photodynamic therapy. Its first in human recurrent high grade glioma results (Science
Translational Medicine 2025) reported safety and only modest cell death. Type I (oxygen
independent) sonosensitizers exist but are preclinical only. There is also a self limiting
hypoxia paradox: ROS generation consumes oxygen and can worsen local hypoxia during
treatment. Our model's optimistic, largely oxygen independent SDT is therefore not
representative of the agent actually in trials. Issue: #336 (oxygen dependent SDT yield).

### Immunotherapy and immunogenic cell death

Immunotherapy dominates the corpus (about 1,685 of 4,830 full text articles), yet our
immune model under represents it, and specifically models only the pro immune side of
ferroptosis. Recent in vivo evidence shows ferroptosis is frequently net immunosuppressive:
ferroptotic cancer cells impede dendritic cell mediated anti tumor immunity (Wiernicki,
Nature Communications 2022, PMID 35760796), extracellular GPX4 released by ferroptotic
cells impairs antitumor immunity via dendritic ZP3 receptors (Cell 2026, PMID 41494530),
and reviews describe ferroptosis as having dual, context dependent immune roles (Tang,
Immunological Reviews 2024, PMID 37424139). Our SDT to RSL3 immune ratio (about 104:1 in
2D, shrinking to about 4:1 in 3D) captures only the amplifying side; the sign itself can
flip. Issue: #337 (immunosuppressive axis of ferroptosis). The magnitude of the ratio was
already addressed by the closed issue #288; #337 addresses the missing mechanism and sign.

### Computational oncology methods

The field solves reaction diffusion PDEs for oxygen, nutrient, and drug fields over
explicit vessel sources with consumption (PhysiCell, CompuCell3D, Chaste, Morpheus), which
reproduces non monotonic, biphasic gradients with high drug, low oxygen pockets. Our edge
distance and Krogh exponential proxies average those away. We already ship a PhysiCell
C-FFI, so benchmarking and interoperation are within reach. Issues: #343 (reaction
diffusion fields), #344 (ODE cross validation), #351 (SBML export).

### Drug PK and transport

`tumor_pk` and `drug_transport` use order of magnitude estimates. Real intratumoral
penetration is shaped by the binding site barrier and interstitial transport, not pure
exponential decay; antibody and nanoparticle delivery face an even steeper barrier and the
EPR effect is far weaker than once assumed. RSL3 has no clinical PK, so a named preclinical
analog (for example imidazole ketone erastin) should anchor the curve. Issues: #334
(PBPK/PK anchor, superseding the closed #316), #335 (penetration and binding site barrier).

### Persisters and resistance evolution

Our persister model captures recovery and a competing rate acquire/revert step, but not
the triggers of entry into the slow cycling state, nor the empirical reversible to
irreversible (epigenetically locked) transition after sustained drug exposure. FSP1 and
HDAC mediated suppression of alternative defenses in persisters is now documented (Science
Advances 2026, PMID 41481741, already cited in the manuscript). Therapy induced senescence
is listed as a primary escape route in `analysis/principle-resistance-tradeoff.md` but is
not simulated at all. Issues: #341 (senescence/SASP), #342 (persister entry and locking).

### Under covered modalities and taxonomy

The manuscript already names a next pass taxonomy scaffold. Executing it makes the
consolidation current. The verified landmark recovery targets include VISION (Lu-177-PSMA
for metastatic castration resistant prostate cancer, NEJM 2021, PMID 34161051) and
PANOVA-3 (tumor treating fields plus chemotherapy for pancreatic cancer, JCO 2025, PMID
40448572), both known absent from the local full text corpus and both distorting the
maturity picture for their mechanisms. Issues: #345 (landmark recovery), #347 (taxonomy
expansion to radioligand, targeted protein degradation, oncolytic, mRNA vaccine, TTFields,
bispecifics, cuproptosis, disulfidptosis, cold atmospheric plasma).

### Calibration and validation (the keystone)

Public datasets exist that we have not touched: GDSC (https://www.cancerrxgene.org/),
DepMap (https://depmap.org/), and the PRISM repurposing screen (Corsello, Nature Cancer
2020, PMID 32613204), all with measured drug responses across hundreds of cell lines,
including ferroptosis relevant compounds. The discipline of the field is held out
prediction, global sensitivity (Sobol) and identifiability analysis, and uncertainty
propagation (Bayesian/ABC) so that claims are reported as intervals. The existing
univariate PRCC work (#134, closed) is a good foundation to extend. Issues: #330 (GDSC/
DepMap/PRISM calibration), #331 (Sobol and identifiability), #332 (Bayesian/ABC intervals),
and #333 (multi size spheroid validation against Browning, eLife 2021, PMID 34842141).

### Corpus and evidence methodology

The keyword tagger is 96% binary evidence-presence precision but only 55% recall, so absence claims are provisional
and rare mechanisms are under counted. MeSH hierarchical expansion and embedding based
semantic retrieval are the standard upgrades; both can be re measured on the existing 100
article gold set. The open access fraction of the full text corpus is about 98.7%, which
may distort mechanism rankings, and this is quantifiable. Issues: #346 (recall via MeSH and
embeddings), #348 (open access bias), #349 (living review pipeline).

### Scientific contribution strategy

To be cited and built upon, the work needs a stable archived release with a DOI, a
reproducibility container, a standard format export for interoperability, transparent
reporting (a model card and assumptions checklist), and a preregistered, falsification
oriented validation plan that splits the corpus synthesis and the simulation methods into
their two natural contributions. Issues: #350 (Zenodo DOI and container), #351 (SBML
export), #352 (model card and reporting standards), #353 (preprint and preregistration
plan), #354 (interactive dashboard).

## Contradictions to our current claims (summary)

1. SDT is oxygen dependent in the clinic (SONALA-001 is Type II), so the hypoxia leg of
   the physical versus pharmacologic asymmetry is the weakest of the three and may close
   for the agent actually in trials. Addressed by #336; the manuscript caveat already
   reflects this and is now backed by a 2025 clinical citation.
2. Ferroptosis is often net immunosuppressive in vivo (extracellular GPX4, suppressive
   myeloid enrichment), so the immune coupling claim is on shakier ground than the model's
   pro immune only treatment implies. Addressed by #337.
3. Hypoxia is not uniformly protective against pharmacologic ferroptosis; iron scavenging
   can increase vulnerability. Addressed by #340.
4. The steady state MUFA assumption (about 0.40) does not hold for acute dosing, where
   SCD1 driven protection has not yet accumulated. Addressed by #339.
5. Real oxygen and drug gradients are biphasic, not monotone exponential. Addressed by
   #343.

## What solidifies our findings

1. The manuscript's de biased caveat that the lead SDT agent is Type II and oxygen
   dependent is directly supported by the SONALA-001 first in human results (Science
   Translational Medicine 2025).
2. The §7.2 caveat that ferroptosis immunogenicity is contested and possibly
   immunosuppressive is strongly supported by 2022 to 2026 evidence (PMID 35760796, PMID
   41494530, PMID 37424139).
3. The dual GPX4/FSP1 dependency that underlies the 1.99x Bliss synergy is consistent with
   FSP1 being a critical backup and with persister ferroptosis suppression via FSP1 and
   HDACs (PMID 41481741).
4. The spheroid zone geometry is already anchored to measured structure (Browning, eLife
   2021, PMID 34842141).

## The issue backlog (epics)

**Validation and calibration (keystone):** #330, #331, #332, #333.
**Drug PK and transport:** #334, #335.
**Model biology gaps:** #336 (oxygen dependent SDT), #337 (immunosuppressive ferroptosis),
then #338 (DHODH/GCH1), #339 (lipid remodeling escape), #340 (dynamic iron/ferritinophagy),
plus #341 (senescence/SASP) and #342 (persister entry and locking).
**Spatial and numerical methods:** #343 (reaction diffusion fields), #344 (ODE cross
validation).
**Corpus and evidence:** #345 (landmark recovery), #346 (recall), #347 (taxonomy
expansion), #348 (open access bias), #349 (living review).
**Scientific contribution and reproducibility:** #350 (Zenodo DOI and container), #351
(SBML export), #352 (model card), #353 (preprint and preregistration), #354 (dashboard).

## Verified key references

- SONALA-001 sonodynamic therapy, first in human recurrent high grade glioma. Science
  Translational Medicine 2025, DOI 10.1126/scitranslmed.ads5813.
- Wiernicki et al. Cancer cells dying from ferroptosis impede dendritic cell mediated anti
  tumor immunity. Nature Communications 2022, PMID 35760796.
- Extracellular GPX4 impairs antitumor immunity via dendritic ZP3 receptors. Cell 2026,
  PMID 41494530.
- Tang et al. Ferroptosis in immunostimulation and immunosuppression. Immunological Reviews
  2024, PMID 37424139.
- FSP1 and histone deacetylases suppress cancer persister cell ferroptosis. Science
  Advances 2026, PMID 41481741.
- VISION: Lutetium-177-PSMA-617 for metastatic castration resistant prostate cancer. NEJM
  2021, PMID 34161051.
- PANOVA-3: Tumor treating fields with gemcitabine and nab-paclitaxel for locally advanced
  pancreatic adenocarcinoma. JCO 2025, PMID 40448572.
- Corsello et al. Discovering the anti-cancer potential of non-oncology drugs by systematic
  viability profiling (PRISM). Nature Cancer 2020, PMID 32613204.
- Browning et al. Quantitative analysis of tumour spheroid structure. eLife 2021, PMID
  34842141.
- GDSC, https://www.cancerrxgene.org/ . DepMap, https://depmap.org/ .
