# News Source Criteria and Authentication Framework

## 1. Purpose

The local corpus (4,830 PubMed-indexed articles) captures peer-reviewed research but misses real-world context that shapes how cancer research translates into practice:

- **Trial outcomes before publication** (6-18 month PubMed lag)
- **Regulatory milestones** (FDA approvals, EMA decisions)
- **Industry pipeline decisions** (trial failures, drug discontinuations)
- **Funding and policy shifts** (NCI priorities, WHO resolutions)
- **Patient-reported perspectives** (treatment access, quality of life)
- **Expert debate and synthesis** (commentary on emerging evidence)

Incorporating this material strengthens the manuscript's role as a comprehensive guide for the scientific community. However, non-peer-reviewed sources carry higher risk of inaccuracy, bias, and misinformation. This framework defines rigorous criteria for evaluating, verifying, and citing news sources.

---

## 2. Source Tiers

### Tier 1 — Institutional and Peer-Adjacent (high trust, cite directly)

Sources with institutional editorial oversight, fact-checking processes, and direct ties to primary research data. Claims from Tier 1 sources can be cited as supporting evidence when verified.

| Source | URL | Type |
|--------|-----|------|
| NIH/NCI Newsroom | cancer.gov/news-events | Government research agency |
| WHO Fact Sheets & News | who.int/news-room | International health authority |
| FDA Drug Approvals | fda.gov/drugs | Regulatory authority |
| ClinicalTrials.gov Announcements | clinicaltrials.gov | Trial registry |
| Nature News & Comment | nature.com/news | Journal publisher news arm |
| Science News | science.org/news | Journal publisher news arm |
| Cell Press News | cell.com/news | Journal publisher news arm |
| The Lancet News | thelancet.com/news | Journal publisher news arm |
| NEJM Journal Watch | jwatch.org | Journal publisher commentary |
| University press offices | (varies) | Institutional communications |
| IARC/GCO | gco.iarc.fr | Cancer statistics authority |

**Tier 1 criteria**: The source is an official communication channel of a government agency, international health organization, major research university, or peer-reviewed journal publisher. The content undergoes editorial review by the publishing institution. Authors or institutions are clearly identified.

### Tier 2 — Science Journalism (medium trust, cite with verification)

Professional journalism outlets with dedicated health/science desks, editorial standards, and named reporters. Claims require cross-referencing against the underlying primary source before citing.

| Source | URL | Focus |
|--------|-----|-------|
| STAT News | statnews.com | Biotech, pharma, health policy |
| The Cancer Letter | cancerletter.com | Oncology policy and clinical trials |
| Endpoints News | endpointsnews.com | Clinical trials and biotech |
| Reuters Health | reuters.com/business/healthcare-pharmaceuticals | Wire service health desk |
| Associated Press Health | apnews.com/health | Wire service health desk |
| Science Daily | sciencedaily.com | Research press release aggregator |
| Medical News Today | medicalnewstoday.com | Health news aggregator |
| Ars Technica Science | arstechnica.com/science | Technology-adjacent science |
| The Conversation (Health) | theconversation.com | Academic expert commentary |
| FiercePharma / FierceBiotech | fiercepharma.com | Industry-focused journalism |

**Tier 2 criteria**: The source employs professional journalists or editors with health/science expertise. Articles have bylines (named authors). The outlet has published corrections when errors are identified. Editorial independence from industry sponsors is maintained (or disclosed when not).

### Tier 3 — Expert Blogs and Commentary (use for context, not facts)

Individual expert voices, advocacy organizations, and institutional blogs. These provide valuable interpretation and context but should never be cited as factual evidence. Use for narrative framing, patient perspectives, and expert opinion.

| Source | URL | Voice |
|--------|-----|-------|
| Derek Lowe — In the Pipeline | science.org/blogs/pipeline | Pharmaceutical chemistry |
| Vinay Prasad (Substack/blog) | vinayprasad.com | Oncology trial criticism |
| ASCO Connection Blog | connection.asco.org | Clinical oncology |
| ACS Research News | cancer.org/research/acs-research-news | Advocacy + research |
| Leukemia & Lymphoma Society | lls.org/news | Patient advocacy |
| Broad Institute Blog | broadinstitute.org/blog | Computational biology |
| MD Anderson Cancerwise | mdanderson.org/cancerwise | Institutional perspective |
| Patient Power | patientpower.info | Patient experience |
| Cancer Research UK Blog | cancerresearchuk.org/about-cancer | UK advocacy + research |
| Institute for Cancer Research Blog | icr.ac.uk/blogs | UK research institution |

**Tier 3 criteria**: The author has verifiable expertise (medical degree, PhD, research position, or relevant professional experience). The source represents a recognized institution or has a track record of responsible commentary. Content is opinion/interpretation, not presented as peer-reviewed evidence.

### Excluded — Do Not Ingest

| Category | Examples | Rationale |
|----------|----------|-----------|
| Social media posts | Twitter/X, Reddit, Facebook, TikTok | No editorial oversight, high misinformation risk |
| Anonymous content | Unsigned articles, anonymous forums | No accountability |
| Alternative medicine sites | Natural News, Mercola, GreenMedInfo | Documented misinformation track records |
| Industry-sponsored "news" without editorial independence | Pharma-funded advertorials, sponsored content without disclosure | Conflict of interest without transparency |
| Preprint commentary without peer review | bioRxiv/medRxiv social media discussions | Preliminary, unverified findings |
| AI-generated content without expert review | Automated health articles, chatbot outputs | No human expert accountability |

---

## 3. Authentication Pipeline

For each candidate news article, apply these five steps sequentially:

### Step 1: Source Verification

| Check | How | Pass/Fail |
|-------|-----|-----------|
| Is the source in Tier 1-3? | Match URL domain against tier lists | Required |
| Is the author identifiable? | Look for byline, author bio, institutional affiliation | Required for Tier 1-2 |
| Is there an editorial process? | Check for editorial board, corrections policy, masthead | Required for Tier 1-2 |
| Is the publication date verifiable? | Check for date stamp, last-updated timestamp | Required |
| Is the content behind a paywall? | Check accessibility | Note (not disqualifying) |

### Step 2: Claim Extraction

Read the article and classify each assertion:

- **FACTUAL**: Contains a specific, verifiable number, statistic, trial name, drug name, or outcome (e.g., "Response rate was 43%", "The FDA approved drug X on date Y")
- **INTERPRETIVE**: Expert opinion, editorial judgment, or contextualization (e.g., "This represents a paradigm shift", "The results are encouraging")
- **SPECULATIVE**: Prediction, hypothesis, or unverified extrapolation (e.g., "This could lead to new treatments within 5 years", "If confirmed, this would change standard of care")

### Step 3: Cross-Reference Verification

For each FACTUAL claim:

1. Search the local corpus (INDEX.jsonl) by keywords (drug name, trial name, author name)
2. Search PubMed API: `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term={query}`
3. Search ClinicalTrials.gov if a trial name is mentioned

**Verification status per claim**:
- **VERIFIED**: Primary source found (PMID or DOI linked)
- **UNVERIFIED**: No primary source found; news article is the only source
- **DISPUTED**: Primary source found but contradicts the news claim
- **SELF-REFERENCING**: The news article IS the primary report (e.g., WHO fact sheet citing its own IARC data)

### Step 4: Credibility Scoring

Compute a 0-100 credibility score:

```
score = tier_weight × (40 × verified_ratio + 30 × author_score + 20 × recency + 10 × cross_citation)
```

| Component | Weight | Computation |
|-----------|--------|-------------|
| `tier_weight` | multiplier | Tier 1: 1.0, Tier 2: 0.8, Tier 3: 0.6 |
| `verified_ratio` | 0-1 | verified_claims / total_factual_claims |
| `author_score` | 0-1 | 1.0 if named + credentialed, 0.7 if named only, 0.3 if anonymous |
| `recency` | 0-1 | 1.0 if <6 months old, 0.8 if <1 year, 0.5 if <3 years, 0.2 if older |
| `cross_citation` | 0-1 | 1.0 if 3+ other trusted sources report same finding, 0.5 if 1-2, 0.0 if unique |

**Score interpretation**:
- 70-100: High confidence — integrate as supporting evidence
- 40-69: Moderate — cite as contextual information with caveats
- 20-39: Low — use for narrative framing only, not factual claims
- 0-19: Exclude from manuscript

### Step 5: Integration Decision

| Tier | Verified ratio | Score range | Integration level |
|------|---------------|-------------|-------------------|
| Tier 1 | >80% | 70+ | Cite as supporting evidence alongside PMIDs |
| Tier 1 | 50-80% | 40-69 | Cite as context: "[News: Source, Date]" |
| Tier 2 | >60% | 50+ | Cite as context with verification note |
| Tier 2 | <60% | <50 | Narrative framing only (no specific claims) |
| Tier 3 | Any | Any | Expert opinion attribution: "According to [expert]..." |
| Any | <20% | <20 | Exclude |

---

## 4. Citation Format

News sources use a distinct citation format to differentiate from peer-reviewed references:

```
[News: Author, "Title", Publication, Date. Verified: PMID:XXXXX]
[News: Author, "Title", Publication, Date. Unverified — no primary source found]
[Commentary: Author, "Title", Blog/Publication, Date. Expert opinion]
```

---

## 5. Example Processing

### Example 1 — Tier 1: WHO Cancer Fact Sheet

**Source**: WHO, "Cancer" fact sheet, who.int/news-room/fact-sheets/detail/cancer
**Date**: Updated 2024
**Tier**: 1 (International health authority)

**Step 1 — Source verification**: WHO is Tier 1. Institutional authorship (WHO editorial team). Editorial process: WHO fact sheets undergo multi-department review. Date verified. Publicly accessible. **PASS.**

**Step 2 — Claim extraction** (21 factual claims identified):
- FACTUAL: "Cancer accounted for nearly 10 million deaths in 2022"
- FACTUAL: "Lung cancer: 2.5 million new cases [in 2022]"
- FACTUAL: "Approximately 38% of cancers can currently be prevented"
- FACTUAL: "Approximately 10% of cancers diagnosed in 2022 were attributed to carcinogenic infections"
- INTERPRETIVE: "Cancer is a leading cause of death worldwide"
- (16 additional factual claims with specific statistics)

**Step 3 — Cross-reference**:
- Claim "10 million deaths in 2022" → Source: IARC Global Cancer Observatory (Ferlay et al., gco.iarc.fr). **VERIFIED (self-referencing to IARC data).**
- Claim "38% preventable" → Source: Fink et al., Nature Medicine 2026. **VERIFIED (named publication).**
- Claim "10% attributed to infections" → Source: IARC. **VERIFIED (self-referencing).**
- Verified ratio: 21/21 = 100% (all claims traceable to named sources)

**Step 4 — Credibility score**:
- tier_weight: 1.0 (Tier 1)
- verified_ratio: 1.0 (21/21)
- author_score: 1.0 (institutional, editorial process)
- recency: 0.8 (updated 2024, <3 years old)
- cross_citation: 1.0 (WHO data widely cited by NCI, IARC, national registries)
- **Score: 1.0 × (40×1.0 + 30×1.0 + 20×0.8 + 10×1.0) = 96/100**

**Step 5 — Integration decision**: Score 96, Tier 1, verified ratio 100%. **INTEGRATE as supporting evidence.** Citation: `[News: WHO, "Cancer Fact Sheet", who.int, 2024. Verified: IARC GCO, Fink et al. Nat Med 2026]`

---

### Example 2 — Tier 2: STAT News on Grail Cancer Test

**Source**: Matthew Herper & Angus Chen, "Key study of Grail's cancer detection test fails", STAT News, Feb 19, 2026
**URL**: statnews.com/2026/02/19/grail-cancer-test-galleri-results/
**Tier**: 2 (Science journalism)

**Step 1 — Source verification**: STAT News is Tier 2. Authors: Matthew Herper (Senior Writer, Medicine) and Angus Chen (Cancer Reporter) — both named with credentials. STAT has a corrections policy and editorial board. Date verified. Partially paywalled. **PASS.**

**Step 2 — Claim extraction**:
- FACTUAL: "Grail sold 185,000 tests in 2025, generating $136.8 million in revenue"
- FACTUAL: "List price for Galleri: $1,000"
- FACTUAL: "Stock price declined 47% in after-hours trading"
- FACTUAL: "Test is not yet FDA-approved despite being commercially available"
- FACTUAL: "[The test] failed to meet its main goal in a giant study being conducted with England's NHS"
- INTERPRETIVE: Implied significance of the failure for the liquid biopsy field

**Step 3 — Cross-reference**:
- Grail/Galleri NHS study → ClinicalTrials.gov NCT identifier exists (NHS-Galleri trial). **VERIFIED (trial registry).**
- Revenue figures → Grail SEC filings / earnings reports. **VERIFIED (public financial data).**
- FDA approval status → fda.gov search shows no Galleri approval. **VERIFIED.**
- Stock price decline → Financial data publicly available. **VERIFIED.**
- Verified ratio: 5/5 = 100%

**Step 4 — Credibility score**:
- tier_weight: 0.8 (Tier 2)
- verified_ratio: 1.0
- author_score: 1.0 (named, credentialed)
- recency: 1.0 (<6 months old)
- cross_citation: 1.0 (widely reported: Reuters, Endpoints, FiercePharma)
- **Score: 0.8 × (40×1.0 + 30×1.0 + 20×1.0 + 10×1.0) = 80/100**

**Step 5 — Integration decision**: Score 80, Tier 2, verified ratio 100%. **CITE as context.** Citation: `[News: Herper & Chen, "Key study of Grail's cancer detection test fails", STAT News, Feb 2026. Verified: NHS-Galleri trial registry, Grail SEC filings]`

---

### Example 3 — Tier 3: ACS Research News on Colorectal Cancer

**Source**: ACS Research News, "Colorectal Cancer Drops in Older Adults and Rises in Younger Ones", cancer.org
**URL**: cancer.org/research/acs-research-news/colorectal-cancer-drops-in-older-adults-and-rises-in-young-ones.html
**Tier**: 3 (Advocacy organization blog)

**Step 1 — Source verification**: ACS Research News is Tier 3 (advocacy + research organization). Authors: ACS research team (Siegel, Wagle, Star, Kratzer, Smith, Jemal — named researchers). Editorial process: ACS internal review. Date verified. Publicly accessible. **PASS.**

**Step 2 — Claim extraction**:
- FACTUAL: "158,850 new CRC cases will be diagnosed in 2026, and 55,230 will die"
- FACTUAL: "Under 50: 3% per year increase; 65+: 2.5% per year decrease"
- FACTUAL: "27% of under-50 patients have distant stage CRC"
- FACTUAL: "Rectal cancer now represents 32% of all cases vs 27% in mid-2000s"
- INTERPRETIVE: "CRC can no longer be called an old person's disease" (Ahmedin Jemal)
- FACTUAL: Published as "Colorectal Cancer Statistics, 2026" in CA: A Cancer Journal for Clinicians

**Step 3 — Cross-reference**:
- All statistics → Siegel et al., "Colorectal Cancer Statistics, 2026", CA Cancer J Clin. **VERIFIED (named peer-reviewed publication).**
- The blog post IS a summary of the peer-reviewed paper by the same authors.
- Verified ratio: 5/5 factual claims = 100% (all traceable to the CA paper)

**Step 4 — Credibility score**:
- tier_weight: 0.6 (Tier 3)
- verified_ratio: 1.0
- author_score: 1.0 (named researchers, MD/PhD credentials)
- recency: 1.0 (<6 months)
- cross_citation: 0.5 (original paper widely cited; blog post itself less so)
- **Score: 0.6 × (40×1.0 + 30×1.0 + 20×1.0 + 10×0.5) = 57/100**

**Step 5 — Integration decision**: Score 57, Tier 3, verified ratio 100%. **CITE as context**, but prefer the underlying CA Cancer J Clin paper for factual claims. Citation: `[Commentary: Siegel et al./ACS, "Colorectal Cancer Drops in Older Adults...", ACS Research News, 2026. Primary source: Siegel et al., CA Cancer J Clin 2026]`

---

## 6. Framework Limitations

1. **The credibility score is a starting point**, not a gold standard. Weights (40/30/20/10) are reasonable defaults but not empirically validated. Adjust based on experience.

2. **Paywalled sources** (STAT News, The Cancer Letter) may have higher-quality content that's harder to verify. Note accessibility status but don't penalize for paywalls.

3. **Self-referencing sources** (WHO citing IARC, NCI citing its own data) receive "verified" status because the institution IS the authority. This is different from independent verification.

4. **Expert commentary** (Tier 3) can contain valuable insights that no peer-reviewed paper captures. The framework intentionally preserves this voice while preventing it from being cited as factual evidence.

5. **The framework evaluates individual articles**, not sources as a whole. A Tier 2 source can publish both excellent and poor articles. Each article is scored independently.

6. **This framework is for the cancer research manuscript specifically**. It does not claim to be a general-purpose news credibility framework.
