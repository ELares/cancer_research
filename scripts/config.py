"""Shared configuration for all fetch/enrichment scripts."""

import os
import time
from pathlib import Path

import requests
from dotenv import load_dotenv

# Load .env from project root
PROJECT_ROOT = Path(__file__).resolve().parent.parent
load_dotenv(PROJECT_ROOT / ".env")

# --- Paths ---
CORPUS_DIR = PROJECT_ROOT / "corpus"
PMID_DIR = CORPUS_DIR / "by-pmid"
ABSTRACT_PMID_DIR = CORPUS_DIR / "abstracts" / "by-pmid"
DOI_LOOKUP = CORPUS_DIR / "by-doi" / "DOI_LOOKUP.jsonl"
INDEX_FILE = CORPUS_DIR / "INDEX.jsonl"
TAGS_DIR = PROJECT_ROOT / "tags"

# --- API Keys ---
NCBI_API_KEY = os.getenv("NCBI_API_KEY", "")
OPENALEX_EMAIL = os.getenv("OPENALEX_EMAIL", "")
OPENALEX_API_KEY = os.getenv("OPENALEX_API_KEY", "")
SEMANTIC_SCHOLAR_API_KEY = os.getenv("SEMANTIC_SCHOLAR_API_KEY", "")
CORE_API_KEY = os.getenv("CORE_API_KEY", "")

# --- Rate Limiters ---

class RateLimiter:
    """Simple token-bucket rate limiter."""

    def __init__(self, requests_per_second: float):
        self.min_interval = 1.0 / requests_per_second
        self.last_request = 0.0

    def wait(self):
        now = time.time()
        elapsed = now - self.last_request
        if elapsed < self.min_interval:
            time.sleep(self.min_interval - elapsed)
        self.last_request = time.time()


# Rate limiters per API
NCBI_RATE = RateLimiter(9 if NCBI_API_KEY else 2.5)  # slightly under limits
OPENALEX_RATE = RateLimiter(9)
PMC_BIOC_RATE = RateLimiter(2)
PUBTATOR_RATE = RateLimiter(2)
ICITE_RATE = RateLimiter(5)


def resilient_get(url: str, params: dict = None, timeout: int = 30, retries: int = 2, rate_limiter: 'RateLimiter | None' = None) -> requests.Response:
    """GET with automatic retry on transient failures (5xx, timeout, connection error)."""
    last_exc = None
    for attempt in range(1 + retries):
        if rate_limiter:
            rate_limiter.wait()
        try:
            resp = requests.get(url, params=params, timeout=timeout)
            if resp.status_code < 500:
                return resp
            # 5xx — retry after backoff
            last_exc = requests.HTTPError(f"HTTP {resp.status_code}")
        except (requests.ConnectionError, requests.Timeout) as e:
            last_exc = e

        if attempt < retries:
            time.sleep(2 ** attempt)  # 1s, 2s backoff

    raise last_exc  # type: ignore


# --- API Base URLs ---
PUBMED_ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
PUBMED_EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
PMC_BIOC = "https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_json"
PUBTATOR_API = "https://www.ncbi.nlm.nih.gov/research/pubtator3-api/publications/export/biocjson"
ICITE_API = "https://icite.od.nih.gov/api/pubs"
OPENALEX_WORKS = "https://api.openalex.org/works"
CLINICALTRIALS_API = "https://clinicaltrials.gov/api/v2/studies"
EUROPEPMC_API = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"

# --- Mechanism & Cancer Type Keywords (for auto-tagging) ---

MECHANISM_KEYWORDS = {
    "ttfields": [
        "tumor treating fields", "tumour treating fields", "ttfields", "optune",
        "alternating electric fields", "tumor-treating fields",
    ],
    "immunotherapy": [
        "immunotherapy", "immune checkpoint", "checkpoint inhibitor",
        "anti-pd-1", "anti-pd-l1", "anti-ctla-4", "pembrolizumab",
        "nivolumab", "ipilimumab", "atezolizumab", "durvalumab",
    ],
    "car-t": [
        "car-t", "car t", "chimeric antigen receptor", "cart cell",
    ],
    "crispr": [
        "crispr", "cas9", "cas12", "cas13", "gene editing", "genome editing",
    ],
    "nanoparticle": [
        "nanoparticle", "nanocarrier", "nanomedicine", "liposome",
        "lipid nanoparticle", "nano-delivery", "quantum dot",
    ],
    "metabolic-targeting": [
        "warburg effect", "metabolic reprogramming", "glycolysis inhibit",
        "cancer metabolism", "metabolic targeting", "glutamine deprivation",
        "oxidative phosphorylation", "metabolic vulnerability",
    ],
    "oncolytic-virus": [
        "oncolytic virus", "oncolytic virotherapy", "viral oncolysis",
        "t-vec", "talimogene", "oncolytic herpes", "oncolytic adenovirus",
    ],
    "mRNA-vaccine": [
        "mrna vaccine", "mrna cancer vaccine", "mrna-based", "messenger rna vaccine",
        "personalized neoantigen", "neoantigen vaccine",
    ],
    "synthetic-lethality": [
        "synthetic lethality", "synthetic lethal", "parp inhibitor",
        "olaparib", "niraparib", "rucaparib", "talazoparib", "brca synthetic",
    ],
    "bioelectric": [
        "bioelectricity", "bioelectric signaling", "bioelectric signalling",
        "membrane potential", "transmembrane potential", "vmem",
        "depolarization cancer", "ion channel cancer",
    ],
    "electrolysis": [
        "electrolysis", "electrochemical therapy", "electrochemical treatment",
        "electrolytic ablation", "echt", "galvanotherapy",
    ],
    "sonodynamic": [
        "sonodynamic", "sonosensitizer", "ultrasound therapy cancer",
        "acoustic therapy cancer",
    ],
    "cold-atmospheric-plasma": [
        "cold atmospheric plasma", "cap therapy", "plasma jet cancer",
        "non-thermal plasma cancer", "atmospheric pressure plasma",
    ],
    "hifu": [
        "hifu", "high intensity focused ultrasound", "focused ultrasound ablation",
        "mrgfus", "mr-guided focused ultrasound",
    ],
    "electrochemical-therapy": [
        "irreversible electroporation", "nanoknife",
        "pulsed electric field", "electroporation cancer",
        "electroporation therapy", "electroporation ablation",
    ],
    "epigenetic": [
        "epigenetic therapy", "dna methylation cancer", "hdac inhibitor",
        "histone deacetylase", "epigenetic reprogramming", "dnmt inhibitor",
        "azacitidine", "decitabine", "vorinostat", "romidepsin",
    ],
    "microbiome": [
        "microbiome cancer", "gut microbiota cancer", "fecal microbiota",
        "microbiome immunotherapy", "intratumoral bacteria",
    ],
    "frequency-therapy": [
        "radiofrequency ablation", "rfa cancer", "microwave ablation",
        "electromagnetic frequency", "resonant frequency cancer",
        "pemf cancer", "pulsed electromagnetic",
    ],
    "antibody-drug-conjugate": [
        "antibody-drug conjugate", "antibody drug conjugate",
        "trastuzumab deruxtecan", "enhertu", "sacituzumab", "brentuximab",
        "ado-trastuzumab", "trastuzumab emtansine", "polatuzumab",
    ],
    "bispecific-antibody": [
        "bispecific antibody", "bispecific t-cell", "bispecific engager",
        "blinatumomab", "teclistamab", "mosunetuzumab", "glofitamab",
    ],
    "radioligand-therapy": [
        "radioligand therapy", "radiopharmaceutical therapy", "radionuclide therapy",
        "targeted radionuclide therapy", "targeted radioligand therapy",
        "peptide receptor radionuclide therapy", "prrt", "psma radioligand",
        "lutetium-177", "actinium-225", "radium-223", "lutathera", "pluvicto",
        "177lu-dotatate", "177lu psma", "225ac psma", "radioiodine therapy",
    ],
    "targeted-protein-degradation": [
        "protac", "proteolysis targeting chimera", "molecular glue",
        "targeted protein degradation", "degrader",
    ],
    "phagocytosis-checkpoint": [
        "cd47", "sirpa", "sirpalpha", "phagocytosis checkpoint",
        "don't eat me signal", "dont eat me signal",
    ],
}

BIOLOGY_PROCESS_KEYWORDS = {
    "autophagy": [
        "autophagy", "autophagic", "lysosomal degradation", "autophagosome",
    ],
    "senescence-sasp": [
        "senescence", "senescent", "sasp", "senolytic",
        "therapy-induced senescence",
    ],
    "tme-stroma": [
        "tumor microenvironment", "tumour microenvironment", "cancer-associated fibroblast",
        "caf", "cafs", "extracellular matrix", "stromal barrier", "ecm remodeling",
    ],
    "cuproptosis": [
        "cuproptosis", "copper ionophore", "elesclomol", "fdx1",
    ],
    "disulfidptosis": [
        "disulfidptosis", "glucose starvation-induced disulfide stress",
    ],
}

PATHWAY_TARGET_KEYWORDS = {
    "dhodh-defense": [
        "dhodh", "dihydroorotate dehydrogenase", "brequinar",
    ],
    "dhcr7-7dhc-axis": [
        "dhcr7", "7-dhc", "7-dehydrocholesterol", "7 dehydrocholesterol",
    ],
    "mboat1-mboat2-axis": [
        "mboat1", "mboat2",
    ],
    "scd-mufa-axis": [
        "scd", "scd1", "stearoyl-coa desaturase", "stearoyl coa desaturase",
        "monounsaturated fatty acid", "mufa enrichment",
    ],
    "fdx1-cuproptosis-axis": [
        "fdx1", "ferredoxin 1", "cuproptosis", "elesclomol", "copper ionophore",
    ],
    "trim25-gpx4-degradation": [
        "trim25", "n6f11", "gpx4 degradation", "selective gpx4 degradation",
    ],
    "cuproptosis-core": [
        "cuproptosis", "copper ionophore", "elesclomol", "fdx1", "lipoylated tca",
    ],
    "disulfidptosis-core": [
        "disulfidptosis", "disulfide stress", "glucose starvation-induced disulfide stress",
        "slc7a11", "xct-dependent disulfide stress",
    ],
}

RADIOLIGAND_TARGET_KEYWORDS = {
    "psma": [
        "psma", "prostate-specific membrane antigen", "prostate specific membrane antigen",
        "vipivotide tetraxetan",
    ],
    "fap": [
        "fap", "fibroblast activation protein", "fibroblast activation protein alpha", "fapi",
    ],
    "sstr": [
        "sstr", "somatostatin receptor", "somatostatin receptor 2", "sst2", "dotatate", "dotatoc",
    ],
    "cea": [
        "cea", "ceacam5", "carcinoembryonic antigen",
    ],
}

RESISTANT_STATE_RULES = {
    "oxphos-dependent-persister": {
        "all_of": [
            ["drug-tolerant persister", "drug tolerant persister", "persister cell", "residual disease"],
            ["oxidative phosphorylation", "oxphos", "mitochondrial respiration"],
        ]
    },
    "nrf2-compensated-ferroptosis-resistant": {
        "all_of": [
            ["nrf2", "antioxidant response"],
            ["ferroptosis resistance", "gpx4 compensation", "gsh homeostasis", "redox compensation"],
        ]
    },
    "slc7a11-high-disulfidptosis-prone": {
        "all_of": [
            ["slc7a11", "xct", "cystine transporter"],
            ["disulfidptosis", "glucose starvation", "disulfide stress"],
        ]
    },
    "therapy-induced-senescence": {
        "all_of": [
            ["therapy-induced senescence", "treatment-induced senescence", "senescent tumor cell"],
            ["drug resistance", "persister", "adaptive resistance", "residual disease"],
        ]
    },
    "stromal-sheltered-immune-excluded": {
        "all_of": [
            ["cancer-associated fibroblast", "tumor microenvironment", "extracellular matrix", "stromal barrier"],
            ["immune exclusion", "immune desert", "t cell exclusion", "stromal-mediated resistance", "stromal shelter"],
        ]
    },
    "epigenetically-plastic": {
        "all_of": [
            ["epigenetic plasticity", "chromatin state", "kdm5", "ezh2", "hdac inhibitor", "dedifferentiation"],
            ["persister", "drug tolerance", "adaptive resistance", "cell state transition"],
        ]
    },
}

CANCER_TYPE_KEYWORDS = {
    "breast": ["breast cancer", "breast neoplasm", "breast carcinoma", "triple-negative breast", "tnbc", "her2-positive breast", "breast tumor", "breast tumour", "mammary carcinoma"],
    "lung": ["lung cancer", "lung neoplasm", "nsclc", "non-small cell lung", "small cell lung", "sclc", "lung adenocarcinoma"],
    "colorectal": ["colorectal cancer", "colon cancer", "rectal cancer", "colorectal neoplasm", "crc"],
    "glioblastoma": ["glioblastoma", "gbm", "glioma", "brain tumor", "brain tumour", "brain cancer"],
    "pancreatic": ["pancreatic cancer", "pancreatic neoplasm", "pancreatic ductal", "pdac"],
    "melanoma": ["melanoma", "cutaneous melanoma", "uveal melanoma", "skin cancer melanoma"],
    "leukemia": ["leukemia", "leukaemia", "acute myeloid leukemia", "chronic myeloid leukemia", "acute lymphoblastic leukemia", "chronic lymphocytic leukemia", "myeloid leukemia", "lymphocytic leukemia", "lymphoblastic leukemia"],
    "lymphoma": ["lymphoma", "hodgkin", "non-hodgkin", "dlbcl", "follicular lymphoma"],
    "prostate": ["prostate cancer", "prostate neoplasm", "prostate carcinoma", "castration-resistant"],
    "ovarian": ["ovarian cancer", "ovarian neoplasm", "ovarian carcinoma", "epithelial ovarian"],
    "liver": ["hepatocellular carcinoma", "hcc", "liver cancer", "hepatic cancer"],
    "gastric": ["gastric cancer", "stomach cancer", "gastric neoplasm", "gastric carcinoma"],
    "cervical": ["cervical cancer", "cervical neoplasm", "cervical carcinoma"],
    "bladder": ["bladder cancer", "urothelial carcinoma", "bladder neoplasm"],
    "kidney": ["renal cell carcinoma", "kidney cancer", "renal cancer", "rcc"],
    "thyroid": ["thyroid cancer", "thyroid carcinoma", "papillary thyroid", "thyroid neoplasm"],
    "esophageal": ["esophageal cancer", "oesophageal cancer", "esophageal carcinoma"],
    "head-and-neck": ["head and neck cancer", "oral cancer", "oropharyngeal", "nasopharyngeal", "hnscc"],
    "sarcoma": [
        "sarcoma", "osteosarcoma", "osteogenic sarcoma", "soft tissue sarcoma",
        "soft-tissue sarcoma", "ewing sarcoma", "ewing's sarcoma",
        "rhabdomyosarcoma", "synovial sarcoma",
    ],
    "myeloma": ["multiple myeloma", "myeloma", "plasma cell myeloma"],
    "mesothelioma": ["mesothelioma", "pleural mesothelioma"],
    "neuroblastoma": ["neuroblastoma"],
}

CANCER_SUBTYPE_KEYWORDS = {
    "osteosarcoma": [
        "osteosarcoma", "osteogenic sarcoma",
    ],
    "ewing-sarcoma": [
        "ewing sarcoma", "ewing's sarcoma", "ewing family tumor", "ewing family tumour",
    ],
    "rhabdomyosarcoma": [
        "rhabdomyosarcoma", "embryonal rhabdomyosarcoma", "alveolar rhabdomyosarcoma",
    ],
    "synovial-sarcoma": [
        "synovial sarcoma",
    ],
    "soft-tissue-sarcoma": [
        "soft tissue sarcoma", "soft-tissue sarcoma",
    ],
}

CANCER_SUBTYPE_ORDER = [
    "osteosarcoma",
    "ewing-sarcoma",
    "rhabdomyosarcoma",
    "synovial-sarcoma",
    "soft-tissue-sarcoma",
]

TISSUE_CATEGORY_ORDER = [
    "epithelial",
    "hematologic",
    "mesenchymal",
    "neuroectodermal",
    "mesothelial",
]

CANCER_TYPE_TO_TISSUE = {
    "breast": "epithelial",
    "lung": "epithelial",
    "colorectal": "epithelial",
    "pancreatic": "epithelial",
    "melanoma": "neuroectodermal",
    "leukemia": "hematologic",
    "lymphoma": "hematologic",
    "prostate": "epithelial",
    "ovarian": "epithelial",
    "liver": "epithelial",
    "gastric": "epithelial",
    "cervical": "epithelial",
    "bladder": "epithelial",
    "kidney": "epithelial",
    "thyroid": "epithelial",
    "esophageal": "epithelial",
    "head-and-neck": "epithelial",
    "sarcoma": "mesenchymal",
    "myeloma": "hematologic",
    "mesothelioma": "mesothelial",
    "glioblastoma": "neuroectodermal",
    "neuroblastoma": "neuroectodermal",
}


def derive_tissue_categories(cancer_types: list[str]) -> list[str]:
    derived = {CANCER_TYPE_TO_TISSUE[c] for c in cancer_types if c in CANCER_TYPE_TO_TISSUE}
    return [t for t in TISSUE_CATEGORY_ORDER if t in derived]


def derive_sarcoma_subtypes(
    matched_subtypes: list[str],
    cancer_types: list[str],
    title_subtypes: list[str] | None = None,
    abstract_subtypes: list[str] | None = None,
) -> list[str]:
    """Only surface sarcoma-family subtypes when the paper looks subtype-focused.

    Title mentions are treated as strong focus signals. Abstract-only mentions are
    accepted only when the article is otherwise sarcoma-focused, which avoids
    promoting broad multi-cancer comparison papers into subtype counts.
    """
    if "sarcoma" not in cancer_types:
        return []
    title_set = set(title_subtypes or [])
    abstract_set = set(abstract_subtypes or [])
    subtype_set = set(matched_subtypes)
    if title_set:
        subtype_set &= title_set
    elif len(cancer_types) == 1:
        subtype_set &= abstract_set
    else:
        return []
    return [subtype for subtype in CANCER_SUBTYPE_ORDER if subtype in subtype_set]

# ---------------------------------------------------------------------------
# Diagnostic-to-Therapy Matching (pilot layer — issue #41)
# ---------------------------------------------------------------------------
# Each chain maps a diagnostic modality through a targetable feature to an
# intervention class.  Matching requires the intervention link PLUS at least
# one of (diagnostic, feature) to reduce false positives from papers that
# discuss only a diagnostic or only a therapy in passing.

DIAGNOSTIC_THERAPY_KEYWORDS = {
    "psma-imaging-to-radioligand": {
        "diagnostic": [
            "psma pet", "psma imaging", "psma scan", "68ga-psma", "psma-11",
            "psma pet/ct", "psma-pet",
        ],
        "feature": [
            "psma expression", "psma-positive", "psma positive",
            "prostate-specific membrane antigen",
        ],
        "intervention": [
            "177lu-psma", "lu-psma", "psma radioligand", "psma-617",
            "vipivotide", "pluvicto", "lutetium-psma",
        ],
    },
    "sstr-imaging-to-prrt": {
        "diagnostic": [
            "sstr scintigraphy", "dotatate pet", "dotatoc pet",
            "68ga-dotatate", "68ga-dotatoc",
            "somatostatin receptor imaging", "sstr pet",
        ],
        "feature": [
            "sstr expression", "sstr2-positive", "sstr2 positive",
            "somatostatin receptor positive",
        ],
        "intervention": [
            "lutathera", "177lu-dotatate", "177lu-dotatoc",
            "peptide receptor radionuclide therapy", "prrt",
        ],
    },
    "pdl1-ihc-to-checkpoint": {
        "diagnostic": [
            "pd-l1 immunohistochemistry", "pd-l1 ihc", "pd-l1 staining",
            "tps score", "tumor proportion score",
            "cps score", "combined positive score",
        ],
        "feature": [
            "pd-l1 positive", "pd-l1 high", "pd-l1 expression",
            "pd-l1-positive",
        ],
        "intervention": [
            "pembrolizumab", "nivolumab", "atezolizumab", "durvalumab",
            "avelumab", "cemiplimab",
        ],
    },
    "tmb-msi-to-immunotherapy": {
        "diagnostic": [
            "tumor mutational burden", "tmb-high", "tmb-h",
            "microsatellite instability", "msi-high", "msi-h",
            "mismatch repair deficient", "dmmr",
        ],
        "feature": [
            "tmb-high", "tmb-h", "msi-h", "msi-high",
            "hypermutated", "mismatch repair deficient",
        ],
        "intervention": [
            "pembrolizumab", "nivolumab", "checkpoint inhibitor",
            "immune checkpoint", "anti-pd-1", "anti-pd1",
        ],
    },
    "neoantigen-profiling-to-mrna-vaccine": {
        "diagnostic": [
            "neoantigen prediction", "neoantigen profiling",
            "neoantigen identification", "neoantigen discovery",
            "whole exome sequencing", "tumor sequencing",
            "mutanome", "immunopeptidome",
        ],
        "feature": [
            "neoantigen", "neo-antigen", "tumor-specific antigen",
            "personalized antigen", "individualized neoantigen",
        ],
        "intervention": [
            "mrna vaccine", "mrna cancer vaccine", "personalized vaccine",
            "individualized mrna", "neoantigen vaccine",
            "autogene cevumeran", "mrna-4157",
        ],
    },
    "oncolytic-susceptibility-to-virotherapy": {
        "diagnostic": [
            "viral receptor expression", "nectin-1 expression",
            "cd46 expression", "coxsackievirus receptor",
            "herpes simplex entry", "oncolytic susceptibility",
        ],
        "feature": [
            "viral entry receptor", "nectin-1", "cd46",
            "interferon deficiency", "interferon-deficient",
        ],
        "intervention": [
            "t-vec", "talimogene", "oncolytic herpes",
            "oncolytic adenovirus", "oncolytic virus therapy",
            "oncolytic vaccinia", "oncolytic reovirus",
        ],
    },
}

DIAGNOSTIC_THERAPY_ORDER = [
    "psma-imaging-to-radioligand",
    "sstr-imaging-to-prrt",
    "pdl1-ihc-to-checkpoint",
    "tmb-msi-to-immunotherapy",
    "neoantigen-profiling-to-mrna-vaccine",
    "oncolytic-susceptibility-to-virotherapy",
]

EVIDENCE_LEVEL_KEYWORDS = {
    "phase3-clinical": [
        "phase 3", "phase iii", "phase-3",
        "pivotal trial", "randomized controlled trial", "randomized clinical trial",
    ],
    "phase2-clinical": ["phase 2", "phase ii", "phase-2"],
    "phase1-clinical": ["phase 1", "phase i ", "phase-1", "first-in-human", "dose-escalation"],
    "clinical-other": [
        "pilot study", "pilot trial", "feasibility study", "feasibility trial",
        "single-arm", "single arm", "retrospective study", "retrospective analysis",
        "retrospective cohort", "retrospective review",
        "case report", "case series", "investigator-initiated", "investigator initiated",
        "real-world study", "real world study", "real-world analysis", "real world analysis",
        "real-world cohort", "real world cohort", "registry study", "registry analysis",
        "clinical experience",
        "reported a case", "single patient",
    ],
    "preclinical-invivo": ["in vivo", "mouse model", "xenograft", "animal model", "murine", "tumor-bearing mice"],
    "preclinical-invitro": ["in vitro", "cell line", "cell culture", "cultured cells"],
    "theoretical": ["computational model", "mathematical model", "simulation", "theoretical framework", "in silico"],
}

# --- News Source Tiers ---
# See analysis/news-source-criteria.md for full framework documentation.
#
# This is a first-pass subset of the 31 sources listed in the criteria doc.
# Sources behind paywalls or those requiring manual URL discovery (e.g.
# university press offices) are omitted here and must be matched manually.
# classify_source() should match the longest path prefix, not the bare domain.
#
# IMPORTANT: Only editorial/news surfaces belong here.  Primary research
# paths (e.g. nature.com/articles, science.org/doi, nejm.org/doi) are
# PubMed-indexed journal content and must NOT be classified as news-tier
# sources.  If a URL matches a journal's primary-research path, it belongs
# in the corpus layer, not the news layer.

SOURCE_TIER_DEFINITIONS = {
    "tier1": {
        "name": "Institutional and peer-adjacent",
        "trust": "high",
        "cite_as": "evidence",
        "verified_threshold": 0.80,
        "path_prefixes": [
            "cancer.gov/news-events",
            "cancer.gov/about-nci",
            "nih.gov/news-events",
            "who.int/news-room",
            "fda.gov/drugs",
            "fda.gov/news-events",
            "clinicaltrials.gov",
            "nature.com/news",
            "science.org/news",
            "cell.com/news",
            "thelancet.com/news",
            "jwatch.org",
            "gco.iarc.fr",
        ],
    },
    "tier2": {
        "name": "Science journalism",
        "trust": "medium",
        "cite_as": "context",
        "verified_threshold": 0.60,
        # Tier 2 outlets are editorially independent; full-domain matching is
        # acceptable because these outlets don't host user-generated content.
        "path_prefixes": [
            "statnews.com",
            "cancerletter.com",
            "endpointsnews.com",
            "reuters.com/business/healthcare-pharmaceuticals",
            "apnews.com/health",
            "sciencedaily.com",
            "medicalnewstoday.com",
            "arstechnica.com/science",
            "theconversation.com",
            "fiercepharma.com",
            "fiercebiotech.com",
        ],
    },
    "tier3": {
        "name": "Expert blogs and commentary",
        "trust": "context-only",
        "cite_as": "opinion",
        "verified_threshold": 0.0,
        "path_prefixes": [
            "connection.asco.org",
            "science.org/blogs/pipeline",
            "vinayprasad.com",
            "cancer.org/research/acs-research-news",
            "lls.org/news",
            "broadinstitute.org/blog",
            "mdanderson.org/cancerwise",
            "icr.ac.uk/blogs",
            "cancerresearchuk.org/about-cancer",
            "patientpower.info",
        ],
    },
}

# --- News Pipeline ---

NEWS_RATE = RateLimiter(2)  # conservative: 2 req/s for news sites

NEWS_DIR = PROJECT_ROOT / "news"

# Regex patterns that indicate a sentence contains a verifiable factual claim.
CLAIM_FACTUAL_MARKERS = [
    r'\d+\.?\d*\s*%',                          # percentages
    r'\$[\d,.]+\s*(?:million|billion|M|B)?',    # dollar amounts
    r'[Pp]hase\s+[I1-3]{1,3}\b',               # trial phases
    r'FDA\s+approv',                            # FDA actions
    r'[Ee]nrolled\s+[\d,]+',                    # enrollment numbers
    r'[\d,]+\s+patients',                       # patient counts
    r'overall\s+survival',                      # clinical endpoints
    r'progression.free\s+survival',             # clinical endpoints
    r'response\s+rate',                         # clinical endpoints
    r'hazard\s+ratio',                          # statistical measures
    r'p\s*[<=]\s*0\.\d+',                       # p-values
    r'median\s+(?:survival|OS|PFS)',            # survival endpoints
    r'five.year\s+survival|5.year\s+survival',  # survival rates
]

# Keywords that suggest a claim's type (checked in priority order).
CLAIM_TYPE_MARKERS = {
    'event':       ['approved', 'announced', 'launched', 'granted', 'designated',
                    'authorized', 'cleared', 'recalled', 'withdrew', 'submitted'],
    'result':      ['showed', 'demonstrated', 'found', 'observed', 'measured',
                    'produced', 'achieved', 'reported a', 'yielded', 'detected'],
    'mechanism':   ['through', 'via', 'pathway', 'mediated', 'mechanism',
                    'inhibit', 'activat', 'regulat', 'modulat', 'target'],
    'opinion':     ['believes', 'argues', 'suggests that', 'according to',
                    'noted that', 'said', 'commented', 'emphasized'],
    'speculation': ['could', 'might', 'may lead', 'potential', 'if confirmed',
                    'promising', 'expected to', 'likely to', 'remains to be seen'],
}
