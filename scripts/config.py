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
        "radioligand therapy", "radiopharmaceutical therapy", "theranostic",
        "theranostics", "lutetium-177", "pluvicto", "lutathera",
        "psma radioligand", "radionuclide therapy",
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
    "sarcoma": ["sarcoma", "osteosarcoma", "soft tissue sarcoma", "ewing sarcoma"],
    "myeloma": ["multiple myeloma", "myeloma", "plasma cell myeloma"],
    "mesothelioma": ["mesothelioma", "pleural mesothelioma"],
    "neuroblastoma": ["neuroblastoma"],
}

EVIDENCE_LEVEL_KEYWORDS = {
    "phase3-clinical": [
        "phase 3", "phase iii", "phase-3",
        "pivotal trial", "randomized controlled trial", "randomized clinical trial",
    ],
    "phase2-clinical": ["phase 2", "phase ii", "phase-2"],
    "phase1-clinical": ["phase 1", "phase i ", "phase-1", "first-in-human", "dose-escalation"],
    "preclinical-invivo": ["in vivo", "mouse model", "xenograft", "animal model", "murine", "tumor-bearing mice"],
    "preclinical-invitro": ["in vitro", "cell line", "cell culture", "cultured cells"],
    "theoretical": ["computational model", "mathematical model", "simulation", "theoretical framework", "in silico"],
}
