#!/usr/bin/env python3
"""Convert v1.md to v1.tex with proper LaTeX formatting."""
import re
import unicodedata
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
MD = ROOT / "article" / "drafts" / "v1.md"
TEX = ROOT / "article" / "drafts" / "v1.tex"

md = MD.read_text()
title = re.search(r'^# (.+)$', md, re.MULTILINE).group(1).strip()

# Build cite map from reference list
cite_map = {}
for line in md.split('\n'):
    m = re.match(r'^(\d+)\. PMID: (\d+) -- (\w+)', line.strip())
    if m:
        pmid = m.group(2)
        first = unicodedata.normalize('NFKD', m.group(3)).encode('ascii', 'ignore').decode().lower()
        ym = re.search(r'\((\d{4})\)', line)
        year = ym.group(1) if ym else '2024'
        cite_map[pmid] = f'{first}{year}_{pmid}'
# Fix known parsing issues
cite_map['31130474'] = 'unknown2019_31130474'
cite_map['29978216'] = 'unknown2018_29978216'

# Extract sections
abstract = re.search(r'## Abstract\n\n(.*?)(?=\n\*\*Keywords)', md, re.DOTALL).group(1).strip()
keywords = re.search(r'\*\*Keywords:\*\*\s*(.+)', md).group(1).strip()
body = re.search(r'## 1\. Introduction(.*?)## References', md, re.DOTALL).group(0).replace('## References','').strip()

# Convert PMID citations to placeholders (will become \cite after text conversion)
def repl(m):
    key = cite_map.get(m.group(1), 'pmid'+m.group(1))
    return f'CITEPLACEHOLDER{{{key}}}'
body = re.sub(r'PMID: (\d+)', repl, body)
abstract = re.sub(r'PMID: (\d+)', repl, abstract)

# Markdown → LaTeX
def cvt(t):
    t = re.sub(r'^#### (\d+\.\d+\.\d+) (.+)$', r'\\subsubsection{\2}', t, flags=re.MULTILINE)
    t = re.sub(r'^## (\d+)\. (.+)$', r'\\section{\2}', t, flags=re.MULTILINE)
    t = re.sub(r'^### (\d+\.\d+) (.+)$', r'\\subsection{\2}', t, flags=re.MULTILINE)
    t = re.sub(r'\*\*(.+?)\*\*', r'\\textbf{\1}', t)
    t = re.sub(r'(?<!\*)\*([^*\n]+?)\*(?!\*)', r'\\textit{\1}', t)
    # Escape special chars BEFORE replacing unicode
    t = t.replace('%', '\\%')
    t = t.replace('&', '\\&')
    t = t.replace('#', '\\#')
    # Escape underscores in prose (word_word patterns) but not in LaTeX commands
    # (LaTeX commands use _ for subscripts inside $ delimiters)
    t = re.sub(r'(?<=\w)_(?=\w)', r'\\_', t)
    # Unicode → LaTeX
    t = t.replace('→', '$\\rightarrow$')
    t = t.replace('×', '$\\times$')
    t = t.replace('~', '$\\sim$')
    t = t.replace('↓', '$\\downarrow$')
    t = t.replace('↑', '$\\uparrow$')
    t = t.replace('—', '---')
    t = t.replace('≥', '$\\geq$')
    t = t.replace('≤', '$\\leq$')
    t = t.replace('≈', '$\\approx$')
    t = t.replace('δ', '$\\delta$')
    t = t.replace('α', '$\\alpha$')
    t = t.replace('µ', '$\\mu$')
    t = t.replace('μ', '$\\mu$')      # U+03BC (Greek mu) — distinct from U+00B5 (micro sign)
    t = t.replace('λ', '$\\lambda$')
    t = re.sub(r'√\(([^)]+)\)', r'$\\sqrt{\1}$', t)  # √(x) → $\sqrt{x}$
    t = t.replace('√', '$\\sqrt{}$')                    # bare √ fallback
    t = t.replace('²', '$^2$')
    t = t.replace('₂', '$_2$')
    t = t.replace('−', '$-$')
    t = t.replace('₀', '$_0$')
    # No blanket brace fixes needed — protection handles it
    return t

body_tex = cvt(body)
abstract_tex = cvt(abstract)

# Now convert CITEPLACEHOLDER{key} → \cite{key}
body_tex = re.sub(r'CITEPLACEHOLDER\{([^}]+)\}', r'\\cite{\1}', body_tex)
abstract_tex = re.sub(r'CITEPLACEHOLDER\{([^}]+)\}', r'\\cite{\1}', abstract_tex)
# Remove square brackets around \cite (markdown artifact: [PMID: X] → [\cite{x}])
body_tex = re.sub(r'\[\\cite\{', r'\\cite{', body_tex)
abstract_tex = re.sub(r'\[\\cite\{', r'\\cite{', abstract_tex)
body_tex = re.sub(r'\\cite\{([^}]+)\}\]', r'\\cite{\1}', body_tex)
abstract_tex = re.sub(r'\\cite\{([^}]+)\}\]', r'\\cite{\1}', abstract_tex)

# Remove markdown horizontal rules, which are invalid in LaTeX body text.
body_tex = re.sub(r'^\s*---\s*$', '', body_tex, flags=re.MULTILINE)

# Replace markdown tables with LaTeX tables
# Find pipe-delimited tables and replace
def replace_table(text, marker, caption, label, headers, rows):
    pattern = r'(?m)^' + re.escape(marker) + r'.*\n(?:^\|.*\n)+'
    match = re.search(pattern, text)
    if not match:
        return text

    h = ' & '.join(f'\\textbf{{{h}}}' for h in headers)
    r = ' \\\\\n'.join(' & '.join(cells) for cells in rows)
    cols = 'l' + 'c' * (len(headers)-1)

    table = f"""\\begin{{table}}[ht]
\\centering
\\caption{{{caption}}}
\\label{{{label}}}
\\begin{{tabular}}{{{cols}}}
\\toprule
{h} \\\\
\\midrule
{r} \\\\
\\bottomrule
\\end{{tabular}}
\\end{{table}}"""
    return text[:match.start()] + table + text[match.end():]

# Simulation table
body_tex = replace_table(body_tex, '| Phenotype |',
    'Monte Carlo ferroptosis simulation (n=1M cells/condition).', 'tab:sim',
    ['Phenotype', 'Control', 'RSL3', 'SDT', 'PDT'],
    [['Glycolytic', '0.00\\%', '0.00\\%', '87.2\\%', '87.2\\%'],
     ['OXPHOS', '0.04\\%', '1.1\\%', '99.9\\%', '99.9\\%'],
     ['Persister (FSP1$\\downarrow$)', '1.2\\%', '42.5\\%', '100.0\\%', '100.0\\%'],
     ['Persister + NRF2', '0.00\\%', '0.05\\%', '99.5\\%', '99.5\\%']])

# Modality table
body_tex = replace_table(body_tex, '| Modality |',
    'Ferroptosis engagement across physical modalities (PubMed, March 2026).', 'tab:mod',
    ['Modality', 'Ferroptosis', 'Ferro+ICD', 'Depth'],
    [['\\textbf{PDT}', '\\textbf{355}', '\\textbf{67}', 'Superficial (mm)'],
     ['\\textbf{SDT}', '\\textbf{121}', '\\textbf{25}', 'Deep (cm)'],
     ['IRE', '15', 'emerging', 'Invasive'],
     ['HIFU', '3', 'minimal', 'Deep (cm)'],
     ['TTFields', '0', '0', 'Surface']])

# Clean leftover pipe tables
body_tex = re.sub(r'\|[-|]+\|', '', body_tex)
body_tex = re.sub(r'^\|.*\|$', '', body_tex, flags=re.MULTILINE)

# Replace figure placeholders
figs = {
    '1': ('fig5_publication_trends', 'Publication volume 2015--2025.'),
    '2': ('fig2_mechanism_heatmap', 'Mechanism $\\times$ cancer type matrix.'),
    '3': ('fig1_ferroptosis_comparison', 'Ferroptosis engagement ($\\chi^2=97.3$, $p<5.9\\times10^{-23}$).'),
    '4': ('fig4_molecular_overlap', 'Molecular pathway engagement (normalized \\%).'),
    '5': ('fig3_literature_disconnect', 'Literature disconnect between communities.'),
    '6': ('fig6_sdt_chain_evidence', 'SDT ferroptosis-ICD chain evidence.'),
    '7': ('fig7_monte_carlo_simulation', 'Monte Carlo simulation (1M cells/condition).'),
    '8': ('fig8_simulation_by_treatment', 'Spatial tumor simulation: depth-kill curves and 2D death heatmaps.'),
    '9': ('fig13_gold_set_eval', 'Evidence tagger performance: gold-set evaluation (100-article stratified sample).'),
    '10': ('fig9_evidence_tiers', 'Evidence tier composition by mechanism.'),
    '11': ('fig14_tissue_mechanism_heatmap', 'Tissue-of-origin $\\times$ mechanism article counts (coverage: 62\\%).'),
    '12': ('fig15_designed_combinations', 'Classification of multi-mechanism articles into designed combinations, co-mentions, and reviews.'),
    '13': ('fig16_weighted_evidence', 'Weighted evidence score by mechanism (tier $\\times$ citation percentile $\\times$ recency).'),
}
def repl_figure(match):
    num = match.group(1)
    if num not in figs:
        return match.group(0)
    fn, cap = figs[num]
    has_description = ':' in match.group(0)
    # Standalone placeholders (with description, on own line) → full figure environment
    # Inline references (no description, inside paragraph) → ref only
    if has_description:
        return f"""\\begin{{figure}}[ht]
\\centering
\\includegraphics[width=\\textwidth]{{../figures/{fn}.pdf}}
\\caption{{{cap}}}
\\label{{fig:{fn}}}
\\end{{figure}}"""
    else:
        return f'(Figure~\\ref{{fig:{fn}}})'

body_tex = re.sub(r'\[FIGURE (\d+)(?::[^\]]*)?\]', repl_figure, body_tex)

# Fail loudly if any figure placeholders survived substitution
leftover = re.findall(r'\[FIGURE \d+(?::[^\]]*)?\]', body_tex)
if leftover:
    print("ERROR: Unhandled figure placeholders in manuscript:")
    for placeholder in leftover:
        print(f"  {placeholder[:80]}")
    print("Add missing entries to the `figs` dict in generate_latex.py.")
    raise SystemExit(1)

latex = f"""\\documentclass[12pt,a4paper]{{article}}
\\usepackage[utf8]{{inputenc}}
\\usepackage[T1]{{fontenc}}
\\usepackage{{amsmath,amssymb}}
\\usepackage{{graphicx}}
\\usepackage{{hyperref}}
\\usepackage{{natbib}}
\\usepackage{{booktabs}}
\\usepackage{{geometry}}
\\usepackage{{setspace}}
\\usepackage{{authblk}}

\\geometry{{margin=1in}}
\\onehalfspacing

\\title{{{title}}}

\\author[1]{{Ezequiel Lares}}
\\affil[1]{{Independent Researcher}}

\\date{{}}

\\begin{{document}}

\\maketitle

\\begin{{abstract}}
{abstract_tex}
\\end{{abstract}}

\\textbf{{Keywords:}} {keywords}

{body_tex}

\\bibliographystyle{{unsrtnat}}
\\bibliography{{../references/bibliography}}

\\end{{document}}
"""

TEX.write_text(latex)
print(f'Written {TEX}: {len(latex)} chars, {latex.count(chr(92)+"cite{")} citations')
