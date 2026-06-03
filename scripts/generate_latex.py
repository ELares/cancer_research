#!/usr/bin/env python3
"""Convert v1.md to v1.tex with proper LaTeX formatting.

Supports book-style structure with Parts, Chapters, Sections, and Subsections.
Document class: report (not book — avoids forced recto chapter starts).
See article/AUTHORING.md for heading conventions.
"""
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
MD = ROOT / "article" / "drafts" / "v1.md"
TEX = ROOT / "article" / "drafts" / "v1.tex"

md = MD.read_text()
title = re.search(r'^# (.+)$', md, re.MULTILINE).group(1).strip()

# Build footnote definition map from Markdown [^label]: text patterns
footnote_defs = {}
for m in re.finditer(r'^\[\^(\w+)\]:\s*(.+)$', md, re.MULTILINE):
    footnote_defs[m.group(1)] = m.group(2).strip()

# Extract sections
abstract = re.search(r'## Abstract\n\n(.*?)(?=\n\*\*Keywords)', md, re.DOTALL).group(1).strip()
keywords = re.search(r'\*\*Keywords:\*\*\s*(.+)', md).group(1).strip()

# Body: everything from the first Part header to just before References.
body_start = re.search(r'^# Part [IVX]+: ', md, re.MULTILINE)
ref_match = re.search(r'^## References', md, re.MULTILINE)
if not body_start or not ref_match:
    raise SystemExit("ERROR: Could not find '# Part ...' or '## References' boundaries in v1.md")
body = md[body_start.start():ref_match.start()].strip()

# Remove footnote definition lines from the body (they'll become \footnote{} inline)
body = re.sub(r'^\[\^\w+\]:\s*.+$', '', body, flags=re.MULTILINE)

# Markdown → LaTeX
def cvt(t):
    # Book-structure headings (report document class)
    t = re.sub(r'^# Part [IVX]+: (.+)$', r'\\part{\1}', t, flags=re.MULTILINE)
    t = re.sub(r'^## Chapter \d+: (.+)$', r'\\chapter{\1}', t, flags=re.MULTILINE)
    t = re.sub(r'^### \d+\.\d+ (.+)$', r'\\section{\1}', t, flags=re.MULTILINE)
    t = re.sub(r'^### (.+)$', r'\\section{\1}', t, flags=re.MULTILINE)  # unnumbered fallback
    t = re.sub(r'^#### \d+\.\d+\.\d+ (.+)$', r'\\subsection{\1}', t, flags=re.MULTILINE)
    t = re.sub(r'^#### (.+)$', r'\\subsection{\1}', t, flags=re.MULTILINE)  # unnumbered fallback
    t = re.sub(r'\*\*(.+?)\*\*', r'\\textbf{\1}', t)
    t = re.sub(r'(?<!\*)\*([^*\n]+?)\*(?!\*)', r'\\textit{\1}', t)
    # Escape special chars BEFORE replacing unicode.
    # $, %, &, # are all special in LaTeX and must be escaped when they
    # appear as literal characters in prose.  The $ escape must run BEFORE
    # the unicode→LaTeX block below (which inserts $...$ math wrappers);
    # escaping first ensures only literal prose $ are hit.
    # NOTE: { and } are NOT escaped here because CITEPLACEHOLDER{key}
    # tokens are still in the text at this point — escaping braces would
    # break the \cite conversion that runs after cvt().  Bare braces in
    # prose will cause LaTeX errors; avoid them in v1.md.
    t = t.replace('$', '\\$')
    t = t.replace('%', '\\%')
    t = t.replace('&', '\\&')
    t = t.replace('#', '\\#')
    # NOTE: underscore escaping is done AFTER all LaTeX conversions
    # (cites, figures, tables) to avoid breaking citation keys and labels.
    # See the escape_prose_underscores() call below cvt().
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

# Convert footnote references [^label] → \footnote{definition text}
def repl_footnote(m):
    label = m.group(1)
    text = footnote_defs.get(label, f'[{label}]')
    # Escape LaTeX special chars in footnote text
    text = text.replace('&', '\\&').replace('%', '\\%').replace('#', '\\#')
    # Handle Unicode that pdflatex can't render directly
    text = text.replace('≤', '$\\leq$').replace('≥', '$\\geq$')
    text = text.replace('×', '$\\times$')
    text = text.replace('\u2009', ' ')  # thin space → regular space
    text = text.replace('—', '---').replace('–', '--')
    text = text.replace('→', '$\\rightarrow$')
    # Accented chars: keep as-is (fontenc T1 handles common Latin accents)
    return f'\\footnote{{{text}}}'

body_tex = re.sub(r'\[\^(\w+)\]', repl_footnote, body_tex)
# Remove any leftover footnote definition lines that survived the earlier cleanup
body_tex = re.sub(r'^\[\^\w+\]:\s*.+$', '', body_tex, flags=re.MULTILINE)

# Remove markdown horizontal rules, which are invalid in LaTeX body text.
body_tex = re.sub(r'^\s*---\s*$', '', body_tex, flags=re.MULTILINE)

# Collapse the runs of blank lines left behind when footnote-definition lines
# and horizontal rules are stripped above (cosmetic; LaTeX already treats any
# blank-line run as a single paragraph break, so this changes no output).
body_tex = re.sub(r'\n{3,}', '\n\n', body_tex)

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

# Hardcoded LaTeX tables intentionally diverge from the markdown source:
# - Simulation table splits SDT/PDT into separate columns (markdown combines
#   them because values are identical).
# - Modality table uses abbreviated headers for column width.
# If the markdown table data changes, these hardcoded versions must be
# updated manually to match.
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
    '3': ('fig1_ferroptosis_comparison', 'Ferroptosis engagement ($\\chi^2=97.3$, $p<5.9\\times10^{-23}$; corpus-derived, subject to tagging and taxonomy uncertainty).'),
    '4': ('fig4_molecular_overlap', 'Molecular pathway engagement (normalized \\%).'),
    '5': ('fig3_literature_disconnect', 'Literature disconnect between communities.'),
    '6': ('fig6_sdt_chain_evidence', 'SDT ferroptosis-ICD chain evidence.'),
    '7': ('fig7_monte_carlo_simulation', 'Monte Carlo simulation (1M cells/condition).'),
    '8': ('fig8_simulation_by_treatment', 'Depth-kill curves: tissue penetration sets modality reach (2D model). (a) Observed tumor kill versus depth across 1 cm of tissue: SDT (ultrasound) stays near 95 to 100\\% throughout, PDT (light) collapses from ${\\sim}93\\%$ at the surface to ${\\sim}0\\%$ by 10 mm (Beer-Lambert attenuation), and RSL3 (systemic drug) is a flat, near-zero baseline at every depth. (b) The driving physics from the model\'s own equations: PDT light decays as $\\exp(-\\mu_{\\mathrm{eff}} z)$ with $\\mu_{\\mathrm{eff}}{=}0.31$/mm ($\\delta{\\approx}3.2$ mm), SDT acoustic as $10^{-\\alpha f z/10}$ with $\\alpha{=}0.7$ dB/cm/MHz at 1 MHz, RSL3 uniform at 100\\%. RSL3 reaches every depth yet kills little, a biochemical limit, not a penetration one. The depth profiles follow well-measured physics (high confidence); absolute kill \\% rests on uncalibrated biochemistry, so the profile shape is the result, not the magnitudes. \\textbf{SDT is modeled as O$_2$-independent, an optimistic upper bound} (Section 7.1).'),
    '9': ('fig13_gold_set_eval', 'Evidence tagger performance: gold-set evaluation (100-article stratified sample).'),
    '10': ('fig9_evidence_tiers', 'Evidence tier composition by mechanism.'),
    '11': ('fig14_tissue_mechanism_heatmap', 'Tissue-of-origin $\\times$ mechanism article counts (coverage: 62\\%).'),
    '12': ('fig15_designed_combinations', 'Classification of multi-mechanism articles into designed combinations, co-mentions, and reviews.'),
    '13': ('fig16_weighted_evidence', 'Weighted evidence score by mechanism (tier $\\times$ citation percentile $\\times$ recency).'),
    '14': ('fig17_damp_heatmap', 'DAMP spatial distribution after immune coupling (O$_2$ gradient $\\lambda$=120$\\mu$m, per-panel scaling---intensity not comparable across panels). SDT covers the full tumor area (139,641 kills, 539 immune kills); RSL3 produces sparse isolated hotspots (163 kills, 2 immune kills).'),
    '15': ('fig18_hypoxia_crosssection', 'Hypoxia cross-section: O$_2$ gradient from blood vessel (left) into tumor core (right). RSL3 efficacy collapses as basal ROS disappears; SDT maintains efficacy via exogenous ROS delivery.'),
    '16': ('fig19_immune_coupling_flow', 'Immune coupling pathway: SDT produces dense kill with high LP overshoot, generating strong DAMP fields and 104$\\times$ more immune kills than RSL3.'),
    '17': ('fig20_stromal_shielding', 'Stromal shielding: CAF-mediated GSH and MUFA supply halves RSL3 kill at the tumor boundary (3.0\\% $\\rightarrow$ 1.5\\%) while barely affecting SDT (96.1\\% $\\rightarrow$ 91.2\\%).'),
    '18': ('fig21_ph_ion_trapping', 'pH-driven ion trapping: acidic tumor core protonates and traps drug molecules, reducing RSL3 kills by 53\\%. SDT is pH-independent (no drug to trap).'),
    '19': ('fig22_decision_flowchart', 'Decision framework: which modality for which clinical context, based on tumor localizability, depth, ferroptosis-prone residual state, and immunocompetence.'),
    '20': ('fig23_prisma_flow', 'PRISMA-inspired corpus construction flow. 10,414 unique PubMed records across 19 mechanism queries; 4,830 full-text articles indexed (803 journals, 2001--2026), 5,584 abstract-only records retained separately. No manual screening or exclusion criteria were applied; this is an automated keyword-driven pipeline, not a formal systematic review.'),
    '21': ('fig24_hypoxia_killcurve', 'Hypoxia kill-collapse (2D model). (a) RSL3 kill collapses from 3.7\\% (normoxic) to ${\\sim}0.1\\%$ (hypoxic) while SDT holds 91.9\\%${\\to}$87.8\\%; (b) the gradient result is flat across O$_2$ penetration length $\\lambda$=80--150$\\mu$m. \\textbf{SDT is modeled as O$_2$-independent, an optimistic upper bound} --- SDT\'s own O$_2$-dependence is contested (Section 7.1), so the direction is more robust than the magnitude of the gap.'),
    '22': ('fig25_bliss_synergy', 'Dual-pathway depletion synergy. (a) RSL3+FSP1i kills 84.1\\%, far above the 42.2\\% Bliss-independent prediction (1.99$\\times$ synergy); (b) pairwise synergy scores (SDT pairs excluded for a 100\\% single-agent ceiling). Drug potencies are estimates; the directional finding (dual-pathway $>$ single) held across the $\\pm$50\\% sensitivity sweep (Section 5).'),
    '23': ('fig26_vulnerability_window', 'The ferroptosis-sensitive treatment window. (a) After chemotherapy, RSL3 kill collapses from 42.4\\% to 1.4\\% by day 3 and to ${\\sim}0$ by day 7 as GPX4 is re-expressed, while SDT holds ${\\sim}100\\%{\\to}99.5\\%$ through day 28; (b) the RSL3 collapse tracks mean GPX4 recovery (twin axis). Defense-recovery half-times (GPX4 3 d, FSP1 7 d, NRF2 5 d, GSH 1 d) are literature-estimated, so the window durations are approximate (medium confidence) until experimentally validated.'),
    '24': ('fig27_resistance_asymmetry', 'The resistance-mechanism asymmetry (2D model, flagship). Under each tumor-microenvironment resistance mechanism, pharmacologic RSL3 collapses while physical SDT holds; each panel uses the same metric its section reports, so the figure and the text agree. (a) Hypoxia (Section 7.1): overall kill, RSL3 falls from 3.7\\% to 0.1\\% (normoxic to hypoxic) while SDT holds 91.9\\% to 87.8\\%. (b) Stromal/CAF (Section 7.3): kill among the CAF-adjacent boundary cells, RSL3 halved from 3.0\\% to 1.5\\% while SDT barely moves (96.1\\% to 91.2\\%). (c) Acidic pH (Section 7.4): ferroptosis kills (an immune-free counter), RSL3 from 163 to 77 (a 53\\% drop) while SDT is unaffected (139,640 to 140,693). (d) Immune/ICD coupling (Section 7.2): SDT produces 521 ICD-driven immune kills versus 5 for RSL3 (104:1 in 2D). Panel (a) is computed without the immune layer (a clean O$_2$-only comparison); panels (b)-(d) share the gradient-O$_2$ plus immune-on baseline the sim runs those mechanisms under (the pH ``neutral\'\' bar reuses the stromal-off run, the only available reference). \\textbf{Confidence tiers differ per panel} (titles): the hypoxia leg is the most contested (SDT is modeled O$_2$-independent, an optimistic upper bound, Section 7.1) and the 2D immune ratio over-extrapolates (${\\sim}4{:}1$ under 3D volumetric dilution). Magnitudes rest on uncalibrated biochemistry; the cross-modality direction, not the numbers, is the result.'),
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

# Escape underscores in prose AFTER all LaTeX conversions (cites, figures,
# tables) are complete. Then un-escape inside \cite{}, \label{}, \ref{},
# and \includegraphics{} commands where underscores are valid.
def escape_prose_underscores(t):
    # Step 1: escape underscores between word characters (e.g., gene_name).
    # Pattern is intentionally narrow: (?<=\w)_(?=\w) catches the common
    # case but misses edge cases like _foo or foo_.  A broader pattern
    # would also match the _ in math-mode subscripts ($_2$, $_0$) inserted
    # by cvt(), breaking them.  The narrow pattern is a trade-off:
    # - Catches: GPX4_activity, SLC7A11_high, file_name
    # - Misses: _italic_ (already handled by bold/italic conversion),
    #   bare _ at word boundaries (rare in scientific prose)
    t = re.sub(r'(?<=\w)_(?=\w)', r'\\_', t)
    # Step 2: un-escape inside LaTeX commands that use underscored keys
    def unescape_braces(m):
        return m.group(0).replace('\\_', '_')
    t = re.sub(r'\\cite\{[^}]+\}', unescape_braces, t)
    t = re.sub(r'\\label\{[^}]+\}', unescape_braces, t)
    t = re.sub(r'\\ref\{[^}]+\}', unescape_braces, t)
    t = re.sub(r'\\includegraphics\[[^\]]*\]\{[^}]+\}', unescape_braces, t)
    return t

body_tex = escape_prose_underscores(body_tex)
abstract_tex = escape_prose_underscores(abstract_tex)

latex = f"""\\documentclass[12pt,a4paper]{{report}}
\\usepackage[utf8]{{inputenc}}
\\usepackage[T1]{{fontenc}}
\\usepackage{{amsmath,amssymb}}
\\usepackage{{graphicx}}
\\usepackage{{hyperref}}
% Citations use inline footnotes — no natbib/bibtex needed
\\usepackage{{booktabs}}
\\usepackage{{geometry}}
\\usepackage{{setspace}}

\\geometry{{margin=1in}}
\\onehalfspacing

\\title{{{title}}}
\\author{{Ezequiel Lares \\\\ Independent Researcher}}
\\date{{}}

\\begin{{document}}

\\maketitle

\\begin{{abstract}}
{abstract_tex}
\\end{{abstract}}

\\noindent\\textbf{{Keywords:}} {keywords}

\\tableofcontents

{body_tex}

% Reference list kept as appendix in v1.md; citations are inline footnotes.

\\end{{document}}
"""

TEX.write_text(latex)
print(f'Written {TEX}: {len(latex)} chars, {latex.count(chr(92)+"cite{")} citations')
