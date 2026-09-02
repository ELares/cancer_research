#!/usr/bin/env python3
r"""Convert v1.md into the book's LaTeX source.

THIS FILE EMITS MEANING, NOT LAYOUT. Every visual decision -- trim size,
typefaces, palette, measure, how a chapter opens, what a callout looks like --
lives in `article/drafts/bookdesign.sty`, and `tests/test_book_design.py` fails
if a layout primitive (\geometry, \fontsize, \vspace, a colour, a font
family) appears here. The split is what makes the design adaptable: a chapter
written next month is set exactly like the twelve that exist, and the book can
be redrawn by editing one style file without touching this pipeline.

What the markdown may contain, and what each construct becomes:

    # Part I: Title          \part      (opener page, standfirst below it)
    ## Chapter 3: Title      \chapter   (opener page, standfirst, drop cap)
    > standfirst text        \standfirst / \partstandfirst when it directly
                             follows a Part/Chapter heading; \pullquote anywhere
                             else
    ::: finding ... :::      a callout box (finding / refusal / numbers)
    1. item / - item         enumerate / itemize
    ```  ...  ```            a terminal block, set verbatim
    `identifier`             \bookcode
    ---  (alone on a line)   a scene break
    [FIGURE n: ...]          the figure, with its caption and pinned number

See article/AUTHORING.md for the authoring rules these mirror.
"""
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
MD = ROOT / "article" / "drafts" / "v1.md"
TEX = ROOT / "article" / "drafts" / "v1.tex"

md = MD.read_text()
full_title = re.search(r'^# (.+)$', md, re.MULTILINE).group(1).strip()
# The book's title and subtitle are the two halves of the manuscript's own
# title line, split rather than retyped so the title page cannot drift from the
# source. A title with no colon simply has no subtitle.
if ':' in full_title:
    title, subtitle = (p.strip() for p in full_title.split(':', 1))
else:
    title, subtitle = full_title, ''

# Fenced code blocks are pulled OUT before any conversion runs and put back
# after, because every step in between escapes LaTeX specials and a shell
# command is not prose. Placeholders are unlikely-to-collide sentinels rather
# than markers made of markdown, which the converters would rewrite.
code_blocks = []


def stash_code_blocks(text):
    def take(m):
        code_blocks.append(m.group(1).rstrip('\n'))
        return f'@@CODEBLOCK{len(code_blocks) - 1}@@'
    return re.sub(r'```[a-z]*\n(.*?)```', take, text, flags=re.DOTALL)


def restore_code_blocks(text):
    def put(m):
        body = code_blocks[int(m.group(1))]
        return ('\\begin{bookterminal}\n\\begin{verbatim}\n'
                + body + '\n\\end{verbatim}\n\\end{bookterminal}')
    return re.sub(r'@@CODEBLOCK(\d+)@@', put, text)


md = stash_code_blocks(md)

# Build footnote definition map from Markdown [^label]: text patterns
footnote_defs = {}
for m in re.finditer(r'^\[\^(\w+)\]:\s*(.+)$', md, re.MULTILINE):
    footnote_defs[m.group(1)] = m.group(2).strip()

# Extract sections
abstract = re.search(r'## Abstract\n\n(.*?)(?=\n\*\*Keywords)', md, re.DOTALL).group(1).strip()
keywords = re.search(r'\*\*Keywords:\*\*\s*(.+)', md).group(1).strip()

# Front matter written as prose in the manuscript. Optional: absent until the
# section is written, and the book simply omits the page rather than failing.
_note = re.search(r'^## A Note On This Edition\n\n(.*?)(?=\n# Part )', md,
                  re.DOTALL | re.MULTILINE)
howto = _note.group(1).strip() if _note else ''

# Body: everything from the first Part header to just before References.
body_start = re.search(r'^# Part [IVX]+: ', md, re.MULTILINE)
ref_match = re.search(r'^## References', md, re.MULTILINE)
if not body_start or not ref_match:
    raise SystemExit("ERROR: Could not find '# Part ...' or '## References' boundaries in v1.md")
body = md[body_start.start():ref_match.start()].strip()

# Remove footnote definition lines from the body (they'll become \footnote{} inline)
body = re.sub(r'^\[\^\w+\]:\s*.+$', '', body, flags=re.MULTILINE)

# Greek letters → LaTeX math commands. pdflatex with inputenc=utf8 does NOT map
# bare Greek code points to glyphs (unlike accented Latin / ± / curly quotes,
# which it does), so an unmapped Greek letter is a hard `! LaTeX Error: Unicode
# character ...` that makes pdflatex exit non-zero and breaks the release-pdf
# build even though a PDF is produced. Mapping every Greek letter with a distinct
# LaTeX command lets the prose (and citation footnotes) use them freely. Capital
# letters whose glyph is identical to a Latin capital (Α, Β, Ε, ...) have no
# \Command and are intentionally omitted; use the Latin letter for those.
_GREEK_TO_LATEX = (
    ('α', 'alpha'), ('β', 'beta'), ('γ', 'gamma'), ('δ', 'delta'),
    ('ε', 'epsilon'), ('ζ', 'zeta'), ('η', 'eta'), ('θ', 'theta'),
    ('ι', 'iota'), ('κ', 'kappa'), ('λ', 'lambda'), ('ν', 'nu'),
    ('ξ', 'xi'), ('π', 'pi'), ('ρ', 'rho'), ('σ', 'sigma'),
    ('τ', 'tau'), ('φ', 'phi'), ('χ', 'chi'), ('ψ', 'psi'), ('ω', 'omega'),
    ('Γ', 'Gamma'), ('Δ', 'Delta'), ('Θ', 'Theta'), ('Λ', 'Lambda'),
    ('Ξ', 'Xi'), ('Π', 'Pi'), ('Σ', 'Sigma'), ('Υ', 'Upsilon'),
    ('Φ', 'Phi'), ('Ψ', 'Psi'), ('Ω', 'Omega'),
)


def map_greek(t):
    """Replace bare Greek code points with `$\\command$`. `μ` (U+03BC) is handled
    separately by the caller (it shares a command with the micro sign U+00B5)."""
    for greek, cmd in _GREEK_TO_LATEX:
        t = t.replace(greek, f'$\\{cmd}$')
    return t


# Markdown → LaTeX
def cvt(t):
    # Inline literals FIRST, before the escapers touch them: a backticked
    # `identifier` was reaching pdflatex as a backtick, which sets an opening
    # quote, so 205 file names and enum variants printed as mismatched quote
    # marks. Nothing warned; the PDF simply had a typographic error in it 205
    # times.
    t = re.sub(r'`([^`\n]+)`', r'\\bookcode{\1}', t)
    # Book-structure headings (report document class)
    t = re.sub(r'^# Part [IVX]+: (.+)$', r'\\part{\1}', t, flags=re.MULTILINE)
    t = re.sub(r'^## Chapter (\d+): (.+)$',
               r'\\booknotesheading{Chapter \1}\n\\chapter{\2}', t,
               flags=re.MULTILINE)
    # Appendices. `## Appendix A: Title` matched NO rule, so all three reached
    # the page as the literal text "## Appendix A: ..." -- a heading printed as
    # its own markup, in the compiled book, for as long as the appendices have
    # existed. \appendix switches the chapter counter to letters, so the
    # lettered section numbers the prose already uses (A.1, B.3) are produced
    # rather than typed twice.
    t = re.sub(r'^## Appendix A: (.+)$', r'\\appendix\n\\chapter{\1}', t,
               flags=re.MULTILINE)
    t = re.sub(r'^## Appendix [B-Z]: (.+)$', r'\\chapter{\1}', t,
               flags=re.MULTILINE)
    t = re.sub(r'^### [A-Z]\.\d+ (.+)$', r'\\section{\1}', t, flags=re.MULTILINE)
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
    # `^` is a superscript command in text mode, so a formula written in prose
    # -- `exp(-alpha*D - beta*D^2)` -- reaches pdflatex as an unclosed group
    # and the document compiles with a group still open. The negative
    # lookbehind spares `[^label]`, which is a footnote reference and is
    # converted later.
    t = re.sub(r'(?<!\[)\^', r'\\textasciicircum{}', t)
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
    # Straight quotation marks set as two apostrophes on both sides, which is
    # the oldest giveaway that a document was typed rather than typeset. Paired
    # left-to-right so an odd one out is left alone rather than guessed at.
    def _curly(m):
        return "``" + m.group(1) + "''"
    t = re.sub(r'"([^"]{1,400})"', _curly, t, flags=re.DOTALL)
    t = t.replace('≥', '$\\geq$')
    t = t.replace('≤', '$\\leq$')
    t = t.replace('≳', '$\\gtrsim$')
    t = t.replace('≲', '$\\lesssim$')
    t = t.replace('≫', '$\\gg$')
    t = t.replace('≪', '$\\ll$')
    t = t.replace('≈', '$\\approx$')
    t = t.replace('†', '$\\dagger$')   # table footnote marker
    t = t.replace('‡', '$\\ddagger$')
    t = t.replace('µ', '$\\mu$')
    t = t.replace('μ', '$\\mu$')      # U+03BC (Greek mu) — distinct from U+00B5 (micro sign)
    t = map_greek(t)                  # all remaining Greek letters (incl. Σ, Π, Δ, ...)
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
    text = text.replace('≲', '$\\lesssim$').replace('≳', '$\\gtrsim$')
    text = text.replace('≫', '$\\gg$').replace('≪', '$\\ll$')
    text = text.replace('×', '$\\times$')
    text = text.replace('\u2009', ' ')  # thin space → regular space
    text = text.replace('—', '---').replace('–', '--')
    text = re.sub(r'"([^"]{1,400})"', lambda m: '``' + m.group(1) + "''",
                  text, flags=re.DOTALL)
    text = text.replace('→', '$\\rightarrow$')
    text = text.replace('µ', '$\\mu$').replace('μ', '$\\mu$')
    text = map_greek(text)  # Greek in citation titles/authors (same hard-error class as the body)
    # Accented chars: keep as-is (fontenc T1 handles common Latin accents)
    return f'\\footnote{{{text}}}'

body_tex = re.sub(r'\[\^(\w+)\]', repl_footnote, body_tex)
# A note mark belongs against the word it annotates. The markdown writes a
# space before `[^ref]`, which reached the page as a gap between the word and
# its superscript -- small, and on every one of several hundred marks.
body_tex = re.sub(r'[ \t]+(\\footnote\{)', r'\1', body_tex)
# Remove any leftover footnote definition lines that survived the earlier cleanup
body_tex = re.sub(r'^\[\^\w+\]:\s*.+$', '', body_tex, flags=re.MULTILINE)

# A rule alone on a line is a scene break -- a shift of subject inside a
# section. It used to be deleted, which silently joined two passages that the
# author had separated.
body_tex = re.sub(r'^\s*---\s*$', r'\\dinkus', body_tex, flags=re.MULTILINE)

# Collapse the runs of blank lines left behind when footnote-definition lines
# and horizontal rules are stripped above (cosmetic; LaTeX already treats any
# blank-line run as a single paragraph break, so this changes no output).
body_tex = re.sub(r'\n{3,}', '\n\n', body_tex)

# ── Structure: standfirsts, callouts, lists, drop caps ──────────────────
# Everything below turns a markdown convention into a SEMANTIC command that
# bookdesign.sty draws. None of it decides how the result looks.

# Fence name -> environment. The label and its colour belong to the
# environment, which lives in the style file; a fence may pass a custom label
# as an optional argument and nothing here decides how it is set.
CALLOUTS = {
    'finding': 'finding',
    'refusal': 'refusal',
    'numbers': 'numbers',
}


def convert_callouts(t):
    """`::: finding [label] ... :::` -> a callout environment.

    An unknown fence name is a hard error rather than a passthrough: a
    misspelled `::: findng` would otherwise print three colons and its label
    into the running text, where it reads as a typesetting accident rather
    than as the author's mistake.
    """
    pattern = re.compile(r'^::: *(\S+) *(.*?)\n(.*?)\n::: *$',
                         re.DOTALL | re.MULTILINE)

    def one(m):
        name, custom, inner = m.group(1), m.group(2).strip(), m.group(3).strip()
        if name not in CALLOUTS:
            raise SystemExit(
                f"ERROR: unknown callout ':::  {name}' in v1.md. "
                f"Known: {', '.join(sorted(CALLOUTS))}. Add it to CALLOUTS "
                "here AND define the environment in bookdesign.sty.")
        env = CALLOUTS[name]
        opt = f'[{custom}]' if custom else ''
        return f'\\begin{{{env}}}{opt}\n{inner}\n\\end{{{env}}}'

    t = pattern.sub(one, t)
    leftover = re.findall(r'^:::.*$', t, flags=re.MULTILINE)
    if leftover:
        raise SystemExit(
            "ERROR: unclosed callout fence in v1.md: " + leftover[0][:70])
    return t


def convert_blockquotes(t):
    """A blockquote is a standfirst when it opens a Part or Chapter, and a
    pull quote anywhere else.

    The position carries the meaning, which is why the two cases are not
    distinguished by two markdown syntaxes: an author writing the sentence
    that goes under a chapter title should not have to remember a second
    spelling for it."""
    lines = t.split('\n')
    out, i = [], 0
    while i < len(lines):
        line = lines[i]
        if line.startswith('> '):
            block, j = [], i
            while j < len(lines) and lines[j].startswith('> '):
                block.append(lines[j][2:].strip())
                j += 1
            text = ' '.join(block).strip()
            prev = next((x for x in reversed(out) if x.strip()), '')
            if prev.startswith('\\part{'):
                # Fold the two into one command: a part page ends with a
                # \newpage, so a standfirst emitted after it opens the next
                # leaf instead of finishing the opener.
                title = prev[len('\\part{'):].rstrip()[:-1]
                for k in range(len(out) - 1, -1, -1):
                    if out[k].strip() == prev.strip():
                        out[k] = f'\\bookpart{{{title}}}{{{text}}}'
                        break
            elif prev.startswith('\\chapter{'):
                out.append(f'\\standfirst{{{text}}}')
            else:
                out.append(f'\\pullquote{{{text}}}')
            i = j
            continue
        out.append(line)
        i += 1
    return '\n'.join(out)


def convert_lists(t):
    """Consecutive numbered or bulleted PARAGRAPHS become one list.

    Items in this manuscript are blank-line separated, so a naive line-based
    grouper sees each as its own list. They were not converted at all before:
    135 numbered items reached the page as ordinary paragraphs beginning with
    a digit and a full stop."""
    blocks = t.split('\n\n')
    out, i = [], 0
    while i < len(blocks):
        b = blocks[i].strip()
        num = re.match(r'^\d+\.\s+(.*)$', b, re.DOTALL)
        bul = re.match(r'^[-*]\s+(.*)$', b, re.DOTALL)
        if num or bul:
            env = 'enumerate' if num else 'itemize'
            items, j = [], i
            while j < len(blocks):
                nb = blocks[j].strip()
                m = (re.match(r'^\d+\.\s+(.*)$', nb, re.DOTALL) if num
                     else re.match(r'^[-*]\s+(.*)$', nb, re.DOTALL))
                if not m:
                    break
                items.append(m.group(1).strip())
                j += 1
            body = '\n'.join(f'\\item {it}' for it in items)
            out.append(f'\\begin{{{env}}}\n{body}\n\\end{{{env}}}')
            i = j
            continue
        out.append(blocks[i])
        i += 1
    return '\n\n'.join(out)


def add_drop_caps(t):
    """Open each chapter's first paragraph with a drop cap.

    Skipped -- deliberately, and silently -- when the paragraph opens with
    anything but a plain capital letter: a bold lead-in, a quantity, a
    quotation. A drop cap on `\\textbf{` would set a backslash three lines
    high, and the alternative (unwrapping the markup to reach the letter)
    would put layout logic in this file."""
    blocks = t.split('\n\n')
    for i, b in enumerate(blocks):
        # The chapter command shares its block with the notes heading emitted
        # just above it, so this looks INSIDE the block rather than at its
        # first characters.
        if not re.match(r'(\\booknotesheading\{[^}]*\}\s*)?\\chapter\{', b.lstrip()):
            continue
        for j in range(i + 1, min(i + 6, len(blocks))):
            nxt = blocks[j].strip()
            # A chapter here always opens with a numbered section, so the
            # first PARAGRAPH is two or three blocks down; the search skips
            # the furniture rather than giving up at the first thing that is
            # not prose.
            if (not nxt or nxt.startswith('\\standfirst{')
                    or nxt.startswith('\\section{')
                    or nxt.startswith('\\booknotesheading{')):
                continue
            m = re.match(r'^([A-Z])([a-z]+)(.*)$', nxt, re.DOTALL)
            if m and len(m.group(2)) >= 2:
                blocks[j] = (f'\\bookopen{{{m.group(1)}}}{{{m.group(2)}}}'
                             f'{{{m.group(3)}}}')
            break
    return '\n\n'.join(blocks)


body_tex = convert_callouts(body_tex)
body_tex = convert_blockquotes(body_tex)
body_tex = convert_lists(body_tex)
body_tex = add_drop_caps(body_tex)

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

    table = f"""\\begin{{table}}[tbp]
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
    '1': ('fig1c_ratio_straddle', 'The pharmacological:physical volume ratio under four class definitions, against the 9.1:1 an earlier keyword method reported. Two of the four exceed it and two fall below, so the DIRECTION of the imbalance survives every reading while the claim that the census understates the earlier case does not. Which definition is chosen decides the direction of that second conclusion, which is why all four are drawn.'),
    '2': ('fig2c_census_volume', 'Census publication volume by year, with the retrieved corpus overlaid on a second axis. The corpus rises 31-fold over the decade the census grows 1.10-fold; separate axes because the two series differ by orders of magnitude, so their shapes are comparable and their heights are not.'),
    '3': ('fig1_ferroptosis_comparison', 'Ferroptosis engagement ($\\chi^2=38.8$, $p=4.7\\times10^{-10}$; corpus-derived, subject to tagging and taxonomy uncertainty).'),
    '4': ('fig4_molecular_overlap', 'Molecular pathway engagement (normalized \\%).'),
    '5': ('fig5c_mechanism_site_matrix', 'Mechanism by anatomical site across the census, coloured by observed over expected under independent marginals rather than by article count -- a count heatmap reproduces the product of the marginals, so its brightest cell is always the largest mechanism crossed with the largest site. Cells whose expectation falls below 20 articles are grey, because a ratio computed on a handful is noise and colouring it would put the loudest colours on the least reliable cells. The concentrated cells recover known clinical concentrations from descriptors alone.'),
    '6': ('fig6_sdt_chain_evidence', 'SDT ferroptosis-ICD chain evidence.'),
    '7': ('fig7_monte_carlo_simulation', 'Monte Carlo simulation (1M cells/condition).'),
    '8': ('fig8_simulation_by_treatment', 'Depth-kill curves: tissue penetration sets modality reach (2D model). (a) Observed tumor kill versus depth across 1 cm of tissue: SDT (ultrasound) stays near 95 to 100\\% throughout, PDT (light) collapses from ${\\sim}93\\%$ at the surface to ${\\sim}0\\%$ by 10 mm (Beer-Lambert attenuation), and RSL3 (systemic drug) is a flat, near-zero baseline at every depth. (b) The driving physics from the model\'s own equations: PDT light decays as $\\exp(-\\mu_{\\mathrm{eff}} z)$ with $\\mu_{\\mathrm{eff}}{=}0.31$/mm ($\\delta{\\approx}3.2$ mm), SDT acoustic as $10^{-\\alpha f z/10}$ with $\\alpha{=}0.7$ dB/cm/MHz at 1 MHz, RSL3 uniform at 100\\%. RSL3 reaches every depth yet kills little, a biochemical limit, not a penetration one. The depth profiles follow well-measured physics (high confidence); absolute kill \\% rests on uncalibrated biochemistry, so the profile shape is the result, not the magnitudes. \\textbf{SDT is modeled as O$_2$-independent, an optimistic upper bound} (Section 7.1).'),
    '9': ('fig13_gold_set_eval', 'Evidence tagger performance: gold-set evaluation (100-article stratified sample).'),
    '10': ('fig9c_design_composition', 'Study-design composition of the census, from NLM publication types and MeSH check tags. The undetermined class is drawn rather than dropped: at 44.5\\% it is the largest, and omitting it would imply the census assigns a design to every record.'),
    '11': ('fig14c_class_by_site', 'Mechanism class by anatomical site, as enrichment against each site\'s own share of site-assigned records. Bold labels mark the sites where the physical and pharmacological classes move in opposite directions, the reading that does not depend on how much a site is written about. The physical class omits radiotherapy, its largest real member.'),
    '12': ('fig15c_mechanism_pairs', 'The ten most frequent mechanism pairs in the census. Counts, not a rate: co-tagging records that two vocabularies appear on one article, and the co-occurrence rate is a property of the labelling instrument rather than of the field.'),
    '13': ('fig16c_trial_share', 'Clinical-trial share against article volume (log axis). HIFU sits above CAR-T on a fraction of the volume, so physical modality is not a maturity class; nanoparticle delivery sits at high volume and 0.49\\%.'),
    '14': ('fig17_damp_heatmap', 'DAMP spatial distribution after immune coupling, baseline run (O$_2$ gradient $\\lambda$=120$\\mu$m, no stromal shielding, no pH gradient; per-panel scaling---intensity not comparable across panels). SDT covers the full tumor area (139,640 kills, 521 immune kills); RSL3 produces sparse isolated hotspots (163 kills, 5 immune kills).'),
    '15': ('fig18_hypoxia_crosssection', 'Hypoxia cross-section: O$_2$ gradient from blood vessel (left) into tumor core (right). RSL3 efficacy collapses as basal ROS disappears; SDT maintains efficacy via exogenous ROS delivery.'),
    '16': ('fig19_immune_coupling_flow', 'Immune coupling pathway: SDT produces dense kill with high LP overshoot, generating strong DAMP fields and 104$\\times$ more immune kills than RSL3.'),
    '17': ('fig20_stromal_shielding', 'Stromal shielding: CAF-mediated GSH and MUFA supply halves RSL3 kill at the tumor boundary (3.0\\% $\\rightarrow$ 1.5\\%) while barely affecting SDT (96.1\\% $\\rightarrow$ 91.2\\%).'),
    '18': ('fig21_ph_ion_trapping', 'pH-driven ion trapping: acidic tumor core protonates and traps drug molecules, reducing RSL3 kills by 53\\%. SDT is pH-independent (no drug to trap).'),
    '19': ('fig22_decision_flowchart', 'Decision framework: which modality for which clinical context, based on tumor localizability, depth, ferroptosis-prone residual state, and immunocompetence.'),
    '20': ('fig23_census_flow', 'Census construction. The PubMed annual baseline (1,334 XML files) filtered to MeSH tree C04 (704 descriptors) union nine adjacent experimental-context descriptors, yielding 4,403,994 MeSH-indexed cancer articles, plus a second stream of 783,271 records recovered by text matching where MeSH indexing has not yet been applied (matcher precision 75.7\\%, recall 95.6\\%): 5,187,265 in total. The two streams are parallel admissions under different rules, not a branch after an exclusion -- nothing is screened out.'),
    '21': ('fig24_hypoxia_killcurve', 'Hypoxia kill-collapse (2D model). (a) RSL3 kill collapses from 3.7\\% (normoxic) to ${\\sim}0.1\\%$ (hypoxic) while SDT holds 91.9\\%${\\to}$87.8\\%; (b) the gradient result is flat across O$_2$ penetration length $\\lambda$=80--150$\\mu$m. \\textbf{SDT is modeled as O$_2$-independent, an optimistic upper bound} --- SDT\'s own O$_2$-dependence is contested (Section 7.1), so the direction is more robust than the magnitude of the gap.'),
    '22': ('fig25_bliss_synergy', 'Dual-pathway depletion synergy. (a) RSL3+FSP1i kills 84.1\\%, far above the 42.2\\% Bliss-independent prediction (1.99$\\times$ synergy); (b) pairwise synergy scores (SDT pairs excluded for a 100\\% single-agent ceiling). Drug potencies are estimates; the directional finding (dual-pathway $>$ single) held across the $\\pm$50\\% sensitivity sweep (Section 5).'),
    '23': ('fig26_vulnerability_window', 'The ferroptosis-sensitive treatment window. (a) After chemotherapy, RSL3 kill collapses from 42.4\\% to 1.4\\% by day 3 and to ${\\sim}0$ by day 7 as GPX4 is re-expressed, while SDT holds ${\\sim}100\\%{\\to}99.5\\%$ through day 28; (b) the RSL3 collapse tracks mean GPX4 recovery (twin axis). Defense-recovery half-times (GPX4 3 d, FSP1 7 d, NRF2 5 d, GSH 1 d) are literature-estimated, so the window durations are approximate (medium confidence) until experimentally validated.'),
    '24': ('fig27_resistance_asymmetry', 'The resistance-mechanism asymmetry (2D model, flagship). Under each tumor-microenvironment resistance mechanism, pharmacologic RSL3 collapses while physical SDT holds; each panel uses the same metric its section reports, so the figure and the text agree. (a) Hypoxia (Section 7.1): overall kill, RSL3 falls from 3.7\\% to 0.1\\% (normoxic to hypoxic) while SDT holds 91.9\\% to 87.8\\% only under the contested O$_2$-independent assumption, so the hypoxic-zone SDT result is best read as the bracket 0\\% to 86.6\\% (it collapses to roughly 0\\% if SDT\'s ROS is instead fully O$_2$-dependent, the regime the lead clinical agent occupies). (b) Stromal/CAF (Section 7.3): kill among the CAF-adjacent boundary cells, RSL3 halved from 3.0\\% to 1.5\\% while SDT barely moves (96.1\\% to 91.2\\%). (c) Acidic pH (Section 7.4): ferroptosis kills (an immune-free counter), RSL3 from 163 to 77 (a 53\\% drop) while SDT is unaffected (139,640 to 140,693). (d) Immune/ICD coupling (Section 7.2): SDT produces 521 ICD-driven immune kills versus 5 for RSL3 (104:1 in 2D). Panel (a) is computed without the immune layer (a clean O$_2$-only comparison); panels (b)-(d) share the gradient-O$_2$ plus immune-on baseline the sim runs those mechanisms under (the pH ``neutral\'\' bar reuses the stromal-off run, the only available reference). \\textbf{Confidence tiers differ per panel} (titles): the hypoxia leg is the most contested (SDT is modeled O$_2$-independent, an optimistic upper bound, Section 7.1) and the 2D immune ratio over-extrapolates (${\\sim}4{:}1$ under 3D volumetric dilution). Magnitudes rest on uncalibrated biochemistry; the cross-modality direction, not the numbers, is the result.'),
    '25': ('fig30_modality_landscape', 'What the engine can be asked, against what the field publishes. Panel (a) is per-mechanism census volume against engine representation; panel (b) changes subject once nothing is absent, because a column emptied to zero invites a reader to take it as a result. The informative split is then APPLICABILITY -- how many mechanisms a run can actually select -- which is a harder number and a smaller one.'),
    '26': ('fig31_modality_panel', 'Every applicable arm against the identical tumour from the identical seed, so a difference between two rows is the arm and nothing else. Bars are coloured by ROUTE rather than ranked, because a chart of kill fractions alone invites exactly the reading the analysis refuses. All but radiation\\textquotesingle s DNA channel are uncalibrated placeholders.'),
    '27': ('fig32_modality_tme', 'What the microenvironment does to every arm, split by cell state and drawn on a SIGNED scale. Red is a loss and blue a gain, and the gains are real: clonal heterogeneity supplies a low-glutathione tail that dies while the average cell resists. Which axis bites depends on the cell state, which is why one phenotype alone reported two axes inert -- a statement about the run, not the biology.'),
    '28': ('fig33_adoptive_barriers', 'The same CAR-T construct against a leukaemia and a solid tumour, on a log axis. No single step looks catastrophic: three barriers leave six per cent, exhaustion removes most of the remainder, and the antigen ceiling is drawn although it contributes nothing here, because it is a cap rather than a coefficient. Every barrier value is an uncalibrated placeholder.'),
    '29': ('fig29_rare_event_resolution', 'How far down the death-rate tail a given sample size can see. Several conditions report exactly 0\\% at one million cells; because the support is unbounded, such a rate is a statement about the sample size rather than about the biology, and it can be pushed down by simulating more cells. The curve is the resolution floor, not a result.'),
    '30': ('fig34_depth_reach', 'How far each modality reaches into tissue. (a) Kill at the surface and at 9.4 mm: photodynamic therapy falls from 92.5\\% to 0.6\\% while megavoltage radiation goes 45.4\\% to 43.7\\%, and the pharmacologic arm does not attenuate at all. (b) Retention as a share of surface kill. Both panels share one order, sorted by kill at depth, and the UNTREATED arm is the TALLEST bar in (b) at 400\\% -- a ratio of two near-zero numbers is not robustness, and the two panels must be read together or neither. Attenuation is fixed physics; every kill magnitude rests on uncalibrated biochemistry.'),
    '31': ('fig35_calibration_verdicts', 'What a published target could settle, per arm. Drawn because the informative rows are the THREE that failed: a thermal-ablation arm whose target admits almost the whole scanned range, checkpoint blockade, whose response band constrains a product of two factors neither of which is identifiable from it, and the antibody-drug-conjugate bystander effect, which has no published target at all. ADMISSIBLE means a parameter reproduces a band, not that the arm is validated -- none of these feeds a number the quantitative chapters of this manuscript report.'),
    '32': ('fig36_fractionation', 'The schedule, and the two things it can be checked against. Panel (a) is the radiation arm\'s external check, and it runs backwards: two schedules a trial reported as not differing imply the $\\alpha/\\beta$ at which they are equivalent, and that value is compared against estimates derived from other data (shaded bands). Prostate lands inside its band at 2.29 Gy; breast crosses at 0.70 Gy, below any plausible tissue, because its shorter arm delivers less total dose and was still not inferior --- a statement about what EQD2 leaves out rather than about the trial. Panels (b) and (c) are DIRECTION-only: the late-responding $\\alpha/\\beta$ is a convention and the reoxygenation half-life is a free parameter no dataset here constrains.'),
    '33': ('fig37_chemotherapy', 'Chemotherapy\'s two structural predictions, neither needing a fitted potency. (a) Dose-response by class in a proliferating population, on a log axis: the phase-nonspecific agent falls away while the phase-specific agents flatten, because raising the dose cannot reach a cell that is not in the phase. (b) The residue each phase-specific class leaves relative to a phase-nonspecific agent at the same dose --- the second pair is the result that was not designed in, because the gap NARROWS when fewer cells are dividing. (c) The dose-density window: the burden ratio between a 21-day and a 14-day schedule at the same total dose, against the Gompertz regrowth rate, with an interior peak and nothing at either end. No dose-response here is fitted --- the CTRPv2 route is access-blocked --- so every absolute value is a placeholder and only the shapes are results.'),
    '34': ('fig38_checkpoint', 'What a ratio can constrain that an absolute response rate cannot. (a) The model\'s antigenicity response saturates with mutational burden, with the two KEYNOTE-158 strata marked; a trial can refute the RATIO between them and cannot refute either height, because the mapping from a model kill to a radiological response is unknown. (b) The region of shape parameters the measured band admits --- 31\\% of the grid scanned --- with the shipped constants marked. (c) The limit, drawn at the same size as the result: the trial reports strata rather than tumours, and the representative-burden choice moves the model\'s answer across most of the target band.'),
    '35': ('fig39_adoptive_escalation', 'What escalating a CAR-T dose buys, and when it buys nothing. (a) A CAR needs a minimum target density to trigger lysis at all, so engagement is a THRESHOLD rather than a gradient --- which is what makes a density failure different in kind from every barrier this module otherwise models. (b) Two tumours with identical barriers and identical poor outcomes: a tenfold dose increase multiplies the delivery-limited one by ten and the density-limited one by less than one per cent. (c) Expansion is driven by the antigen it then consumes, so peak expansion tracks burden. Nothing here is fitted, and the separation in (b) is a prediction about an experiment rather than a reproduction of one.'),
    '36': ('fig40_oncolytic_bind', 'The oncolytic double bind, and the condition that decides it. (a) The two arms as separate terms: virus surviving to spread falls with immune competence, anti-tumour priming rises with it, and neither is non-monotonic. (b) Their sum is, but only sometimes --- with weak priming the best outcome is at full suppression, with strong priming it is in the interior. (c) Where the switch happens, with the saturating regime in grey because there the outcome has hit its ceiling and the optimum no longer reflects the trade-off. The crossover\'s EXISTENCE is the claim; its position is in units this model invented and nothing here constrains it.'),
    '37': ('fig41_adc_loading', 'An optimum that falls out of somebody else\'s measurement. (a) The two clearance ratios measured by Hamblett 2004, with the piecewise interpolation through both and the single power law that misses one of them by a third. (b) Delivered payload per unit antibody dose --- loading divided by clearance --- peaking at a drug-antibody ratio of four, which nothing here was tuned to produce and which is where the field settled for this payload class. (c) The in-vitro and in-vivo orderings disagreeing: potency rises monotonically with loading in a dish where nothing is cleared, while delivery peaks and falls in an animal. Delivered payload is not efficacy, and one conjugate is not all conjugates.'),
    '38': ('fig42_ablation_sleeve', 'Where a thermal ablation fails, and where electroporation does not fail. (a) Blood flowing in a vessel holds nearby tissue near body temperature, so the peak temperature rises with distance from the VESSEL rather than from the applicator. (b) The surviving sleeve as a radius rather than a coverage percentage --- shrinking as the applicator gets hotter and never closing --- against electroporation\'s zero, which is structural rather than measured because it is not a thermal modality. The coolest applicator fails everywhere rather than leaving a sleeve, and is marked as such. (c) Why the size is not a result: the cooling length stands in for vessel calibre and flow, which this layer does not represent.'),
    '39': ('fig43_sonodynamic_frequency', 'An optimum the model finds, and a scaling an independent study contradicts. (a) Three frequency dependences act at once and two oppose: focusing favours going high, while attenuation and the mechanical index\'s own inverse square root both favour going low, so the product peaks in between. (b) The disagreement --- this model\'s optimum falls threefold across 50 to 150 mm, crossing the near-flat 750\\,kHz an independent full-wave study reports at all three depths. No attenuation coefficient makes a one-over-depth law flat, so it is a refutation rather than a tuning gap. (c) Why this arm is structurally unlike the dose-response arms: inertial cavitation is a THRESHOLD, so below the line no exposure time produces anything. The threshold height is a parameter, not a measurement.'),
    '40': ('fig44_pdt_fluence_rate', 'Photodynamic therapy consumes the oxygen it needs, so the same light delivered faster does less. (a) The two opposing monotonic effects separately: oxygen depletion punishes speed, sensitizer clearance punishes slowness. (b) Their product, an interior optimum whose position moves with the drug\'s half-life --- a statement about drugs rather than doses, which is why it is the registered prediction. (c) The refusal: that position scales as roughly the square root of the critical fluence rate, measured nowhere in this project, so the ordering is the result and the milliwatt figures restate an assumption.'),
    '41': ('fig45_radiation_oer', 'Radiation in the spatial engine: an oxygen effect measured from a population rather than restated from a formula. (a) The same single-fraction dose against an oxygenated rim and a hypoxic core. (b) The dose-modifying factor against the steepness of the oxygen gradient --- a parameter the Alper--Howard-Flanders formula has no term for --- moving from 1.24 to 2.55 and peaking INTERIOR to the swept range, against the 2.86 the formula restated returns regardless. (c) Why every value is a lower bound: one gradient length sets the rim and the core together, so a gradient steep enough to make the core anoxic leaves the rim hypoxic too, and the model can never present the fully-oxic-versus-anoxic pair the published band was measured on.'),
    '42': ('fig46_oncolytic_percolation', 'Oncolytic spread is a threshold, and a closed-form front speed cannot express one. (a) The front either makes a local cluster or spans the tumour, with almost nothing between; crosses mark runs where no enterable cell existed near the seed at all --- a failure to START rather than to cross. (b) The permissive fraction needed to cross, against transmission probability, converging on the 26-neighbour site-percolation constant as transmission approaches certainty: the signature of site--bond percolation. (c) The Fisher--KPP speed the engine already contained is positive on every one of these conditions.'),
    '43': ('fig47_adc_bystander_reach', 'A bystander payload multiplies a conjugate\'s reach most where the conjugate reaches least. Every point is a matched pair of runs differing only in the linker, so the bystander\'s contribution is measured rather than supplied. (a) is the prediction and the condition on it together: the advantage falls with penetration at a payload reach of two or three cells and is flat at one, because a payload reaching only its immediate neighbours lands on cells the conjugate reached anyway. (b) is the column that points the other way --- the absolute count mostly rises --- reported beside the ratio rather than instead of it. (c) is why the arm needed a distance at all: the scalar the point model returns has no term for penetration.'),
    '44': ('fig48_chemo_decomposition', 'Two effects push a chemotherapy kill to the rim, and they have opposite remedies. (a) is the control: with neither term active the kill is flat with depth, which is what makes the other columns terms rather than artifacts of the zoning. (b) Delivery dominates on these parameters, taking the core kill to a fraction of a percent in every class --- the penetration length was chosen, not measured, and a longer one shifts the balance. (c) is the half better delivery cannot fix, and the check that the cycle arm is really the cell cycle: an alkylator damages DNA whatever the cell is doing and keeps most of its core kill, while the phase-specific agents lose most of theirs.'),
    '45': ('fig49_cart_independence', 'Multiplying the barrier fractions is right only when the barriers are independent. The spatial run replaces the infiltration scalar with a depth field of the same mean, so the two models differ in distribution and nothing else. (b) is the control and the finding at once: at zero antigen--depth correlation they agree to under one per cent, which is what makes the rest a finding rather than a bug --- the first version of this measurement double-counted the infiltration barrier and its control failed everywhere. (c) is a bucket the product has no slot for: cells that failed both ways, which improving trafficking would not rescue.'),
    '46': ('fig50_checkpoint_priming', 'Checkpoint blockade cannot start a response, and a fold-change cannot say whether one matters. (a) and (b) are the same runs: the fold-change is identical for both active treatments, because a multiplier does not know what it multiplies --- so the quantity the point model reports is the one that cannot tell the two cases apart --- while the share differs thirty-one fold. The untreated arm gains nothing at any blockade strength. (c) is the warning: the treatment producing four times the danger signal has the smaller immune share, because it already killed almost everything directly.'),
    '47': ('fig51_ablation_superposition', 'Heat sinks do not take turns, and a one-vessel radius cannot say so. Both models in (a) run on the same cells and the same vessels, differing only in whether the non-nearest vessels cool. The agreement at wide spacing in (b) is a control: there the extra cooling factors are negligible and the analytic model is simply right, which is what licenses reading the divergence elsewhere as physics. The shaded tail is excluded rather than read as renewed agreement --- there the ablation fails everywhere, and a ratio across a saturated arm is not a comparison. (c) is the mechanism: solid curves count both vessels, dotted count only the nearer.'),
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
        # PIN THE PRINTED NUMBER TO THE MANUSCRIPT NUMBER. LaTeX numbers floats
        # in document order, so the four Chapter 6 figures printed as 13-16
        # while the prose beside them said 25-28 -- the deliverable was true of
        # the markdown reading and false of the compiled PDF, and adding any
        # float shifted every later figure's number. `setcounter` to N-1 makes
        # the environment's own increment land on N, so `\\ref` and the literal
        # prose agree.
        return f"""\\begin{{figure}}[tbp]
\\setcounter{{figure}}{{{int(num) - 1}}}
\\bookgraphic{{../figures/{fn}.pdf}}
\\caption{{\\figlabel{{{num}}}{cap}}}
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
    # \bookgraphic wraps the image now, and it takes the same kind of argument:
    # a PATH, where an escaped underscore is not an underscore. It compiled
    # anyway, which is worse than failing -- the escape survived into a
    # filename and only graphicx's own expansion rescued it.
    t = re.sub(r'\\bookgraphic\{[^}]+\}', unescape_braces, t)
    return t

body_tex = escape_prose_underscores(body_tex)
abstract_tex = escape_prose_underscores(abstract_tex)

# The tables the manuscript sets are data, not prose: smaller, and in the
# book's table face. The command is defined in the style file.
body_tex = body_tex.replace('\\begin{tabular}', '\\booktablesize\\begin{tabular}')

howto_tex = escape_prose_underscores(cvt(howto)) if howto else ''
# Front matter sits inside an unnumbered chapter, so its own headings are
# unnumbered too -- numbered they would read 0.1, 0.2, which is what a section
# counter does before the first chapter starts.
howto_tex = howto_tex.replace('\\section{', '\\section*{')
howto_tex = convert_callouts(howto_tex)
howto_tex = convert_blockquotes(howto_tex)
howto_tex = convert_lists(howto_tex)
howto_tex = re.sub(r'\[\^(\w+)\]', repl_footnote, howto_tex)

# Front matter is arranged by the style file and WORDED here, from the
# manuscript. The colophon is the one place the pipeline describes itself,
# because a reader holding a generated book should be told it is generated.
colophon = (
    'This edition was set from \\bookcode{article/drafts/v1.md} by '
    '\\bookcode{scripts/generate\\_latex.py} and the book design in '
    '\\bookcode{article/drafts/bookdesign.sty}. Every figure in it is drawn by '
    'a committed script from committed data; every number it reports can be '
    'recomputed from the repository it comes with.\n\n'
    'The work is open. Source, data, simulation code and the issue history are '
    'at \\bookcode{github.com/ELares/cancer\\_research}.\n\n'
    'This is a research document. It reports simulations and literature '
    'measurements, not clinical evidence, and nothing in it is medical '
    'advice.\n\n'
    '\\textit{Keywords:} ' + cvt(keywords)
)

latex = f"""\\documentclass[11pt,twoside,openright]{{report}}
% The design, and the ONLY package this document loads: everything about how
% the book looks -- including which packages that requires -- is in that one
% file, and this document declares structure and lets the style draw it.
\\usepackage[notes=end]{{bookdesign}}
% FLAT, PINNED FIGURE NUMBERS. `report.cls` sets \\thefigure to
% \\thechapter.\\arabic{{figure}}, so the per-float \\setcounter produced 6.25
% where the prose said 25, and NOT ONE of the 29 citations resolved -- the
% deliverable was true of the markdown and false of the compiled PDF.
\\renewcommand{{\\thefigure}}{{\\arabic{{figure}}}}

\\hypersetup{{pdftitle={{{full_title}}}, pdfauthor={{Ezequiel Lares}}}}

\\begin{{document}}

\\frontmatter@bookish
\\bookhalftitle{{{title}}}
\\booktitlepage{{{title}}}{{{subtitle}}}{{Ezequiel Lares}}{{Independent researcher}}
\\bookcolophon{{{colophon}}}

\\pagenumbering{{roman}}

\\chapter*{{In Brief}}
\\markboth{{In Brief}}{{In Brief}}
\\addcontentsline{{toc}}{{chapter}}{{In Brief}}
{abstract_tex}

\\cleardoublepage
\\tableofcontents

{('\\chapter*{A Note On This Edition}' + chr(10) +
  '\\markboth{A Note On This Edition}{A Note On This Edition}' + chr(10) +
  '\\addcontentsline{toc}{chapter}{A Note On This Edition}' + chr(10) +
  howto_tex) if howto_tex else ''}

\\cleardoublepage
\\pagenumbering{{arabic}}

{body_tex}

\\printbooknotes

% Reference list kept as appendix in v1.md; citations are inline notes.

\\end{{document}}
"""

# `\frontmatter` belongs to book.cls and this document is a report; the marker
# above is removed rather than left to fail, so the intent stays readable in
# the source of this generator without emitting a command that does not exist.
latex = latex.replace('\\frontmatter@bookish\n', '')

latex = restore_code_blocks(latex)

TEX.write_text(latex)
print(f'Written {TEX}: {len(latex)} chars, {latex.count(chr(92)+"cite{")} citations')
