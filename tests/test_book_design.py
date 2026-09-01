"""The book design stays in one file, and the pipeline stays out of it.

WHY THIS FILE EXISTS
--------------------
The manuscript was reset from a LaTeX report into a 6x9 trade book. A design is
only worth having if the NEXT chapter inherits it, and the way that stops being
true is always the same: one hand-placed `\\vspace` to fix one page, then a
`\\fontsize` for one caption, and six months later the layout lives in forty
places and nobody can change the trim size.

So the split is a rule, not an intention. `article/drafts/bookdesign.sty` holds
every visual decision. `scripts/generate_latex.py` emits meaning -- \\chapter,
\\standfirst, \\begin{finding} -- and is forbidden from naming a length, a
size, a colour or a font. These guards are what makes "the design is
adaptable" a checkable claim rather than a hope.

WHAT THEY DO NOT CHECK
----------------------
Whether the book looks good. Nothing here can tell you that, and a green suite
is not a substitute for compiling it and turning the pages.
"""
import re
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
STY = REPO / "article" / "drafts" / "bookdesign.sty"
GEN = REPO / "scripts" / "generate_latex.py"
TEX = REPO / "article" / "drafts" / "v1.tex"
MD = REPO / "article" / "drafts" / "v1.md"
AUTHORING = REPO / "article" / "AUTHORING.md"
WORKFLOW = REPO / ".github" / "workflows" / "release-pdf.yml"

# The release build installs exactly these Debian texlive sets. Every package
# the style file loads was checked against their file lists before it was used;
# a package outside them compiles here and fails there, which is the worst
# place to find out.
ALLOWED_PACKAGES = {
    # texlive-latex-base
    "geometry", "fontenc", "inputenc", "graphicx", "fancyhdr", "hyperref",
    "mathpazo", "avant", "pifont", "amsmath", "amssymb",
    # texlive-latex-recommended
    "xcolor", "microtype", "letterspace", "ragged2e", "setspace", "booktabs",
    "caption", "etoolbox",
    # texlive-latex-extra
    "titlesec", "titletoc", "lettrine", "enumitem", "needspace", "endnotes",
}

# Layout primitives. A generator that emits any of these has started making
# design decisions, which is the failure this file exists to catch.
FORBIDDEN_IN_GENERATOR = [
    (r"\\geometry\b", "page geometry"),
    (r"\\fontsize\b", "a type size"),
    (r"\\linespread\b", "leading"),
    (r"\\vspace\b", "vertical space"),
    (r"\\hspace\b", "horizontal space"),
    (r"\\textwidth\b", "the measure"),
    (r"\\centering\b", "alignment"),
    (r"\\rule\{", "a rule"),
    (r"\\color\b", "a colour"),
    (r"\\bfseries\b|\\itshape\b|\\small\b|\\footnotesize\b", "a font switch"),
    (r"\\usepackage\{(?!bookdesign)", "a package other than the design"),
    (r"documentclass\[[^\]]*a4paper", "a paper size"),
]


def test_the_generator_emits_no_layout():
    src = GEN.read_text()
    # The prose that explains the rule necessarily names the commands it
    # forbids, so only CODE is scanned: the module docstring and comment lines
    # are removed first.
    src = re.sub(r'r?""".*?"""', "", src, count=1, flags=re.DOTALL)
    code = "\n".join(l for l in src.split("\n") if not l.lstrip().startswith("#"))
    found = [(pat, what) for pat, what in FORBIDDEN_IN_GENERATOR
             if re.search(pat, code)]
    assert not found, (
        "scripts/generate_latex.py has started making layout decisions: "
        + "; ".join(f"{what} ({pat})" for pat, what in found)
        + ". Move it into article/drafts/bookdesign.sty, which is the only "
        "place the book's design is allowed to live.")


def test_every_package_the_design_loads_is_in_the_release_build():
    loaded = set(re.findall(r"\\RequirePackage(?:\[[^\]]*\])?\{([^}]+)\}",
                            STY.read_text()))
    names = {n.strip() for group in loaded for n in group.split(",")}
    unknown = names - ALLOWED_PACKAGES
    assert not unknown, (
        f"bookdesign.sty loads {sorted(unknown)}, which are not recorded as "
        "present in the release build's texlive package set "
        "(.github/workflows/release-pdf.yml). Check the package is in one of "
        "those Debian sets and add it to ALLOWED_PACKAGES, or install the set "
        "that carries it in the workflow -- do not find out from a red build.")


def test_the_style_file_is_what_the_document_loads():
    tex = TEX.read_text()
    assert re.search(r"\\usepackage(\[[^\]]*\])?\{bookdesign\}", tex), (
        "v1.tex does not load bookdesign.sty; the book has no design")
    assert STY.exists()
    # It has to sit beside v1.tex, because that is the directory pdflatex runs
    # in and a .sty anywhere else is simply not found.
    assert STY.parent == TEX.parent, (
        f"bookdesign.sty is in {STY.parent} and pdflatex runs in {TEX.parent}")


def test_the_document_is_a_book_and_not_a_report():
    """The three decisions that make it one, pinned so a regeneration cannot
    quietly undo them."""
    sty = STY.read_text()
    assert "paperwidth=6in" in sty and "paperheight=9in" in sty, (
        "the trim size is no longer 6x9; if that was deliberate, change this "
        "guard deliberately too")
    assert "\\RequirePackage[sc]{mathpazo}" in sty, "the body face has changed"
    tex = TEX.read_text()
    assert "twoside" in tex, "a book is set two-sided"


def test_every_callout_the_manuscript_uses_has_an_environment():
    """A fence name with no environment is a build error, and a fence the
    style file defines but nobody uses is dead design."""
    md = MD.read_text()
    used = set(re.findall(r"^::: *(\S+)", md, flags=re.MULTILINE))
    defined = set(re.findall(r"\\newenvironment\{(finding|refusal|numbers)\}",
                             STY.read_text()))
    missing = used - defined
    assert not missing, (
        f"v1.md uses callouts {sorted(missing)} that bookdesign.sty does not "
        "define")
    # Documented, because an author cannot use a construct they cannot find.
    doc = AUTHORING.read_text()
    for name in defined:
        assert f"::: {name}" in doc, (
            f"the {name} callout is not documented in article/AUTHORING.md")


def test_no_fenced_line_is_wider_than_the_measure():
    """A verbatim line cannot wrap: a long one runs off the page.

    This is the one authoring rule the design cannot enforce for itself, and
    it shipped broken -- a command line of 120 characters ran into the margin
    and past the trim.
    """
    LIMIT = 78
    long_lines = []
    for block in re.findall(r"```[a-z]*\n(.*?)```", MD.read_text(), re.DOTALL):
        for line in block.split("\n"):
            if len(line) > LIMIT:
                long_lines.append(f"{len(line)} chars: {line[:60]}...")
    assert not long_lines, (
        f"fenced lines longer than {LIMIT} characters run off the 6x9 page:\n  "
        + "\n  ".join(long_lines))


def test_the_standfirsts_are_attached_to_headings():
    """A blockquote under a heading is a standfirst; anywhere else it is a
    pull quote. Both are legitimate, so this checks the ones that exist are
    where the author meant them -- every part and chapter has one, which is
    what makes the openers work."""
    md = MD.read_text()
    headings = re.findall(r"^(# Part [IVX]+: .+|## Chapter \d+: .+)$", md,
                          flags=re.MULTILINE)
    missing = []
    for h in headings:
        after = md.split(h, 1)[1].lstrip("\n")
        if not after.startswith("> "):
            missing.append(h)
    assert not missing, (
        "these openers have no standfirst, so the page under the title is "
        f"empty: {missing}. Add one as a blockquote directly under the "
        "heading (see article/AUTHORING.md).")


def test_the_workflow_still_installs_what_the_design_needs():
    """The design's dependencies and the build's package list are two files
    that have to agree; this is the pair that fails silently."""
    wf = WORKFLOW.read_text()
    for pkg in ("texlive-latex-base", "texlive-latex-recommended",
                "texlive-latex-extra", "texlive-fonts-recommended"):
        assert pkg in wf, (
            f"{pkg} is no longer installed by release-pdf.yml, and "
            "bookdesign.sty depends on it")
