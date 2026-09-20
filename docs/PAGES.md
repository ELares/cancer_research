# GitHub Pages research atlas

The public site at <https://elares.github.io/cancer_research/> is served directly
from `main:/docs`. The front page is plain HTML, CSS and JavaScript: there is no
build framework, Python download, third-party font, tracking script or backend
required to explore it. `.nojekyll` keeps these files as ordinary static assets.

## What is available

- **Evidence atlas:** searchable mechanisms, sorting, publication trends with
  accessible count tables, and cancer-site enrichment. Charts cover the indexed
  census, not the larger indexed-plus-text-recovered collection. Trial shares
  describe publication labels, not treatment efficacy. Descriptor breadth and
  overlapping tags limit comparisons.
- **Model checks:** the biological-model sampling screens, the failed additional
  geometry challenge, and the pending independent-assay validation are distinct.
  The geometry switcher shows each run's actual results; its diagrams are
  explicitly two-dimensional schematics of seven-dimensional targets.
- **Roadmap and resources:** links to current methods, failed results, missing
  inputs, the manuscript, the model card and contribution instructions.
- **Optional Python dashboard:** `dashboard.html` loads stlite only after the
  visitor clicks Launch. It provides the historical 4,830-record archive,
  census aggregates and a historical read-only simulation report, explicitly
  distinguished from current calibration studies. It requires a
  network connection for its pinned stlite runtime and repository files; a
  first launch may take 30–60 seconds. It does not run the compiled simulator.

Navigation and source links work without JavaScript. Interactive panels explain
when their local snapshot is unavailable and provide a retry action. The site
supports keyboard operation, small screens and reduced-motion preferences.

## Update the snapshot

The small `docs/assets/research-data.json` is generated from six committed
sources. `scripts/build_pages_data.py` records their exact paths and SHA256
hashes, validates denominators and run completeness, and writes atomically.
It does not run any simulations or change frozen study artifacts.

```bash
python3 scripts/build_pages_data.py
python3 scripts/build_pages_data.py --check
```

`--check` fails for a missing or stale snapshot and never writes. Source changes
must include a regenerated snapshot. When adding a new study, explicitly select
its source, review its interpretation, update the narrative in `index.html` and
`site.js`, and update `SNAPSHOT_DATE` to the date of that content review. The date
is a reviewed research snapshot, not the visitor's date or a claim of live data.
Do not rewrite previous scientific results to refresh the site.

The headline collection count, indexed denominator, mechanism rows and run
metrics come from the snapshot. The project interpretation is curated prose,
so the Pages contract tests also pin the cited studies and outcomes. Review
both when changing which study the page presents.

## Preview and verify

```bash
python3 -m http.server 8000 --directory docs
# Open http://localhost:8000/
python3 -m unittest discover -s tests -p 'test_pages_*.py' -v
node --check docs/assets/site.js
node --check docs/assets/dashboard.js
node --test tests/test_pages_*.cjs
```

Use an HTTP server, not `file://`, because browsers restrict local JSON fetches.
All asset links are relative and work under the `/cancer_research/` project path.
For a matching subpath preview, serve a parent directory containing a
`cancer_research` symlink to `docs`.

For a browser review, check desktop, tablet, 390px and 320px layouts; keyboard
focus and Escape in the sources dialog; search/no-results/clear; sorting and
mechanism selection (including the last row after expanding the list); geometry
switching; JSON failure/retry; and no-JavaScript
navigation. Run an accessibility audit in the main and dialog states. Check that
no external resource is downloaded from the atlas or before launching the
optional dashboard. The dashboard uses files from `main`; to test a pending
Python change, intercept those raw GitHub requests with the local files in a
browser harness or use a temporary test branch URL, without committing that URL.
Check source-download HTTP errors and network failures, then retry: incomplete
files must not launch a partial dashboard or remain labelled as loading.

The dedicated Pages workflow checks the snapshot, site contracts, optional
dashboard file mounts, JavaScript syntax, and download/retry behavior through
Node's built-in test runner. The regular Python suite includes the Python tests
too. After edits, refresh
`MANIFEST.sha256` using the instructions in `CONTRIBUTING.md`. A PR changes the
previewable files; the public site updates only after merge into `main`.
