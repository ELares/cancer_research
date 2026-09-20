"""Public-site contracts: working links, bounded loading and honest study scope."""

from collections import Counter
from html.parser import HTMLParser
import json
from pathlib import Path
import unittest
from urllib.parse import unquote, urlparse

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"


class Document(HTMLParser):
    def __init__(self, path):
        super().__init__()
        self.path = path
        self.tags = []
        self.feed(path.read_text(encoding="utf-8"))

    def handle_starttag(self, tag, attrs):
        self.tags.append((tag, dict(attrs)))

    def ids(self):
        return [attrs["id"] for _, attrs in self.tags if "id" in attrs]


class PagesSiteTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pages = [Document(DOCS / name) for name in ("index.html", "dashboard.html")]
        cls.snapshot = json.loads((DOCS / "assets/research-data.json").read_text())

    def test_internal_assets_fragments_and_repository_source_links_exist(self):
        for page in self.pages:
            ids = page.ids()
            self.assertEqual(len(ids), len(set(ids)), f"Duplicate IDs in {page.path.name}")
            for tag, attrs in page.tags:
                for key in ("href", "src"):
                    value = attrs.get(key)
                    if not value:
                        continue
                    with self.subTest(page=page.path.name, url=value):
                        url = urlparse(value)
                        if url.netloc == "github.com" and url.path.startswith("/ELares/cancer_research/blob/main/"):
                            path = unquote(url.path.removeprefix("/ELares/cancer_research/blob/main/"))
                            self.assertTrue((ROOT / path).is_file(), path)
                        elif not url.scheme and not url.netloc:
                            if not url.path and url.fragment:
                                self.assertIn(unquote(url.fragment), ids)
                            elif url.path:
                                self.assertFalse(url.path.startswith("/"), "Project Pages needs relative asset paths")
                                self.assertTrue((page.path.parent / unquote(url.path)).exists(), value)

    def test_each_page_has_one_main_heading_and_labelled_inputs(self):
        for page in self.pages:
            with self.subTest(page=page.path.name):
                self.assertEqual(sum(tag == "h1" for tag, _ in page.tags), 1)
                self.assertEqual(sum(tag == "main" for tag, _ in page.tags), 1)
                labels = {attrs["for"] for tag, attrs in page.tags if tag == "label" and "for" in attrs}
                for tag, attrs in page.tags:
                    if tag in ("input", "select"):
                        self.assertTrue(attrs.get("id") in labels or attrs.get("aria-label"))
                self.assertTrue(any(tag == "html" and attrs.get("lang") == "en" for tag, attrs in page.tags))

    def test_first_render_is_local_and_the_atlas_stays_small(self):
        resources = {DOCS / "index.html", DOCS / "assets/research-data.json"}
        for page in self.pages:
            for tag, attrs in page.tags:
                if tag not in ("script", "link", "img", "iframe"):
                    continue
                target = attrs.get("src") or attrs.get("href")
                if target:
                    self.assertFalse(urlparse(target).scheme or target.startswith("//"), target)
                    if page.path.name == "index.html":
                        resources.add(DOCS / target)
        self.assertLess(sum(path.stat().st_size for path in resources), 200_000,
                        "Review the initial-load budget before adding large payloads")
        css = (DOCS / "assets/site.css").read_text()
        self.assertNotIn("@import", css)
        self.assertNotIn("https://", css)
        self.assertTrue((DOCS / ".nojekyll").is_file())

    def test_curated_study_text_still_matches_the_selected_historical_results(self):
        """A regenerated count must not silently contradict the adjacent prose.

        Deliberately pin this *presentation's* selected outcomes. Selecting a
        later study requires reviewing its interpretation, not rewriting the old
        scientific result or merely replacing a visible number.
        """
        sampling, challenge = self.snapshot["sampling"], self.snapshot["challenge"]
        self.assertTrue(sampling["passed"])
        self.assertEqual([run["passed"] for run in sampling["runs"]], [True] * 3)
        self.assertFalse(challenge["passed"])
        self.assertEqual((challenge["passed_runs"], challenge["total_runs"]), (7, 9))
        failures = {group["id"]: Counter(check for run in group["runs"] for check in run["failed_checks"])
                    for group in challenge["groups"]}
        self.assertEqual(failures, {"rotated_box": Counter({"maximum_weight": 2, "region_masses": 1}),
                                    "annular_cylinder": Counter(), "unequal_balls": Counter()})
        controls = challenge["controls"]
        self.assertEqual(controls["oracle_passed"], controls["oracle_total"])
        self.assertEqual(controls["negative_failed_truth"], controls["negative_total"])
        self.assertEqual(controls["negative_passed_usual"], controls["negative_total"])
        source_paths = {source["path"] for source in self.snapshot["sources"]}
        self.assertIn("analysis/calibration/proposal-coverage-challenges.json", source_paths)
        html = (DOCS / "index.html").read_text()
        self.assertIn("Overall study: failed", html)
        self.assertIn("Independent validation pending", html)
        self.assertIn("not how effective a treatment is", html)

    def test_census_archive_and_trial_publications_remain_distinct(self):
        census = self.snapshot["census"]
        self.assertEqual(census["archive_records"], 4830)
        self.assertGreater(census["collection_total"], census["indexed"])
        for page in self.pages:
            html = page.path.read_text().lower()
            self.assertIn("record-level browsing of the census is not offered", html)
            self.assertIn("4,830", html)
        js = (DOCS / "assets/site.js").read_text()
        self.assertIn("trial-labelled records", js)
        self.assertNotIn("${number(row.trials)} trials", js)
        self.assertIn("data-stat=\"trial-classifiable-share\"", (DOCS / "index.html").read_text())

    def test_the_pages_workflow_covers_all_snapshot_inputs(self):
        workflow = (ROOT / ".github/workflows/pages-check.yml").read_text()
        self.assertIn("build_pages_data.py --check", workflow)
        self.assertIn("test_pages_*.py", workflow)
        # Both PR and main-push events must cover each input family, including
        # future changes to a source that don't touch the site itself.
        for pattern in ("'docs/**'", "'analysis/**'", "'corpus/INDEX.jsonl'",
                        "'scripts/build_pages_data.py'"):
            self.assertEqual(workflow.count(pattern), 2, pattern)


if __name__ == "__main__":
    unittest.main()
