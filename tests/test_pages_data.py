"""The public snapshot preserves source results and refuses incomplete data."""

import contextlib
import copy
import hashlib
import io
import json
from pathlib import Path
import shutil
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import build_pages_data as pages  # noqa: E402


def source(path):
    return json.loads((ROOT / path).read_text(encoding="utf-8"))


@contextlib.contextmanager
def copied_sources():
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        for path in pages.SOURCE_PATHS:
            target = root / path
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(ROOT / path, target)
        yield root


def quiet_main(arguments):
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
        return pages.main(arguments)


class PagesDataTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.snapshot = pages.build_snapshot(ROOT)

    def test_committed_snapshot_is_fresh_with_exact_source_provenance(self):
        self.assertEqual((ROOT / pages.OUTPUT).read_bytes(), pages.encode_snapshot(self.snapshot))
        self.assertEqual(self.snapshot["schema_version"], 1)
        self.assertEqual(self.snapshot["snapshot_date"], "2026-09-21")
        self.assertEqual({entry["path"] for entry in self.snapshot["sources"]}, set(pages.SOURCE_PATHS))
        for entry in self.snapshot["sources"]:
            with self.subTest(path=entry["path"]):
                self.assertEqual(entry["sha256"], hashlib.sha256((ROOT / entry["path"]).read_bytes()).hexdigest())
                self.assertEqual(entry["url"], pages.REPOSITORY + entry["path"])

    def test_census_keeps_both_denominators_and_separate_archive(self):
        design = source("analysis/census-evidence-design.json")
        atlas = source("analysis/atlas-site-coverage.json")
        census = self.snapshot["census"]
        self.assertEqual(census["indexed"], design["census"])
        self.assertEqual(census["trials"], design["classes"]["trial"])
        self.assertEqual(census["classifiable"], design["classifiable"])
        self.assertEqual(census["undetermined"], design["classes"]["undetermined"])
        self.assertEqual(census["share_of_census"], round(100 * census["trials"] / census["indexed"], 2))
        self.assertEqual(census["share_of_classifiable"], round(100 * census["trials"] / census["classifiable"], 2))
        self.assertEqual(census["recovered"], atlas["excluded_streams"]["text_matched_no_mesh"])
        self.assertEqual(census["collection_total"], census["indexed"] + census["recovered"])
        archive = [json.loads(line) for line in (ROOT / "corpus/INDEX.jsonl").read_text().splitlines() if line.strip()]
        self.assertEqual(census["archive_records"], len(archive))
        self.assertLess(census["archive_records"], census["indexed"])

    def test_mechanisms_preserve_counts_rankings_years_and_caveats(self):
        profile = source("analysis/census-mechanism-profile.json")
        actual = {row["id"]: row for row in self.snapshot["mechanisms"]}
        self.assertEqual(set(actual), set(profile["count"]))
        for row in profile["rows"]:
            key = row["mechanism"]
            with self.subTest(mechanism=key):
                view = actual[key]
                self.assertEqual([view["articles"], view["trials"], view["trial_share"], view["growth"]],
                                 [row["census"], row["trials"], row["trial_share"], row["growth"]])
                self.assertEqual(view["years"], sorted([int(year), n] for year, n in profile["by_year"][key].items()))
                self.assertEqual([(item["id"], item["articles"], item["enrichment"]) for item in view["top_sites"]],
                                 [(item["site"], item["n"], item["enrichment"]) for item in row["top_sites"]])
                self.assertEqual([(item["id"], item["articles"]) for item in view["top_partners"]],
                                 [(item["mechanism"], item["n"]) for item in row["top_partners"]])
        self.assertEqual([row["trial_share"] for row in self.snapshot["mechanisms"]],
                         sorted((row["trial_share"] for row in actual.values()), reverse=True))
        self.assertIn("broad descriptor", actual["sonodynamic"]["label"])
        self.assertIn("broader than sonodynamic therapy", actual["sonodynamic"]["note"])
        self.assertEqual(actual["car-t"]["label"], "CAR T-cell / adoptive therapy")
        self.assertIn("Immunotherapy, Adoptive", actual["car-t"]["note"])
        self.assertIn("beyond CAR T-cell therapy", actual["car-t"]["note"])
        self.assertEqual(actual["epigenetic"]["label"], "Epigenetic processes & therapies")
        self.assertIn("DNA Methylation", actual["epigenetic"]["note"])
        self.assertIn("basic research", actual["epigenetic"]["note"])
        self.assertIsNone(actual["mrna-vaccine"]["growth"])
        self.assertEqual({row["id"]: row["articles"] for row in self.snapshot["sites"]}, profile["site_totals"])

    def test_sampling_and_challenge_keep_successes_and_failures(self):
        sampling = source("analysis/calibration/joint-resample-sampling.json")
        challenge = source("analysis/calibration/proposal-coverage-challenges.json")
        self.assertEqual(self.snapshot["sampling"]["passed"], sampling["stability"]["passed"])
        for actual, expected in zip(self.snapshot["sampling"]["runs"], sampling["runs"]):
            self.assertEqual(actual["seed"], expected["seed"])
            for view, original in (("accepted", "n_accepted"), ("attempts", "n_attempts"), ("ess", "ess")):
                self.assertEqual(actual[view], expected["importance"][original])
        summary = self.snapshot["challenge"]
        self.assertEqual(summary["passed"], challenge["passed"])
        self.assertEqual(summary["passed_runs"], sum(run["assessment"]["passed"] for run in challenge["positive_runs"]))
        self.assertEqual(summary["total_runs"], len(challenge["positive_runs"]))
        self.assertEqual((summary["passed_runs"], summary["total_runs"]), (7, 9))
        self.assertEqual(summary["thresholds"], challenge["specification"]["gates"])
        actual_runs = {run["seed"]: run for group in summary["groups"] for run in group["runs"]}
        for original in challenge["positive_runs"]:
            assessment = original["assessment"]
            actual = actual_runs[original["seed"]]
            self.assertEqual(actual["passed"], assessment["passed"])
            self.assertEqual(actual["region_error"], max(assessment["absolute_region_mass_errors"].values()))
            self.assertEqual(actual["max_weight"], assessment["importance"]["max_normalized_weight"])
            expected_failed = sorted(key for group in ("usual_checks", "truth_checks")
                                     for key, passed in assessment[group].items() if not passed)
            self.assertEqual(actual["failed_checks"], expected_failed)
        self.assertEqual(summary["controls"]["oracle_passed"], len(challenge["oracle_controls"]))
        self.assertEqual(summary["controls"]["negative_failed_truth"], len(challenge["negative_controls"]))
        self.assertEqual(summary["controls"]["negative_passed_usual"], len(challenge["negative_controls"]))

    def test_incomplete_and_contradictory_reports_fail_closed(self):
        profile = source("analysis/census-mechanism-profile.json")
        challenge = source("analysis/calibration/proposal-coverage-challenges.json")
        bad_profiles = []
        bad = copy.deepcopy(profile)
        bad["rows"].pop()
        bad_profiles.append(bad)
        bad = copy.deepcopy(profile)
        bad["by_year"]["car-t"].pop("2025")
        bad_profiles.append(bad)
        bad = copy.deepcopy(profile)
        bad["rows"][0]["trial_share"] = 99
        bad_profiles.append(bad)
        for index, bad in enumerate(bad_profiles):
            with self.subTest(profile=index), self.assertRaises(ValueError):
                pages.mechanism_summary(bad)
        for missing in ("positive_runs", "negative_controls", "oracle_controls"):
            bad = copy.deepcopy(challenge)
            bad[missing].pop()
            with self.subTest(missing=missing), self.assertRaises(ValueError):
                pages.challenge_summary(bad)
        bad = copy.deepcopy(challenge)
        bad["passed"] = True
        with self.assertRaises(ValueError):
            pages.challenge_summary(bad)

    def test_invalid_types_ranges_and_nonfinite_values_are_rejected(self):
        sampling = source("analysis/calibration/joint-resample-sampling.json")
        for field, value in (("n_accepted", True), ("n_accepted", -1), ("n_attempts", 1),
                             ("ess", float("nan")), ("ess", float("inf")), ("ess", 9000)):
            bad = copy.deepcopy(sampling)
            bad["runs"][0]["importance"][field] = value
            with self.subTest(field=field, value=value), self.assertRaises(ValueError):
                pages.sampling_summary(bad)
        bad = copy.deepcopy(sampling)
        bad["stability"]["passed"] = "false"
        with self.assertRaises(ValueError):
            pages.sampling_summary(bad)
        for raw in ('{"value": NaN}', '{"value": Infinity}', '{"nested": [1e999]}', '{"value": 1, "value": 2}'):
            with self.subTest(raw=raw), self.assertRaises(ValueError):
                pages.parse_json(raw)

    def test_check_detects_missing_stale_and_changed_sources_without_writing(self):
        with copied_sources() as root:
            output = root / pages.OUTPUT
            self.assertEqual(quiet_main(["--root", str(root), "--check"]), 1)
            self.assertFalse(output.parent.exists())
            self.assertEqual(quiet_main(["--root", str(root)]), 0)
            original = output.read_bytes()
            self.assertEqual(quiet_main(["--root", str(root), "--check"]), 0)
            input_path = root / pages.SOURCE_PATHS[0]
            input_path.write_bytes(input_path.read_bytes() + b"\n")
            self.assertEqual(quiet_main(["--root", str(root), "--check"]), 1)
            self.assertEqual(output.read_bytes(), original)
            self.assertEqual(quiet_main(["--root", str(root)]), 0)
            self.assertNotEqual(output.read_bytes(), original)
            self.assertEqual(quiet_main(["--root", str(root), "--check"]), 0)

    def test_invalid_inputs_never_replace_or_create_an_output(self):
        with copied_sources() as root:
            output = root / pages.OUTPUT
            bad_source = root / pages.SOURCE_PATHS[0]
            bad_source.write_text("{}", encoding="utf-8")
            self.assertEqual(quiet_main(["--root", str(root)]), 1)
            self.assertFalse(output.parent.exists())
            output.parent.mkdir(parents=True)
            output.write_bytes(b"keep existing output\n")
            self.assertEqual(quiet_main(["--root", str(root)]), 1)
            self.assertEqual(output.read_bytes(), b"keep existing output\n")
            self.assertEqual(list(output.parent.iterdir()), [output])

    def test_archive_identity_and_cross_source_census_are_validated(self):
        with copied_sources() as root:
            archive = root / "corpus/INDEX.jsonl"
            original = archive.read_bytes()
            archive.write_bytes(original + original.splitlines()[0] + b"\n")
            with self.assertRaisesRegex(ValueError, "duplicate archive"):
                pages.build_snapshot(root)
            archive.write_bytes(original)
            atlas_path = root / "analysis/atlas-site-coverage.json"
            atlas = json.loads(atlas_path.read_text())
            atlas["census"] -= 1
            atlas_path.write_text(json.dumps(atlas), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "census sources disagree"):
                pages.build_snapshot(root)


if __name__ == "__main__":
    unittest.main()
