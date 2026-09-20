#!/usr/bin/env python3
"""Build the small Pages snapshot from committed reports; never run analyses.

Percent shares use 0–100 units; growth is the source's 2015–2025 multiplier
(null where its baseline is too small). Source hashes cover the exact bytes.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = Path("docs/assets/research-data.json")
SNAPSHOT_DATE = "2026-09-20"
REPOSITORY = "https://github.com/ELares/cancer_research/blob/main/"
SOURCE_PATHS = (
    "analysis/census-evidence-design.json",
    "analysis/census-mechanism-profile.json",
    "analysis/atlas-site-coverage.json",
    "analysis/calibration/joint-resample-sampling.json",
    "analysis/calibration/proposal-coverage-challenges.json",
    "corpus/INDEX.jsonl",
)
MECHANISM_LABELS = {
    "antibody-drug-conjugate": "Antibody–drug conjugates",
    "bispecific-antibody": "Bispecific antibodies",
    "car-t": "CAR T-cell / adoptive therapy",
    "crispr": "CRISPR gene editing",
    "electrochemical-therapy": "Electrochemical therapy",
    "epigenetic": "Epigenetic processes & therapies",
    "hifu": "High-intensity focused ultrasound",
    "immunotherapy": "Immunotherapy",
    "metabolic-targeting": "Metabolic targeting",
    "microbiome": "Microbiome",
    "mrna-vaccine": "mRNA vaccines",
    "nanoparticle": "Nanoparticles",
    "oncolytic-virus": "Oncolytic viruses",
    "phagocytosis-checkpoint": "Phagocytosis checkpoints",
    "sonodynamic": "Therapeutic ultrasound (broad descriptor)",
    "synthetic-lethality": "Synthetic lethality",
}
SITE_LABELS = {
    "brain/CNS": "Brain / central nervous system",
    "cervix/uterus": "Cervix / uterus",
    "head and neck": "Head and neck",
    "leukaemia": "Leukemia",
    "oesophagus": "Esophagus",
    "skin/melanoma": "Skin / melanoma",
}
CHALLENGE_LABELS = {
    "rotated_box": "Rotated box",
    "annular_cylinder": "Annular cylinder",
    "unequal_balls": "Unequal balls",
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def count(value, name, minimum=0):
    require(type(value) is int and value >= minimum, f"{name}: expected integer >= {minimum}")
    return value


def number(value, name, minimum=0, maximum=None):
    require(type(value) in (int, float) and math.isfinite(value), f"{name}: expected finite number")
    require(value >= minimum and (maximum is None or value <= maximum), f"{name}: out of range")
    return value


def boolean(value, name):
    require(type(value) is bool, f"{name}: expected boolean")
    return value


def mapping(value, name):
    require(isinstance(value, dict) and bool(value), f"{name}: expected nonempty object")
    return value


def sequence(value, name):
    require(isinstance(value, list) and bool(value), f"{name}: expected nonempty array")
    return value


def checks(value, name):
    return {key: boolean(item, f"{name}.{key}") for key, item in mapping(value, name).items()}


def unique_object(pairs):
    result = {}
    for key, value in pairs:
        require(key not in result, f"duplicate JSON key: {key}")
        result[key] = value
    return result


def reject_nonfinite(value):
    if isinstance(value, dict):
        for item in value.values():
            reject_nonfinite(item)
    elif isinstance(value, list):
        for item in value:
            reject_nonfinite(item)
    elif isinstance(value, float):
        require(math.isfinite(value), "nonfinite JSON number")


def parse_json(raw):
    value = json.loads(raw, object_pairs_hook=unique_object)
    reject_nonfinite(value)
    return mapping(value, "source")


def importance_summary(run, name):
    importance = mapping(run["importance"], f"{name}.importance")
    accepted = count(importance["n_accepted"], f"{name}.accepted")
    attempts = count(importance["n_attempts"], f"{name}.attempts", 1)
    require(accepted <= attempts, f"{name}: accepted exceeds attempts")
    require(boolean(importance["assessable"], name), f"{name}: unassessable importance summary")
    ess = number(importance["ess"], f"{name}.ess", maximum=accepted + 1e-8)
    return {"accepted": accepted, "attempts": attempts, "ess": ess}


def census_summary(design, profile, atlas, archive_records):
    indexed = count(design["census"], "census", 1)
    classifiable = count(design["classifiable"], "classifiable", 1)
    classes = {key: count(value, f"classes.{key}")
               for key, value in mapping(design["classes"], "classes").items()}
    trials, undetermined = classes["trial"], classes["undetermined"]
    require(sum(classes.values()) == indexed, "evidence classes do not sum to census")
    require(classifiable + undetermined == indexed and trials <= classifiable,
            "inconsistent evidence denominators")
    require(count(profile["census"], "profile.census") == indexed
            and count(atlas["census"], "atlas.census") == indexed, "census sources disagree")
    recovered = count(atlas["excluded_streams"]["text_matched_no_mesh"], "recovered")
    return {"indexed": indexed, "trials": trials, "classifiable": classifiable,
            "undetermined": undetermined, "share_of_census": round(100 * trials / indexed, 2),
            "share_of_classifiable": round(100 * trials / classifiable, 2),
            "archive_records": archive_records, "recovered": recovered,
            "collection_total": indexed + recovered}


def mechanism_summary(profile):
    counts = mapping(profile["count"], "mechanism counts")
    require(set(counts) == set(MECHANISM_LABELS), "mechanisms missing or unknown; update labels explicitly")
    rows = sequence(profile["rows"], "mechanism rows")
    require(len(rows) == len(counts) and {row["mechanism"] for row in rows} == set(counts),
            "mechanism rows missing or duplicated")
    for field in ("trials", "by_year", "by_site", "partners"):
        require(set(mapping(profile[field], field)) == set(counts), f"incomplete {field}")
    site_counts = {site: count(value, f"site_totals.{site}")
                   for site, value in mapping(profile["site_totals"], "site totals").items()}
    require(all(value <= profile["census"] for value in site_counts.values()), "site count exceeds census")
    sites = [{"id": site, "label": SITE_LABELS.get(site, site.capitalize()), "articles": value}
             for site, value in sorted(site_counts.items(), key=lambda item: (-item[1], item[0]))]
    site_labels = {site["id"]: site["label"] for site in sites}
    start_year = count(profile["start_year"], "start_year", 1)
    end_year = count(profile["end_year"], "end_year", start_year + 1)
    output = []
    for row in rows:
        key = row["mechanism"]
        articles = count(row["census"], f"{key}.articles", 1)
        trials = count(row["trials"], f"{key}.trials")
        require(articles == count(counts[key], key) and articles <= profile["census"], f"{key}: invalid article total")
        require(trials == count(profile["trials"][key], key) and trials <= articles, f"{key}: invalid trial total")
        share = number(row["trial_share"], f"{key}.trial_share", maximum=100)
        require(share == round(100 * trials / articles, 2), f"{key}: trial share disagrees with counts")
        years = []
        for year, value in mapping(profile["by_year"][key], f"{key}.years").items():
            require(year.isdigit() and len(year) == 4, f"{key}: invalid year")
            years.append([int(year), count(value, f"{key}.{year}")])
        years.sort()
        require(sum(value for _, value in years) == articles, f"{key}: incomplete annual counts")
        by_year = dict(years)
        start, end = count(row["start"], key), count(row["end"], key)
        require(start == by_year.get(start_year, 0) and end == by_year.get(end_year, 0), f"{key}: growth years disagree")
        growth = row["growth"]
        if growth is not None:
            number(growth, f"{key}.growth")
        require(growth == (round(end / start, 2) if start >= 30 else None), f"{key}: invalid growth")
        top_sites = []
        for site in sequence(row["top_sites"], f"{key}.top_sites"):
            sid = site["site"]
            require(sid in site_counts, f"{key}: unknown site {sid}")
            n = count(site["n"], f"{key}.{sid}")
            require(n == count(profile["by_site"][key][sid], sid)
                    and n <= min(articles, site_counts[sid]), f"{key}: invalid site count")
            top_sites.append({"id": sid, "label": site_labels[sid], "articles": n,
                              "enrichment": number(site["enrichment"], f"{key}.{sid}.enrichment")})
        top_partners = []
        for partner in sequence(row["top_partners"], f"{key}.top_partners"):
            pid = partner["mechanism"]
            require(pid in counts and pid != key, f"{key}: unknown or self partner")
            n = count(partner["n"], f"{key}.{pid}")
            require(n == count(profile["partners"][key][pid], pid)
                    and n <= min(articles, counts[pid]), f"{key}: invalid partner count")
            top_partners.append({"id": pid, "label": MECHANISM_LABELS[pid], "articles": n})
        item = {"id": key, "label": MECHANISM_LABELS[key], "articles": articles,
                "trials": trials, "trial_share": share, "growth": growth,
                "top_sites": top_sites, "top_partners": top_partners, "years": years}
        if key == "sonodynamic":
            item["note"] = ("The source id is sonodynamic, but its Ultrasonic Therapy descriptor "
                            "is broader than sonodynamic therapy. These are broad descriptor counts.")
        elif key == "car-t":
            item["note"] = ("The mapping includes Receptors, Chimeric Antigen and the broader "
                            "Immunotherapy, Adoptive descriptor. Counts include adoptive therapies "
                            "beyond CAR T-cell therapy.")
        elif key == "epigenetic":
            item["note"] = ("The mapping includes therapy descriptors and DNA Methylation, which "
                            "also tags basic methylation biology. Counts include basic research "
                            "as well as epigenetic treatments.")
        output.append(item)
    output.sort(key=lambda row: (-row["trial_share"], row["id"]))
    return output, sites, {"start_year": start_year, "end_year": end_year}


def sampling_summary(source):
    require(type(source["schema_version"]) is int and source["schema_version"] == 2,
            "unsupported sampling schema")
    runs = sequence(source["runs"], "sampling runs")
    require(len(runs) == 3, "sampling requires all three prespecified runs")
    output = []
    for run in runs:
        seed = count(run["seed"], "sampling seed", 1)
        summary = importance_summary(run, f"sampling.{seed}")
        require(summary["attempts"] == count(source["plan"]["production_attempts"], "sampling plan", 1),
                "sampling attempts disagree with plan")
        passed = boolean(run["adequacy"]["passed"], "sampling adequacy")
        require(passed == all(checks(run["adequacy"]["checks"], "adequacy checks").values()),
                "sampling adequacy flag disagrees with checks")
        output.append({"seed": seed, **summary, "passed": passed})
    require(len({run["seed"] for run in output}) == len(output), "duplicate sampling seeds")
    stability = source["stability"]
    passed = boolean(stability["passed"], "sampling passed")
    stability_checks = checks(stability["checks"], "stability checks")
    require(passed == all(stability_checks.values()), "sampling passed disagrees with checks")
    require(stability_checks["all_runs_adequate"] == all(run["passed"] for run in output),
            "sampling run adequacy disagrees with stability")
    thresholds = {key: number(value, f"sampling gate {key}")
                  for key, value in mapping(source["gates"], "sampling gates").items()}
    return {"passed": passed, "runs": sorted(output, key=lambda run: run["seed"]), "thresholds": thresholds}


def challenge_runs(source, field, seeds_field):
    runs = sequence(source[field], field)
    seeds = mapping(source["specification"][seeds_field], seeds_field)
    require(set(seeds) == set(CHALLENGE_LABELS), f"{field}: unexpected or missing fixture")
    expected = set()
    for fixture, values in seeds.items():
        values = sequence(values, f"{seeds_field}.{fixture}")
        require(len(values) == 3 and len(values) == len(set(values)), f"{seeds_field}: missing or duplicate seeds")
        expected.update((fixture, count(value, "seed", 1)) for value in values)
    observed = set()
    output = []
    for run in runs:
        fixture, seed = run["fixture"], count(run["seed"], "challenge seed", 1)
        require((fixture, seed) not in observed, f"{field}: duplicate run")
        observed.add((fixture, seed))
        assessment = run["assessment"]
        summary = importance_summary(assessment, f"{fixture}.{seed}")
        require(summary["attempts"] == count(source["specification"]["production_attempts"], "challenge plan", 1),
                "challenge attempts disagree with plan")
        usual = checks(assessment["usual_checks"], "usual checks")
        truth = checks(assessment["truth_checks"], "truth checks")
        failed = sorted(key for group in (usual, truth) for key, passed in group.items() if not passed)
        if field == "positive_runs" and not boolean(assessment["pilot_reached_epsilon"], "pilot reached"):
            failed.append("pilot_reached_epsilon")
        passed = boolean(assessment["passed"], "challenge passed")
        require(passed == (not failed), "challenge passed disagrees with checks")
        region_errors = mapping(assessment["absolute_region_mass_errors"], "region errors")
        output.append({"fixture": fixture, "seed": seed, **summary,
                       "max_weight": number(assessment["importance"]["max_normalized_weight"], "max weight", maximum=1),
                       "region_error": max(number(value, "region error", maximum=1) for value in region_errors.values()),
                       "passed": passed, "failed_checks": failed,
                       "usual_passed": all(usual.values()), "truth_passed": all(truth.values())})
    require(observed == expected, f"{field}: missing or unexpected prespecified runs")
    return output


def challenge_summary(source):
    require(type(source["schema_version"]) is int and source["schema_version"] == 1,
            "unsupported challenge schema")
    learned = challenge_runs(source, "positive_runs", "learned_seeds")
    oracles = challenge_runs(source, "oracle_controls", "oracle_seeds")
    negatives = challenge_runs(source, "negative_controls", "negative_seeds")
    expected_checks = {
        "all_nine_learned_runs_pass": all(run["passed"] for run in learned),
        "all_nine_full_target_oracles_pass": all(run["passed"] for run in oracles),
        "all_nine_support_holes_pass_usual_screens": all(run["usual_passed"] for run in negatives),
        "all_nine_support_holes_fail_truth": all(not run["truth_passed"] for run in negatives),
    }
    require(checks(source["checks"], "challenge checks") == expected_checks, "challenge summary disagrees with runs")
    passed = boolean(source["passed"], "challenge passed")
    require(passed == all(expected_checks.values()), "challenge passed disagrees with summary")
    groups = []
    for fixture, label in CHALLENGE_LABELS.items():
        runs = [{key: value for key, value in run.items()
                 if key not in {"fixture", "usual_passed", "truth_passed"}}
                for run in learned if run["fixture"] == fixture]
        groups.append({"id": fixture, "label": label, "passed_runs": sum(run["passed"] for run in runs),
                       "total_runs": len(runs), "runs": sorted(runs, key=lambda run: run["seed"])})
    return {"passed": passed, "passed_runs": sum(run["passed"] for run in learned),
            "total_runs": len(learned), "groups": groups,
            "controls": {"oracle_passed": sum(run["passed"] for run in oracles), "oracle_total": len(oracles),
                         "negative_failed_truth": sum(not run["truth_passed"] for run in negatives),
                         "negative_passed_usual": sum(run["usual_passed"] for run in negatives),
                         "negative_total": len(negatives)},
            "thresholds": {key: number(value, f"challenge gate {key}")
                           for key, value in mapping(source["specification"]["gates"], "challenge gates").items()}}


def build_snapshot(root=ROOT):
    """Read and validate every source before returning one complete snapshot."""
    raw = {path: (Path(root) / path).read_bytes() for path in SOURCE_PATHS}
    design, profile, atlas, sampling, challenge = [parse_json(raw[path]) for path in SOURCE_PATHS[:-1]]
    records = [parse_json(line) for line in raw["corpus/INDEX.jsonl"].splitlines() if line.strip()]
    require(bool(records), "empty archive index")
    pmids = [record["pmid"] for record in records]
    require(all(isinstance(pmid, str) and pmid.isdigit() for pmid in pmids), "invalid archive PMID")
    require(len(set(pmids)) == len(pmids), "duplicate archive PMID")
    census = census_summary(design, profile, atlas, len(records))
    mechanisms, sites, growth_period = mechanism_summary(profile)
    return {"schema_version": 1, "snapshot_date": SNAPSHOT_DATE,
            "sources": [{"path": path, "sha256": hashlib.sha256(raw[path]).hexdigest(),
                         "url": REPOSITORY + path} for path in SOURCE_PATHS],
            "census": census, "mechanisms": mechanisms, "sites": sites, "growth_period": growth_period,
            "sampling": sampling_summary(sampling), "challenge": challenge_summary(challenge)}


def encode_snapshot(snapshot):
    return (json.dumps(snapshot, ensure_ascii=False, sort_keys=True, separators=(",", ":"),
                       allow_nan=False) + "\n").encode("utf-8")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="check freshness without writing")
    parser.add_argument("--root", type=Path, default=ROOT, help="repository root")
    args = parser.parse_args(argv)
    output = args.root / OUTPUT
    temporary = None
    try:
        encoded = encode_snapshot(build_snapshot(args.root))
        if args.check:
            require(output.is_file() and output.read_bytes() == encoded,
                    f"{OUTPUT} is missing or stale; run python scripts/build_pages_data.py")
            print(f"Fresh: {OUTPUT}")
        else:
            output.parent.mkdir(parents=True, exist_ok=True)
            with tempfile.NamedTemporaryFile(dir=output.parent, prefix=".research-data-", delete=False) as stream:
                temporary = Path(stream.name)
                stream.write(encoded)
            os.replace(temporary, output)
            print(f"Wrote {OUTPUT} ({len(encoded):,} bytes)")
        return 0
    except (OSError, ValueError, KeyError, TypeError, AttributeError) as error:
        print(f"Pages data error: {error}", file=sys.stderr)
        return 1
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


if __name__ == "__main__":
    raise SystemExit(main())
