#!/usr/bin/env python3
"""Summarize complete smoke and study runs into results/verified.json.

Each directory holds C1/, C2/, C3/ and C4/ with the summary.json that run.py
writes (C1 with --include-next-reaction). Only complete, passed runs are recorded;
the record keeps the numbers that VALIDATION.md cites, not host paths or
timings.

`--publications DIR` records the smoke runs of C5 to C10 instead, in their own
file (results/publications-smoke.json): DIR holds C5/ to C10/, and the record
keeps each chapter's check names and verdicts and the package versions, no
numbers, paths or timings. The study numbers of those chapters are the tables
under results/C5/ to results/C10/.
"""

import argparse
import json
from pathlib import Path

CHAPTERS = ["C1", "C2", "C3", "C4"]
PUBLICATIONS = ["C5", "C6", "C7", "C8", "C9", "C10"]


def fail(field, message):
    raise ValueError(f"record.{field}: {message}")


def pick(row, keys):
    return {k: row[k] for k in keys}


def c1(s):
    policies = sorted({a["policy"] for a in s["arms"]})
    if policies != ["direct", "next_reaction"]:
        fail("C1.policies", "record C1 with --include-next-reaction")
    return dict(
        arms=[
            pick(
                a,
                [
                    "model",
                    "policy",
                    "sampler",
                    "request_schema",
                    "route",
                    "seed",
                    "rates",
                    "truth",
                    "training_loss",
                    "truth_training_loss",
                    "best_is_start",
                    "holdout_loss",
                    "truth_holdout_loss",
                    "parameter_relative_error",
                    "independent_moment_max_standard_errors",
                ],
            )
            for a in s["arms"]
        ]
    )


def c2(s):
    audit = s["audit"]
    return dict(
        audit=dict(
            pick(
                audit,
                [
                    "gradient_error",
                    "integration_error",
                    "tightened_error",
                    "full_observation_rank",
                    "terminal_product_twin",
                ],
            ),
            condition_numbers={
                k: v["condition_number"]
                for k, v in audit["local_identifiability"].items()
            },
        ),
        arms=[
            pick(
                a,
                [
                    "strategy",
                    "rates",
                    "native_loss",
                    "independent_loss",
                    "parameter_relative_error",
                    "native_solves",
                    "native_failures",
                    "stochastic_vs_ode_normalized_rmse",
                    "validation_sampler",
                    "validation_request_schema",
                ],
            )
            for a in s["arms"]
        ],
    )


def c3(s):
    return dict(
        pick(
            s,
            [
                "pareto_points",
                "occupied_niches",
                "failed_evaluations",
                "clipped_descriptors",
                "generator_max_error",
            ],
        ),
        validation=[
            dict(
                pick(
                    v,
                    [
                        "origin",
                        "training_niche",
                        "heldout_niche",
                        "niche_changed",
                        "descriptor_distance",
                        "training_oscillatory_fraction",
                        "heldout_oscillatory_fraction",
                    ],
                ),
                training_descriptors=v["training"]["descriptors"],
            )
            for v in s["validation"]
        ],
        backend_control=[
            pick(r, ["max_standard_errors", "crn_sampler", "reference_backend"])
            for r in s["backend_control"]["rows"]
        ],
    )


def c4(s):
    return dict(
        arms=[
            pick(
                a,
                [
                    "kernel",
                    "fitted_rate",
                    "true_rate",
                    "training_loss",
                    "best_is_start",
                    "exact_law_optimal_rate",
                    "exact_law_optimal_loss",
                    "exact_heldout_prediction_mse",
                    "heldout_law_max_cdf_distance",
                    "dkw_bound",
                ],
            )
            for a in s["arms"]
        ]
    )


def load(folder, name):
    path = folder / name / "summary.json"
    if not path.is_file():
        fail(f"{name}.summary", f"{path.name} missing for {name}")
    s = json.loads(path.read_text())
    if s.get("state") != "complete" or s.get("passed") is not True:
        fail(f"{name}.state", f"{name} did not complete with passed checks")
    return s


def chapter(folder, name):
    s = load(folder, name)
    return s, dict(
        passed=True,
        checks=[c["name"] for c in s["checks"]],
        budgets=s["budgets"],
        results={"C1": c1, "C2": c2, "C3": c3, "C4": c4}[name](s),
    )


def publications(folder):
    record = dict(
        schema_version="crn-craft.optimization-record/0.1",
        scope="Smoke runs of the publication tutorials: which checks passed.",
        chapters={},
    )
    versions = {}
    for name in PUBLICATIONS:
        summary = load(folder, name)
        if summary["preset"] != "smoke":
            fail(f"{name}.preset", f"{folder.name}/{name} is not a smoke run")
        if not summary["checks"]:
            fail(f"{name}.checks", f"{name} ran no checks")
        for package, version in summary["versions"].items():
            if versions.setdefault(package, version) != version:
                fail(f"{name}.versions", f"runs used different {package} versions")
        record["chapters"][name] = [
            pick(c, ["name", "passed"]) for c in summary["checks"]
        ]
    record["versions"] = versions
    return record


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--smoke", type=Path)
    p.add_argument("--study", type=Path)
    p.add_argument("--publications", type=Path)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    if args.publications is not None:
        if args.smoke or args.study:
            fail("arguments", "--publications excludes --smoke and --study")
        record = publications(args.publications)
        args.output.write_text(json.dumps(record, indent=2, allow_nan=False) + "\n")
        return
    if args.smoke is None or args.study is None:
        fail("arguments", "give --smoke and --study, or --publications")
    record = dict(
        schema_version="crn-craft.optimization-record/0.1",
        scope="Executed teaching workflows, not optimizer-performance claims.",
        presets={},
    )
    versions = None
    for preset, folder in [("smoke", args.smoke), ("study", args.study)]:
        record["presets"][preset] = {}
        for name in CHAPTERS:
            summary, row = chapter(folder, name)
            if summary["preset"] != preset:
                fail(f"{name}.preset", f"{folder.name}/{name} is not a {preset} run")
            if versions is None:
                versions = summary["versions"]
            elif summary["versions"] != versions:
                fail(f"{name}.versions", "runs used different package versions")
            record["presets"][preset][name] = row
    record["versions"] = versions
    args.output.write_text(json.dumps(record, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
