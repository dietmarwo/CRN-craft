"""Reproduce the fixed-mtDNA specialization benchmark on one host."""

import argparse
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import tempfile
from time import perf_counter

import numpy as np

from benchmark_mtdna import GRID, fixed_cross_check, run_java
from check_mtdna import cdf_distance, read_model, require, rules_files

TRACKERS = [{"observable": name} for name in ("neutral_mutants", "advantaged_mutants")]
STATIC = "chacha8-u53/rules-static-direct/1"
FINAL = "chacha8-u53-ziggurat/rules-mtdna-direct/1"


def run_crn(binary, folder, rows, model, label, cells, workers, kind, phases=False):
    paths = rules_files(
        folder,
        model,
        n=cells,
        threads=workers,
        times=GRID,
        method=kind,
        seed=17,
        trackers=TRACKERS,
        keep=True,
    )
    args = [str(binary), *([] if phases else ["rules"]), *map(str, paths)]
    start = perf_counter()
    run = subprocess.run(args, capture_output=True, check=False, timeout=900)
    process = perf_counter() - start
    require(run.returncode == 0, run.stderr.decode())
    result = json.loads(run.stdout)
    elapsed = perf_counter() - start
    require(result["termination"] == "completed", str(result["diagnostics"]))
    sampler = result["sampler"]
    expected = {"static": STATIC, "final": FINAL, "phases": FINAL}
    require(sampler == expected.get(label, "chacha8-u53/rules-mtdna-direct/1"), sampler)
    row = {
        "case": label,
        "method": sampler,
        "cells": cells,
        "threads": workers,
        "seconds": elapsed,
        "process_seconds": process,
        "decode_seconds": elapsed - process,
        "events": result["events"],
        "stopped": result["stopped_replicates"],
        "output_bytes": len(run.stdout),
    }
    if phases:
        row["phases"] = json.loads(run.stderr)
    rows.append(row)
    print(json.dumps(row), flush=True)
    return result


def compare_crn(base, final):
    # Compare each of six outcomes with a two-sample DKW union bound.
    # Report bins include the never-takeover mass.
    keys = ["N", "wild_type", "neutral_mutants", "advantaged_mutants", "mutants"]
    bound = 2 * np.sqrt(np.log(4 * 6 / 0.001) / (2 * len(base["paths"])))
    distances = {}
    for key in keys:
        columns = [
            [p["values"][-1][r["observables"].index(key)] for p in r["paths"]]
            for r in (base, final)
        ]
        distances[key] = cdf_distance(*(np.asarray(v) for v in columns))
    stops = [
        np.asarray(
            [
                np.inf if p["termination"] != "stopped" else p["stop_time"]
                for p in r["replicates"]
            ]
        )
        for r in (base, final)
    ]
    distances["first_passage"] = cdf_distance(*stops)
    return {
        "familywise_alpha": 0.001,
        "bound": float(bound),
        "distances": distances,
        "passed": bool(max(distances.values()) <= bound),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("static", "fused", "lazy", "final", "phases", "java", "output"):
        parser.add_argument(
            f"--{name}", type=Path, required=name not in {"fused", "lazy", "phases"}
        )
    parser.add_argument("--java-classpath", required=True)
    parser.add_argument("--java-probe-classpath")
    parser.add_argument("--full-repeats", type=int, default=2)
    args = parser.parse_args()
    rows = []
    cpu_max = Path("/sys/fs/cgroup/cpu.max")
    evidence = {
        "platform": platform.platform(),
        "cpu_max": cpu_max.read_text().strip() if cpu_max.exists() else None,
        "binary_sha256": {
            name: hashlib.sha256(path.read_bytes()).hexdigest()
            for name in ("static", "fused", "lazy", "final", "phases")
            if (path := getattr(args, name)) is not None
        },
        "rows": rows,
    }

    def save():
        args.output.write_text(json.dumps(evidence, indent=2) + "\n")

    model = read_model("mtdna-three-class.rules.json")
    model["stop"] = "N > 0 && mutants >= 0.9 * N"
    with tempfile.TemporaryDirectory() as directory:
        folder = Path(directory)
        for repeat in range(args.full_repeats):
            # Reverse order on even-numbered observations.
            order = ["java", "final", "static"][:: (-1 if repeat % 2 else 1)]
            work = {}
            for label in order:
                if label == "java":
                    work[label] = run_java(
                        args, folder, rows, "fixedPop", "ssa", 10_000, 8
                    )
                else:
                    kind = "static_direct" if label == "static" else "mtdna_direct"
                    work[label] = run_crn(
                        getattr(args, label),
                        folder,
                        rows,
                        model,
                        label,
                        10_000,
                        8,
                        kind,
                    )
                save()
            if repeat == 0:
                evidence["crn_law_cross_check"] = compare_crn(
                    work["static"], work["final"]
                )
                evidence["java_cross_check"] = fixed_cross_check(
                    work["final"], *work["java"]
                )
                save()
                require(
                    evidence["crn_law_cross_check"]["passed"], "CRN two-sample gate"
                )
                require(
                    evidence["java_cross_check"]["passed"], "Java first-passage gate"
                )

        for workers in (1, 8):
            for label in ("static", "fused", "lazy", "final"):
                if (binary := getattr(args, label)) is not None:
                    kind = "static_direct" if label == "static" else "mtdna_direct"
                    run_crn(binary, folder, rows, model, label, 512, workers, kind)
                    save()
            run_java(args, folder, rows, "fixedPop", "ssa", 512, workers)
            save()

        if args.phases is not None:
            run_crn(
                args.phases,
                folder,
                rows,
                model,
                "phases",
                10_000,
                8,
                "mtdna_direct",
                True,
            )
            save()
        if args.java_probe_classpath:
            args.java_classpath = args.java_probe_classpath
            run_java(args, folder, rows, "fixedPop", "ssa", 10_000, 8)
            save()


if __name__ == "__main__":
    main()
