"""Time compiled Static Direct event loops end to end from Python.

Declared before measuring (docs/decisions/compiled-rules-loop-python-2026-11-06.md):
on the 1000-cell three-class workload of docs/evidence/mtdna-performance-analysis.md
(takeover stop, three years, monthly reports, seed 17, two threads, two
count trackers), `run.simulate(compiled=loop)` must be at least 1.5x faster
than interpreted `run.simulate()` under `method="static_direct"`, by the
median of three runs in alternating order. The row `static_direct` uses
inverse waits and per-event trackers, the meaning of `"static_direct"` in
rules-run/0.2-0.5; the row `ziggurat + sparse trackers` uses the options
that rules-run/0.6 makes the default. Times cover the whole
`simulate()` call, including the subprocess and parsing its output, and
exclude the build. Every timed pair must be byte-identical.

The build is timed twice: into an empty cache and target directory, then
again from the cache.

Usage: python validation/benchmark_compiled_loop.py --output timing.json
[--source CHECKOUT]. Run it under the host's CPU lock. Timing is
host-specific.
"""

import argparse
import json
import platform
from pathlib import Path
import sys
import tempfile
import time

from crn_craft import rules

ROOT = Path(__file__).resolve().parents[1]
TARGET = 1.5
GRID = [24.0 * day for day in range(0, 1081, 30)]
TRACK = [{"observable": "neutral_mutants"}, {"observable": "advantaged_mutants"}]
# Explicit options, so the rows keep their meaning whatever the request
# version's defaults are.
METHODS = {
    "static_direct": {
        "kind": "static_direct",
        "waiting_time": "inverse",
        "sparse_trackers": False,
    },
    "ziggurat + sparse trackers": {
        "kind": "static_direct",
        "waiting_time": "ziggurat",
        "sparse_trackers": True,
    },
}


def require(condition, message):
    if not condition:
        raise SystemExit(f"benchmark_compiled_loop: {message}")


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--source", type=Path, help="checkout; default crates.io")
    parser.add_argument("--cells", type=int, default=1000)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    model = json.loads(rules.load(ROOT / "models/rules/mtdna-three-class.rules.json"))
    model["stop"] = "N > 0 && mutants >= 0.9 * N"
    report = {
        "platform": platform.platform(),
        "python": sys.version.split()[0],
        "crn_craft": rules._native.__version__,
        "source": "crates.io" if args.source is None else "checkout",
        "declared_target": f"compiled >= {TARGET}x interpreted static_direct",
        "cells": args.cells,
        "rows": [],
    }
    with tempfile.TemporaryDirectory() as cache:
        build = {"cache_dir": cache, "source": args.source, "progress": False}
        start = time.perf_counter()
        loop = rules.compile_loop(model, **build)
        report["build_first_s"] = time.perf_counter() - start
        report["cargo_first_s"] = loop.build_seconds
        start = time.perf_counter()
        again = rules.compile_loop(model, **build)
        report["build_cached_s"] = time.perf_counter() - start
        require(again.cached, "the second build was not cached")
        print(
            f"build: first {report['build_first_s']:.1f} s "
            f"(cargo {loop.build_seconds:.1f} s), cached {report['build_cached_s']:.3f} s",
            flush=True,
        )
        for label, method in METHODS.items():
            run = rules.ensemble(
                model,
                end=25_920.0,
                report_times=GRID,
                replicates=args.cells,
                seed=17,
                threads=2,
                method=method,
                trackers=TRACK,
            )
            payloads = {}
            for repeat in range(args.repeats):
                order = [("interpreted", None), ("compiled", loop)]
                for route, compiled in order[:: -1 if repeat % 2 else 1]:
                    start = time.perf_counter()
                    result = run.simulate(compiled=compiled)
                    seconds = time.perf_counter() - start
                    meta = result.metadata
                    payloads.setdefault(
                        route, (meta, result.records, result.rule_events)
                    )
                    row = {
                        "method": label,
                        "route": route,
                        "repeat": repeat,
                        "seconds": seconds,
                        "events": meta["events"],
                        "ns_per_event": 1e9 * seconds / meta["events"],
                        "sampler": meta["sampler"],
                    }
                    report["rows"].append(row)
                    print(json.dumps(row), flush=True)
            (meta, records, counts), other = (
                payloads["interpreted"],
                payloads["compiled"],
            )
            require(
                meta == other[0]
                and records.equals(other[1])
                and counts.equals(other[2]),
                f"{label}: compiled and interpreted results differ",
            )
            # Untimed: metadata and record columns byte for byte.
            loop.verify(run)
    summary = {}
    for label in METHODS:
        medians = {}
        for route in ("interpreted", "compiled"):
            rows = [
                r
                for r in report["rows"]
                if r["method"] == label and r["route"] == route
            ]
            seconds = sorted(r["seconds"] for r in rows)[len(rows) // 2]
            medians[route] = {
                "median_seconds": seconds,
                "ns_per_event": 1e9 * seconds / rows[0]["events"],
            }
        medians["speedup"] = (
            medians["interpreted"]["median_seconds"]
            / medians["compiled"]["median_seconds"]
        )
        summary[label] = medians
        print(f"{label}: {medians['speedup']:.2f}x", flush=True)
    report["summary"] = summary
    report["target_met"] = summary["static_direct"]["speedup"] >= TARGET
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"target_met": report["target_met"]}))


if __name__ == "__main__":
    main()
