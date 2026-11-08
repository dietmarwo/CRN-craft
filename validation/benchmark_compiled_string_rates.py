"""Time a compiled event loop with string rates end to end from Python.

Declared before measuring (2026-11-08): on the SIR model as `rules.Model`
writes it, with the string rates `"beta * S * I / N"` and `"gamma * I"`
(whole-rule propensities), scaled to N = 10,000 (S = 9,990 and I = 10 at
time 0, beta = 0.3, gamma = 0.1), 2,000 replicates to time 400 with
reports every 10 time units, seed 19, two threads and no trackers,
`run.simulate(compiled=loop)` must be at least 1.5x faster than interpreted
`run.simulate()` under Static Direct, by the median of three runs in
alternating order, in both rows: inverse waits, and Ziggurat waits (the
rules-run/0.6 default). Times cover the whole `simulate()` call, including
the subprocess and parsing its output, and exclude the build. Every timed
cell is checked byte for byte with `CompiledLoop.verify` afterwards
(untimed).

Usage: python validation/benchmark_compiled_string_rates.py --output
timing.json [--source CHECKOUT] [--cache-dir DIR]. Run it under the host's
CPU lock. Timing is host-specific.
"""

import argparse
import json
from pathlib import Path
import platform
import sys
import time

from crn_craft import rules

TARGET = 1.5
METHODS = {
    "inverse": {"kind": "static_direct", "waiting_time": "inverse"},
    "ziggurat": {"kind": "static_direct", "waiting_time": "ziggurat"},
}


def sir(population):
    """The document `rules.Model` writes for the reported SIR model."""
    return {
        "schema": "crn-craft.rules-model/0.3",
        "name": "sir",
        "parameters": {"beta": 0.3, "gamma": 0.1, "N": float(population)},
        "types": [{"name": "S"}, {"name": "I"}, {"name": "R"}],
        "initial": [
            {"type": "S", "count": population - 10},
            {"type": "I", "count": 10},
        ],
        "observables": {s: {"type": s} for s in ("S", "I", "R")},
        "rules": [
            {
                "name": "infection",
                "reactants": [
                    {"label": "x", "type": "S"},
                    {"label": "y", "type": "I"},
                ],
                "propensity": "beta * S * I / N",
                "products": [{"type": "I"}, {"type": "I"}],
            },
            {
                "name": "recovery",
                "reactants": [{"label": "x", "type": "I"}],
                "propensity": "gamma * I",
                "products": [{"type": "R"}],
            },
        ],
    }


def require(condition, message):
    if not condition:
        raise SystemExit(f"benchmark_compiled_string_rates: {message}")


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--source", type=Path, help="checkout; default crates.io")
    parser.add_argument("--cache-dir", type=Path)
    parser.add_argument("--replicates", type=int, default=2000)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    model = sir(10_000)
    report = {
        "platform": platform.platform(),
        "python": sys.version.split()[0],
        "crn_craft": rules._native.__version__,
        "source": "crates.io" if args.source is None else "checkout",
        "declared_target": f"compiled >= {TARGET}x interpreted in every row",
        "replicates": args.replicates,
        "rows": [],
    }
    loop = rules.compile_loop(
        model, source=args.source, cache_dir=args.cache_dir, progress=False
    )
    report["cached_build"] = loop.cached
    for label, method in METHODS.items():
        run = rules.ensemble(
            model,
            end=400.0,
            report_times=[10.0 * i for i in range(41)],
            replicates=args.replicates,
            seed=19,
            threads=2,
            method=method,
        )
        for repeat in range(args.repeats):
            order = [("interpreted", None), ("compiled", loop)]
            for route, compiled in order[:: -1 if repeat % 2 else 1]:
                start = time.perf_counter()
                result = run.simulate(compiled=compiled)
                seconds = time.perf_counter() - start
                meta = result.metadata
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
        # Untimed: metadata and record columns byte for byte.
        require(loop.verify(run).metadata["events"] > 0, f"{label}: verify")
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
    report["target_met"] = all(s["speedup"] >= TARGET for s in summary.values())
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"target_met": report["target_met"]}))


if __name__ == "__main__":
    main()
