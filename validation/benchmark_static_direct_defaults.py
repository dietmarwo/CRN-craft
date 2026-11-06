"""Static Direct: the rules-run/0.6 default against the rules-run/0.5 default.

The workload is the 1000-cell three-class case of the methods section of
examples/mtdna/performance.py: the takeover stop, three years, monthly
reports, seed 17, two threads and two count trackers. The rules-run/0.5
request selects Static Direct without options (inverse waits, per-event
trackers: what Python's ensemble(method="static_direct") ran in 0.3.0); the
rules-run/0.6 request is the one that call writes now (Ziggurat waits,
sparse trackers). The two run in alternating order, 0.5 first, and each is
summarized by its median.

Run it under the CPU lock with the wheel built from the tree under test:

    python validation/benchmark_static_direct_defaults.py --output out.json
"""

import argparse
import hashlib
import json
import os
import platform
from pathlib import Path
import subprocess
import sys
import time

import numpy as np
from crn_craft import rules

ROOT = Path(__file__).resolve().parents[1]
TARGET = 1.10  # declared before the first measurement
CELLS = 1000
GRID = 24.0 * np.arange(0, 1081, 30)
TRACK = [{"observable": "neutral_mutants"}, {"observable": "advantaged_mutants"}]
SAMPLERS = {
    "0.5 default": "chacha8-u53/rules-static-direct/1",
    "0.6 default": "chacha8-u53-ziggurat/rules-static-direct-sparse-trackers/1",
}


def model_data():
    model = json.loads(rules.load(ROOT / "models/rules/mtdna-three-class.rules.json"))
    model["stop"] = "N > 0 && mutants >= 0.9 * N"
    return json.dumps(model).encode()


def released_request(data):
    """A rules-run/0.5 Static Direct request without options."""
    return {
        "schema": "crn-craft.rules-run/0.5",
        "model": {
            "path": "model.rules.json",
            "sha256": hashlib.sha256(data).hexdigest(),
        },
        "time": {"end": 25_920.0, "report_times": [float(t) for t in GRID]},
        "ensemble": {"replicates": CELLS, "root_seed": 17},
        "outputs": {
            "statistics": ["mean", "variance"],
            "keep_paths": 0,
            "keep_instances": 0,
        },
        "budget": {
            "threads": 2,
            "max_events_per_replicate": 1_000_000,
            "max_sampling_steps_per_replicate": 8_000_000,
            "max_instances": 100_000,
            "max_output_bytes": 64 << 20,
        },
        "trackers": TRACK,
        "method": {"kind": "static_direct"},
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--repeats", type=int, default=3)
    args = parser.parse_args()
    data = model_data()
    runs = {
        "0.5 default": rules.prepare(data, released_request(data)),
        "0.6 default": rules.ensemble(
            data,
            end=25_920.0,
            report_times=GRID,
            replicates=CELLS,
            seed=17,
            threads=2,
            method="static_direct",
            trackers=TRACK,
        ),
    }
    rows = []
    for repeat in range(args.repeats):
        for label, run in runs.items():
            start = time.perf_counter()
            result = run.simulate()
            seconds = time.perf_counter() - start
            meta = result.metadata
            if meta["sampler"] != SAMPLERS[label] or meta["termination"] != "completed":
                raise RuntimeError(f"{label}: {meta['sampler']} {meta['termination']}")
            rows.append(
                {
                    "repeat": repeat,
                    "case": label,
                    "schema": meta["schema"],
                    "sampler": meta["sampler"],
                    "seconds": seconds,
                    "events": meta["events"],
                    "ns_per_event": 1e9 * seconds / meta["events"],
                    "stopped": meta["stopped_replicates"],
                }
            )
            print(
                f"{repeat} {label}: {seconds:6.2f} s {rows[-1]['ns_per_event']:5.1f} "
                f"ns/event {meta['events']:,} events",
                flush=True,
            )
    summary = {}
    for label in runs:
        mine = [r for r in rows if r["case"] == label]
        summary[label] = {
            "median_seconds": float(np.median([r["seconds"] for r in mine])),
            "median_ns_per_event": float(np.median([r["ns_per_event"] for r in mine])),
            "events": mine[0]["events"],
        }
    old, new = summary["0.5 default"], summary["0.6 default"]
    speedup = old["median_ns_per_event"] / new["median_ns_per_event"]
    report = {
        "workload": (
            "examples/mtdna/performance.py methods section: three-class model, "
            "takeover stop, 1000 cells, 25920 h, monthly reports, seed 17, "
            "two threads, two count trackers"
        ),
        "order": "alternating, 0.5 default first",
        "metric": "median wall ns per event; wall seconds reported as well",
        "target_speedup": TARGET,
        "host": {
            "platform": platform.platform(),
            "cpus": len(os.sched_getaffinity(0)),
            "python": sys.version.split()[0],
        },
        "commit": subprocess.run(
            ["git", "-C", str(ROOT), "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            check=False,
        ).stdout.strip(),
        "rows": rows,
        "summary": summary,
        "speedup_ns_per_event": speedup,
        "speedup_seconds": old["median_seconds"] / new["median_seconds"],
        "met": speedup >= TARGET,
    }
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"speedup": round(speedup, 3), "met": report["met"]}))


if __name__ == "__main__":
    main()
