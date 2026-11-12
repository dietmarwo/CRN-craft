"""Added time and memory of quantile statistics (task M4e; no target).

Declared in `validation/merge/RESULTS.md`: 4,000 cells of
`models/rules/mtdna-three-class-cme.rules.json` (6 observables) under
`static_direct` on the 37-point grid of `validation/check_compiled_loop.py`
(every 30 days to 1,080 days, in hours), 2 threads, seed 1, with and
without the quantiles 0.025, 0.5 and 0.975. Each configuration runs in fresh
processes (three, alternating with the other configuration); each process
prepares once, simulates once untimed, then five timed times, and reports
its peak resident memory. The quantile buffer is 4,000 x 37 x 6 x 8 bytes.

    flock /home/claude/craft5/cpu.lock python validation/merge/measure_quantiles.py \\
        --output quantile-cost.json

Run it with the interpreter of the wheel under test.
"""

import argparse
import json
from pathlib import Path
import platform
import statistics
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
CHILD = r"""
import json, resource, sys, time
from crn_craft import rules
quantiles = json.loads(sys.argv[2])
model = rules.load(sys.argv[1])
run = rules.ensemble(
    model,
    end=25920.0,
    report_times=[24.0 * day for day in range(0, 1081, 30)],
    replicates=4000,
    seed=1,
    threads=2,
    method="static_direct",
    quantiles=quantiles or None,
)
run.simulate()
times = []
for _ in range(5):
    start = time.perf_counter()
    result = run.simulate()
    times.append(time.perf_counter() - start)
stats = result.statistics
print(json.dumps({
    "seconds": times,
    "peak_rss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
    "events": result.metadata["events"],
    "cells": len(stats["mean"]) * len(stats["mean"].columns),
    "quantiles": list(stats.get("quantiles", {})),
}))
"""


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    model = ROOT / "models/rules/mtdna-three-class-cme.rules.json"
    configurations = {"none": [], "quantiles": [0.025, 0.5, 0.975]}
    samples = {name: [] for name in configurations}
    for _ in range(3):
        for name, quantiles in configurations.items():
            done = subprocess.run(
                [sys.executable, "-c", CHILD, str(model), json.dumps(quantiles)],
                capture_output=True,
                text=True,
                check=True,
            )
            samples[name].append(json.loads(done.stdout))
    summary = {}
    for name, runs in samples.items():
        times = [t for run in runs for t in run["seconds"]]
        summary[name] = {
            "median_seconds": statistics.median(times),
            "min_seconds": min(times),
            "max_seconds": max(times),
            "peak_rss_kib": [run["peak_rss_kib"] for run in runs],
            "events": runs[0]["events"],
            "quantiles": runs[0]["quantiles"],
        }
    if summary["none"]["events"] != summary["quantiles"]["events"]:
        raise SystemExit("measure_quantiles: the two configurations ran differently")
    buffer = 4000 * 37 * 6 * 8
    report = {
        "schema": "crn-craft.merge-quantile-cost/0.1",
        "host": {"machine": platform.machine(), "python": platform.python_version()},
        "replicates": 4000,
        "report_times": 37,
        "observables": 6,
        "buffer_bytes": buffer,
        "sort_column_bytes": 4000 * 8,
        "configurations": summary,
        "added_median_seconds": summary["quantiles"]["median_seconds"]
        - summary["none"]["median_seconds"],
        "added_peak_rss_kib": statistics.median(summary["quantiles"]["peak_rss_kib"])
        - statistics.median(summary["none"]["peak_rss_kib"]),
        "samples": samples,
    }
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(
        json.dumps(
            {key: report[key] for key in ("added_median_seconds", "added_peak_rss_kib")}
            | {name: summary[name]["median_seconds"] for name in summary}
        )
    )


if __name__ == "__main__":
    main()
