"""Time hybrid with independent against dependent continuous variables.

Declared before measuring in validation/hybrid/RESULTS.md (section E): model
C4 of validation/check_hybrid.py (an eight-compartment linear chain, all
variables independent, drives arrivals; deaths at rate 0.1), 4,000
replicates to T = 12, seed 6301, two threads, rtol 1e-6 and atol 1e-9.
`rules.ensemble(...)` plus `simulate()` is timed for `independent: true`
and for `false`, three runs each in alternating order; the median with
`true` must be at least 2x faster. The arms' results differ within the
tolerances, so they are reported, not compared.

Usage: python validation/benchmark_independent_ode.py --output timing.json.
Run it under the host's CPU lock. Timing is host-specific.
"""

import argparse
import json
from pathlib import Path
import platform
import sys
import time

from check_hybrid import HYBRID, column, pulse
from crn_craft import rules
import numpy as np

TARGET = 2.0


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--replicates", type=int, default=4000)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = {
        "platform": platform.platform(),
        "python": sys.version.split()[0],
        "crn_craft": rules._native.__version__,
        "declared_target": f"independent >= {TARGET}x faster than dependent",
        "model": "C4: eight-compartment pulse, rules-model/0.4",
        "replicates": args.replicates,
        "rows": [],
    }
    arms = (("independent", True), ("dependent", False))
    for repeat in range(args.repeats):
        for label, flag in arms[:: -1 if repeat % 2 else 1]:
            start = time.perf_counter()
            run = rules.ensemble(
                pulse(flag),
                end=12.0,
                report_times=[12.0],
                replicates=args.replicates,
                seed=6301,
                threads=2,
                method=HYBRID,
                keep_paths=args.replicates,
            )
            result = run.simulate()
            seconds = time.perf_counter() - start
            meta = result.metadata
            row = {
                "arm": label,
                "repeat": repeat,
                "seconds": seconds,
                "events": meta["events"],
                "completed": meta["completed_replicates"],
                "mean_count": float(np.mean(column(meta, "A")[:, -1])),
            }
            report["rows"].append(row)
            print(json.dumps(row), flush=True)
    medians = {}
    for label, _ in arms:
        seconds = sorted(r["seconds"] for r in report["rows"] if r["arm"] == label)
        medians[label] = seconds[len(seconds) // 2]
    report["median_seconds"] = medians
    report["speedup"] = medians["dependent"] / medians["independent"]
    report["target_met"] = report["speedup"] >= TARGET
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(
        json.dumps({"speedup": report["speedup"], "target_met": report["target_met"]})
    )


if __name__ == "__main__":
    main()
