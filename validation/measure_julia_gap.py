"""Paired timing of exact SSA engines and a projected Julia gap.

Build `cargo build --release --example bench_ssa` at two commits and copy each
binary aside, then:

    python validation/measure_julia_gap.py pair --before OLD --after NEW \
        --output pair.json
    python validation/measure_julia_gap.py project pair.json

`pair` alternates parent and candidate binaries per round on W1-W7 of the
acceleration report and records every timed repeat, the event count and the
result digest. `project` multiplies the report's CRN prepared medians (another
host) by the measured paired speedups and divides by Julia's fastest median.
That projection is not a Julia measurement; see
docs/decisions/julia-gap-2026-11-02.md.
"""

import argparse
import json
import math
import statistics
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
EVIDENCE = ROOT / "validation" / "ssa-acceleration"
WORKLOADS = {
    "W1": (100.0, 64),
    "W2": (100.0, 64),
    "W3": (100.0, 64),
    "W4": (100.0, 64),
    "W5": (20.0, 16),
    "W6": (20.0, 16),
    "W7": (20.0, 16),
}
METHODS = [
    ("legacy_cached", "inverse"),
    ("cached_direct", "ziggurat"),
    ("direct_cr", "ziggurat"),
    ("sorting_direct", "ziggurat"),
    ("rssa_cr", "ziggurat"),
]
REPORT_ENGINE = {
    "legacy_cached": "legacy_cached/prepared",
    "cached_direct": "cached_direct/ziggurat/prepared",
    "direct_cr": "direct_cr/ziggurat/prepared",
    "sorting_direct": "sorting_direct/ziggurat/prepared",
    "rssa_cr": "rssa_cr/ziggurat/prepared",
}


def bench(binary, workload, method, wait, repeats):
    horizon, paths = WORKLOADS[workload]
    model = EVIDENCE / "models" / f"{workload}.model.json"
    command = [binary, model, horizon, paths, method, wait, repeats]
    output = subprocess.run(
        [str(part) for part in command], capture_output=True, text=True, check=True
    ).stdout
    return json.loads(output)


def pair(args):
    methods = [m for m in METHODS if not args.methods or m[0] in args.methods]
    report = {"protocol": "docs/decisions/julia-gap-2026-11-02.md", "rows": {}}
    for workload in args.workloads:
        for method, wait in methods:
            row = {"before": [], "after": [], "digests": {}, "events": {}}
            for round_index in range(args.rounds):
                order = [("before", args.before), ("after", args.after)]
                if round_index % 2:
                    order.reverse()
                for label, binary in order:
                    result = bench(binary, workload, method, wait, args.repeats)
                    row[label].append(result["seconds"])
                    row["digests"].setdefault(label, set()).add(result["digest"])
                    row["events"][label] = result["events"]
            medians = {
                k: [statistics.median(r) for r in row[k]] for k in ("before", "after")
            }
            ratios = [
                b / a for b, a in zip(medians["before"], medians["after"], strict=True)
            ]
            row["digests"] = {k: sorted(v) for k, v in row["digests"].items()}
            row["identical"] = len(set().union(*map(set, row["digests"].values()))) == 1
            row["ratios"] = ratios
            row["speedup"] = statistics.median(ratios)
            row["rounds_faster"] = sum(r > 1.0 for r in ratios)
            key = f"{workload}/{method}/{wait}"
            report["rows"][key] = row
            print(
                f"{key:30s} before {statistics.median(medians['before']):.4f} "
                f"after {statistics.median(medians['after']):.4f} "
                f"speedup {row['speedup']:.3f} ({row['rounds_faster']}/{args.rounds}) "
                + ("identical" if row["identical"] else "CHANGED"),
                flush=True,
            )
    if args.output:
        args.output.write_text(json.dumps(report, indent=1) + "\n")


def project(args):
    measured = json.loads(args.pair.read_text())["rows"]
    comparison = json.loads((EVIDENCE / "julia-comparison.json").read_text())
    rows, gaps = [], []
    for workload in WORKLOADS:
        timings = comparison["timings"][workload]
        julia = min(
            (v["median_seconds"], k)
            for k, v in timings.items()
            if k.startswith("Julia/")
        )
        best = None
        for method, wait in METHODS:
            row = measured.get(f"{workload}/{method}/{wait}")
            if row is None:
                continue
            seconds = timings[REPORT_ENGINE[method]]["median_seconds"] / row["speedup"]
            if best is None or seconds < best[0]:
                best = (seconds, method, row["speedup"])
        if best is None:
            continue
        gap = best[0] / julia[0]
        gaps.append(gap)
        rows.append(
            f"| {workload} | {best[1]} | {best[2]:.2f}× | {best[0]:.4f} | "
            f"{julia[1]} | {julia[0]:.4f} | {gap:.2f}× |"
        )
    print(
        "| Workload | Fastest projected CRN | Measured speedup | Projected seconds "
        "| Fastest Julia | Julia seconds | Projected CRN / Julia |"
    )
    print("|---|---|---|---|---|---|---|")
    print("\n".join(rows))
    if gaps:
        geometric = math.exp(sum(map(math.log, gaps)) / len(gaps))
        print(f"\nGeometric mean projected gap: {geometric:.2f}×")


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    commands = parser.add_subparsers(dest="command", required=True)
    paired = commands.add_parser("pair")
    paired.add_argument("--before", type=Path, required=True)
    paired.add_argument("--after", type=Path, required=True)
    paired.add_argument("--workloads", nargs="+", default=list(WORKLOADS))
    paired.add_argument("--methods", nargs="+")
    paired.add_argument("--repeats", type=int, default=3)
    paired.add_argument("--rounds", type=int, default=5)
    paired.add_argument("--output", type=Path)
    projected = commands.add_parser("project")
    projected.add_argument("pair", type=Path)
    args = parser.parse_args()
    {"pair": pair, "project": project}[args.command](args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
