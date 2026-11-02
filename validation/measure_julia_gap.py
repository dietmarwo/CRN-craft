"""Paired timing of exact SSA engines and a projected Julia gap.

Build `cargo build --release --example bench_ssa` at two commits and copy each
binary aside, then:

    python validation/measure_julia_gap.py chain base=OLD a01=NEW1 a02=NEW2 \
        --output chain.json
    python validation/measure_julia_gap.py project chain.json

`chain` times every binary in each round on W1-W7 of the acceleration report,
rotating the order between rounds, and records every timed repeat, the event
count and the result digests. Each commit is paired with its parent within
the same rounds; the first-to-last pair gives the cumulative speedup.
`project` multiplies the report's CRN prepared medians (another host) by the
cumulative speedups and divides by Julia's fastest median.
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


def chain(args):
    names = [b.split("=", 1)[0] for b in args.binaries]
    paths = {b.split("=", 1)[0]: Path(b.split("=", 1)[1]) for b in args.binaries}
    methods = [m for m in METHODS if not args.methods or m[0] in args.methods]
    report = {"protocol": "docs/decisions/julia-gap-2026-11-02.md", "order": names}
    report["rows"] = {}
    for workload in args.workloads:
        for method, wait in methods:
            seconds = {name: [] for name in names}
            digests = {name: set() for name in names}
            numerical = {name: set() for name in names}
            events = {}
            for round_index in range(args.rounds):
                # Rotate, and reverse every other round, so no commit always runs first.
                order = (
                    names[round_index % len(names) :]
                    + names[: round_index % len(names)]
                )
                if round_index % 2:
                    order.reverse()
                for name in order:
                    result = bench(paths[name], workload, method, wait, args.repeats)
                    seconds[name].append(result["seconds"])
                    digests[name].add(result["digest"])
                    numerical[name].add(
                        result.get("numerical_digest", result["digest"])
                    )
                    events[name] = result["events"]
            medians = {k: [statistics.median(r) for r in v] for k, v in seconds.items()}
            steps = {}
            for parent, child in [
                *zip(names, names[1:], strict=False),
                (names[0], names[-1]),
            ]:
                ratios = [
                    b / a for b, a in zip(medians[parent], medians[child], strict=True)
                ]
                steps[f"{parent}->{child}"] = {
                    "ratios": ratios,
                    "speedup": statistics.median(ratios),
                    "rounds_faster": sum(r > 1.0 for r in ratios),
                    "identical": digests[parent] == digests[child]
                    and len(digests[child]) == 1,
                    "numerically_identical": numerical[parent] == numerical[child]
                    and len(numerical[child]) == 1,
                }
            key = f"{workload}/{method}/{wait}"
            report["rows"][key] = {
                "seconds": seconds,
                "events": events,
                "digests": {k: sorted(v) for k, v in digests.items()},
                "numerical": {k: sorted(v) for k, v in numerical.items()},
                "steps": steps,
            }
            cells = []
            for label, step in steps.items():
                flag = (
                    ""
                    if step["identical"]
                    else "~"
                    if step["numerically_identical"]
                    else "!"
                )
                cells.append(
                    f"{label} {step['speedup']:.3f}({step['rounds_faster']}){flag}"
                )
            print(f"{key:30s} " + " ".join(cells), flush=True)
            if args.output:
                args.output.write_text(json.dumps(report, indent=1) + "\n")


def project(args):
    report = json.loads(args.chain.read_text())
    first, last = report["order"][0], report["order"][-1]
    measured = {
        key: {"speedup": row["steps"][f"{first}->{last}"]["speedup"]}
        for key, row in report["rows"].items()
    }
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
    chained = commands.add_parser("chain")
    chained.add_argument(
        "binaries", nargs="+", help="NAME=PATH, oldest commit first (at least two)"
    )
    chained.add_argument("--workloads", nargs="+", default=list(WORKLOADS))
    chained.add_argument("--methods", nargs="+")
    chained.add_argument("--repeats", type=int, default=3)
    chained.add_argument("--rounds", type=int, default=5)
    chained.add_argument("--output", type=Path)
    projected = commands.add_parser("project")
    projected.add_argument("chain", type=Path)
    args = parser.parse_args()
    {"chain": chain, "project": project}[args.command](args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
