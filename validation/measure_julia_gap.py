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
import re
import statistics
import subprocess
import sys
import tempfile
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


def bench(binary, workload, method, wait, repeats, fast=False):
    horizon, paths = WORKLOADS[workload]
    model = EVIDENCE / "models" / f"{workload}.model.json"
    command = [binary, model, horizon, paths, method, wait, repeats]
    if fast:
        command.append("fast")
    output = subprocess.run(
        [str(part) for part in command], capture_output=True, text=True, check=True
    ).stdout
    return json.loads(output)


def chain(args):
    names = [b.split("=", 1)[0] for b in args.binaries]
    # NAME=PATH@fast times the same binary with method.fast (ssa-run/0.2).
    paths = {
        b.split("=", 1)[0]: Path(b.split("=", 1)[1].split("@")[0])
        for b in args.binaries
    }
    fast = {b.split("=", 1)[0]: b.endswith("@fast") for b in args.binaries}
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
                    result = bench(
                        paths[name], workload, method, wait, args.repeats, fast[name]
                    )
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


def instructions(args):
    """Callgrind instruction counts inside the path functions, first four paths."""
    report = {}
    for spec in args.binaries:
        name, binary = spec.split("=", 1)
        for workload in WORKLOADS:
            horizon = WORKLOADS[workload][0]
            model = EVIDENCE / "models" / f"{workload}.model.json"
            for method, wait in METHODS:
                with tempfile.TemporaryDirectory() as temp:
                    out = Path(temp) / "callgrind.out"
                    command = [
                        "valgrind",
                        "--tool=callgrind",
                        "--toggle-collect=crn_craft::ssa::path*",
                        "--toggle-collect=crn_craft::simulation::execute_inner::*",
                        f"--callgrind-out-file={out}",
                        binary,
                        model,
                        horizon,
                        4,
                        method,
                        wait,
                        1,
                    ]
                    subprocess.run(
                        [str(part) for part in command], capture_output=True, check=True
                    )
                    summary = re.search(r"^summary: (\d+)", out.read_text(), re.M)
                key = f"{workload}/{method}/{wait}"
                report.setdefault(name, {})[key] = int(summary.group(1))
                print(name, key, report[name][key], flush=True)
    args.output.write_text(json.dumps(report, indent=1) + "\n")


# First-round target rows, as declared in docs/decisions/julia-gap-2026-11-02.md.
NEW = ["cached_direct", "direct_cr", "sorting_direct", "rssa_cr"]
TREE = ["direct_cr", "sorting_direct", "rssa_cr"]
DECLARED = {
    "a01": (list(WORKLOADS), NEW),
    "a02": (list(WORKLOADS), ["legacy_cached", *NEW]),
    "a03": (["W3"], ["legacy_cached", *NEW]),
    "a04": (["W4", "W5", "W6", "W7"], ["rssa_cr"]),
    "a05": (list(WORKLOADS), ["direct_cr", "rssa_cr"]),
    "a06": (list(WORKLOADS), TREE),
    "a07": (list(WORKLOADS), ["legacy_cached", "cached_direct"]),
    "a09": (list(WORKLOADS), TREE),
    "a10": (list(WORKLOADS), ["direct_cr", "rssa_cr"]),
}


def decide(args):
    report = json.loads(args.chain.read_text())
    rows, order = report["rows"], report["order"]
    counts = json.loads(args.instructions.read_text()) if args.instructions else None
    # Five more rounds for a row that failed only the no-regression rule.
    extra = {}
    for path in args.remeasure:
        more = json.loads(path.read_text())
        for key, row in more["rows"].items():
            for label, step in row["steps"].items():
                extra[(label, key)] = step["ratios"]
    for parent, child in zip(order, order[1:], strict=False):
        step = f"{parent}->{child}"
        if counts is None:
            workloads, methods = DECLARED.get(child, ([], []))
            targets = [
                k
                for k in rows
                if k.split("/")[0] in workloads and k.split("/")[1] in methods
            ]
        else:
            targets = [k for k in rows if counts[child][k] <= 0.95 * counts[parent][k]]
        steps = {k: dict(rows[k]["steps"][step]) for k in rows}
        for k in rows:
            if k not in targets and (step, k) in extra:
                ratios = steps[k]["ratios"] + extra[(step, k)]
                steps[k].update(ratios=ratios, speedup=statistics.median(ratios))
        median = (
            statistics.median(steps[k]["speedup"] for k in targets) if targets else 0.0
        )
        weak = [k for k in targets if steps[k]["rounds_faster"] < 4]
        regress = [k for k in rows if k not in targets and steps[k]["speedup"] < 0.97]
        changed = [k for k in rows if not steps[k]["numerically_identical"]]
        if child in args.path_changing:
            changed = []  # decided by the statistical gate instead
        performance = child not in args.no_target
        accepted = (
            not regress
            and not changed
            and (not performance or (targets and median >= 1.03 and not weak))
        )
        print(
            f"{step:10s} targets {len(targets):2d} median {median:.3f} "
            f"weak {len(weak)} regressions {len(regress)} changed {len(changed)} "
            f"-> {'ACCEPT' if accepted else 'REJECT'}"
        )
        for k in targets:
            print(
                f"   target {k:30s} {steps[k]['speedup']:.3f} ({steps[k]['rounds_faster']})"
            )
        for k in regress:
            ratios = [round(r, 3) for r in steps[k]["ratios"]]
            print(f"   regression {k:26s} {steps[k]['speedup']:.3f} {ratios}")
        for k in changed:
            print(f"   numbers changed {k}")


def project(args):
    report = json.loads(args.chain.read_text())
    first, last = report["order"][0], report["order"][-1]
    measured = {
        key: {"speedup": row["steps"][f"{first}->{last}"]["speedup"]}
        for key, row in report["rows"].items()
    }
    comparison = json.loads(args.reference.read_text())
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
        "binaries",
        nargs="+",
        help="NAME=PATH or NAME=PATH@fast, oldest commit first (at least two)",
    )
    chained.add_argument("--workloads", nargs="+", default=list(WORKLOADS))
    chained.add_argument("--methods", nargs="+")
    chained.add_argument("--repeats", type=int, default=3)
    chained.add_argument("--rounds", type=int, default=5)
    chained.add_argument("--output", type=Path)
    counted = commands.add_parser("instructions")
    counted.add_argument("binaries", nargs="+", help="NAME=PATH")
    counted.add_argument("--output", type=Path, required=True)
    decided = commands.add_parser("decide")
    decided.add_argument("chain", type=Path)
    decided.add_argument(
        "--instructions", type=Path, help="B re-test: targets from instruction counts"
    )
    decided.add_argument(
        "--no-target", nargs="*", default=[], help="non-performance commits (fixes)"
    )
    decided.add_argument(
        "--path-changing", nargs="*", default=[], help="commits gated statistically"
    )
    decided.add_argument(
        "--remeasure", nargs="*", type=Path, default=[], help="extra chain reports"
    )
    projected = commands.add_parser("project")
    projected.add_argument("chain", type=Path)
    projected.add_argument(
        "--reference",
        type=Path,
        default=EVIDENCE / "julia-comparison.json",
        help="same-host CRN/Julia timing report whose CRN medians are scaled",
    )
    args = parser.parse_args()
    commands = {
        "chain": chain,
        "instructions": instructions,
        "decide": decide,
        "project": project,
    }
    commands[args.command](args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
