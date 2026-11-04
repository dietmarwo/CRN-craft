"""Time five rules workloads on two CLI binaries and compare their outputs.

Used for the structural iterations in
docs/decisions/generic-rules-structural-2026-11-03.md. Each workload runs
REPEATS times per binary in alternating order. The report keeps every
timing, the median, the events of both sides, nanoseconds per event and
whether the two outputs are byte-identical. Timing is host-specific.
"""

import argparse
import hashlib
import json
import platform
import subprocess
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
GRID = [float(t) for t in range(0, 25921, 720)]
TRACK = [{"observable": "neutral_mutants"}, {"observable": "advantaged_mutants"}]
FAST = {"kind": "static_direct", "waiting_time": "ziggurat", "sparse_trackers": True}


def three_class():
    path = ROOT / "models/rules/mtdna-three-class.rules.json"
    model = json.loads(path.read_text())
    model["stop"] = "N > 0 && mutants >= 0.9 * N"
    return model


def sparse():
    model = {
        "schema": "crn-craft.rules-model/0.1",
        "name": "sparse_48",
        "types": [],
        "initial": [],
        "rules": [],
    }
    for i in range(48):
        name = f"X{i}"
        model["types"].append({"name": name})
        model["initial"].append({"type": name, "count": 8})
        model["rules"] += [
            {
                "name": f"birth{i}",
                "reactants": [],
                "rate": "2",
                "products": [{"type": name}],
            },
            {
                "name": f"death{i}",
                "reactants": [{"type": name, "label": "x"}],
                "rate": "0.2",
                "products": [],
            },
        ]
    return model


# label: (model, request schema, method, replicates, end, report times, trackers)
WORK = {
    "mt/generic": (
        three_class,
        "crn-craft.rules-run/0.4",
        FAST,
        128,
        25920.0,
        GRID,
        TRACK,
    ),
    "mt/default": (
        three_class,
        "crn-craft.rules-run/0.2",
        {"kind": "static_direct"},
        128,
        25920.0,
        GRID,
        TRACK,
    ),
    "mt/special": (
        three_class,
        "crn-craft.rules-run/0.3",
        {"kind": "mtdna_direct"},
        128,
        25920.0,
        GRID,
        TRACK,
    ),
    "sparse/default": (
        sparse,
        "crn-craft.rules-run/0.2",
        {"kind": "static_direct"},
        64,
        100.0,
        [0.0, 50.0, 100.0],
        (),
    ),
    "sparse/next": (
        sparse,
        "crn-craft.rules-run/0.2",
        {"kind": "next_reaction"},
        64,
        100.0,
        [0.0, 50.0, 100.0],
        (),
    ),
}


def request(data, schema, method, n, end, times, trackers):
    run = {
        "schema": schema,
        "model": {
            "path": "m.rules.json",
            "sha256": hashlib.sha256(data).hexdigest(),
        },
        "time": {"end": end, "report_times": times},
        "ensemble": {"replicates": n, "root_seed": 17},
        "outputs": {
            "statistics": ["mean", "variance"],
            "keep_paths": n,
            "keep_instances": 0,
        },
        "budget": {
            "threads": 1,
            "max_events_per_replicate": 10_000_000,
            "max_sampling_steps_per_replicate": 100_000_000,
            "max_instances": 100_000,
            "max_output_bytes": 1 << 30,
        },
        "method": method,
    }
    if trackers:
        run["trackers"] = list(trackers)
    return run


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference", type=Path, help="reference crn-craft binary")
    parser.add_argument("candidate", type=Path, help="candidate crn-craft binary")
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    binaries = {"ref": args.reference.resolve(), "new": args.candidate.resolve()}
    report = {
        "schema": "crn-craft.structural-probe/0.1",
        "platform": platform.platform(),
        "repeats": args.repeats,
        "rows": {},
    }
    with tempfile.TemporaryDirectory() as tmp:
        folder = Path(tmp)
        for label, (make, schema, method, n, end, times, trackers) in WORK.items():
            data = json.dumps(make()).encode()
            (folder / "m.rules.json").write_bytes(data)
            run = request(data, schema, method, n, end, times, trackers)
            (folder / "r.json").write_text(json.dumps(run))
            timings = {"ref": [], "new": []}
            outputs = {}
            for repeat in range(args.repeats):
                order = ("ref", "new") if repeat % 2 == 0 else ("new", "ref")
                for side in order:
                    start = time.perf_counter()
                    done = subprocess.run(
                        [str(binaries[side]), "rules", "m.rules.json", "r.json"],
                        cwd=folder,
                        capture_output=True,
                        check=True,
                    )
                    timings[side].append(time.perf_counter() - start)
                    outputs[side] = done.stdout
            median = {k: sorted(v)[args.repeats // 2] for k, v in timings.items()}
            events = {k: json.loads(v)["events"] for k, v in outputs.items()}
            per_event = {k: 1e9 * median[k] / events[k] for k in median}
            report["rows"][label] = {
                "seconds": timings,
                "median": median,
                "events": events,
                "ns_per_event": per_event,
                "identical": outputs["ref"] == outputs["new"],
            }
            print(
                f"{label:15s} {per_event['ref']:.0f} -> {per_event['new']:.0f} ns/event "
                f"(x{per_event['ref'] / per_event['new']:.2f}) "
                f"identical={outputs['ref'] == outputs['new']}",
                flush=True,
            )
    rows = report["rows"]
    report["gap_generic_vs_special"] = (
        rows["mt/generic"]["ns_per_event"]["new"]
        / rows["mt/special"]["ns_per_event"]["new"]
    )
    print(f"generic/special gap: {report['gap_generic_vs_special']:.2f}")
    args.output.write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
