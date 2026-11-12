"""Frozen dynamic API benchmark. Exact payload identities precede comparisons.

Since 0.7.0 (M8a) the frozen network models W1-W4 are converted with
`crn-craft convert` and run as `crn-craft.run/1.0` Static Direct requests by
`crn-craft run MODEL REQUEST` (default selector and Ziggurat waits), so both
binaries under comparison must be 0.7.0 or later. The 0.6.0 version timed the
core of `simulate-cached` from the archive's provenance file; archives of the
1.0 contract have none, so the time is the wall time of the whole process
(`process_seconds`, which includes the process start and printing the result).
"""

import argparse
import hashlib
import json
from pathlib import Path
import statistics
import subprocess
import tempfile
import time


def models():
    def base(name, species, reactions, parameters=None):
        return {
            "schema": "crn-craft.model/0.1",
            "name": name,
            "convention": "stochastic_counts",
            "omega": 1.0,
            "parameters": parameters or {},
            "species": [{"name": s, "initial_count": x} for s, x in species],
            "reactions": reactions,
        }

    def rx(n, a, b, k):
        return {
            "name": n,
            "reactants": a,
            "products": b,
            "rate": {"kind": "mass_action", "coefficient": k},
        }

    yield (
        "two_channel",
        base(
            "TwoChannel",
            [("X", 0)],
            [rx("birth", {}, {"X": 1}, 100.0), rx("death", {"X": 1}, {}, 1.0)],
        ),
    )
    species = [(f"X{i}", 0) for i in range(33)]
    reactions = []
    for s, _ in species:
        reactions += [
            rx("birth_" + s, {}, {s: 1}, 3.0),
            rx("death_" + s, {s: 1}, {}, 0.2),
        ]
    yield "sparse_66", base("Sparse66", species, reactions)
    root = Path(__file__).resolve().parents[1]
    yield (
        "circuit_12",
        json.loads((root / "models/examples/repressilator.model.json").read_text()),
    )
    species = [(f"X{i}", 20) for i in range(8)]
    both = {s: 1 for s, _ in species}
    reactions = []
    for i in range(33):
        reactions.append(rx(f"birth{i}", {}, both, 3.0))
        reactions.append(
            {
                "name": f"death{i}",
                "reactants": both,
                "products": {},
                "rate": {
                    "kind": "michaelis_menten",
                    "substrate": "X0",
                    "vmax": 3.0,
                    "km": 20.0,
                },
            }
        )
    yield "dense_66", base("Dense66", species, reactions)


def request_for(model_bytes, threads, horizon, selector=None, waiting_time=None):
    method = {"kind": "static_direct"}
    if selector is not None:
        method["selector"] = selector
    if waiting_time is not None:
        method["waiting_time"] = waiting_time
    return {
        "schema": "crn-craft.run/1.0",
        "model": {
            "path": "model.json",
            "sha256": hashlib.sha256(model_bytes).hexdigest(),
        },
        "method": method,
        "time": {"end": horizon, "report_times": [0.0, horizon]},
        "ensemble": {"replicates": 64, "root_seed": 4242},
        "outputs": {
            "statistics": ["mean", "variance"],
            "keep_paths": 0,
            "keep_instances": 0,
        },
        "budget": {
            "threads": threads,
            "max_events_per_replicate": 10000000,
            "max_sampling_steps_per_replicate": 1000000000,
            "max_instances": 100000,
            "max_output_bytes": 10000000,
        },
    }


def convert(binary, model, folder):
    """The converted `crn-craft.model/1.0` bytes of a network model dictionary."""
    (folder / "network.json").write_bytes((json.dumps(model) + "\n").encode())
    p = subprocess.run(
        [str(binary), "convert", str(folder / "network.json")],
        capture_output=True,
        check=False,
        timeout=60,
    )
    if p.returncode:
        raise RuntimeError(p.stdout.decode() + p.stderr.decode())
    return p.stdout


def run_case(
    binary,
    command,
    model,
    repetitions,
    threads=1,
    horizon=100.0,
    selector=None,
    waiting_time=None,
):
    """Time `repetitions` runs of `crn-craft COMMAND` (`run`) on a network model."""
    rows = []
    with tempfile.TemporaryDirectory(prefix="crn-bench-") as tmp:
        folder = Path(tmp)
        m = convert(binary, model, folder)
        (folder / "model.json").write_bytes(m)
        run = request_for(m, threads, horizon, selector, waiting_time)
        (folder / "run.json").write_text(json.dumps(run))
        for _ in range(repetitions):
            start = time.perf_counter()
            p = subprocess.run(
                [
                    str(binary),
                    command,
                    str(folder / "model.json"),
                    str(folder / "run.json"),
                ],
                capture_output=True,
                check=False,
                timeout=180,
            )
            elapsed = time.perf_counter() - start
            if p.returncode:
                raise RuntimeError(p.stdout.decode() + p.stderr.decode())
            result = json.loads(p.stdout)
            result.pop("run_sha256")
            digest = hashlib.sha256(
                json.dumps(result, sort_keys=True, separators=(",", ":")).encode()
            ).hexdigest()
            rows.append(
                {
                    "process_seconds": elapsed,
                    "events": result["events"],
                    "numerical_sha256": digest,
                }
            )
    return {
        "median_process_seconds": statistics.median(r["process_seconds"] for r in rows),
        "process_ns_per_event": 1e9
        * statistics.median(r["process_seconds"] / r["events"] for r in rows),
        "rows": rows,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--before", type=Path, required=True)
    ap.add_argument("--after", type=Path)
    ap.add_argument("--before-command", default="run")
    ap.add_argument("--after-command", default="run")
    ap.add_argument("--repetitions", type=int, default=7)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--only", choices=[n for n, _ in models()])
    ap.add_argument("--threads", type=int, default=1)
    args = ap.parse_args()
    report = {
        "schema": "crn-craft.performance/0.2",
        "seed": 4242,
        "replicates": 64,
        "end": 100.0,
        "threads": args.threads,
        "workloads": {},
    }
    for name, model in models():
        if args.only and name != args.only:
            continue
        # Pair alternating before/after calls to reduce time-order bias.
        before = []
        after = []
        for _ in range(args.repetitions):
            before.extend(
                run_case(
                    args.before.resolve(), args.before_command, model, 1, args.threads
                )["rows"]
            )
            if args.after:
                after.extend(
                    run_case(
                        args.after.resolve(), args.after_command, model, 1, args.threads
                    )["rows"]
                )

        def summary(rows):
            return {
                "process_range_seconds": [
                    min(r["process_seconds"] for r in rows),
                    max(r["process_seconds"] for r in rows),
                ],
                "median_process_seconds": statistics.median(
                    r["process_seconds"] for r in rows
                ),
                "process_ns_per_event": 1e9
                * statistics.median(r["process_seconds"] / r["events"] for r in rows),
                "rows": rows,
            }

        entry = {"before": summary(before)}
        if after:
            entry["after"] = summary(after)
            entry["identical"] = (
                len({r["numerical_sha256"] for r in before + after}) == 1
            )
            entry["paired_speedup_median"] = statistics.median(
                b["process_seconds"] / a["process_seconds"]
                for b, a in zip(before, after, strict=True)
            )
            if not entry["identical"]:
                raise RuntimeError(f"{name}: numerical payload changed")
        report["workloads"][name] = entry
        args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
        print(name, entry.get("paired_speedup_median"), flush=True)


if __name__ == "__main__":
    main()
