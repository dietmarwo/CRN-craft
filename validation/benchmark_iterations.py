#!/usr/bin/env python3
"""Paired timings for P11–P20, as declared in docs/campaign/performance-p11-p20.md.

Before/after binaries alternate on frozen workloads. For exact iterations every
numerical digest must match before a timing ratio counts. The decision follows
the declaration: median paired speedup >= 1.03 on a target workload with a
bootstrap 95% interval above 1.0, and no W1–W4 interval entirely below 0.97.
"""

import argparse
import hashlib
import json
from pathlib import Path
import random
import statistics
import subprocess
import tempfile
import time

from benchmark_campaign import models

SEED = 4242
REPLICATES = 64
FROZEN = ["W1", "W2", "W3", "W4"]


def encode(value):
    return (json.dumps(value) + "\n").encode()


def ssa_workload(model_name, horizon, saves, keep):
    model = dict(models())[model_name]
    return {
        "kind": "ssa",
        "model": model,
        "horizon": horizon,
        "saves": saves,
        "keep": keep,
    }


def fsp_birth_death():
    return {
        "schema": "crn-craft.model/0.1",
        "name": "FspCost",
        "convention": "stochastic_counts",
        "omega": 1.0,
        "parameters": {},
        "species": [{"name": "X", "initial_count": 100}],
        "reactions": [
            {
                "name": "birth",
                "reactants": {},
                "products": {"X": 1},
                "rate": {"kind": "mass_action", "coefficient": 100.0},
            },
            {
                "name": "death",
                "reactants": {"X": 1},
                "products": {},
                "rate": {"kind": "mass_action", "coefficient": 1.0},
            },
        ],
    }


def fsp_two_species():
    reactions = []
    for name in ["A", "B"]:
        reactions += [
            {
                "name": "birth_" + name,
                "reactants": {},
                "products": {name: 1},
                "rate": {"kind": "mass_action", "coefficient": 20.0},
            },
            {
                "name": "death_" + name,
                "reactants": {name: 1},
                "products": {},
                "rate": {"kind": "mass_action", "coefficient": 0.5},
            },
        ]
    return {
        "schema": "crn-craft.model/0.1",
        "name": "FspTwoSpecies",
        "convention": "stochastic_counts",
        "omega": 1.0,
        "parameters": {},
        "species": [
            {"name": "A", "initial_count": 0},
            {"name": "B", "initial_count": 0},
        ],
        "reactions": reactions,
    }


def tau_model():
    return {
        "schema": "crn-craft.model/0.1",
        "name": "TauLargeMean",
        "convention": "stochastic_counts",
        "omega": 1.0,
        "parameters": {},
        "species": [{"name": "X", "initial_count": 10000}],
        "reactions": [
            {
                "name": "birth",
                "reactants": {},
                "products": {"X": 1},
                "rate": {"kind": "mass_action", "coefficient": 10000.0},
            },
            {
                "name": "death",
                "reactants": {"X": 1},
                "products": {},
                "rate": {"kind": "mass_action", "coefficient": 1.0},
            },
        ],
    }


DENSE = [round(0.05 * i, 10) for i in range(2001)]
WORKLOADS = {
    "W1": (ssa_workload("two_channel", 100.0, [0.0, 100.0], 0), 7),
    "W2": (ssa_workload("sparse_66", 100.0, [0.0, 100.0], 0), 7),
    "W3": (ssa_workload("circuit_12", 100.0, [0.0, 100.0], 0), 7),
    "W4": (ssa_workload("dense_66", 100.0, [0.0, 100.0], 0), 7),
    "W5": (ssa_workload("circuit_12", 100.0, DENSE, 0), 7),
    "W6": (ssa_workload("circuit_12", 100.0, DENSE, REPLICATES), 7),
    "W7": (ssa_workload("two_channel", 0.1, [0.0, 0.1], 0), 21),
    "W8a": (
        {
            "kind": "fsp",
            "model": fsp_birth_death(),
            "method": {
                "kind": "fsp",
                "bounds": {"X": 400},
                "max_states": 401,
                "tolerance": 1e-12,
            },
            "horizon": 40.0,
        },
        7,
    ),
    "W8b": (
        {
            "kind": "fsp",
            "model": fsp_two_species(),
            "method": {
                "kind": "fsp",
                "bounds": {"A": 99, "B": 99},
                "max_states": 10000,
                "tolerance": 1e-12,
            },
            "horizon": 5.0,
        },
        5,
    ),
    "W9": ({"kind": "tau", "model": tau_model(), "epsilon": 0.05, "horizon": 5.0}, 5),
}


def prepare(folder, workload):
    """Write the model and run/request; return the command arguments after the binary."""
    m = encode(workload["model"])
    (folder / "model.json").write_bytes(m)
    reference = {"path": "model.json", "sha256": hashlib.sha256(m).hexdigest()}
    if workload["kind"] == "fsp":
        run = {
            "schema": "crn-craft.run/0.1",
            "model": reference,
            "method": workload["method"],
            "time": {
                "end": workload["horizon"],
                "save_times": [0.0, workload["horizon"]],
            },
            "outputs": {"statistics": ["mean"], "keep_trajectories": 0},
            "budget": {"threads": 1, "max_output_bytes": 100000000},
        }
        (folder / "run.json").write_text(json.dumps(run))
        return ["fsp", str(folder / "model.json"), str(folder / "run.json")]
    if workload["kind"] == "tau":
        run = {
            "schema": "crn-craft.run/0.1",
            "model": reference,
            "method": {"kind": "ssa", "aggregator": "direct"},
            "ensemble": {"replicates": 256, "root_seed": SEED},
            "time": {
                "end": workload["horizon"],
                "save_times": [0.0, workload["horizon"]],
            },
            "outputs": {"statistics": ["mean", "variance"], "keep_trajectories": 256},
            "budget": {
                "threads": 1,
                "max_events_per_replicate": 100000000,
                "max_output_bytes": 100000000,
            },
        }
        request = {
            "schema": "crn-craft.tau-request/0.1",
            "epsilon": workload["epsilon"],
            "max_leaps_per_replicate": 1000000,
            "reference_run": run,
        }
        (folder / "tau.json").write_text(json.dumps(request))
        return ["tau", str(folder / "model.json"), str(folder / "tau.json")]
    run = {
        "schema": "crn-craft.run/0.1",
        "model": reference,
        "method": {"kind": "ssa", "aggregator": "direct"},
        "ensemble": {"replicates": REPLICATES, "root_seed": SEED},
        "time": {"end": workload["horizon"], "save_times": workload["saves"]},
        "outputs": {
            "statistics": ["mean", "variance"],
            "keep_trajectories": workload["keep"],
        },
        "budget": {
            "threads": 1,
            "max_events_per_replicate": 10000000,
            "max_output_bytes": 200000000,
        },
    }
    (folder / "run.json").write_text(json.dumps(run))
    return [str(folder / "model.json"), str(folder / "run.json")]


def run_once(binary, command, workload, folder, label):
    arguments = prepare(folder, workload)
    if workload["kind"] == "ssa":
        output = folder / label
        arguments = [command, *arguments, str(output)]
    start = time.perf_counter()
    proc = subprocess.run(
        [str(binary), *arguments], capture_output=True, check=False, timeout=900
    )
    elapsed = time.perf_counter() - start
    if workload["kind"] == "ssa":
        if proc.returncode:
            raise RuntimeError(proc.stdout.decode() + proc.stderr.decode())
        result = (output / "result.json").read_bytes()
        provenance = json.loads((output / "provenance.json").read_text())
        return {
            "process_seconds": elapsed,
            "core_seconds": provenance["elapsed_seconds"],
            "events": json.loads(result)["events"],
            "numerical_sha256": hashlib.sha256(result).hexdigest(),
        }
    if proc.returncode:
        raise RuntimeError(proc.stdout.decode() + proc.stderr.decode())
    payload = json.loads(proc.stdout)
    # Tau results carry no timer; FSP has none either. Digest the full payload.
    return {
        "process_seconds": elapsed,
        "core_seconds": elapsed,
        "events": payload.get("events", payload.get("work")),
        "numerical_sha256": hashlib.sha256(proc.stdout).hexdigest(),
    }


def interval(ratios, resamples=20000, seed=20261101):
    rng = random.Random(seed)
    medians = sorted(
        statistics.median(rng.choices(ratios, k=len(ratios))) for _ in range(resamples)
    )
    return [medians[int(0.025 * resamples)], medians[int(0.975 * resamples) - 1]]


def measure(name, before, after, command, pairs, exact):
    workload, default_pairs = WORKLOADS[name]
    pairs = pairs or default_pairs
    rows = {"before": [], "after": []}
    with tempfile.TemporaryDirectory(prefix="crn-iter-") as temp:
        folder = Path(temp)
        for pair in range(pairs):
            for label, binary in [("before", before), ("after", after)]:
                rows[label].append(
                    run_once(binary, command, workload, folder, f"{label}{pair}")
                )
    digests = {r["numerical_sha256"] for r in rows["before"] + rows["after"]}
    ratios = [
        b["core_seconds"] / a["core_seconds"]
        for b, a in zip(rows["before"], rows["after"], strict=True)
    ]
    entry = {
        "pairs": pairs,
        "identical": len(digests) == 1,
        "events": rows["before"][0]["events"],
        "before_median_core_seconds": statistics.median(
            r["core_seconds"] for r in rows["before"]
        ),
        "after_median_core_seconds": statistics.median(
            r["core_seconds"] for r in rows["after"]
        ),
        "paired_speedup_median": statistics.median(ratios),
        "paired_speedup_interval_95": interval(ratios),
        "rows": rows,
    }
    if workload["kind"] == "ssa":
        entry["before_ns_per_event"] = 1e9 * statistics.median(
            r["core_seconds"] / r["events"] for r in rows["before"]
        )
        entry["after_ns_per_event"] = 1e9 * statistics.median(
            r["core_seconds"] / r["events"] for r in rows["after"]
        )
    if exact and not entry["identical"]:
        raise RuntimeError(f"{name}: numerical payload changed")
    return entry


def decide(report, targets):
    workloads = report["workloads"]
    gained = [
        t
        for t in targets
        if workloads[t]["paired_speedup_median"] >= 1.03
        and workloads[t]["paired_speedup_interval_95"][0] > 1.0
    ]
    regressed = [
        w
        for w in FROZEN
        if w in workloads and workloads[w]["paired_speedup_interval_95"][1] < 0.97
    ]
    return {
        "targets": targets,
        "targets_gained": gained,
        "frozen_regressions": regressed,
        "accepted": bool(gained) and not regressed,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--iteration", required=True)
    ap.add_argument("--before", type=Path, required=True)
    ap.add_argument("--after", type=Path, required=True)
    ap.add_argument("--command", default="simulate-cached")
    ap.add_argument("--workloads", nargs="+", required=True, choices=list(WORKLOADS))
    ap.add_argument("--targets", nargs="+", required=True, choices=list(WORKLOADS))
    ap.add_argument("--pairs", type=int)
    ap.add_argument("--approximate", action="store_true")
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    report = {
        "schema": "crn-craft.performance-iteration/0.1",
        "iteration": args.iteration,
        "declaration": "docs/campaign/performance-p11-p20.md",
        "command": args.command,
        "seed": SEED,
        "exact": not args.approximate,
        "workloads": {},
    }
    for name in args.workloads:
        report["workloads"][name] = measure(
            name,
            args.before.resolve(),
            args.after.resolve(),
            args.command,
            args.pairs,
            not args.approximate,
        )
        entry = report["workloads"][name]
        print(
            name,
            f"{entry['paired_speedup_median']:.3f}",
            [round(x, 3) for x in entry["paired_speedup_interval_95"]],
            "identical" if entry["identical"] else "changed",
            flush=True,
        )
        args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    report["decision"] = decide(report, args.targets)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report["decision"]))


if __name__ == "__main__":
    main()
