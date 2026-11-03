"""Replay and alternating benchmark for the Static Direct optimization campaign."""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory
from time import perf_counter

from benchmark_mtdna import GRID, sparse
from check_mtdna import read_model, require, rules_files

TRACKERS = [{"observable": x} for x in ("neutral_mutants", "advantaged_mutants")]
FAST = {"kind": "static_direct", "waiting_time": "ziggurat", "sparse_trackers": True}


def mt(stop=True):
    model = read_model("mtdna-three-class.rules.json")
    if stop:
        model["stop"] = "N > 0 && mutants >= 0.9 * N"
    return model


def invoke(binary, files):
    start = perf_counter()
    result = subprocess.run(
        [binary, "rules", *files], capture_output=True, check=False, timeout=300
    )
    return perf_counter() - start, result


def replay(folder, baseline, candidate):
    model = read_model("immigration-death.rules.json")
    stiff = read_model("immigration-death.rules.json")
    stiff["parameters"]["k_in"] = 1e-30
    stiff["initial"] = [{"type": "X", "count": 10}]
    aggregate = read_model("aggregation-constant.rules.json")
    aggregate["rules"][0]["products"][0]["attributes"]["size"] = 1
    overflow = read_model("immigration-death.rules.json")
    overflow["observables"]["cold"] = {"type": "X", "weight": "1e308"}
    undefined = read_model("immigration-death.rules.json")
    undefined["observables"]["cold"] = {"type": "X", "weight": "2"}
    undefined["stop"] = "if(X > 0, 0/0, 0)"
    small = read_model("mtdna-three-class-cme.rules.json")
    small["stop"] = "N > 0 && mutants >= 0.5 * N"
    options = [
        ("mt", mt(False), {}),
        ("stop", mt(), {}),
        ("trackers", mt(), {"trackers": TRACKERS}),
        (
            "predicate",
            mt(),
            {
                "trackers": [
                    {
                        "observable": "mutants",
                        "first_passage": "N > 0 && mutants >= 0.5 * N",
                    }
                ]
            },
        ),
        ("log", mt(), {"log": {"max_entries": 4}}),
        ("sparse", sparse(), {}),
        ("aggregate", aggregate, {}),
        ("immigration", model, {}),
        ("next", mt(), {"method": "next_reaction"}),
        (
            "thinning",
            read_model("mtdna-three-class-atp.rules.json"),
            {"method": "thinning"},
        ),
        ("tau", mt(), {"method": {"kind": "tau", "epsilon": 0.01}}),
        ("binomial", mt(), {"method": {"kind": "binomial", "step": 1.0}}),
        ("small-stop", small, {"trackers": TRACKERS}),
        ("cold-overflow", overflow, {}),
        ("stop-diagnostic", undefined, {}),
        ("resolution", stiff, {}),
        ("instance", sparse(), {"max_instances": 1}),
    ]
    rows = []
    for label, fixture, extra in options:
        for workers in (1, 8):
            settings = {
                "method": "static_direct",
                "n": 8,
                "times": [0, 10, 100],
                "threads": workers,
            } | extra
            files = rules_files(folder, fixture, **settings)
            _, old = invoke(baseline, files)
            _, new = invoke(candidate, files)
            require(
                old.returncode == new.returncode
                and old.stdout == new.stdout
                and old.stderr == new.stderr,
                f"replay: {label}/{workers}",
            )
            payload = json.loads(new.stdout)
            expected = {
                "resolution": "failed",
                "instance": "instance_limit",
                "cold-overflow": "failed",
                "stop-diagnostic": "failed",
            }.get(label, "completed")
            require(
                payload["termination"] == expected,
                f"fixture {label}: {payload['termination']}",
            )
            rows.append(
                {
                    "case": label,
                    "workers": workers,
                    "sha256": hashlib.sha256(new.stdout).hexdigest(),
                }
            )
    return rows


def timing(folder, left, right, stage, selected=None):
    rows = []
    for model, label, cells, times in [
        (mt(), "mtdna", 512, GRID),
        (sparse(), "sparse", 64, [0, 100]),
    ]:
        if selected is not None and label != selected:
            continue
        for workers in (1, 8) if stage == "s4" else (1,):
            method_a, method_b = "static_direct", "static_direct"
            comparison = right
            if stage == "s4":
                method_a, method_b = (
                    (FAST, "mtdna_direct")
                    if label == "mtdna"
                    else ("static_direct", "next_reaction")
                )
            settings = {
                "n": cells,
                "times": times,
                "threads": workers,
                "seed": 17,
                "trackers": TRACKERS if label == "mtdna" else (),
                "keep": label == "mtdna",
            }
            if stage == "s4" and label == "sparse":
                comparison = left
            for repeat in range(3 if stage == "s4" else 1):
                order = [(left, method_a, "left"), (comparison, method_b, "right")]
                for binary, method, side in order[:: (-1 if repeat % 2 else 1)]:
                    elapsed, result = invoke(
                        binary, rules_files(folder, model, method=method, **settings)
                    )
                    require(result.returncode == 0, result.stdout[:1000])
                    output = json.loads(result.stdout)
                    row = {
                        "case": label,
                        "workers": workers,
                        "repeat": repeat,
                        "side": side,
                        "sampler": output["sampler"],
                        "events": output["events"],
                        "seconds": elapsed,
                    }
                    rows.append(row)
                    print(json.dumps(row), flush=True)
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("s1", "s2", "s4", "replay"))
    parser.add_argument("left", type=Path)
    parser.add_argument("right", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--case", choices=("mtdna", "sparse"))
    args = parser.parse_args()
    with TemporaryDirectory() as tmp:
        folder = Path(tmp)
        rows = (
            replay(folder, args.left, args.right)
            if args.stage == "replay"
            else timing(folder, args.left, args.right, args.stage, args.case)
        )
    args.output.write_text(json.dumps(rows, indent=2) + "\n")


if __name__ == "__main__":
    main()
