"""Same-host mtDNA benchmarks. Java inputs and outputs stay outside the patch."""

import argparse
import hashlib
import json
from pathlib import Path
import platform
import re
import subprocess
import tempfile
import time

import numpy as np

from check_mtdna import cdf_distance, read_model, require, rules_run

GRID = list(range(0, 25921, 720))


def run_rules(binary, folder, rows, model, label, n, workers, method, **kwargs):
    start = time.perf_counter()
    code, result = rules_run(
        binary,
        folder,
        model,
        method=method,
        n=n,
        threads=workers,
        times=kwargs.pop("times", GRID),
        keep=kwargs.pop("keep", False),
        **kwargs,
    )
    seconds = time.perf_counter() - start
    require(
        code == 0 and result["termination"] == "completed",
        str(result.get("diagnostics")),
    )
    events = result["events"]
    row = {
        "case": label,
        "method": method or "tree",
        "cells": n,
        "threads": workers,
        "seconds": seconds,
        "events": events,
        "events_per_second": events / seconds,
        "stopped": result["stopped_replicates"],
    }
    rows.append(row)
    print(json.dumps(row), flush=True)
    return result


def run_java(args, folder, rows, model, method, n):
    output = folder / "java.txt"
    command = [
        str(args.java),
        "-Xmx3g",
        "-Djava.awt.headless=true",
        "-XX:ActiveProcessorCount=20",
        "-cp",
        args.java_classpath,
        "replViaTrans.replViaTrans",
        model,
        "0.00003",
        "0.2",
        "1.5",
        str(n),
        "17",
        str(output),
        method,
    ]
    start = time.perf_counter()
    result = subprocess.run(
        command, capture_output=True, text=True, check=False, timeout=900
    )
    seconds = time.perf_counter() - start
    require(result.returncode == 0, result.stderr)
    text = output.read_text()
    keys = [
        "nTakenOver",
        "nCOXneg",
        "swMutationEvents",
        "smMutationEvents",
        "swMutantsVanished",
        "smMutantsVanished",
        "nReplications",
    ]
    totals = {k: int(re.search(rf"^{k}:\s+(\d+)", text, re.M)[1]) for k in keys}
    passages = []
    for values in re.findall(r"^ATP\d+;(.*)$", text, re.M):
        # Variable ATP is printed as integer. Its zero bins can begin when
        # ATP is positive but below one, so use these bins descriptively.
        passage = next((i for i, x in enumerate(values.split(";")) if x == "0"), None)
        passages.append(passage)
    require(len(passages) == n, f"Java returned {len(passages)} paths, expected {n}")
    row = {
        "case": model,
        "method": "java_" + method,
        "cells": n,
        "threads": 20,
        "seconds": seconds,
        "totals": totals,
        "output_bytes": output.stat().st_size,
    }
    rows.append(row)
    print(json.dumps(row), flush=True)
    return totals, np.array([np.inf if x is None else x for x in passages])


def fixed_cross_check(result, java, bins):
    n = len(result["replicates"])
    stops = np.array(
        [
            np.inf
            if r["termination"] != "stopped"
            else np.searchsorted(GRID, r["stop_time"], side="left")
            for r in result["replicates"]
        ]
    )
    counts = np.sum([r["rule_events"] for r in result["replicates"]], axis=0)
    rules = dict(zip(result["rules"], map(int, counts), strict=True))
    # DKW for two independent empirical passage distributions, including
    # the never-stopped mass and Java's report-time observation granularity.
    bound = 2 * np.sqrt(np.log(4 / 0.001) / (2 * n))
    distance = cdf_distance(stops, bins)
    return {
        "takeover_crn": result["stopped_replicates"],
        "takeover_java": java["nTakenOver"],
        "report_binned_passage_cdf_distance": distance,
        "bound": float(bound),
        "passed": bool(distance <= bound),
        "mutations_crn": [rules["mutation_mt"], rules["mutation_mtAdv"]],
        "mutations_java": [java["swMutationEvents"], java["smMutationEvents"]],
        "clone_extinctions_java": [
            java["swMutantsVanished"],
            java["smMutantsVanished"],
        ],
        "class_extinctions_crn": [
            sum(ob["trackers"][i]["extinctions"] for ob in result["observations"])
            for i in range(2)
        ],
    }


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
        model["rules"].extend(
            [
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
        )
    return model


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--java", type=Path)
    parser.add_argument("--java-classpath")
    parser.add_argument("--java-source", type=Path)
    args = parser.parse_args()
    quota = Path("/sys/fs/cgroup/cpu.max")
    report = {
        "platform": platform.platform(),
        "crn_binary_sha256": hashlib.sha256(args.binary.read_bytes()).hexdigest(),
        "cpu_quota": quota.read_text().strip() if quota.exists() else None,
        "rows": [],
        "complete": False,
    }
    if args.java_source:
        report["java_source_sha256"] = hashlib.sha256(
            args.java_source.read_bytes()
        ).hexdigest()
    rows = report["rows"]

    def save():
        args.output.write_text(json.dumps(report, indent=2) + "\n")

    try:
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)
            model = read_model("mtdna-three-class.rules.json")
            for workers in (1, 20):
                for method in (None, "static_direct", "next_reaction"):
                    run_rules(
                        args.binary,
                        folder,
                        rows,
                        model,
                        "fixed_no_stop",
                        200,
                        workers,
                        method,
                    )
                    save()
            for method in ("static_direct", "next_reaction"):
                run_rules(
                    args.binary,
                    folder,
                    rows,
                    sparse(),
                    "sparse_48",
                    64,
                    1,
                    method,
                    times=[0.0, 100.0],
                )
                save()
            from crn_craft import rules

            best = min(
                (
                    r
                    for r in rows
                    if r["case"] == "fixed_no_stop" and r["threads"] == 20
                ),
                key=lambda r: r["seconds"],
            )["method"]
            model_bytes = rules.load(
                Path(__file__).resolve().parents[1]
                / "models/rules/mtdna-three-class.rules.json"
            )
            hashes = []
            for threads in (20, 1):
                prepared = rules.ensemble(
                    model_bytes,
                    end=GRID[-1],
                    report_times=GRID,
                    replicates=10_000,
                    seed=17,
                    threads=threads,
                    method=best,
                )
                start = time.perf_counter()
                out = prepared.simulate()
                seconds = time.perf_counter() - start
                require(
                    out.metadata["termination"] == "completed",
                    "Python full ensemble incomplete",
                )
                signature = hashlib.sha256()
                signature.update(
                    json.dumps(
                        {k: v for k, v in out.metadata.items() if k != "run_sha256"},
                        sort_keys=True,
                    ).encode()
                )
                for key, values in sorted(out.records.items()):
                    signature.update(key.encode())
                    signature.update(values.tobytes())
                hashes.append(signature.hexdigest())
                row = {
                    "case": "fixed_no_stop_python",
                    "method": best,
                    "cells": 10000,
                    "threads": threads,
                    "seconds": seconds,
                    "events": out.metadata["events"],
                    "events_per_second": out.metadata["events"] / seconds,
                    "numerical_sha256": hashes[-1],
                }
                rows.append(row)
                print(json.dumps(row), flush=True)
                save()
            require(hashes[0] == hashes[1], "full 1/20 Python replay failed")
            report["full_python_replay"] = True

            small = read_model("mtdna-three-class-cme.rules.json")
            for label, options, callback in (
                ("bare", {}, None),
                (
                    "tracker_log",
                    {
                        "trackers": [{"observable": "mutants"}],
                        "event_log": {"max_entries": 1024},
                    },
                    None,
                ),
                ("live_callback", {}, lambda _id, _event: None),
            ):
                prepared = rules.ensemble(
                    small,
                    end=10.0,
                    report_times=[0.0, 10.0],
                    replicates=64,
                    seed=17,
                    threads=1,
                    method="static_direct",
                    **options,
                )
                observations = []
                for _ in range(3):
                    start = time.perf_counter()
                    out = prepared.simulate(callback=callback)
                    observations.append(time.perf_counter() - start)
                row = {
                    "case": "observer_overhead",
                    "method": label,
                    "cells": 64,
                    "threads": 1,
                    "seconds": float(np.median(observations)),
                    "observations_seconds": observations,
                    "events": out.metadata["events"],
                }
                rows.append(row)
                print(json.dumps(row), flush=True)
                save()
            if args.java:
                model["stop"] = "N > 0 && mutants >= 0.9 * N"
                model["observables"]["neutral_mutants"] = {
                    "type": "mtDNA",
                    "where": "x.clone == 1",
                }
                crn = run_rules(
                    args.binary,
                    folder,
                    rows,
                    model,
                    "fixed_cox_stop",
                    10000,
                    20,
                    "static_direct",
                    trackers=[
                        {"observable": "neutral_mutants"},
                        {"observable": "advantaged_mutants"},
                    ],
                    seed=17,
                    keep=True,
                )
                save()
                totals, bins = run_java(args, folder, rows, "fixedPop", "ssa", 10000)
                report["fixed_cross_check"] = fixed_cross_check(crn, totals, bins)
                save()
                run_java(args, folder, rows, "fixedPop", "hourly", 10000)
                save()
                variable = read_model("mtdna-three-class-atp.rules.json")
                crn = run_rules(
                    args.binary,
                    folder,
                    rows,
                    variable,
                    "variable_atp",
                    512,
                    20,
                    "thinning",
                    seed=17,
                    keep=True,
                )
                totals, _ = run_java(args, folder, rows, "variablePop", "ssa", 512)
                counts = np.sum([r["rule_events"] for r in crn["replicates"]], axis=0)
                report["atp_cross_check"] = {
                    "cells": 512,
                    "takeover_crn": crn["stopped_replicates"],
                    "takeover_java": totals["nTakenOver"],
                    "mutations_crn": {
                        k: int(v)
                        for k, v in zip(crn["rules"], counts, strict=True)
                        if k.startswith("mutation")
                    },
                    "mutations_java": [
                        totals["swMutationEvents"],
                        totals["smMutationEvents"],
                    ],
                    "boundary": "Java freezes propensities over its waiting interval, CRN uses continuous thinning.",
                }
            report["complete"] = True
    finally:
        save()


if __name__ == "__main__":
    main()
