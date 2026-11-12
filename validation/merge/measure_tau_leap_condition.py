"""Cost of the rules `tau` leap condition on W5 and W6 (M5c, recorded).

`tau` checks every channel's rate after one firing of every enabled channel
before a leap, O(channels²) per attempted leap. This runs the converted W5
(510 channels) and W6 (992 channels) with `tau` (ε = 0.05) and with
`static_direct`, as declared in validation/merge/RESULTS.md (M5c): end 20,
report times 0 and 20, 16 replicates, root seed 4242, 1 thread, the wall time
of the `crn-craft rules` process, medians of 3.
"""

import argparse
import hashlib
import json
from pathlib import Path
import statistics
import subprocess
import tempfile
import time

ROOT = Path(__file__).resolve().parents[2]
MODELS = ROOT / "validation/ssa-acceleration/models"
METHODS = {
    "tau": {"kind": "tau", "epsilon": 0.05},
    "static_direct": {"kind": "static_direct"},
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--repetitions", type=int, default=3)
    parser.add_argument("--replicates", type=int, default=16)
    args = parser.parse_args()
    binary = str(args.binary.resolve())
    report = dict(
        schema="crn-craft.tau-leap-condition-cost/0.1",
        declared_in="validation/merge/RESULTS.md (M5c declaration)",
        replicates=args.replicates,
        repetitions=args.repetitions,
        rows=[],
    )
    with tempfile.TemporaryDirectory(prefix="crn-tau-cost-") as temp:
        folder = Path(temp)
        for name in ["W5", "W6"]:
            converted = subprocess.run(
                [binary, "convert", str(MODELS / f"{name}.model.json")],
                capture_output=True,
                check=True,
            ).stdout
            (folder / "model.json").write_bytes(converted)
            channels = len(json.loads(converted)["rules"])
            for label, method in METHODS.items():
                request = {
                    "schema": "crn-craft.run/1.0",
                    "model": {
                        "path": "model.json",
                        "sha256": hashlib.sha256(converted).hexdigest(),
                    },
                    "method": method,
                    "ensemble": {"replicates": args.replicates, "root_seed": 4242},
                    "time": {"end": 20.0, "report_times": [0.0, 20.0]},
                    "outputs": {
                        "statistics": ["mean"],
                        "keep_paths": 0,
                        "keep_instances": 0,
                    },
                    "budget": {
                        "threads": 1,
                        "max_events_per_replicate": 1000000000,
                        "max_sampling_steps_per_replicate": 1000000000000,
                        "max_instances": 1000000,
                        "max_output_bytes": 100000000,
                    },
                }
                (folder / "request.json").write_text(json.dumps(request))
                walls = []
                result = None
                for _ in range(args.repetitions):
                    start = time.perf_counter()
                    proc = subprocess.run(
                        [
                            binary,
                            "rules",
                            str(folder / "model.json"),
                            str(folder / "request.json"),
                        ],
                        capture_output=True,
                        check=False,
                    )
                    walls.append(time.perf_counter() - start)
                    if proc.returncode:
                        raise RuntimeError(proc.stdout.decode()[:2000])
                    result = json.loads(proc.stdout)
                row = dict(
                    model=name,
                    channels=channels,
                    method=label,
                    request=method,
                    sampler=result["sampler"],
                    termination=result["termination"],
                    events_per_replicate=result["events"] / args.replicates,
                    rng_words_per_replicate=sum(
                        r["rng_words"] for r in result["replicates"]
                    )
                    / args.replicates,
                    wall_seconds=statistics.median(walls),
                    wall_samples=walls,
                )
                report["rows"].append(row)
                args.output.write_text(json.dumps(report, indent=2) + "\n")
                print(name, label, row["wall_seconds"], flush=True)
    for name in ["W5", "W6"]:
        by = {r["method"]: r for r in report["rows"] if r["model"] == name}
        by["tau"]["tau_over_static_direct"] = (
            by["tau"]["wall_seconds"] / by["static_direct"]["wall_seconds"]
        )
    args.output.write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
