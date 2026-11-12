"""P10 sparse FSP allocation control, process cost and exact payload identity.

Since 0.7.0 (M8a) the model is converted with `crn-craft convert` and solved as
a `crn-craft.run/1.0` request with method `fsp` by `crn-craft run`; both
binaries must be 0.7.0 or later.
"""

import argparse
import hashlib
import json
from pathlib import Path
import statistics
import subprocess
import tempfile
import time

from benchmark_campaign import convert


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--before", type=Path, required=True)
    ap.add_argument("--after", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    rows = []
    with tempfile.TemporaryDirectory(prefix="crn-fsp-cost-") as temp:
        folder = Path(temp)
        model = {
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
        m = convert(args.before.resolve(), model, folder)
        (folder / "model.json").write_bytes(m)
        run = {
            "schema": "crn-craft.run/1.0",
            "model": {"path": "model.json", "sha256": hashlib.sha256(m).hexdigest()},
            "method": {
                "kind": "fsp",
                "bounds": {"X": 400},
                "max_states": 401,
                "tolerance": 1e-12,
            },
            "time": {"end": 40.0, "report_times": [0.0, 40.0]},
            "budget": {"max_output_bytes": 1000000},
        }
        (folder / "run.json").write_text(json.dumps(run))
        digests = []
        for _ in range(7):
            pair = {}
            for label, binary in [("before", args.before), ("after", args.after)]:
                start = time.perf_counter()
                proc = subprocess.run(
                    [
                        str(binary.resolve()),
                        "run",
                        str(folder / "model.json"),
                        str(folder / "run.json"),
                    ],
                    capture_output=True,
                    check=False,
                    timeout=120,
                )
                elapsed = time.perf_counter() - start
                if proc.returncode:
                    raise RuntimeError(proc.stdout.decode() + proc.stderr.decode())
                digests.append(hashlib.sha256(proc.stdout).hexdigest())
                pair[label] = elapsed
            rows.append(pair)
    report = {
        "schema": "crn-craft.fsp-performance/0.2",
        "identical": len(set(digests)) == 1,
        "numerical_sha256": digests[0],
        "rows": rows,
        "paired_median_speedup": statistics.median(
            r["before"] / r["after"] for r in rows
        ),
        "scope": "process timings, 401 reachable states, horizon 40, identical full FSP result including work/error fields",
    }
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    if not report["identical"]:
        raise RuntimeError("FSP payload changed")
    print(json.dumps({"identical": True, "speedup": report["paired_median_speedup"]}))


if __name__ == "__main__":
    main()
