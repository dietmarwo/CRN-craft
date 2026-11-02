"""Step-one costs and prepared/one-shot equivalence on frozen W1–W4."""

import argparse
import concurrent.futures
import json
from pathlib import Path
import statistics
import subprocess
import tempfile
import time

import crn_craft as crn
from benchmark_campaign import models


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--profiler", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    report = {"schema": "crn-craft.ssa-profile/0.1", "workloads": {}}
    for name, model in models():
        m = (json.dumps(model) + "\n").encode()
        run = crn.run_spec(
            m,
            end=100.0,
            save_times=[0.0, 100.0],
            replicates=64,
            seed=4242,
            threads=1,
            max_events=10**9,
            max_output_bytes=200_000_000,
        )
        with tempfile.TemporaryDirectory(prefix="crn-profile-") as tmp:
            folder = Path(tmp)
            (folder / "model.json").write_bytes(m)
            (folder / "run.json").write_text(json.dumps(run))
            core = json.loads(
                subprocess.check_output(
                    [
                        str(args.profiler.resolve()),
                        str(folder / "model.json"),
                        str(folder / "run.json"),
                    ]
                )
            )
        prepared = crn.prepare(m, run, cached=True)
        reference = crn.simulate(m, run, cached=True)
        rows = {"one_shot": [], "prepared": []}
        for i in range(5):
            order = ["one_shot", "prepared"] if i % 2 == 0 else ["prepared", "one_shot"]
            for kind in order:
                t = time.perf_counter()
                value = (
                    prepared.simulate()
                    if kind == "prepared"
                    else crn.simulate(m, run, cached=True)
                )
                rows[kind].append(time.perf_counter() - t)
                if value != reference:
                    raise RuntimeError("prepared identity differs")
        # Concurrent reuse is checked after timing, on a small fixed request.
        small = crn.run_spec(m, end=1.0, save_times=[0.0, 1.0], replicates=4, keep=4)
        shared = crn.prepare(m, small, cached=True)
        expected = shared.simulate()
        with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
            values = list(
                pool.map(lambda _, shared=shared: shared.simulate(), range(3))
            )
        if any(v != expected for v in values):
            raise RuntimeError("concurrent prepared reuse differs")
        report["workloads"][name] = {
            "components": core,
            "python_seconds": rows,
            "medians": {k: statistics.median(v) for k, v in rows.items()},
            "prepared_replay": True,
            "concurrent_replay": True,
        }
        args.output.write_text(json.dumps(report, indent=2) + "\n")
        print(name, report["workloads"][name]["medians"], flush=True)
    report["passed"] = True
    args.output.write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
