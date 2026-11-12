"""Time the parameter benchmark worker before and after M5g (declared in RESULTS.md).

Each route runs in its own worker process: `init` with a frozen case, 128
training evaluations at log-uniform points (the same for every route), one
held-out evaluation at the true values. The report holds the medians over
repetitions of the replies' bind and simulation seconds and of the wall time
of the 129 calls, and asserts that the routes of one engine return identical
moments.
"""

import argparse
import json
import platform
from pathlib import Path
import statistics
import subprocess
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
FOLDER = ROOT / "validation/parameter-optimization"
EVALUATIONS = 128


class Worker:
    def __init__(self, binary):
        self.process = subprocess.Popen(
            [str(binary)],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            bufsize=1,
        )

    def call(self, **request):
        self.process.stdin.write(json.dumps(request, allow_nan=False) + "\n")
        self.process.stdin.flush()
        line = self.process.stdout.readline()
        if not line:
            raise RuntimeError("worker stopped: " + self.process.stderr.read())
        reply = json.loads(line)
        if not reply["ok"]:
            raise RuntimeError(reply["error"])
        return reply["result"]

    def close(self):
        self.process.stdin.close()
        self.process.wait(timeout=60)
        errors = self.process.stderr.read()
        self.process.stdout.close()
        self.process.stderr.close()
        if self.process.returncode or errors:
            raise RuntimeError("worker failed: " + errors)


def points(row, index):
    rng = np.random.default_rng(20261112 + index)
    low, high = np.log(row["lower"]), np.log(row["upper"])
    u = rng.random((EVALUATIONS, len(low)))
    return np.exp(low + u * (high - low)).tolist()


def run_route(binary, engine, model, row, values):
    worker = Worker(binary)
    try:
        worker.call(
            command="init",
            engine=engine,
            case=row["case"],
            model=model,
            times=row["times"],
            paths=row["training_paths"],
            heldout_paths=row["heldout_paths"],
            seed=row["training_seed"],
        )
        bind = simulation = 0.0
        moments = []
        start = time.perf_counter()
        for phase, value in [("train", v) for v in values] + [
            ("heldout", row["true_values"])
        ]:
            reply = worker.call(command="eval", values=value, phase=phase)
            bind += reply["bind_seconds"]
            simulation += reply["simulation_seconds"]
            moments.append(reply["moments"])
        wall = time.perf_counter() - start
    finally:
        worker.close()
    return dict(bind=bind, simulation=simulation, wall=wall), moments


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--old-dynamic", type=Path, required=True)
    parser.add_argument("--old-generated", type=Path, required=True)
    parser.add_argument("--new", type=Path, required=True)
    parser.add_argument("--repetitions", type=int, default=3)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    routes = [
        ("old", "fresh", args.old_dynamic),
        ("old", "reuse", args.old_dynamic),
        ("old", "generated", args.old_generated),
        ("new", "fresh", args.new),
        ("new", "reuse", args.new),
    ]
    manifest = json.loads((FOLDER / "fixtures.json").read_text())
    report = dict(
        schema="crn-craft.parameter-worker-speed/0.1",
        declared_in="validation/merge/RESULTS.md (M5g declaration)",
        python=sys.version.split()[0],
        platform=platform.platform(),
        evaluations=EVALUATIONS + 1,
        repetitions=args.repetitions,
        cases=[],
    )
    for index, row in enumerate(manifest["rows"]):
        model = json.loads(
            (FOLDER / "models" / f"{row['case']}.model.json").read_text()
        )
        values = points(row, index)
        samples = {f"{e}/{r}": [] for e, r, _ in routes}
        moments = {}
        for repetition in range(args.repetitions):
            shift = (index + repetition) % len(routes)
            for engine, route, binary in routes[shift:] + routes[:shift]:
                timing, seen = run_route(binary, route, model, row, values)
                samples[f"{engine}/{route}"].append(timing)
                key = f"{engine}/{route}"
                if key in moments and moments[key] != seen:
                    raise RuntimeError(f"{row['case']} {key}: moments changed")
                moments[key] = seen
                print(row["case"], key, repetition, timing, flush=True)
        for engine in ("old", "new"):
            keys = [k for k in moments if k.startswith(engine + "/")]
            if any(moments[k] != moments[keys[0]] for k in keys):
                raise RuntimeError(f"{row['case']}: {engine} routes differ")
        report["cases"].append(
            dict(
                case=row["case"],
                routes={
                    key: {
                        field: statistics.median(t[field] for t in timings)
                        for field in ("bind", "simulation", "wall")
                    }
                    | {"wall_samples": [t["wall"] for t in timings]}
                    for key, timings in samples.items()
                },
                identical_moments_within_engine=True,
            )
        )
        args.output.write_text(json.dumps(report, indent=2) + "\n")
    for case in report["cases"]:
        r = case["routes"]
        case["ratios"] = {
            "new_fresh/old_fresh": r["new/fresh"]["wall"] / r["old/fresh"]["wall"],
            "new_reuse/old_reuse": r["new/reuse"]["wall"] / r["old/reuse"]["wall"],
            "new_reuse/old_generated": r["new/reuse"]["wall"]
            / r["old/generated"]["wall"],
        }
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps([(c["case"], c["ratios"]) for c in report["cases"]], indent=1))


if __name__ == "__main__":
    main()
