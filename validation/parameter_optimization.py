"""Complete fixed-topology optimization comparison, including AOT compilation."""

import argparse
import hashlib
import importlib.metadata
import json
import math
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import time
import numpy as np
import rebop
import fcmaes_rust
from parameter_optimization_models import ROOT, FOLDER, exact_moments, fixtures, loss

ENGINES = [
    "crn_fresh_fast",
    "crn_reuse_fast",
    "crn_generated_fast",
    "rebop_expression",
    "rebop_rebuild_lma",
    "rebop_macro",
]
WORKER = FOLDER / "worker/Cargo.toml"


def source_hashes():
    paths = list((ROOT / "crates/crn-craft/src").rglob("*.rs"))
    paths += list((FOLDER / "worker/src").glob("*.rs"))
    paths += [
        ROOT / "Cargo.lock",
        ROOT / "Cargo.toml",
        WORKER,
        WORKER.with_name("Cargo.lock"),
        FOLDER / "fixtures.json",
        ROOT / "validation/parameter_optimization.py",
        ROOT / "validation/parameter_optimization_models.py",
        ROOT / "validation/generate_parameter_workers.py",
        ROOT / "validation/check_parameter_interface.py",
        ROOT / "python/crn_craft/__init__.py",
        ROOT / "bindings/python/src/lib.rs",
    ]
    return {
        str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in sorted(paths)
    }


def save(path, payload):
    path.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n")


def build(target):
    env = os.environ.copy()
    env["CARGO_TARGET_DIR"] = str(target)
    env["CARGO_BUILD_JOBS"] = "4"
    common = [
        "cargo",
        "build",
        "--release",
        "--locked",
        "--offline",
        "--manifest-path",
        str(WORKER),
    ]
    # Prime model-independent dependencies and retain first/failed build logs.
    with (FOLDER / "evidence/build-prime.log").open("w") as log:
        subprocess.run(
            common + ["--features", "generated,rebop-macro"],
            env=env,
            cwd=ROOT,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=True,
        )
    result = dict(
        note="Dependency-cache-warm package rebuilds; whole three-model bundle plus arithmetic guard. No clean-machine installation claim.",
        workers={},
    )
    binaries = target / "parameter-worker-binaries"
    binaries.mkdir(exist_ok=True)
    for name, feature in [
        ("dynamic", None),
        ("generated", "generated"),
        ("rebop_macro", "rebop-macro"),
    ]:
        subprocess.run(
            [
                "cargo",
                "clean",
                "--release",
                "--manifest-path",
                str(WORKER),
                "-p",
                "crn-parameter-bench-worker",
            ],
            env=env,
            cwd=ROOT,
            capture_output=True,
            check=True,
        )
        start = time.perf_counter()
        with (FOLDER / f"evidence/build-{name}.log").open("w") as log:
            subprocess.run(
                common + (["--features", feature] if feature else []),
                env=env,
                cwd=ROOT,
                stdout=log,
                stderr=subprocess.STDOUT,
                check=True,
            )
        seconds = time.perf_counter() - start
        destination = binaries / name
        shutil.copyfile(target / "release/crn-parameter-bench-worker", destination)
        destination.chmod(0o755)
        result["workers"][name] = dict(
            binary=str(destination),
            compilation_seconds=seconds,
            sha256=hashlib.sha256(destination.read_bytes()).hexdigest(),
            feature=feature,
        )
        save(FOLDER / "evidence/builds.json", result)
        print("build", name, seconds, flush=True)
    return result


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
            raise RuntimeError("Worker stopped: " + self.process.stderr.read())
        reply = json.loads(line)
        if not reply["ok"]:
            raise RuntimeError(reply["error"])
        return reply["result"]

    def close(self):
        self.process.stdin.close()
        self.process.wait(timeout=30)
        errors = self.process.stderr.read()
        if self.process.returncode or errors:
            raise RuntimeError("Worker failed: " + errors)
        self.process.stdout.close()
        self.process.stderr.close()


def path_seeds(seed, paths):
    return [
        int(
            np.random.default_rng(seed + i).integers(
                np.iinfo(np.uint64).max, dtype=np.uint64
            )
        )
        for i in range(paths)
    ]


def reference_model(model, values, expression):
    names = sorted(model["parameters"])
    parameters = dict(zip(names, values, strict=True))
    simulation = rebop.Gillespie()
    for r in model["reactions"]:
        a = [s for s, n in r["reactants"].items() for _ in range(n)]
        b = [s for s, n in r["products"].items() for _ in range(n)]
        coefficient = r["rate"]["coefficient"]
        divisor = math.prod(math.factorial(n) for n in r["reactants"].values())
        if expression:
            terms = [coefficient if isinstance(coefficient, str) else repr(coefficient)]
            for s, n in r["reactants"].items():
                terms += [s if i == 0 else f"({s}-{i})" for i in range(n)]
            rate = "*".join(terms) + f"/{divisor}"
        else:
            rate = (
                parameters[coefficient] if isinstance(coefficient, str) else coefficient
            ) / divisor
        simulation.add_reaction(rate, a, b)
    return simulation


class PythonRebop:
    def __init__(self, fixture, row, expression):
        self.fixture, self.row, self.expression = fixture, row, expression
        self.initial = {
            s["name"]: s["initial_count"] for s in fixture["model"]["species"]
        }
        self.simulation = (
            reference_model(fixture["model"], row["true_values"], True)
            if expression
            else None
        )

    def evaluate(self, values, phase):
        row = self.row
        paths = row["training_paths"] if phase == "train" else row["heldout_paths"]
        seed = row["training_seed"] if phase == "train" else row["heldout_seed"]
        simulation = (
            self.simulation
            if self.expression
            else reference_model(self.fixture["model"], values, False)
        )
        parameters = dict(zip(row["parameter_names"], values, strict=True))
        states = []
        for i in range(paths):
            result = simulation.run(
                self.initial,
                row["times"][-1],
                len(row["times"]) - 1,
                params=parameters if self.expression else {},
                rng=seed + i,
                sparse=self.fixture["case"] == "sparse33",
            )
            if not np.array_equal(result.time.values, row["times"]):
                raise RuntimeError("rebop save grid mismatch")
            states.append(np.stack([result[s].values for s in self.initial], axis=1))
        states = np.stack(states)
        if not np.all(states >= 0):
            raise RuntimeError("Negative rebop population")
        return dict(
            mean=states.mean(axis=0).tolist(),
            variance=states.var(axis=0, ddof=1).tolist(),
        )


def start_worker(binary, engine, fixture, row):
    worker = Worker(binary)
    worker.call(
        command="init",
        engine=engine,
        case=fixture["case"],
        model=fixture["model"],
        times=row["times"],
        paths=row["training_paths"],
        heldout_paths=row["heldout_paths"],
        seed=row["training_seed"],
    )
    return worker


def one_run(engine, fixture, row, optimizer_seed, builds, generation):
    start_total = time.perf_counter()
    if engine.startswith("crn_") or engine == "rebop_macro":
        variant = (
            "generated"
            if engine == "crn_generated_fast"
            else "rebop_macro"
            if engine == "rebop_macro"
            else "dynamic"
        )
        mode = dict(
            crn_fresh_fast="fresh",
            crn_reuse_fast="reuse",
            crn_generated_fast="generated",
            rebop_macro="rebop_macro",
        )[engine]
        worker = start_worker(builds["workers"][variant]["binary"], mode, fixture, row)
        seeds = {
            phase: path_seeds(row[key], row[count])
            for phase, key, count in [
                ("train", "training_seed", "training_paths"),
                ("heldout", "heldout_seed", "heldout_paths"),
            ]
        }

        def evaluate(values, phase):
            reply = worker.call(
                command="eval",
                values=list(map(float, values)),
                phase=phase,
                path_seeds=seeds[phase],
            )
            return reply["moments"], {
                k: reply[k] for k in ["bind_seconds", "simulation_seconds"]
            }
    else:
        variant, worker = None, None
        reference = PythonRebop(fixture, row, engine == "rebop_expression")

        def evaluate(values, phase):
            start = time.perf_counter()
            moments = reference.evaluate(values, phase)
            return moments, dict(
                bind_seconds=None, simulation_seconds=time.perf_counter() - start
            )

    startup = time.perf_counter() - start_total
    trace = []

    def objective(x):
        begin = time.perf_counter()
        values = np.exp(x)
        actual, components = evaluate(values, "train")
        value = loss(actual, row["target"])
        trace.append(
            dict(
                x=list(map(float, x)),
                values=values.tolist(),
                loss=value,
                wall_seconds=time.perf_counter() - begin,
                **components,
            )
        )
        return value

    start = time.perf_counter()
    try:
        optimum = fcmaes_rust.optimize_bite(
            objective,
            np.log(row["guess"]).tolist(),
            np.log(row["lower"]).tolist(),
            np.log(row["upper"]).tolist(),
            seed=optimizer_seed,
            max_evaluations=128,
        )
        optimization = time.perf_counter() - start
        values = np.exp(optimum[0])
        begin = time.perf_counter()
        heldout, _ = evaluate(values, "heldout")
        heldout_seconds = time.perf_counter() - begin
        exact = exact_moments(fixture, values)
        if worker:
            worker.close()
            worker = None
        elapsed = time.perf_counter() - start_total
        compilation = (
            builds["workers"][variant]["compilation_seconds"] if variant else 0.0
        )
        generation_seconds = (
            sum(
                v
                for k, v in generation.items()
                if k.endswith("/crn" if variant == "generated" else "/rebop")
            )
            if variant in ["generated", "rebop_macro"]
            else 0.0
        )
        return dict(
            engine=engine,
            case=fixture["case"],
            optimizer_seed=optimizer_seed,
            passed=True,
            requested_evaluations=128,
            actual_evaluations=len(trace),
            x=list(map(float, optimum[0])),
            values=values.tolist(),
            training_loss=float(optimum[1]),
            heldout_loss=loss(heldout, row["target"]),
            exact_moment_loss=loss(exact, row["target"]),
            relative_parameter_errors=np.abs(
                values / row["true_values"] - 1.0
            ).tolist(),
            compilation_seconds=compilation,
            generation_seconds=generation_seconds,
            startup_seconds=startup,
            optimization_seconds=optimization,
            objective_seconds=sum(t["wall_seconds"] for t in trace),
            heldout_seconds=heldout_seconds,
            measured_run_seconds=elapsed,
            total_with_compilation_seconds=elapsed + compilation + generation_seconds,
            trace=trace,
        )
    finally:
        if worker:
            worker.close()


def run():
    validation = json.loads((FOLDER / "evidence/parameter-validation.json").read_text())
    gates = json.loads((FOLDER / "evidence/final/summary.json").read_text())
    if (
        not validation["passed"]
        or not gates.get("complete")
        or not gates["passed"]
        or validation["source_sha256"] != source_hashes()
    ):
        raise RuntimeError("Final gates/source validation required before timing")
    builds = json.loads((FOLDER / "evidence/builds.json").read_text())
    for row in builds["workers"].values():
        if (
            hashlib.sha256(Path(row["binary"]).read_bytes()).hexdigest()
            != row["sha256"]
        ):
            raise RuntimeError("Worker binary changed")
    generation = json.loads((FOLDER / "evidence/generation.json").read_text())
    manifest = json.loads((FOLDER / "fixtures.json").read_text())
    report = dict(
        schema="crn-craft.parameter-optimization-comparison/0.1",
        complete=False,
        source_sha256=source_hashes(),
        builds=builds,
        versions={
            name: importlib.metadata.version(name)
            for name in ["crn-craft", "rebop", "fcmaes-rust", "numpy", "scipy"]
        },
        environment=dict(
            platform=platform.platform(),
            python=sys.version,
            cpu=next(
                line.split(":", 1)[1].strip()
                for line in Path("/proc/cpuinfo").read_text().splitlines()
                if line.startswith("model name")
            ),
            affinity=len(os.sched_getaffinity(0)),
            # cgroup v2 only; absent on cgroup v1 hosts and in a v2 root.
            quota=(
                Path("/sys/fs/cgroup/cpu.max").read_text().strip()
                if Path("/sys/fs/cgroup/cpu.max").exists()
                else None
            ),
        ),
        rows=[],
    )
    output = FOLDER / "evidence/optimization.json"
    fixtures_by_name = {f["case"]: f for f in fixtures()}
    for case_index, row in enumerate(manifest["rows"]):
        for seed_index, seed in enumerate(manifest["optimizer_seeds"]):
            offset = (case_index + seed_index) % len(ENGINES)
            order = ENGINES[offset:] + ENGINES[:offset]
            if seed_index % 2:
                order = order[::-1]
            for engine in order:
                print("start", row["case"], seed, engine, flush=True)
                try:
                    result = one_run(
                        engine,
                        fixtures_by_name[row["case"]],
                        row,
                        seed,
                        builds,
                        generation,
                    )
                except Exception as exc:
                    result = dict(
                        engine=engine,
                        case=row["case"],
                        optimizer_seed=seed,
                        passed=False,
                        error=repr(exc),
                    )
                report["rows"].append(result)
                save(output, report)
                print(
                    "done",
                    engine,
                    result.get("optimization_seconds"),
                    result.get("exact_moment_loss"),
                    result.get("error"),
                    flush=True,
                )
    # Whole optimizer traces must agree across the three CRN routes, excluding timings.
    identity = []
    for row in manifest["rows"]:
        for seed in manifest["optimizer_seeds"]:
            choices = [
                r
                for r in report["rows"]
                if r["case"] == row["case"]
                and r["optimizer_seed"] == seed
                and r["engine"].startswith("crn_")
            ]
            traces = [
                [(v["x"], v["loss"]) for v in r.get("trace", [])] for r in choices
            ]
            equal = (
                len(choices) == 3
                and all(r["passed"] for r in choices)
                and all(t == traces[0] for t in traces)
            )
            identity.append(dict(case=row["case"], seed=seed, passed=equal))
    report["trace_identities"] = identity
    report["complete"] = True
    report["passed"] = all(r["passed"] for r in report["rows"]) and all(
        r["passed"] for r in identity
    )
    save(output, report)
    if not report["passed"]:
        raise RuntimeError("Failed optimization rows retained")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--build", action="store_true")
    p.add_argument("--run", action="store_true")
    p.add_argument("--target", type=Path, default=ROOT / "target/parameter-bench")
    args = p.parse_args()
    if args.build:
        build(args.target.resolve())
    if args.run:
        run()


if __name__ == "__main__":
    main()
