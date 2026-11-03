"""Predeclared seven-workload comparison; validation precedes five timing rounds."""

import argparse
import hashlib
import gzip
import itertools
import json
import math
import os
from pathlib import Path
import platform
import random
import statistics
import sys
import subprocess
import time
import crn_craft as crn
import numpy as np
from scipy import stats
from acceleration_workloads import FOLDER, workloads, encode
from compare_julia import (
    Julia,
    ALGORITHMS,
    to_dsl,
    probe_native,
    cdf_error,
    digest,
    require,
)

ROOT = Path(__file__).resolve().parents[1]
VARIANTS = list(
    itertools.product(
        ["direct", "cached_direct", "direct_cr", "sorting_direct", "rssa_cr"],
        ["inverse", "ziggurat"],
    )
)
RUST = ["legacy_direct", "legacy_cached"] + [f"{a}/{w}" for a, w in VARIANTS]
ENGINES = RUST + [f"Julia/{a}" for a in ALGORITHMS]
TESTS, N, ALPHA = 6946, 512, 0.01


def rust_job(engine, model, end, paths, seed, keep, boundary):
    kwargs = dict(
        end=end,
        save_times=[0.0, end],
        replicates=paths,
        seed=seed,
        threads=1,
        keep=paths if keep else 0,
        max_events=10**9,
        max_output_bytes=200_000_000,
    )
    if engine.startswith("legacy_"):
        request = crn.run_spec(model, **kwargs)
        cached = engine == "legacy_cached"
        if boundary == "prepared":
            return crn.prepare(model, request, cached=cached).simulate
        return lambda: crn.simulate(model, request, cached=cached)
    method, wait = engine.split("/")
    request = crn.ssa_spec(
        model, aggregator=method, waiting_time=wait, fast=False, **kwargs
    )
    if boundary == "prepared":
        return crn.prepare_ssa(model, request).simulate
    return lambda: crn.simulate_ssa(model, request)


def payload(engine, result):
    if engine.startswith("Julia/"):
        return {k: v for k, v in result.items() if k != "seconds"}
    return {
        k: result[k]
        for k in ["events", "statistics", "trajectories", "replicates", "termination"]
    }


class Baseline:
    def __init__(self, package):
        env = os.environ.copy()
        env["PYTHONPATH"] = str(package)
        self.process = subprocess.Popen(
            [sys.executable, str(ROOT / "validation/gap10_python_reference.py")],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            text=True,
            env=env,
        )

    def call(self, **request):
        self.process.stdin.write(json.dumps(request) + "\n")
        self.process.stdin.flush()
        line = self.process.stdout.readline()
        require(bool(line), "baseline process stopped")
        return json.loads(line)

    def close(self):
        self.process.stdin.close()
        self.process.wait(timeout=30)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--julia", type=Path, required=True)
    p.add_argument("--binary", type=Path, required=True)
    p.add_argument("--depot", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--validation-only", action="store_true")
    p.add_argument("--validated", type=Path)
    p.add_argument("--baseline-package", type=Path, required=True)
    args = p.parse_args()
    project = ROOT / "validation/ssa-acceleration/julia"
    manifest = json.loads((FOLDER / "manifest.json").read_text())
    data = list(workloads())
    require(
        TESTS
        == sum(
            len(m["species"])
            * (len(ENGINES) if name in ["W1", "W2"] else len(ENGINES) - 1)
            for name, m, _, _ in data
        ),
        "comparison count changed",
    )
    source = (
        "\n\n".join(
            to_dsl(name, model).replace("initial=", f"horizon={end!r}, initial=")
            for name, model, end, _ in data
        )
        + "\n"
    )
    generated = project / "generated_models.jl"
    generated.write_text(source)
    models = {}
    for (name, model, _end, _paths), entry in zip(data, manifest["rows"], strict=True):
        m = encode(model)
        require(
            hashlib.sha256(m).hexdigest() == entry["sha256"]
            and m == (FOLDER / entry["path"]).read_bytes(),
            "frozen model changed",
        )
        models[name] = m
    report = dict(
        schema="crn-craft.julia-gap10-comparison/0.1",
        protocol="docs/decisions/julia-gap10-2026-11-02.md",
        workloads=manifest,
        engines=ENGINES,
        rounds=5,
        timing_seed=4242,
        validation_seed="777 + 1009 * engine index",
        validation={},
        probes={},
        timings={},
        first_runs={},
        preparation_seconds={},
        source_sha256={},
        environment=dict(
            platform=platform.platform(),
            python=sys.version,
            crn_craft=crn.__version__,
            numpy=np.__version__,
            cpu_model=next(
                line.split(":", 1)[1].strip()
                for line in Path("/proc/cpuinfo").read_text().splitlines()
                if line.startswith("model name")
            ),
            cpu_affinity=len(os.sched_getaffinity(0)),
            cpu_quota=Path("/sys/fs/cgroup/cpu.max").read_text().strip(),
        ),
    )
    sources = [
        Path(__file__).resolve(),
        ROOT / "validation/acceleration_workloads.py",
        ROOT / "validation/compare_julia.py",
        ROOT / "validation/compare_julia_accelerated.py",
        ROOT / "validation/gap10_python_reference.py",
        project / "reference_server.jl",
        project / "Project.toml",
        project / "Manifest.toml",
        generated,
    ]
    sources += sorted((ROOT / "crates/crn-craft/src").rglob("*.rs"))
    sources += [
        ROOT / "Cargo.lock",
        ROOT / "Cargo.toml",
        ROOT / "bindings/python/Cargo.toml",
        ROOT / "bindings/python/src/lib.rs",
        ROOT / "python/crn_craft/__init__.py",
    ]
    report["source_sha256"] = {
        str(v.relative_to(ROOT)): hashlib.sha256(v.read_bytes()).hexdigest()
        for v in sources
    }

    def save():
        args.output.write_text(json.dumps(report, indent=2) + "\n")

    baseline = Baseline(args.baseline_package)
    julia = Julia(
        args.julia.resolve(),
        project,
        args.depot,
        args.output.with_suffix(".stderr.log"),
    )
    try:
        start = time.perf_counter()
        report["julia_setup"] = julia.call(command="load", path=str(generated))
        report["julia_process_to_ready_seconds"] = time.perf_counter() - start
        save()
        print("Julia ready", flush=True)
        if args.validated:
            prior = json.loads(args.validated.read_text())
            require(
                prior["passed"]
                and prior["source_sha256"] == report["source_sha256"]
                and prior["workloads"] == manifest,
                "validation identity mismatch",
            )
            report["validation"] = prior["validation"]
            report["probes"] = prior["probes"]
            report["first_runs"] = prior["first_runs"]
            report["validation_evidence_sha256"] = hashlib.sha256(
                args.validated.read_bytes()
            ).hexdigest()
        else:
            rng = random.Random(1729)
            for name, model, end, _paths in data:
                ns = len(model["species"])
                states = [[v] * ns for v in [0, 1, 2, 10, 100]]
                if name != "W4":
                    states += [
                        [rng.randrange(101) for _ in range(ns)] for _ in range(3)
                    ]
                ours = probe_native(args.binary.resolve(), model, states)
                theirs = julia.call(command="probe", name=name, states=states)
                a, b = np.array(ours["rates"]), np.array(theirs["rates"])
                passed = ours["changes"] == theirs["changes"] and bool(
                    np.allclose(a, b, rtol=2e-12, atol=0)
                )
                report["probes"][name] = dict(
                    passed=passed,
                    comparisons=a.size,
                    max_relative_error=float(
                        np.max(np.abs(a - b) / np.maximum(a, 1e-300))
                    ),
                )
                save()
                require(passed, "propensity/stoichiometry probe failure")
                samples = {}
                report["first_runs"][name] = {}
                for engine in ENGINES:
                    # Distinct seeds: engines sharing a stream would make the
                    # two-sample gates compare correlated (or identical) samples.
                    seed = 777 + 1009 * ENGINES.index(engine)
                    if engine.startswith("Julia/"):
                        algorithm = engine.split("/")[1]
                        first = julia.call(
                            command="run",
                            name=name,
                            algorithm=algorithm,
                            paths=2,
                            seed=123,
                            keep=False,
                        )
                        report["first_runs"][name][engine] = first["seconds"]
                        result = julia.call(
                            command="run",
                            name=name,
                            algorithm=algorithm,
                            paths=N,
                            seed=seed,
                            keep=True,
                        )
                        final = np.array(result["final_counts"], dtype=np.int64)
                    else:
                        result = rust_job(
                            engine, models[name], end, N, seed, True, "prepared"
                        )()
                        require(
                            result["termination"] == "completed",
                            f"{name}/{engine} incomplete",
                        )
                        final = np.array(
                            [v["counts"][-1] for v in result["trajectories"]],
                            dtype=np.int64,
                        )
                    require(
                        final.shape == (N, ns) and np.all(final >= 0),
                        "invalid sample shape/counts",
                    )
                    if name == "W4":
                        require(np.all(final == final[:, :1]), "W4 invariant")
                    if name in ["W5", "W6"]:
                        require(
                            np.all(final.sum(axis=1) == 1024), "conversion conservation"
                        )
                    if name == "W7":
                        require(
                            np.all(
                                (final[:, 0::3] + final[:, 2::3]).sum(axis=1) == 32 * 28
                            ),
                            "A+C conservation",
                        )
                        require(
                            np.all(final[:, 1::3] + final[:, 2::3] == 16),
                            "B+C compartment conservation",
                        )
                    samples[engine] = final
                    print(name, engine, "validated samples", flush=True)
                rows = []
                reference = samples["legacy_cached"]
                for engine, values in samples.items():
                    if name not in ["W1", "W2"] and engine == "legacy_cached":
                        continue
                    for j in range(ns):
                        if name in ["W1", "W2"]:
                            mu = (
                                100 * (1 - math.exp(-100))
                                if name == "W1"
                                else 15 * (1 - math.exp(-20))
                            )
                            bound = math.sqrt(math.log(2 * TESTS / ALPHA) / (2 * N))
                            error = cdf_error(values[:, j], mu)
                        else:
                            bound = 2 * math.sqrt(math.log(4 * TESTS / ALPHA) / (2 * N))
                            error = float(
                                stats.ks_2samp(reference[:, j], values[:, j]).statistic
                            )
                        rows.append(
                            dict(
                                engine=engine,
                                species_index=j,
                                cdf_error=error,
                                cdf_bound=bound,
                                passed=error <= bound,
                            )
                        )
                sample_path = args.output.parent / f"{name}.samples.json.gz"
                sample_bytes = json.dumps(
                    {k: v.tolist() for k, v in samples.items()}, separators=(",", ":")
                ).encode()
                sample_path.write_bytes(gzip.compress(sample_bytes, mtime=0))
                report["validation"][name] = dict(
                    paths_per_engine=N,
                    rows=rows,
                    sample_sha256={k: digest(v.tolist()) for k, v in samples.items()},
                    sample_file=sample_path.name,
                    sample_file_sha256=hashlib.sha256(
                        sample_path.read_bytes()
                    ).hexdigest(),
                    passed=all(r["passed"] for r in rows),
                )
                save()
                require(
                    report["validation"][name]["passed"],
                    name + " declared distribution gate failed",
                )
            require(
                sum(len(v["rows"]) for v in report["validation"].values()) == TESTS,
                "validation count",
            )
        report["passed"] = all(v["passed"] for v in report["validation"].values())
        save()
        if args.validation_only:
            return
        for name, _model, end, paths in data:
            report.setdefault("baseline_setup", {})[name] = baseline.call(
                command="prepare",
                model_path=str(FOLDER / f"{name}.model.json"),
                end=end,
                paths=paths,
            )
            jobs = {}
            replays = {}
            report["timings"][name] = {}
            report["preparation_seconds"][name] = {}
            for engine in RUST:
                for boundary in ["one_shot", "prepared"]:
                    key = engine + "/" + boundary
                    start = time.perf_counter()
                    jobs[key] = rust_job(
                        engine, models[name], end, paths, 4242, False, boundary
                    )
                    report["preparation_seconds"][name][key] = (
                        time.perf_counter() - start
                    )
                    result = jobs[key]()
                    require(result["termination"] == "completed", key + " incomplete")
                    replays[key] = digest(payload(engine, result))
                    if boundary == "prepared":
                        require(
                            replays[key] == replays[engine + "/one_shot"],
                            "boundary replay differs",
                        )
                    report["timings"][name][key] = dict(
                        seconds=[],
                        events=result["events"],
                        work=result.get("work"),
                        numerical_sha256=replays[key],
                    )
            for engine in RUST:
                key = "Baseline/" + engine + "/prepared"
                result = baseline.call(command="run", engine=engine)["value"]
                require(result["termination"] == "completed", "baseline incomplete")
                replays[key] = digest(payload(engine, result))

                # Work counters may fall; numerical results and complete records must match.
                def without_work(value):
                    value = json.loads(json.dumps(value))
                    for rec in value.get("replicates", []):
                        rec.pop("work", None)
                    return value

                own = jobs[engine + "/prepared"]()
                require(
                    without_work(payload(engine, result))
                    == without_work(payload(engine, own)),
                    f"{name}/{engine} changed numerical paths",
                )
                report["timings"][name][key] = dict(
                    seconds=[], events=result["events"], numerical_sha256=replays[key]
                )
            for algorithm in ALGORITHMS:
                key = "Julia/" + algorithm
                result = julia.call(
                    command="run",
                    name=name,
                    algorithm=algorithm,
                    paths=paths,
                    seed=4242,
                    keep=False,
                )
                replays[key] = digest(payload(key, result))
                report["timings"][name][key] = dict(
                    seconds=[], numerical_sha256=replays[key]
                )
            keys = list(report["timings"][name])
            # Deterministic shuffled initial order, rotated/reversed each round.
            random.Random(4242).shuffle(keys)
            for round_index in range(5):
                offset = (round_index * 7) % len(keys)
                order = keys[offset:] + keys[:offset]
                if round_index % 2:
                    order = list(reversed(order))
                for key in order:
                    if key.startswith("Baseline/"):
                        engine = key[len("Baseline/") :].rsplit("/", 1)[0]
                        reply = baseline.call(command="run", engine=engine)
                        result = reply["value"]
                        elapsed = reply["seconds"]
                    elif key.startswith("Julia/"):
                        result = julia.call(
                            command="run",
                            name=name,
                            algorithm=key.split("/")[1],
                            paths=paths,
                            seed=4242,
                            keep=False,
                        )
                        elapsed = result["seconds"]
                        engine = key
                    else:
                        engine = key.rsplit("/", 1)[0]
                        start = time.perf_counter()
                        result = jobs[key]()
                        elapsed = time.perf_counter() - start
                        require(
                            result["termination"] == "completed", key + " incomplete"
                        )
                    require(
                        digest(payload(engine, result)) == replays[key],
                        key + " replay changed",
                    )
                    report["timings"][name][key]["seconds"].append(elapsed)
                save()
                print(name, "timing round", round_index + 1, flush=True)
            for row in report["timings"][name].values():
                row.update(
                    median_seconds=statistics.median(row["seconds"]),
                    min_seconds=min(row["seconds"]),
                    max_seconds=max(row["seconds"]),
                )
            save()
        report["complete"] = True
        save()
    except Exception as error:
        report["passed"] = False
        report["error"] = str(error)
        save()
        raise
    finally:
        julia.close()
        baseline.close()


if __name__ == "__main__":
    main()
