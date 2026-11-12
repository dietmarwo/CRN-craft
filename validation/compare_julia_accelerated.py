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
import time
import crn_craft as crn
import numpy as np
from scipy import stats
from acceleration_workloads import FOLDER, workloads, encode
from compare_julia import (
    Julia,
    ALGORITHMS,
    SELECTORS,
    crn_job,
    crn_payload,
    final_counts,
    to_dsl,
    probe_native,
    cdf_error,
    digest,
    require,
)

ROOT = Path(__file__).resolve().parents[1]
VARIANTS = list(itertools.product(SELECTORS, ["inverse", "ziggurat"]))
# linear/inverse is the reference engine of the two-sample gates (it replaces
# the retired cached Direct); cached Direct and `fast` no longer exist.
REFERENCE = "linear/inverse"
RUST = [f"{a}/{w}" for a, w in VARIANTS]
ENGINES = RUST + [f"Julia/{a}" for a in ALGORITHMS]
# 6946 before the retired legacy and cached engines were removed; same alpha.
TESTS, N, ALPHA = 5218, 512, 0.01


def rust_job(engine, model, end, paths, seed, keep, boundary):
    selector, wait = engine.split("/")
    job = crn_job(model, end, paths, seed, keep, selector, wait)
    if boundary == "prepared":
        return job.simulate
    return lambda: crn_job(model, end, paths, seed, keep, selector, wait).simulate()


def payload(engine, result):
    if engine.startswith("Julia/"):
        return {k: v for k, v in result.items() if k != "seconds"}
    return crn_payload(result)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--julia", type=Path, required=True)
    p.add_argument("--binary", type=Path, required=True)
    p.add_argument("--depot", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--validation-only", action="store_true")
    p.add_argument("--validated", type=Path)
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
        models[name] = model
    report = dict(
        schema="crn-craft.julia-comparison/0.2",
        protocol="docs/decisions/ssa-acceleration-2026-11-02.md",
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
        project / "reference_server.jl",
        project / "Project.toml",
        project / "Manifest.toml",
        generated,
    ]
    sources += sorted((ROOT / "crates/crn-craft/src").rglob("*.rs"))
    sources += [
        ROOT / "Cargo.lock",
        ROOT / "bindings/python/src/lib.rs",
        ROOT / "python/crn_craft/__init__.py",
    ]
    report["source_sha256"] = {
        str(v.relative_to(ROOT)): hashlib.sha256(v.read_bytes()).hexdigest()
        for v in sources
    }

    def save():
        args.output.write_text(json.dumps(report, indent=2) + "\n")

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
                            result.metadata["termination"] == "completed",
                            f"{name}/{engine} incomplete",
                        )
                        final = final_counts(result, models[name])
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
                reference = samples[REFERENCE]
                for engine, values in samples.items():
                    if name not in ["W1", "W2"] and engine == REFERENCE:
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
                    require(
                        result.metadata["termination"] == "completed",
                        key + " incomplete",
                    )
                    replays[key] = digest(payload(engine, result))
                    if boundary == "prepared":
                        require(
                            replays[key] == replays[engine + "/one_shot"],
                            "boundary replay differs",
                        )
                    report["timings"][name][key] = dict(
                        seconds=[],
                        events=int(result.metadata["events"]),
                        numerical_sha256=replays[key],
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
                    if key.startswith("Julia/"):
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
                            result.metadata["termination"] == "completed",
                            key + " incomplete",
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


if __name__ == "__main__":
    main()
