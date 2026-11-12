"""Measured Catalyst/JumpProcesses comparison on frozen W1-W4; see declared protocol.

Run with the installed current CRN-craft wheel and an instantiated validation/julia
environment. Numerical engines are unmodified. No Julia installation is automatic.
"""

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import random
import re
import statistics
import subprocess
import sys
import tempfile
import time

import crn_craft._legacy as crn
import numpy as np
from scipy import stats

from benchmark_campaign import models

ALGORITHMS = ["Direct", "SortingDirect", "DirectCR", "NRM", "RSSACR"]
NAMES = dict(zip([n for n, _ in models()], ["W1", "W2", "W3", "W4"], strict=True))
PATHS, CHECK_PATHS, SEED, ALPHA, TESTS = 64, 512, 4242, 0.01, 239


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def encode(model):
    return (json.dumps(model) + "\n").encode()


def value(model, v):
    return float(model["parameters"][v] if isinstance(v, str) else v)


def to_dsl(name, model):
    names = [s["name"] for s in model["species"]]
    require(
        all(re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*", s) for s in names), "unsafe name"
    )

    def side(part):
        return " + ".join(f"{n}{s}" if n != 1 else s for s, n in part.items()) or "0"

    lines = [f'MODELS["{name}"] = @reaction_network begin']
    for r in model["reactions"]:
        rate = r["rate"]
        kind = rate["kind"]
        arrow = "-->" if kind == "mass_action" else "=>"
        if kind == "mass_action":
            law = repr(value(model, rate["coefficient"]))
        elif kind == "hill_repression":
            v, b, k, n = (value(model, rate[x]) for x in ("vmax", "basal", "k", "n"))
            law = f"({b!r} + {v!r}/(1 + ({rate['regulator']}/{k!r})^{n!r}))"
            require(not r["reactants"], "unexpected Hill availability gate")
        elif kind == "michaelis_menten":
            require(name == "W4", "gate omission is authorized only for frozen W4")
            v, k = (value(model, rate[x]) for x in ("vmax", "km"))
            s = rate["substrate"]
            law = f"({v!r}*{s}/({k!r}+{s}))"
        else:
            raise ValueError(kind)
        lines.append(f"    {law}, {side(r['reactants'])} {arrow} {side(r['products'])}")
    lines.append("end")
    initial = [s["initial_count"] for s in model["species"]]
    lines.append(
        f'MODEL_INPUTS["{name}"] = (names={json.dumps(names)}, initial={initial})'
    )
    return "\n".join(lines)


class Julia:
    def __init__(self, binary, project, depot, error_log):
        env = os.environ.copy()
        env.update(JULIA_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
        if depot:
            env["JULIA_DEPOT_PATH"] = str(depot)
        self.log = error_log.open("w")
        self.process = subprocess.Popen(
            [
                str(binary),
                "--startup-file=no",
                "--threads=1",
                f"--project={project}",
                str(project / "reference_server.jl"),
            ],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=self.log,
            text=True,
            bufsize=1,
            env=env,
        )

    def call(self, **kwargs):
        self.process.stdin.write(json.dumps(kwargs) + "\n")
        self.process.stdin.flush()
        line = self.process.stdout.readline()
        require(bool(line), "Julia stopped; inspect the retained stderr log")
        reply = json.loads(line)
        require(reply["ok"], reply.get("error", "Julia error"))
        return reply["result"]

    def close(self):
        self.process.stdin.close()
        self.process.wait(timeout=30)
        self.log.close()


def run_crn(model, paths, seed, cached, keep=False):
    m = encode(model)
    request = crn.run_spec(
        m,
        end=100.0,
        save_times=[0.0, 100.0],
        replicates=paths,
        seed=seed,
        threads=1,
        keep=paths if keep else 0,
        max_events=10**9,
        max_output_bytes=200_000_000,
    )
    start = time.perf_counter()
    result = crn.simulate(m, request, cached=cached)
    seconds = time.perf_counter() - start
    require(result["termination"] == "completed", "incomplete CRN run")
    return result, seconds


def probe_native(binary, model, states):
    rows = []
    with tempfile.TemporaryDirectory(prefix="crn-julia-probe-") as temp:
        d = Path(temp)
        (d / "model.json").write_bytes(encode(model))
        names = [s["name"] for s in model["species"]]
        for state in states:
            (d / "state.json").write_text(
                json.dumps(
                    {
                        "schema": "crn-craft.state/0.1",
                        "counts": dict(zip(names, state, strict=True)),
                    }
                )
            )
            p = subprocess.run(
                [str(binary), "analyze", str(d / "model.json"), str(d / "state.json")],
                capture_output=True,
                check=True,
            )
            a = json.loads(p.stdout)["analysis"]
            rows.append(a["propensities"])
    changes = []
    for r in model["reactions"]:
        changes.append(
            [r["products"].get(s, 0) - r["reactants"].get(s, 0) for s in names]
        )
    return {"rates": rows, "changes": changes}


def cdf_error(samples, mu):
    support = np.arange(max(int(np.max(samples)) + 2, int(mu + 15 * math.sqrt(mu)) + 2))
    empirical = np.searchsorted(np.sort(samples), support, side="right") / len(samples)
    return float(np.max(np.abs(empirical - stats.poisson.cdf(support, mu))))


def digest(payload):
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--julia", type=Path, required=True)
    p.add_argument("--binary", type=Path, required=True)
    p.add_argument("--depot", type=Path)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--rounds", type=int, default=5)
    p.add_argument("--smoke", action="store_true")
    args = p.parse_args()
    root = Path(__file__).resolve().parents[1]
    project = root / "validation/julia"
    fixtures = {NAMES[n]: m for n, m in models()}
    dimer = json.loads(json.dumps(fixtures["W1"]))
    dimer["name"] = "DimerConvention"
    dimer["reactions"] = [
        {
            "name": "dimer",
            "reactants": {"X": 2},
            "products": {},
            "rate": {"kind": "mass_action", "coefficient": 0.25},
        }
    ]
    fixtures["dimer"] = dimer
    source = (
        "# Generated from the frozen CRN-craft workload definitions.\n"
        + "\n\n".join(to_dsl(n, m) for n, m in fixtures.items())
        + "\n"
    )
    generated = project / "generated_models.jl"
    generated.write_text(source)
    report = {
        "schema": "crn-craft.julia-comparison/0.1",
        "base_commit": "baed0e4590f93e4c1dcc77cd7f039a2604a1594c",
        "smoke": args.smoke,
        "protocol": "docs/decisions/julia-comparison-2026-11-02.md",
        "environment": {
            "platform": platform.platform(),
            "python": sys.version,
            "crn_craft": crn.__version__,
            "numpy": np.__version__,
            "cpu_model": next(
                line.split(":", 1)[1].strip()
                for line in Path("/proc/cpuinfo").read_text().splitlines()
                if line.startswith("model name")
            ),
            "cpu_affinity": len(os.sched_getaffinity(0)),
            "cpu_quota": Path("/sys/fs/cgroup/cpu.max").read_text().strip(),
        },
        "paths": 2 if args.smoke else PATHS,
        "timing_seed": SEED,
        "validation_seed": 777,
        "horizon": 100.0,
        "rounds": args.rounds,
        "source_sha256": {
            str(path.relative_to(root)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in [
                Path(__file__).resolve(),
                project / "reference_server.jl",
                project / "Project.toml",
                project / "Manifest.toml",
                generated,
                root / "validation/benchmark_campaign.py",
            ]
        },
        "model_sha256": {
            n: hashlib.sha256(encode(m)).hexdigest() for n, m in fixtures.items()
        },
        "probes": {},
        "validation": {},
        "timings": {},
        "first_runs": {},
    }

    def save():
        args.output.parent.mkdir(parents=True, exist_ok=True)
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
        print("Julia ready", report["julia_setup"], flush=True)
        rng = random.Random(1729)
        for name, model in fixtures.items():
            ns = len(model["species"])
            states = [[v] * ns for v in [0, 1, 2, 10, 20, 100, 1000]]
            if name != "W4":
                states += [[rng.randrange(101) for _ in range(ns)] for _ in range(16)]
            ours = probe_native(args.binary.resolve(), model, states)
            theirs = julia.call(command="probe", name=name, states=states)
            a, b = np.array(ours["rates"]), np.array(theirs["rates"])
            require(
                ours["changes"] == theirs["changes"], name + " stoichiometry mismatch"
            )
            require(
                np.allclose(a, b, rtol=2e-12, atol=0), name + " propensity mismatch"
            )
            report["probes"][name] = {
                "passed": True,
                "rate_comparisons": a.size,
                "maximum_relative_error": float(
                    np.max(np.abs(a - b) / np.maximum(a, 1e-300))
                ),
                "stoichiometry_identical": True,
            }
        save()
        check_paths = 16 if args.smoke else CHECK_PATHS
        for name, model in fixtures.items():
            if name == "dimer":
                continue
            ours, _ = run_crn(model, check_paths, 777, True, keep=True)
            a = np.array([t["counts"][-1] for t in ours["trajectories"]])
            all_samples = {"CRN-cached": a}
            report["first_runs"][name] = {}
            for algorithm in ALGORITHMS:
                first = julia.call(
                    command="run",
                    name=name,
                    algorithm=algorithm,
                    paths=2,
                    seed=123,
                    keep=False,
                )
                report["first_runs"][name][algorithm] = first["seconds"]
                result = julia.call(
                    command="run",
                    name=name,
                    algorithm=algorithm,
                    paths=check_paths,
                    seed=777,
                    keep=True,
                )
                b = np.array(result["final_counts"], dtype=np.int64)
                require(b.shape == a.shape, "sample shape")
                require(np.all(b >= 0), "negative counts")
                all_samples[algorithm] = b
                print(name, algorithm, "validation complete", flush=True)
            mu = (
                100 * (1 - math.exp(-100))
                if name == "W1"
                else 15 * (1 - math.exp(-20))
                if name == "W2"
                else None
            )
            rows = []
            if mu is not None:
                bound = math.sqrt(math.log(2 * TESTS / ALPHA) / (2 * check_paths))
                for engine, samples in all_samples.items():
                    for j in range(a.shape[1]):
                        error = cdf_error(samples[:, j], mu)
                        rows.append(
                            {
                                "engine": engine,
                                "species_index": j,
                                "cdf_error": error,
                                "cdf_bound": bound,
                                "mean_z": float(
                                    (samples[:, j].mean() - mu)
                                    / math.sqrt(mu / check_paths)
                                ),
                                "passed": error <= bound,
                            }
                        )
            else:
                bound = 2 * math.sqrt(math.log(4 * TESTS / ALPHA) / (2 * check_paths))
                for engine, samples in all_samples.items():
                    if name == "W4":
                        require(np.all(samples == samples[:, :1]), "W4 invariant")
                    if engine == "CRN-cached":
                        continue
                    for j in range(1 if name == "W4" else a.shape[1]):
                        error = float(stats.ks_2samp(a[:, j], samples[:, j]).statistic)
                        rows.append(
                            {
                                "engine": engine,
                                "species_index": j,
                                "cdf_error": error,
                                "cdf_bound": bound,
                                "passed": error <= bound,
                            }
                        )
            report["validation"][name] = {
                "paths_per_engine": check_paths,
                "rows": rows,
                "sample_sha256": {
                    k: digest(v.tolist()) for k, v in all_samples.items()
                },
                "final_samples": {k: v.tolist() for k, v in all_samples.items()},
                "passed": all(r["passed"] for r in rows),
            }
            save()
            require(
                report["validation"][name]["passed"], name + " distribution gate failed"
            )
        # All loading, JIT and independent validation precede these serial timings.
        engines = ["CRN-cached", "CRN-direct", *ALGORITHMS]
        for name, model in fixtures.items():
            if name == "dimer":
                continue
            report["timings"][name] = {e: [] for e in engines}
            paths = 2 if args.smoke else PATHS
            for rep in range(args.rounds):
                order = engines[rep % len(engines) :] + engines[: rep % len(engines)]
                for engine in order:
                    if engine.startswith("CRN-"):
                        result, elapsed = run_crn(
                            model, paths, SEED, engine == "CRN-cached"
                        )
                        row = {
                            "seconds": elapsed,
                            "events": result["events"],
                            "digest": digest(result),
                        }
                    else:
                        result = julia.call(
                            command="run",
                            name=name,
                            algorithm=engine,
                            paths=paths,
                            seed=SEED,
                            keep=False,
                        )
                        elapsed = result.pop("seconds")
                        row = {"seconds": elapsed, "digest": digest(result)}
                    report["timings"][name][engine].append(row)
                save()
                print(name, "timing round", rep + 1, flush=True)
            for engine, rows in report["timings"][name].items():
                require(
                    len({r["digest"] for r in rows}) == 1,
                    name + " " + engine + " replay differs",
                )
            require(
                report["timings"][name]["CRN-cached"][0]["digest"]
                == report["timings"][name]["CRN-direct"][0]["digest"],
                "native identity differs",
            )
        report["summary"] = {
            name: {
                engine: {
                    "median_seconds": statistics.median(r["seconds"] for r in rows),
                    "min_seconds": min(r["seconds"] for r in rows),
                    "max_seconds": max(r["seconds"] for r in rows),
                }
                for engine, rows in engines.items()
            }
            for name, engines in report["timings"].items()
        }
        report["passed"] = True
        save()
        print(json.dumps(report["summary"], indent=2), flush=True)
    except Exception as exc:
        report["passed"] = False
        report["failure"] = str(exc)
        save()
        raise
    finally:
        julia.close()


if __name__ == "__main__":
    main()
