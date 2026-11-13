#!/usr/bin/env python3
"""CRN-craft 0.7.0 against rebop's Python (dynamic) API, with compiled loops.

Declared in `validation/merge/RESULTS.md`, "rebop on 0.7.0, for the
comparison page: declaration". No speed target. Workloads W1 (`two_channel`)
and W2 (`sparse_66`) of `benchmark_campaign.models()`, the two that
`compare_references.py` measures: 64 paths over t in [0, 100], reports at 0
and 100, one thread, independent random streams in the two simulators.

Rebop has two routes: the Python API that builds a system dynamically, and a
Rust macro that compiles a system at build time. This script measures the
Python/dynamic route only.

Per workload the script

1. checks agreement: 512 paths from each simulator, two-sample
   Kolmogorov-Smirnov tests and mean z-scores on final counts (up to three
   species), and a pooled chi-square test of each simulator against the exact
   Poisson law at t = 100; a p-value below 0.01 / 8 is a failed test;
2. times, interleaved with alternating order, three repeats each: the
   CRN-craft default request through `crn.ensemble(...).simulate()`, the same
   request with a compiled Static Direct loop (`crn.compile_loop`, compile
   time reported separately), the same two with two threads, and 64 calls of
   `rebop.Gillespie.run`.

The report schema is `crn-craft.rebop-comparison/0.1`. `--only W1` measures
one workload; `--merge A.json B.json` merges parts by workload.

    python validation/merge/compare_rebop.py --output rebop.json [--only W1]
"""

import argparse
import importlib.metadata
import json
import math
import os
from pathlib import Path
import platform
import statistics
import subprocess
import sys
import time

import numpy as np
from scipy import stats

import crn_craft as crn

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "validation"))

from benchmark_campaign import models  # noqa: E402
from compare_references import reference  # noqa: E402

SCHEMA = "crn-craft.rebop-comparison/0.1"
NAMES = {"two_channel": "W1", "sparse_66": "W2"}
HORIZON, PATHS, SEED, CHECK_PATHS, REPEATS = 100.0, 64, 4242, 512, 3
CRN_CHECK_SEED, REBOP_CHECK_SEED = 778, 9000
TESTS = 8  # 4 Kolmogorov-Smirnov and 4 chi-square tests in all
THRESHOLD = 0.01 / TESTS


def ensemble(model_bytes, paths, seed, keep, threads=1, compiled=None):
    return crn.ensemble(
        model_bytes,
        end=HORIZON,
        report_times=[0.0, HORIZON],
        replicates=paths,
        seed=seed,
        threads=threads,
        keep_paths=keep,
        max_events=10**9,
        max_sampling_steps=10**10,
        max_output_bytes=200_000_000,
        method={"kind": "static_direct"},
        quantiles=(),
        compiled=compiled,
    )


def exact_poisson_mean(name):
    if name == "two_channel":
        return 100.0 * (1 - math.exp(-100.0))
    if name == "sparse_66":
        return 15.0 * (1 - math.exp(-20.0))
    return None


def rebop_run(reb, initial, seed):
    return reb.run(initial, HORIZON, 1, rng=seed)


def timed(function):
    start = time.perf_counter()
    function()
    return time.perf_counter() - start


def agreement(name, model, mb, reb, initial):
    names = [s["name"] for s in model["species"]]
    ours = ensemble(mb, CHECK_PATHS, CRN_CHECK_SEED, CHECK_PATHS).simulate()
    a = (
        ours.paths.sel(variable=names)
        .isel(time=-1)
        .to_numpy()
        .astype(float)
        .reshape(CHECK_PATHS, len(names))
    )
    b = np.array(
        [
            [
                float(rebop_run(reb, initial, REBOP_CHECK_SEED + i)[n].values[-1])
                for n in names
            ]
            for i in range(CHECK_PATHS)
        ]
    )
    report = {"paths_each": CHECK_PATHS, "two_sample": [], "tests": []}
    for j, species in enumerate(names[:3]):
        se = math.sqrt(a[:, j].var(ddof=1) / len(a) + b[:, j].var(ddof=1) / len(b))
        p = float(stats.ks_2samp(a[:, j], b[:, j]).pvalue)
        report["two_sample"].append(
            {
                "species": species,
                "crn_mean": float(a[:, j].mean()),
                "rebop_mean": float(b[:, j].mean()),
                "mean_z": float((a[:, j].mean() - b[:, j].mean()) / (se or 1.0)),
                "ks_pvalue": p,
            }
        )
        report["tests"].append({"test": f"ks {species}", "pvalue": p})
    mu = exact_poisson_mean(name)
    if mu is not None:
        law = stats.poisson(mu)
        for label, x in [("crn_craft", a), ("rebop", b)]:
            z = (x.mean(axis=0) - mu) / math.sqrt(mu / len(x))
            pooled = x.ravel().astype(int)
            support = np.arange(int(mu + 12 * math.sqrt(mu)))
            observed = np.bincount(pooled, minlength=len(support))[: len(support)]
            expected = law.pmf(support) * len(pooled)
            keep = expected > 5
            chi2 = ((observed[keep] - expected[keep]) ** 2 / expected[keep]).sum()
            p = float(stats.chi2.sf(chi2, keep.sum() - 1))
            report.setdefault("exact_poisson", {})[label] = {
                "mean": mu,
                "rms_mean_z_over_species": float(np.sqrt((z**2).mean())),
                "pooled_chi2_pvalue": p,
                "variance_over_mean": float(x.var(axis=0, ddof=1).mean() / mu),
            }
            report["tests"].append({"test": f"chi2 {label}", "pvalue": p})
    report["threshold"] = THRESHOLD
    for t in report["tests"]:
        t["failed"] = t["pvalue"] < THRESHOLD
    report["passed"] = not any(t["failed"] for t in report["tests"])
    return report


def host():
    cpu = ""
    try:
        for line in Path("/proc/cpuinfo").read_text().splitlines():
            if line.startswith("model name"):
                cpu = line.split(":", 1)[1].strip()
                break
    except OSError:
        pass
    return {
        "cpu": cpu,
        "cpu_count": os.cpu_count(),
        "kernel": platform.release(),
        "python": sys.version.split()[0],
    }


def measure(name, model, build_dir):
    label = NAMES[name]
    mb = crn.convert((json.dumps(model) + "\n").encode())
    initial = {s["name"]: s["initial_count"] for s in model["species"]}
    reb = reference(model)

    # Compiled loop: accepted or recorded as not accepted.
    compiled, compile_info = None, {}
    start = time.perf_counter()
    try:
        compiled = crn.compile_loop(mb, source=ROOT, cache_dir=build_dir)
        compile_info = {
            "accepted": True,
            "compile_seconds": time.perf_counter() - start,
            "cargo_build_seconds": compiled.build_seconds,
            "cached": compiled.cached,
        }
    except Exception as error:  # a refusal is a result
        compile_info = {"accepted": False, "error": f"{type(error).__name__}: {error}"}

    # Warm-ups, and the identity of the compiled and interpreted results.
    interpreted = ensemble(mb, PATHS, SEED, PATHS).simulate()
    identity = None
    if compiled:
        loop = ensemble(mb, PATHS, SEED, PATHS, compiled=compiled).simulate()
        identity = bool(
            np.array_equal(interpreted.paths.to_numpy(), loop.paths.to_numpy())
            and interpreted.metadata["events"] == loop.metadata["events"]
        )
    for i in range(PATHS):
        rebop_run(reb, initial, SEED + i)
    crn_events = int(interpreted.metadata["events"])
    rebop_events = sum(
        len(reb.run(initial, HORIZON, 0, rng=SEED + i)["time"]) - 1
        for i in range(PATHS)
    )

    check = agreement(name, model, mb, reb, initial)

    def run_rebop():
        for i in range(PATHS):
            rebop_run(reb, initial, SEED + i)

    rows = {
        "crn_default": lambda: ensemble(mb, PATHS, SEED, 0).simulate(),
        "crn_default_threads2": lambda: ensemble(mb, PATHS, SEED, 0, 2).simulate(),
    }
    if compiled:
        rows["crn_compiled"] = lambda: ensemble(
            mb, PATHS, SEED, 0, compiled=compiled
        ).simulate()
        rows["crn_compiled_threads2"] = lambda: ensemble(
            mb, PATHS, SEED, 0, 2, compiled=compiled
        ).simulate()
    seconds = {k: [] for k in [*rows, "rebop_dynamic"]}
    for repeat in range(REPEATS):
        order = [*rows, "rebop_dynamic"]
        if repeat % 2:
            order = ["rebop_dynamic", *rows]
        for key in order:
            seconds[key].append(
                timed(run_rebop if key == "rebop_dynamic" else rows[key])
            )
    median = {k: statistics.median(v) for k, v in seconds.items()}
    ratios = {
        f"rebop_over_{k}": median["rebop_dynamic"] / median[k]
        for k in ("crn_default", "crn_compiled")
        if k in median
    }
    return label, {
        "model": name,
        "species": len(model["species"]),
        "reactions": len(model["reactions"]),
        "seconds": seconds,
        "median_seconds": median,
        "ratios_one_thread": ratios,
        "crn_events": crn_events,
        "rebop_events": rebop_events,
        "ns_per_event": {
            **{k: 1e9 * median[k] / crn_events for k in median if k != "rebop_dynamic"},
            "rebop_dynamic": 1e9 * median["rebop_dynamic"] / rebop_events,
        },
        "compile_loop": compile_info,
        "compiled_identical_to_interpreter": identity,
        "agreement": check,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--only", choices=list(NAMES.values()))
    parser.add_argument("--build-dir", type=Path, help="cache of compiled loops")
    parser.add_argument("--merge", type=Path, nargs="+", help="merge part reports")
    args = parser.parse_args()
    if args.merge:
        parts = [json.loads(p.read_text()) for p in args.merge]
        merged = {k: v for k, v in parts[0].items() if k != "workloads"}
        merged["workloads"] = {}
        for part in parts:
            other = {k: v for k, v in part.items() if k != "workloads"}
            if other != merged and {k: v for k, v in other.items() if k != "part"} != {
                k: v for k, v in merged.items() if k != "part"
            }:
                raise SystemExit("parts differ in fields other than workloads")
            merged["workloads"].update(part["workloads"])
        merged.pop("part", None)
        args.output.write_text(json.dumps(merged, indent=2, sort_keys=True) + "\n")
        return
    report = {
        "schema": SCHEMA,
        "versions": {
            "crn_craft_python": crn.__version__,
            "rebop": importlib.metadata.version("rebop"),
            "numpy": np.__version__,
            "python": sys.version.split()[0],
            "rustc": subprocess.run(
                ["rustc", "--version"], capture_output=True, text=True
            ).stdout.strip(),
        },
        "host": host(),
        "work": f"{PATHS} paths, t in [0, {HORIZON}], output at 0 and {HORIZON}",
        "route": "rebop Python/dynamic API (rebop.Gillespie); not the Rust macro",
        "repeats": REPEATS,
        "agreement_threshold": THRESHOLD,
        "no_speed_target": True,
        "workloads": {},
    }
    if args.only:
        report["part"] = args.only
    build_dir = args.build_dir
    for name, model in list(models())[:2]:
        if args.only and NAMES[name] != args.only:
            continue
        label, entry = measure(name, model, build_dir)
        report["workloads"][label] = entry
        args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
        print(label, entry["median_seconds"], entry["ratios_one_thread"], flush=True)


if __name__ == "__main__":
    main()
