#!/usr/bin/env python3
"""CRN-craft against GillesPy2 on the frozen SSA workloads W1-W4.

Each workload is 64 paths over t in [0, 100] with output at 0 and 100, on one
thread, with independent random streams in the two simulators. The script

1. converts each CRN-craft model to a GillesPy2 model (mass action as
   GillesPy2 mass action; Hill and Michaelis-Menten laws as propensity
   expressions; CRN-craft's reactant gates as x / (x + 1e-300), which is
   exactly 0 at x = 0 and 1 otherwise);
2. checks that both simulate the same process: 512 paths each, two-sample
   Kolmogorov-Smirnov tests and mean z-scores on final counts, and for W1 and
   W2, whose laws are Poisson at t = 100, a comparison of each simulator
   with the exact law;
3. times, interleaved: the CRN-craft CLI (`simulate`, `simulate-cached`;
   process wall and core time), the CRN-craft Python API (`crn.simulate`,
   plain and cached; and `crn.simulate_ssa` with `fast=True` for the
   cached Direct, DirectCR, SortingDirect and RSSACR methods with Ziggurat
   waits; call wall time), GillesPy2's compiled SSACSolver (compiled once,
   compile time reported separately) and its NumPySSASolver.

`fast=True` is CRN-craft's JumpProcesses-style policy (Xoshiro256++,
running totals without loss or resolution checks, batched stop polls). Like
GillesPy2, it does not check every positive channel against the 53-bit
selection resolution, so it is the like-for-like setting; the checked policy
(`fast` off) is timed as well. The script passes `fast` explicitly, so its
requests mean the same before and after `fast` became the default
(ssa-run/0.3). The fast cached Direct method is the closest
algorithmic match to GillesPy2's SSACSolver (Direct method); it is also
validated against GillesPy2 below.

W4: every reaction changes all eight species by the same amount, so they
stay equal and CRN-craft's reactant gates on the Michaelis-Menten deaths are
redundant. The timed GillesPy2 model omits them (the validation checks that
the species stay equal); the gated C++ model is timed as well and reported.

Requires the installed `crn_craft` wheel and validation/reference-requirements.txt.

    python validation/compare_gillespy2.py --binary target/release/crn-craft \
        --output gillespy2.json --figure gillespy2.png
"""

import argparse
import hashlib
import json
import math
from pathlib import Path
import statistics
import subprocess
import sys
import tempfile
import time

import gillespy2 as gp
import numpy as np
from scipy import stats

import crn_craft._legacy as crn
from benchmark_campaign import models

NAMES = {"two_channel": "W1", "sparse_66": "W2", "circuit_12": "W3", "dense_66": "W4"}
HORIZON, PATHS, SEED, CHECK_PATHS = 100.0, 64, 4242, 512
FAST = ["cached_direct", "direct_cr", "sorting_direct", "rssa_cr"]


def value(model, x):
    return model["parameters"][x] if isinstance(x, str) else x


def to_gillespy2(model, gates=True):
    g = gp.Model(name=model["name"], population=True, volume=1.0)
    species = {
        s["name"]: gp.Species(
            name=s["name"], initial_value=s["initial_count"], mode="discrete"
        )
        for s in model["species"]
    }
    g.add_species(list(species.values()))
    for i, r in enumerate(model["reactions"]):
        if any(n != 1 for n in r["reactants"].values()):
            raise ValueError("W1-W4 have unit reactant coefficients only")
        rate = r["rate"]
        reactants = {species[s]: n for s, n in r["reactants"].items()}
        products = {species[s]: n for s, n in r["products"].items()}
        if rate["kind"] == "mass_action":
            k = gp.Parameter(name=f"k{i}", expression=value(model, rate["coefficient"]))
            g.add_parameter(k)
            g.add_reaction(
                gp.Reaction(
                    name=r["name"], reactants=reactants, products=products, rate=k
                )
            )
            continue
        if rate["kind"] in ("hill_repression", "hill_activation"):
            x = rate["regulator"]
            vmax, basal, k, n = (
                value(model, rate[f]) for f in ("vmax", "basal", "k", "n")
            )
            ratio = f"pow({x}/{k!r}, {n!r})"
            law = (
                f"({basal!r} + {vmax!r}/(1 + {ratio}))"
                if rate["kind"] == "hill_repression"
                else f"({basal!r} + {vmax!r}*{ratio}/(1 + {ratio}))"
            )
        elif rate["kind"] == "michaelis_menten":
            s = rate["substrate"]
            vmax, km = value(model, rate["vmax"]), value(model, rate["km"])
            law = f"({vmax!r}*{s}/({km!r} + {s}))"
        else:
            raise ValueError(rate["kind"])
        terms = [law]
        if gates:
            terms += [f"({s}/({s} + 1e-300))" for s in r["reactants"]]
        g.add_reaction(
            gp.Reaction(
                name=r["name"],
                reactants=reactants,
                products=products,
                propensity_function="*".join(terms),
            )
        )
    g.timespan(np.array([0.0, HORIZON]))
    return g


def encoded(model):
    return (json.dumps(model) + "\n").encode()


def run_spec(model, paths, seed, keep):
    return crn.run_spec(
        encoded(model),
        end=HORIZON,
        save_times=[0.0, HORIZON],
        replicates=paths,
        seed=seed,
        keep=keep,
        max_events=10**9,
        max_output_bytes=200_000_000,
    )


def fast_spec(model, paths, seed, keep, aggregator):
    return crn.ssa_spec(
        encoded(model),
        aggregator=aggregator,
        waiting_time="ziggurat",
        fast=True,
        end=HORIZON,
        save_times=[0.0, HORIZON],
        replicates=paths,
        seed=seed,
        keep=keep,
        max_events=10**9,
        max_sampling_steps=10**10,
        max_output_bytes=200_000_000,
    )


def exact_poisson(name):
    """Mean of the Poisson law of every species at t = 100, if known."""
    if name == "two_channel":
        return 100.0 * (1 - math.exp(-100.0))
    if name == "sparse_66":
        return 15.0 * (1 - math.exp(-20.0))
    return None


def validate(name, model, gmodel, solver, fast=False):
    if fast:
        spec = fast_spec(model, CHECK_PATHS, 778, CHECK_PATHS, "cached_direct")
        ours = crn.simulate_ssa(encoded(model), spec)
    else:
        spec = run_spec(model, CHECK_PATHS, 777, CHECK_PATHS)
        ours = crn.simulate(encoded(model), spec, cached=True)
    theirs = gmodel.run(
        solver=solver, number_of_trajectories=CHECK_PATHS, seed=779 if fast else 777
    )
    names = [s["name"] for s in model["species"]]
    a = np.array([p["counts"][-1] for p in ours["trajectories"]], dtype=float)
    b = np.array([[t[n][-1] for n in names] for t in theirs], dtype=float)
    report = {"paths_each": CHECK_PATHS, "two_sample": []}
    for j, species in enumerate(names[:3]):
        se = math.sqrt(a[:, j].var(ddof=1) / len(a) + b[:, j].var(ddof=1) / len(b))
        report["two_sample"].append(
            {
                "species": species,
                "crn_mean": float(a[:, j].mean()),
                "gillespy2_mean": float(b[:, j].mean()),
                "mean_z": float((a[:, j].mean() - b[:, j].mean()) / (se or 1.0)),
                "ks_pvalue": float(stats.ks_2samp(a[:, j], b[:, j]).pvalue),
            }
        )
    mu = exact_poisson(name)
    if mu is not None:
        law = stats.poisson(mu)
        for label, x in [("crn_craft", a), ("gillespy2", b)]:
            z = (x.mean(axis=0) - mu) / math.sqrt(mu / len(x))
            pooled = x.ravel().astype(int)
            support = np.arange(int(mu + 12 * math.sqrt(mu)))
            observed = np.bincount(pooled, minlength=len(support))[: len(support)]
            expected = law.pmf(support) * len(pooled)
            keep = expected > 5
            chi2 = ((observed[keep] - expected[keep]) ** 2 / expected[keep]).sum()
            report.setdefault("exact_poisson", {})[label] = {
                "mean": mu,
                "rms_mean_z_over_species": float(np.sqrt((z**2).mean())),
                "pooled_chi2_pvalue": float(stats.chi2.sf(chi2, keep.sum() - 1)),
                "variance_over_mean": float(x.var(axis=0, ddof=1).mean() / mu),
            }
    if name == "dense_66":
        report["species_stay_equal"] = {
            "crn_craft": bool((a == a[:, :1]).all()),
            "gillespy2": bool((b == b[:, :1]).all()),
        }
    return report


def cli_once(binary, command, model, folder, label):
    m = encoded(model)
    (folder / "model.json").write_bytes(m)
    run = run_spec(model, PATHS, SEED, 0)
    run["model"] = {"path": "model.json", "sha256": hashlib.sha256(m).hexdigest()}
    (folder / "run.json").write_text(json.dumps(run))
    out = folder / label
    start = time.perf_counter()
    proc = subprocess.run(
        [str(binary), command, str(folder / "model.json"), str(folder / "run.json")]
        + [str(out)],
        capture_output=True,
        check=False,
    )
    wall = time.perf_counter() - start
    if proc.returncode:
        raise RuntimeError(proc.stdout.decode())
    core = json.loads((out / "provenance.json").read_text())["elapsed_seconds"]
    events = json.loads((out / "result.json").read_text())["events"]
    return wall, core, events


def timed(function, *args, **kwargs):
    start = time.perf_counter()
    function(*args, **kwargs)
    return time.perf_counter() - start


def summary(values):
    return {"median": statistics.median(values), "runs": values}


def plot(report, path):
    """Time per event of the three compiled solvers, from a report."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    labels = {
        "W1": "W1 two channels",
        "W2": "W2 sparse 66",
        "W3": "W3 repressilator",
        "W4": "W4 dense 66",
    }
    series = [
        ("python_fast_cached_direct", "CRN-craft fast cached Direct", "#2a78d6"),
        ("python_cached", "CRN-craft simulate-cached (checked)", "#8fb5e6"),
        ("python_simulate", "CRN-craft simulate", "#eb6834"),
        ("gillespy2_ssa_cpp", "GillesPy2 SSACSolver (C++)", "#1baf7a"),
    ]
    names = list(report["workloads"])
    fig, axis = plt.subplots(figsize=(9, 3.8), facecolor="#fcfcfb")
    axis.set_facecolor("#fcfcfb")
    width = 0.2

    def events(n, key):
        entry = report["workloads"][n]
        if key.startswith("python_fast_"):
            return entry["fast_events"][key.removeprefix("python_fast_")]
        return entry["crn_events"]

    for i, (key, label, color) in enumerate(series):
        values = [
            1e9 * report["workloads"][n][key]["median"] / events(n, key) for n in names
        ]
        x = [j + (i - 1.5) * width for j in range(len(names))]
        bars = axis.bar(x, values, width * 0.92, color=color, label=label)
        axis.bar_label(bars, fmt="%.0f", fontsize=8, color="#52514e", padding=2)
    axis.set_xticks(range(len(names)), [labels[n] for n in names])
    axis.set_ylabel("ns per event (median of 5)")
    axis.set_title(
        "Time per SSA event, 64 paths, one thread (lower is faster)",
        fontsize=11,
        fontweight="bold",
    )
    for side in ("top", "right"):
        axis.spines[side].set_visible(False)
    axis.grid(axis="y", color="#e4e3df", linewidth=0.8)
    axis.set_axisbelow(True)
    axis.legend(frameon=False, loc="upper left")
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--numpy-repeats", type=int, default=3)
    parser.add_argument("--figure", type=Path, help="also plot time per event")
    parser.add_argument(
        "--plot-only",
        action="store_true",
        help="plot an existing --output report instead of measuring",
    )
    args = parser.parse_args()
    if args.plot_only:
        plot(json.loads(args.output.read_text()), args.figure)
        return
    if args.binary is None:
        parser.error("--binary is required unless --plot-only is given")
    binary = args.binary.resolve()
    version = subprocess.run([str(binary), "--version"], capture_output=True)
    report = {
        "schema": "crn-craft.gillespy2-comparison/0.2",
        "versions": {
            "crn_craft_cli": json.loads(version.stdout)["version"],
            "crn_craft_python": crn.__version__,
            "gillespy2": gp.__version__,
            "python": sys.version.split()[0],
            "cxx": subprocess.run(
                ["g++", "--version"], capture_output=True, text=True
            ).stdout.splitlines()[0],
        },
        "work": f"{PATHS} paths, t in [0, {HORIZON}], output at 0 and {HORIZON}, one thread",
        "workloads": {},
    }
    with tempfile.TemporaryDirectory(prefix="crn-gillespy2-") as temp:
        folder = Path(temp)
        for name, model in models():
            label = NAMES[name]
            redundant_gates = name == "dense_66"
            gmodel = to_gillespy2(model, gates=not redundant_gates)
            start = time.perf_counter()
            csolver = gp.SSACSolver(model=gmodel)
            compile_seconds = time.perf_counter() - start
            nsolver = gp.NumPySSASolver(model=gmodel)
            check = validate(name, model, gmodel, csolver)
            fast_check = validate(name, model, gmodel, csolver, fast=True)
            mb, plain = encoded(model), run_spec(model, PATHS, SEED, 0)
            crn.simulate(mb, plain)  # warm-up
            fast = {m: fast_spec(model, PATHS, SEED, 0, m) for m in FAST}
            fast_events = {m: crn.simulate_ssa(mb, fast[m])["events"] for m in FAST}
            gmodel.run(solver=csolver, number_of_trajectories=PATHS, seed=SEED)
            rows = {
                k: []
                for k in [
                    "cli_simulate_process",
                    "cli_simulate_core",
                    "cli_cached_process",
                    "cli_cached_core",
                    "python_simulate",
                    "python_cached",
                    *[f"python_fast_{m}" for m in FAST],
                    "gillespy2_ssa_cpp",
                ]
            }
            events = None
            for repeat in range(args.repeats):
                for command, key in [
                    ("simulate", "cli_simulate"),
                    ("simulate-cached", "cli_cached"),
                ]:
                    wall, core, events = cli_once(
                        binary, command, model, folder, f"{label}-{key}-{repeat}"
                    )
                    rows[key + "_process"].append(wall)
                    rows[key + "_core"].append(core)
                rows["python_simulate"].append(timed(crn.simulate, mb, plain))
                rows["python_cached"].append(
                    timed(crn.simulate, mb, plain, cached=True)
                )
                for m in FAST:
                    rows[f"python_fast_{m}"].append(
                        timed(crn.simulate_ssa, mb, fast[m])
                    )
                rows["gillespy2_ssa_cpp"].append(
                    timed(
                        gmodel.run,
                        solver=csolver,
                        number_of_trajectories=PATHS,
                        seed=SEED + repeat,
                    )
                )
            entry = {k: summary(v) for k, v in rows.items()}
            entry["gillespy2_numpy"] = summary(
                [
                    timed(
                        gmodel.run,
                        solver=nsolver,
                        number_of_trajectories=PATHS,
                        seed=SEED + repeat,
                    )
                    for repeat in range(args.numpy_repeats)
                ]
            )
            if redundant_gates:
                gated = to_gillespy2(model, gates=True)
                gsolver = gp.SSACSolver(model=gated)
                entry["gillespy2_ssa_cpp_with_gates"] = summary(
                    [
                        timed(
                            gated.run,
                            solver=gsolver,
                            number_of_trajectories=PATHS,
                            seed=SEED + repeat,
                        )
                        for repeat in range(3)
                    ]
                )
            entry["gillespy2_ssa_cpp_compile_seconds"] = compile_seconds
            entry["gillespy2_gates_omitted_as_redundant"] = redundant_gates
            entry["crn_events"] = events
            entry["fast_events"] = fast_events
            entry["validation"] = check
            entry["fast_validation"] = fast_check
            report["workloads"][label] = entry
            medians = {
                k: round(v["median"], 4)
                for k, v in entry.items()
                if isinstance(v, dict) and "median" in v
            }
            print(label, events, medians, flush=True)
            args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    if args.figure:
        plot(report, args.figure)


if __name__ == "__main__":
    main()
