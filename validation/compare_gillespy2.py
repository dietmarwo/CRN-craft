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
3. times, interleaved: the CRN-craft CLI (`crn-craft run`, process wall
   time), the CRN-craft Python API (`crn.ensemble(...).simulate()` for the
   default Static Direct configuration, for the checked configuration
   (linear selector, inverse waits) and for the `linear`,
   `composition_rejection`, `sorting` and `rssa_cr` selectors with Ziggurat
   waits; call wall time including the preparation), GillesPy2's compiled
   SSACSolver (compiled once, compile time reported separately) and its
   NumPySSASolver.

Since 0.7.0 (M8a) the script runs the unified engine: each network model is
converted with `crn_craft.convert` and run as a `crn-craft.run/1.0` Static
Direct request. The 0.6.0 rows `simulate-cached` and the `fast` policy
(cached Direct, xoshiro256++, no loss or resolution checks) no longer exist,
and neither does the check of the cached Direct method against GillesPy2; the
selectors with Ziggurat waits are the like-for-like setting now, and the
default configuration is validated against GillesPy2 as well. `--crn-only`
leaves GillesPy2 out (validation and timings), which is the smoke form.

W4: every reaction changes all eight species by the same amount, so they
stay equal and CRN-craft's reactant gates on the Michaelis-Menten deaths are
redundant. The timed GillesPy2 model omits them (the validation checks that
the species stay equal); the gated C++ model is timed as well and reported.

Requires the installed `crn_craft` wheel and validation/reference-requirements.txt.

    python validation/compare_gillespy2.py --binary target/release/crn-craft \
        --output gillespy2.json --figure gillespy2.png
"""

import argparse
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

import crn_craft as crn
from benchmark_campaign import models

NAMES = {"two_channel": "W1", "sparse_66": "W2", "circuit_12": "W3", "dense_66": "W4"}
HORIZON, PATHS, SEED, CHECK_PATHS = 100.0, 64, 4242, 512
SELECTORS = ["linear", "composition_rejection", "sorting", "rssa_cr"]


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


def run_spec(model_bytes, paths, seed, keep, selector=None, waiting_time=None):
    """A prepared Static Direct run of the converted model (`model_bytes`)."""
    method = {"kind": "static_direct"}
    if waiting_time is not None:
        method["waiting_time"] = waiting_time
    return crn.ensemble(
        model_bytes,
        end=HORIZON,
        report_times=[0.0, HORIZON],
        replicates=paths,
        seed=seed,
        keep_paths=keep,
        max_events=10**9,
        max_sampling_steps=10**10,
        max_output_bytes=200_000_000,
        method=method,
        selector=selector,
        quantiles=(),
    )


def checked(model_bytes, paths, seed, keep):
    """The checked configuration of 0.6.0: linear scan, inverse waits."""
    return run_spec(model_bytes, paths, seed, keep, "linear", "inverse")


def final_counts(result, names):
    """Final counts of every replicate, one column per species, in `names` order."""
    return result.paths.sel(variable=names).isel(time=-1).to_numpy().astype(float)


def exact_poisson(name):
    """Mean of the Poisson law of every species at t = 100, if known."""
    if name == "two_channel":
        return 100.0 * (1 - math.exp(-100.0))
    if name == "sparse_66":
        return 15.0 * (1 - math.exp(-20.0))
    return None


def validate(name, model, gmodel, solver, ziggurat=False):
    mb = crn.convert(encoded(model))
    names = [s["name"] for s in model["species"]]
    if ziggurat:
        ours = run_spec(mb, CHECK_PATHS, 778, CHECK_PATHS).simulate()
    else:
        ours = checked(mb, CHECK_PATHS, 777, CHECK_PATHS).simulate()
    theirs = gmodel.run(
        solver=solver, number_of_trajectories=CHECK_PATHS, seed=779 if ziggurat else 777
    )
    a = final_counts(ours, names)
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


def cli_once(binary, model_bytes, folder, request):
    (folder / "model.json").write_bytes(model_bytes)
    (folder / "run.json").write_bytes(crn.rules.model_bytes(request))
    start = time.perf_counter()
    proc = subprocess.run(
        [str(binary), "run", str(folder / "model.json"), str(folder / "run.json")],
        capture_output=True,
        check=False,
    )
    wall = time.perf_counter() - start
    if proc.returncode:
        raise RuntimeError(proc.stdout.decode())
    return wall, json.loads(proc.stdout)["events"]


def timed(function, *args, **kwargs):
    start = time.perf_counter()
    function(*args, **kwargs)
    return time.perf_counter() - start


def summary(values):
    return {"median": statistics.median(values), "runs": values}


def plot(report, path):
    """Time per event of the compiled solvers, from a report."""
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
        ("python_default", "CRN-craft Static Direct (default)", "#2a78d6"),
        ("python_checked", "CRN-craft linear, inverse waits", "#8fb5e6"),
        ("python_composition_rejection", "CRN-craft composition-rejection", "#eb6834"),
        ("gillespy2_ssa_cpp", "GillesPy2 SSACSolver (C++)", "#1baf7a"),
    ]
    names = [
        n for n in report["workloads"] if "gillespy2_ssa_cpp" in report["workloads"][n]
    ]
    fig, axis = plt.subplots(figsize=(9, 3.8), facecolor="#fcfcfb")
    axis.set_facecolor("#fcfcfb")
    width = 0.2

    def events(n, key):
        entry = report["workloads"][n]
        if key == "python_checked" or key == "gillespy2_ssa_cpp":
            return entry["crn_events"]
        return entry["event_counts"][key.removeprefix("python_")]

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
    parser.add_argument(
        "--crn-only",
        action="store_true",
        help="smoke form: time and run only the CRN-craft side, no GillesPy2",
    )
    parser.add_argument(
        "--only", choices=list(NAMES.values()), help="measure one workload only"
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
        "schema": "crn-craft.gillespy2-comparison/0.3",
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
        "crn_only": args.crn_only,
        "workloads": {},
    }
    with tempfile.TemporaryDirectory(prefix="crn-gillespy2-") as temp:
        folder = Path(temp)
        for name, model in models():
            label = NAMES[name]
            if args.only and label != args.only:
                continue
            redundant_gates = name == "dense_66"
            mb = crn.convert(encoded(model))
            variants = {
                "default": lambda mb=mb: run_spec(mb, PATHS, SEED, 0),
                "checked": lambda mb=mb: checked(mb, PATHS, SEED, 0),
                **{
                    s: lambda mb=mb, s=s: run_spec(mb, PATHS, SEED, 0, s, "ziggurat")
                    for s in SELECTORS
                },
            }
            warm = variants["default"]().simulate()  # warm-up
            events = {
                k: int(make().simulate().metadata["events"])
                for k, make in variants.items()
            }
            rows = {
                k: []
                for k in [
                    "cli_run_process",
                    *[f"python_{k}" for k in variants],
                    *([] if args.crn_only else ["gillespy2_ssa_cpp"]),
                ]
            }
            request = json.loads(json.dumps(warm.request))
            request["model"]["path"] = "model.json"
            if not args.crn_only:
                gmodel = to_gillespy2(model, gates=not redundant_gates)
                start = time.perf_counter()
                csolver = gp.SSACSolver(model=gmodel)
                compile_seconds = time.perf_counter() - start
                nsolver = gp.NumPySSASolver(model=gmodel)
                check = validate(name, model, gmodel, csolver)
                ziggurat_check = validate(name, model, gmodel, csolver, ziggurat=True)
                gmodel.run(solver=csolver, number_of_trajectories=PATHS, seed=SEED)
            cli_events = None
            for repeat in range(args.repeats):
                wall, cli_events = cli_once(binary, mb, folder, request)
                rows["cli_run_process"].append(wall)
                for key, make in variants.items():
                    rows[f"python_{key}"].append(
                        timed(lambda make=make: make().simulate())
                    )
                if not args.crn_only:
                    rows["gillespy2_ssa_cpp"].append(
                        timed(
                            gmodel.run,
                            solver=csolver,
                            number_of_trajectories=PATHS,
                            seed=SEED + repeat,
                        )
                    )
            entry = {k: summary(v) for k, v in rows.items()}
            if not args.crn_only:
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
                entry["validation"] = check
                entry["ziggurat_validation"] = ziggurat_check
            entry["crn_events"] = events["checked"]
            entry["cli_events"] = cli_events
            entry["event_counts"] = events
            report["workloads"][label] = entry
            medians = {
                k: round(v["median"], 4)
                for k, v in entry.items()
                if isinstance(v, dict) and "median" in v
            }
            print(label, events["checked"], medians, flush=True)
            args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    if args.figure and not args.crn_only:
        plot(report, args.figure)


if __name__ == "__main__":
    main()
