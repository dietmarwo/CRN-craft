#!/usr/bin/env python3
"""CRN-craft against JumpProcesses.jl (through Catalyst.jl) on W1-W4.

The frozen SSA workloads are written as Catalyst reaction networks
(mass action with `-->`; Hill and Michaelis-Menten laws as complete
propensities with `=>`, which is CRN-craft's custom-rate semantics), and
validation/julia/bench.jl times 64 trajectories over t in [0, 100], saving at
0 and 100, serially, with several JumpProcesses aggregators. The first solve
of each aggregator, which includes Julia's compilation, is reported apart.

CRN-craft is timed in the same session through its Python API (plain and
cached), matching docs/comparison-gillespy2.md. Before the numbers are used,
512 Direct trajectories from Julia are compared with 512 CRN-craft paths
(two-sample Kolmogorov-Smirnov and mean z on final counts) and, for W1 and
W2, with the exact Poisson laws.

    python validation/compare_julia.py --julia julia --output julia.json \
        --figure julia.png
"""

import argparse
import json
import math
from pathlib import Path
import statistics
import subprocess
import tempfile
import time

import numpy as np
from scipy import stats

import crn_craft as crn
from benchmark_campaign import models

ROOT = Path(__file__).resolve().parent
NAMES = {"two_channel": "W1", "sparse_66": "W2", "circuit_12": "W3", "dense_66": "W4"}
HORIZON, PATHS, SEED, CHECK_PATHS = 100.0, 64, 4242, 512


def value(model, x):
    return float(model["parameters"][x] if isinstance(x, str) else x)


def side(terms):
    if not terms:
        return "0"
    return " + ".join(s if n == 1 else f"{n}{s}" for s, n in terms.items())


def catalyst_source(label, model):
    """A Catalyst @reaction_network with numeric rates."""
    lines = []
    for r in model["reactions"]:
        rate = r["rate"]
        left, right = side(r["reactants"]), side(r["products"])
        if rate["kind"] == "mass_action":
            # Catalyst's mass action divides by stoichiometric factorials,
            # as CRN-craft's binomial convention does.
            lines.append(f"{value(model, rate['coefficient'])!r}, {left} --> {right}")
            continue
        if rate["kind"] in ("hill_repression", "hill_activation"):
            function = "hillr" if rate["kind"] == "hill_repression" else "hill"
            law = (
                f"{value(model, rate['basal'])!r} + {function}({rate['regulator']}, "
                f"{value(model, rate['vmax'])!r}, {value(model, rate['k'])!r}, "
                f"{value(model, rate['n'])!r})"
            )
        elif rate["kind"] == "michaelis_menten":
            law = (
                f"mm({rate['substrate']}, {value(model, rate['vmax'])!r}, "
                f"{value(model, rate['km'])!r})"
            )
        else:
            raise ValueError(rate["kind"])
        # `=>` uses the rate as the complete propensity, without a mass-action
        # factor. CRN-craft also gates custom laws on available reactants;
        # in W3 (no reactants) and W4 (species stay equal, the law vanishes at
        # zero substrate) the gate never changes a propensity.
        lines.append(f"({law}), {left} => {right}")
    body = "\n        ".join(lines)
    u0 = ", ".join(f":{s['name']} => {s['initial_count']}" for s in model["species"])
    names = ", ".join(f'"{s["name"]}"' for s in model["species"])
    return (
        f'    ("{label}", @reaction_network {label} begin\n        {body}\n    end,\n'
        f"     [{u0}], [{names}]),\n"
    )


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


def timed(function, *args, **kwargs):
    start = time.perf_counter()
    function(*args, **kwargs)
    return time.perf_counter() - start


def exact_poisson(name):
    if name == "two_channel":
        return 100.0 * (1 - math.exp(-100.0))
    if name == "sparse_66":
        return 15.0 * (1 - math.exp(-20.0))
    return None


def validate(name, model, julia_counts):
    names = [s["name"] for s in model["species"]]
    ours = crn.simulate(
        encoded(model), run_spec(model, CHECK_PATHS, 777, CHECK_PATHS), cached=True
    )
    a = np.array([p["counts"][-1] for p in ours["trajectories"]], dtype=float)
    b = np.array([julia_counts[n] for n in names], dtype=float).T
    report = {"paths_each": CHECK_PATHS, "two_sample": []}
    for j, species in enumerate(names[:3]):
        se = math.sqrt(a[:, j].var(ddof=1) / len(a) + b[:, j].var(ddof=1) / len(b))
        report["two_sample"].append(
            {
                "species": species,
                "crn_mean": float(a[:, j].mean()),
                "julia_mean": float(b[:, j].mean()),
                "mean_z": float((a[:, j].mean() - b[:, j].mean()) / (se or 1.0)),
                "ks_pvalue": float(stats.ks_2samp(a[:, j], b[:, j]).pvalue),
            }
        )
    mu = exact_poisson(name)
    if mu is not None:
        law = stats.poisson(mu)
        for label, x in [("crn_craft", a), ("julia", b)]:
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
            "julia": bool((b == b[:, :1]).all()),
        }
    return report


def plot(report, path):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    labels = {
        "W1": "W1 two channels",
        "W2": "W2 sparse 66",
        "W3": "W3 repressilator",
        "W4": "W4 dense 66",
    }
    names = list(report["workloads"])
    best = {
        n: min(
            report["workloads"][n]["julia"].items(),
            key=lambda item: statistics.median(item[1]["runs"]),
        )
        for n in names
    }
    series = [
        ("CRN-craft simulate-cached", "#2a78d6", lambda n: "python_cached"),
        ("CRN-craft simulate", "#eb6834", lambda n: "python_simulate"),
    ]
    fig, axis = plt.subplots(figsize=(9, 3.8), facecolor="#fcfcfb")
    axis.set_facecolor("#fcfcfb")
    width = 0.26

    def per_event(n, seconds):
        return 1e9 * seconds / report["workloads"][n]["crn_events"]

    columns = [
        (
            label,
            color,
            [per_event(n, report["workloads"][n][key(n)]["median"]) for n in names],
        )
        for label, color, key in series
    ]
    columns.append(
        (
            "JumpProcesses.jl, best aggregator",
            "#1baf7a",
            [per_event(n, statistics.median(best[n][1]["runs"])) for n in names],
        )
    )
    for i, (label, color, values) in enumerate(columns):
        x = [j + (i - 1) * width for j in range(len(names))]
        bars = axis.bar(x, values, width * 0.92, color=color, label=label)
        axis.bar_label(bars, fmt="%.0f", fontsize=8, color="#52514e", padding=2)
    axis.set_xticks(
        range(len(names)), [f"{labels[n]}\n(Julia: {best[n][0]})" for n in names]
    )
    axis.set_ylabel("ns per event (median of 5)")
    axis.set_title(
        "Time per SSA event, 64 paths, one thread (lower is faster)",
        fontsize=11,
        fontweight="bold",
    )
    for spine in ("top", "right"):
        axis.spines[spine].set_visible(False)
    axis.grid(axis="y", color="#e4e3df", linewidth=0.8)
    axis.set_axisbelow(True)
    axis.legend(frameon=False, loc="upper left")
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--julia", default="julia")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--figure", type=Path)
    parser.add_argument("--plot-only", action="store_true")
    args = parser.parse_args()
    if args.plot_only:
        plot(json.loads(args.output.read_text()), args.figure)
        return
    workloads = list(models())
    with tempfile.TemporaryDirectory(prefix="crn-julia-") as temp:
        folder = Path(temp)
        source = "const WORKLOADS = [\n" + "".join(
            catalyst_source(NAMES[name], model) for name, model in workloads
        )
        (folder / "models.jl").write_text(source + "]\n")
        start = time.perf_counter()
        subprocess.run(
            [
                args.julia,
                f"--project={ROOT / 'julia'}",
                str(ROOT / "julia/bench.jl"),
                str(folder / "models.jl"),
                str(folder / "julia.json"),
            ],
            env={**__import__("os").environ, "REPEATS": str(args.repeats)},
            check=True,
        )
        julia_wall = time.perf_counter() - start
        julia = json.loads((folder / "julia.json").read_text())
        julia_models = source
    version = subprocess.run(
        [args.julia, "--version"], capture_output=True, text=True
    ).stdout.strip()
    report = {
        "schema": "crn-craft.julia-comparison/0.1",
        "versions": {
            "crn_craft_python": crn.__version__,
            "julia": julia["julia"],
            "julia_cli": version,
            **julia["packages"],
        },
        "work": f"{PATHS} paths, t in [0, {HORIZON}], output at 0 and {HORIZON}, one thread",
        "julia_process_seconds": julia_wall,
        "catalyst_models": julia_models,
        "workloads": {},
    }
    for name, model in workloads:
        label = NAMES[name]
        mb, plain = encoded(model), run_spec(model, PATHS, SEED, 0)
        crn.simulate(mb, plain)  # warm-up
        events = crn.simulate(mb, plain)["events"]
        rows = {"python_simulate": [], "python_cached": []}
        for _ in range(args.repeats):
            rows["python_simulate"].append(timed(crn.simulate, mb, plain))
            rows["python_cached"].append(timed(crn.simulate, mb, plain, cached=True))
        entry = {
            k: {"median": statistics.median(v), "runs": v} for k, v in rows.items()
        }
        jw = julia["workloads"][label]
        entry["julia"] = {
            k: v for k, v in jw.items() if k != "final_counts" and "runs" in v
        }
        entry["crn_events"] = events
        entry["validation"] = validate(name, model, jw["final_counts"])
        report["workloads"][label] = entry
        print(
            label,
            events,
            {k: round(v["median"], 4) for k, v in entry.items() if "median" in v},
            {
                k: round(statistics.median(v["runs"]), 4)
                for k, v in entry["julia"].items()
            },
            flush=True,
        )
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    if args.figure:
        plot(report, args.figure)


if __name__ == "__main__":
    main()
