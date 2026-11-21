"""Measure neutral heteroplasmy variance under mitochondrial fusion and fission."""

import argparse
import json
from pathlib import Path

import crn_craft as crn
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
M = json.loads((ROOT / "manifest.json").read_text())
NAMES = ("WS", "MS", "WF", "MF")
# One row per channel in Equations 1–9, with 0 for empty output.
NETWORK = (
    ("fusion_WW", ["WS", "WS"], ["WF", "WF"]),
    ("fusion_MM", ["MS", "MS"], ["MF", "MF"]),
    ("fusion_WF_WS", ["WF", "WS"], ["WF", "WF"]),
    ("fusion_MF_MS", ["MF", "MS"], ["MF", "MF"]),
    ("fusion_WF_MS", ["WF", "MS"], ["WF", "MF"]),
    ("fusion_MF_WS", ["MF", "WS"], ["MF", "WF"]),
    ("fusion_WM", ["WS", "MS"], ["WF", "MF"]),
)


def build(counts, *, mu, beta, gamma, kappa, b, delta):
    mod = crn.Model("mtDNA_network")
    mod.add_species(dict(zip(NAMES, map(int, counts), strict=True)))
    mod.add_parameter(
        {
            "mu": mu,
            "beta": beta,
            "gamma": gamma,
            "kappa": kappa,
            "feedback": b,
            "delta": delta,
        }
    )
    for name, source, product in NETWORK:
        mod.add_reaction(crn.mass_action("gamma"), source, product, name=name)
    birth = crn.mass_action("max(0, mu + feedback*(kappa - WS - WF - delta*(MS + MF)))")
    for allele in ("W", "M"):
        s, f = f"{allele}S", f"{allele}F"
        mod.add_reaction(crn.mass_action("beta"), [f], [s], name=f"fission_{allele}")
        mod.add_reaction(crn.mass_action("mu"), [s], [], name=f"mitophagy_{allele}")
        for species in (s, f):
            mod.add_reaction(birth, [species], [f, f], name=f"replicate_{species}")
    return mod


def parameters(row):
    counts = [int(row[f"{x.lower()}_init"]) for x in NAMES]
    return counts, {
        p: float(row[p]) for p in ("mu", "beta", "gamma", "kappa", "b", "delta")
    }


def predicted_slope(counts, mu):
    n = sum(counts)
    h = (counts[1] + counts[3]) / n
    fs = sum(counts[:2]) / n
    return 2 * mu * fs * h * (1 - h) / n


def sample(row, hours, n, seed, threads):
    initial, params = parameters(row)
    model = build(initial, **params)
    result = crn.ensemble(
        model,
        end=float(hours[-1]),
        report_times=list(hours),
        seed=seed,
        replicates=n,
        threads=threads,
        method="static_direct",
        keep_paths="all",
        strict=True,
        max_output_bytes=128 * 2**20,
    ).simulate()
    paths = result.paths.sel(variable=list(NAMES)).values
    totals = paths.sum(axis=-1)
    mutant = paths[..., 1] + paths[..., 3]
    h = np.divide(mutant, totals, out=np.full_like(totals, np.nan), where=totals > 0)
    influence = (h - np.nanmean(h, axis=0)) ** 2
    w = np.asarray(hours) / np.dot(hours, hours)
    path_slope = influence @ w
    prediction = predicted_slope(initial, params["mu"])
    rng = np.random.default_rng(seed + 40000)
    # Bootstrap the centered squared deviations by whole paths. Retain the
    # centering correction from the resampled means as well.
    weights = rng.multinomial(n, np.full(n, 1 / n), size=600)
    finite = np.isfinite(h)
    sum_weights = weights @ finite.astype(float)
    average = np.divide(weights @ np.nan_to_num(h), sum_weights, where=sum_weights > 0)
    raw_second = (weights @ np.nan_to_num(h * h)) / sum_weights
    bootstrap = ((raw_second - average * average) @ w) / prediction
    lo, hi = np.quantile(bootstrap, [0.01 / 6, 1 - 0.01 / 6])
    trace = pd.DataFrame(
        {
            "hour": hours,
            "variance": np.nanvar(h, axis=0),
            "mean_copies": totals.mean(axis=0),
            "mean_heteroplasmy": np.nanmean(h, axis=0),
            "empty": (~finite).sum(axis=0),
            "fixed_fraction": np.mean((h == 0) | (h == 1), axis=0),
        }
    )
    summary = {
        "n": n,
        "initial_copies": sum(initial),
        "singleton_fraction": sum(initial[:2]) / sum(initial),
        "equation13": prediction,
        "slope": float(np.nanmean(path_slope)),
        "ratio": float(np.nanmean(path_slope) / prediction),
        "ratio_low": float(lo),
        "ratio_high": float(hi),
        "within_band": bool(lo >= 0.8 and hi <= 1.2),
    }
    return trace, summary


def compare_archive(sweep):
    path = ROOT / "reference-variance.csv"
    if not path.exists():
        return sweep
    rows = []
    for (panel, index), data in pd.read_csv(path).groupby(["panel", "index"]):
        series = data.query("0 <= hour <= 100")
        slopes = np.dot(series.hour, series.variance) / np.dot(series.hour, series.hour)
        rows.append((panel, index, slopes))
    return sweep.merge(
        pd.DataFrame(rows, columns=["panel", "index", "archived_slope_0_100h"]),
        how="left",
    )


def campaign(size, threads, output):
    params = pd.read_csv(ROOT / "conditions.csv")
    selected = M["primary"]
    if size == "quick":
        selected = ["H:60"]
    elif size == "draft":
        selected += ["H:55", "H:27", "H:82", "G:4"]
    summaries, traces = [], []
    for i, row in params.iterrows():
        key = f"{row.panel}:{row['index']}"
        init, settings = parameters(row)
        common = {"panel": row.panel, "index": int(row["index"]), "condition": key}
        if size != "full" and key not in selected:
            summaries.append(
                {**common, "n": 0, "equation13": predicted_slope(init, settings["mu"])}
            )
            continue
        hours = np.array(M["early_hours"])
        n = 256 if size == "quick" else 2000
        if size == "quick":
            hours = hours / 8
        if size == "full" and key not in M["primary"] and settings["beta"] >= 100:
            # Prespecified effort for the very fast switching control sweep.
            hours *= min(1.0, 100.0 / settings["beta"])
            n = 1024
        file = (
            output
            / "condition-checkpoints"
            / f"{row.panel}-{row['index']}-{size}-{n}.csv"
        )
        note = file.with_suffix(".json")
        file.parent.mkdir(parents=True, exist_ok=True)
        if file.exists() and note.exists():
            trace = pd.read_csv(file)
            row_summary = json.loads(note.read_text())
        else:
            trace, row_summary = sample(row, hours, n, M["seed"] + i, threads)
            tmp = file.with_suffix(".tmp")
            trace.to_csv(tmp, index=False)
            tmp.replace(file)
            tmp = note.with_suffix(".tmp")
            tmp.write_text(json.dumps(row_summary, sort_keys=True))
            tmp.replace(note)
        summaries.append({**common, **row_summary, "end_hour": float(hours[-1])})
        print(
            f"{key} n={n} at {hours[-1]:g} h, slope ratio {row_summary['ratio']:.3f}",
            flush=True,
        )
        traces.append(trace.assign(**common))
    if size != "quick":
        row = params[(params.panel == "H") & (params["index"] == 60)].iloc[0]
        trace, _ = sample(
            row,
            np.array([0.0, 1.0, 2.0, 5.0, 10.0, 25.0, 50.0]),
            96,
            M["seed"] + 701,
            threads,
        )
        trace.to_csv(output / "nominal-long.csv", index=False)
    return compare_archive(pd.DataFrame(summaries)), pd.concat(
        traces, ignore_index=True
    )


def plot(output, folder):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    table = pd.read_csv(output / "slopes.csv")
    trace = pd.read_csv(output / "variance.csv")
    fig, axes = plt.subplots(2, 2, figsize=(9, 6.5), layout="constrained")
    representatives = {"E": [0, 3], "F": [0, 4], "G": [0, 4, 8], "H": [5, 60, 115]}
    names = {
        "E": "Turnover",
        "F": "Copy number",
        "G": "Mutant sensing",
        "H": "Network ratio",
    }
    for ax, (panel, indices) in zip(axes.flat, representatives.items(), strict=True):
        for index in indices:
            condition = f"{panel}:{index}"
            group = trace[trace.condition == condition]
            if group.empty:
                continue
            slope = table.loc[table.condition == condition, "equation13"].iloc[0]
            line = ax.plot(
                group.hour, group.variance, marker="o", ms=3, label=condition
            )[0]
            ax.plot(
                group.hour,
                group.hour * slope,
                color=line.get_color(),
                linestyle="--",
                alpha=0.8,
            )
        ax.set(
            title=f"{panel} · {names[panel]}",
            xlabel="Time (hours)",
            ylabel="Variance of mutant fraction",
        )
        if ax.lines:
            ax.legend(frameon=False, fontsize=8)
    folder.mkdir(parents=True, exist_ok=True)
    fig.savefig(folder / "variance.png", dpi=155)
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--size", choices=("quick", "draft", "full"), default="quick")
    p.add_argument("--output", type=Path)
    p.add_argument("--threads", type=int, default=1)
    p.add_argument("--figures", type=Path)
    args = p.parse_args()
    args.output = args.output or ROOT / (
        "results" if args.figures else f"runs/{args.size}"
    )
    if args.figures:
        plot(args.output, args.figures)
        return
    args.output.mkdir(parents=True, exist_ok=True)
    table, trace = campaign(args.size, args.threads, args.output)
    table.to_csv(args.output / "slopes.csv", index=False)
    trace.to_csv(args.output / "variance.csv", index=False)
    print(
        table.query("n > 0")[
            ["condition", "ratio", "ratio_low", "ratio_high"]
        ].to_string(index=False)
    )


if __name__ == "__main__":
    main()
