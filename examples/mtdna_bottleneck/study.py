"""Reconstruct copy number and heteroplasmy variance in a developing mouse."""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import differential_evolution

from model import accept_distance, chain, ensemble_chain, forecast_sample, predictions

ROOT = Path(__file__).resolve().parent
M = json.loads((ROOT / "manifest.json").read_text())
# Uniform priors: replication and degradation rates in [0, 1] per hour,
# initial copies in (0, 1e6], quiescent phase lengths in [0, 50] days.
# Phase six has no replication; by default its degradation is fixed to
# zero as well (--free-final-degradation frees it, as in the appendix's
# prior). Search bounds for starting points merely put the initial
# candidates near measurable trajectories.
FREE = np.array([0, 1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 13, 14, 15])
FREE_FINAL = np.array([0, 1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14, 15])
START_BOUNDS = [
    (5e4, 2.9e5),
    (0.052, 0.095),
    (0.026, 0.069),
    (-0.012, 0.025),
    (-0.009, 0.015),
    (-0.004, 0.009),
    *[(0.01, 0.82)] * 5,
    (1, 15),
    (5, 22),
    (15, 50),
]
SCALES = np.array(
    [
        6500,
        0.0018,
        0.005,
        0.0028,
        0.0016,
        0.001,
        0,
        0.04,
        0.04,
        0.045,
        0.045,
        0.045,
        0,
        1.0,
        1.8,
        2.0,
    ]
)


def data():
    t = pd.read_csv(ROOT / "reference-table1.csv").sort_values("day")
    kinds = (t.observable == "normalized_variance").to_numpy(dtype=np.int64)
    raw = t.value.to_numpy()
    return t, t.day.to_numpy(), kinds, np.where(kinds, raw, np.log(raw))


def expand(values):
    theta = np.zeros(16)
    theta[FREE if len(values) == len(FREE) else FREE_FINAL] = values
    return theta


def start(age, kind, observed, seed, free_final=False):
    def cost(values):
        return accept_distance(expand(values), age, kind, observed)

    bounds = list(START_BOUNDS)
    if free_final:
        bounds.insert(11, (0.0, 0.002))
    result = differential_evolution(
        cost, bounds, maxiter=250, popsize=9, rng=seed, polish=False
    )
    if result.fun > M["figure3_threshold"]:
        raise RuntimeError(
            f"Could not initialize an ABC chain at the declared threshold: {result.fun}"
        )
    return expand(result.x), float(result.fun)


def ess(x):
    x = np.asarray(x) - np.mean(x)
    if np.dot(x, x) == 0:
        return 0.0
    power = np.abs(np.fft.rfft(x, 2 * len(x))) ** 2
    corr = np.fft.irfft(power)[: len(x)]
    corr /= corr[0]
    pairs = corr[1 : len(x) - 1 : 2] + corr[2 : len(x) : 2]
    first = np.flatnonzero(pairs <= 0)
    length = first[0] if len(first) else len(pairs)
    return min(len(x), len(x) / (1 + 2 * pairs[:length].sum()))


def hb_comparison(draws, seed):
    cells = pd.read_csv(ROOT / "reference-hb.csv", dtype={"litter": str})
    table = pd.read_csv(ROOT / "reference-table2.csv", dtype={"litter": str})
    rows = []
    rng = np.random.default_rng(seed)
    for record in table.itertuples(index=False):
        h = cells.loc[cells.group == record.group, "heteroplasmy"].to_numpy()
        mean, var = h.mean(), h.var(ddof=1)
        ratio = var / (mean * (1 - mean))
        pred = np.array(
            [predictions(x, np.array([float(record.day)]))[0][0, 1] for x in draws]
        )
        band = np.quantile(pred, [0.025, 0.5, 0.975])
        variance_sample = []
        for i in rng.integers(len(draws), size=400):
            parameter = draws[i]
            moments = predictions(parameter, np.array([float(record.day)]))[0][0]
            cell_h = forecast_sample(parameter[0], mean, *moments, len(h), rng)
            h_bar = np.nanmean(cell_h)
            if 0 < h_bar < 1:
                variance_sample.append(
                    np.nanvar(cell_h, ddof=1) / (h_bar * (1 - h_bar))
                )
        low, high = np.quantile(variance_sample, [0.025, 0.975])
        rows.append(
            {
                "group": record.group,
                "litter": record.litter,
                "day": record.day,
                "n": len(h),
                "holdout": record.litter in M["hb_holdout_litters"],
                "mean": mean,
                "variance": var,
                "normalized_variance": ratio,
                "reference_mean": record.mean,
                "reference_variance": record.variance,
                "reference_normalized": record.normalized_variance,
                "prediction_low": band[0],
                "prediction_median": band[1],
                "prediction_high": band[2],
                "expected_sample_low": low,
                "expected_sample_high": high,
                "in_sample_band": low <= ratio <= high,
            }
        )
    return pd.DataFrame(rows)


def run(size, output, free_final=False):
    _, age, kinds, observed = data()
    width = 16
    sweeps = {"quick": 200, "draft": 15000, "full": 360000}[size]
    count = 1 if size == "quick" else M["abc_chains"]
    changes = SCALES.copy()
    if free_final:
        changes[12] = 0.0005
    indices = np.flatnonzero(changes)
    ages = np.unique(np.r_[np.linspace(0, 70, 141), 100.0, M["reference_ages"]])
    draw_rows, diag = [], []
    for i in range(count):
        checkpoint = output / f"ensemble-{sweeps}-{i}.npz"
        if checkpoint.exists():
            with np.load(checkpoint, allow_pickle=False) as saved:
                states, distances, accepted, first = (
                    saved[k] for k in ("states", "distances", "accepted", "first")
                )
        else:
            initial, first = start(age, kinds, observed, 61 + i, free_final)
            spread, _, _ = chain(
                initial,
                160000 if size != "quick" else 2000,
                changes,
                indices,
                age,
                kinds,
                observed,
                float(M["figure3_threshold"]),
                M["seed"] + 100 * i,
            )
            population = spread[np.linspace(0, len(spread) - 1, width, dtype=int)]
            states, distances, accepted = ensemble_chain(
                population,
                sweeps,
                changes,
                indices,
                age,
                kinds,
                observed,
                float(M["figure3_threshold"]),
                M["seed"] + 100 * i + 1,
            )
            temp = checkpoint.with_name(checkpoint.stem + ".tmp.npz")
            np.savez_compressed(
                temp, states=states, distances=distances, accepted=accepted, first=first
            )
            temp.replace(checkpoint)
        states, distances = states[len(states) // 4 :], distances[len(states) // 4 :]
        if np.any(distances > M["figure3_threshold"]):
            raise RuntimeError("Chain includes a state outside the threshold")
        selected = np.array([10.0, 35.0, 70.0, 100.0])
        values = np.asarray(
            [predictions(x, selected)[0] for x in states.reshape(-1, len(changes))]
        )
        # Report ESS using the variance of each ensemble mean relative to
        # the marginal variance over walkers, with temporal correlation.
        series = values.reshape(len(states), width, -1).mean(axis=1)
        marginal = values.reshape(-1, series.shape[1]).var(axis=0)
        ratios = np.divide(
            marginal,
            series.var(axis=0),
            out=np.zeros_like(marginal),
            where=series.var(axis=0) > 0,
        )
        effective = [ess(v) * ratio for v, ratio in zip(series.T, ratios, strict=True)]
        diag.append(
            {
                "chain": i,
                "sweeps": sweeps,
                "walkers": width,
                "accepted": int(accepted),
                "acceptance": int(accepted) / (sweeps * width),
                "initial_distance": float(first),
                "min_distance": float(distances.min()),
                "max_distance": float(distances.max()),
                "ess_min_observable": min(effective),
            }
        )
        rng = np.random.default_rng(M["seed"] + i)
        flat = states.reshape(-1, len(changes))
        pick = rng.choice(len(flat), size=min(1600, len(flat)), replace=False)
        draw_rows.extend(flat[pick])
        print(json.dumps(diag[-1]), flush=True)
    samples = np.array(draw_rows)
    model = np.array([predictions(theta, ages)[0] for theta in samples])
    model[:, :, 0] = np.exp(model[:, :, 0])
    quantile = np.quantile(model, [0.025, 0.25, 0.5, 0.75, 0.975], axis=0)
    rows = [
        (day, name, *quantile[:, i, j])
        for i, day in enumerate(ages)
        for j, name in enumerate(("copies", "normalized_variance"))
    ]
    trajectories = pd.DataFrame(
        rows, columns=["day", "observable", "q025", "q25", "median", "q75", "q975"]
    )
    trajectories.to_csv(output / "trajectories.csv", index=False)
    reference = pd.read_csv(ROOT / "reference-figure3.csv", dtype={"day": float}).merge(
        trajectories
    )
    reference["center_in_95"] = (reference.q025 <= reference.center) & (
        reference.center <= reference.q975
    )
    reference["difference"] = np.where(
        reference.observable == "copies",
        np.log(reference["median"] / reference.center),
        reference["median"] - reference.center,
    )
    reference.to_csv(output / "figure3-comparison.csv", index=False)
    diagnostic = pd.DataFrame(diag)
    diagnostic.to_csv(output / "abc.csv", index=False)
    hb_comparison(samples, M["seed"] + 99).to_csv(
        output / "hb-comparison.csv", index=False
    )
    print(diagnostic.to_string(index=False))


def plot(output, folder):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    sample = pd.read_csv(output / "trajectories.csv")
    reference = pd.read_csv(ROOT / "reference-table1.csv")
    paper_curve = pd.read_csv(ROOT / "reference-figure3.csv")
    folder.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 1, figsize=(8, 6.4), layout="constrained", sharex=True)
    for ax, observable, ylabel in zip(
        axes,
        ["copies", "normalized_variance"],
        ["Mean copies", "Normalized variance"],
        strict=True,
    ):
        s = sample[sample.observable == observable]
        r = reference[reference.observable == observable]
        p = paper_curve[paper_curve.observable == observable]
        ax.fill_between(
            s.day,
            s.q025,
            s.q975,
            color="#93b8d3",
            alpha=0.5,
            label="Reconstruction, 95% band",
        )
        ax.plot(s.day, s["median"], color="#234a74", label="Reconstruction, median")
        ax.scatter(
            p.day,
            p.center,
            marker="x",
            color="#ad563e",
            label="Published curve, digitized",
        )
        ax.scatter(
            r.day, r.value, s=11, color=".27", alpha=0.7, label="Table 1 measurements"
        )
        ax.set(ylabel=ylabel)
    axes[0].set_yscale("log")
    axes[0].legend(frameon=False, fontsize=8, ncol=2)
    axes[-1].set_xlabel("Days post-conception")
    fig.savefig(folder / "bottleneck.png", dpi=155)
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--size", choices=("quick", "draft", "full"), default="quick")
    p.add_argument("--output", type=Path)
    p.add_argument("--figures", type=Path)
    p.add_argument(
        "--free-final-degradation",
        action="store_true",
        help="sensitivity run: degradation in phase six free on [0, 1] per hour",
    )
    a = p.parse_args()
    suffix = "-free-final" if a.free_final_degradation else ""
    a.output = a.output or ROOT / ("results" if a.figures else f"runs/{a.size}{suffix}")
    if a.figures:
        plot(a.output, a.figures)
        return
    a.output.mkdir(parents=True, exist_ok=True)
    run(a.size, a.output, a.free_final_degradation)


if __name__ == "__main__":
    main()
