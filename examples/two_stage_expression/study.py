"""Reproduce the protein-count comparisons and generate their figures."""

import argparse
import json
from pathlib import Path

import crn_craft as crn
import numpy as np
import pandas as pd
from scipy.sparse.linalg import expm_multiply

from model import (
    build,
    comparison,
    joint_generator,
    moments,
    negative_binomial,
    protein_pmf,
)

ROOT = Path(__file__).resolve().parent
M = json.loads((ROOT / "manifest.json").read_text())


def stationary_case(a, b, gamma):
    horizon = 32 * max(1.0, 1 / gamma)
    coarse = protein_pmf(a, b, gamma, horizon, 2048, 10)
    long = protein_pmf(a, b, gamma, 2 * horizon, 2048, 10)
    fine = protein_pmf(a, b, gamma, 2 * horizon, 4096, 14)
    x = np.arange(len(fine))
    mu = np.dot(x, fine)
    fano = np.dot((x - mu) ** 2, fine) / mu
    true_mean, true_fano = moments(a, b, gamma)
    short_gap = np.max(np.abs(np.cumsum(coarse) - np.cumsum(long)))
    refinement = np.max(np.abs(np.cumsum(long) - np.cumsum(fine)[:2048]))
    tail = fine[2048:].sum()
    q = negative_binomial(a, b, len(fine))
    kl, tv = comparison(fine, q)
    passed = (
        short_gap < M["cdf_double_horizon_tolerance"]
        and refinement < M["numerical_tail_target"]
        and tail < M["numerical_tail_target"]
        and abs(mu / true_mean - 1) < M["stationary_relative_moment_tolerance"]
        and abs(fano / true_fano - 1) < M["stationary_relative_moment_tolerance"]
    )
    return (
        fine,
        q,
        {
            "horizon": horizon,
            "mean": mu,
            "mean_exact": true_mean,
            "fano": fano,
            "fano_exact": true_fano,
            "horizon_cdf_gap": short_gap,
            "refinement_cdf_gap": refinement,
            "tail_beyond_2047": tail,
            "distribution_tail_gt150": fine[151:].sum(),
            "kl": kl,
            "tv": tv,
            "numerical_checks": passed,
        },
    )


def native_check():
    a, b, gamma = 0.7, 1.2, 2.0
    bounded = build(a, b, gamma).fsp(1.0, [1.0], {"M": 5, "P": 8}, tolerance=1e-12)
    meta = bounded.metadata
    obs = meta["snapshots"][0]
    q = joint_generator(a, b, gamma, 5, 8)
    init = np.zeros(q.shape[0])
    init[0] = 1
    projected = expm_multiply(q, init)
    # CRN-craft orders labels lexicographically and states by discovery.
    ids = [m * 9 + p for m, p in meta["states"]]
    gap = np.max(np.abs(projected[ids] - obs["probabilities"]))
    sink_gap = abs(projected[-1] - obs["physical_leakage_estimate"])
    return {"max_gap": float(gap), "sink_gap": float(sink_gap)}


def narrow_panel_fsp():
    """Compare an entire production panel with CRN-craft's bounded solver."""
    bound = {"M": 50, "P": 106}
    result = build(20, 2.5, 1).fsp(3, [3], bound, max_states=5500, tolerance=5e-10)
    meta = result.metadata
    snap = meta["snapshots"][0]
    proteins = np.asarray(meta["states"])[:, 1]
    projected = np.bincount(proteins, weights=snap["probabilities"])
    exact = protein_pmf(20, 2.5, 1, 3.0)
    return {
        "max_marginal_cdf_gap": float(
            np.max(np.abs(np.cumsum(projected) - np.cumsum(exact)[: len(projected)]))
        ),
        "physical_leakage": snap["physical_leakage_estimate"],
        "series_tail_upper_bound": snap["series_tail_upper_bound"],
    }


def stochastic_case(a, b, gamma, horizon, n, seed, threads):
    run = crn.ensemble(
        build(a, b, gamma),
        end=horizon,
        report_times=[horizon],
        replicates=n,
        threads=threads,
        seed=seed,
        method="static_direct",
        keep_paths="all",
        strict=True,
    )
    sample = run.simulate().paths.sel(variable="P").values.ravel().astype(int)
    return np.bincount(sample)


def run(size, output, threads):
    check = native_check()
    if max(check.values()) > 1e-9:
        raise RuntimeError("Native FSP differs from the independent matrix exponential")
    if size == "quick":
        pc = protein_pmf(0.7, 1.2, 2, 1.0, 256)
        hist = stochastic_case(0.7, 1.2, 2, 1.0, 1000, M["seed"], threads)
        gap = np.max(np.abs(np.cumsum(hist) / hist.sum() - np.cumsum(pc)[: len(hist)]))
        print(json.dumps({**check, "sample_cdf_gap": float(gap)}, sort_keys=True))
        return
    rows, profiles = [], []
    main = M["gamma_figure_bc"]
    gamma_values = (
        main
        if size == "draft"
        else sorted(set(main + [10**x for x in M["log10_gamma_figure_d"]]))
    )
    for i, (a, b) in enumerate(M["regimes"]):
        for j, gamma in enumerate(gamma_values):
            p, q, row = stationary_case(a, b, gamma)
            row.update(a=a, b=b, gamma=gamma, n=0)
            if gamma in main:
                hist = stochastic_case(
                    a, b, gamma, row["horizon"], 10000, M["seed"] + i * 100 + j, threads
                )
                empirical = np.pad(hist / hist.sum(), (0, max(0, len(p) - len(hist))))[
                    : len(p)
                ]
                gap = abs(np.cumsum(empirical) - np.cumsum(p)).max()
                bound = np.sqrt(np.log(2 * 4 / 0.01) / (2 * hist.sum()))
                row.update(
                    n=int(hist.sum()),
                    cdf_gap=gap,
                    dkw_bound=bound,
                    within_dkw=gap <= bound,
                )
                profiles.extend(
                    (a, b, gamma, k, p[k], q[k], empirical[k]) for k in range(151)
                )
                profiles.append(
                    (
                        a,
                        b,
                        gamma,
                        151,
                        p[151:].sum(),
                        q[151:].sum(),
                        empirical[151:].sum(),
                    )
                )
            rows.append(row)
            print(
                f"{a:g},{b:g},{gamma:g}: KL {row['kl']:.6g}, fano {row['fano']:.5g}",
                flush=True,
            )
    pd.DataFrame(rows).to_csv(output / "divergence.csv", index=False)
    pd.DataFrame(
        profiles,
        columns=[
            "a",
            "b",
            "gamma",
            "count",
            "reaction_pmf",
            "negative_binomial",
            "ssa",
        ],
    ).to_csv(output / "histograms.csv", index=False)
    check["figure1b_gamma1_age3"] = narrow_panel_fsp()
    (output / "native-check.json").write_text(json.dumps(check, indent=2) + "\n")


def plot(output, folder):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    full = pd.read_csv(output / "histograms.csv")
    above = full[full["count"] == 151].set_index(["a", "gamma"]).reaction_pmf
    hist = full.query("count <= 150")
    divergence = pd.read_csv(output / "divergence.csv").sort_values("gamma")
    folder.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.0), layout="constrained")
    for ax, ((a, b), group) in zip(
        axes, hist.groupby(["a", "b"], sort=False), strict=True
    ):
        for gamma, s in group.groupby("gamma"):
            ax.plot(
                s["count"],
                s.reaction_pmf,
                label=f"Reaction law, γ={gamma:g}, {above[a, gamma]:.1%} above 150",
            )
            ax.scatter(s["count"].iloc[::5], s.ssa.iloc[::5], s=5, alpha=0.55)
        ax.plot(
            s["count"],
            s.negative_binomial,
            color="0.35",
            linestyle="--",
            label=f"Negative binomial, {1 - s.negative_binomial.sum():.1%} above 150",
        )
        ax.set(
            title=f"a={a:g}, b={b:g}, mean 50",
            xlabel="Protein copies",
            ylabel="Probability",
        )
        ax.legend(frameon=False, fontsize=8)
    fig.savefig(folder / "same-mean.png", dpi=155)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(7.1, 4.2), layout="constrained")
    archive = pd.read_csv(ROOT / "reference-figure1d.csv")
    sampled = pd.read_csv(output / "review-figure1d-conventions.csv")
    for (a, b), series in divergence.groupby(["a", "b"]):
        line = ax.plot(
            series.gamma,
            series.kl / np.log(2),
            marker="o",
            label=f"a={a:g}, b={b:g}, exact law",
        )[0]
        draws = sampled[sampled.a == a]
        ax.plot(
            draws.gamma,
            draws.bits_sampled_10000_mean,
            linestyle="--",
            color=line.get_color(),
            label=f"a={a:g}, 10,000 simulated draws",
        )
        values = archive[(archive.a == a) & (archive.b == b)]
        ax.scatter(
            values.gamma,
            values.published_kl_display,
            facecolors="none",
            edgecolors=line.get_color(),
            label=f"Published, a={a:g}",
            s=38,
        )
    ax.set(
        xscale="log",
        xlabel="mRNA-to-protein decay-rate ratio γ",
        ylabel="KL divergence (bits)",
    )
    ax.legend(frameon=False)
    fig.savefig(folder / "divergence.png", dpi=155)
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--size", choices=("quick", "draft", "full"), default="quick")
    p.add_argument("--threads", type=int, default=1)
    p.add_argument("--output", type=Path)
    p.add_argument("--figures", type=Path)
    a = p.parse_args()
    a.output = a.output or ROOT / ("results" if a.figures else f"runs/{a.size}")
    if a.figures:
        plot(a.output, a.figures)
        return
    a.output.mkdir(parents=True, exist_ok=True)
    run(a.size, a.output, a.threads)


if __name__ == "__main__":
    main()
