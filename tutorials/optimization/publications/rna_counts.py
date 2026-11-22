"""Fit the RNA-count distributions from Poovathingal and Gunawan, case 1."""

import crn_craft as crn
import numpy as np
from scipy.special import ndtr

from common import checked, write
from .search import PACKAGES as PACKAGES
from .search import PENALTY, check, save_plot, search

DOI = "10.1186/1471-2105-11-414"
NAMES = ("k1", "k2", "k3")
TRUTH = (0.0277, 0.1667, 0.4)
TABLE = {
    10: [0.0200, 0.1029, 0.3643],
    20: [0.0371, 0.2124, 0.5263],
    100: [0.0273, 0.1702, 0.4121],
    10_000: [0.0276, 0.1659, 0.4102],
}
TIMES = np.arange(0.5, 75.01, 0.5)
SIGMA = 0.5  # N(0, 0.25) with 0.25 the variance.


def model():
    net = crn.Model("rna_count_transcription")
    net.add_species({"silent": 1, "active": 0, "rna": 0})
    for name, value in zip(NAMES, TRUTH, strict=True):
        net.add_parameter(name, value)
    net.add_reaction("k1", ["silent"], ["active"], name="activate")
    net.add_reaction("k2", ["active"], ["silent"], name="deactivate")
    net.add_reaction("k3", ["active"], ["active", "rna"], name="transcribe")
    return net


def paths(net, rates, size, seed, times):
    source = crn.ensemble(
        net,
        end=float(times[-1]),
        report_times=times.tolist(),
        replicates=size,
        seed=seed,
        keep_paths="all",
        quantiles=(),
        method="static_direct",
        strict=True,
    )
    result = checked(
        source.bind_parameters(dict(zip(NAMES, rates, strict=True))).simulate()
    )
    return np.asarray(result.paths.sel(variable="rna"), dtype=np.int32)


def observed(net, rates, size, seed, times):
    rng = np.random.default_rng(seed)
    latent = paths(net, rates, size, seed, times)
    return latent + np.rint(rng.normal(0, SIGMA, latent.shape)).astype(np.int32)


class CdfScore:
    """Equation 7 with plotting positions from Equation 9 for sparse samples."""

    def __init__(self, observations, *, population=False, midrank=False):
        self.midrank = midrank
        self.low = int(observations.min())
        self.high = int(observations.max())
        size, times = observations.shape
        # Collapse repeated thresholds without changing their rank weights.
        coefficients = np.zeros((3, times, self.high - self.low + 1))
        for i, column in enumerate(observations.T):
            x = np.sort(column)
            if midrank:
                x, counts = np.unique(x, return_counts=True)
                rank = (np.cumsum(counts) - 0.5 * counts) / size
            elif population:
                x, counts = np.unique(x, return_counts=True)
                rank = np.cumsum(counts)[:-1] / size
                x = x[:-1]
            else:
                rank = (np.arange(size - 1) + 0.5) / size
                x = x[:-1]
            w = size / (rank * (1 - rank))
            for j in range(3):
                np.add.at(coefficients[j, i], x - self.low, w * rank**j)
        self.w, self.wf, self.wf2 = coefficients

    def __call__(self, predicted):
        predicted = np.asarray(predicted)
        return float(np.sum(self.w * predicted**2 - 2 * self.wf * predicted + self.wf2))

    def from_counts(self, counts):
        support = np.arange(self.low, self.high + 1)
        cdf = []
        for x in counts.T:
            ordered = np.sort(x)
            right = np.searchsorted(ordered, support, side="right")
            if self.midrank:
                left = np.searchsorted(ordered, support, side="left")
                right = 0.5 * (left + right)
            cdf.append(right / len(x))
        cdf = np.array(cdf)
        return self(cdf)


def distributions(graph, rates, times):
    solution = graph.solve(
        float(times[-1]),
        times.tolist(),
        parameters=dict(zip(NAMES, rates, strict=True)),
        tolerance=1e-10,
    )
    if solution.snapshots[-1].physical_leakage > 1e-6:
        raise ValueError("RNA cutoff lost more than 1e-6 of the law")
    return np.array([solution.marginal("rna", at=i) for i in range(len(times))])


def objective(graph, data, times):
    count_grid = np.arange(int(graph.states[:, graph.species.index("rna")].max()) + 1)
    readout = ndtr(
        (np.arange(data.low, data.high + 1)[None, :] + 0.5 - count_grid[:, None])
        / SIGMA
    )
    if data.midrank:
        left = ndtr(
            (np.arange(data.low, data.high + 1)[None, :] - 0.5 - count_grid[:, None])
            / SIGMA
        )
        readout = 0.5 * (left + readout)
    return lambda rates: data(distributions(graph, rates, times) @ readout)


def options(parser):
    parser.add_argument(
        "--source-budget",
        action="store_true",
        help="Run the publication-scale DE/SSA arm for m=100.",
    )
    parser.add_argument("--source-evaluations", type=int, default=120_030)
    parser.add_argument(
        "--published-rank",
        action="store_true",
        help="Use the individually ranked Equation 9 thresholds in every new fit.",
    )
    parser.add_argument(
        "--evaluations",
        type=int,
        help="Override the faster search evaluations per start.",
    )
    parser.add_argument(
        "--datasets",
        type=int,
        default=8,
        help="Regenerated 100-trajectory datasets in the study (default 8).",
    )
    parser.add_argument(
        "--starts",
        type=int,
        default=2,
        help="Optimizer starts per 100-trajectory dataset in the study (default 2).",
    )


def run(args, report):
    """Regenerate the RNA data, check the FSP against SSA, and recover rates."""
    net = model()
    write(args.output / "model.json", net.to_dict())
    fast = args.preset == "smoke"
    seed = 5193 + args.seed_offset
    samples = (100,) if fast else (10, 20, 100, 10_000)
    # 300 calls (ten DE generations) left 11 of 16 delivered starts above the
    # score of the generating rates on the same data; 3,000 reached it on all
    # four datasets of the review rerun.
    budget = args.evaluations or (40 if fast else 3000)
    times = TIMES if not fast else TIMES[::10].copy()
    if fast:
        times[-1] = TIMES[-1]
    graph = crn.prepare_fsp(net, {"silent": 1, "active": 1, "rna": 510})
    chosen = np.array([10.0, 30.0, 75.0])
    independent = observed(net, TRUTH, 5000 if fast else 10_000, seed + 300, chosen)
    pmf = distributions(graph, TRUTH, chosen)
    counts = np.arange(-2, 121)
    noised = pmf @ ndtr(
        (counts[None, :] + 0.5 - np.arange(pmf.shape[1])[:, None]) / SIGMA
    )
    empirical = np.array(
        [(independent[:, j, None] <= counts).mean(axis=0) for j in range(len(chosen))]
    )
    gap = float(np.max(np.abs(noised - empirical)))
    dkw = float(np.sqrt(np.log(2 * len(chosen) / 0.01) / (2 * len(independent))))
    check(report, "noise-adjusted-FSP-and-SSA", gap <= dkw, cdf_gap=gap, dkw_99=dkw)
    example = CdfScore(np.array([[0], [1], [2]]))
    target = sum(
        3 * (r - y) ** 2 / (r * (1 - r)) for r, y in [(1 / 6, 0.2), (0.5, 0.6)]
    )
    check(report, "two-rank-score", abs(example([[0.2, 0.6, 1.0]]) - target) < 1e-12)
    bounds = [(0.0, 5.0)] * 3
    results = []
    reproducible = {}
    generated = 0
    for size in samples:
        for data_index in range(1 if fast or size != 100 else args.datasets):
            data_seed = seed + size + data_index * 20000
            observations = observed(net, TRUTH, size, data_seed, times)
            generated += size
            if data_index == 0:
                reproducible[str(size)] = observations.tolist()
            problem = CdfScore(
                observations,
                population=size >= 10_000,
                midrank=not args.published_rank and size < 10_000,
            )
            f = objective(graph, problem, times)
            for retry in range(1 if fast or size != 100 else args.starts):
                arm = f"m{size}-data{data_index}-de{retry}"
                # Unit-box cubing is a bijection onto [0, 5] and allocates
                # more early proposals near the lower end of the rate box.
                fit = search(
                    args,
                    arm,
                    "de",
                    lambda u, f=f: f(5 * u**3),
                    [(0.0, 1.0)] * 3,
                    data_seed + 91 + retry,
                    budget,
                    pop=10 if fast else 30,
                )
                fitted = 5 * fit.best**3
                at_truth = f(TRUTH)
                held = problem.from_counts(
                    observed(
                        net,
                        fitted,
                        1500 if fast else 4000,
                        data_seed + 421 + retry,
                        times,
                    )
                )
                results.append(
                    dict(
                        cells=size,
                        data_seed=data_seed,
                        rates=fitted,
                        score=fit.value,
                        held_ssa_score=held,
                        truth_score=at_truth,
                        # False means the search stopped short of the optimum.
                        reached_truth_score=fit.value <= at_truth,
                        paper=TABLE[size],
                        evaluations=len(fit.candidates),
                        failures=fit.failures,
                    )
                )
    if args.source_budget:
        observations = observed(net, TRUTH, 100, seed + 100, TIMES)
        problem = CdfScore(observations, midrank=not args.published_rank)
        prediction_seed = seed + 200
        fit = search(
            args,
            "source-ssa",
            "de",
            lambda k: problem.from_counts(
                observed(net, k, 10_000, prediction_seed, TIMES)
            ),
            bounds,
            seed + 201,
            args.source_evaluations,
            pop=30,
        )
        held = problem.from_counts(observed(net, fit.best, 10_000, seed + 301, TIMES))
        report["source_ssa"] = dict(
            rates=fit.best,
            score=fit.value,
            independent_score=held,
            truth_score=problem.from_counts(
                observed(net, TRUTH, 10_000, prediction_seed, TIMES)
            ),
            sparse_score="published ranks" if args.published_rank else "midranks",
            candidate_paths=10_000 * len(fit.candidates),
            evaluations=len(fit.candidates),
        )
    # Figure 2A: ten trajectories, first rate at truth, the score of the fits.
    contour_times = times
    small = observed(net, TRUTH, 10, seed + 10, contour_times)
    f = objective(
        graph, CdfScore(small, midrank=not args.published_rank), contour_times
    )
    grid = np.linspace(0.1, 1.0, 5 if fast else 31)
    cost = np.array([[f((TRUTH[0], a, b)) for a in grid] for b in grid])
    report.update(
        doi=DOI,
        fits=results,
        truth=TRUTH,
        time_points=len(times),
        generated_paths=generated,
        search_coordinates="rates = 5*u^3",
        sparse_score="published ranks" if args.published_rank else "discrete midranks",
        starts_reaching_truth_score=sum(r["reached_truth_score"] for r in results),
        measurement="rounded normal, signed observations",
        source_sample="regenerated",
    )
    check(report, "search-ran", all(r["score"] < PENALTY for r in results))
    write(args.output / "observations.json", dict(times=times, sizes=reproducible))
    write(
        args.output / "cdf-contour.json",
        dict(grid=grid, loss=cost, k1=TRUTH[0], source_seed=seed + 10),
    )

    def draw(fig):
        ax = fig.subplots()
        fig.colorbar(
            ax.contourf(grid, grid, np.log10(np.maximum(cost, 1e-12))),
            ax=ax,
            label="log10 CDF score",
        )
        ax.scatter(
            [TRUTH[1]],
            [TRUTH[2]],
            color="white",
            edgecolors="black",
            label="generating rates",
        )
        ax.set(xlabel="Switching-off rate", ylabel="RNA production rate")
        ax.legend()

    save_plot(args, "cdf-contour.png", draw)


def draw(tables, folder):
    """Draw C5-contour.png and C5-recovery.png from results/C5 tables."""
    import csv
    from pathlib import Path

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    tables, folder = Path(tables), Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    with (tables / "recovery.csv").open(newline="") as stream:
        fits = list(csv.DictReader(stream))

    def rates(row):
        return float(row["k2"]), float(row["k3"])

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.4), layout="constrained")
    for ax, cells in zip(axes, (10, 100), strict=True):
        with (tables / f"contour-m{cells}.csv").open(newline="") as stream:
            grid = np.array(
                [[float(v) for v in row.values()] for row in csv.DictReader(stream)]
            )
        k2, k3 = np.unique(grid[:, 0]), np.unique(grid[:, 1])
        score = grid[:, 2].reshape(len(k3), len(k2))
        shade = ax.contourf(
            k2, k3, np.log10(score / score.min()), levels=12, cmap="Blues_r"
        )
        fig.colorbar(shade, ax=ax, label="log10(score / grid minimum)")
        i, j = np.unravel_index(np.argmin(score), score.shape)
        ax.plot(
            k2[j], k3[i], "*", ms=13, color="white", mec="black", label="grid minimum"
        )
        ax.plot(
            TRUTH[1],
            TRUTH[2],
            "o",
            ms=9,
            color="white",
            mec="black",
            label="generating rates",
        )
        ax.plot(
            *TABLE[cells][1:],
            "P",
            ms=10,
            color="#eb6834",
            mec="black",
            label="paper's fit, Table 1 (own data)",
        )
        for row in fits:
            if (row["run"], row["cells"], row["data_index"]) == (
                "review-3000",
                str(cells),
                "0",
            ) and 0.1 <= min(rates(row)) <= max(rates(row)) <= 1:
                ax.plot(
                    *rates(row),
                    "D",
                    ms=8,
                    color="#1baf7a",
                    mec="black",
                    label="our fit, 3,000 calls (k1 free)",
                )
        ax.set(
            xlim=(0.1, 1),
            ylim=(0.1, 1),
            xlabel="k2, switching off (1/min)",
            ylabel="k3, RNA production (1/min)",
        )
        ax.set_title(
            f"{cells} regenerated trajectories, k1 = 0.0277 fixed", fontsize=10
        )
        ax.legend(fontsize=7, loc="lower right")
    fig.suptitle(
        "Midrank CDF score over (k2, k3), the setting of Figure 2A", fontsize=10
    )
    fig.savefig(folder / "C5-contour.png", dpi=155)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.4, 4.6), layout="constrained")
    styles = (
        (
            "delivered-300",
            "published",
            "s",
            "#eb6834",
            "none",
            "300 calls, individual ranks",
        ),
        ("delivered-300", "midrank", "o", "#2a78d6", "none", "300 calls, midranks"),
        ("review-3000", "midrank", "D", "#1baf7a", "#1baf7a", "3,000 calls, midranks"),
    )
    for run, variant, marker, color, face, label in styles:
        chosen = [
            rates(row)
            for row in fits
            if (row["run"], row["sparse_score"], row["cells"]) == (run, variant, "100")
        ]
        ax.plot(
            *np.array(chosen).T,
            marker,
            ms=7,
            color=color,
            mfc=face,
            ls="none",
            label=f"ours, {label} ({len(chosen)})",
        )
    line = np.array([0.02, 6])
    ax.plot(
        line,
        line * TRUTH[2] / TRUTH[1],
        "--",
        color="#52514e",
        lw=1,
        label="k3/k2 of the generating rates",
    )
    ax.plot(
        TRUTH[1],
        TRUTH[2],
        "*",
        ms=14,
        color="white",
        mec="black",
        label="generating rates",
    )
    ax.plot(
        *TABLE[100][1:],
        "P",
        ms=10,
        color="black",
        label="paper, Table 1, 100 trajectories",
    )
    ax.set(
        xscale="log",
        yscale="log",
        xlim=(0.02, 6),
        ylim=(0.04, 6),
        xlabel="fitted k2 (1/min)",
        ylabel="fitted k3 (1/min)",
    )
    ax.grid(alpha=0.3, lw=0.5)
    ax.legend(fontsize=7, loc="upper left")
    ax.set_title(
        "Fits to regenerated 100-trajectory datasets (upper bound 5)", fontsize=10
    )
    fig.savefig(folder / "C5-recovery.png", dpi=155)
    plt.close(fig)
