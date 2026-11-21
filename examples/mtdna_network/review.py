"""Review runs: interval precision, the H:65 row and the archived variance series.

Each run was declared before it ran (question, size, seed and decision rule);
the rules are repeated in the guide page and in results/review-*.csv.
"""

import argparse
from pathlib import Path

import crn_craft as crn
import numpy as np
import pandas as pd
from study import M, NAMES, ROOT, build, parameters, predicted_slope

CONDITIONS = pd.read_csv(ROOT / "conditions.csv")
CONDITIONS["condition"] = CONDITIONS.panel + ":" + CONDITIONS["index"].astype(str)
EARLY = np.array(M["early_days"], dtype=float)


def heteroplasmy(condition, days, paths, seed):
    row = CONDITIONS.loc[CONDITIONS.condition == condition].iloc[0]
    initial, params = parameters(row)
    result = crn.ensemble(
        build(initial, **params),
        end=float(days[-1]),
        report_times=list(days),
        seed=seed,
        replicates=paths,
        threads=1,
        method="static_direct",
        keep_paths="all",
        strict=True,
        max_output_bytes=1024 * 2**20,
    ).simulate()
    counts = result.paths.sel(variable=list(NAMES)).values.astype(float)
    total = counts.sum(axis=-1)
    if (total == 0).any():
        raise SystemExit(f"{condition}: a path lost all copies")
    h = (counts[..., 1] + counts[..., 3]) / total
    return h, total, predicted_slope(initial, params["mu"])


def resample(h, statistic, replicates, seed):
    """Whole-path bootstrap: each replicate draws paths with replacement."""
    rng = np.random.default_rng(seed)
    out = np.empty(replicates)
    for b in range(replicates):
        out[b] = statistic(h[rng.integers(0, len(h), len(h))])
    return out


def early_ratio(days, prediction):
    # Ensemble variance (divisor n) regressed on time through the origin.
    return lambda h: days @ h.var(axis=0) / (days @ days) / prediction


def primary(scale):
    """R2: the three primary rows with 20,000 bootstrap replicates."""
    rows = []
    for condition in M["primary"]:
        i = int(CONDITIONS.index[CONDITIONS.condition == condition][0])
        paths = 2000 // scale
        h, _, prediction = heteroplasmy(condition, EARLY, paths, M["seed"] + i)
        ratio = early_ratio(EARLY, prediction)
        tail = 0.01 / 6
        lo, hi = np.quantile(resample(h, ratio, 20000 // scale, 1), [tail, 1 - tail])
        small = [
            np.quantile(resample(h, ratio, 600, s), [tail, 1 - tail])
            for s in range(2, 7)
        ]
        rows.append(
            {
                "condition": condition,
                "paths": paths,
                "seed": M["seed"] + i,
                "ratio": ratio(h),
                "ratio_low": lo,
                "ratio_high": hi,
                "replicates": 20000 // scale,
                "low_600_min": min(s[0] for s in small),
                "low_600_max": max(s[0] for s in small),
                "high_600_min": min(s[1] for s in small),
                "high_600_max": max(s[1] for s in small),
                "within_band": bool(lo >= 0.8 and hi <= 1.2),
            }
        )
    return pd.DataFrame(rows)


def h65(scale):
    """R3: H:65 again, same scaled window, fresh seed, 2,000 paths."""
    beta = float(CONDITIONS.loc[CONDITIONS.condition == "H:65", "beta"].iloc[0])
    days = EARLY * min(1.0, 100.0 / beta)
    seed = M["seed"] + 90065
    h, _, prediction = heteroplasmy("H:65", days, 2000 // scale, seed)
    ratio = early_ratio(days, prediction)
    lo, hi = np.quantile(resample(h, ratio, 20000 // scale, 1), [0.005, 0.995])
    return pd.DataFrame(
        [
            {
                "condition": "H:65",
                "paths": len(h),
                "seed": seed,
                "end_day": days[-1],
                "equation13": prediction,
                "ratio": ratio(h),
                "ratio_low": lo,
                "ratio_high": hi,
                "contains_one": bool(lo <= 1 <= hi),
            }
        ]
    )


def archive(scale):
    """R4: fresh runs on the archive's 10-day grid against the archived series."""
    reference = pd.read_csv(ROOT / "reference-variance.csv")
    plans = (("H:55", 500.0, 2000, 290555), ("H:60", 100.0, 400, 290560))
    rows, traces = [], []
    for condition, end, paths, seed in plans:
        if scale > 1:
            end, paths = end / 10, paths // scale
        days = np.arange(0.0, end + 5, 10.0)
        h, total, prediction = heteroplasmy(condition, days, paths, seed)
        panel, index = condition.split(":")
        series = reference[
            (reference.panel == panel) & (reference["index"] == int(index))
        ]
        series = series[series.day <= end + 1e-9]
        if not np.allclose(series.day, days):
            raise SystemExit(f"{condition}: archive grid differs")

        # The authors' estimator: mean of the finite-difference gradient of the
        # sample variance (divisor n - 1) on the 10-day grid.
        def slope(hh):
            return np.mean(np.gradient(hh.var(axis=0, ddof=1), 10.0))

        archived = np.mean(np.gradient(series.variance.to_numpy(), 10.0))
        boot = resample(h, slope, 2000, 1)
        lo, hi = np.quantile(boot, [0.005, 0.995])
        fresh = slope(h)
        fixed = np.mean((h == 0) | (h == 1), axis=0)
        rows.append(
            {
                "condition": condition,
                "paths": paths,
                "seed": seed,
                "end_day": end,
                "equation13": prediction,
                "archived_slope": archived,
                "fresh_slope": fresh,
                "fresh_over_archive": fresh / archived,
                "low_over_archive": lo / archived,
                "high_over_archive": hi / archived,
                "archive_over_equation13": archived / prediction,
                "fresh_over_equation13": fresh / prediction,
                "low_over_equation13": lo / prediction,
                "high_over_equation13": hi / prediction,
                "fixed_fraction_end": fixed[-1],
            }
        )
        traces.append(
            pd.DataFrame(
                {
                    "condition": condition,
                    "day": days,
                    "fresh_variance": h.var(axis=0, ddof=1),
                    "archived_variance": series.variance.to_numpy(),
                    "equation13_variance": prediction * days,
                    "fresh_mean_copies": total.mean(axis=0),
                    "fresh_fixed_fraction": fixed,
                }
            )
        )
    return pd.DataFrame(rows), pd.concat(traces, ignore_index=True)


def plot(output, folder):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    trace = pd.read_csv(output / "review-archive-traces.csv")
    table = pd.read_csv(output / "review-archive.csv").set_index("condition")
    titles = {"H:55": "H:55, network rates x 0.01", "H:60": "H:60, nominal"}
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.6), layout="constrained")
    for ax, (condition, group) in zip(axes, trace.groupby("condition"), strict=True):
        n = int(table.loc[condition, "paths"])
        ax.plot(group.day, group.archived_variance, "o", ms=3, label="archived, 10,000")
        ax.plot(
            group.day,
            group.fresh_variance,
            "s",
            ms=3,
            mfc="none",
            label=f"fresh, {n:,}",
        )
        ax.plot(group.day, group.equation13_variance, "k--", lw=1, label="Equation 13")
        ax.set(
            title=titles[condition],
            xlabel="Time (days)",
            ylabel="Variance of mutant fraction",
        )
        ax.legend(frameon=False, fontsize=8)
    folder.mkdir(parents=True, exist_ok=True)
    fig.savefig(folder / "long-variance.png", dpi=155)
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--size", choices=("quick", "full"), default="quick")
    p.add_argument("--output", type=Path)
    p.add_argument("--figures", type=Path)
    args = p.parse_args()
    if args.figures:
        plot(args.output or ROOT / "results", args.figures)
        return
    scale = 1 if args.size == "full" else 20
    output = args.output or ROOT / f"runs/review-{args.size}"
    output.mkdir(parents=True, exist_ok=True)
    primary(scale).to_csv(output / "review-primary-intervals.csv", index=False)
    h65(scale).to_csv(output / "review-h65.csv", index=False)
    table, trace = archive(scale)
    table.to_csv(output / "review-archive.csv", index=False)
    trace.to_csv(output / "review-archive-traces.csv", index=False)
    print(f"review results written to {output}")


if __name__ == "__main__":
    main()
