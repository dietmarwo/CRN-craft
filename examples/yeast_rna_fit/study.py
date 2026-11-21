"""Reproduce the mean versus histogram fits and predict independent RNA measurements."""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from fitting import (
    Data,
    exact_moments,
    fit_nonspatial,
    nonspatial_laws,
    observations,
    refine_nonspatial,
    spatial_grid,
    spatial_laws,
)
from model import KEYS, ROOT, parameters
from sites import prediction

M = json.loads((ROOT / "manifest.json").read_text())
ANCHORS = json.loads((ROOT / "source-parameters.json").read_text())["parameters"]


def curve(x, label, data, cap, step):
    p = nonspatial_laws(x, data, cap, step)
    result, profiles = [], []
    for dose, law in p.items():
        for minute in sorted(
            set([0, *[t for d, _, t, *_ in data.groups if d == dose]])
        ):
            row = law[minute]
            result.extend(
                (label, dose, minute, k, value) for k, value in enumerate(row[:201])
            )
            result.append((label, dose, minute, 201, 1 - row[:201].sum()))
            xgrid = np.arange(len(row))
            mean = np.dot(xgrid, row)
            sd = np.sqrt(max(0.0, np.dot((xgrid - mean) ** 2, row)))
            profiles.append((label, dose, minute, mean, sd, float(row[4:].sum())))
    return result, profiles


def fitting_run(args):
    preset = M["presets"][args.size].copy()
    if args.budget:
        preset["evaluations"] = args.budget
    frame = observations()
    if args.size == "quick":
        frame = frame[(frame.dose == 0.2) & frame.minute.isin([0, 8, 15, 25, 45, 55])]
    anchor = parameters(ANCHORS["NonSpatial-Distributions"])
    data = Data(frame)
    table, totals, rows, profiles, report = [], [], [], [], {}
    dose_labels = ["all"]
    if args.cross_validate:
        dose_labels.extend([str(d) for d in data.dose])
    for k, part in enumerate(dose_labels):
        mask = np.ones(len(frame), dtype=bool)
        if part != "all":
            held = M["holdouts"][part]
            mask = frame.replicate != held
        train = Data(frame[mask])
        test = Data(frame[~mask])
        for kind in ("mean", "histogram"):
            path = (
                args.output
                / f"fit-{args.size}-{part}-{kind}-{preset['evaluations']}-{int(args.budget_check)}.json"
            )
            if path.exists():
                note = json.loads(path.read_text())
                x = parameters(note["parameters"])
                restarts = note["restarts"]
            else:
                x, _, restarts, _ = fit_nonspatial(
                    train,
                    anchor,
                    preset,
                    M["seed"] + k * 1000,
                    workers=args.workers,
                    kind=kind,
                    extra=args.budget_check,
                )
                note = {
                    "parameters": dict(zip(KEYS, x, strict=True)),
                    "restarts": restarts,
                }
                temp = path.with_suffix(".tmp")
                temp.write_text(json.dumps(note, allow_nan=False, indent=2) + "\n")
                temp.replace(path)
            table.extend({"partition": part, "objective": kind, **r} for r in restarts)
            final_cap = 360
            high = nonspatial_laws(x, train, final_cap, 1.5)
            while (
                max(max(1 - v.sum(-1)) for v in high.values()) > 1e-5
                and final_cap < 2880
            ):
                final_cap *= 2
                high = nonspatial_laws(x, train, final_cap, 1.5)
            low = nonspatial_laws(x, train, max(360, final_cap // 2), 3.0)
            loss, low_loss = train.counts_cost(high), train.counts_cost(low)
            if not np.isfinite(loss):
                raise RuntimeError(
                    f"Refined {part}/{kind} fit has impossible occupied bins"
                )
            totals.append(
                {
                    "partition": part,
                    "objective": kind,
                    "training_cells": train.n(),
                    "mean_cost": refine_nonspatial(x, train, "mean", final_cap, 1.5),
                    "final_cap": final_cap,
                    "histogram_nll_per_cell": loss / train.n(),
                    "coarse_refined_nats_per_cell": abs(low_loss - loss) / train.n(),
                    "max_physical_leakage_estimate": max(
                        max(1 - v.sum(-1)) for v in high.values()
                    ),
                }
            )
            if part == "all":
                mh, moments = curve(x, kind, data, final_cap, 1.5)
                rows += mh
                profiles += moments
                report[kind] = note["parameters"]
            else:
                check = test.counts_cost(
                    nonspatial_laws(x, test, final_cap, 1.5), per_rep=True
                )
                if isinstance(check, float):
                    raise RuntimeError("Holdout has impossible occupied bins")
                for (dose, rep), (score, n) in check.items():
                    totals[-1].update(
                        held_dose=dose,
                        held_rep=rep,
                        held_n=n,
                        held_nll_per_cell=score / n,
                    )
            print(part, kind, totals[-1], flush=True)
        pd.DataFrame(totals).to_csv(args.output / "scores.csv", index=False)
        pd.DataFrame(table).to_csv(args.output / "restarts.csv", index=False)
    for flavor in ("Means", "Distributions"):
        x = parameters(ANCHORS[f"NonSpatial-{flavor}"])
        if flavor == "Distributions":
            mh, moments = curve(x, f"source-{flavor.lower()}", data, 400, 1.5)
            rows += mh
            profiles += moments
        # Exact moment equations handle the high-turnover mean-only anchor.
        means = exact_moments(x, data)
        for dose, series in means.items():
            for minute, value in enumerate(series[:, 0]):
                rows.append(
                    (f"source-{flavor.lower()}-moment", dose, minute, -1, value)
                )
    pd.DataFrame(
        rows, columns=["fit", "dose", "minute", "count", "probability"]
    ).to_csv(args.output / "histograms.csv", index=False)
    pd.DataFrame(
        profiles, columns=["fit", "dose", "minute", "mean", "sd", "on_fraction"]
    ).to_csv(args.output / "model-moments.csv", index=False)
    (args.output / "parameters.json").write_text(json.dumps(report, indent=2) + "\n")


def spatial_run(args):
    """Evaluate the source spatial model at its mature-count fit parameters."""
    data = Data(observations(), spatial_model=True)
    x = parameters(ANCHORS["Spatial-Distributions"])
    if args.spatial_fit:
        preset = M["presets"][args.size]
        budget = args.budget or preset["evaluations"]
        budgets = [budget] * preset["restarts"]
        if args.budget_check:
            budgets.append(4 * budget)
        fits, fitted_law = spatial_grid(
            data, x, args.output, min(args.workers, 4), budgets, M["seed"]
        )
        (args.output / "spatial-parameters.json").write_text(
            json.dumps(
                {mode: dict(zip(KEYS, p, strict=True)) for mode, p in fits.items()},
                indent=2,
            )
            + "\n"
        )
        x = fits["histogram"]
        pp = fitted_law
    else:
        pp = spatial_laws(x, data, 28, 305, 4e-7)
    rows = []
    for dose, predicted in pp.items():
        for minute, grid in predicted.items():
            rows.append((dose, minute, 1 - grid.sum(), float(grid.min())))
    score = data.counts_cost(pp)
    reference = pd.read_csv(ROOT / "reference-site-activity.csv")
    predicted, tail = prediction(x, reference)
    activity = predicted[predicted.intensity_bin == 0][
        ["dose", "minute", "model_probability"]
    ].copy()
    activity["fresh_active_fraction"] = 1 - activity.pop("model_probability")
    reference.merge(activity, on=["dose", "minute"], how="left").to_csv(
        args.output / "site-activity-comparison.csv", index=False
    )
    predicted.to_csv(args.output / "sites.csv", index=False)
    pd.DataFrame(
        rows, columns=["dose", "minute", "leakage_estimate", "minimum_mass"]
    ).to_csv(args.output / "spatial-accounting.csv", index=False)
    meta = {
        "parameters": dict(zip(KEYS, x, strict=True)),
        "mature_nll_per_cell": score / data.n(),
        "site_max_leakage": tail,
        "total_usable_cells": data.n(),
    }
    (args.output / "spatial.json").write_text(json.dumps(meta, indent=2) + "\n")
    print(json.dumps(meta), flush=True)


def plot(output, folder):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    model = pd.read_csv(output / "histograms.csv")
    obs = observations()
    folder.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(
        2, 4, figsize=(12, 5.4), sharex=True, sharey=True, layout="constrained"
    )
    for i, dose in enumerate((0.2, 0.4)):
        for ax, minute in zip(axes[i], (8, 15, 30, 45), strict=True):
            data = (
                obs[(obs.dose == dose) & (obs.minute == minute)]
                .groupby("total")
                .cells.sum()
            )
            ax.bar(
                data.index,
                data / data.sum(),
                color="#cbd8e1",
                width=1,
                label="Measured",
            )
            for fit, color in (("mean", "#b14731"), ("histogram", "#234a74")):
                pred = model[
                    (model.dose == dose)
                    & (model.minute == minute)
                    & (model.fit == fit)
                    & (model["count"] >= 0)
                ]
                ax.plot(
                    pred["count"], pred.probability, label=fit.capitalize(), color=color
                )
            ax.set(xlim=(-1, 125), ylim=(0, 0.24), title=f"{dose} M, {minute} min")
            if i == 1:
                ax.set_xlabel("Total mature RNA")
        axes[i, 0].set_ylabel("Probability")
    axes[0, 0].legend(frameon=False, fontsize=8)
    fig.savefig(folder / "rna-fits.png", dpi=155)
    plt.close(fig)
    summary = pd.read_csv(output / "model-moments.csv")
    fig, axes = plt.subplots(2, 3, figsize=(10, 5.7), layout="constrained")
    for i, dose in enumerate((0.2, 0.4)):
        relevant = obs[obs.dose == dose]
        for (_rep, minute), rows in relevant.groupby(["replicate", "minute"]):
            p = rows.cells / rows.cells.sum()
            mean = p @ rows.total
            sd = np.sqrt(p @ (rows.total - mean) ** 2)
            for ax, value in zip(
                axes[i], (mean, sd, p[rows.total > 3].sum()), strict=True
            ):
                ax.scatter(minute, value, color=".5", alpha=0.65, s=10)
        for ax, (column, title) in zip(
            axes[i],
            (
                ("mean", "Mean count"),
                ("sd", "Standard deviation"),
                ("on_fraction", "Fraction above three"),
            ),
            strict=True,
        ):
            for kind, color in (("mean", "#b14731"), ("histogram", "#234a74")):
                model_values = summary[(summary.dose == dose) & (summary.fit == kind)]
                ax.plot(
                    model_values.minute,
                    model_values[column],
                    color=color,
                    label=kind.capitalize(),
                )
            ax.set(title=f"{dose} M · {title}", xlabel="Time (minutes)")
    axes[0, 0].legend(frameon=False, fontsize=8)
    fig.savefig(folder / "rna-summaries.png", dpi=155)
    plt.close(fig)
    if not (output / "sites.csv").exists():
        return
    data = pd.read_csv(ROOT / "reference-site-activity.csv")
    dist = pd.read_csv(ROOT / "reference-site-distribution.csv")
    predictions = pd.read_csv(output / "sites.csv")
    fig, axes = plt.subplots(2, 2, figsize=(9.2, 6.4), layout="constrained")
    for i, dose in enumerate((0.2, 0.4)):
        p = predictions[predictions.dose == dose]
        ts = data[data.dose == dose]
        for name, subset in ts.groupby("series"):
            if name == "published-model":
                axes[i, 0].plot(
                    subset.minute,
                    subset.active_fraction,
                    color=".50",
                    linestyle="--",
                    label="Authors’ model",
                )
            else:
                axes[i, 0].scatter(
                    subset.minute,
                    subset.active_fraction,
                    s=13,
                    alpha=0.65,
                    label="Recorded replicates"
                    if name == sorted(ts.series.unique())[0]
                    else None,
                )
        line = p[p.intensity_bin == 0]
        axes[i, 0].plot(
            line.minute,
            1 - line.model_probability,
            color="#234a74",
            label="Frozen prediction",
        )
        axes[i, 0].set(
            xlabel="Time (minutes)", ylabel="Fraction active", title=f"{dose} M NaCl"
        )

        part = dist[(dist.dose == dose) & (dist.upper_intensity >= 3)]
        for series, u in part.groupby("series"):
            color, dash = (
                (".5", "--") if series == "published-model" else ("#7a9daa", "-")
            )
            axes[i, 1].step(
                u.upper_intensity - 0.5,
                u.probability,
                color=color,
                linestyle=dash,
                alpha=0.7,
                where="mid",
            )
        weights = obs[obs.dose == dose].groupby("minute").cells.sum()
        law = p[p.intensity_bin > 0].copy()
        law["weight"] = law.minute.map(weights)
        law["weighted"] = law.weight * law.model_probability
        averaged = law.groupby("intensity_bin").weighted.sum() / weights.sum()
        positive = averaged > 1e-6
        axes[i, 1].plot(
            averaged.index[positive] + 0.5,
            averaged[positive],
            color="#234a74",
            label="Frozen prediction",
        )
        axes[i, 1].set(
            xlabel="Mature-RNA intensity equivalents",
            ylabel="Probability per cell",
            title="Pooled sampling times",
            xlim=(1.4, 20),
            yscale="log",
            ylim=(1e-5, 0.2),
        )
    axes[0, 0].legend(frameon=False, fontsize=8)
    axes[0, 1].legend(frameon=False, fontsize=8)
    fig.savefig(folder / "transcription-sites.png", dpi=155)
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--size", choices=("quick", "draft", "full"), default="quick")
    p.add_argument("--workers", type=int, default=1)
    p.add_argument("--output", type=Path)
    p.add_argument("--figures", type=Path)
    p.add_argument("--spatial", action="store_true")
    p.add_argument("--spatial-only", action="store_true")
    p.add_argument(
        "--spatial-fit",
        action="store_true",
        help="Refit initiation and cytoplasmic decay on a verified sparse grid",
    )
    p.add_argument("--cross-validate", action="store_true")
    p.add_argument("--budget-check", action="store_true")
    p.add_argument("--budget", type=int)
    a = p.parse_args()
    a.output = a.output or ROOT / ("results" if a.figures else f"runs/{a.size}")
    if a.figures:
        plot(a.output, a.figures)
        return
    a.output.mkdir(parents=True, exist_ok=True)
    if not a.spatial_only:
        fitting_run(a)
    if a.spatial or a.spatial_only or a.spatial_fit or a.size == "full":
        spatial_run(a)


if __name__ == "__main__":
    main()
