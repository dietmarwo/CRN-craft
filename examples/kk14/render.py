"""Generate the figures and tables from saved numerical results."""

from __future__ import annotations

import csv
from collections import defaultdict

import numpy as np

from campaign import save_csv
from reference import HERE, MECHANISMS

NAMES = {
    "drift": "Neutral drift",
    "size": "Faster copying of shorter mtDNA",
    "transcription": "Transcription-linked replication",
}
COLORS = {"drift": "#244d88", "size": "#a74d0d", "transcription": "#166452"}
STYLE = {"drift": "o", "size": "s", "transcription": "^"}


def read_csv(path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def _axes():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.titleweight": "bold",
            "axes.labelcolor": "#253346",
            "text.color": "#182638",
            "axes.edgecolor": "#8797a9",
            "xtick.color": "#253346",
            "ytick.color": "#253346",
            "grid.color": "#e1e7ee",
            "savefig.dpi": 170,
        }
    )
    return plt


def _validation(out):
    return [r for r in read_csv(out / "comparison.csv") if r["cohort"] == "validation"]


def reference_differences(out, dest=None):
    comparison = _validation(out)
    rates = read_csv(out / "rates.csv")
    result = []
    for r in read_csv(HERE / "reference.csv"):
        label = "size-half-life-2d" if r["figure"] == "7" else "primary"
        is_mutation = r["metric"] == "mutation"
        value_column = r["metric"]
        if (
            r["article"] == "PNAS 2014"
            and r["panel"] == "B"
            and r["mechanism"] == "drift"
        ):
            value_column = "all_types"
        # The old drift plot reports all types. Include the alternative
        # count convention as a second difference in the same row.
        low_column = "fit_low" if is_mutation else value_column + "_low"
        high_column = "fit_high" if is_mutation else value_column + "_high"
        case = next(
            (
                c
                for c in (rates if is_mutation else comparison)
                if (c["label"], c["mechanism"], c["years"])
                == (label, r["mechanism"], r["years"])
            ),
            None,
        )
        if case is None or not case.get(value_column):
            continue
        current = float(case[value_column])
        paper = float(r["value"])
        result.append(
            {
                **r,
                "statistic_used": value_column,
                "simulation": current,
                "simulation_low": case.get(low_column),
                "simulation_high": case.get(high_column),
                "difference": current - paper,
                "mutant_only_difference": (
                    float(case["mutant_types"]) - paper
                    if value_column == "all_types"
                    else None
                ),
            }
        )
    save_csv(dest or out / "reference-differences.csv", result)
    return result


def comparison_table(out):
    rows = _validation(out)
    lines = [
        "| Age | Mechanism | Affected / cells | Prevalence (95% Wilson) | Mutant lineages (95% bootstrap) | Target band |",
        "| ---: | --- | ---: | --- | --- | --- |",
    ]
    for row in rows:
        if row["label"] != "primary":
            continue
        prevalence = "—"
        if row.get("prevalence"):
            low, high = (100 * float(row[k]) for k in ("wilson_low", "wilson_high"))
            prevalence = f"{float(row['prevalence']):.2%} ({low:.2f}–{high:.2f}%)"
        clones = "—"
        if row.get("mutant_types"):
            clones = (
                f"{float(row['mutant_types']):.2f} "
                f"({float(row['mutant_types_low']):.2f}–"
                f"{float(row['mutant_types_high']):.2f})"
            )
        lines.append(
            f"| {row['years']} y | {NAMES[row['mechanism']]} | "
            f"{row['cox_negative']} / {row['assigned']} | {prevalence} | "
            f"{clones} | {row['target_band']} |"
        )
    return "\n".join(lines) + "\n"


def tissue_figures(out, images):
    plt = _axes()
    images.mkdir(parents=True, exist_ok=True)
    rows = [r for r in _validation(out) if r["label"] == "primary"]
    if not rows:
        raise ValueError("the results contain no primary validation rows")
    published = [
        r for r in read_csv(HERE / "reference.csv") if r["article"] == "PNAS 2014"
    ]
    ages = sorted({int(r["years"]) for r in rows})
    labels = {(r["mechanism"], r["years"]): r for r in rows}

    fig, (rates, types) = plt.subplots(1, 2, figsize=(10.8, 4.6), layout="constrained")
    for model in MECHANISMS:
        vals = [labels[model, str(age)] for age in ages if (model, str(age)) in labels]
        if not vals:
            continue
        color, marker = COLORS[model], STYLE[model]
        x = [int(r["years"]) for r in vals]
        rates.plot(
            x,
            [float(r["mutation"]) for r in vals],
            color=color,
            marker=marker,
            label=NAMES[model],
        )
        mean = [r for r in vals if r.get("mutant_types")]
        if mean:
            xx = [int(r["years"]) for r in mean]
            m = np.array([float(r["mutant_types"]) for r in mean])
            left = m - [float(r["mutant_types_low"]) for r in mean]
            right = [float(r["mutant_types_high"]) for r in mean] - m
            types.errorbar(
                xx,
                m,
                yerr=np.maximum(0, [left, right]),
                color=color,
                marker=marker,
                capsize=2,
            )
        for axis, panel in ((rates, "A"), (types, "B")):
            paper = [
                r
                for r in published
                if r["panel"] == panel
                and r["mechanism"] == model
                and int(r["years"]) in x
            ]
            axis.scatter(
                [int(r["years"]) for r in paper],
                [float(r["value"]) for r in paper],
                facecolors="none",
                edgecolors=color,
                marker=marker,
                s=65,
                linewidths=1.0,
            )
    rates.set(
        yscale="log",
        title="A. Fit the affected fraction",
        ylabel="Mutation probability per daughter",
    )
    types.set(
        title="B. Predict lineage diversity",
        ylabel="Mutant lineages per affected cell",
        ylim=(0, None),
    )
    for ax in (rates, types):
        ax.set(xticks=ages, xlabel="Virtual lifespan (years)")
        ax.grid(axis="y", linewidth=0.7)
        ax.set_axisbelow(True)
    rates.legend(frameon=False, fontsize=8, loc="upper right")
    cases = sorted({int(r["assigned"]) for r in rows})
    n_label = f"{cases[0]:,}" if len(cases) == 1 else f"{cases[0]:,}–{cases[-1]:,}"
    fig.suptitle(
        f"mtDNA mechanisms at a 10% COX-negative target · {n_label} held-out cells per case"
    )
    types.text(
        0.98,
        0.96,
        "Open symbols show published values",
        fontsize=8,
        ha="right",
        va="top",
        transform=types.transAxes,
    )
    fig.savefig(images / "mechanisms.png")
    plt.close(fig)

    histogram = out / "diversity-3y.csv"
    if histogram.is_file():
        grouped = defaultdict(list)
        for r in read_csv(histogram):
            if r["label"] == "primary":
                grouped[r["mechanism"]].append((int(r["types"]), int(r["cells"])))
        if grouped:
            fig, axes = plt.subplots(
                3, 1, figsize=(8, 6.2), sharex=True, layout="constrained"
            )
            for name, ax in zip(MECHANISMS, axes, strict=True):
                values = sorted(grouped[name])
                n = sum(count for _, count in values)
                if n:
                    ax.bar(
                        [k for k, _ in values],
                        [count / n for _, count in values],
                        color=COLORS[name],
                        width=0.83,
                    )
                ax.set_title(f"{NAMES[name]} · {n:,} affected", loc="left", fontsize=10)
                ax.grid(axis="y", linewidth=0.7)
                ax.set_axisbelow(True)
            axes[-1].set_xlabel("Distinct mutant lineages in one affected cell")
            axes[1].set_ylabel("Fraction of affected cells")
            fig.suptitle("Lineage diversity at 3 years")
            fig.savefig(images / "diversity-3y.png")
            plt.close(fig)


def single_cell_figure(atp, delayed, images):
    plt = _axes()
    images.mkdir(parents=True, exist_ok=True)
    fig, (left, right) = plt.subplots(1, 2, figsize=(10.5, 4.05), layout="constrained")
    t, w, m, a = atp.T
    left.plot(t / 24, w, color="#244d88", label="Wild type")
    left.plot(t / 24, m, color="#a74d0d", label="Mutants")
    left.plot(t[a >= 1] / 24, a[a >= 1], color="#166452", label="ATP")
    left.axvline(t[-1] / 24, color="#8797a9", lw=1, ls=":")
    left.set(
        title="ATP-limited competition",
        xlabel="Days",
        ylabel="Copies or ATP units",
        yscale="log",
    )
    h, w, bw, m, bm = delayed.T
    right.plot(h / 8760, w + bw, color="#244d88", label="Wild type")
    right.plot(h / 8760, m + bm, color="#a74d0d", label="Mutants")
    right.set(
        title="Fixed-delay competition",
        xlabel="Years (365 days)",
        ylabel="mtDNA copies",
        yscale="log",
    )
    for axis in (left, right):
        axis.grid(axis="y", linewidth=0.7)
        axis.legend(frameon=False, fontsize=8)
    fig.savefig(images / "single-cell.png")
    plt.close(fig)
