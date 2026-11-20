"""Plot results using the tables produced by this example."""

from __future__ import annotations

import numpy as np

from fitting import read_csv
from model import THRESHOLDS

COLORS = ("#234d82", "#148068", "#a04716", "#6e5294")
NAMES = {
    "mutation": "Mutation probability `p`",
    "advantage": "Selection advantage `s`",
    "advantageous": "Advantageous fraction `f`",
}


def _number(name, value):
    value = float(value)
    if name == "mutation":
        mantissa, exponent = f"{value:.2e}".split("e")
        return f"{mantissa}e-{int(exponent[1:])}"
    return f"{value:.3f}"


def recovery_table(source):
    """The chapter's table of reconstructed and published recovery results."""
    lines = [
        "| Regime | Parameter | Truth | Search box | Mean (CV) | 2.5–97.5% | Published mean (CV) | Published 2.5–97.5% | Truth inside: here / paper |",
        "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |",
    ]
    for row in read_csv(source / "parameter-differences.csv"):
        name = row["parameter"]

        def n(key, row=row, name=name):
            return _number(name, row[key])

        inside = (
            float(row["paper_p2_5"]) <= float(row["truth"]) <= float(row["paper_p97_5"])
        )
        here = row["truth_in_interval"] == "True"
        lines.append(
            f"| {row['regime'].capitalize()} | {NAMES[name]} | {n('truth')} | "
            f"{n('box_low')}–{n('box_high')} | "
            f"{n('mean')} ({float(row['cv_percent']):.1f}%) | "
            f"{n('p2_5')}–{n('p97_5')} | "
            f"{n('paper_mean')} ({float(row['paper_cv_percent']):.1f}%) | "
            f"{n('paper_p2_5')}–{n('paper_p97_5')} | "
            f"{'yes' if here else 'no'} / {'yes' if inside else 'no'} |"
        )
    return "\n".join(lines) + "\n"


def _plt():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {
            "font.size": 10,
            "font.family": "DejaVu Sans",
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.titleweight": "bold",
            "axes.edgecolor": "#8693a1",
            "grid.color": "#e2e8ef",
            "savefig.dpi": 170,
        }
    )
    return plt


def figures(source, images):
    plt = _plt()
    images.mkdir(parents=True, exist_ok=True)
    data = read_csv(source / "fig6-curves.csv")
    paper = read_csv(source / "fig6-differences.csv")
    fig, panels = plt.subplots(
        1, 2, figsize=(10.7, 4.3), layout="constrained", sharey=True
    )
    for ax, panel, backend, title in zip(
        panels,
        ("A", "B"),
        ("ATP 1 h", "discrete"),
        ("A. Hourly ATP model", "B. Fixed-population Moran model"),
        strict=True,
    ):
        ax.set(title=title, xlabel="Age (days)", xlim=(0, 1080), ylim=(0, 6.0))
        ax.grid(axis="y", linewidth=0.7)
        for color, k in zip(COLORS, THRESHOLDS, strict=True):
            rows = [
                r
                for r in data
                if r["backend"] == backend
                and r["case"] == "baseline"
                and int(r["threshold_percent"]) == k
            ]
            refs = [
                r
                for r in paper
                if r["panel"] == panel
                and int(r["threshold_percent"]) == k
                and int(r["age_days"]) % 90 == 0
            ]
            ax.plot(
                [float(r["age_days"]) for r in rows],
                [100 * float(r["fraction"]) for r in rows],
                color=color,
                label=f">{k}%",
            )
            ax.scatter(
                [int(r["age_days"]) for r in refs],
                [100 * float(r["value"]) for r in refs],
                color=color,
                facecolors="white",
                edgecolors=color,
                s=18,
                lw=0.8,
                zorder=3,
            )
    panels[0].set_ylabel("Cells exceeding the mutant threshold (%)")
    panels[1].legend(title="Mutant load", frameon=False, loc="upper left", fontsize=9)
    fig.suptitle(
        "Cross-sectional mtDNA snapshots · lines reproduce the models, circles read the paper"
    )
    fig.savefig(images / "fig6-snapshots.png")
    plt.close(fig)

    curve = read_csv(source / "fig5-passage.csv")
    fig, (ax, censor) = plt.subplots(1, 2, figsize=(10.7, 4.25), layout="constrained")
    labels = list(
        dict.fromkeys(r["label"] for r in curve if r["label"].startswith("s="))
    )
    for color, label in zip(COLORS, labels, strict=True):
        rows = [r for r in curve if r["label"] == label]
        days = np.array([float(r["elapsed_days"]) for r in rows])
        cdf = np.array([float(r["conditional_cdf"]) for r in rows])
        mass = np.diff(cdf)
        mids = (days[1:] + days[:-1]) / 2
        ax.plot(mids, mass, label=label, color=color)
        if label == "s=0.1":
            i = np.searchsorted(days, 1080)
            censor.plot(days, cdf, color=color, label="All successful paths")
            censor.plot(
                days[: i + 1],
                cdf[: i + 1] / cdf[i],
                color="#a04716",
                ls="--",
                label="Takeover observed by day 1,080",
            )
            censor.axvline(1080, color="#8693a1", lw=1, ls=":")
    ax.set(
        title="A. 30-day bins conditional on takeover",
        xlabel="Days since one founder",
        ylabel="Probability in a 30-day bin",
        xlim=(120, 1850),
        ylim=(0, None),
    )
    censor.set(
        title="B. Bias from a short observation window",
        xlabel="Days since one founder",
        ylabel="Conditional cumulative probability",
        xlim=(400, 2000),
        ylim=(0, 1.03),
    )
    for panel in (ax, censor):
        panel.grid(axis="y", linewidth=0.7)
        panel.legend(frameon=False, fontsize=8)
    fig.suptitle("Time to 90% mutant abundance in the Moran chain")
    fig.savefig(images / "fig5-takeover.png")
    plt.close(fig)
