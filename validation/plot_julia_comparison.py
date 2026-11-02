"""Render the measured Julia comparison without rounding the underlying data."""

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def main():
    root = Path(__file__).resolve().parents[1]
    report = json.loads(
        (root / "validation/campaign/julia-comparison.json").read_text()
    )
    if not report["passed"] or report["smoke"]:
        raise RuntimeError("A passing full comparison is required")
    engines = [
        "CRN-cached",
        "CRN-direct",
        "Direct",
        "SortingDirect",
        "DirectCR",
        "NRM",
        "RSSACR",
    ]
    labels = [
        "CRN cached Direct",
        "CRN Direct",
        "Julia Direct",
        "Julia SortingDirect",
        "Julia DirectCR",
        "Julia NRM",
        "Julia RSSACR",
    ]
    titles = [
        "W1 · Immigration–death · 2 reactions",
        "W2 · Sparse · 66 reactions",
        "W3 · Repressilator · 12 reactions",
        "W4 · Dense · 66 reactions",
    ]
    plt.rcParams.update({"font.size": 10, "svg.fonttype": "none"})
    fig, axes = plt.subplots(2, 2, figsize=(13, 8.5))
    for ax, name, title in zip(axes.flat, report["summary"], titles, strict=True):
        rows = report["summary"][name]
        medians = np.array([rows[e]["median_seconds"] for e in engines])
        low = np.array([rows[e]["min_seconds"] for e in engines])
        high = np.array([rows[e]["max_seconds"] for e in engines])
        y = np.arange(len(engines))
        ax.barh(y, medians, color=["#183b66", "#8199b3", *["#19857d"] * 5])
        ax.errorbar(
            medians,
            y,
            xerr=np.array([medians - low, high - medians]),
            fmt="none",
            ecolor="#243746",
            elinewidth=1,
            capsize=2.5,
        )
        for i, value in enumerate(medians):
            ax.text(high[i] + high.max() * 0.025, i, f"{value:.3f}", va="center")
        ax.set_yticks(y, labels)
        ax.invert_yaxis()
        ax.set_xlim(0, high.max() * 1.23)
        ax.set_xlabel("Seconds for 64 paths")
        ax.set_title(title, loc="left", fontsize=11, fontweight="bold", pad=13)
        ax.grid(axis="x", color="#e1e6eb", linewidth=0.6)
        ax.set_axisbelow(True)
        ax.tick_params(axis="both", length=0)
        for spine in ax.spines.values():
            spine.set_visible(False)
    fig.suptitle(
        "Exact SSA: warm ensemble time", x=0.04, ha="left", fontsize=20, weight="bold"
    )
    fig.text(
        0.04,
        0.925,
        "t = 100 · one numerical worker · median and range of five runs · lower is faster",
        color="#46596b",
    )
    fig.text(
        0.04,
        0.025,
        "2 Nov 2026 · Julia 1.13.3 / Catalyst 17.0.0 / JumpProcesses 10.0.1"
        " · shared AMD EPYC host · panel scales differ",
        fontsize=9,
        color="#46596b",
    )
    fig.tight_layout(rect=(0.025, 0.055, 0.99, 0.90), h_pad=3.0, w_pad=3.0)
    for suffix in ["png", "svg"]:
        path = root / f"docs/images/julia-comparison-warm.{suffix}"
        fig.savefig(path, dpi=180)
        if suffix == "svg":
            path.write_text(
                "\n".join(line.rstrip() for line in path.read_text().splitlines())
                + "\n"
            )


if __name__ == "__main__":
    main()
