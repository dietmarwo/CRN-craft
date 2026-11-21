"""Markdown tables of the guide page, rendered from the stored result files.

``validation/check_two_stage_expression.py`` compares the page with these.
"""

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
LN2 = np.log(2)


def _row(cells):
    return "| " + " | ".join(cells) + " |"


def same_mean_table(results=RESULTS):
    hist = pd.read_csv(results / "histograms.csv")
    div = pd.read_csv(results / "divergence.csv")
    lines = [
        _row(
            (
                "a, b",
                "γ",
                "Mean mRNA a/γ",
                "Most likely count",
                "P(0)",
                "Mass above 150: reaction law / NB / sample",
                "Sample CDF gap",
                "DKW bound",
            )
        ),
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for (a, b, gamma), s in hist.groupby(["a", "b", "gamma"], sort=False):
        s = s.sort_values("count")
        body = s[s["count"] <= 150]
        tail = s[s["count"] == 151].iloc[0]
        d = div[(div.a == a) & np.isclose(div.gamma, gamma)].iloc[0]
        lines.append(
            _row(
                (
                    f"{a:g}, {b:g}",
                    f"{gamma:g}",
                    f"{a / gamma:g}",
                    f"{int(body['count'].iloc[body.reaction_pmf.argmax()])}",
                    f"{body.reaction_pmf.iloc[0]:.4f}",
                    f"{tail.reaction_pmf:.4f} / {tail.negative_binomial:.4f}"
                    f" / {tail.ssa:.4f}",
                    f"{d.cdf_gap:.4f}",
                    f"{d.dkw_bound:.4f}",
                )
            )
        )
    return "\n".join(lines) + "\n"


def divergence_table(results=RESULTS):
    review = pd.read_csv(results / "review-figure1d-conventions.csv")
    lines = [
        _row(
            (
                "γ",
                "a=20: published",
                "a=20: bits",
                "a=20: bits, 10,000 draws",
                "a=0.5: published",
                "a=0.5: bits",
                "a=0.5: bits, 10,000 draws",
            )
        ),
        "| ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    narrow = review[review.a == 20].reset_index(drop=True)
    broad = review[review.a == 0.5].reset_index(drop=True)
    for n, w in zip(narrow.itertuples(), broad.itertuples(), strict=True):
        cells = [f"{n.gamma:.3g}"]
        for r in (n, w):
            cells += [
                f"{r.published:.3f}",
                f"{r.bits_reaction_vs_nb:.3f}",
                f"{r.bits_sampled_10000_mean:.3f} ± {r.bits_sampled_10000_sd:.3f}",
            ]
        lines.append(_row(cells))
    return "\n".join(lines) + "\n"


CANDIDATES = (
    ("nats_reaction_vs_nb", "Reaction law vs NB, nats, 0–150 + tail (shipped)"),
    ("nats_nb_vs_reaction", "Reversed direction, nats"),
    ("nats_no_tail_bin", "No tail bin, nats"),
    ("nats_renormalized_0_150", "0–150 renormalized, nats"),
    ("nats_full_support", "All counts, nats"),
    ("nats_tail_above_100", "Tail bin above 100, nats"),
    ("nats_tail_above_300", "Tail bin above 300, nats"),
    ("nats_width5_bins", "Bins of width 5, nats"),
    ("bits_reaction_vs_nb", "Shipped convention in bits"),
    ("bits_sampled_1000_mean", "Bits, 1,000 simulated draws"),
    ("bits_sampled_10000_mean", "Bits, 10,000 simulated draws"),
)


def convention_table(results=RESULTS):
    review = pd.read_csv(results / "review-figure1d-conventions.csv")
    large = review.published >= 0.1
    lines = [
        _row(
            (
                "Convention",
                "Published / candidate, 8 values ≥ 0.1",
                "Largest gap, 18 values",
            )
        ),
        "| --- | ---: | ---: |",
    ]
    for column, label in CANDIDATES:
        ratio = review.published[large] / review[column][large]
        gap = (review.published - review[column]).abs().max()
        lines.append(
            _row((label, f"{ratio.min():.2f}–{ratio.max():.2f}", f"{gap:.3f}"))
        )
    return "\n".join(lines) + "\n"


TABLES = {
    "two-stage-same-mean": same_mean_table,
    "two-stage-divergence": divergence_table,
    "two-stage-conventions": convention_table,
}
