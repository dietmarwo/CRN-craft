"""Reproduce the forward figures and compare the results with the paper."""

from __future__ import annotations

import numpy as np

from fitting import CAPS, STEP, TOLERANCE, axes_for, from_x, read_csv, save_csv
from model import (
    AGES,
    BASE,
    HERE,
    REGIMES,
    THRESHOLDS,
    Model,
    cohort,
    forward,
    passage,
    tails,
)

PAPER_MODES = {0.10: 990, 0.20: 540, 0.30: 390, 0.50: 240}
ATP_CASES = ((1.0, 100_000, 712648), (0.5, 30_000, 821937))


def forward_evidence(out, arrays, atp):
    """Write comparable series from computed arrays and the digitized source."""
    report, checks = [], []
    for case in REGIMES:
        exact = arrays[f"{case}-discrete"]
        fast = arrays[f"{case}-split"]
        finer = arrays[f"{case}-fine"]
        expected = tails(exact[:, :5])
        for backend, data in (("split", fast), ("finer", finer)):
            difference = tails(data[:, :5]) - expected
            gap = float(np.max(np.abs(difference)))
            checks.append(
                {
                    "regime": case,
                    "comparison": backend + "_vs_discrete",
                    "max_curve_difference": gap,
                    "max_leaked_mass": float(data[:, 5].max()),
                    "max_negative_mass": float(data[:, 6].max()),
                    "reference_leakage": float(exact[:, 5].max()),
                    "accepted": gap <= TOLERANCE,
                }
            )
        for backend, data in (("discrete", exact), ("split", fast), ("finer", finer)):
            for age, values, row in zip(AGES, tails(data[:, :5]), data, strict=True):
                for threshold, value in zip(THRESHOLDS, values, strict=True):
                    report.append(
                        {
                            "case": case,
                            "backend": backend,
                            "age_days": int(age),
                            "threshold_percent": threshold,
                            "fraction": value,
                            "leakage": row[5],
                            "negative_mass": row[6],
                        }
                    )
    for step, count, seed in ATP_CASES:
        bins, stops = atp[step]
        for age, values, stopped in zip(AGES, tails(bins) / count, stops, strict=True):
            for threshold, value in zip(THRESHOLDS, values, strict=True):
                report.append(
                    {
                        "case": "baseline",
                        "backend": f"ATP {step:g} h",
                        "age_days": int(age),
                        "threshold_percent": threshold,
                        "fraction": value,
                        "cells": count,
                        "stops": stopped,
                        "seed": seed,
                    }
                )

    lookup = {
        (r["case"], r["backend"], r["age_days"], r["threshold_percent"]): r["fraction"]
        for r in report
    }
    pairs = []
    for row in read_csv(HERE / "paper-figure6.csv"):
        name = "ATP 1 h" if row["panel"] == "A" else "discrete"
        key = ("baseline", name, int(row["age_days"]), int(row["threshold_percent"]))
        pred = lookup[key]
        pairs.append(
            {
                **row,
                "model": name,
                "prediction": pred,
                "signed_difference": pred - float(row["value"]),
            }
        )
    for panel in ("A", "B"):
        d = np.array([r["signed_difference"] for r in pairs if r["panel"] == panel])
        checks.append(
            {
                "regime": "baseline",
                "comparison": f"Figure 6 {panel}",
                "max_curve_difference": float(np.max(np.abs(d))),
                "rmse": float(np.sqrt(np.mean(d * d))),
                "ssq": float(d @ d),
            }
        )
    (n1, n2) = (ATP_CASES[0][1], ATP_CASES[1][1])
    hourly, half = tails(atp[1.0][0]) / n1, tails(atp[0.5][0]) / n2
    d = hourly - half
    # Two-proportion z with the pooled proportion, where it is not 0 or 1.
    pooled = (hourly * n1 + half * n2) / (n1 + n2)
    varied = (pooled > 0) & (pooled < 1)
    se = np.sqrt(pooled[varied] * (1 - pooled[varied]) * (1 / n1 + 1 / n2))
    z = d[varied] / se
    checks.append(
        {
            "regime": "baseline",
            "comparison": "ATP mesh and Monte Carlo",
            "max_curve_difference": float(np.max(np.abs(d))),
            "rmse": float(np.sqrt(np.mean(d * d))),
            "ssq": float(np.sum(d * d)),
            "compared_points": int(varied.sum()),
            "max_abs_z": float(np.max(np.abs(z))),
        }
    )
    for name in REGIMES:
        for corner, values in arrays.items():
            if not corner.startswith(f"{name}-corner-"):
                continue
            split = values[1]
            difference = tails(split[:, :5]) - tails(values[0][:, :5])
            gap = float(np.max(np.abs(difference)))
            checks.append(
                {
                    "regime": name,
                    "comparison": corner.removeprefix(f"{name}-")
                    + " split_vs_discrete",
                    "max_curve_difference": gap,
                    "max_leaked_mass": float(split[:, 5].max()),
                    "max_negative_mass": float(split[:, 6].max()),
                    "reference_leakage": float(values[0][:, 5].max()),
                    "accepted": gap <= TOLERANCE,
                }
            )
    if "baseline-discrete-complete" in arrays:
        data = arrays["baseline-discrete-complete"]
        d = tails(data[:, :5]) - tails(arrays["baseline-discrete"][:, :5])
        checks.append(
            {
                "regime": "baseline",
                "comparison": "neutral cap vs complete discrete",
                "max_curve_difference": float(np.max(np.abs(d))),
                "max_leaked_mass": float(data[:, 5].max()),
                "max_negative_mass": float(data[:, 6].max()),
                "accepted": bool(np.max(np.abs(d)) <= TOLERANCE),
            }
        )
    save_csv(out / "fig6-curves.csv", report)
    save_csv(out / "fig6-differences.csv", pairs)
    save_csv(out / "forward-checks.csv", checks)
    return checks


def compute_forward(out, threads, complete=False):
    arrays = {}
    for name, model in REGIMES.items():
        cap = CAPS[name]
        print(f"{name}: discrete and two implicit resolutions", flush=True)
        arrays[f"{name}-discrete"] = forward(model, cap=cap)
        arrays[f"{name}-split"] = forward(model, cap=cap, step=STEP)
        arrays[f"{name}-fine"] = forward(model, cap=cap, step=STEP / 2)
        # Two corners of the search box with the largest p and s, at both
        # ends of the advantageous fraction (the largest neutral load).
        axes = axes_for(name)
        for label, f in (("high-f", axes[2][-1]), ("low-f", axes[2][0])):
            corner = from_x(model, [axes[0][-1], axes[1][-1], f])
            print(f"{name}: corner {label}, discrete and implicit", flush=True)
            arrays[f"{name}-corner-{label}"] = (
                forward(corner, cap=cap),
                forward(corner, cap=cap, step=STEP),
            )
    if complete:
        arrays["baseline-discrete-complete"] = forward(BASE)
    atp = {}
    for step, count, seed in ATP_CASES:
        print(f"ATP: {count:,} cells at {step:g} h", flush=True)
        atp[step] = cohort(BASE, count, seed, threads=threads, atp_step=step)
    return forward_evidence(out, arrays, atp)


def passage_evidence(out, fitted=()):
    rows, summary = [], []
    choices = [(f"s={s:g}", s) for s in PAPER_MODES] + list(fitted)
    for label, s in choices:
        data = passage(Model(0, s, 1))
        elapsed, cdf = data["days"], data["cdf"]
        gap_max = float(np.max(np.abs(data["absorption"].sum(axis=1) - 1)))
        if gap_max > 1e-9 or np.min(np.diff(cdf)) < -1e-11:
            raise RuntimeError("invalid first-passage probability curve")
        mode_start = int(elapsed[np.argmax(np.diff(cdf))])
        summary.append(
            {
                "label": label,
                "advantage": s,
                "copies": 1000,
                "target": data["target"],
                "success_probability": data["success"],
                "median_days": data["quantiles_days"]["p50"],
                "q2_5_days": data["quantiles_days"]["p2.5"],
                "q97_5_days": data["quantiles_days"]["p97.5"],
                "mode_bin_lower_days": mode_start,
                "paper_mode_bin_lower_days": PAPER_MODES.get(s),
                "conditional_fraction_by_1080": float(np.interp(1080, elapsed, cdf)),
                "no_takeover_yet_given_success": data["right_tail"],
                "max_mass_error": gap_max,
            }
        )
        for age, v, mass in zip(elapsed, cdf, data["absorption"], strict=True):
            rows.append(
                {
                    "label": label,
                    "advantage": s,
                    "elapsed_days": age,
                    "conditional_cdf": v,
                    "loss": mass[0],
                    "takeover": mass[1],
                    "transient": mass[2],
                }
            )
    save_csv(out / "fig5-passage.csv", rows)
    save_csv(out / "takeover-summary.csv", summary)
    return summary


def paper_table_differences(out):
    lookup = {
        (r["regime"], r["parameter"]): r
        for r in read_csv(HERE / "paper-tables.csv")
        if r["table"] in ("2", "3")
    }
    data = []
    for r in read_csv(out / "bootstrap.csv"):
        ref = lookup[r["regime"], r["parameter"]]
        data.append(
            {
                **r,
                "paper_table": ref["table"],
                "paper_cohorts": ref["cohorts"],
                "paper_cells": ref["cells"],
                "paper_mean": ref["mean"],
                "paper_cv_percent": ref["cv_percent"],
                "paper_p2_5": ref["p2_5"],
                "paper_p97_5": ref["p97_5"],
                "mean_difference": float(r["mean"]) - float(ref["mean"]),
            }
        )
    save_csv(out / "parameter-differences.csv", data)
    return data
