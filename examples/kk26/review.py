#!/usr/bin/env python3
"""Checks of the recovery experiment: optimizer budget, a wider search box,
Model-1 data and independent cells per age.

Each command writes one ``review-*.csv`` file. The questions, runs and
decision rules were fixed before the runs; guide page 18b gives them.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import replace
from pathlib import Path

import numpy as np
from scipy.optimize import Bounds

from fitting import (
    PARAMETERS,
    ROOT_SEED,
    STORED,
    code_hash,
    fit,
    fit_summary,
    from_x,
    prepare_grid,
    read_csv,
    save_csv,
    seed_for,
    to_x,
)
from model import AGES, REGIMES, cohort, forward, tails

LOCAL = Path(__file__).resolve().parent / "results" / "local"
# Wide boxes in (log10 p, s, f) and the neutral cap they need. The strong
# s range follows the 0..2 fitting range of the authors' companion preprint.
WIDE = {
    "baseline": ([np.log10(2e-6), 0.01, 0.05], [np.log10(3e-5), 0.5, 1.0], 400),
    "strong": ([np.log10(8e-6), 0.10, 0.02], [np.log10(1e-4), 2.0, 1.0], 500),
}
WIDE_COHORTS = (("baseline", 0), ("baseline", 8), ("strong", 0), ("strong", 2))


def saved_cohorts(regime, cache, threads):
    """The saved fit rows and their cohorts, regenerated from the saved seeds."""
    rows = [r for r in read_csv(STORED / "fits.csv") if r["regime"] == regime]
    grid, _ = prepare_grid(regime, cache / "grids", threads=threads)
    work = cache / "cohorts"
    work.mkdir(parents=True, exist_ok=True)
    data = []
    for row in rows:
        seed = int(row["cohort_seed"])
        cells = int(row["cells"])
        file = work / f"{regime}-{cells}-{seed:016x}-{code_hash()[:12]}.npz"
        if file.is_file():
            with np.load(file, allow_pickle=False) as z:
                counts = z["counts"]
        else:
            counts = cohort(grid.model, cells, seed, threads=threads)[0]
            np.savez_compressed(file, counts=counts)
        data.append((row, counts))
    return grid, data


def budget(out, cache, threads):
    result = []
    for regime in REGIMES:
        grid, data = saved_cohorts(regime, cache, threads)
        for row, counts in data:
            i = int(row["cohort"])
            large = fit(grid, counts, int(row["optimizer_seed"]), 18000)
            other = fit(
                grid, counts, seed_for(ROOT_SEED, regime + "-de-review", i), 1800
            )
            entry = {"regime": regime, "cohort": i}
            for tag, values in (
                ("saved", row),
                ("budget18000", large),
                ("seed2", other),
            ):
                for name in (*PARAMETERS, "ssq"):
                    entry[f"{tag}_{name}"] = float(values[name])
            result.append(entry)
    save_csv(out / "review-budget.csv", result)
    for tag in ("budget18000", "seed2"):
        worst = {
            name: max(
                abs(r[f"{tag}_{name}"] - r[f"saved_{name}"])
                / (1 if name == "advantageous" else r[f"saved_{name}"])
                for r in result
            )
            for name in PARAMETERS
        }
        print(tag, "largest change (relative for p and s, absolute for f):", worst)


def wide_box(out, cache, threads, only=None):
    """DE on the direct split solver (15-day step) over a wide box."""
    from fcmaes.decpp import minimize

    result = []
    for regime, i in WIDE_COHORTS:
        if only and (regime, i) not in only:
            continue
        low, high, cap = WIDE[regime]
        grid, data = saved_cohorts(regime, cache, threads)
        row, counts = data[i]
        observed = tails(counts) / counts.sum(axis=1, keepdims=True)
        model = REGIMES[regime]
        memo_file = cache / f"wide-{regime}-{i}-{code_hash()[:12]}.jsonl"
        memo = {}
        if memo_file.is_file():
            with memo_file.open() as fh:
                memo = dict(json.loads(line) for line in fh)

        def curves(x, step, model=model, cap=cap):
            return tails(forward(from_x(model, x), cap=cap, step=step)[:, :5])

        with memo_file.open("a") as fh:

            def objective(x, observed=observed, memo=memo, fh=fh):
                key = np.asarray(x, dtype=float).tobytes().hex()
                if key not in memo:
                    memo[key] = float(np.sum((curves(x, 15.0) - observed) ** 2))
                    fh.write(json.dumps([key, memo[key]]) + "\n")
                    fh.flush()
                return memo[key]

            best = minimize(
                objective,
                bounds=Bounds(low, high),
                popsize=32,
                max_evaluations=960,
                rg=np.random.default_rng(seed_for(ROOT_SEED, regime + "-wide", i)),
                workers=1,
            )
        local = to_x(replace(model, **{name: float(row[name]) for name in PARAMETERS}))
        result.append(
            {
                "regime": regime,
                "cohort": i,
                "box_low": [float(10 ** low[0]), low[1], low[2]],
                "box_high": [float(10 ** high[0]), high[1], high[2]],
                "neutral_cap": cap,
                "evaluations": int(best.nfev),
                "wide_mutation": float(10 ** best.x[0]),
                "wide_advantage": float(best.x[1]),
                "wide_advantageous": float(best.x[2]),
                "saved_mutation": float(row["mutation"]),
                "saved_advantage": float(row["advantage"]),
                "saved_advantageous": float(row["advantageous"]),
                "wide_ssq_5day": float(np.sum((curves(best.x, 5.0) - observed) ** 2)),
                "saved_ssq_5day": float(np.sum((curves(local, 5.0) - observed) ** 2)),
            }
        )
        print(result[-1], flush=True)
    previous = []
    if (out / "review-wide-box.csv").is_file():
        done = {(r["regime"], r["cohort"]) for r in result}
        previous = [
            r
            for r in read_csv(out / "review-wide-box.csv")
            if (r["regime"], int(r["cohort"])) not in done
        ]
    rows = sorted(
        previous + result,
        key=lambda r: [(g, int(c)) for g, c in WIDE_COHORTS].index(
            (r["regime"], int(r["cohort"]))
        ),
    )
    save_csv(out / "review-wide-box.csv", rows)


def atp_fit(out, cache, threads):
    """Fit the saved 100,000-cell hourly ATP cohort with the Moran grid."""
    curves = read_csv(STORED / "fig6-curves.csv")
    cells = 100_000
    exceed = np.array(
        [
            [
                float(r["fraction"])
                for r in curves
                if r["backend"] == "ATP 1 h" and float(r["age_days"]) == age
            ]
            for age in AGES
        ]
    )
    exceed = np.rint(exceed * cells).astype(int)
    counts = np.column_stack(
        [cells - exceed[:, 0], -np.diff(exceed, axis=1), exceed[:, 3]]
    )
    grid, _ = prepare_grid("baseline", cache / "grids", threads=threads)
    rows = [
        {"seed_index": k, **fit(grid, counts, seed_for(ROOT_SEED, "atp-fit", k), 1800)}
        for k in range(4)
    ]
    save_csv(out / "review-atp-fit.csv", rows)
    print(rows[0])


def design(out, cache, threads):
    """Independent cells at every age instead of one cohort followed over time."""
    fits, summary = [], []
    for regime in REGIMES:
        grid, _ = prepare_grid(regime, cache / "grids", threads=threads)
        current = []
        for i in range(40):
            counts = np.empty((len(AGES), 5), dtype=np.int64)
            for j in range(len(AGES)):
                seed = seed_for(ROOT_SEED, f"{regime}-independent-{i}", j)
                counts[j] = cohort(
                    grid.model, 10_000, seed, AGES[: j + 1], threads=threads
                )[0][-1]
            row = {
                "regime": regime,
                "cohort": i,
                **fit(grid, counts, seed_for(ROOT_SEED, regime + "-de", i), 1800),
            }
            current.append(row)
            fits.append(row)
        summary.extend(fit_summary(regime, current, grid.bounds))
    save_csv(out / "review-design-fits.csv", fits)
    save_csv(out / "review-design.csv", summary)
    for row in summary:
        print(
            row["regime"],
            row["parameter"],
            row["cv_percent"],
            row["p2_5"],
            row["p97_5"],
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("check", choices=("budget", "wide-box", "atp-fit", "design"))
    parser.add_argument("--results", type=Path, default=LOCAL / "review")
    parser.add_argument("--cache", type=Path, default=LOCAL / "review-cache")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument(
        "--cohort",
        action="append",
        help="wide-box only: regime:index, for example strong:2 (repeatable)",
    )
    args = parser.parse_args()
    args.results.mkdir(parents=True, exist_ok=True)
    args.cache.mkdir(parents=True, exist_ok=True)
    if args.check == "budget":
        budget(args.results, args.cache, args.threads)
    elif args.check == "wide-box":
        only = None
        if args.cohort:
            only = {(c.split(":")[0], int(c.split(":")[1])) for c in args.cohort}
        wide_box(args.results, args.cache, args.threads, only)
    elif args.check == "atp-fit":
        atp_fit(args.results, args.cache, args.threads)
    else:
        design(args.results, args.cache, args.threads)


if __name__ == "__main__":
    main()
