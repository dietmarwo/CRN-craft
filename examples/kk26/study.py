#!/usr/bin/env python3
"""Run the mtDNA inverse-model tutorial, bootstrap and figures."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from evidence import compute_forward, paper_table_differences, passage_evidence
from fitting import (
    ROOT_SEED,
    STORED,
    campaign,
    fit,
    load_counts,
    prepare_grid,
    read_csv,
    seed_for,
)
from model import HERE, REGIMES
from render import figures

DEFAULT = HERE / "results" / "local"


def positive(n):
    value = int(n)
    if value < 1:
        raise argparse.ArgumentTypeError("must be positive")
    return value


def main():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="command", required=True)
    boot = sub.add_parser(
        "bootstrap", help="generate synthetic cohorts and fit them independently"
    )
    boot.add_argument("--size", choices=("quick", "full"), default="quick")
    boot.add_argument("--regime", choices=REGIMES, action="append")
    boot.add_argument("--cohorts", type=positive)
    boot.add_argument("--cells", type=positive)
    boot.add_argument("--budget", type=positive)
    boot.add_argument("--seed", type=int, default=ROOT_SEED)
    boot.add_argument("--threads", type=positive, default=1)
    boot.add_argument("--results", type=Path)
    boot.add_argument("--cache", type=Path)

    single = sub.add_parser("fit", help="fit a CSV with age_days and five bin counts")
    single.add_argument("counts", type=Path)
    single.add_argument("--regime", choices=REGIMES, default="baseline")
    single.add_argument("--budget", type=positive, default=1800)
    single.add_argument("--seed", type=int, default=ROOT_SEED)
    single.add_argument("--threads", type=positive, default=1)
    single.add_argument("--cache", type=Path, default=DEFAULT / "grids")

    evidence = sub.add_parser(
        "evidence", help="regenerate forward and first-passage checks"
    )
    evidence.add_argument("--threads", type=positive, default=1)
    evidence.add_argument("--results", type=Path, default=DEFAULT)
    evidence.add_argument("--complete-grid", action="store_true")
    evidence.add_argument("--images", type=Path)

    images = sub.add_parser(
        "figures", help="render the saved Figure 5 and Figure 6 data"
    )
    images.add_argument("--results", type=Path, default=STORED)
    images.add_argument("--images", type=Path, default=DEFAULT / "figures")
    args = p.parse_args()

    if args.command == "fit":
        ages, counts = load_counts(args.counts)
        grid, _ = prepare_grid(args.regime, args.cache, ages, threads=args.threads)
        row = fit(grid, counts, seed_for(args.seed, "input-fit", 0), args.budget)
        print(json.dumps(row, indent=2))
        return
    if args.command == "bootstrap":
        large = args.size == "full"
        reps = args.cohorts or (40 if large else 4)
        cells = args.cells or (10_000 if large else 2_000)
        budget = args.budget or (1800 if large else 640)
        labels = args.regime or (tuple(REGIMES) if large else ("baseline",))
        if len(set(labels)) != len(labels):
            p.error("each regime may be selected once")
        out = args.results or DEFAULT / args.size
        cache = args.cache or DEFAULT / "grids"
        campaign(out, cache, labels, reps, cells, budget, args.seed, args.threads)
        paper_table_differences(out)
        print(f"Results written to {out}")
        return
    if args.command == "evidence":
        args.results.mkdir(parents=True, exist_ok=True)
        compute_forward(args.results, args.threads, args.complete_grid)
        table = args.results / "bootstrap.csv"
        fitted = (
            [
                (f"{r['regime']} fitted", float(r["mean"]))
                for r in read_csv(table)
                if r["parameter"] == "advantage"
            ]
            if table.is_file()
            else ()
        )
        passage_evidence(args.results, fitted)
        figures(args.results, args.images or args.results / "figures")
        return
    figures(args.results, args.images)


if __name__ == "__main__":
    main()
