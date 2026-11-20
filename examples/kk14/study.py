#!/usr/bin/env python3
"""Run the mitochondrial DNA mechanism comparison and regenerate its figures."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from campaign import AGES, BASE_SEED, STORED, TIERS, compare, save_csv
from reference import (
    HERE,
    MECHANISMS,
    atp_equilibrium,
    size_capacity,
    solve_atp,
    solve_delay,
)
from render import (
    comparison_table,
    reference_differences,
    single_cell_figure,
    tissue_figures,
)


def single_cell(out, images=None):
    out.mkdir(parents=True, exist_ok=True)
    atp = solve_atp()
    size = solve_delay()
    other = {
        "mesh_0.5_h": solve_delay(step=0.5),
        "mesh_0.125_h": solve_delay(step=0.125),
        "initiation_0.1": solve_delay(initiation=0.1),
        "half_life_2_days": solve_delay(half_life=48, max_years=10),
        "wild_replication_10_h": solve_delay(wild_delay=10, max_years=50),
        "mutant_size_0.7": solve_delay(step=0.1, mutant_size=0.7),
        "identical_delays": solve_delay(mutant_size=1, max_years=3),
    }
    for path, data, cols in (
        ("atp-trajectory.csv", atp, "hour,wild,mutants,ATP"),
        (
            "size-trajectory.csv",
            size,
            "hour,wild_free,wild_busy,mutant_free,mutant_busy",
        ),
    ):
        np.savetxt(
            out / path, data, header=cols, delimiter=",", comments="", fmt="%.12g"
        )
    delta, initial_atp, scale = atp_equilibrium()
    result = {
        "atp_collapse_days": float(atp[-1, 0] / 24),
        "atp_mutant_wild_ratio_at_collapse": float(atp[-1, 2] / atp[-1, 1]),
        "atp_equilibrium": initial_atp,
        "atp_scale": scale,
        "atp_method": {
            "solver": "Radau",
            "rtol": 1e-9,
            "atol": 1e-11,
            "max_step_h": 12,
        },
        "size_capacity": size_capacity(1000, 240, 2, 1),
        "size_extinction_hours": float(size[-1, 0]),
        "size_extinction_years": float(size[-1, 0] / 8760),
        "size_wild_busy_equilibrium": math.expm1(2 * delta),
        "size_mutant_busy_equilibrium": math.expm1(delta),
        "size_final_mutant_busy_fraction": float(
            size[-1, 4] / (size[-1, 3] + size[-1, 4])
        ),
        "size_checks_years": {
            name: float(points[-1, 0] / 8760)
            for name, points in other.items()
            if name != "identical_delays"
        },
        "same_size_final_ratio": float(
            (other["identical_delays"][-1, 1:3].sum())
            / (other["identical_delays"][-1, 3:].sum())
        ),
    }
    (out / "single-cell.json").write_text(json.dumps(result, indent=2) + "\n")
    if images is not None:
        single_cell_figure(atp, size, images)
    return result


def _positive(text):
    value = int(text)
    if value < 1:
        raise argparse.ArgumentTypeError("must be at least one")
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    one = sub.add_parser("single-cell", help="verify both deterministic references")
    one.add_argument("--size", choices=TIERS, default="quick")
    one.add_argument("--results", type=Path)
    one.add_argument("--images", type=Path)
    together = sub.add_parser(
        "compare", help="calibrate and validate the tissue comparison"
    )
    together.add_argument("--size", choices=TIERS, default="quick")
    together.add_argument("--ages", nargs="+", choices=AGES, type=int)
    together.add_argument(
        "--mechanisms", nargs="+", choices=MECHANISMS, default=MECHANISMS
    )
    together.add_argument("--threads", type=_positive, default=1)
    together.add_argument("--results", type=Path)
    together.add_argument("--images", type=Path)
    together.add_argument("--resume", action="store_true")
    together.add_argument("--seed", type=int, default=BASE_SEED)
    together.add_argument("--sensitivities", action="store_true")
    together.add_argument("--from-rates", type=Path)
    together.add_argument("--validation-cells", type=_positive)
    figures = sub.add_parser(
        "figures", help="regenerate the figures and comparison table"
    )
    figures.add_argument("--results", type=Path, default=STORED)
    figures.add_argument("--images", type=Path)
    figures.add_argument("--table", type=Path)
    native = sub.add_parser(
        "native", help="exercise the CRN-craft chain and ATP models"
    )
    native.add_argument("--threads", type=_positive, default=2)
    native.add_argument("--results", type=Path, default=HERE / "results" / "native")
    args = parser.parse_args()
    if args.command == "native":
        from native import atp_comparison, marked_hybrid, staged_comparison

        args.results.mkdir(parents=True, exist_ok=True)
        marks = marked_hybrid(threads=args.threads)
        (args.results / "native-marks.json").write_text(
            json.dumps(marks, indent=2) + "\n"
        )
        save_csv(args.results / "native-atp.csv", atp_comparison(threads=args.threads))
        save_csv(
            args.results / "native-size.csv", staged_comparison(threads=args.threads)
        )
        print(f"Native results written to {args.results}")
        return

    out = args.results or HERE / "results" / args.size
    if args.command == "single-cell":
        result = single_cell(out, args.images or out / "figures")
        print(
            f"ATP zero at {result['atp_collapse_days']:.3f} days; wild type <1 at "
            f"{result['size_extinction_years']:.3f} years."
        )
        return
    if args.command == "compare":
        if args.from_rates is not None and not args.from_rates.is_file():
            parser.error("--from-rates must name an existing rate table")
        ages = tuple(sorted(set(args.ages or ((3,) if args.size == "quick" else AGES))))
        methods = tuple(name for name in MECHANISMS if name in args.mechanisms)
        manifest = compare(
            out,
            args.size,
            ages,
            args.threads,
            args.resume,
            args.seed,
            args.sensitivities or args.size == "full",
            args.from_rates,
            args.validation_cells,
            methods,
        )
        print(f"Three-year criterion: {manifest['three_year_gate']}")
    if not (out / "comparison.csv").is_file():
        parser.error(f"no comparison.csv in {out}")
    tissue_figures(out, args.images or out / "figures")
    reference_differences(out)
    path = (
        args.table
        if args.command == "figures" and args.table
        else out / "comparison-table.md"
    )
    path.write_text(comparison_table(out))


if __name__ == "__main__":
    main()
