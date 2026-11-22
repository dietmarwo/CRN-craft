"""Compare MODE on the published virtual-reactor benchmark."""

import csv
from itertools import product

import numpy as np
from scipy.integrate import solve_ivp

from common import write
from .concentrations import Solver, document
from .search import PACKAGES as PACKAGES
from .search import check, save_plot, search
from .sources import data_option, fixtures

DOI = "10.1002/cmtd.202000051"
BOUNDS = ((0.5, 2.0), (1.0, 5.0), (0.1, 0.5), (30.0, 120.0))
HV_POINT = np.array([-2957.0, 10.7])
MOLAR = np.array([159.09, 71.12, 210.21, 210.21, 261.33])
MODEL = document(
    ("A", "B", "P", "Q", "R"),
    ("a0", "b0", "ka", "kb", "kc", "kd"),
    ("a0", "b0", "0", "0", "0"),
    ("-v1-v2", "-v1-v2-v3-v4", "v1-v3", "v2-v4", "v3+v4"),
    assignments=(
        ("v1", "ka*A*B"),
        ("v2", "kb*A*B"),
        ("v3", "kc*P*B"),
        ("v4", "kd*Q*B"),
    ),
)


def inputs(design):
    _, equivalents, a0, temp = design
    delta = 1 / (temp + 273.71) - 1 / (90 + 273.71)
    kinetic = (
        0.6
        * np.array([57.9, 2.70, 0.865, 1.63])
        * np.exp(-np.array([33.3, 35.3, 38.9, 44.8]) * delta / 0.008314)
    )
    return (a0, equivalents * a0, *kinetic)


def outlets_reference(design, *, loose=False):
    a0, b0, ka, kb, kc, kd = inputs(design)

    def rhs(_, state):
        a, b, p, q, _ = state
        if loose and a < 1e-6 * a0:
            a = 0.0
        if loose and b < 1e-6 * b0:
            b = 0.0
        consumption = (ka + kb) * a * b
        return [
            -consumption,
            -consumption - kc * b * p - kd * b * q,
            ka * a * b - kc * b * p,
            kb * a * b - kd * b * q,
            kc * b * p + kd * b * q,
        ]

    options = {} if loose else dict(rtol=2e-11, atol=1e-13)
    solution = solve_ivp(rhs, [0.0, design[0]], [a0, b0, 0.0, 0.0, 0.0], **options)
    if not solution.success:
        raise RuntimeError(solution.message)
    return solution.y[:, -1]


def metrics(residence, outlets):
    outlet = np.maximum(outlets, 0.0)
    product_mass = MOLAR[2] * outlet[2]
    other_mass = MOLAR @ outlet - product_mass
    sty = max(60 * product_mass / residence, 1e-6)
    waste = (
        1000.0
        if np.isclose(outlet[2], 0.0)
        else min(1000.0, (789.0 + other_mass) / product_mass)
    )
    return np.array([-sty, waste])


def hypervolume(points, point=HV_POINT):
    """Area of the union of minimizing rectangles anchored at point."""
    ordered = points[np.argsort(points[:, 0], kind="stable")]
    previous = point[1]
    area = 0.0
    for x, y in ordered:
        if x < point[0] and y < previous:
            area += (point[0] - x) * (previous - y)
            previous = y
    return float(area)


def options(parser):
    data_option(parser)


def run(args, report):
    """Verify the virtual reactor, then run the published evaluation budget."""
    sources = fixtures(args.data_dir, "C7")
    with sources["snar_results.csv"].open(newline="") as archive:
        archived = list(csv.DictReader(archive))
    noise = sorted(set(row["noise_level"] for row in archived))
    check(
        report,
        "archived-reference-noise",
        noise == ["0"],
        rows=len(archived),
        noise=noise,
    )
    selected = [
        r for r in archived if r["noise_level"] == "0" and r["batch_size"] == "1"
    ]
    groups = {}
    for row in selected:
        key = (
            row["strategy_name"],
            row["transform_name"],
            row["sty_tolerance"],
            row["e_factor_tolerance"],
            row["expression"],
            row["maximize"],
        )
        groups.setdefault(key, []).append(float(row["terminal_hypervolume"]))
    archived_scores = [
        dict(
            strategy=k[0],
            transform=k[1],
            sty_tolerance=k[2],
            e_factor_tolerance=k[3],
            expression=k[4],
            maximize=k[5],
            runs=len(v),
            median=np.median(v),
            quantiles=np.quantile(v, [0.25, 0.75]),
        )
        for k, v in groups.items()
    ]
    rng = np.random.default_rng(4984)
    lows, highs = np.asarray(BOUNDS).T
    designs = np.array([*product(*BOUNDS), *rng.uniform(lows, highs, (8, 4))])
    rows, results, series = [], [], []
    starts = 2 if args.preset == "smoke" else 20
    with Solver(args.ode_binary, args.output / "model.json", MODEL) as engine:

        def outlets(design, *, tight=False):
            state, _ = engine.batch(
                [inputs(design)],
                [0.0, design[0]],
                **(dict(rtol=2e-10, atol=1e-12) if tight else {}),
            )
            return state[0, -1]

        for design in designs:
            native = outlets(design)
            fine = outlets_reference(design)
            baseline = outlets_reference(design, loose=True)
            rows.append(
                dict(
                    design=design,
                    outputs=metrics(design[0], native),
                    native_minus_tight=native - fine,
                    reference_minus_tight=baseline - fine,
                    output_gap=metrics(design[0], native) - metrics(design[0], fine),
                    source_output_gap=metrics(design[0], baseline)
                    - metrics(design[0], fine),
                )
            )
        state_err = float(np.max(np.abs([row["native_minus_tight"] for row in rows])))
        score_err = float(np.max(np.abs([row["output_gap"] for row in rows])))
        check(
            report,
            "reactor-outlet-parity",
            state_err < 1e-6 and score_err < 0.02,
            state_max_abs=state_err,
            objective_max_abs=score_err,
        )
        example = np.array([[1.0, 3.0], [2.0, 2.0], [3.0, 4.0], [6.0, 0.0]])
        check(
            report,
            "hypervolume-hand-fixture",
            hypervolume(example, np.array([5.0, 5.0])) == 11.0,
        )
        for name in ("mode", "random"):
            for seed in range(starts):
                fit = search(
                    args,
                    f"{name}-{seed}",
                    name,
                    lambda x: metrics(x[0], outlets(x)),
                    BOUNDS,
                    9971 + seed + args.seed_offset,
                    50,
                    pop=10,
                    observables=2,
                )
                hv = [hypervolume(fit.scores[:n]) for n in range(1, 51)]
                series.append((name, hv))
                results.append(
                    dict(
                        algorithm=name,
                        seed=seed,
                        hypervolume=hv[-1],
                        evaluations=len(fit.scores),
                        failures=fit.failures,
                    )
                )
        check(
            report,
            "exact-50-evaluations",
            len(results) == 2 * starts
            and all(r["evaluations"] == 50 and not r["failures"] for r in results),
        )
        physical = np.array([r["outputs"] for r in rows])
        rescore = max(
            np.max(np.abs(metrics(d[0], outlets(d, tight=True)) - v))
            for d, v in zip(designs, physical, strict=True)
        )
        report["tight_objective_difference"] = rescore
    report.update(
        doi=DOI,
        reference_point=HV_POINT,
        stages=results,
        designs=rows,
        archived_groups=archived_scores,
        terminal={
            n: dict(
                median=np.median(
                    [r["hypervolume"] for r in results if r["algorithm"] == n]
                ),
                q25_q75=np.quantile(
                    [r["hypervolume"] for r in results if r["algorithm"] == n],
                    [0.25, 0.75],
                ),
            )
            for n in ("mode", "random")
        },
    )
    write(args.output / "hypervolume.json", dict(series=series))
    with sources["pareto_front_snar.csv"].open(newline="") as file:
        reader = csv.DictReader(file)
        next(reader)  # Summit's type row.
        front = np.array([[float(q["sty"]), float(q["e_factor"])] for q in reader])

    def draw(fig):
        front_axis, trace = fig.subplots(1, 2)
        front_axis.scatter(
            front[:, 0], front[:, 1], c="#d4d4d4", s=3, label="Summit sample"
        )
        match = []
        for i in range(starts):
            with (args.output / f"mode-{i}.csv").open(newline="") as h:
                match.extend(
                    [
                        [-float(r["f0"]), float(r["f1"])]
                        for r in csv.DictReader(h)
                        if not r["error"]
                    ]
                )
        match = np.array(match)
        front_axis.scatter(match[:, 0], match[:, 1], s=7, label="MODE")
        front_axis.set(xlabel="Space–time yield", ylabel="E-factor", ylim=(0, 40))
        front_axis.legend(fontsize=7)
        calls = np.arange(1, 51)
        for name in ("mode", "random"):
            values = np.array([y for n, y in series if n == name])
            lo, med, hi = np.quantile(values, [0.25, 0.5, 0.75], axis=0)
            trace.plot(calls, med, label=name)
            trace.fill_between(calls, lo, hi, alpha=0.2)
        trace.set(xlabel="Evaluations", ylabel="Hypervolume")
        trace.legend()

    save_plot(args, "pareto-and-hypervolume.png", draw)
