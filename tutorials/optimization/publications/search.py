"""Evaluate fcmaes ask/tell populations with exact candidate accounting."""

from concurrent.futures import ThreadPoolExecutor
from contextlib import nullcontext
import csv
import json
from dataclasses import dataclass
from pathlib import Path
import time

from fcmaes import cmaes, de, mode
from loguru import logger
import numpy as np
from scipy.optimize import Bounds

from common import arm_manifest, require

PACKAGES = ("crn-craft", "fcmaes", "numpy", "scipy", "matplotlib")
PENALTY = 1e50
for module_name in ("fcmaes.de", "fcmaes.cmaes", "fcmaes.mode"):
    logger.disable(module_name)


@dataclass
class Search:
    candidates: np.ndarray
    scores: np.ndarray
    best: np.ndarray
    value: float
    failures: int


def options_resume(parser):
    parser.add_argument(
        "--resume-from",
        type=Path,
        help="Replay an earlier output's saved populations and extend each run.",
    )


def search(
    args,
    name,
    algorithm,
    function,
    bounds,
    seed,
    budget,
    *,
    pop=10,
    x0=None,
    checkpoints=(),
    observables=1,
    context=None,
    sigma=0.3,
):
    """Keep proposals in evaluation order. A final partial batch is scored only."""
    rng = np.random.default_rng(seed)
    bounds = Bounds(*np.asarray(bounds, dtype=float).T)
    pop = min(pop, budget)
    require(pop >= 4, "search needs at least four evaluations")
    if algorithm == "de":
        opt = de.DE(len(bounds.lb), bounds=bounds, popsize=pop, rg=rng)
    elif algorithm == "cmaes":
        start = rng.uniform(bounds.lb, bounds.ub) if x0 is None else np.asarray(x0)
        opt = cmaes.Cmaes(
            bounds=bounds,
            x0=start,
            popsize=pop,
            rg=rng,
            input_sigma=sigma,
            randn=lambda *shape: rng.standard_normal(shape),
            repair_z=True,
        )
    elif algorithm == "mode":
        opt = mode.MODE(2, 0, bounds, popsize=pop, rg=rng, nsga_update=True)
    elif algorithm == "random":
        opt = None
    else:
        raise ValueError(f"unknown optimizer {algorithm}")

    def evaluate(x):
        try:
            score = np.atleast_1d(function(x)).astype(float)
            if score.shape != (observables,) or not np.isfinite(score).all():
                return np.full(observables, PENALTY), "nonfinite score"
        except (ArithmeticError, RuntimeError, ValueError) as error:
            return np.full(observables, PENALTY), str(error)
        return score, ""

    path = Path(args.output) / f"{name}.csv"
    previous = []
    resume = getattr(args, "resume_from", None)
    if resume is not None and (resume / path.name).is_file():
        meta = json.loads((resume / f"{name}.run.json").read_text())
        require(
            meta["seed"] == seed
            and meta["algorithm"] == algorithm
            and meta["population"] == pop
            and meta.get("context") == context
            and meta.get("sigma", 0.3) == sigma,
            "saved run has a different optimizer or objective context",
        )
        with (resume / path.name).open(newline="") as stream:
            previous = list(csv.DictReader(stream))
        require(
            len(previous) == meta["actual_evaluations"], "incomplete saved run trace"
        )
    X, Y = [], []
    failures = 0
    start = time.perf_counter()
    pool = ThreadPoolExecutor(args.threads) if args.threads > 1 else nullcontext()
    with pool as executor, path.open("w", newline="") as output:
        writer = csv.writer(output)
        writer.writerow(
            [
                "call",
                *[f"x{i}" for i in range(len(bounds.lb))],
                *[f"f{i}" for i in range(observables)],
                "error",
            ]
        )
        while len(X) < budget:
            points = (
                rng.uniform(bounds.lb, bounds.ub, (pop, len(bounds.lb)))
                if opt is None
                else np.asarray(opt.ask()).copy()
            )
            points = points[: budget - len(X)]
            mapped = map if executor is None else executor.map
            reused = min(len(points), max(0, len(previous) - len(X)))
            old = []
            for j, point in enumerate(points[:reused]):
                row = previous[len(X) + j]
                stored = np.array([float(row[f"x{k}"]) for k in range(len(point))])
                require(
                    int(row["call"]) == len(X) + j + 1
                    and np.array_equal(point, stored),
                    "saved coordinates differ from replayed proposals",
                )
                old.append(
                    (
                        np.array([float(row[f"f{k}"]) for k in range(observables)]),
                        row["error"],
                    )
                )
            fresh = list(mapped(evaluate, points[reused:]))
            for point, (value, error) in zip(points, old + fresh, strict=True):
                X.append(point)
                Y.append(value)
                failures += bool(error)
                writer.writerow([len(X), *point, *value, error])
            output.flush()
            if opt is not None and len(points) == pop:
                opt.tell(
                    np.asarray(Y[-pop:]).reshape((pop, observables))
                    if observables > 1
                    else np.asarray(Y[-pop:]).ravel()
                )
    X, Y = np.asarray(X), np.asarray(Y)
    i = int(np.argmin(Y[:, 0]))
    arm_manifest(
        args.output,
        name,
        budget,
        len(X),
        time.perf_counter() - start,
        seed,
        workers=args.threads,
        formulation="scalar" if observables == 1 else "pareto",
        algorithm=algorithm,
        backend="numpy" if opt is None else "fcmaes",
        population=pop,
        failures=failures,
        trace=path.name,
        context=context,
        sigma=sigma,
        replayed=min(budget, len(previous)),
        checkpoints={n: float(np.min(Y[:n, 0])) for n in checkpoints if n <= len(X)},
    )
    return Search(X, Y, X[i], float(Y[i, 0]), failures)


def check(report, name, condition, **details):
    report["checks"].append(dict(name=name, passed=bool(condition), **details))


def save_plot(args, filename, draw):
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure

    figure = Figure(figsize=(8, 4.5), layout="constrained")
    FigureCanvasAgg(figure)
    draw(figure)
    figure.savefig(Path(args.output) / filename, dpi=130)
