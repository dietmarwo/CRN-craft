"""C4: infer an aggregation kernel using the native dynamic-species engine."""

import hashlib
import json
import math
import time
import numpy as np
from scipy.linalg import expm
from scipy.optimize import minimize_scalar
import fcmaes_rust
from common import ROOT, arm_manifest, call_json, require, write


START, LOWER, UPPER = 0.6, 0.05, 3.0  # rate search start and bounds


def partitions(n, maximum=None):
    if n == 0:
        yield ()
    else:
        for k in range(min(n, maximum if maximum is not None else n), 0, -1):
            for rest in partitions(n - k, k):
                yield (k, *rest)


def exact_law(monomers, rate, times, multiplicative=False):
    """Independent finite CTMC on integer partitions, using a matrix exponential."""
    states = list(partitions(monomers))
    index = {state: i for i, state in enumerate(states)}
    q = np.zeros((len(states), len(states)))
    for state, i in index.items():
        for a in range(len(state)):
            for b in range(a + 1, len(state)):
                target = [x for j, x in enumerate(state) if j not in (a, b)]
                target.append(state[a] + state[b])
                value = rate * (state[a] * state[b] if multiplicative else 1)
                q[i, index[tuple(sorted(target, reverse=True))]] += value
                q[i, i] -= value
    require(
        np.max(np.abs(q.sum(axis=1))) < 1e-10, "generator does not conserve probability"
    )
    p0 = np.zeros(len(states))
    p0[index[(1,) * monomers]] = 1
    rows = []
    for t in times:
        p = p0 @ expm(q * t)
        require(
            np.min(p) > -1e-12 and abs(p.sum() - 1) < 1e-10, "invalid exact CTMC law"
        )
        rows.append(
            [
                np.cumsum(
                    np.bincount(
                        [f(s) for s in states], weights=p, minlength=monomers + 1
                    )
                )
                for f in [len, max]
            ]
        )
    return np.asarray(rows)


def exact_optimum(monomers, target, times, multiplicative, lower, upper):
    """Best rate of a kernel against the exact training law: the population
    reference that the finite-path search approximates."""
    answer = minimize_scalar(
        lambda x: float(
            np.mean(
                (exact_law(monomers, np.exp(x), times, multiplicative) - target) ** 2
            )
        ),
        bounds=(lower, upper),
        method="bounded",
        options=dict(xatol=1e-10),
    )
    require(answer.success, "exact-law reference optimization failed")
    return float(np.exp(answer.x)), float(answer.fun)


def empirical(result, monomers):
    names = result["observables"]
    paths = np.asarray([p["values"] for p in result["paths"]])
    require(
        np.all(paths[:, :, names.index("mass")] == monomers), "aggregation lost mass"
    )
    return np.asarray(
        [
            [
                np.mean(
                    paths[:, t, names.index(name), None] <= np.arange(monomers + 1),
                    axis=0,
                )
                for name in ["clusters", "largest"]
            ]
            for t in range(paths.shape[1])
        ]
    )


def run(args, report):
    """Fit a scalar rate, then quantify wrong-kernel predictions at held-out times."""
    budget = 16 if args.preset == "smoke" else 64
    paths = 64 if args.preset == "smoke" else 256
    validation_paths = 1024 if args.preset == "smoke" else 8192
    monomers, truth = 6, 1.2
    train_times, test_times = [0.15, 0.4], [0.25, 0.7]
    target = exact_law(monomers, truth, train_times)
    report["budgets"] = dict(
        requested_calls=budget,
        training_paths=paths,
        validation_paths=validation_paths,
        max_events_per_path=10,
        max_instances=16,
        seconds_per_call=20,
        monomers=monomers,
        training_seed=1001,
        validation_seed=9001,
        family_alpha=0.001,
        marginal_checks=4,
    )
    counter = 0
    work = []

    def simulate(rate, multiplicative, times, n, seed):
        nonlocal counter
        tag = f"native-{counter:04d}"
        counter += 1
        folder = args.output / tag
        folder.mkdir()
        filename = (
            "aggregation-multiplicative.rules.json"
            if multiplicative
            else "aggregation-constant.rules.json"
        )
        model = json.loads((ROOT / "models/rules" / filename).read_text())
        model["parameters"].update(k=float(rate), monomers=monomers)
        model_path = folder / "model.json"
        write(model_path, model)
        request = dict(
            schema="crn-craft.rules-run/0.1",
            model=dict(
                path="model.json",
                sha256=hashlib.sha256(model_path.read_bytes()).hexdigest(),
            ),
            time=dict(end=times[-1], report_times=times),
            ensemble=dict(replicates=n, root_seed=seed),
            outputs=dict(
                statistics=["mean", "variance"], keep_paths=n, keep_instances=0
            ),
            budget=dict(
                max_events_per_replicate=10,
                max_sampling_steps_per_replicate=10000,
                max_instances=16,
                threads=args.threads,
                max_output_bytes=16 * 1024 * 1024,
            ),
        )
        write(folder / "request.json", request)
        result = call_json(
            [args.binary.resolve(), "rules", model_path, folder / "request.json"], 20
        )
        write(folder / "result.json", result)
        require(result["termination"] == "completed", "rules run incomplete")
        work.append(
            dict(
                path=tag,
                rate=rate,
                multiplicative=multiplicative,
                paths=n,
                events=result["events"],
                seed=seed,
            )
        )
        report["native_runs"] = work
        write(args.output / "summary.json", report)
        return empirical(result, monomers)

    epsilon = math.sqrt(math.log(2 * 4 / 0.001) / (2 * validation_paths))
    for multiplicative in [False, True]:
        calls = []
        start = time.perf_counter()

        def objective(x, multiplicative=multiplicative, calls=calls):
            rate = float(np.exp(x[0]))
            cdf = simulate(rate, multiplicative, train_times, paths, 1001)
            value = float(np.mean((cdf - target) ** 2))
            calls.append(dict(rate=rate, loss=value))
            return value

        answer = fcmaes_rust.optimize_bite(
            objective,
            [math.log(START)],
            [math.log(LOWER)],
            [math.log(UPPER)],
            seed=713,
            max_evaluations=budget,
        )
        rate = float(np.exp(answer[0][0]))
        reference_rate, reference_loss = exact_optimum(
            monomers,
            target,
            train_times,
            multiplicative,
            math.log(LOWER),
            math.log(UPPER),
        )
        seconds = time.perf_counter() - start
        held = simulate(rate, multiplicative, test_times, validation_paths, 9001)
        law = exact_law(monomers, rate, test_times, multiplicative)
        distance = float(np.max(np.abs(held - law)))
        require(
            distance <= epsilon,
            "held-out rules paths fail the predeclared DKW family bound",
        )
        row = dict(
            kernel="multiplicative" if multiplicative else "constant",
            fitted_rate=rate,
            true_rate=truth,
            training_loss=answer[1],
            best_is_start=bool(math.isclose(rate, START, rel_tol=1e-12)),
            exact_law_optimal_rate=reference_rate,
            exact_law_optimal_loss=reference_loss,
            actual_evaluations=len(calls),
            exact_heldout_prediction_mse=float(
                np.mean((law - exact_law(monomers, truth, test_times)) ** 2)
            ),
            heldout_law_max_cdf_distance=distance,
            dkw_bound=epsilon,
            search_seconds=seconds,
            trace=calls,
        )
        report["arms"].append(row)
        report["checks"].append(
            dict(name=row["kernel"] + "-mass-and-exact-law", passed=True)
        )
        arm_manifest(
            args.output,
            "rules-" + row["kernel"],
            budget,
            len(calls),
            seconds,
            713,
            args.threads,
        )
    report["interpretation"] = (
        "Kernel prediction errors are reported, not a guaranteed model-selection result at this pilot budget."
    )
