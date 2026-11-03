"""C2: nonlinear BDF sensitivity fitting, independent ODE checks and SSA transfer."""

import copy
import json
import time
import numpy as np
from scipy.integrate import solve_ivp
import crn_craft as crn
import fcmaes_rust
import sqp_craft as sqp
from common import ROOT, arm_manifest, call_json, checked, digest, require, write

NAMES = ["kon", "kcat"]
TIMES = np.linspace(0, 4, 9)
SCALE = np.array([50, 10, 10, 50])


def independent(rates, times=TIMES):
    kon, kcat = rates

    def rhs(_t, y):
        s, e, es, _p = y
        bind, release, cat = kon * s * e, 4 * es, kcat * es
        return [-bind + release, -bind + release + cat, bind - release - cat, cat]

    result = solve_ivp(
        rhs,
        [0, times[-1]],
        [50, 10, 0, 0],
        t_eval=times,
        method="DOP853",
        rtol=1e-11,
        atol=1e-12,
    )
    require(result.success, "independent DOP853 solve failed")
    return result.y.T


class Objective:
    def __init__(self, binary, model, folder, target, cap, rtol=1e-9):
        self.binary, self.model, self.folder = binary, model, folder
        self.target, self.cap, self.rtol = target, cap, rtol
        self.calls, self.attempts, self.cache = 0, 0, None
        self.best, self.best_value = None, float("inf")
        folder.mkdir()
        write(
            folder / "request.json",
            dict(
                schema="crn-craft.ode-request/0.1",
                times=TIMES.tolist(),
                rtol=rtol,
                atol=rtol * 0.01,
            ),
        )

    def evaluate(self, x):
        if self.cache is not None and np.array_equal(x, self.cache[0]):
            return self.cache[1:]
        self.attempts += 1
        if self.calls >= self.cap:
            raise sqp.UndefinedPointError("predeclared native-solve budget exhausted")
        model = copy.deepcopy(self.model)
        model["parameters"].update(dict(zip(NAMES, np.exp(x).tolist(), strict=True)))
        write(self.folder / "model.json", model)
        self.calls += 1
        write(self.folder / f"{self.calls:04d}.model.json", model)
        result = call_json(
            [
                self.binary.resolve(),
                self.folder / "model.json",
                self.folder / "request.json",
            ],
            30,
        )
        write(self.folder / f"{self.calls:04d}.result.json", result)
        order = [result["parameters"].index(n) for n in NAMES]
        values = np.asarray(result["values"])
        sensitivity = np.asarray(result["sensitivities"])[order].transpose(
            1, 2, 0
        ) * np.exp(x)
        residual = (values - self.target) / SCALE
        value = float(np.mean(residual**2))
        gradient = (
            2 * np.einsum("ts,tsp->p", residual / SCALE, sensitivity) / residual.size
        )
        self.cache = (np.asarray(x).copy(), value, gradient, values, sensitivity)
        if value < self.best_value:
            self.best, self.best_value = np.asarray(x).copy(), value
        return self.cache[1:]

    def functions(self, x):
        return self.evaluate(x)[0], []

    def derivatives(self, x):
        return self.evaluate(x)[1].tolist(), []


def run(args, report):
    """Fit stiff nonlinear enzyme kinetics and separate accuracy, derivatives and identification."""
    model = json.loads(
        (ROOT / "models/examples/enzyme-mechanism.model.json").read_text()
    )
    model["parameters"].update(kon=1.0, koff=4.0, kcat=0.8)
    truth = np.array([1.0, 0.8])
    target = independent(truth)
    budget = 48 if args.preset == "smoke" else 160
    report["budgets"] = dict(
        native_solves_per_strategy=budget,
        seconds_per_solve=30,
        validation_paths=256 if args.preset == "smoke" else 4096,
    )
    report["ode_binary_sha256"] = digest(args.ode_binary)
    write(args.output / "model.json", model)
    write(args.output / "target.json", target)
    audit = Objective(args.ode_binary, model, args.output / "audit", target, 30)
    point = np.log([0.8, 1.1])
    value, gradient, values, jac = audit.evaluate(point)
    h = 1e-4
    numeric = np.array(
        [
            (audit.evaluate(point + h * e)[0] - audit.evaluate(point - h * e)[0])
            / (2 * h)
            for e in np.eye(2)
        ]
    )
    gradient_error = float(np.max(np.abs(gradient - numeric)))
    integration_error = float(
        np.max(np.abs(values - independent(np.exp(point))) / SCALE)
    )
    tight = Objective(
        args.ode_binary, model, args.output / "tight-audit", target, 1, rtol=1e-11
    )
    tightened_error = float(
        np.max(np.abs(tight.evaluate(point)[2] - independent(np.exp(point))) / SCALE)
    )
    require(
        gradient_error < 2e-5,
        "native log-coordinate gradient disagrees with central differences",
    )
    require(
        integration_error < 2e-6 and tightened_error < 2e-7,
        "BDF disagrees with independent ODE integration",
    )
    full_rank = int(np.linalg.matrix_rank((jac / SCALE[None, :, None]).reshape(-1, 2)))
    terminal_rank = int(np.linalg.matrix_rank(jac[-1, 3:4, :]))
    require(
        full_rank == 2 and terminal_rank == 1,
        "identifiability rank demonstration failed",
    )
    report["audit"] = dict(
        gradient_error=gradient_error,
        integration_error=integration_error,
        tightened_error=tightened_error,
        full_observation_rank=full_rank,
        terminal_product_rank=terminal_rank,
        native_solves=audit.calls + tight.calls,
    )
    report["checks"].append(
        dict(name="independent-integration-log-gradient-and-rank", passed=True)
    )
    lower, upper = np.log([0.2, 0.1]), np.log([3.0, 2.0])
    spec = sqp.ProblemSpec(
        lower=lower.tolist(),
        upper=upper.tolist(),
        constraint_lower=[],
        constraint_upper=[],
        jacobian_indptr=[0],
        jacobian_indices=[],
    )
    for strategy in ["global", "multistart", "global-local"]:
        start = time.perf_counter()
        status = []
        objective = Objective(
            args.ode_binary, model, args.output / strategy, target, budget
        )
        if strategy != "multistart":
            global_budget = budget if strategy == "global" else budget // 3

            def fitness(x, objective=objective):
                try:
                    return objective.evaluate(x)[0]
                except sqp.UndefinedPointError:
                    return float("inf")

            answer = fcmaes_rust.optimize_bite(
                fitness, point, lower, upper, seed=7101, max_evaluations=global_budget
            )
            status.append(
                dict(method="BiteOpt", requested=global_budget, actual=answer[2])
            )
        if strategy != "global":
            starts = (
                [np.log([0.4, 0.3]), np.log([2.0, 1.5])]
                if strategy == "multistart"
                else [objective.best]
            )
            for i, initial in enumerate(starts):
                objective.cap = (i + 1) * budget // len(starts)
                local = sqp.solve(
                    spec,
                    objective,
                    initial.tolist(),
                    method="sqp",
                    method_options=sqp.SqpOptions(
                        qp_solver="active-set",
                        max_iterations=40,
                        stationarity_tolerance=1e-7,
                    ),
                )
                status.append(dict(method="SQP", status=local.status))
        require(objective.best is not None, "strategy found no finite candidate")
        rates = np.exp(objective.best)
        search_seconds = time.perf_counter() - start
        native_calls = objective.calls
        request = crn.ssa_spec(
            model,
            aggregator="cached_direct",
            end=4,
            save_times=TIMES.tolist(),
            replicates=report["budgets"]["validation_paths"],
            seed=9103,
            threads=args.threads,
            keep=2,
            max_events=100000,
        )
        template = crn.prepare_parameters(model, request)
        params = dict(model["parameters"])
        params.update(dict(zip(NAMES, rates.tolist(), strict=True)))
        stochastic = checked(template.simulate(params))
        write(args.output / f"{strategy}.ssa.json", stochastic)
        mean = np.asarray(stochastic["statistics"]["mean"])
        require(
            np.max(np.abs(mean[:, 0] + mean[:, 2] + mean[:, 3] - 50)) < 1e-9,
            "substrate conservation",
        )
        require(
            np.max(np.abs(mean[:, 1] + mean[:, 2] - 10)) < 1e-9, "enzyme conservation"
        )
        row = dict(
            strategy=strategy,
            rates=rates,
            native_loss=objective.best_value,
            independent_loss=float(
                np.mean(((independent(rates) - target) / SCALE) ** 2)
            ),
            native_solves=native_calls,
            requested_native_solve_cap=budget,
            status=status,
            search_seconds=search_seconds,
            parameter_relative_error=float(
                np.linalg.norm(rates - truth) / np.linalg.norm(truth)
            ),
            stochastic_vs_ode_normalized_rmse=float(
                np.sqrt(np.mean(((mean - independent(rates)) / SCALE) ** 2))
            ),
            validation_paths=stochastic["completed_replicates"],
            validation_events=stochastic["events"],
        )
        report["arms"].append(row)
        arm_manifest(args.output, strategy, budget, native_calls, search_seconds, 7101)
        report["checks"].append(
            dict(name=strategy + "-finite-fit-and-SSA-conservation", passed=True)
        )
        write(args.output / "summary.json", report)
