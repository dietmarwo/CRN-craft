"""Synthetic global-to-local calibration using native BDF forward sensitivities."""

import argparse
from dataclasses import asdict
import json
from pathlib import Path
import subprocess
import tempfile
import numpy as np
import fcmaes_rust
import sqp_craft


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--binary", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    root = Path(__file__).resolve().parents[2]
    model = json.loads(
        (root / "models/examples/immigration-death.model.json").read_text()
    )
    times = np.linspace(0.0, 4.0, 9)
    target = 3.0 / 0.8 * (-np.expm1(-0.8 * times))
    with tempfile.TemporaryDirectory(prefix="crn-ode-fit-") as tmp:
        folder = Path(tmp)
        request = folder / "ode.json"
        request.write_text(
            json.dumps(
                {
                    "schema": "crn-craft.ode-request/0.1",
                    "times": times.tolist(),
                    "rtol": 1e-10,
                    "atol": 1e-12,
                }
            )
        )

        class Objective:
            def __init__(self):
                self.cached = None
                self.value = None
                self.gradient = None
                self.calls = 0

            def functions(self, x):
                if self.cached != list(x):
                    rates = np.exp(x)
                    model["parameters"].update(
                        {"birth": float(rates[0]), "death": float(rates[1])}
                    )
                    mp = folder / "model.json"
                    mp.write_text(json.dumps(model))
                    proc = subprocess.run(
                        [str(args.binary.resolve()), str(mp), str(request)],
                        capture_output=True,
                        check=False,
                        timeout=30,
                    )
                    if proc.returncode:
                        raise sqp_craft.UndefinedPointError(proc.stderr.decode())
                    result = json.loads(proc.stdout)
                    actual = np.asarray(result["values"])[:, 0]
                    residual = actual - target
                    self.value = float(residual @ residual)
                    sensitivity = np.asarray(result["sensitivities"])[:, :, 0]
                    self.gradient = (2 * sensitivity @ residual * rates).tolist()
                    self.cached = list(x)
                    self.calls += 1
                return self.value, []

            def derivatives(self, x):
                self.functions(x)
                return self.gradient, []

        objective = Objective()
        global_result = fcmaes_rust.optimize_bite(
            lambda x: objective.functions(x)[0],
            [0.5, -0.5],
            [-1.0, -2.0],
            [2.0, 1.0],
            seed=7001,
            max_evaluations=64,
        )
        start = list(global_result[0])
        spec = sqp_craft.ProblemSpec(
            lower=[-1.0, -2.0],
            upper=[2.0, 1.0],
            constraint_lower=[],
            constraint_upper=[],
            jacobian_indptr=[0],
            jacobian_indices=[],
        )
        derivative_report = sqp_craft.check_derivatives(
            spec, objective, start, absolute_tolerance=1e-5, relative_tolerance=1e-4
        )
        local = sqp_craft.solve(
            spec,
            objective,
            start,
            method="sqp",
            method_options=sqp_craft.SqpOptions(
                qp_solver="active-set", max_iterations=60, stationarity_tolerance=1e-7
            ),
        )
        value, _ = objective.functions(local.candidate)
        gradient, _ = objective.derivatives(local.candidate)
        fitted = np.exp(local.candidate)
        independent = 3.0 / 0.8 * (-np.expm1(-0.8 * times))
        check_values = fitted[0] / fitted[1] * (-np.expm1(-fitted[1] * times))
        independent_loss = float(np.sum((check_values - independent) ** 2))
        passed = (
            value < 1e-10
            and independent_loss < 1e-10
            and max(abs(g) for g in gradient) < 1e-5
        )
        report = {
            "schema": "crn-craft.global-local-fit/0.1",
            "passed": passed,
            "versions": {
                "fcmaes-rust": fcmaes_rust.__version__,
                "sqp-craft": sqp_craft.__version__,
                "diffsol": "0.18.0",
            },
            "synthetic_truth": [3.0, 0.8],
            "global": {
                "x": start,
                "loss": float(global_result[1]),
                "evaluations": int(global_result[2]),
            },
            "local": asdict(local),
            "fitted_rates": fitted.tolist(),
            "native_loss": value,
            "independent_analytic_loss": independent_loss,
            "gradient_inf_norm": max(abs(g) for g in gradient),
            "derivative_check": asdict(derivative_report),
            "actual_native_solves_including_checks": objective.calls,
        }
        args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
        if not passed:
            raise RuntimeError("fit gate failed; retained report")
        print(
            json.dumps(
                {
                    "passed": passed,
                    "global_loss": float(global_result[1]),
                    "local_loss": value,
                    "native_solves": objective.calls,
                }
            )
        )


if __name__ == "__main__":
    main()
