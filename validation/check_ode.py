"""Independent linear, nonlinear and forward sensitivity checks of diffsol."""

import argparse
import json
import math
from pathlib import Path
import subprocess
import tempfile
import numpy as np
from scipy.integrate import solve_ivp


def require(value, message):
    if not value:
        raise RuntimeError(message)


def rhs(model, y):
    names = [s["name"] for s in model["species"]]
    x = dict(zip(names, y, strict=True))
    p = model["parameters"]
    out = np.zeros(len(y))

    def resolve(k):
        return p[k] if isinstance(k, str) else k

    for rx in model["reactions"]:
        rate = rx["rate"]
        kind = rate["kind"]
        if kind == "mass_action":
            a = resolve(rate["coefficient"]) * math.prod(
                x[s] ** n / math.factorial(n) for s, n in rx["reactants"].items()
            )
        elif kind == "michaelis_menten":
            z = x[rate["substrate"]]
            a = resolve(rate["vmax"]) * z / (resolve(rate["km"]) + z)
        else:
            z = x[rate["regulator"]]
            power = resolve(rate["n"])
            k = resolve(rate["k"])
            f = z**power / (k**power + z**power)
            a = resolve(rate["basal"]) + resolve(rate["vmax"]) * (
                1 - f if kind == "hill_repression" else f
            )
        out += a * np.array(
            [rx["products"].get(s, 0) - rx["reactants"].get(s, 0) for s in names]
        )
    return out


def reference(model, times):
    y0 = [s["initial_count"] for s in model["species"]]
    return solve_ivp(
        lambda t, y, model=model: rhs(model, y),
        [0.0, times[-1]],
        y0,
        method="DOP853",
        t_eval=times,
        rtol=2e-13,
        atol=1e-13,
    ).y.T


def scaled_error(actual, expected):
    actual = np.asarray(actual, dtype=float)
    expected = np.asarray(expected, dtype=float)
    return float(np.max(np.abs(actual - expected) / np.maximum(1.0, np.abs(expected))))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--binary", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    root = Path(__file__).resolve().parents[1]
    fixtures = root / "integrations/diffsol/fixtures"
    rows = []
    with tempfile.TemporaryDirectory(prefix="crn-ode-") as tmp:
        folder = Path(tmp)

        def run(model, times, rtol=1e-10, atol=1e-12):
            """Exit code and parsed stdout result or one-line stderr diagnostic."""
            m = folder / "model.json"
            r = folder / "request.json"
            m.write_text(json.dumps(model))
            r.write_text(
                json.dumps(
                    {
                        "schema": "crn-craft.ode-request/0.1",
                        "times": times,
                        "rtol": rtol,
                        "atol": atol,
                    }
                )
            )
            proc = subprocess.run(
                [str(args.binary.resolve()), str(m), str(r)],
                capture_output=True,
                check=False,
                timeout=60,
            )
            if proc.returncode == 0:
                return 0, json.loads(proc.stdout)
            lines = proc.stderr.decode().splitlines()
            # Anything but one JSON diagnostic line (for example a panic) fails its row.
            try:
                outcome = json.loads(lines[0]) if len(lines) == 1 else None
            except ValueError:
                outcome = None
            if proc.stdout or not isinstance(outcome, dict):
                return proc.returncode, {"malformed_failure_output": True}
            return proc.returncode, outcome

        times = [0.0, 0.5, 1.0, 2.0, 5.0]
        # The original 16 predeclared comparisons, unchanged.
        for name in ["immigration-death", "enzyme-mechanism", "repressilator"]:
            model = json.loads((root / f"models/legacy/{name}.model.json").read_text())
            code, result = run(model, times)
            require(code == 0, json.dumps(result))
            y0 = [s["initial_count"] for s in model["species"]]
            ref = solve_ivp(
                lambda t, y, model=model: rhs(model, y),
                [0.0, times[-1]],
                y0,
                method="DOP853",
                t_eval=times,
                rtol=2e-13,
                atol=1e-13,
            )
            actual = np.array(result["values"])
            expected = ref.y.T
            error = float(
                np.max(np.abs(actual - expected) / np.maximum(1.0, np.abs(expected)))
            )
            rows.append(
                {
                    "model": name,
                    "check": "independent_DOP853",
                    "error": error,
                    "passed": error <= 1e-5,
                }
            )
            for j, parameter in enumerate(result["parameters"]):
                original = model["parameters"][parameter]
                h = 1e-4 * max(1.0, original)
                if original <= h:
                    continue
                model["parameters"][parameter] = original + h
                a = solve_ivp(
                    lambda t, y, model=model: rhs(model, y),
                    [0.0, times[-1]],
                    y0,
                    method="DOP853",
                    t_eval=times,
                    rtol=2e-13,
                    atol=1e-13,
                ).y.T
                model["parameters"][parameter] = original - h
                b = solve_ivp(
                    lambda t, y, model=model: rhs(model, y),
                    [0.0, times[-1]],
                    y0,
                    method="DOP853",
                    t_eval=times,
                    rtol=2e-13,
                    atol=1e-13,
                ).y.T
                model["parameters"][parameter] = original
                numeric = (a - b) / (2 * h)
                sensitivity = np.array(result["sensitivities"][j])
                error = float(
                    np.max(
                        np.abs(sensitivity - numeric) / np.maximum(1.0, np.abs(numeric))
                    )
                )
                rows.append(
                    {
                        "model": name,
                        "parameter": parameter,
                        "check": "forward_sensitivity",
                        "error": error,
                        "passed": error <= 1e-4,
                    }
                )
            if name == "immigration-death":
                birth = model["parameters"]["birth"]
                death = model["parameters"]["death"]
                exact = np.array(
                    [birth / death * (-math.expm1(-death * t)) for t in times]
                )
                rows.append(
                    {
                        "model": name,
                        "check": "analytic_linear",
                        "error": float(np.max(np.abs(actual[:, 0] - exact))),
                        "passed": bool(np.max(np.abs(actual[:, 0] - exact)) <= 1e-7),
                    }
                )

        # Order-2 mass action, Hill activation/repression and a mixed X + 2Y binding
        # with a moving regulator: DOP853 values and relative-step central differences
        # of every raw parameter (the same 1e-5 and 1e-4 thresholds as above).
        fixture_times = [0.0, 0.5, 1.0, 2.0, 5.0, 10.0]
        for name in ["dimerization", "hill-constant", "hill-dynamic"]:
            model = json.loads((fixtures / f"{name}.model.json").read_text())
            code, result = run(model, fixture_times)
            if code != 0:
                rows.append(
                    {
                        "model": name,
                        "check": "solve_completed",
                        "exit": code,
                        "code": result.get("code"),
                        "passed": False,
                    }
                )
                continue
            actual = np.array(result["values"])
            error = scaled_error(actual, reference(model, fixture_times))
            rows.append(
                {
                    "model": name,
                    "check": "independent_DOP853",
                    "error": error,
                    "passed": error <= 1e-5,
                }
            )
            require(
                result["parameters"] == sorted(model["parameters"]),
                "sensitivity parameter order",
            )
            for j, parameter in enumerate(result["parameters"]):
                original = model["parameters"][parameter]
                h = 1e-4 * original
                model["parameters"][parameter] = original + h
                a = reference(model, fixture_times)
                model["parameters"][parameter] = original - h
                b = reference(model, fixture_times)
                model["parameters"][parameter] = original
                error = scaled_error(result["sensitivities"][j], (a - b) / (2 * h))
                rows.append(
                    {
                        "model": name,
                        "parameter": parameter,
                        "check": "forward_sensitivity_relative_step",
                        "error": error,
                        "passed": error <= 1e-4,
                    }
                )
            if name == "dimerization":
                # dA/dt = -c A^2 (rate c A^2/2!, two A per event); A(1) = 50.
                c = model["parameters"]["c"]
                a0 = model["species"][0]["initial_count"]
                t = np.array(fixture_times)
                exact = a0 / (1 + c * a0 * t)
                ds = -(a0**2) * t / (1 + c * a0 * t) ** 2
                checks = {
                    "analytic_dimerization_A1_equals_50": abs(actual[2, 0] - 50.0),
                    "analytic_dimerization_values": scaled_error(
                        actual, np.stack([exact, (a0 - exact) / 2], axis=1)
                    ),
                    "analytic_dimerization_sensitivity": scaled_error(
                        result["sensitivities"][0], np.stack([ds, -ds / 2], axis=1)
                    ),
                }
                for check, error in checks.items():
                    rows.append(
                        {
                            "model": name,
                            "check": check,
                            "error": float(error),
                            "passed": bool(error <= 1e-6),
                        }
                    )
            if name == "hill-constant":
                p = model["parameters"]
                v, b, k, n, g = (p[q] for q in ["v", "b", "K", "n", "g"])
                regulator = model["species"][0]["initial_count"]
                f = 1 / (1 + (k / regulator) ** n)
                t = np.array(fixture_times)
                shape = -np.expm1(-g * t) / g
                order = result["parameters"]
                for i, kind, sign, fraction in [
                    (1, "activation", 1.0, f),
                    (2, "repression", -1.0, 1 - f),
                ]:
                    a = b + v * fraction
                    exact = {
                        "v": fraction * shape,
                        "b": shape,
                        "K": sign * v * (-n / k * f * (1 - f)) * shape,
                        "n": sign * v * f * (1 - f) * math.log(regulator / k) * shape,
                        "g": -a / g * shape + a * t * np.exp(-g * t) / g,
                    }
                    sens = np.array(result["sensitivities"])[:, :, i]
                    error = max(
                        scaled_error(actual[:, i], a * shape),
                        *(
                            scaled_error(sens[order.index(q)], exact[q])
                            for q in sorted(exact)
                        ),
                    )
                    rows.append(
                        {
                            "model": name,
                            "check": f"analytic_hill_{kind}",
                            "error": error,
                            "passed": error <= 1e-6,
                        }
                    )

        # Decay through zero: the value band is one error weight, atol + rtol * A0.
        decay = {
            "schema": "crn-craft.model/0.1",
            "name": "decay",
            "convention": "stochastic_counts",
            "omega": 1.0,
            "parameters": {"k": 1.0},
            "species": [{"name": "A", "initial_count": 10}],
            "reactions": [
                {
                    "name": "d",
                    "reactants": {"A": 1},
                    "products": {},
                    "rate": {"kind": "mass_action", "coefficient": "k"},
                }
            ],
        }
        for horizon in [[0.0, 20.0, 40.0], [0.0, 100.0, 1000.0], [0.0, 1000.0, 1e4]]:
            for rtol, atol in [(1e-10, 1e-12), (1e-6, 1e-8)]:
                code, result = run(decay, horizon, rtol, atol)
                error = (
                    float(
                        np.max(
                            np.abs(
                                np.array(result["values"])[:, 0]
                                - 10 * np.exp(-np.array(horizon))
                            )
                        )
                    )
                    if code == 0
                    else None
                )
                rows.append(
                    {
                        "model": "decay",
                        "check": "analytic_decay_through_zero",
                        "times": horizon,
                        "rtol": rtol,
                        "atol": atol,
                        "exit": code,
                        "error": error,
                        "passed": code == 0 and error <= atol + rtol * 10,
                    }
                )
        # Long horizons that formerly failed on BDF trial states below zero.
        for name, horizon, conserved in [
            ("enzyme-mechanism", 1000.0, [["E", "ES"], ["S", "ES", "P"]]),
            ("michaelis-menten", 200.0, [["S", "P"]]),
        ]:
            model = json.loads((root / f"models/legacy/{name}.model.json").read_text())
            grid = [0.0, horizon / 4, horizon / 2, horizon]
            code, result = run(model, grid)
            error = drift = None
            names = [s["name"] for s in model["species"]]
            y0 = np.array([s["initial_count"] for s in model["species"]], dtype=float)
            if code == 0:
                actual = np.array(result["values"])
                drift = max(
                    float(
                        np.max(
                            np.abs(
                                actual[:, [names.index(s) for s in group]].sum(axis=1)
                                - y0[[names.index(s) for s in group]].sum()
                            )
                        )
                    )
                    for group in conserved
                )
                error = scaled_error(actual, reference(model, grid))
            rows.append(
                {
                    "model": name,
                    "check": "long_horizon_DOP853_and_conservation",
                    "horizon": horizon,
                    "exit": code,
                    "error": error,
                    "conservation_drift": drift,
                    "passed": code == 0 and error <= 1e-5 and drift <= 1e-8,
                }
            )

        # Failure contract: one JSON line on stderr, exit 1, a specific code/message.
        immigration = json.loads(
            (root / "models/legacy/immigration-death.model.json").read_text()
        )
        limited = json.loads(json.dumps(immigration))
        limited["species"][0]["population_limit"] = 5
        chain = {
            "schema": "crn-craft.model/0.1",
            "name": "chain",
            "convention": "stochastic_counts",
            "omega": 1.0,
            "parameters": {"k": 1.0},
            "species": [
                {"name": "A", "initial_count": 100},
                {"name": "B", "initial_count": 0, "population_limit": 30},
            ],
            "reactions": [
                {
                    "name": "ab",
                    "reactants": {"A": 1},
                    "products": {"B": 1},
                    "rate": {"kind": "mass_action", "coefficient": "k"},
                },
                {
                    "name": "b",
                    "reactants": {"B": 1},
                    "products": {},
                    "rate": {"kind": "mass_action", "coefficient": "k"},
                },
            ],
        }
        singular = {
            "schema": "crn-craft.model/0.1",
            "name": "singular",
            "convention": "stochastic_counts",
            "omega": 1.0,
            "parameters": {"v": 2.0},
            "species": [
                {"name": "X", "initial_count": 0},
                {"name": "Y", "initial_count": 1},
            ],
            "reactions": [
                {
                    "name": "h",
                    "reactants": {},
                    "products": {"Y": 1},
                    "rate": {
                        "kind": "hill_activation",
                        "regulator": "X",
                        "vmax": "v",
                        "basal": 0.0,
                        "k": 5.0,
                        "n": 0.5,
                    },
                }
            ],
        }
        for check, model, grid, code_expected, text in [
            (
                "population_limit_at_save_time",
                limited,
                [0.0, 1.0, 2.0],
                "CRN-E044",
                "population_limit 5",
            ),
            (
                "population_limit_between_save_times",
                chain,
                [0.0, 10.0],
                "CRN-E044",
                "species B",
            ),
            (
                "huge_time_no_panic",
                immigration,
                [0.0, 1.7976931348623157e308],
                "CRN-E080",
                "",
            ),
            (
                "kinetics_diagnostic_retained",
                singular,
                [0.0, 1.0],
                "CRN-E080",
                "singular at zero",
            ),
        ]:
            code, outcome = run(model, grid, 1e-8, 1e-10)
            rows.append(
                {
                    "check": check,
                    "exit": code,
                    "code": outcome.get("code"),
                    "passed": code == 1
                    and outcome.get("success") is False
                    and outcome.get("code") == code_expected
                    and text in outcome.get("message", ""),
                }
            )
    report = {
        "schema": "crn-craft.ode-validation/0.1",
        "passed": all(r["passed"] for r in rows),
        "rows": rows,
    }
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    require(report["passed"], "ODE numerical gate failed; retained report")
    print(json.dumps({"passed": True, "checks": len(rows)}))


if __name__ == "__main__":
    main()
