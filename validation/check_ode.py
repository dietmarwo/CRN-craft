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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--binary", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    root = Path(__file__).resolve().parents[1]
    rows = []
    with tempfile.TemporaryDirectory(prefix="crn-ode-") as tmp:
        folder = Path(tmp)
        r = folder / "ode.json"
        times = [0.0, 0.5, 1.0, 2.0, 5.0]
        r.write_text(
            json.dumps(
                {
                    "schema": "crn-craft.ode-request/0.1",
                    "times": times,
                    "rtol": 1e-10,
                    "atol": 1e-12,
                }
            )
        )
        for name in ["immigration-death", "enzyme-mechanism", "repressilator"]:
            model = json.loads(
                (root / f"models/examples/{name}.model.json").read_text()
            )
            m = folder / "model.json"
            m.write_text(json.dumps(model))
            proc = subprocess.run(
                [str(args.binary.resolve()), str(m), str(r)],
                capture_output=True,
                check=False,
                timeout=60,
            )
            require(proc.returncode == 0, proc.stderr.decode())
            result = json.loads(proc.stdout)
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
