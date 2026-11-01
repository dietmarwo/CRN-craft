"""ODE prioritization with a mandatory stochastic allocation for every topology."""

import argparse
import json
from pathlib import Path
import subprocess
import tempfile
import numpy as np
import crn_craft as crn
from crn_craft.observables import oscillation


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ode", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    root = Path(__file__).resolve().parents[2]
    rows = []
    with tempfile.TemporaryDirectory() as folder:
        folder = Path(folder)
        for name in ["repressilator", "vilar-noise"]:
            model = (root / f"models/examples/{name}.model.json").read_bytes()
            (folder / "model.json").write_bytes(model)
            times = np.linspace(0.0, 200.0, 401).tolist()
            (folder / "ode.json").write_text(
                json.dumps(
                    {
                        "schema": "crn-craft.ode-request/0.1",
                        "times": times,
                        "rtol": 1e-8,
                        "atol": 1e-10,
                    }
                )
            )
            p = subprocess.run(
                [
                    str(args.ode.resolve()),
                    str(folder / "model.json"),
                    str(folder / "ode.json"),
                ],
                capture_output=True,
                check=False,
                timeout=120,
            )
            ode = json.loads(p.stdout) if not p.returncode else None
            diagnostic = (
                oscillation(times, np.asarray(ode["values"])[:, 0], burn_fraction=0.5)
                if ode
                else {"classification": "domain_or_solver_error"}
            )
            # No hard filter: even an apparently stationary ODE gets held-out SSA paths.
            request = crn.run_spec(
                model,
                end=200.0,
                save_times=times,
                replicates=8,
                seed=20261301,
                keep=8,
                max_events=10000000,
            )
            result = crn.simulate(model, request)
            if result["termination"] != "completed":
                raise RuntimeError("stochastic screen incomplete")
            score = [
                oscillation(times, np.asarray(x["counts"])[:, 0], burn_fraction=0.5)
                for x in result["trajectories"]
            ]
            rows.append(
                {
                    "model": name,
                    "ode_diagnostic": diagnostic,
                    "stochastic_allocation": 8,
                    "stochastic_diagnostics": score,
                    "events": result["events"],
                    "rejected_by_ode": False,
                    "scope": "measured controls, not population-wide detector sensitivity",
                }
            )
    args.output.write_text(
        json.dumps(
            {
                "schema": "crn-craft.screen-study/0.1",
                "passed": all(
                    r["stochastic_allocation"] > 0 and not r["rejected_by_ode"]
                    for r in rows
                ),
                "rows": rows,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
