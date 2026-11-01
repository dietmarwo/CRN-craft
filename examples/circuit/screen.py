"""ODE prioritization with a mandatory stochastic allocation for every topology.

Declared gate (screen-study/0.2), fixed before the study runs:
1. Every topology receives exactly ALLOCATION completed SSA paths, each with a
   diagnostic, whatever its ODE classification.
2. Every ODE solve completes (no domain or solver error).
3. No topology is rejected because of its ODE classification.
4. The positive control (the three-gene ring) is classified oscillatory by the
   ODE diagnostic of its screened species.
5. At least one topology is not classified oscillatory by the ODE (the Vilar
   noise-driven control) and still receives the full allocation, so the run
   exercises the case the mandatory allocation exists for.
The gate is also applied to perturbed copies of the measured rows; each
perturbation must fail it, or the script fails.
"""

import argparse
import copy
import json
from pathlib import Path
import subprocess
import tempfile
import numpy as np
import crn_craft as crn
from crn_craft.observables import oscillation

ALLOCATION = 8
# (model, screened species, role); species are named, never positional.
TOPOLOGIES = [
    ("repressilator", "m1", "positive_control"),
    ("vilar-noise", "A", "noise_driven_control"),
]


def gate(rows):
    """Criteria 1-5 of the module docstring, each reported separately."""
    criteria = {
        "full_allocation_for_every_topology": all(
            r["stochastic_allocation"] == ALLOCATION
            and r["completed_paths"] == ALLOCATION
            and len(r["stochastic_diagnostics"]) == ALLOCATION
            for r in rows
        ),
        "every_ode_solve_completed": all(
            r["ode_diagnostic"]["classification"] != "domain_or_solver_error"
            for r in rows
        ),
        "no_ode_rejection": not any(r["rejected_by_ode"] for r in rows),
        "positive_control_ode_oscillatory": all(
            r["ode_diagnostic"]["classification"] == "oscillatory"
            for r in rows
            if r["role"] == "positive_control"
        )
        and any(r["role"] == "positive_control" for r in rows),
        "ode_negative_topology_still_allocated": any(
            r["ode_diagnostic"]["classification"] != "oscillatory"
            and r["completed_paths"] == ALLOCATION
            for r in rows
        ),
    }
    return all(criteria.values()), criteria


def perturbations(rows):
    """Copies of the measured rows that each violate one criterion."""

    def edit(change):
        broken = copy.deepcopy(rows)
        change(broken)
        return broken

    def drop_path(rs):
        rs[1]["stochastic_diagnostics"].pop()
        rs[1]["completed_paths"] -= 1

    def ode_failure(rs):
        rs[0]["ode_diagnostic"] = {"classification": "domain_or_solver_error"}

    def reject(rs):
        rs[1]["rejected_by_ode"] = True

    def ambiguous_control(rs):
        rs[0]["ode_diagnostic"]["classification"] = "ambiguous"

    def all_ode_positive(rs):
        for r in rs:
            r["ode_diagnostic"]["classification"] = "oscillatory"

    return {
        "one_path_missing": edit(drop_path),
        "ode_solver_failure": edit(ode_failure),
        "ode_rejection": edit(reject),
        "positive_control_not_oscillatory": edit(ambiguous_control),
        "no_ode_negative_topology": edit(all_ode_positive),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ode", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    # 400 leaves about four ~45-unit ring periods after the 50% burn-in; the
    # detector needs three cycles.
    ap.add_argument("--horizon", type=float, default=400.0)
    args = ap.parse_args()
    root = Path(__file__).resolve().parents[2]
    rows = []
    times = np.linspace(0.0, args.horizon, int(round(2 * args.horizon)) + 1).tolist()
    with tempfile.TemporaryDirectory() as folder:
        folder = Path(folder)
        for name, species, role in TOPOLOGIES:
            model = (root / f"models/examples/{name}.model.json").read_bytes()
            index = [s["name"] for s in json.loads(model)["species"]].index(species)
            (folder / "model.json").write_bytes(model)
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
                oscillation(
                    times, np.asarray(ode["values"])[:, index], burn_fraction=0.5
                )
                if ode
                else {"classification": "domain_or_solver_error"}
            )
            # No hard filter: even an apparently stationary ODE gets held-out SSA paths.
            request = crn.run_spec(
                model,
                end=args.horizon,
                save_times=times,
                replicates=ALLOCATION,
                seed=20261301,
                keep=ALLOCATION,
                max_events=10000000,
            )
            result = crn.simulate(model, request)
            completed = (
                len(result["trajectories"])
                if result["termination"] == "completed"
                else 0
            )
            score = [
                oscillation(times, np.asarray(x["counts"])[:, index], burn_fraction=0.5)
                for x in result["trajectories"]
            ]
            rows.append(
                {
                    "model": name,
                    "species": species,
                    "role": role,
                    "ode_diagnostic": diagnostic,
                    "stochastic_allocation": ALLOCATION,
                    "completed_paths": completed,
                    "stochastic_diagnostics": score,
                    "events": result["events"],
                    "rejected_by_ode": False,
                    "scope": "measured controls, not population-wide detector sensitivity",
                }
            )
    passed, criteria = gate(rows)
    self_test = {key: gate(broken)[0] for key, broken in perturbations(rows).items()}
    report = {
        "schema": "crn-craft.screen-study/0.2",
        "horizon": args.horizon,
        "allocation": ALLOCATION,
        "passed": passed and not any(self_test.values()),
        "gate": criteria,
        "gate_fails_on_perturbed_rows": {k: not v for k, v in self_test.items()},
        "rows": rows,
    }
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    if not report["passed"]:
        raise RuntimeError("screen gate failed; retained report")
    print(json.dumps({"passed": True, "gate": criteria}))


if __name__ == "__main__":
    main()
