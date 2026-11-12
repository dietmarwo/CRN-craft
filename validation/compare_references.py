"""Same-host dynamic-API reference comparison; no RNG/path-identity claim."""

import argparse
import importlib.metadata
import json
import math
from pathlib import Path
import statistics
import time

import numpy as np
import rebop
from benchmark_campaign import models, run_case


def reference(model):
    simulation = rebop.Gillespie()
    for r in model["reactions"]:
        a = [s for s, n in r["reactants"].items() for _ in range(n)]
        b = [s for s, n in r["products"].items() for _ in range(n)]
        rate = r["rate"]
        if rate["kind"] != "mass_action":
            raise ValueError("elementary benchmark only")
        k = rate["coefficient"]
        k = model["parameters"][k] if isinstance(k, str) else k
        # rebop uses falling factorials rather than binomial factors.
        k /= math.prod(math.factorial(n) for n in r["reactants"].values())
        simulation.add_reaction(k, a, b)
    return simulation


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--binary", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--command", default="run")
    p.add_argument("--rebop-only", action="store_true")
    args = p.parse_args()
    report = {
        "schema": "crn-craft.reference-cost/0.1",
        "versions": {"rebop": importlib.metadata.version("rebop")},
        "replicates": 64,
        "horizon": 100.0,
        "seed_policy": "CRN root seed 4242 streams; references per-replicate seeds 4242+i; independent samples",
        "rows": [],
        "unavailable": [
            {
                "name": "Julia JumpProcesses/FiniteStateProjection",
                "reason": "Julia runtime not installed; no measured comparison",
            }
        ],
    }
    try:
        import gillespy2 as gp

        report["versions"]["gillespy2"] = gp.__version__
    except ImportError:
        gp = None
    if args.rebop_only:
        gp = None
    for name, model in list(models())[:2]:
        reb = reference(model)

        def run_reb(reb=reb, model=model):
            for i in range(64):
                reb.run(
                    {s["name"]: s["initial_count"] for s in model["species"]},
                    100.0,
                    1,
                    rng=4242 + i,
                )

        run_reb()  # warmup
        values = []
        direct = []
        gilles = []
        gmodel = None
        if gp:
            gmodel = gp.Model(name="reference", volume=1.0)
            gs = {
                s["name"]: gp.Species(
                    name=s["name"], initial_value=s["initial_count"], mode="discrete"
                )
                for s in model["species"]
            }
            gmodel.add_species(list(gs.values()))
            for i, r in enumerate(model["reactions"]):
                k = gp.Parameter(name=f"k{i}", expression=r["rate"]["coefficient"])
                gmodel.add_parameter(k)
                gmodel.add_reaction(
                    gp.Reaction(
                        name=r["name"],
                        reactants={gs[s]: n for s, n in r["reactants"].items()},
                        products={gs[s]: n for s, n in r["products"].items()},
                        rate=k,
                    )
                )
            gmodel.timespan(np.array([0.0, 100.0]))
        for _ in range(5):
            d = run_case(args.binary.resolve(), args.command, model, 1)
            direct.append(d["median_process_seconds"])
            t = time.perf_counter()
            run_reb()
            values.append(time.perf_counter() - t)
            if gp:
                t = time.perf_counter()
                gmodel.run(
                    solver=gp.NumPySSASolver, number_of_trajectories=64, seed=4242
                )
                gilles.append(time.perf_counter() - t)
        report["rows"].append(
            {
                "workload": name,
                "CRN_process_seconds": direct,
                "rebop_dynamic_python_seconds": values,
                "rebop_over_CRN": statistics.median(values) / statistics.median(direct),
                "gillespy2_NumPySSA_seconds": gilles,
                "work": "64 complete trajectories, [0,100] output, independent reference event counts not instrumented",
                "comparison_limit": "reference includes Python/xarray overhead; CRN is the wall time of a crn-craft run process (start, run with statistics, printing)",
            }
        )
        args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
        print(name, report["rows"][-1]["rebop_over_CRN"], flush=True)


if __name__ == "__main__":
    main()
