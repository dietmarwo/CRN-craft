"""The same 64 cells under every method that runs the three-class model."""

import json
import time

from crn_craft import rules

model = json.loads(rules.load("models/rules/mtdna-three-class.rules.json"))
model["stop"] = "N > 0 && mutants >= 0.9 * N"
methods = {
    "tree": None,
    "static_direct": "static_direct",
    "static_direct, ziggurat": {"kind": "static_direct", "waiting_time": "ziggurat"},
    "next_reaction": "next_reaction",
    "tau (epsilon 0.03)": {"kind": "tau", "epsilon": 0.03},
    "binomial (1 h steps)": {"kind": "binomial", "step": 1.0},
    "mtdna_direct (opt-in)": "mtdna_direct",
}
for label, method in methods.items():
    run = rules.ensemble(
        model,
        end=25_920.0,
        report_times=[0.0, 12_960.0, 25_920.0],
        replicates=64,
        seed=5,
        threads=2,
        method=method,
    )
    start = time.perf_counter()
    result = run.simulate()
    seconds = time.perf_counter() - start
    meta = result.metadata
    stopped = result.records["termination"] == "stopped"
    final = result.statistics["mean"].iloc[-1]  # the means at the last report time
    fraction = final["mutants"] / final["N"]
    print(
        f"{label:24s} {seconds:6.2f} s  {meta['events']:>11,} events  "
        f"taken over {stopped.sum():2d}  final mutant fraction {fraction:.3f}  "
        f"{meta['sampler']}"
    )
