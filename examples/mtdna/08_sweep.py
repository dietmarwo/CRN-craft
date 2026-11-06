"""A parameter sweep that reuses one prepared template."""

import json
import time

from crn_craft import rules

model = json.loads(rules.load("models/rules/mtdna-three-class.rules.json"))
model["stop"] = "N > 0 && mutants >= 0.9 * N"
template = rules.ensemble(
    model,
    end=25_920.0,
    report_times=[0.0, 25_920.0],
    replicates=400,
    seed=8,
    threads=2,
    method="mtdna_direct",
)
start = time.perf_counter()
print("takeover fraction after three years")
print("   mu \\ sm " + "".join(f"{sm:8.2f}" for sm in (1.0, 1.25, 1.5, 2.0)))
for mu in (1e-5, 3e-5, 1e-4):
    row = []
    for sm in (1.0, 1.25, 1.5, 2.0):
        result = template.bind_parameters({"mu": mu, "sm": sm}).simulate()
        row.append((result.records["termination"] == "stopped").mean())
    print(f"{mu:9.0e}  " + "".join(f"{v:8.3f}" for v in row))
print(f"12 candidates in {time.perf_counter() - start:.1f} s")
# The same seed gives every candidate the same random streams (common
# random numbers), so differences between cells of the table are not noise
# from different seeds alone.
