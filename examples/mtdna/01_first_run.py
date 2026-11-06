"""Three-class model: 200 cells for three years with the takeover stop."""

import json
import time

import numpy as np
from crn_craft import rules

model = json.loads(rules.load("models/rules/mtdna-three-class.rules.json"))
model["stop"] = "N > 0 && mutants >= 0.9 * N"
days = np.arange(0, 1081, 90)
run = rules.ensemble(
    model,
    end=25_920.0,
    report_times=24.0 * days,
    replicates=200,
    seed=1,
    threads=2,
    method="mtdna_direct",
)
start = time.perf_counter()
result = run.simulate()
seconds = time.perf_counter() - start

meta = result.metadata
records = result.records  # one row per cell
stopped = records["termination"] == "stopped"
print(meta["sampler"], meta["termination"])
print(f"{meta['events']:,} events in {seconds:.2f} s")
print(f"taken over: {stopped.sum()} of {len(records)} cells")
print(records[["termination", "stop_time", "events"]].head(3))
first = records.loc[stopped, "stop_time"].nsmallest(5)
print("takeover days (first five):", first.to_numpy() / 24)
mean = result.statistics["mean"]  # one row per report time
fraction = mean["mutants"] / mean["N"]
for hours, value in fraction.items():
    print(f"day {hours / 24:5.0f}: mean mutant fraction {value:.4f}")
