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
records = result.records
stopped = records["termination"] == rules.TERMINATIONS.index("stopped")
print(meta["sampler"], meta["termination"])
print(f"{meta['events']:,} events in {seconds:.2f} s")
print(f"taken over: {stopped.sum()} of {len(stopped)} cells")
print("takeover days (first five):", np.sort(records["stop_time"][stopped])[:5] / 24)
names = meta["observables"]
mean = result.statistics["mean"]
fraction = mean[:, names.index("mutants")] / mean[:, names.index("N")]
for day, value in zip(days, fraction, strict=True):
    print(f"day {day:5d}: mean mutant fraction {value:.4f}")
