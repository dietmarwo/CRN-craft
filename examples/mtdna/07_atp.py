"""Variable population with the ATP balance (exact thinning)."""

import json

import numpy as np
from crn_craft import rules

model = json.loads(rules.load("models/rules/mtdna-three-class-atp.rules.json"))
print("ATP block:", model["atp"])
run = rules.ensemble(
    model,
    end=25_920.0,
    report_times=24.0 * np.arange(0, 1081, 180),
    replicates=64,
    seed=7,
    threads=2,
    method="thinning",
)
result = run.simulate()
meta = result.metadata
names = meta["observables"]
print("observables:", names)
print("sampler:", meta["sampler"], "stopped:", meta["stopped_replicates"])
mean = result.statistics["mean"]
for day, row in zip(np.arange(0, 1081, 180), mean, strict=True):
    print(
        f"day {day:4d}: N {row[names.index('N')]:7.1f}  "
        f"mutants {row[names.index('mutants')]:6.1f}  ATP {row[names.index('ATP')]:7.2f}"
    )
