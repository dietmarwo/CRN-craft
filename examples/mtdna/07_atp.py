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
print("observables:", meta["observables"])
print("sampler:", meta["sampler"], "stopped:", meta["stopped_replicates"])
mean = result.statistics["mean"]
for hours, row in mean.iterrows():
    print(
        f"day {hours / 24:4.0f}: N {row['N']:7.1f}  "
        f"mutants {row['mutants']:6.1f}  ATP {row['ATP']:7.2f}"
    )
