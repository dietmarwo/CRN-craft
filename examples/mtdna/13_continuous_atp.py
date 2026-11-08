"""Variable population with ATP as a continuous variable (hybrid method)."""

import json

import numpy as np
from crn_craft import rules

# The ATP model as an ordinary rules-model/0.3 file: ATP is a continuous
# variable with its ODE written out, and replication reads it.
data = rules.load("models/rules/mtdna-three-class-atp-continuous.rules.json")
model = json.loads(data)
print("continuous:", model["continuous"])
print("replication_wt rate:", model["rules"][1]["rate"])
run = rules.ensemble(
    data,
    end=25_920.0,
    report_times=24.0 * np.arange(0, 1081, 180),
    replicates=64,
    seed=7,
    threads=2,
    method="hybrid",  # rtol 1e-6, atol 1e-9 unless a method dict says otherwise
)
result = run.simulate()
meta = result.metadata
print("columns:", list(result.statistics["mean"].columns))
print("continuous columns:", meta["continuous"])
print("sampler:", meta["sampler"], "method:", meta["method"])
stopped = result.records["termination"] == "stopped"
print("stopped (ATP exhausted):", stopped.sum())
mean = result.statistics["mean"]
for hours, row in mean.iterrows():
    print(
        f"day {hours / 24:4.0f}: N {row['N']:7.1f}  "
        f"mutants {row['mutants']:6.1f}  ATP {row['ATP']:7.2f}"
    )

# The same ODE, written with the builder for a single wild-type class.
cell = rules.Model("atp_wild_type")
cell.add_parameter({"f": 0.1, "v1": 0.2, "v2": 0.01, "c": 1.3034153747670567})
cell.add_species("W", 1000)
cell.add_ode("ATP", 450.0, "f * W - v1 * ATP - v2 * W", lower_bound=0)
cell.add_reaction("c / (c + ATP) * W", ["W"], ["W", "W"], name="replication")
cell.add_reaction(0.0028881132523331052, ["W"], [], name="degradation")
small = rules.ensemble(
    cell,
    end=2_400.0,
    report_times=[0.0, 1_200.0, 2_400.0],
    replicates=16,
    seed=3,
    threads=2,
    method={"kind": "hybrid", "rtol": 1e-7, "atol": 1e-9},
).simulate()
print("builder model method:", small.metadata["method"])
print(small.statistics["mean"].round(2))
