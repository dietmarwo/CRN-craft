"""Per-lineage fates with the clone model (tree method)."""

import json

import numpy as np
from crn_craft import rules

# Every new mutant is its own clone (attribute `clone` from new_id()).
model = json.loads(rules.load("models/rules/mtdna-fixed.rules.json"))
model["parameters"].update({"K": 100, "w0": 100, "mu": 0.01, "delta": 0.05})
run = rules.ensemble(
    model,
    end=400.0,
    report_times=[0.0, 100.0, 200.0, 400.0],
    replicates=20,
    seed=6,
    threads=2,
    keep_instances=3,
    trackers=[{"observable": "mutant_clones"}],
    event_log={"max_entries": 5000, "rules": []},
    max_output_bytes=1 << 30,
)
result = run.simulate()
print("sampler:", result.metadata["sampler"])
mean = result.statistics["mean"]
names = result.metadata["observables"]
print("mean number of distinct mutant clones:", mean[:, names.index("mutant_clones")])
print(
    "peak instances (clones + wild type), max over cells:",
    result.records["peak_instances"].max(),
)
print("instances created per cell:", result.records["created_instances"].mean())
omitted = sum(o["omitted_log_entries"] for o in result.observations)
print("log entries omitted:", omitted)  # 0: every event is in the log

# Lineage lifetimes from the full event log: replay each cell's net changes
# per clone id (reactants are consumed, products created).
born, died, peak = {}, {}, {}
for observation in result.observations:
    counts = {}
    for entry in observation["log"]:
        change = {}
        for sign, records in ((-1, entry["reactants"]), (1, entry["products"])):
            for record in records:
                clone, _, wt = record["attributes"]
                if wt == 0:  # mutant clones only
                    key = (observation["id"], clone)
                    change[key] = change.get(key, 0) + sign * record["count"]
        for key, delta in change.items():
            counts[key] = counts.get(key, 0) + delta
            born.setdefault(key, entry["time"])
            peak[key] = max(peak.get(key, 0), counts[key])
            if counts[key] == 0:
                died[key] = entry["time"]
life = np.array([died[k] - born[k] for k in born if k in died])
print(f"{len(born)} clones born, {len(life)} lost within the horizon")
print("median lifetime of lost clones (h):", np.median(life))
print("largest clone ever (molecules):", max(peak.values()))
print("final instances of cell 0:", result.metadata["final_instances"][0])
