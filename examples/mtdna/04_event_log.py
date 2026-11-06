"""Exact event records: a bounded event log for every cell."""

import json
from collections import Counter

import numpy as np
from crn_craft import rules

model = json.loads(rules.load("models/rules/mtdna-three-class.rules.json"))
run = rules.ensemble(
    model,
    end=25_920.0,
    report_times=[0.0, 25_920.0],
    replicates=100,
    seed=4,
    threads=2,
    method="static_direct",
    event_log={"max_entries": 50, "rules": ["mutation_mt", "mutation_mtAdv"]},
)
result = run.simulate()
# Log entries carry a rule index; result.rule_names turns it into a name.
total = Counter()
for observation in result.observations:
    for entry in observation["log"]:
        total[result.rule_names[entry["rule"]]] += 1
print("logged mutation events:", dict(total))
print("omitted (log full):", sum(o["omitted_log_entries"] for o in result.observations))
columns = result.event_columns(result.observations[0]["id"])
print(
    "cell",
    result.observations[0]["id"],
    "mutation days:",
    np.round(columns["time"] / 24, 1),
)
observables = result.metadata["observables"]
print("N after each logged event:", columns["observables"][:, observables.index("N")])
