"""Exact event records: the event log and live callbacks."""

import json
from collections import Counter

import numpy as np
from crn_craft import rules

model = json.loads(rules.load("models/rules/mtdna-three-class.rules.json"))
names = [rule["name"] for rule in model["rules"]]

# 1. A bounded log of selected rules for every cell.
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
total = Counter()
for observation in result.observations:
    for entry in observation["log"]:
        total[names[entry["rule"]]] += 1
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

# 2. A live callback sees every committed event in order; returning False
#    stops that replicate. Callbacks hold the GIL, so use them on few cells.
seen = Counter()


def watch(replicate, event):
    seen[event["rule"]] += 1
    # Stop a cell as soon as an advantaged mutant arises.
    return event["rule"] != names.index("mutation_mtAdv")


run = rules.ensemble(
    model,
    end=25_920.0,
    report_times=[0.0, 25_920.0],
    replicates=20,
    seed=4,
    threads=1,
    method="static_direct",
)
result = run.simulate(callback=watch)
stopped = result.records["termination"] == rules.TERMINATIONS.index("stopped")
print("events seen by the callback:", sum(seen.values()))
print("cells stopped at their first advantaged mutation:", stopped.sum())
print("their stop days:", np.round(result.records["stop_time"][stopped] / 24, 1))
