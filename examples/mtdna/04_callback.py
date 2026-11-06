"""Exact event records: a live callback that stops cells."""

import json
from collections import Counter

from crn_craft import rules

model = json.loads(rules.load("models/rules/mtdna-three-class.rules.json"))
run = rules.ensemble(
    model,
    end=25_920.0,
    report_times=[0.0, 25_920.0],
    replicates=20,
    seed=4,
    threads=1,
    method="static_direct",
)
seen = Counter()


def watch(replicate, event):
    """Count events by rule; stop a cell at its first advantaged mutant."""
    name = run.rule_names[event["rule"]]
    seen[name] += 1
    return name != "mutation_mtAdv"


# Callbacks hold the GIL, so use them on few cells.
result = run.simulate(callback=watch)
stopped = result.records["termination"] == "stopped"
print("events seen by the callback:", sum(seen.values()))
print("cells stopped at their first advantaged mutation:", stopped.sum())
days = result.records.loc[stopped, "stop_time"] / 24
print("their stop days:", days.round(1).to_numpy())
