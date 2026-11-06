"""Exact event records: a bounded event log for every cell."""

import json

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
# Every logged event, indexed by (replicate, entry); rule_name names the rule.
log = result.event_columns()
counts = log.groupby("rule_name", observed=True).size()
print("logged mutation events:", counts.to_dict())
print("omitted (log full):", sum(o["omitted_log_entries"] for o in result.observations))
cell = log.loc[0]  # the log of cell 0
print(cell[["time", "rule", "rule_name", "N", "mutants"]])
print("cell 0 mutation days:", (cell["time"] / 24).round(1).tolist())
