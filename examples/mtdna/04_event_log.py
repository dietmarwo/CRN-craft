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
# Every logged event as a DataArray: cells x entries x observables after the
# event, NaN-padded; time, rule and rule_name are coordinates of each entry.
log = result.event_columns()
names = log["rule_name"].to_series().dropna()  # padding has no name
print("logged mutation events:", names.value_counts().sort_index().to_dict())
print("omitted (log full):", sum(o["omitted_log_entries"] for o in result.observations))
table = result.event_frame(0)  # the log of cell 0 as a DataFrame
print(table[["time", "rule", "rule_name", "N", "mutants"]])
days = result.event_columns(0)["time"] / 24  # cell 0, without padding
print("cell 0 mutation days:", days.round(1).values.tolist())
