"""Budgets, wall time and cancellation."""

import json

from crn_craft import rules

model = json.loads(rules.load("models/rules/mtdna-three-class.rules.json"))


def summary(result):
    terms = sorted(set(result.records["termination"]))
    codes = [d["code"] for d in result.metadata["diagnostics"]]
    return result.metadata["termination"], terms, codes


# Too few events for three years: every cell ends with event_limit and the
# ensemble withholds its statistics.
small = rules.ensemble(
    model,
    end=25_920.0,
    report_times=[0.0, 25_920.0],
    replicates=4,
    seed=1,
    max_events=10_000,
    method="mtdna_direct",
)
result = small.simulate()
print(summary(result), "statistics:", result.statistics)

# A wall-time budget stops the ensemble cooperatively.
slow = rules.ensemble(
    model,
    end=25_920.0,
    report_times=[0.0, 25_920.0],
    replicates=20_000,
    seed=1,
    threads=2,
    wall_time_s=1.0,
    method="mtdna_direct",
)
print(summary(slow.simulate()))
