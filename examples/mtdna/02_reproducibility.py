"""Seeds, threads and reproducibility."""

import json

import numpy as np
from crn_craft import rules

model = json.loads(rules.load("models/rules/mtdna-three-class.rules.json"))
model["stop"] = "N > 0 && mutants >= 0.9 * N"


def run(seed, threads):
    ensemble = rules.ensemble(
        model,
        end=25_920.0,
        report_times=[0.0, 8640.0, 17_280.0, 25_920.0],
        replicates=64,
        seed=seed,
        threads=threads,
        method="mtdna_direct",
        keep_paths=64,
    )
    return ensemble.simulate()


one, two, other = run(7, 1), run(7, 2), run(8, 2)
same = all(np.array_equal(one.records[k], two.records[k]) for k in one.records)
print("1 vs 2 threads, same seed, identical records:", same)
print(
    "same statistics:", np.array_equal(one.statistics["mean"], two.statistics["mean"])
)
print(
    "seed 7 vs 8, identical stop times:",
    np.array_equal(one.records["stop_time"], other.records["stop_time"]),
)
# Replicate i always uses the same stream, whatever the ensemble size.
small = rules.ensemble(
    model,
    end=25_920.0,
    report_times=[0.0, 8640.0, 17_280.0, 25_920.0],
    replicates=8,
    seed=7,
    threads=1,
    method="mtdna_direct",
).simulate()
print(
    "first 8 of 64 equal an 8-cell run:",
    np.array_equal(small.records["events"], one.records["events"][:8]),
)
